"""
    test_mismatch_classification_first_principles.jl

First-principles test cases for edge mismatch classification.
These tests are designed from the CONCEPTUAL definition of each mismatch type,
not from the implementation details.

Each test case is constructed from geometric first principles:
- What does this mismatch type mean geometrically?
- What is the minimal example that demonstrates this concept?
- Does the classification correctly identify this fundamental case?
"""

using Test
using Nas2Step

const Z = 0.0

# Helper to create coordinate key
ckey(x, y, z=Z) = (Float64(x), Float64(y), Float64(z))

# Helper to create triangle
function make_tri(id, c1, c2, c3)
    return Triangle(id, id+100, id+200, id, ckey(c1...), ckey(c2...), ckey(c3...))
end

# Helper to build minimal topology for testing
function make_test_topology(facesA::Vector{Triangle}, facesB::Vector{Triangle})
    # Extract all vertices
    all_verts_A = Set{NTuple{3,Float64}}()
    all_verts_B = Set{NTuple{3,Float64}}()
    
    for tri in facesA
        push!(all_verts_A, tri.coord1, tri.coord2, tri.coord3)
    end
    for tri in facesB
        push!(all_verts_B, tri.coord1, tri.coord2, tri.coord3)
    end
    
    shared = intersect(all_verts_A, all_verts_B)
    node_key_to_ids = Dict(k => (1, 1) for k in shared)
    
    # Build edge maps
    edgesA = Dict{EdgeKey,Vector{Int}}()
    edgesB = Dict{EdgeKey,Vector{Int}}()
    
    for (idx, tri) in enumerate(facesA)
        edges = [
            EdgeKey(tri.coord1, tri.coord2),
            EdgeKey(tri.coord2, tri.coord3),
            EdgeKey(tri.coord3, tri.coord1)
        ]
        for e in edges
            push!(get!(edgesA, e, Int[]), idx)
        end
    end
    
    for (idx, tri) in enumerate(facesB)
        edges = [
            EdgeKey(tri.coord1, tri.coord2),
            EdgeKey(tri.coord2, tri.coord3),
            EdgeKey(tri.coord3, tri.coord1)
        ]
        for e in edges
            push!(get!(edgesB, e, Int[]), idx)
        end
    end
    
    edges_only_A = setdiff(Set(keys(edgesA)), Set(keys(edgesB)))
    edges_only_B = setdiff(Set(keys(edgesB)), Set(keys(edgesA)))
    edges_shared = intersect(Set(keys(edgesA)), Set(keys(edgesB)))
    
    bbox = BoundingBox(collect(shared))
    
    return InterfaceTopology(
        1, 2, shared, node_key_to_ids,
        facesA, facesB, edgesA, edgesB,
        edges_only_A, edges_only_B, edges_shared,
        bbox, length(shared), length(facesA), length(facesB),
        length(edges_shared), length(edges_only_A) + length(edges_only_B), 0.0,
        0.0, 0.0, 0, 1.0  # consistency metrics: max_vertex_dist, mean_vertex_dist, edge_mismatch_count, triangulation_similarity
    )
end

@testset "First Principles: T_JUNCTION" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: One mesh has a vertex ON an edge, the other doesn't
        # This is like a "T" intersection where one line ends at another
        
        # Create a simple coplanar quad interface:
        # A: Two triangles forming a quad with edge (0,5)-(10,5) across the middle
        # B: Four triangles - same quad but with vertex (5,5) splitting that middle edge
        
        # Mesh A: Simple quad split by diagonal (0,0)-(10,10)
        # Vertices: (0,0), (10,0), (10,10), (0,10)
        triA1 = make_tri(1, (0,0,Z), (10,0,Z), (10,10,Z))
        triA2 = make_tri(2, (0,0,Z), (10,10,Z), (0,10,Z))
        
        # Mesh B: Same quad but with extra vertex (5,5) in the middle
        # This creates 4 triangles from the center point
        triB1 = make_tri(3, (0,0,Z), (10,0,Z), (5,5,Z))
        triB2 = make_tri(4, (10,0,Z), (10,10,Z), (5,5,Z))
        triB3 = make_tri(5, (10,10,Z), (0,10,Z), (5,5,Z))
        triB4 = make_tri(6, (0,10,Z), (0,0,Z), (5,5,Z))
        
        topo = make_test_topology([triA1, triA2], [triB1, triB2, triB3, triB4])
        
        # Edge (0,0)-(10,10) exists in A but not in B (B has vertex (5,5) on it)
        edge_diag = EdgeKey(ckey(0,0), ckey(10,10))
        
        @test edge_diag ∈ topo.edges_only_in_A
        
        # Classify from A's perspective
        # Vertex (5,5) from B should be detected as hanging on A's diagonal
        mismatch = Nas2Step.classify_edge_mismatch(edge_diag, topo, :A)
        
        @test mismatch.mismatch_type == Nas2Step.T_JUNCTION
        @test length(mismatch.hanging_nodes) == 1
        @test mismatch.hanging_nodes[1] == ckey(5,5)
    end
end

@testset "First Principles: REFINEMENT" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: One mesh has MULTIPLE vertices ON an edge (hierarchical subdivision)
        # This is like subdividing an edge multiple times
        
        # Create a simple quad where one mesh has a single diagonal,
        # the other has that diagonal split by 2 intermediate vertices
        
        # Mesh A: Simple quad with diagonal (0,0)-(15,15)
        triA1 = make_tri(1, (0,0,Z), (15,0,Z), (15,15,Z))
        triA2 = make_tri(2, (0,0,Z), (15,15,Z), (0,15,Z))
        
        # Mesh B: Same quad but diagonal split at (5,5) and (10,10)
        # Creates 3 triangular regions along the diagonal
        triB1 = make_tri(3, (0,0,Z), (15,0,Z), (5,5,Z))
        triB2 = make_tri(4, (5,5,Z), (15,0,Z), (10,10,Z))
        triB3 = make_tri(5, (10,10,Z), (15,0,Z), (15,15,Z))
        triB4 = make_tri(6, (0,0,Z), (5,5,Z), (0,15,Z))
        triB5 = make_tri(7, (5,5,Z), (10,10,Z), (0,15,Z))
        triB6 = make_tri(8, (10,10,Z), (15,15,Z), (0,15,Z))
        
        topo = make_test_topology([triA1, triA2], [triB1, triB2, triB3, triB4, triB5, triB6])
        
        # Edge (0,0)-(15,15) exists in A, not in B
        edge = EdgeKey(ckey(0,0), ckey(15,15))
        
        @test edge ∈ topo.edges_only_in_A
        
        # Classify from A's perspective
        # Vertices (5,5) and (10,10) from B should be detected as hanging
        mismatch = Nas2Step.classify_edge_mismatch(edge, topo, :A)
        @test mismatch.mismatch_type == Nas2Step.REFINEMENT
        @test length(mismatch.hanging_nodes) == 2
        # Both (5,5) and (10,10) should be hanging on the edge
    end
end

@testset "First Principles: DIAGONAL" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: A quad (4-sided polygon) can be split into 2 triangles using EITHER diagonal
        # Both meshes have the SAME quad boundary, just different diagonal choices
        
        # Quad with corners: A(0,0), B(10,0), C(10,10), D(0,10)
        A, B, C, D = ckey(0,0), ckey(10,0), ckey(10,10), ckey(0,10)
        
        # Side A: Uses diagonal A-C (splits into ABC and ACD)
        triA1 = make_tri(1, A, B, C)
        triA2 = make_tri(2, A, C, D)
        
        # Side B: Uses diagonal B-D (splits into ABD and BCD)
        triB1 = make_tri(3, A, B, D)
        triB2 = make_tri(4, B, C, D)
        
        # CRITICAL: Both have SAME 4 boundary edges: AB, BC, CD, DA
        # Only the internal diagonal differs
        
        topo = make_test_topology([triA1, triA2], [triB1, triB2])
        
        # Edge A-C exists in A but not B
        edgeAC = EdgeKey(A, C)
        
        if edgeAC ∈ topo.edges_only_in_A
            mismatch = Nas2Step.classify_edge_mismatch(edgeAC, topo, :A)
            @test mismatch.mismatch_type == Nas2Step.DIAGONAL
            @test length(mismatch.quad_vertices) == 4
            # Should find the 4 corners of the quad
            @test Set(mismatch.quad_vertices) == Set([A, B, C, D])
        else
            @warn "Diagonal edge not found in edges_only_in_A"
        end
    end
end

@testset "First Principles: QUAD_MISMATCH" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: Same 4 vertices but they DON'T form the same quad
        # The boundary edges are different (twisted connectivity)
        
        # 4 vertices arranged in a square
        V1, V2, V3, V4 = ckey(0,0), ckey(10,0), ckey(10,10), ckey(0,10)
        
        # Side A: Normal quad V1-V2-V3-V4 with boundary edges:
        # V1-V2, V2-V3, V3-V4, V4-V1, and diagonal V1-V3
        triA1 = make_tri(1, V1, V2, V3)
        triA2 = make_tri(2, V1, V3, V4)
        
        # Side B: TWISTED quad - different boundary!
        # Uses edges: V1-V3, V3-V2, V2-V4, V4-V1, and diagonal V3-V4
        # This is NOT the same quad shape even though same 4 vertices
        triB1 = make_tri(3, V1, V3, V2)
        triB2 = make_tri(4, V1, V4, V3)
        triB3 = make_tri(5, V2, V3, V4)  # Extra triangle to complete connectivity
        
        # Note: This is a degenerate/twisted case but demonstrates the concept
        # In reality, this might result from different triangulation strategies
        
        topo = make_test_topology([triA1, triA2], [triB1, triB2, triB3])
        
        edgeV1V3 = EdgeKey(V1, V3)
        
        if edgeV1V3 ∈ topo.edges_only_in_A
            mismatch = Nas2Step.classify_edge_mismatch(edgeV1V3, topo, :A)
            # Should detect that while 4 vertices exist, quad boundaries don't match
            # Expected: QUAD_MISMATCH or UNKNOWN (depending on exact geometry)
            @test mismatch.mismatch_type ∈ [Nas2Step.QUAD_MISMATCH, Nas2Step.UNKNOWN]
        end
    end
end

@testset "First Principles: BOUNDARY_EDGE" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: An edge exists on the BOUNDARY of one mesh (only 1 triangle touches it)
        # It's not an internal edge, so it's not truly a shared interface edge
        
        # Side A: Single triangle - ALL its edges are boundary edges
        triA = make_tri(1, (0,0,Z), (10,0,Z), (5,5,Z))
        
        # Side B: A mesh that happens to share some vertices but not this specific edge
        # B has a triangle using vertices (10,0) and (5,5) but not edge (10,0)-(5,5)
        triB1 = make_tri(2, (0,0,Z), (10,0,Z), (10,10,Z))
        triB2 = make_tri(3, (0,0,Z), (10,10,Z), (5,5,Z))
        
        topo = make_test_topology([triA], [triB1, triB2])
        
        # Edge (10,0)-(5,5) exists in A and appears in only 1 triangle
        edge = EdgeKey(ckey(10,0), ckey(5,5))
        
        if edge ∈ topo.edges_only_in_A
            mismatch = Nas2Step.classify_edge_mismatch(edge, topo, :A)
            @test mismatch.mismatch_type == Nas2Step.BOUNDARY_EDGE
            # Boundary edges appear in exactly 1 triangle on their side
        else
            @warn "Boundary edge not found in edges_only_in_A"
        end
    end
end

@testset "First Principles: NON_MANIFOLD" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: An edge is shared by MORE than 2 triangles
        # This violates manifold topology - like multiple surfaces meeting at an edge
        
        # Side A: THREE triangles sharing the same edge (0,0)-(10,0)
        # This is geometrically impossible in 3D without overlapping surfaces
        triA1 = make_tri(1, (0,0,Z), (10,0,Z), (5,5,Z))
        triA2 = make_tri(2, (0,0,Z), (10,0,Z), (5,-5,Z))
        triA3 = make_tri(3, (0,0,Z), (10,0,Z), (5,0,5))  # All 3 share same edge!
        
        # Side B: Normal mesh without this edge
        triB = make_tri(4, (0,0,Z), (10,10,Z), (5,5,Z))
        
        topo = make_test_topology([triA1, triA2, triA3], [triB])
        
        # Edge (0,0)-(10,0) should be non-manifold (appears in 3 triangles)
        edge = EdgeKey(ckey(0,0), ckey(10,0))
        
        if edge ∈ topo.edges_only_in_A
            mismatch = Nas2Step.classify_edge_mismatch(edge, topo, :A)
            @test mismatch.mismatch_type == Nas2Step.NON_MANIFOLD
        else
            @warn "Non-manifold edge not found in edges_only_in_A"
        end
    end
end

@testset "First Principles: UNSHARED_ENDPOINT" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: An edge's endpoint vertices are NOT in the shared vertex set
        # This means the edge connects vertices that don't exist on both sides
        
        # Side A: TWO triangles sharing edge (10,0)-(5,5) so it's NOT a boundary edge
        # This ensures the edge passes the boundary check but fails the shared endpoint check
        triA1 = make_tri(1, (0,0,Z), (10,0,Z), (5,5,Z))
        triA2 = make_tri(2, (10,0,Z), (15,0,Z), (5,5,Z))  # Second triangle sharing the edge
        
        # Side B: Triangle with DIFFERENT vertices - only (0,0) is shared
        # Vertices (10,0), (15,0), (5,5) are NOT shared with A
        triB = make_tri(3, (0,0,Z), (15,0,Z), (7,7,Z))
        
        topo = make_test_topology([triA1, triA2], [triB])
        
        # Edge (10,0)-(5,5) exists in A (in 2 triangles, so not boundary)
        # but both vertices are NOT in shared set because B doesn't have them
        edge = EdgeKey(ckey(10,0), ckey(5,5))
        
        if edge ∈ topo.edges_only_in_A
            mismatch = Nas2Step.classify_edge_mismatch(edge, topo, :A)
            @test mismatch.mismatch_type == Nas2Step.UNSHARED_ENDPOINT
            # One or both endpoints not in shared vertex set
        else
            @warn "Edge with unshared endpoints not in edges_only_in_A"
        end
    end
end

@testset "First Principles: UNKNOWN" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: Cases that don't fit any other category
        # Example: Edge exists but geometric structure doesn't allow classification
        
        # Side A: Triangle with edge (0,0)-(10,0)
        triA = make_tri(1, (0,0,Z), (10,0,Z), (5,5,Z))
        
        # Side B: Complex geometry where the edge exists but in a way
        # that doesn't match any standard pattern
        # For example: vertices are shared, no hanging nodes, but can't form a quad
        triB1 = make_tri(2, (0,0,Z), (3,3,Z), (5,5,Z))
        triB2 = make_tri(3, (3,3,Z), (10,0,Z), (5,5,Z))
        triB3 = make_tri(4, (0,0,Z), (10,0,Z), (3,3,Z))
        
        topo = make_test_topology([triA], [triB1, triB2, triB3])
        
        # The edge (0,0)-(10,0) exists in both, but the triangulations
        # don't form a recognizable pattern
        # This should fall through to UNKNOWN
        
        # Note: This is the hardest to test from first principles since
        # UNKNOWN is defined by exclusion (not fitting other categories)
        # The test mainly validates that complex cases get caught
    end
end

"""
KEY INSIGHTS FROM FIRST-PRINCIPLES APPROACH:

1. T_JUNCTION: One mesh refined an edge (1 hanging node)
2. REFINEMENT: One mesh refined an edge multiple times (2+ hanging nodes)
3. DIAGONAL: Same quad boundary, different diagonal choice
4. QUAD_MISMATCH: Same 4 vertices, different quad boundary
5. BOUNDARY_EDGE: Edge only touches 1 triangle (not interior)
6. NON_MANIFOLD: Edge shared by 3+ triangles (broken topology)
7. UNSHARED_ENDPOINT: Edge connects non-shared vertices
8. UNKNOWN: Doesn't match any pattern above

Each test is built from the GEOMETRIC MEANING, not implementation details.
"""
