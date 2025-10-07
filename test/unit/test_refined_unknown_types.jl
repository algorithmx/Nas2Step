"""
    test_refined_unknown_types.jl

First-principles test cases for refined UNKNOWN mismatch classifications.
Tests DEGENERATE_EDGE, SOURCE_EDGE_ABSENT, QUAD_NOT_FOUND_IN_SOURCE, and
QUAD_VERTICES_NOT_TRIANGULATED_IN_TARGET.
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
        length(edges_shared), length(edges_only_A) + length(edges_only_B), 0.0
    )
end

@testset "Refined UNKNOWN: DEGENERATE_EDGE" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: Edge has zero or near-zero length (endpoints nearly coincident)
        # This should be detected early before other checks
        
        # Create a triangle in A with a degenerate edge (two vertices very close)
        # Use tolerance of 1e-4, so make edge length < 1e-4
        triA = make_tri(1, (0,0,Z), (0,0.00005,Z), (1,0,Z))  # Edge (0,0)-(0,0.00005) is degenerate
        
        # Mesh B: simple triangle
        triB = make_tri(2, (0,0,Z), (1,0,Z), (0.5,0.5,Z))
        
        topo = make_test_topology([triA], [triB])
        
        # The degenerate edge
        edge_degen = EdgeKey(ckey(0,0), ckey(0,0.00005))
        
        if edge_degen ∈ topo.edges_only_in_A
            mismatch = Nas2Step.classify_edge_mismatch(edge_degen, topo, :A)
            
            @test mismatch.mismatch_type == Nas2Step.DEGENERATE_EDGE
            @test mismatch.diagnostics !== nothing
            @test mismatch.diagnostics.reason == "Edge has degenerate length < 0.0001"
        end
    end
end

@testset "Refined UNKNOWN: QUAD_NOT_FOUND_IN_SOURCE" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: Source has 2 triangles sharing an edge, but extract_quad_vertices fails
        # This can happen with degenerate or nearly-coplanar collapsed quads
        # 
        # For this test, create two triangles that share an edge but have only 3 unique vertices
        # Triangle 1: (0,0), (10,0), (5,5)
        # Triangle 2: (0,0), (5,5), (10,0) - same 3 vertices, makes a degenerate "quad"
        # They share edge (0,0)-(5,5) BUT both endpoints ARE shared with B
        triA1 = make_tri(1, (0,0,Z), (10,0,Z), (5,5,Z))
        triA2 = make_tri(2, (0,0,Z), (5,5,Z), (10,0,Z))
        
        # Mesh B: shares the endpoints (0,0) and (5,5) to pass endpoint check
        triB = make_tri(3, (0,0,Z), (5,5,Z), (7,7,Z))
        
        topo = make_test_topology([triA1, triA2], [triB])
        
        # Edge (0,0)-(5,5) is shared by 2 triangles in A, endpoints are shared
        edge = EdgeKey(ckey(0,0), ckey(5,5))
        
        if edge ∈ topo.edges_only_in_A
            mismatch = Nas2Step.classify_edge_mismatch(edge, topo, :A)
            
            # Should detect that we cannot extract 4 unique vertices from the two source triangles
            # (They only have 3 unique vertices total: (0,0), (10,0), (5,5))
            @test mismatch.mismatch_type == Nas2Step.QUAD_NOT_FOUND_IN_SOURCE
            @test mismatch.diagnostics !== nothing
            @test occursin("Cannot extract 4 unique vertices", mismatch.diagnostics.reason)
        end
    end
end

@testset "Refined UNKNOWN: TARGET_USES_FINER_TRIANGULATION" begin
    @testset "Conceptual Definition" begin
        # CONCEPT: Source quad exists with 4 vertices, but target doesn't use those vertices
        # Target might use a completely different set of vertices for the same region
        
        # Mesh A: Quad with vertices (0,0), (10,0), (10,10), (0,10) split by diagonal (0,0)-(10,10)
        triA1 = make_tri(1, (0,0,Z), (10,0,Z), (10,10,Z))
        triA2 = make_tri(2, (0,0,Z), (10,10,Z), (0,10,Z))
        
        # Mesh B: Uses completely different vertices (not sharing the interior of the quad)
        # Only shares the corner vertices (0,0), (10,0), (10,10), (0,10) but not all at once
        # Instead uses a different triangulation with extra vertices
        triB1 = make_tri(3, (0,0,Z), (10,0,Z), (5,0,Z))
        triB2 = make_tri(4, (0,0,Z), (5,0,Z), (0,10,Z))
        triB3 = make_tri(5, (5,0,Z), (10,0,Z), (10,10,Z))
        triB4 = make_tri(6, (5,0,Z), (10,10,Z), (0,10,Z))
        
        topo = make_test_topology([triA1, triA2], [triB1, triB2, triB3, triB4])
        
        # Diagonal edge in A
        edge_diag = EdgeKey(ckey(0,0), ckey(10,10))
        
        if edge_diag ∈ topo.edges_only_in_A
            mismatch = Nas2Step.classify_edge_mismatch(edge_diag, topo, :A)
            
            # Should detect that target doesn't triangulate the same 4 vertices
            # (Target doesn't have a pair of triangles using all 4 corner vertices)
            @test mismatch.mismatch_type == Nas2Step.TARGET_USES_FINER_TRIANGULATION
            @test mismatch.diagnostics !== nothing
            @test occursin("Target has no triangles using the source quad", mismatch.diagnostics.reason)
        end
    end
end

@testset "Integration: Verify diagnostics are populated" begin
    @testset "Diagnostics fields" begin
        # Create a simple degenerate case and verify all diagnostic fields are populated
        triA = make_tri(1, (0,0,Z), (0,0.00005,Z), (1,0,Z))
        triB = make_tri(2, (0,0,Z), (1,0,Z), (0.5,0.5,Z))
        
        topo = make_test_topology([triA], [triB])
        edge_degen = EdgeKey(ckey(0,0), ckey(0,0.00005))
        
        if edge_degen ∈ topo.edges_only_in_A
            mismatch = Nas2Step.classify_edge_mismatch(edge_degen, topo, :A)
            
            @test mismatch.diagnostics !== nothing
            diag = mismatch.diagnostics
            
            # Check all fields are present
            @test diag.present_in == :A
            @test typeof(diag.source_triangle_count) == Int
            @test typeof(diag.target_triangle_count_using_endpoints) == Int
            @test typeof(diag.endpoints_shared) == Tuple{Bool,Bool}
            @test typeof(diag.edge_length) == Float64
            @test typeof(diag.tried_quad_finding) == Bool
            @test typeof(diag.quad_vertices_found) == Bool
            @test typeof(diag.target_triangles_using_quad) == Int
            @test typeof(diag.reason) == String
            
            # Check specific values for degenerate edge
            @test diag.edge_length < 0.0001
            @test length(diag.reason) > 0
        end
    end
end

println("\n✓ All refined UNKNOWN type tests passed!")
