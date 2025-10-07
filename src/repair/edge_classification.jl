# Edge mismatch classification for interface repair
# Phase 1.2: Classify edge mismatches into repair types

# Note: geometric_utilities.jl is included by the parent module

# ============================================================================
# Edge mismatch types and classification
# ============================================================================

@enum MismatchType begin
    T_JUNCTION          # One region has intermediate node, other doesn't
    DIAGONAL            # Same VALID quad, different triangulation (strict check)
    REFINEMENT          # Hierarchical mesh refinement difference
    QUAD_MISMATCH       # Same 4 vertices but different quad boundary topology
    BOUNDARY_EDGE       # Edge on mesh boundary, not true interface edge
    NON_MANIFOLD        # Edge shared by != 2 triangles (topology issue)
    UNSHARED_ENDPOINT   # One or both endpoints not in shared vertex set
    DEGENERATE_EDGE     # Edge has zero or near-zero length
    SOURCE_EDGE_ABSENT  # Edge not found in any source triangle (data anomaly)
    QUAD_NOT_FOUND_IN_SOURCE           # Cannot extract 4 unique vertices from source
    TARGET_USES_FINER_TRIANGULATION    # Target boundary uses additional vertices beyond source quad vertices
    UNKNOWN             # Cannot classify (catch-all)
end

"""
Detailed diagnostics for UNKNOWN mismatch cases.
"""
struct UnknownDiagnostics
    # Edge and topology context
    present_in::Symbol                 # :A or :B (which side has the edge)
    source_triangle_count::Int         # multiplicity in source
    target_triangle_count_using_endpoints::Int
    endpoints_shared::Tuple{Bool,Bool} # (node1_shared, node2_shared)
    edge_length::Float64
    
    # Quad-finding context
    tried_quad_finding::Bool
    quad_vertices_found::Bool
    target_triangles_using_quad::Int
    boundary_matched::Union{Bool,Nothing}
    
    # Reason summary
    reason::String
end

"""
Edge mismatch with complete geometric context for repair planning.
"""
struct EdgeMismatch
    edge_key::EdgeKey
    mismatch_type::MismatchType
    
    # Context
    present_in::Symbol  # :A_only or :B_only
    should_be_in::Int   # Which PID needs the edge inserted
    
    # Geometric details
    hanging_nodes::Vector{NTuple{3,Float64}}  # nodes lying on this edge
    affected_triangles::Vector{Int}            # triangle indices that need splitting
    
    # Quad retriangulation (for DIAGONAL type)
    quad_vertices::Vector{NTuple{3,Float64}}   # 4 vertices forming the quad (empty if not diagonal)
    triangles_to_replace::Vector{Int}          # 2 triangle indices to replace (empty if not diagonal)
    
    # Repair complexity (0.0 = simple, 1.0 = complex)
    complexity_score::Float64
    
    # Quality assessment
    min_affected_triangle_quality::Float64
    repair_feasible::Bool
    
    # Diagnostics (optional, populated for UNKNOWN and refined unknown types)
    diagnostics::Union{UnknownDiagnostics,Nothing}
end

# Backward-compatibility convenience constructor (without diagnostics)
# Older tests/builders construct EdgeMismatch without the optional diagnostics field.
# This method forwards to the full constructor with `nothing` diagnostics.
function EdgeMismatch(
    edge_key::EdgeKey,
    mismatch_type::MismatchType,
    present_in::Symbol,
    should_be_in::Int,
    hanging_nodes::Vector{NTuple{3,Float64}},
    affected_triangles::Vector{Int},
    quad_vertices::Vector{NTuple{3,Float64}},
    triangles_to_replace::Vector{Int},
    complexity_score::Float64,
    min_affected_triangle_quality::Float64,
    repair_feasible::Bool
)
    return EdgeMismatch(
        edge_key,
        mismatch_type,
        present_in,
        should_be_in,
        hanging_nodes,
        affected_triangles,
        quad_vertices,
        triangles_to_replace,
        complexity_score,
        min_affected_triangle_quality,
        repair_feasible,
        nothing,
    )
end

"""
Complete classification report for an interface.
"""
struct InterfaceClassification
    topology::InterfaceTopology
    
    # Classified mismatches
    mismatches_A::Vector{EdgeMismatch}  # Edges missing in A
    mismatches_B::Vector{EdgeMismatch}  # Edges missing in B
    
    # Statistics by type
    t_junction_count::Int
    diagonal_count::Int
    refinement_count::Int
    quad_mismatch_count::Int
    boundary_edge_count::Int
    non_manifold_count::Int
    unshared_endpoint_count::Int
    degenerate_edge_count::Int
    source_edge_absent_count::Int
    quad_not_found_in_source_count::Int
    target_uses_finer_triangulation_count::Int
    unknown_count::Int
    
    # Repair assessment
    total_feasible::Int
    total_infeasible::Int
    average_complexity::Float64
end

# Backward-compatibility constructor (legacy short signature)
# InterfaceClassification(topology, mismatches_A, mismatches_B,
#   t_junction_count, diagonal_count, refinement_count, quad_mismatch_count,
#   boundary_edge_count, non_manifold_count, average_complexity)
function InterfaceClassification(
    topology::InterfaceTopology,
    mismatches_A::Vector{EdgeMismatch},
    mismatches_B::Vector{EdgeMismatch},
    t_junction_count::Int,
    diagonal_count::Int,
    refinement_count::Int,
    quad_mismatch_count::Int,
    boundary_edge_count::Int,
    non_manifold_count::Int,
    average_complexity::Float64
)
    all_mismatches = vcat(mismatches_A, mismatches_B)
    feasible_count = count(m -> m.repair_feasible, all_mismatches)
    infeasible_count = length(all_mismatches) - feasible_count
    # Fill missing categories with zeros for legacy construction
    unshared_endpoint_count = 0
    degenerate_edge_count = 0
    source_edge_absent_count = 0
    quad_not_found_in_source_count = 0
    target_uses_finer_triangulation_count = 0
    unknown_count = 0
    return InterfaceClassification(
        topology,
        mismatches_A,
        mismatches_B,
        t_junction_count,
        diagonal_count,
        refinement_count,
        quad_mismatch_count,
        boundary_edge_count,
        non_manifold_count,
        unshared_endpoint_count,
        degenerate_edge_count,
        source_edge_absent_count,
        quad_not_found_in_source_count,
        target_uses_finer_triangulation_count,
        unknown_count,
        feasible_count,
        infeasible_count,
        average_complexity,
    )
end

# ============================================================================
# Geometric utilities
# ============================================================================

"""
    point_on_segment(p, a, b; tol=1e-4)

Check if point p lies on segment from a to b (within tolerance).
Returns (is_on_segment, parameter_t, distance²).
"""
function point_on_segment(p::NTuple{3,Float64}, 
                         a::NTuple{3,Float64}, 
                         b::NTuple{3,Float64}; 
                         tol::Real=1e-4)
    
    ab = (b[1] - a[1], b[2] - a[2], b[3] - a[3])
    ap = (p[1] - a[1], p[2] - a[2], p[3] - a[3])
    
    # Dot products
    ab_dot_ab = ab[1]^2 + ab[2]^2 + ab[3]^2
    ap_dot_ab = ap[1]*ab[1] + ap[2]*ab[2] + ap[3]*ab[3]
    
    if ab_dot_ab < (tol^2)
        # Degenerate edge
        dist2 = ap[1]^2 + ap[2]^2 + ap[3]^2
        return (dist2 <= tol^2, 0.0, dist2)
    end
    
    # Parameter along segment
    t = clamp(ap_dot_ab / ab_dot_ab, 0.0, 1.0)
    
    # Closest point on segment
    proj = (a[1] + t*ab[1], a[2] + t*ab[2], a[3] + t*ab[3])
    
    # Distance²
    dp = (p[1] - proj[1], p[2] - proj[2], p[3] - proj[3])
    dist2 = dp[1]^2 + dp[2]^2 + dp[3]^2
    
    # Check if interior point (not endpoint)
    is_interior = (t > 1e-6) && (t < 1.0 - 1e-6)
    is_on = (dist2 <= tol^2) && is_interior
    
    return (is_on, t, dist2)
end

"""
    find_hanging_nodes_on_edge(edge, nodes, tol=1e-4)

Find all nodes that lie on the given edge (excluding endpoints).
"""
function find_hanging_nodes_on_edge(edge::EdgeKey, 
                                   nodes::Set{NTuple{3,Float64}}; 
                                   tol::Real=1e-4)::Vector{NTuple{3,Float64}}
    
    hanging = NTuple{3,Float64}[]
    
    for node in nodes
        # Skip if node is an endpoint
        # NOTE: EdgeKey coordinates are rounded to 4 digits during topology construction,
        # but Triangle coordinates are not. We need to check both:
        # 1) Exact match (for perfect coordinate matches)
        # 2) Coordinate rounding match (for coordinates that round to the same value)
        
        # Helper: round coordinate to 4 digits (matching EdgeKey rounding)
        round_coord(c) = (round(c[1], digits=4), round(c[2], digits=4), round(c[3], digits=4))
        
        is_endpoint = (node == edge.node1) || (node == edge.node2) ||
                     (round_coord(node) == edge.node1) || (round_coord(node) == edge.node2)
        
        if is_endpoint
            continue
        end
        
        is_on, t, _ = point_on_segment(node, edge.node1, edge.node2, tol=tol)
        if is_on
            push!(hanging, node)
        end
    end
    
    # Sort by parameter t for consistent ordering
    if !isempty(hanging)
        sort!(hanging, by = n -> begin
            _, t, _ = point_on_segment(n, edge.node1, edge.node2, tol=tol)
            t
        end)
    end
    
    return hanging
end

"""
    compute_triangle_quality(tri::Triangle)

Compute quality metric for a triangle (min angle / 60°).
Returns 1.0 for equilateral, lower for poor quality.
"""
function compute_triangle_quality(tri::Triangle)::Float64
    # Compute edge lengths
    e1 = sqrt((tri.coord2[1] - tri.coord1[1])^2 + 
              (tri.coord2[2] - tri.coord1[2])^2 + 
              (tri.coord2[3] - tri.coord1[3])^2)
    e2 = sqrt((tri.coord3[1] - tri.coord2[1])^2 + 
              (tri.coord3[2] - tri.coord2[2])^2 + 
              (tri.coord3[3] - tri.coord2[3])^2)
    e3 = sqrt((tri.coord1[1] - tri.coord3[1])^2 + 
              (tri.coord1[2] - tri.coord3[2])^2 + 
              (tri.coord1[3] - tri.coord3[3])^2)
    
    # Compute angles using law of cosines
    function angle_from_sides(a, b, c)
        cos_angle = (a^2 + b^2 - c^2) / (2*a*b + 1e-14)
        cos_angle = clamp(cos_angle, -1.0, 1.0)
        return acos(cos_angle) * 180.0 / π
    end
    
    α1 = angle_from_sides(e1, e3, e2)
    α2 = angle_from_sides(e1, e2, e3)
    α3 = angle_from_sides(e2, e3, e1)
    
    min_angle = min(α1, α2, α3)
    
    # Normalize to [0, 1] where 1.0 = 60° (equilateral)
    return min_angle / 60.0
end

# ============================================================================
# Helper functions for tolerance-based comparisons
# ============================================================================

"""
    are_nodes_equal(node1, node2; tol=1e-4)

Check if two nodes are equal within a given tolerance.
"""
function are_nodes_equal(node1::NTuple{3,Float64}, 
                         node2::NTuple{3,Float64}; 
                         tol::Real=1e-4)::Bool
    dx = node1[1] - node2[1]
    dy = node1[2] - node2[2]
    dz = node1[3] - node2[3]
    dist2 = dx*dx + dy*dy + dz*dz
    return dist2 <= tol*tol
end

"""
    get_triangle_nodes(tri::Triangle)

Extract the three nodes from a triangle as a vector.
"""
function get_triangle_nodes(tri::Triangle)::Vector{NTuple{3,Float64}}
    return [tri.coord1, tri.coord2, tri.coord3]
end

"""
    triangle_has_node(tri::Triangle, node; tol=1e-4)

Check if a triangle contains a node within a given tolerance.
"""
function triangle_has_node(tri::Triangle, 
                          node::NTuple{3,Float64}; 
                          tol::Real=1e-4)::Bool
    return are_nodes_equal(tri.coord1, node, tol=tol) ||
           are_nodes_equal(tri.coord2, node, tol=tol) ||
           are_nodes_equal(tri.coord3, node, tol=tol)
end

# ============================================================================
# Quad validation and boundary edge extraction
# ============================================================================

"""
    extract_quad_boundary_edges(tri1, tri2, diagonal; tol=1e-4)

Extract the 4 boundary edges forming the quad perimeter (excluding the diagonal).
Two triangles sharing a diagonal form a quad with 4 boundary edges.

Returns: Vector{EdgeKey} of the 4 boundary edges, or empty if invalid.
"""
function extract_quad_boundary_edges(
    tri1::Triangle,
    tri2::Triangle,
    diagonal::EdgeKey;
    tol::Real=1e-4
)::Vector{EdgeKey}
    # Get all edges from both triangles
    function get_triangle_edges(tri::Triangle)::Vector{EdgeKey}
        return [
            EdgeKey(tri.coord1, tri.coord2),
            EdgeKey(tri.coord2, tri.coord3),
            EdgeKey(tri.coord3, tri.coord1)
        ]
    end
    
    edges1 = get_triangle_edges(tri1)
    edges2 = get_triangle_edges(tri2)
    all_edges = vcat(edges1, edges2)
    
    # Find edges that appear only once (boundary edges)
    # The diagonal appears in both triangles, so it will appear twice
    edge_counts = Dict{EdgeKey,Int}()
    for edge in all_edges
        # Normalize edge orientation for counting
        normalized = edge
        edge_counts[normalized] = get(edge_counts, normalized, 0) + 1
    end
    
    # Boundary edges appear exactly once
    boundary_edges = EdgeKey[]
    for (edge, count) in edge_counts
        if count == 1
            push!(boundary_edges, edge)
        end
    end
    
    # A valid quad should have exactly 4 boundary edges
    if length(boundary_edges) != 4
        return EdgeKey[]
    end
    
    return boundary_edges
end

"""
    edges_match(edges1, edges2; tol=1e-4)

Check if two sets of edges match (same edges, possibly different order).
"""
function edges_match(
    edges1::Vector{EdgeKey},
    edges2::Vector{EdgeKey};
    tol::Real=1e-4
)::Bool
    if length(edges1) != length(edges2)
        return false
    end
    
    # Check if every edge in edges1 has a match in edges2
    for e1 in edges1
        found = false
        for e2 in edges2
            # Check if edges match (either direction)
            if (are_nodes_equal(e1.node1, e2.node1, tol=tol) && 
                are_nodes_equal(e1.node2, e2.node2, tol=tol)) ||
               (are_nodes_equal(e1.node1, e2.node2, tol=tol) && 
                are_nodes_equal(e1.node2, e2.node1, tol=tol))
                found = true
                break
            end
        end
        if !found
            return false
        end
    end
    
    return true
end

# ============================================================================
# Quad finding for diagonal mismatches - SOURCE-FIRST APPROACH
# ============================================================================
# The key insight: When an edge exists in mesh A but not in mesh B, we should:
# 1. Find the quad structure in mesh A (where the edge EXISTS as a diagonal)
# 2. Extract the 4 vertices forming this quad
# 3. Find how those same 4 vertices are triangulated in mesh B
# This is much more robust than searching for a quad structure in mesh B where
# the edge doesn't exist and the triangulation is different!
# ============================================================================

"""
    find_triangles_with_edge(edge, faces; tol=1e-4)

Find all triangles in a face list that contain a specific edge (both vertices).
Returns indices of triangles that have both edge.node1 and edge.node2 as vertices.
"""
function find_triangles_with_edge(
    edge::EdgeKey,
    faces::Vector{Triangle};
    tol::Real=1e-4
)::Vector{Int}
    triangles_with_edge = Int[]
    
    for (idx, tri) in enumerate(faces)
        has_node1 = triangle_has_node(tri, edge.node1, tol=tol)
        has_node2 = triangle_has_node(tri, edge.node2, tol=tol)
        
        if has_node1 && has_node2
            push!(triangles_with_edge, idx)
        end
    end
    
    return triangles_with_edge
end

"""
    extract_quad_vertices(tri1, tri2; tol=1e-4)

Extract the 4 unique vertices from two triangles that form a quad.
Returns a vector of 4 vertices, or empty if triangles don't form a proper quad.
"""
function extract_quad_vertices(
    tri1::Triangle,
    tri2::Triangle;
    tol::Real=1e-4
)::Vector{NTuple{3,Float64}}
    # Get all 6 vertices from both triangles
    verts1 = get_triangle_nodes(tri1)
    verts2 = get_triangle_nodes(tri2)
    
    # Helper: Check if a node is in a list using tolerance-based comparison
    function node_in_list(node::NTuple{3,Float64}, nodes::Vector{NTuple{3,Float64}})::Bool
        for n in nodes
            if are_nodes_equal(node, n, tol=tol)
                return true
            end
        end
        return false
    end
    
    # Find unique nodes
    unique_nodes = NTuple{3,Float64}[]
    for node in vcat(verts1, verts2)
        if !node_in_list(node, unique_nodes)
            push!(unique_nodes, node)
        end
    end
    
    # A valid quad should have exactly 4 unique vertices
    if length(unique_nodes) != 4
        return NTuple{3,Float64}[]
    end
    
    return unique_nodes
end

"""
    find_triangles_using_vertices(quad_vertices, faces; tol=1e-4)

Find all triangles in a face list that use ONLY vertices from the given set.
This is used to find how a quad's 4 vertices are triangulated in the target mesh.
"""
function find_triangles_using_vertices(
    quad_vertices::Vector{NTuple{3,Float64}},
    faces::Vector{Triangle};
    tol::Real=1e-4
)::Vector{Int}
    if length(quad_vertices) != 4
        return Int[]
    end
    
    # Helper: Check if a node is in the quad vertices
    function is_quad_vertex(node::NTuple{3,Float64})::Bool
        for qv in quad_vertices
            if are_nodes_equal(node, qv, tol=tol)
                return true
            end
        end
        return false
    end
    
    triangles_using_verts = Int[]
    
    for (idx, tri) in enumerate(faces)
        tri_nodes = get_triangle_nodes(tri)
        
        # Check if all 3 vertices of this triangle are in the quad vertex set
        all_in_quad = all(node -> is_quad_vertex(node), tri_nodes)
        
        if all_in_quad
            push!(triangles_using_verts, idx)
        end
    end
    
    return triangles_using_verts
end

"""
    find_quad_for_diagonal(edge, topology, present_in, target_pid; tol=1e-4)

Find the quad (4 vertices) that contains the given edge as a diagonal.
Also identifies the 2 triangles in the target side that need to be replaced.

**SOURCE-FIRST APPROACH:**
- First finds the quad in the SOURCE mesh (where the edge EXISTS)
- Extracts the 4 vertices forming this quad
- Then finds how those vertices are triangulated in the TARGET mesh
- This is much more robust than searching in the target where the edge is missing!

Returns: (quad_vertices, triangles_to_replace) or (empty, empty) if not found.
"""
function find_quad_for_diagonal(
    edge::EdgeKey,
    topology::InterfaceTopology,
    present_in::Symbol,  # :A or :B - which side HAS the edge
    target_pid::Int;     # PID where edge needs to be inserted
    tol::Real=1e-4,
    debug::Bool=false
)::Tuple{Vector{NTuple{3,Float64}}, Vector{Int}}

    # Determine source and target meshes
    # Source = where edge EXISTS (as diagonal)
    # Target = where edge is MISSING (needs to be added)
    if present_in == :A
        source_faces = topology.faces_A
        target_faces = topology.faces_B
        source_pid = topology.pidA
    else  # :B
        source_faces = topology.faces_B
        target_faces = topology.faces_A
        source_pid = topology.pidB
    end
    
    if debug
        println("\n  [DEBUG] Finding quad for diagonal edge (SOURCE-FIRST):")
        println("    Edge: $(edge.node1) → $(edge.node2)")
        println("    Present in: $present_in (PID $source_pid)")
        println("    Target PID: $target_pid")
        println("    Source faces: $(length(source_faces))")
        println("    Target faces: $(length(target_faces))")
        println("    Tolerance: $tol")
    end
    
    # STEP 1: Find the quad in SOURCE mesh where edge EXISTS
    # The diagonal edge should be shared by exactly 2 triangles in the source mesh
    source_triangles_with_edge = find_triangles_with_edge(edge, source_faces, tol=tol)
    
    if debug
        println("    Source triangles containing edge: $(length(source_triangles_with_edge))")
    end
    
    if length(source_triangles_with_edge) != 2
        if debug
            println("    [FAIL] Expected exactly 2 triangles sharing the edge in source, got $(length(source_triangles_with_edge))")
        end
        return (NTuple{3,Float64}[], Int[])
    end
    
    # STEP 2: Extract the 4 vertices forming the quad
    tri1_idx = source_triangles_with_edge[1]
    tri2_idx = source_triangles_with_edge[2]
    tri1 = source_faces[tri1_idx]
    tri2 = source_faces[tri2_idx]
    
    quad_vertices = extract_quad_vertices(tri1, tri2, tol=tol)
    
    if isempty(quad_vertices)
        if debug
            println("    [FAIL] Could not extract 4 unique vertices from source triangles")
        end
        return (NTuple{3,Float64}[], Int[])
    end
    
    if debug
        println("    [SUCCESS] Found quad in source mesh with 4 vertices")
        println("    Quad vertices:")
        for (i, v) in enumerate(quad_vertices)
            println("      V$i: $v")
        end
    end
    
    # STEP 3: Find how these SAME 4 vertices are triangulated in TARGET mesh
    target_triangles = find_triangles_using_vertices(quad_vertices, target_faces, tol=tol)
    
    if debug
        println("    Target triangles using these 4 vertices: $(length(target_triangles))")
    end
    
    if isempty(target_triangles)
        if debug
            println("    [FAIL] No triangles in target mesh use these 4 vertices")
        end
        return (NTuple{3,Float64}[], Int[])
    end
    
    # Typically should be exactly 2 triangles (different diagonal)
    if length(target_triangles) != 2
        if debug
            println("    [WARN] Expected 2 target triangles, got $(length(target_triangles)) (may be OK)")
        end
    end
    
    # STEP 4 (NEW): Verify boundary edges match between source and target
    # Extract boundary edges from source quad
    source_boundary = extract_quad_boundary_edges(tri1, tri2, edge, tol=tol)
    
    if isempty(source_boundary)
        if debug
            println("    [FAIL] Could not extract 4 boundary edges from source quad")
        end
        return (NTuple{3,Float64}[], Int[])
    end
    
    # Extract boundary edges from target quad (if we have exactly 2 triangles)
    if length(target_triangles) == 2
        target_tri1 = target_faces[target_triangles[1]]
        target_tri2 = target_faces[target_triangles[2]]
        
        # Find the diagonal in target (different from source diagonal)
        target_edges1 = [
            EdgeKey(target_tri1.coord1, target_tri1.coord2),
            EdgeKey(target_tri1.coord2, target_tri1.coord3),
            EdgeKey(target_tri1.coord3, target_tri1.coord1)
        ]
        target_edges2 = [
            EdgeKey(target_tri2.coord1, target_tri2.coord2),
            EdgeKey(target_tri2.coord2, target_tri2.coord3),
            EdgeKey(target_tri2.coord3, target_tri2.coord1)
        ]
        
        # Find common edge (the target diagonal)
        target_diagonal = nothing
        for e1 in target_edges1
            for e2 in target_edges2
                if (are_nodes_equal(e1.node1, e2.node1, tol=tol) && 
                    are_nodes_equal(e1.node2, e2.node2, tol=tol)) ||
                   (are_nodes_equal(e1.node1, e2.node2, tol=tol) && 
                    are_nodes_equal(e1.node2, e2.node1, tol=tol))
                    target_diagonal = e1
                    break
                end
            end
            if target_diagonal !== nothing
                break
            end
        end
        
        if target_diagonal === nothing
            if debug
                println("    [FAIL] Could not find diagonal in target triangulation")
            end
            return (NTuple{3,Float64}[], Int[])
        end
        
        target_boundary = extract_quad_boundary_edges(target_tri1, target_tri2, target_diagonal, tol=tol)
        
        if isempty(target_boundary)
            if debug
                println("    [FAIL] Could not extract 4 boundary edges from target quad")
            end
            return (NTuple{3,Float64}[], Int[])
        end
        
        # Check if boundary edges match
        if !edges_match(source_boundary, target_boundary, tol=tol)
            if debug
                println("    [QUAD_MISMATCH] Boundary edges don't match between source and target")
                println("    Source has $(length(source_boundary)) boundary edges")
                println("    Target has $(length(target_boundary)) boundary edges")
                println("    → This is a QUAD_MISMATCH, not a valid DIAGONAL")
            end
            # Return special marker: quad vertices but empty triangles_to_replace
            # The caller will detect this and classify as QUAD_MISMATCH
            return (quad_vertices, Int[])
        end
    end
    
    if debug
        println("    [SUCCESS] Boundary edges match - valid DIAGONAL mismatch")
        println("    Found target triangulation: triangles $target_triangles")
    end
    
    return (quad_vertices, target_triangles)
end

# ============================================================================
# Edge mismatch classification
# ============================================================================

"""
    classify_edge_mismatch(edge, topology, present_in; tol=1e-4, debug=false)

Classify a single edge mismatch and determine repair strategy.
"""
function classify_edge_mismatch(edge::EdgeKey, 
                               topology::InterfaceTopology,
                               present_in::Symbol;  # :A or :B
                               tol::Real=1e-4,
                               debug::Bool=false)::EdgeMismatch
    
    # Determine which side has the edge and which needs it
    if present_in == :A
        source_faces = topology.faces_A
        target_faces = topology.faces_B
        target_pid = topology.pidB
        should_be_in = topology.pidB
    else  # :B
        source_faces = topology.faces_B
        target_faces = topology.faces_A
        target_pid = topology.pidA
        should_be_in = topology.pidA
    end
    
    # EARLY CHECK 0: Check if edge is degenerate (zero or near-zero length)
    dx = edge.node2[1] - edge.node1[1]
    dy = edge.node2[2] - edge.node1[2]
    dz = edge.node2[3] - edge.node1[3]
    edge_length = sqrt(dx*dx + dy*dy + dz*dz)
    
    if edge_length < tol
        if debug
            println("  [DEGENERATE_EDGE DETECTED]")
            println("    Edge length: $edge_length (< $tol)")
            println("    → Edge has zero or near-zero length")
        end
        
        return EdgeMismatch(
            edge, DEGENERATE_EDGE, 
            present_in == :A ? :B_only : :A_only, should_be_in,
            NTuple{3,Float64}[], Int[], NTuple{3,Float64}[], Int[],
            0.95, 1.0, false,
            UnknownDiagnostics(present_in, 0, 0, (false, false), edge_length,
                              false, false, 0, nothing, "Edge has degenerate length < $tol")
        )
    end
    
    # Extract all vertices from the TARGET side to search for hanging nodes
    # Hanging nodes are vertices in the target mesh that lie ON the source edge
    # Include both target face vertices and shared node keys to catch topology nodes
    # that may not be part of the target interface faces in this minimal context
    target_nodes = Set{NTuple{3,Float64}}()
    for tri in target_faces
        push!(target_nodes, tri.coord1, tri.coord2, tri.coord3)
    end
    # Also consider shared node keys (already rounded) as candidate nodes
    # Tests may inject hanging nodes via shared_node_keys only
    for k in topology.shared_node_keys
        push!(target_nodes, k)
    end
    
    if debug
        println("  [DEBUG] Searching for hanging nodes:")
        println("    Edge: $(edge.node1) -> $(edge.node2)")
        println("    Present in: $present_in")
        println("    Edge length: $edge_length")
        println("    Target side has $(length(target_nodes)) vertices")
        println("    Target side vertices: $target_nodes")
    end
    
    # Find nodes from target side that lie on this edge
    hanging = find_hanging_nodes_on_edge(edge, target_nodes, tol=tol)
    
    if debug
        println("    Found $(length(hanging)) hanging nodes: $hanging")
    end
    
    # Classify based on hanging nodes
    mismatch_type = UNKNOWN
    complexity = 0.5
    quad_vertices = NTuple{3,Float64}[]
    triangles_to_replace = Int[]
    
    if length(hanging) == 1
        # Single hanging node → T-junction (most common)
        mismatch_type = T_JUNCTION
        complexity = 0.2
    elseif length(hanging) > 1
        # Multiple hanging nodes → Refinement level difference
        mismatch_type = REFINEMENT
        complexity = 0.6 + 0.1 * length(hanging)
    elseif length(hanging) == 0
        # No hanging nodes → Check various edge topology scenarios
        
        # EARLY CHECK 1: How many triangles share this edge in source?
        source_triangles_with_edge = find_triangles_with_edge(edge, source_faces, tol=tol)
        
        if length(source_triangles_with_edge) == 1
            # Edge appears in only 1 triangle → BOUNDARY_EDGE
            mismatch_type = BOUNDARY_EDGE
            complexity = 0.85
            
            if debug
                println("  [BOUNDARY_EDGE DETECTED]")
                println("    Edge appears in only 1 triangle in source mesh")
                println("    → This is a boundary edge, not an internal interface edge")
            end
        elseif length(source_triangles_with_edge) > 2
            # Edge shared by >2 triangles → NON_MANIFOLD
            mismatch_type = NON_MANIFOLD
            complexity = 0.95
            
            if debug
                println("  [NON_MANIFOLD DETECTED]")
                println("    Edge shared by $(length(source_triangles_with_edge)) triangles (expected 2)")
                println("    → Topology is non-manifold, cannot repair")
            end
        elseif length(source_triangles_with_edge) == 0
            # Edge doesn't exist in any triangle - data anomaly
            mismatch_type = SOURCE_EDGE_ABSENT
            complexity = 0.95
            
            if debug
                println("  [SOURCE_EDGE_ABSENT DETECTED]")
                println("    Edge doesn't appear in any source triangles")
                println("    → Data anomaly: edge key exists but not used by source geometry")
            end
        else
            # length == 2: Normal case, proceed with diagonal detection
            # Extract vertices from both sides
            vertices_A = extract_boundary_vertices(topology.faces_A)
            vertices_B = extract_boundary_vertices(topology.faces_B)
            shared_vertices = compute_shared_vertices(vertices_A, vertices_B, tol=tol)
            
            # Check if both endpoints are in the shared vertex set
            endpoint1_shared = is_vertex_in_set(edge.node1, shared_vertices, tol=tol)
            endpoint2_shared = is_vertex_in_set(edge.node2, shared_vertices, tol=tol)
            
            # EARLY CHECK 2: Are both endpoints shared?
            if !endpoint1_shared || !endpoint2_shared
                # One or both endpoints not shared → UNSHARED_ENDPOINT
                mismatch_type = UNSHARED_ENDPOINT
                complexity = 0.95
                
                if debug
                    println("  [UNSHARED_ENDPOINT DETECTED]")
                    println("    Endpoint 1 shared: $endpoint1_shared")
                    println("    Endpoint 2 shared: $endpoint2_shared")
                    println("    → Cannot be repaired - vertices not properly shared")
                end
            # Proceed only if both endpoints are shared
            elseif endpoint1_shared && endpoint2_shared
                # Both endpoints are shared → POTENTIAL diagonal candidate
                # Need to verify by finding the quad and checking boundary edges
                quad_vertices, triangles_to_replace = find_quad_for_diagonal(edge, topology, present_in, should_be_in, tol=tol, debug=debug)
                
                if !isempty(quad_vertices) && !isempty(triangles_to_replace)
                    # Successfully found quad with matching boundary edges → legitimate DIAGONAL
                    mismatch_type = DIAGONAL
                    complexity = 0.4
                elseif !isempty(quad_vertices) && isempty(triangles_to_replace)
                    # Found quad vertices but boundary edges don't match → QUAD_MISMATCH
                    mismatch_type = QUAD_MISMATCH
                    complexity = 0.8
                    
                    if debug
                        println("  [QUAD_MISMATCH DETECTED]")
                        println("    Same 4 vertices but different quad boundary topology")
                        println("    → Cannot repair with simple diagonal retriangulation")
                    end
                else
                    # Quad finding failed → Need to determine why
                    # We know: both endpoints shared, source has 2 triangles with the edge
                    
                    # Try to extract quad vertices from source triangles
                    source_tris_with_edge = find_triangles_with_edge(edge, source_faces, tol=tol)
                    if length(source_tris_with_edge) == 2
                        tri1 = source_faces[source_tris_with_edge[1]]
                        tri2 = source_faces[source_tris_with_edge[2]]
                        test_quad_verts = extract_quad_vertices(tri1, tri2, tol=tol)
                        
                        if isempty(test_quad_verts)
                            # Cannot extract 4 unique vertices from source
                            mismatch_type = QUAD_NOT_FOUND_IN_SOURCE
                            complexity = 0.75
                            
                            if debug
                                println("  [QUAD_NOT_FOUND_IN_SOURCE DETECTED]")
                                println("    Source has 2 triangles sharing the edge but cannot extract 4 unique vertices")
                                println("    → Likely degenerate or collapsed quad in source")
                            end
                        else
                            # Quad vertices found in source, check if target uses them
                            target_tris_using_quad = find_triangles_using_vertices(test_quad_verts, target_faces, tol=tol)
                            
                            if isempty(target_tris_using_quad)
                                # Target doesn't triangulate using only these 4 vertices
                                # Target likely uses additional vertices (finer triangulation)
                                mismatch_type = TARGET_USES_FINER_TRIANGULATION
                                complexity = 0.8
                                
                                if debug
                                    println("  [TARGET_USES_FINER_TRIANGULATION DETECTED]")
                                    println("    Found 4 vertices in source quad: $test_quad_verts")
                                    println("    But target has no triangles using ONLY these 4 vertices")
                                    println("    → Target boundary mesh uses additional vertices beyond the quad set")
                                    println("    → This indicates target has finer/different triangulation with extra vertices")
                                end
                            else
                                # Target has triangles using these vertices but quad finding still failed
                                # This could be due to boundary mismatch, non-planar quad, etc.
                                mismatch_type = UNKNOWN
                                complexity = 0.7
                                
                                if debug
                                    println("  [UNKNOWN - COMPLEX QUAD FAILURE]")
                                    println("    Source quad vertices found: $test_quad_verts")
                                    println("    Target has $(length(target_tris_using_quad)) triangles using these vertices")
                                    println("    But diagonal check failed (likely boundary mismatch or non-2-triangle configuration)")
                                end
                            end
                        end
                    else
                        # Should not reach here as we already checked length == 2
                        mismatch_type = UNKNOWN
                        complexity = 0.7
                        
                        if debug
                            println("  [UNKNOWN - UNEXPECTED STATE]")
                            println("    Source triangle count unexpected: $(length(source_tris_with_edge))")
                        end
                    end
                end
            end
        end
    end
    
    # Find affected triangles in target side
    affected_triangle_indices = Int[]
    min_quality = 1.0
    
    for (idx, tri) in enumerate(target_faces)
        # Check if this triangle's vertices include both endpoints of the edge
        has_node1 = triangle_has_node(tri, edge.node1, tol=tol)
        has_node2 = triangle_has_node(tri, edge.node2, tol=tol)
        
        if has_node1 && has_node2
            push!(affected_triangle_indices, idx)
            quality = compute_triangle_quality(tri)
            min_quality = min(min_quality, quality)
        end
    end
    
    # Assess repair feasibility
    # Infeasible if: no affected triangles, or quality already very poor
    repair_feasible = !isempty(affected_triangle_indices) && min_quality > 0.05
    
    # Adjust complexity based on quality
    if min_quality < 0.2
        complexity += 0.3  # Poor quality makes repair harder
    end
    
    complexity = clamp(complexity, 0.0, 1.0)
    
    # Create diagnostics for UNKNOWN and refined unknown types
    diagnostics = nothing
    if mismatch_type in [UNKNOWN, SOURCE_EDGE_ABSENT, QUAD_NOT_FOUND_IN_SOURCE, TARGET_USES_FINER_TRIANGULATION, DEGENERATE_EDGE]
        # Gather diagnostic information
        source_triangle_count = length(find_triangles_with_edge(edge, source_faces, tol=tol))
        target_triangle_count = length(affected_triangle_indices)
        
        # Extract vertices from both sides for endpoint sharing check
        vertices_A = extract_boundary_vertices(topology.faces_A)
        vertices_B = extract_boundary_vertices(topology.faces_B)
        shared_vertices = compute_shared_vertices(vertices_A, vertices_B, tol=tol)
        endpoint1_shared = is_vertex_in_set(edge.node1, shared_vertices, tol=tol)
        endpoint2_shared = is_vertex_in_set(edge.node2, shared_vertices, tol=tol)
        
        tried_quad = source_triangle_count == 2 && endpoint1_shared && endpoint2_shared
        quad_found = !isempty(quad_vertices)
        target_using_quad = if quad_found
            length(find_triangles_using_vertices(quad_vertices, target_faces, tol=tol))
        else
            0
        end
        
        boundary_match = if quad_found && !isempty(triangles_to_replace)
            true
        elseif quad_found && isempty(triangles_to_replace)
            false
        else
            nothing
        end
        
        reason = if mismatch_type == DEGENERATE_EDGE
            "Edge has degenerate length < $tol"
        elseif mismatch_type == SOURCE_EDGE_ABSENT
            "Edge not found in any source triangle"
        elseif mismatch_type == QUAD_NOT_FOUND_IN_SOURCE
            "Cannot extract 4 unique vertices from 2 source triangles"
        elseif mismatch_type == TARGET_USES_FINER_TRIANGULATION
            "Target has no triangles using the source quad's 4 vertices"
        else
            "Complex failure: endpoints shared, source has $source_triangle_count triangles, quad_found=$quad_found, target_using_quad=$target_using_quad"
        end
        
        diagnostics = UnknownDiagnostics(
            present_in, source_triangle_count, target_triangle_count,
            (endpoint1_shared, endpoint2_shared), edge_length,
            tried_quad, quad_found, target_using_quad, boundary_match,
            reason
        )
    end
    
    return EdgeMismatch(
        edge,
        mismatch_type,
        present_in == :A ? :B_only : :A_only,
        should_be_in,
        hanging,
        affected_triangle_indices,
        quad_vertices,
        triangles_to_replace,
        complexity,
        min_quality,
        repair_feasible,
        diagnostics
    )
end

"""
    classify_interface_mismatches(topology; tol=1e-4, debug=false, debug_samples=3)

Classify all edge mismatches for an interface.
Returns complete classification report.

If debug=true, will print detailed diagnostics for the first debug_samples mismatches.
"""
function classify_interface_mismatches(topology::InterfaceTopology;
                                      tol::Real=1e-4,
                                      debug::Bool=false,
                                      debug_samples::Int=3,
                                      verbose::Bool=true)::InterfaceClassification
    
    if verbose
        println("Classifying edge mismatches for PID=$(topology.pidA) ↔ PID=$(topology.pidB)...")
        if debug
            println("  [DEBUG MODE ENABLED] Will show detailed diagnostics for first $debug_samples mismatches")
        end
    end
    
    # Classify edges missing in B (present only in A)
    mismatches_B = EdgeMismatch[]
    sample_count_B = 0
    for edge in topology.edges_only_in_A
        enable_debug = debug && sample_count_B < debug_samples
        if enable_debug
            sample_count_B += 1
            println("\n  === Analyzing mismatch $(sample_count_B) (edge missing in B) ===")
        end
        mismatch = classify_edge_mismatch(edge, topology, :A, tol=tol, debug=enable_debug)
        push!(mismatches_B, mismatch)
    end
    
    # Classify edges missing in A (present only in B)
    mismatches_A = EdgeMismatch[]
    sample_count_A = 0
    for edge in topology.edges_only_in_B
        enable_debug = debug && sample_count_A < debug_samples
        if enable_debug
            sample_count_A += 1
            println("\n  === Analyzing mismatch $(sample_count_A) (edge missing in A) ===")
        end
        mismatch = classify_edge_mismatch(edge, topology, :B, tol=tol, debug=enable_debug)
        push!(mismatches_A, mismatch)
    end
    
    # Count by type
    all_mismatches = vcat(mismatches_A, mismatches_B)
    t_junction_count = count(m -> m.mismatch_type == T_JUNCTION, all_mismatches)
    diagonal_count = count(m -> m.mismatch_type == DIAGONAL, all_mismatches)
    refinement_count = count(m -> m.mismatch_type == REFINEMENT, all_mismatches)
    quad_mismatch_count = count(m -> m.mismatch_type == QUAD_MISMATCH, all_mismatches)
    boundary_edge_count = count(m -> m.mismatch_type == BOUNDARY_EDGE, all_mismatches)
    non_manifold_count = count(m -> m.mismatch_type == NON_MANIFOLD, all_mismatches)
    unshared_endpoint_count = count(m -> m.mismatch_type == UNSHARED_ENDPOINT, all_mismatches)
    degenerate_edge_count = count(m -> m.mismatch_type == DEGENERATE_EDGE, all_mismatches)
    source_edge_absent_count = count(m -> m.mismatch_type == SOURCE_EDGE_ABSENT, all_mismatches)
    quad_not_found_in_source_count = count(m -> m.mismatch_type == QUAD_NOT_FOUND_IN_SOURCE, all_mismatches)
    target_uses_finer_triangulation_count = count(m -> m.mismatch_type == TARGET_USES_FINER_TRIANGULATION, all_mismatches)
    unknown_count = count(m -> m.mismatch_type == UNKNOWN, all_mismatches)
    
    # Repair assessment
    feasible_count = count(m -> m.repair_feasible, all_mismatches)
    infeasible_count = length(all_mismatches) - feasible_count
    
    avg_complexity = if !isempty(all_mismatches)
        sum(m.complexity_score for m in all_mismatches) / length(all_mismatches)
    else
        0.0
    end
    
    if verbose
        println("  Classification complete:")
        println("    T-junctions: $t_junction_count")
        println("    Diagonal mismatches: $diagonal_count")
        println("    Refinement differences: $refinement_count")
        println("    Quad mismatches: $quad_mismatch_count")
        println("    Boundary edges: $boundary_edge_count")
        println("    Non-manifold: $non_manifold_count")
        println("    Unshared endpoints: $unshared_endpoint_count")
        println("    Degenerate edges: $degenerate_edge_count")
        println("    Source edge absent: $source_edge_absent_count")
        println("    Quad not found in source: $quad_not_found_in_source_count")
        println("    Target uses finer triangulation: $target_uses_finer_triangulation_count")
        println("    Unknown: $unknown_count")
        println("    Feasible for repair: $feasible_count / $(length(all_mismatches))")
        println("    Average complexity: $(round(avg_complexity, digits=2))")
    end
    
    return InterfaceClassification(
        topology,
        mismatches_A,
        mismatches_B,
        t_junction_count,
        diagonal_count,
        refinement_count,
        quad_mismatch_count,
        boundary_edge_count,
        non_manifold_count,
        unshared_endpoint_count,
        degenerate_edge_count,
        source_edge_absent_count,
        quad_not_found_in_source_count,
        target_uses_finer_triangulation_count,
        unknown_count,
        feasible_count,
        infeasible_count,
        avg_complexity
    )
end

# ============================================================================
# Bidirectional edge classification (strategic improvement)
# ============================================================================

"""
    classify_interface_mismatches_bidirectional(topology; tol=1e-4, verbose=true)

Classify interface mismatches from BOTH directional perspectives and merge results.

The classification process is inherently directional:
- T-junctions are detected by finding hanging nodes on the "source" edge in the "target" mesh
- Diagonal mismatches depend on which triangulation is considered "source" vs "target"
- Quality metrics are computed from the perspective of which mesh will be modified

By classifying from both perspectives, we:
1. Discover mismatches that might only be visible from one direction
2. Get more complete quality and feasibility assessments
3. Avoid missing repair opportunities due to perspective bias

Returns:
- classification_merged: Combined classification with all unique mismatches
- classification_A_perspective: Classification treating A as source
- classification_B_perspective: Classification treating B as source
- comparison_metrics: Dictionary with comparison statistics
"""
function classify_interface_mismatches_bidirectional(
    topology::InterfaceTopology;
    tol::Real=1e-4,
    verbose::Bool=true
)::Tuple{InterfaceClassification, InterfaceClassification, InterfaceClassification, Dict{String,Any}}
    
    if verbose
        println("\n" * "="^70)
        println("BIDIRECTIONAL MISMATCH CLASSIFICATION")
        println("="^70)
        println("Strategy: Classify from both A→B and B→A perspectives")
        println("          to discover all possible mismatches")
        println()
    end
    
    # -------------------------------------------------------------------------
    # Perspective 1: A as source, B as target (original direction)
    # -------------------------------------------------------------------------
    if verbose
        println("-"^70)
        println("PERSPECTIVE 1: A (PID=$(topology.pidA)) as source, B (PID=$(topology.pidB)) as target")
        println("-"^70)
    end
    
    classification_AB = classify_interface_mismatches(topology, tol=tol, debug=false)
    
    if verbose
        println("  Result:")
        println("    Edges missing in A: $(length(classification_AB.mismatches_A))")
        println("    Edges missing in B: $(length(classification_AB.mismatches_B))")
        println("    Total mismatches: $(length(classification_AB.mismatches_A) + length(classification_AB.mismatches_B))")
    end
    
    # -------------------------------------------------------------------------
    # Perspective 2: B as source, A as target (reversed direction)
    # -------------------------------------------------------------------------
    if verbose
        println("\n" * "-"^70)
        println("PERSPECTIVE 2: B (PID=$(topology.pidB)) as source, A (PID=$(topology.pidA)) as target")
        println("-"^70)
    end
    
    # Create a swapped topology to reverse the perspective
    topology_swapped = InterfaceTopology(
        topology.pidB,                 # Swap: B becomes pidA
        topology.pidA,                 # Swap: A becomes pidB
        topology.shared_node_keys,     # Shared nodes remain the same
        topology.node_key_to_ids,      # Node mappings remain the same
        topology.faces_B,              # Swap: B's faces become faces_A
        topology.faces_A,              # Swap: A's faces become faces_B
        topology.edges_B,              # Swap edge dictionaries
        topology.edges_A,
        topology.edges_only_in_B,      # Swap: edges only in B become "only in A"
        topology.edges_only_in_A,      # Swap: edges only in A become "only in B"
        topology.edges_shared,         # Shared edges remain the same
        topology.interface_bbox,       # Interface bbox remains the same
        topology.total_shared_nodes,   # Stats
        topology.total_faces_B,        # Swap face counts
        topology.total_faces_A,
        topology.total_edges_B,        # Swap edge counts
        topology.total_edges_A,
        topology.conformity_ratio      # Conformity ratio remains the same
    )
    
    classification_BA = classify_interface_mismatches(topology_swapped, tol=tol, debug=false)
    
    if verbose
        println("  Result:")
        println("    Edges missing in B (from this perspective): $(length(classification_BA.mismatches_A))")
        println("    Edges missing in A (from this perspective): $(length(classification_BA.mismatches_B))")
        println("    Total mismatches: $(length(classification_BA.mismatches_A) + length(classification_BA.mismatches_B))")
    end
    
    # -------------------------------------------------------------------------
    # Merge results: Combine unique mismatches from both perspectives
    # -------------------------------------------------------------------------
    if verbose
        println("\n" * "="^70)
        println("MERGING RESULTS")
        println("="^70)
    end
    
    # For the merged result, we need to unswap the classifications from perspective 2
    # classification_BA has:
    #   - mismatches_A: edges missing in the "new A" (which is actually original B)
    #   - mismatches_B: edges missing in the "new B" (which is actually original A)
    # So we need to swap them back:
    mismatches_A_from_BA = classification_BA.mismatches_B  # These are for original A
    mismatches_B_from_BA = classification_BA.mismatches_A  # These are for original B
    
    # Merge mismatch lists, removing duplicates based on edge coordinates
    function merge_mismatches(list1::Vector{EdgeMismatch}, list2::Vector{EdgeMismatch})::Vector{EdgeMismatch}
        merged = copy(list1)
        
        for m2 in list2
            # Check if this edge already exists in merged list
            edge_key2 = m2.edge_key
            found = false
            
            for m1 in merged
                edge_key1 = m1.edge_key
                # Check if edges are the same (either direction)
                if ((edge_key1.node1 == edge_key2.node1 && edge_key1.node2 == edge_key2.node2) ||
                    (edge_key1.node1 == edge_key2.node2 && edge_key1.node2 == edge_key2.node1))
                    found = true
                    break
                end
            end
            
            # If not found, add it
            if !found
                push!(merged, m2)
            end
        end
        
        return merged
    end
    
    mismatches_A_merged = merge_mismatches(classification_AB.mismatches_A, mismatches_A_from_BA)
    mismatches_B_merged = merge_mismatches(classification_AB.mismatches_B, mismatches_B_from_BA)
    
    # Recompute statistics for merged classification
    all_mismatches = vcat(mismatches_A_merged, mismatches_B_merged)
    t_junction_count = count(m -> m.mismatch_type == T_JUNCTION, all_mismatches)
    diagonal_count = count(m -> m.mismatch_type == DIAGONAL, all_mismatches)
    refinement_count = count(m -> m.mismatch_type == REFINEMENT, all_mismatches)
    quad_mismatch_count = count(m -> m.mismatch_type == QUAD_MISMATCH, all_mismatches)
    boundary_edge_count = count(m -> m.mismatch_type == BOUNDARY_EDGE, all_mismatches)
    non_manifold_count = count(m -> m.mismatch_type == NON_MANIFOLD, all_mismatches)
    unshared_endpoint_count = count(m -> m.mismatch_type == UNSHARED_ENDPOINT, all_mismatches)
    degenerate_edge_count = count(m -> m.mismatch_type == DEGENERATE_EDGE, all_mismatches)
    source_edge_absent_count = count(m -> m.mismatch_type == SOURCE_EDGE_ABSENT, all_mismatches)
    quad_not_found_in_source_count = count(m -> m.mismatch_type == QUAD_NOT_FOUND_IN_SOURCE, all_mismatches)
    target_uses_finer_triangulation_count = count(m -> m.mismatch_type == TARGET_USES_FINER_TRIANGULATION, all_mismatches)
    unknown_count = count(m -> m.mismatch_type == UNKNOWN, all_mismatches)
    feasible_count = count(m -> m.repair_feasible, all_mismatches)
    infeasible_count = length(all_mismatches) - feasible_count
    avg_complexity = isempty(all_mismatches) ? 0.0 : sum(m.complexity_score for m in all_mismatches) / length(all_mismatches)
    
    classification_merged = InterfaceClassification(
        topology,
        mismatches_A_merged,
        mismatches_B_merged,
        t_junction_count,
        diagonal_count,
        refinement_count,
        quad_mismatch_count,
        boundary_edge_count,
        non_manifold_count,
        unshared_endpoint_count,
        degenerate_edge_count,
        source_edge_absent_count,
        quad_not_found_in_source_count,
        target_uses_finer_triangulation_count,
        unknown_count,
        feasible_count,
        infeasible_count,
        avg_complexity
    )
    
    # Compute comparison metrics
    orig_A_count = length(classification_AB.mismatches_A)
    orig_B_count = length(classification_AB.mismatches_B)
    new_A_count = length(mismatches_A_from_BA)
    new_B_count = length(mismatches_B_from_BA)
    merged_A_count = length(mismatches_A_merged)
    merged_B_count = length(mismatches_B_merged)
    
    added_A = merged_A_count - orig_A_count
    added_B = merged_B_count - orig_B_count
    
    comparison_metrics = Dict(
        "perspective_AB" => Dict(
            "mismatches_A" => orig_A_count,
            "mismatches_B" => orig_B_count,
            "total" => orig_A_count + orig_B_count
        ),
        "perspective_BA" => Dict(
            "mismatches_A" => new_A_count,
            "mismatches_B" => new_B_count,
            "total" => new_A_count + new_B_count
        ),
        "merged" => Dict(
            "mismatches_A" => merged_A_count,
            "mismatches_B" => merged_B_count,
            "total" => merged_A_count + merged_B_count
        ),
        "discovered" => Dict(
            "additional_in_A" => added_A,
            "additional_in_B" => added_B,
            "total_discovered" => added_A + added_B
        )
    )
    
    if verbose
        println("\nComparison:")
        println("  Perspective A→B:")
        println("    Edges missing in A: $orig_A_count")
        println("    Edges missing in B: $orig_B_count")
        println("  Perspective B→A:")
        println("    Edges missing in A: $new_A_count")
        println("    Edges missing in B: $new_B_count")
        println("\n  Merged result:")
        println("    Edges missing in A: $merged_A_count (+$added_A discovered)")
        println("    Edges missing in B: $merged_B_count (+$added_B discovered)")
        println("    Total mismatches: $(merged_A_count + merged_B_count)")
        
        if added_A + added_B > 0
            println("\n  ✓ Bidirectional analysis discovered $(added_A + added_B) additional mismatches!")
        else
            println("\n  → No additional mismatches found (both perspectives agree)")
        end
        
        println("="^70)
    end
    
    return (classification_merged, classification_AB, classification_BA, comparison_metrics)
end

"""
    export_classification_json(classification, output_file)

Export classification report to JSON.
"""
function export_classification_json(classification::InterfaceClassification, output_file::String)
    
    function mismatch_to_dict(m::EdgeMismatch)
        d = Dict(
            "edge" => [[m.edge_key.node1...], [m.edge_key.node2...]],
            "type" => string(m.mismatch_type),
            "present_in" => string(m.present_in),
            "should_be_in_pid" => m.should_be_in,
            "hanging_nodes" => [[n...] for n in m.hanging_nodes],
            "affected_triangles" => m.affected_triangles,
            "quad_vertices" => [[v...] for v in m.quad_vertices],
            "triangles_to_replace" => m.triangles_to_replace,
            "complexity_score" => round(m.complexity_score, digits=3),
            "min_quality" => round(m.min_affected_triangle_quality, digits=3),
            "repair_feasible" => m.repair_feasible
        )
        
        # Add diagnostics if present
        if m.diagnostics !== nothing
            d["diagnostics"] = Dict(
                "reason" => m.diagnostics.reason,
                "source_triangle_count" => m.diagnostics.source_triangle_count,
                "target_triangle_count" => m.diagnostics.target_triangle_count_using_endpoints,
                "endpoints_shared" => m.diagnostics.endpoints_shared,
                "edge_length" => round(m.diagnostics.edge_length, digits=6),
                "tried_quad_finding" => m.diagnostics.tried_quad_finding,
                "quad_vertices_found" => m.diagnostics.quad_vertices_found,
                "target_triangles_using_quad" => m.diagnostics.target_triangles_using_quad
            )
        end
        
        return d
    end
    
    report = Dict(
        "interface" => Dict(
            "pidA" => classification.topology.pidA,
            "pidB" => classification.topology.pidB
        ),
        "summary" => Dict(
            "total_mismatches" => length(classification.mismatches_A) + length(classification.mismatches_B),
            "mismatches_in_A" => length(classification.mismatches_A),
            "mismatches_in_B" => length(classification.mismatches_B),
            "by_type" => Dict(
                "t_junction" => classification.t_junction_count,
                "diagonal" => classification.diagonal_count,
                "refinement" => classification.refinement_count,
                "quad_mismatch" => classification.quad_mismatch_count,
                "boundary_edge" => classification.boundary_edge_count,
                "non_manifold" => classification.non_manifold_count,
                "unshared_endpoint" => classification.unshared_endpoint_count,
                "degenerate_edge" => classification.degenerate_edge_count,
                "source_edge_absent" => classification.source_edge_absent_count,
                "quad_not_found_in_source" => classification.quad_not_found_in_source_count,
                "target_uses_finer_triangulation" => classification.target_uses_finer_triangulation_count,
                "unknown" => classification.unknown_count
            ),
            "repair_assessment" => Dict(
                "feasible" => classification.total_feasible,
                "infeasible" => classification.total_infeasible,
                "average_complexity" => round(classification.average_complexity, digits=3)
            )
        ),
        "mismatches_missing_in_A" => [mismatch_to_dict(m) for m in classification.mismatches_A[1:min(100, length(classification.mismatches_A))]],
        "mismatches_missing_in_B" => [mismatch_to_dict(m) for m in classification.mismatches_B[1:min(100, length(classification.mismatches_B))]]
    )
    
    open(output_file, "w") do io
        # Reuse write_json from interface_topology
        write_json(io, report, 0)
    end
    
    println("Classification report exported to: $output_file")
    return output_file
end
