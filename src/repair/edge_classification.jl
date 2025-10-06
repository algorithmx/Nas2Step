# Edge mismatch classification for interface repair
# Phase 1.2: Classify edge mismatches into repair types

# ============================================================================
# Edge mismatch types and classification
# ============================================================================

@enum MismatchType begin
    T_JUNCTION      # One region has intermediate node, other doesn't
    DIAGONAL        # Same quad, different triangulation
    REFINEMENT      # Hierarchical mesh refinement difference
    UNKNOWN         # Cannot classify
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
    unknown_count::Int
    
    # Repair assessment
    total_feasible::Int
    total_infeasible::Int
    average_complexity::Float64
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
        if node == edge.node1 || node == edge.node2
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
# Quad finding for diagonal mismatches
# ============================================================================

"""
    find_quad_for_diagonal(edge, topology, target_pid)

Find the quad (4 vertices) that contains the given edge as a diagonal.
Also identifies the 2 triangles in the target side that need to be replaced.

Returns: (quad_vertices, triangles_to_replace) or (empty, empty) if not found.
"""
function find_quad_for_diagonal(
    edge::EdgeKey,
    topology::InterfaceTopology,
    target_pid::Int;
    rd=x->round(x, digits=8)
)::Tuple{Vector{NTuple{3,Float64}}, Vector{Int}}

    coord2tuple(c) = (rd(c[1]), rd(c[2]), rd(c[3]))
    triangle_to_set(tri) = Set([coord2tuple(tri.coord1), coord2tuple(tri.coord2), coord2tuple(tri.coord3)])

    # Get target side faces
    target_faces = target_pid == topology.pidA ? topology.faces_A : topology.faces_B
    
    # The edge has 2 endpoints - these are opposite corners of the quad
    corner1 = coord2tuple(edge.node1)
    corner2 = coord2tuple(edge.node2)
    
    # Find triangles that have corner1 as a vertex
    triangles_with_corner1 = Int[
        idx for (idx, tri) in enumerate(target_faces)
            if corner1 in triangle_to_set(tri)
    ]

    # Find triangles that have corner2 as a vertex
    triangles_with_corner2 = Int[
        idx for (idx, tri) in enumerate(target_faces)
            if corner2 in triangle_to_set(tri)
    ]
    
    if isempty(triangles_with_corner1) || isempty(triangles_with_corner2)
        # No triangles found with one of the corners; cannot form a quad
        @warn("No triangles found with corners of edge; cannot form quad.")
        return (NTuple{3,Float64}[], Int[])
    end

    # For a quad with diagonal from corner1 to corner2:
    # - There should be 2 triangles sharing an edge (the OTHER diagonal)
    # - These 2 triangles together form the quad
    # - Find triangles that share an edge and connect corner1 and corner2
    
    for idx1 in triangles_with_corner1
        for idx2 in triangles_with_corner2
            if idx1 == idx2
                @warn("Triangles $idx1 and $idx2 are the same; skipping.")
                continue  # Same triangle can't form a quad
            end
            
            tri1 = target_faces[idx1]
            tri2 = target_faces[idx2]
            
            # Get vertices of both triangles
            verts1 = triangle_to_set(tri1)
            verts2 = triangle_to_set(tri2)
            
            # Triangles should share exactly two vertices (a common edge)
            shared = intersect(verts1, verts2)
            if length(shared) != 2
                @warn("Triangles $idx1 and $idx2 do not share an edge.")
                continue
            end

            # The union of vertices should yield exactly 4 unique quad vertices
            all_verts = union(verts1, verts2)
            if length(all_verts) != 4
                @warn("Triangles $idx1 and $idx2 do not form a quad (found $(length(all_verts)) unique vertices).")
                continue
            end

            # The two corners (edge endpoints) must both be vertices of this quad
            if !(corner1 in all_verts && corner2 in all_verts)
                @warn("Triangles $idx1 and $idx2 do not contain both edge endpoints.")
                continue
            end

            # Ensure the shared edge is NOT the mismatch diagonal itself.
            # We want the triangles that share the opposite diagonal (the other edge).
            if Set(shared) == Set([corner1, corner2])
                # These triangles would already share the diagonal edge; skip.
                @warn("Triangles $idx1 and $idx2 share the mismatch diagonal; skipping.")
                continue
            end

            # Found a valid quad and the two target triangles forming the other diagonal
            quad_verts = collect(all_verts)
            return (quad_verts, [idx1, idx2])
        end
    end
    
    # Couldn't find a valid quad
    return (NTuple{3,Float64}[], Int[])
end

# ============================================================================
# Edge mismatch classification
# ============================================================================

"""
    classify_edge_mismatch(edge, topology, present_in; tol=1e-4)

Classify a single edge mismatch and determine repair strategy.
"""
function classify_edge_mismatch(edge::EdgeKey, 
                               topology::InterfaceTopology,
                               present_in::Symbol;  # :A or :B
                               tol::Real=1e-4)::EdgeMismatch
    
    # Determine which side has the edge and which needs it
    if present_in == :A
        source_faces = topology.faces_A
        target_faces = topology.faces_B
        source_nodes = topology.shared_node_keys
        target_pid = topology.pidB
        should_be_in = topology.pidB
    else  # :B
        source_faces = topology.faces_B
        target_faces = topology.faces_A
        source_nodes = topology.shared_node_keys
        target_pid = topology.pidA
        should_be_in = topology.pidA
    end
    
    # Find nodes from target side that lie on this edge
    hanging = find_hanging_nodes_on_edge(edge, source_nodes, tol=tol)
    
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
        # No hanging nodes → Likely diagonal mismatch
        mismatch_type = DIAGONAL
        complexity = 0.4
        
        # Find the quad for diagonal mismatches
        quad_vertices, triangles_to_replace = find_quad_for_diagonal(edge, topology, should_be_in)
    end
    
    # Find affected triangles in target side
    affected_triangle_indices = Int[]
    min_quality = 1.0
    
    for (idx, tri) in enumerate(target_faces)
        # Check if this triangle's vertices include both endpoints of the edge
        tri_keys = Set([
            (round(tri.coord1[1], digits=8), round(tri.coord1[2], digits=8), round(tri.coord1[3], digits=8)),
            (round(tri.coord2[1], digits=8), round(tri.coord2[2], digits=8), round(tri.coord2[3], digits=8)),
            (round(tri.coord3[1], digits=8), round(tri.coord3[2], digits=8), round(tri.coord3[3], digits=8))
        ])
        
        if edge.node1 in tri_keys && edge.node2 in tri_keys
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
        repair_feasible
    )
end

"""
    classify_interface_mismatches(topology; tol=1e-4)

Classify all edge mismatches for an interface.
Returns complete classification report.
"""
function classify_interface_mismatches(topology::InterfaceTopology;
                                      tol::Real=1e-4)::InterfaceClassification
    
    println("Classifying edge mismatches for PID=$(topology.pidA) ↔ PID=$(topology.pidB)...")
    
    # Classify edges missing in B (present only in A)
    mismatches_B = EdgeMismatch[]
    for edge in topology.edges_only_in_A
        mismatch = classify_edge_mismatch(edge, topology, :A, tol=tol)
        push!(mismatches_B, mismatch)
    end
    
    # Classify edges missing in A (present only in B)
    mismatches_A = EdgeMismatch[]
    for edge in topology.edges_only_in_B
        mismatch = classify_edge_mismatch(edge, topology, :B, tol=tol)
        push!(mismatches_A, mismatch)
    end
    
    # Count by type
    all_mismatches = vcat(mismatches_A, mismatches_B)
    t_junction_count = count(m -> m.mismatch_type == T_JUNCTION, all_mismatches)
    diagonal_count = count(m -> m.mismatch_type == DIAGONAL, all_mismatches)
    refinement_count = count(m -> m.mismatch_type == REFINEMENT, all_mismatches)
    unknown_count = count(m -> m.mismatch_type == UNKNOWN, all_mismatches)
    
    # Repair assessment
    feasible_count = count(m -> m.repair_feasible, all_mismatches)
    infeasible_count = length(all_mismatches) - feasible_count
    
    avg_complexity = if !isempty(all_mismatches)
        sum(m.complexity_score for m in all_mismatches) / length(all_mismatches)
    else
        0.0
    end
    
    println("  Classification complete:")
    println("    T-junctions: $t_junction_count")
    println("    Diagonal mismatches: $diagonal_count")
    println("    Refinement differences: $refinement_count")
    println("    Unknown: $unknown_count")
    println("    Feasible for repair: $feasible_count / $(length(all_mismatches))")
    println("    Average complexity: $(round(avg_complexity, digits=2))")
    
    return InterfaceClassification(
        topology,
        mismatches_A,
        mismatches_B,
        t_junction_count,
        diagonal_count,
        refinement_count,
        unknown_count,
        feasible_count,
        infeasible_count,
        avg_complexity
    )
end

"""
    export_classification_json(classification, output_file)

Export classification report to JSON.
"""
function export_classification_json(classification::InterfaceClassification, output_file::String)
    
    function mismatch_to_dict(m::EdgeMismatch)
        return Dict(
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
