# Boundary constraint identification for safe interface repair
# Phase 1.3: Identify boundaries that must NOT be modified

# Import CoordinateKeys module for consistent EdgeKey handling
using .CoordinateKeys: create_edge_key_int, EdgeKeyInt, convert_to_int, convert_to_float

# ============================================================================
# Boundary constraint structures
# ============================================================================

"""
Constraints that must be preserved during interface repair.
"""
struct BoundaryConstraints
    pidA::Int
    pidB::Int
    
    # External boundary edges (interface with other PIDs)
    pidA_external_edges::Set{EdgeKey}  # edges on boundary with PID ≠ B
    pidB_external_edges::Set{EdgeKey}
    
    # Corner nodes (junction of 3+ PIDs)
    corner_nodes::Set{NTuple{3,Float64}}
    
    # Feature edges (high curvature or user-specified)
    feature_edges::Set{EdgeKey}
    
    # Lock status
    locked_nodes::Set{NTuple{3,Float64}}  # CANNOT move/delete
    locked_edges::Set{EdgeKey}             # CANNOT modify
    
    # Statistics
    total_external_edges_A::Int
    total_external_edges_B::Int
    total_corner_nodes::Int
    total_locked_edges::Int
end

# ============================================================================
# Constraint identification
# ============================================================================

"""
    build_boundary_constraints(nas_file, pidA, pidB; tol=1e-4, feature_angle=30.0)

Identify all boundary constraints for the interface between pidA and pidB.
"""
function build_boundary_constraints(nas_file::String, pidA::Int, pidB::Int;
                                   tol::Real=1e-4, feature_angle::Real=30.0)::BoundaryConstraints
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    
    try
        gmsh.open(nas_file)
        
        # Get all nodes and coordinates
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        coords = Dict{Int,NTuple{3,Float64}}()
        for (i, tag) in enumerate(node_tags)
            idx = (i-1)*3
            coords[Int(tag)] = (node_coords[idx+1], node_coords[idx+2], node_coords[idx+3])
        end
        
        # Helper: coordinate key
        function ckey(p::NTuple{3,Float64})
            return (round(p[1]; digits=4), round(p[2]; digits=4), round(p[3]; digits=4))
        end
        
        # Collect all PIDs and their tetrahedra
        region_tets = Dict{Int,Vector{Tuple{Int,NTuple{4,Int}}}}()
        for (dim, tag) in gmsh.model.getEntities(3)
            etypes, etags, enodes = gmsh.model.mesh.getElements(dim, tag)
            for (etype, etag_vec, enode_vec) in zip(etypes, etags, enodes)
                if etype != 4; continue; end
                for i in 1:length(etag_vec)
                    eid = Int(etag_vec[i])
                    base = (i-1)*4
                    nd = (Int(enode_vec[base+1]), Int(enode_vec[base+2]), 
                          Int(enode_vec[base+3]), Int(enode_vec[base+4]))
                    push!(get!(region_tets, tag, Vector{Tuple{Int,NTuple{4,Int}}}()), (eid, nd))
                end
            end
        end
        
        all_pids = collect(keys(region_tets))
        
        # Build face incidence for all PIDs
        face_inc = Dict{NTuple{3,Int},Dict{Int,Int}}()
        tet_faces = ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        
        for (pid, tets) in region_tets
            for (_, nd) in tets
                n = (nd[1], nd[2], nd[3], nd[4])
                for f in tet_faces
                    tri = (n[f[1]], n[f[2]], n[f[3]])
                    tri_sorted = tuple(sort!(collect(tri))...)
                    d = get!(face_inc, tri_sorted, Dict{Int,Int}())
                    d[pid] = get(d, pid, 0) + 1
                end
            end
        end
        
        # Get boundary faces per PID
        boundary_faces = Dict{Int,Vector{NTuple{3,Int}}}()
        for (face, pid_counts) in face_inc
            for (pid, cnt) in pid_counts
                if isodd(cnt)
                    push!(get!(boundary_faces, pid, NTuple{3,Int}[]), face)
                end
            end
        end
        
        # Build node→PID membership
        node_to_pids = Dict{NTuple{3,Float64},Set{Int}}()
        for (pid, tets) in region_tets
            for (_, nd) in tets
                for nid in nd
                    k = ckey(coords[nid])
                    push!(get!(node_to_pids, k, Set{Int}()), pid)
                end
            end
        end
        
        # Find corner nodes (shared by 3+ PIDs)
        corner_nodes = Set{NTuple{3,Float64}}()
        for (node_key, pids) in node_to_pids
            if length(pids) >= 3
                push!(corner_nodes, node_key)
            end
        end
        
        # Build edge→PIDs map from boundary faces
        edge_to_pids = Dict{EdgeKey,Set{Int}}()
        for (pid, faces) in boundary_faces
            for face in faces
                # Use factory function for consistent EdgeKey creation
                c1 = coords[face[1]]
                c2 = coords[face[2]]
                c3 = coords[face[3]]

                ek1 = create_edge_key_int(c1, c2)
                ek2 = create_edge_key_int(c1, c3)
                ek3 = create_edge_key_int(c2, c3)

                push!(get!(edge_to_pids, ek1, Set{Int}()), pid)
                push!(get!(edge_to_pids, ek2, Set{Int}()), pid)
                push!(get!(edge_to_pids, ek3, Set{Int}()), pid)
            end
        end
        
        # External edges for pidA (touches PID ≠ pidB)
        external_edges_A = Set{EdgeKey}()
        for (edge, pids) in edge_to_pids
            if pidA in pids
                other_pids = setdiff(pids, Set([pidA]))
                if !isempty(other_pids) && pidB ∉ other_pids
                    push!(external_edges_A, edge)
                end
            end
        end
        
        # External edges for pidB (touches PID ≠ pidA)
        external_edges_B = Set{EdgeKey}()
        for (edge, pids) in edge_to_pids
            if pidB in pids
                other_pids = setdiff(pids, Set([pidB]))
                if !isempty(other_pids) && pidA ∉ other_pids
                    push!(external_edges_B, edge)
                end
            end
        end
        
        # Detect feature edges (high curvature) - simplified version
        # For now, feature edges are edges shared by exactly 2 triangles with angle > threshold
        feature_edges = Set{EdgeKey}()
        # TODO: Implement proper feature edge detection based on dihedral angles
        
        # Locked nodes: corner nodes
        locked_nodes = copy(corner_nodes)
        
        # Locked edges: external edges + feature edges
        locked_edges = union(external_edges_A, external_edges_B, feature_edges)
        
        return BoundaryConstraints(
            pidA, pidB,
            external_edges_A,
            external_edges_B,
            corner_nodes,
            feature_edges,
            locked_nodes,
            locked_edges,
            length(external_edges_A),
            length(external_edges_B),
            length(corner_nodes),
            length(locked_edges)
        )
        
    finally
        gmsh.finalize()
    end
end

"""
    check_constraint_violations(edge, constraints)

Check if modifying an edge would violate any constraints.
Returns (has_violation, violation_reasons).
"""
function check_constraint_violations(edge::EdgeKey, constraints::BoundaryConstraints)
    violations = String[]
    
    # Check if edge is locked
    if edge in constraints.locked_edges
        push!(violations, "Edge is locked (external or feature)")
    end
    
    # Check if endpoints are locked
    # Convert EdgeKey integer coordinates to compare with locked_nodes
    # Using user's preferred approach: convert to integer and compare element-wise
    edge_int1 = edge.node1
    edge_int2 = edge.node2

    # Convert locked_nodes to integer coordinates for comparison
    # locked_nodes contains float coordinates, so we convert them to integer keys
    for locked_node in constraints.locked_nodes
        locked_node_int = convert_to_int(locked_node)
        if locked_node_int == edge_int1
            push!(violations, "Endpoint node1 is locked (corner node)")
        end
        if locked_node_int == edge_int2
            push!(violations, "Endpoint node2 is locked (corner node)")
        end
    end
    
    has_violation = !isempty(violations)
    return (has_violation, violations)
end

"""
    export_constraints_json(constraints, output_file)

Export boundary constraints to JSON for inspection.
"""
function export_constraints_json(constraints::BoundaryConstraints, output_file::String)
    report = Dict(
        "interface" => Dict(
            "pidA" => constraints.pidA,
            "pidB" => constraints.pidB
        ),
        "summary" => Dict(
            "external_edges_A" => constraints.total_external_edges_A,
            "external_edges_B" => constraints.total_external_edges_B,
            "corner_nodes" => constraints.total_corner_nodes,
            "total_locked_edges" => constraints.total_locked_edges
        ),
        "details" => Dict(
            "corner_nodes_sample" => [[n...] for n in collect(constraints.corner_nodes)[1:min(20, length(constraints.corner_nodes))]],
            "external_edges_A_sample" => [[[convert_to_float(e.node1)...], [convert_to_float(e.node2)...]] for e in collect(constraints.pidA_external_edges)[1:min(20, length(constraints.pidA_external_edges))]],
            "external_edges_B_sample" => [[[convert_to_float(e.node1)...], [convert_to_float(e.node2)...]] for e in collect(constraints.pidB_external_edges)[1:min(20, length(constraints.pidB_external_edges))]]
        )
    )
    
    open(output_file, "w") do io
        write_json(io, report, 0)
    end
    
    println("Boundary constraints exported to: $output_file")
    return output_file
end
