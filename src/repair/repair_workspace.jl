"""
    repair_workspace.jl

Phase 3: Surgical Mesh Repair Execution - Workspace Management

Provides transactional workspace for safe mesh modification with checkpoint/rollback support.
Tracks all modifications to enable atomic repairs with failure recovery.
"""

using LinearAlgebra
using Printf

"""
    ModificationType

Types of modifications that can be applied to a mesh.
"""
@enum ModificationType begin
    FACE_DELETION
    FACE_ADDITION
    NODE_ADDITION
    CONNECTIVITY_UPDATE
end

"""
    MeshModification

Record of a single modification to the mesh.
Used for rollback if repair fails.
"""
struct MeshModification
    mod_type::ModificationType
    timestamp::Int  # Sequential modification ID
    
    # For face operations
    pid::Union{Int,Nothing}
    face_index::Union{Int,Nothing}  # Index in original face list
    face_nodes::Union{Vector{Int},Nothing}  # Node IDs
    
    # For node operations
    node_id::Union{Int,Nothing}
    node_coords::Union{NTuple{3,Float64},Nothing}
    
    # Metadata
    description::String
end

"""
    RepairWorkspace

Transactional workspace for mesh repair operations.
Maintains checkpoint state and modification history for rollback support.
"""
mutable struct RepairWorkspace
    # Original mesh file
    original_file::String
    
    # Working copies (mutable)
    working_faces::Dict{Int, Vector{Vector{Int}}}  # PID => face list (node IDs)
    working_nodes::Dict{Int, NTuple{3,Float64}}    # Node ID => coordinates
    max_node_id::Int
    
    # Transaction tracking
    modifications::Vector{MeshModification}
    checkpoint_mod_count::Int  # Number of mods at last checkpoint
    transaction_active::Bool
    
    # Statistics
    faces_added::Int
    faces_deleted::Int
    nodes_added::Int
end

function RepairWorkspace(mesh_file::String)
        # Load mesh data using gmsh
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        
        try
            gmsh.open(mesh_file)
            
            # Get nodes
            node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
            working_nodes = Dict{Int, NTuple{3,Float64}}()
            for (i, tag) in enumerate(node_tags)
                idx = (i-1)*3
                working_nodes[Int(tag)] = (node_coords[idx+1], node_coords[idx+2], node_coords[idx+3])
            end
            
            max_node_id = maximum(keys(working_nodes))
            
            # Get boundary faces per PID
            region_tets = Dict{Int,Vector{Tuple{Int,NTuple{4,Int}}}}()
            for (dim, tag) in gmsh.model.getEntities(3)
                etypes, etags, enodes = gmsh.model.mesh.getElements(dim, tag)
                for (etype, etag_vec, enode_vec) in zip(etypes, etags, enodes)
                    if etype != 4; continue; end
                    for i in 1:length(etag_vec)
                        base = (i-1)*4
                        nd = (Int(enode_vec[base+1]), Int(enode_vec[base+2]), 
                              Int(enode_vec[base+3]), Int(enode_vec[base+4]))
                        push!(get!(region_tets, tag, Vector{Tuple{Int,NTuple{4,Int}}}()), (0, nd))
                    end
                end
            end
            
            # Build face incidence
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
            
            # Extract boundary faces (odd incidence)
            working_faces = Dict{Int, Vector{Vector{Int}}}()
            for (face, pid_counts) in face_inc
                for (pid, cnt) in pid_counts
                    if isodd(cnt)
                        face_vec = [face[1], face[2], face[3]]
                        push!(get!(working_faces, pid, Vector{Vector{Int}}()), face_vec)
                    end
                end
            end
            
            ws = RepairWorkspace(
                mesh_file,
                working_faces,
                working_nodes,
                max_node_id,
                MeshModification[],
                0,
                false,
                0, 0, 0
            )
            
            return ws
        finally
            gmsh.finalize()
        end
end

"""
    create_checkpoint!(ws::RepairWorkspace)

Create a checkpoint at the current state.
Can be rolled back to this point if needed.
"""
function create_checkpoint!(ws::RepairWorkspace)
    ws.checkpoint_mod_count = length(ws.modifications)
    @info "Checkpoint created at modification #$(ws.checkpoint_mod_count)"
end

"""
    begin_transaction!(ws::RepairWorkspace)

Begin a new transaction for a series of modifications.
"""
function begin_transaction!(ws::RepairWorkspace)
    if ws.transaction_active
        error("Transaction already active!")
    end
    ws.transaction_active = true
    create_checkpoint!(ws)
    @info "Transaction started"
end

"""
    commit_transaction!(ws::RepairWorkspace)

Commit the current transaction, making all modifications permanent.
"""
function commit_transaction!(ws::RepairWorkspace)
    if !ws.transaction_active
        error("No active transaction to commit!")
    end
    ws.transaction_active = false
    create_checkpoint!(ws)  # New checkpoint after commit
    @info "Transaction committed ($(length(ws.modifications)) total modifications)"
end

"""
    rollback_transaction!(ws::RepairWorkspace)

Rollback all modifications since the last checkpoint.
"""
function rollback_transaction!(ws::RepairWorkspace)
    if !ws.transaction_active
        error("No active transaction to rollback!")
    end
    
    # Count modifications to rollback
    mods_to_undo = length(ws.modifications) - ws.checkpoint_mod_count
    @info "Rolling back $mods_to_undo modifications..."
    
    # Undo modifications in reverse order
    for i in length(ws.modifications):-1:(ws.checkpoint_mod_count + 1)
        mod = ws.modifications[i]
        undo_modification!(ws, mod)
    end
    
    # Truncate modification list
    resize!(ws.modifications, ws.checkpoint_mod_count)
    
    # Restore max_node_id to reflect current nodes after rollback
    if isempty(ws.working_nodes)
        ws.max_node_id = 0
    else
        ws.max_node_id = maximum(keys(ws.working_nodes))
    end

    ws.transaction_active = false
    @info "Transaction rolled back successfully"
end

"""
    undo_modification!(ws::RepairWorkspace, mod::MeshModification)

Undo a single modification.
"""
function undo_modification!(ws::RepairWorkspace, mod::MeshModification)
    if mod.mod_type == FACE_DELETION
        # Re-add the deleted face
        if mod.face_index !== nothing && mod.face_nodes !== nothing && mod.pid !== nothing
            # Insert back at original position
            faces = ws.working_faces[mod.pid]
            if mod.face_index <= length(faces)
                insert!(faces, mod.face_index, mod.face_nodes)
            else
                push!(faces, mod.face_nodes)
            end
            ws.faces_deleted -= 1
        end
        
    elseif mod.mod_type == FACE_ADDITION
        # Remove the added face (should be at the end)
        if mod.pid !== nothing
            faces = ws.working_faces[mod.pid]
            if !isempty(faces)
                pop!(faces)
                ws.faces_added -= 1
            end
        end
        
    elseif mod.mod_type == NODE_ADDITION
        # Remove the added node
        if mod.node_id !== nothing
            delete!(ws.working_nodes, mod.node_id)
            ws.nodes_added -= 1
        end
    end
end

"""
    delete_face!(ws::RepairWorkspace, pid::Int, face_index::Int)

Delete a face from the working mesh.
Records modification for rollback support.
"""
function delete_face!(ws::RepairWorkspace, pid::Int, face_index::Int)
    faces = ws.working_faces[pid]
    if face_index < 1 || face_index > length(faces)
        error("Face index $face_index out of bounds for PID $pid")
    end
    
    face_nodes = faces[face_index]
    
    # Record modification
    mod = MeshModification(
        FACE_DELETION,
        length(ws.modifications) + 1,
        pid,
        face_index,
        copy(face_nodes),
        nothing, nothing,
        "Delete face #$face_index from PID $pid"
    )
    push!(ws.modifications, mod)
    
    # Delete the face
    deleteat!(faces, face_index)
    ws.faces_deleted += 1
end

"""
    add_face!(ws::RepairWorkspace, pid::Int, node_ids::Vector{Int})

Add a new face to the working mesh.
Records modification for rollback support.
"""
function add_face!(ws::RepairWorkspace, pid::Int, node_ids::Vector{Int})
    # Validate nodes exist
    for node_id in node_ids
        if !haskey(ws.working_nodes, node_id)
            error("Node $node_id does not exist in workspace")
        end
    end
    
    # Record modification
    mod = MeshModification(
        FACE_ADDITION,
        length(ws.modifications) + 1,
        pid,
        nothing,
        copy(node_ids),
        nothing, nothing,
        "Add face with nodes $node_ids to PID $pid"
    )
    push!(ws.modifications, mod)
    
    # Add the face
    push!(ws.working_faces[pid], node_ids)
    ws.faces_added += 1
end

"""
    add_node!(ws::RepairWorkspace, coords::NTuple{3,Float64}) -> Int

Add a new node to the working mesh.
Returns the new node ID.
"""
function add_node!(ws::RepairWorkspace, coords::NTuple{3,Float64})
    # Generate new node ID
    ws.max_node_id += 1
    new_node_id = ws.max_node_id
    
    # Record modification
    mod = MeshModification(
        NODE_ADDITION,
        length(ws.modifications) + 1,
        nothing,
        nothing,
        nothing,
        new_node_id,
        coords,
        "Add node $new_node_id at $coords"
    )
    push!(ws.modifications, mod)
    
    # Add the node
    ws.working_nodes[new_node_id] = coords
    ws.nodes_added += 1
    
    return new_node_id
end

"""
    get_face_by_nodes(ws::RepairWorkspace, pid::Int, node_coords::Vector{NTuple{3,Float64}}) -> Union{Int,Nothing}

Find face index in PID by matching node coordinates.
Returns face index or nothing if not found.
"""
function get_face_by_nodes(ws::RepairWorkspace, pid::Int, node_coords::Vector{NTuple{3,Float64}})
    faces = ws.working_faces[pid]
    
    for (face_idx, face_nodes) in enumerate(faces)
        if length(face_nodes) != length(node_coords)
            continue
        end
        
        # Get coordinates of face nodes
        face_coords = [ws.working_nodes[nid] for nid in face_nodes]
        
        # Check if coordinates match (order-independent)
        match = true
        for coord in node_coords
            if !any(c -> all(abs.(c .- coord) .< 1e-6), face_coords)
                match = false
                break
            end
        end
        
        if match
            return face_idx
        end
    end
    
    return nothing
end

"""
    get_node_id_by_coords(ws::RepairWorkspace, coords::NTuple{3,Float64}) -> Union{Int,Nothing}

Find node ID by coordinates.
Returns node ID or nothing if not found.
"""
function get_node_id_by_coords(ws::RepairWorkspace, coords::NTuple{3,Float64})
    for (node_id, node_coords) in ws.working_nodes
        if all(abs.(node_coords .- coords) .< 1e-6)
            return node_id
        end
    end
    return nothing
end

"""
    export_modified_mesh(ws::RepairWorkspace)

Export the modified mesh as a new mesh-like object compatible with write_nastran.
"""
function export_modified_mesh(ws::RepairWorkspace)
    # Create new mesh with updated data
    new_mesh = deepcopy(ws.original_mesh)
    
    # Update faces
    new_mesh.all_pid_surfaces = ws.working_faces
    
    # Update nodes
    new_mesh.nodes = ws.working_nodes
    
    # Rebuild face_to_pid map
    new_mesh.face_to_pid = Dict{Vector{Int}, Int}()
    for (pid, faces) in ws.working_faces
        for face in faces
            new_mesh.face_to_pid[face] = pid
        end
    end
    
    return new_mesh
end

"""
    print_workspace_stats(ws::RepairWorkspace)

Print statistics about the workspace state.
"""
function print_workspace_stats(ws::RepairWorkspace)
    println("\n" * "="^70)
    println("Repair Workspace Statistics")
    println("="^70)
    println("Modifications:   $(length(ws.modifications))")
    println("Checkpoint at:   Modification #$(ws.checkpoint_mod_count)")
    println("Transaction:     $(ws.transaction_active ? "Active" : "Inactive")")
    println()
    println("Changes:")
    println("  Faces added:   $(ws.faces_added)")
    println("  Faces deleted: $(ws.faces_deleted)")
    println("  Nodes added:   $(ws.nodes_added)")
    println()
    println("Current state:")
    total_faces = sum(length(faces) for faces in values(ws.working_faces))
    println("  Total faces:   $total_faces")
    println("  Total nodes:   $(length(ws.working_nodes))")
    println("  Max node ID:   $(ws.max_node_id)")
    println("="^70)
end

# ============================================================================
# Phase 6: Mesh Replacement Operations
# ============================================================================

"""
    get_or_create_node!(ws::RepairWorkspace, coords::NTuple{3,Float64}; tol::Real=1e-6) -> Int

Get existing node ID by coordinates, or create a new node if not found.
Uses tolerance-based matching to handle coordinate rounding.
Returns the node ID (existing or newly created).
"""
function get_or_create_node!(ws::RepairWorkspace, coords::NTuple{3,Float64}; tol::Real=1e-6)
    # Try to find existing node within tolerance
    for (node_id, node_coords) in ws.working_nodes
        dx = node_coords[1] - coords[1]
        dy = node_coords[2] - coords[2]
        dz = node_coords[3] - coords[3]
        dist2 = dx*dx + dy*dy + dz*dz
        if dist2 <= tol*tol
            return node_id
        end
    end
    
    # Node not found, create new one
    return add_node!(ws, coords)
end

"""
    delete_interface_face!(ws::RepairWorkspace, pid::Int, face_index::Int) -> Bool

Delete an interface face from the workspace.
This is a wrapper around delete_face! with additional validation.
Returns true on success, false on failure.
"""
function delete_interface_face!(ws::RepairWorkspace, pid::Int, face_index::Int)
    if !haskey(ws.working_faces, pid)
        @warn "PID $pid not found in workspace"
        return false
    end
    
    faces = ws.working_faces[pid]
    if face_index < 1 || face_index > length(faces)
        @warn "Face index $face_index out of bounds for PID $pid (has $(length(faces)) faces)"
        return false
    end
    
    delete_face!(ws, pid, face_index)
    return true
end

"""
    verify_mesh_integrity(ws::RepairWorkspace, pidA::Int, pidB::Int) -> Bool

Verify mesh integrity after interface replacement.
Checks:
- All faces have valid node references
- All nodes exist in workspace
- No degenerate triangles
- Basic manifoldness (each edge shared by at most 2 triangles)

Returns true if mesh is valid, false otherwise.
"""
function verify_mesh_integrity(ws::RepairWorkspace, pidA::Int, pidB::Int)
    for pid in [pidA, pidB]
        if !haskey(ws.working_faces, pid)
            @error "PID $pid not found in workspace"
            return false
        end
        
        faces = ws.working_faces[pid]
        
        # Check each face
        for (face_idx, face_nodes) in enumerate(faces)
            # Check face has exactly 3 nodes
            if length(face_nodes) != 3
                @error "Face $face_idx in PID $pid has $(length(face_nodes)) nodes (expected 3)"
                return false
            end
            
            # Check all nodes exist
            for node_id in face_nodes
                if !haskey(ws.working_nodes, node_id)
                    @error "Node $node_id in face $face_idx of PID $pid does not exist"
                    return false
                end
            end
            
            # Check for degenerate triangle (repeated nodes)
            if length(unique(face_nodes)) != 3
                @error "Face $face_idx in PID $pid has degenerate triangle (repeated nodes)"
                return false
            end
            
            # Check triangle is not zero-area
            coords = [ws.working_nodes[nid] for nid in face_nodes]
            v1 = (coords[2][1] - coords[1][1], coords[2][2] - coords[1][2], coords[2][3] - coords[1][3])
            v2 = (coords[3][1] - coords[1][1], coords[3][2] - coords[1][2], coords[3][3] - coords[1][3])
            
            nx = v1[2] * v2[3] - v1[3] * v2[2]
            ny = v1[3] * v2[1] - v1[1] * v2[3]
            nz = v1[1] * v2[2] - v1[2] * v2[1]
            
            area = sqrt(nx^2 + ny^2 + nz^2) / 2.0
            
            if area < 1e-12
                @error "Face $face_idx in PID $pid has zero area"
                return false
            end
        end
    end
    
    # ========================================================================
    # CRITICAL: Non-Manifold Edge Analysis
    # ========================================================================
    # A manifold edge is incident to exactly 1 (boundary) or 2 (interior) faces.
    # Non-manifold edges (>2 faces) indicate topology problems that prevent
    # downstream operations like remeshing, offset, or boolean operations.
    #
    # Classification:
    # - Boundary edges (1 face): Normal for surface boundaries
    # - Manifold edges (2 faces): Normal interior edges
    # - Non-manifold edges (>2 faces): PROBLEMS - T-junctions, overlaps, etc.
    #
    # For interface repair:
    # - Interface edges: Expected to transition from 1→2 or 2→2 (conforming)
    # - Non-manifold at interface: Often caused by misaligned meshes
    # - Non-manifold within PID: Indicates mesh quality issues
    
    edge_pid_incidence = Dict{Tuple{Int,Int}, Dict{Int,Int}}()
    
    for pid in (pidA, pidB)
        faces = ws.working_faces[pid]
        for face_nodes in faces
            # Create edges (sorted pairs)
            edges = [
                (min(face_nodes[1], face_nodes[2]), max(face_nodes[1], face_nodes[2])),
                (min(face_nodes[2], face_nodes[3]), max(face_nodes[2], face_nodes[3])),
                (min(face_nodes[3], face_nodes[1]), max(face_nodes[3], face_nodes[1]))
            ]
            
            for edge in edges
                d = get!(edge_pid_incidence, edge, Dict{Int,Int}())
                d[pid] = get(d, pid, 0) + 1
            end
        end
    end
    
    # Analyze edge topology
    total_edges = length(edge_pid_incidence)
    boundary_edges = 0    # Incident to 1 face
    manifold_edges = 0    # Incident to 2 faces  
    nonmanifold_edges = 0 # Incident to >2 faces
    
    # Categorize non-manifold edges by severity and location
    nm_within_A = 0      # Non-manifold entirely within PID A
    nm_within_B = 0      # Non-manifold entirely within PID B
    nm_interface = 0     # Non-manifold at interface (both PIDs involved)
    max_incidence = 0
    worst_edge = nothing
    
    incidence_histogram = Dict{Int,Int}()  # incidence -> count
    
    for (edge, perpid) in edge_pid_incidence
        a_cnt = get(perpid, pidA, 0)
        b_cnt = get(perpid, pidB, 0)
        total = a_cnt + b_cnt
        
        if total == 1
            boundary_edges += 1
        elseif total == 2
            manifold_edges += 1
        else  # total > 2: non-manifold
            nonmanifold_edges += 1
            incidence_histogram[total] = get(incidence_histogram, total, 0) + 1
            
            # Categorize by location
            if a_cnt > 0 && b_cnt > 0
                nm_interface += 1
            elseif a_cnt > 0
                nm_within_A += 1
            else
                nm_within_B += 1
            end
            
            # Track worst case
            if total > max_incidence
                max_incidence = total
                worst_edge = (edge, a_cnt, b_cnt)
            end
        end
    end
    
    # Report non-manifold topology
    if nonmanifold_edges > 0
        pct = round(100 * nonmanifold_edges / total_edges, digits=2)
        
        println("\n  ⚠ Non-Manifold Topology Detected:")
        println("    Total edges analyzed: $total_edges")
        println("    └─ Boundary (1 face): $boundary_edges")
        println("    └─ Manifold (2 faces): $manifold_edges")
        println("    └─ NON-MANIFOLD (>2 faces): $nonmanifold_edges ($pct%)")
        
        println("\n  Non-Manifold Breakdown:")
        println("    Interface (both PIDs): $nm_interface  ← Interface misalignment")
        println("    Within PID $pidA only:  $nm_within_A  ← Internal mesh issues")
        println("    Within PID $pidB only:  $nm_within_B  ← Internal mesh issues")
        
        # Show incidence distribution
        if !isempty(incidence_histogram)
            sorted_inc = sort(collect(incidence_histogram))
            println("\n  Incidence Distribution:")
            for (inc, count) in sorted_inc
                println("    $inc faces → $count edges")
            end
        end
        
        # Show worst case
        if worst_edge !== nothing && max_incidence > 4
            e, ac, bc = worst_edge
            println("\n  ⚠ Worst case: Edge $e with $max_incidence faces (PID A: $ac, PID B: $bc)")
        end
        
        # Actionable guidance
        println("\n  Impact: Non-manifold edges may cause:")
        println("    • Remeshing failures")
        println("    • Invalid boolean operations")
        println("    • Incorrect volume calculations")
        println("    • Visualization artifacts")
        
        if nm_interface > nonmanifold_edges * 0.8
            println("\n  ✓ Most non-manifold edges at interface - expected during repair")
        elseif nm_within_A > 0 || nm_within_B > 0
            println("\n  ✗ Non-manifold edges within PIDs - indicates input mesh quality issues")
        end
        
        println()  # Blank line for readability
    end
    
    return true
end

"""
    map_nodes_to_pid(triangle::Triangle, node_mapping::Dict{NTuple{3,Float64}, Union{Int, Nothing}}, 
                     ws::RepairWorkspace) -> Union{Vector{Int}, Nothing}

Map unified mesh triangle coordinates to workspace node IDs.
Uses the node_mapping from UnifiedInterfaceMesh to find existing nodes.
Creates new nodes if needed.

Returns:
- Vector{Int}: Node IDs for the triangle [node1_id, node2_id, node3_id]
- Nothing: If mapping fails
"""
function map_nodes_to_pid(triangle::Triangle, 
                          node_mapping::Dict{NTuple{3,Float64}, Union{Int, Nothing}},
                          ws::RepairWorkspace)
    node_coords = [triangle.coord1, triangle.coord2, triangle.coord3]
    node_ids = Int[]
    
    for coord in node_coords
        # First try the node mapping
        if haskey(node_mapping, coord)
            mapped_id = node_mapping[coord]
            if mapped_id !== nothing
                push!(node_ids, mapped_id)
                continue
            end
        end
        
        # Fallback: get or create node in workspace
        node_id = get_or_create_node!(ws, coord)
        push!(node_ids, node_id)
    end
    
    if length(node_ids) != 3
        @error "Failed to map all triangle nodes"
        return nothing
    end
    
    return node_ids
end
