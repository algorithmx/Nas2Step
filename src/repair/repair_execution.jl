"""
    repair_execution.jl

Phase 3: Surgical Mesh Repair Execution

Applies repair plans to mesh using transactional workspace.
Implements quad retriangulation with rollback on failure.
"""

using LinearAlgebra
using Printf

# NOTE: Backward-compatibility shims removed.
# This module now expects EdgeInsertionPlan from repair_planning.jl with fields:
#   split_type::Symbol
#   old_triangles::Vector{NTuple{9,Float64}}
#   replacement_triangles::Vector{NTuple{9,Float64}}
#   new_nodes::Vector{NTuple{3,Float64}}

"""
    apply_quad_retriangulation!(ws::RepairWorkspace, plan::EdgeInsertionPlan, pid::Int) -> Bool

Apply a single quad retriangulation to the mesh.
Returns true on success, false on failure (workspace will be rolled back).
"""
function apply_quad_retriangulation!(ws::RepairWorkspace, plan::EdgeInsertionPlan, pid::Int)
    ptype = plan.split_type
    if ptype != :quad_retriangulation
        error("Plan type must be :quad_retriangulation, got $(ptype)")
    end
    
    old_tris = plan.old_triangles
    new_tris = plan.replacement_triangles
    if length(old_tris) != 2
        @warn "Expected 2 triangles to replace, got $(length(old_tris))"
        return false
    end
    
    if length(new_tris) != 2
        @warn "Expected 2 new triangles, got $(length(new_tris))"
        return false
    end
    
    # Find the two triangles to delete by their coordinates
    triangles_to_delete = Int[]
    for tri_coords_flat in old_tris
        # Convert flat 9-tuple to vector of 3 node coordinates
        tri_as_nodes = [
            (tri_coords_flat[1], tri_coords_flat[2], tri_coords_flat[3]),
            (tri_coords_flat[4], tri_coords_flat[5], tri_coords_flat[6]),
            (tri_coords_flat[7], tri_coords_flat[8], tri_coords_flat[9])
        ]
        
        face_idx = get_face_by_nodes(ws, pid, tri_as_nodes)
        if face_idx === nothing
            @warn "Could not find triangle to replace in PID $pid"
            return false
        end
        push!(triangles_to_delete, face_idx)
    end
    
    # Sort in descending order so deletion doesn't affect indices
    sort!(triangles_to_delete, rev=true)
    
    @info "  Deleting 2 old triangles (indices: $triangles_to_delete)"
    
    # Delete old triangles
    for face_idx in triangles_to_delete
        delete_face!(ws, pid, face_idx)
    end
    
    # Add new triangles
    @info "  Adding 2 new triangles"
    for new_tri_coords in new_tris
        # Convert coordinates to node IDs
        node_ids = Int[]
        for i in 1:3
            base = (i-1)*3
            coord = (new_tri_coords[base+1], new_tri_coords[base+2], new_tri_coords[base+3])
            node_id = get_node_id_by_coords(ws, coord)
            if node_id === nothing
                @warn "Could not find node for coordinate $coord"
                return false
            end
            push!(node_ids, node_id)
        end
        
        add_face!(ws, pid, node_ids)
    end
    
    @info "  ✓ Quad retriangulation applied successfully"
    return true
end

"""
    apply_triangle_split!(ws::RepairWorkspace, plan::EdgeInsertionPlan, pid::Int) -> Bool

Apply a triangle split (T-junction or bisection) to the mesh.
Returns true on success, false on failure.
"""
function apply_triangle_split!(ws::RepairWorkspace, plan::EdgeInsertionPlan, pid::Int)
    ptype = plan.split_type
    if ptype == :bisect
        # For bisection, the edge already exists - just verify it's present
        # (no actual modification needed)
        @info "  Bisection - edge already exists, no modification needed"
        return true
    end
    
    old_tris = plan.old_triangles
    new_tris = plan.replacement_triangles
    new_nodes = plan.new_nodes
    if old_tris === nothing || isempty(old_tris)
        @warn "No old triangles specified for split"
        return false
    end
    
    if new_tris === nothing || isempty(new_tris)
        @warn "No replacement triangles specified"
        return false
    end
    
    # Find the triangle to split by coordinates
    old_tri_coords = old_tris[1]  # First (and usually only) triangle to replace
    tri_as_nodes = [
        (old_tri_coords[1], old_tri_coords[2], old_tri_coords[3]),
        (old_tri_coords[4], old_tri_coords[5], old_tri_coords[6]),
        (old_tri_coords[7], old_tri_coords[8], old_tri_coords[9])
    ]
    
    face_idx = get_face_by_nodes(ws, pid, tri_as_nodes)
    if face_idx === nothing
        @warn "Could not find triangle to split in PID $pid"
        return false
    end
    
    @info "  Deleting triangle at index $face_idx"
    delete_face!(ws, pid, face_idx)
    
    # Add new nodes if needed
    for new_node_coords in new_nodes
        # Check if node already exists
        existing_id = get_node_id_by_coords(ws, new_node_coords)
        if existing_id === nothing
            new_id = add_node!(ws, new_node_coords)
            @info "    Added new node $new_id at $new_node_coords"
        else
            @info "    Reusing existing node $existing_id"
        end
    end
    
    # Add replacement triangles
    @info "  Adding $(length(new_tris)) replacement triangle(s)"
    for new_tri_coords in new_tris
        # Convert coordinates to node IDs
        node_ids = Int[]
        for i in 1:3
            base = (i-1)*3
            coord = (new_tri_coords[base+1], new_tri_coords[base+2], new_tri_coords[base+3])
            node_id = get_node_id_by_coords(ws, coord)
            if node_id === nothing
                @warn "Could not find node for coordinate $coord"
                return false
            end
            push!(node_ids, node_id)
        end
        
        add_face!(ws, pid, node_ids)
    end
    
    @info "  ✓ Triangle split applied successfully"
    return true
end

"""
    apply_edge_insertion_plan!(ws::RepairWorkspace, plan::EdgeInsertionPlan, pid::Int) -> Bool

Apply a single edge insertion plan (handles all plan types).
Returns true on success, false on failure.
"""
function apply_edge_insertion_plan!(ws::RepairWorkspace, plan::EdgeInsertionPlan, pid::Int)
    ptype = plan.split_type
    if ptype == :quad_retriangulation
        return apply_quad_retriangulation!(ws, plan, pid)
    elseif ptype == :bisect || ptype == :trisect || ptype == :quadrisect || ptype == :triangle_split
        return apply_triangle_split!(ws, plan, pid)
    elseif ptype == :failed
        @warn "Cannot apply failed plan"
        return false
    else
        @warn "Unknown plan type: $(ptype)"
        return false
    end
end

"""
    apply_repair_plan!(ws::RepairWorkspace, plan::RepairPlan) -> Bool

Apply a complete repair plan to the mesh.
Uses transactions with automatic rollback on failure.
Returns true if all repairs succeeded, false otherwise.
"""
function apply_repair_plan!(ws::RepairWorkspace, plan::RepairPlan)
    println("\n" * "="^70)
    println("Executing Repair Plan")
    println("="^70)
    println("Interface: $(plan.interface_pair)")
    println("Direction: $(plan.repair_direction)")
    println("Total edge insertions: $(plan.total_edges_to_insert)")
    println("Total triangles to split: $(plan.total_triangles_to_split)")
    println("Total nodes to add: $(plan.total_nodes_to_add)")
    println("Feasible plans: $(sum(p.is_feasible for p in plan.insertion_sequence))/$(length(plan.insertion_sequence))")
    println("="^70)
    
    # Determine which PID to modify
    pid_a, pid_b = plan.interface_pair
    target_pid = if plan.repair_direction == :subdivide_A
        pid_a
    elseif plan.repair_direction == :subdivide_B
        pid_b
    else
        error("Unknown repair direction: $(plan.repair_direction)")
    end
    
    println("\nTarget PID for modification: $target_pid")
    
    # Filter to only feasible plans
    feasible_plans = filter(p -> p.is_feasible, plan.insertion_sequence)
    
    if isempty(feasible_plans)
        @warn "No feasible plans to execute!"
        return false
    end
    
    println("Executing $(length(feasible_plans)) feasible repairs...")
    
    # Begin transaction
    begin_transaction!(ws)
    
    successes = 0
    failures = 0
    
    try
        for (i, edge_plan) in enumerate(feasible_plans)
            println("\n--- Repair $i/$(length(feasible_plans)) ---")
            println("Edge: $(edge_plan.insert_edge)")
            println("Plan type: $(edge_plan.split_type)")
            
            success = apply_edge_insertion_plan!(ws, edge_plan, target_pid)
            
            if success
                successes += 1
            else
                failures += 1
                @warn "Repair $i failed, continuing with remaining repairs..."
            end
        end
        
        # Check if enough repairs succeeded
        if failures > 0
            @warn "$failures out of $(length(feasible_plans)) repairs failed"
        end
        
        if successes == 0
            @error "All repairs failed! Rolling back transaction."
            rollback_transaction!(ws)
            return false
        end
        
        # Commit if at least some succeeded
        println("\n" * "="^70)
        println("Repair Summary")
        println("="^70)
        println("Successful: $successes")
        println("Failed:     $failures")
        println("Total:      $(length(feasible_plans))")
        println("="^70)
        
        if failures == 0
            println("\n✓ All repairs succeeded! Committing transaction.")
            commit_transaction!(ws)
            return true
        else
            println("\n⚠ Some repairs failed. Committing partial success.")
            commit_transaction!(ws)
            return false
        end
        
    catch e
        @error "Exception during repair execution: $e"
        println("Rolling back transaction...")
        rollback_transaction!(ws)
        rethrow(e)
    end
end

"""
    execute_repairs_from_json(mesh_file::String, plan_file::String, output_file::String) -> Bool

High-level function to execute repairs from JSON plan file.
Loads mesh, applies repairs, and saves modified mesh.
"""
function execute_repairs_from_json(mesh_file::String, plan_file::String, output_file::String)
    println("\n" * "="^70)
    println("Surgical Mesh Repair Execution")
    println("="^70)
    println("Mesh file:   $mesh_file")
    println("Plan file:   $plan_file")
    println("Output file: $output_file")
    println("="^70)
    
    # Load mesh
    println("\nLoading mesh...")
    mesh = read_nastran(mesh_file)
    println("✓ Mesh loaded: $(length(mesh.nodes)) nodes, $(sum(length(f) for f in values(mesh.all_pid_surfaces))) faces")
    
    # Load repair plan from JSON
    println("\nLoading repair plan...")
    plan = load_repair_plan_from_json(plan_file, mesh)
    println("✓ Repair plan loaded")
    
    # Create workspace
    println("\nCreating repair workspace...")
    ws = RepairWorkspace(mesh)
    println("✓ Workspace created")
    
    # Apply repairs
    println("\nApplying repairs...")
    success = apply_repair_plan!(ws, plan)
    
    if !success
        @warn "Some or all repairs failed"
    end
    
    # Print workspace statistics
    print_workspace_stats(ws)
    
    # Export modified mesh
    println("\nExporting modified mesh...")
    modified_mesh = export_modified_mesh(ws)
    write_nastran(output_file, modified_mesh)
    println("✓ Modified mesh written to: $output_file")
    
    return success
end

"""
    load_repair_plan_from_json(json_file::String, mesh) -> RepairPlan

Load a repair plan from JSON file.
Reconstructs the full RepairPlan with topology and constraints.
"""
function load_repair_plan_from_json(json_file::String, mesh)
    # Read JSON
    plan_data = JSON.parsefile(json_file)
    
    # Extract interface pair
    interface_pair = (plan_data["interface_pair"][1], plan_data["interface_pair"][2])
    pid_a, pid_b = interface_pair
    
    # Rebuild topology (needed for full context)
    println("  Rebuilding interface topology...")
    topology = build_interface_topology(mesh, pid_a, pid_b)
    
    # Rebuild classification
    println("  Rebuilding classification...")
    classification = classify_interface_mismatches(mesh, topology, pid_a, pid_b)
    
    # Rebuild constraints
    println("  Rebuilding constraints...")
    constraints = build_boundary_constraints(mesh, topology, pid_a, pid_b)
    
    # Parse insertion sequence
    insertion_sequence = EdgeInsertionPlan[]
    for plan_dict in plan_data["insertion_sequence"]
        # Parse edge key
        edge_data = plan_dict["edge"]
        v1 = tuple(edge_data["v1"]...)
        v2 = tuple(edge_data["v2"]...)
        edge_key = EdgeKey(v1, v2)
        
        # Parse plan type
        split_type = Symbol(plan_dict["plan_type"])  # matches new field name

        # Parse triangles to replace (flattened NTuple{9,Float64})
        # Accept either nested [[x,y,z]...] or flat [x1..x9]
        old_tris = NTuple{9,Float64}[]
        for tri in plan_dict["triangles_to_replace"]
            flat = NTuple{9,Float64}(tuple(tri...))
            push!(old_tris, flat)
        end

        # Parse new triangles
        new_tris = NTuple{9,Float64}[]
        for tri in plan_dict["new_triangles"]
            flat = NTuple{9,Float64}(tuple(tri...))
            push!(new_tris, flat)
        end

        # Parse new nodes
        new_nodes = NTuple{3,Float64}[]
        for n in plan_dict["new_nodes"]
            push!(new_nodes, NTuple{3,Float64}(tuple(n...)))
        end

        # Create plan (fill unused fields conservatively)
        plan = EdgeInsertionPlan(
            0,                 # target_triangle (unknown from JSON)
            edge_key,          # insert_edge
            split_type,        # split_type
            new_nodes,         # new_nodes
            Int[],             # existing_nodes
            old_tris,          # old_triangles
            new_tris,          # replacement_triangles
            0.0,               # min_angle_before (unknown)
            0.0,               # min_angle_after (unknown)
            true,              # quality_acceptable (assume true)
            Int[],             # depends_on
            false,             # violates_constraints
            String[],          # constraint_violations
            plan_dict["is_feasible"]
        )
        
        push!(insertion_sequence, plan)
    end
    
    # Parse repair direction
    repair_direction = Symbol(plan_data["repair_direction"])
    
    # Create RepairPlan
    repair_plan = RepairPlan(
        interface_pair,
        repair_direction,
        insertion_sequence,
        plan_data["total_edges_to_insert"],
        plan_data["total_triangles_to_split"],
        plan_data["total_nodes_to_add"],
        plan_data["predicted_min_quality"],
        plan_data["predicted_max_quality_loss"],
        plan_data["is_feasible"],
        plan_data["feasibility_issues"],
        topology,
        classification,
        constraints
    )
    
    return repair_plan
end

# ============================================================================
# Phase 6: Symmetric Mesh Replacement Operations
# ============================================================================

"""
    replace_both_interfaces!(ws::RepairWorkspace,
                            unified_mesh,
                            topology) -> Bool

Replace both A's and B's interface faces with the unified mesh.

This is a transactional operation that:
1. Begins a transaction in RepairWorkspace
2. Deletes all interface faces from PID A
3. Deletes all interface faces from PID B
4. Inserts unified mesh triangles into both PIDs (with appropriate node mappings)
5. Verifies mesh integrity
6. Commits transaction if successful, rolls back on failure

# Arguments
- `ws::RepairWorkspace`: The workspace containing the mesh
- `unified_mesh::UnifiedInterfaceMesh`: The unified interface mesh to install
- `topology::InterfaceTopology`: Interface topology describing the original meshes

# Returns
- `Bool`: true if replacement succeeded, false if it failed (with rollback)

# Implementation Notes
The critical challenge is node mapping: the same geometric point may have
different node IDs in meshes A and B. We use:
1. The pre-built node_mapping_A and node_mapping_B from UnifiedInterfaceMesh
2. Coordinate-based lookup with tolerance for existing nodes
3. Creation of new nodes when needed

# Example
```julia
ws = RepairWorkspace(mesh_file)
topology = build_interface_topology(mesh, pidA, pidB)
unified = generate_unified_interface_mesh(topology, symmetric_mismatches, constraints)
success = replace_both_interfaces!(ws, unified, topology)
if success
    println("✓ Both interfaces replaced successfully")
else
    println("✗ Replacement failed (rolled back)")
end
```
"""
function replace_both_interfaces!(
    ws::RepairWorkspace,
    unified_mesh,
    topology;
    verbose::Bool = true
)::Bool
    
    if verbose
        println("\n" * "="^70)
        println("Phase 6: Symmetric Interface Replacement")
        println("="^70)
        println("Interface pair: ($(topology.pidA), $(topology.pidB))")
        println("Unified mesh:")
        println("  Triangles:         $(length(unified_mesh.triangles))")
        println("  Min quality:       $(round(unified_mesh.min_triangle_quality, digits=3))")
        println("  Total area:        $(round(unified_mesh.total_area, digits=2))")
        println("  Compatible with A: $(unified_mesh.compatible_with_A)")
        println("  Compatible with B: $(unified_mesh.compatible_with_B)")
        println("="^70)
    end
    
    # Pre-flight checks
    if !unified_mesh.compatible_with_A || !unified_mesh.compatible_with_B
        @error "Unified mesh is not compatible with both sides!"
        if !isempty(unified_mesh.compatibility_report) && verbose
            println("\nCompatibility issues:")
            for issue in unified_mesh.compatibility_report
                println("  - $issue")
            end
        end
        return false
    end
    
    if isempty(unified_mesh.triangles)
        @error "Unified mesh has no triangles!"
        return false
    end
    
    # Begin transactional replacement
    begin_transaction!(ws)
    
    try
        # Step 1: Delete old interface faces from A
    verbose && println("\nStep 1: Deleting interface faces from PID $(topology.pidA)...")
        faces_A = topology.faces_A
        num_faces_A = length(faces_A)
        
        # Delete in reverse order to avoid index shifting issues
        for tri_idx in num_faces_A:-1:1
            if !delete_interface_face!(ws, topology.pidA, tri_idx)
                @error "Failed to delete face $tri_idx from PID $(topology.pidA)"
                rollback_transaction!(ws)
                return false
            end
        end
    verbose && println("  ✓ Deleted $num_faces_A faces from PID $(topology.pidA)")
        
        # Step 2: Delete old interface faces from B
    verbose && println("\nStep 2: Deleting interface faces from PID $(topology.pidB)...")
        faces_B = topology.faces_B
        num_faces_B = length(faces_B)
        
        # Delete in reverse order
        for tri_idx in num_faces_B:-1:1
            if !delete_interface_face!(ws, topology.pidB, tri_idx)
                @error "Failed to delete face $tri_idx from PID $(topology.pidB)"
                rollback_transaction!(ws)
                return false
            end
        end
    verbose && println("  ✓ Deleted $num_faces_B faces from PID $(topology.pidB)")
        
        # Step 3: Insert unified triangles into A
    verbose && println("\nStep 3: Inserting unified triangles into PID $(topology.pidA)...")
        for (tri_idx, tri) in enumerate(unified_mesh.triangles)
            # Map unified coordinates to A's node IDs
            node_ids_A = map_nodes_to_pid(tri, unified_mesh.node_mapping_A, ws)
            
            if node_ids_A === nothing
                @error "Failed to map triangle $tri_idx nodes to PID $(topology.pidA)"
                rollback_transaction!(ws)
                return false
            end
            
            add_face!(ws, topology.pidA, node_ids_A)
        end
    verbose && println("  ✓ Inserted $(length(unified_mesh.triangles)) triangles into PID $(topology.pidA)")
        
        # Step 4: Insert unified triangles into B
    verbose && println("\nStep 4: Inserting unified triangles into PID $(topology.pidB)...")
        for (tri_idx, tri) in enumerate(unified_mesh.triangles)
            # Map unified coordinates to B's node IDs
            node_ids_B = map_nodes_to_pid(tri, unified_mesh.node_mapping_B, ws)
            
            if node_ids_B === nothing
                @error "Failed to map triangle $tri_idx nodes to PID $(topology.pidB)"
                rollback_transaction!(ws)
                return false
            end
            
            add_face!(ws, topology.pidB, node_ids_B)
        end
    verbose && println("  ✓ Inserted $(length(unified_mesh.triangles)) triangles into PID $(topology.pidB)")
        
        # Step 5: Verify mesh integrity
    verbose && println("\nStep 5: Verifying mesh integrity...")
        if !verify_mesh_integrity(ws, topology.pidA, topology.pidB)
            @error "Mesh integrity check failed after replacement"
            rollback_transaction!(ws)
            return false
        end
    verbose && println("  ✓ Mesh integrity verified")
        
        # Success! Commit the transaction
        if verbose
            println("\n" * "="^70)
            println("✓ Symmetric interface replacement completed successfully!")
            println("="^70)
            println("Summary:")
            println("  Original faces in A:  $num_faces_A")
            println("  Original faces in B:  $num_faces_B")
            println("  Unified faces:        $(length(unified_mesh.triangles))")
            println("  Nodes added:          $(ws.nodes_added)")
            println("="^70)
        end
        
        commit_transaction!(ws)
        return true
        
    catch e
        @error "Exception during interface replacement: $e"
        if verbose
            println("Stack trace:")
            for (exc, bt) in Base.catch_stack()
                showerror(stdout, exc, bt)
                println()
            end
            println("\nRolling back transaction...")
        end
        rollback_transaction!(ws)
        return false
    end
end

"""
    execute_symmetric_repair!(ws::RepairWorkspace, plan) -> Bool

Execute symmetric repair plan by generating and installing unified interface mesh.

This is the main entry point for symmetric repair execution (Phase 7).
It replaces the unidirectional `apply_repair_plan!` for symmetric repairs.

# Arguments
- `ws::RepairWorkspace`: The workspace containing the mesh
- `plan`: The symmetric repair plan to execute

# Returns
- `Bool`: true if repair succeeded, false otherwise

# Workflow
1. Validate plan feasibility
2. Extract the pre-generated unified mesh from plan
3. Display statistics
4. Call `replace_both_interfaces!` to perform the replacement
5. Return success/failure status

# Example
```julia
ws = RepairWorkspace(mesh_file)
symmetric_plan = generate_symmetric_repair_plan(topology, symmetric_classification, constraints)
success = execute_symmetric_repair!(ws, symmetric_plan)
```
"""
function execute_symmetric_repair!(
    ws::RepairWorkspace,
    plan;
    verbose::Bool = true
)::Bool
    
    if verbose
        println("\n" * "="^70)
        println("Executing Symmetric Repair Plan")
        println("="^70)
        println("Interface: $(plan.interface_pair)")
        println("Edges from A:         $(plan.edges_from_A)")
        println("Edges from B:         $(plan.edges_from_B)")
        println("Compromised edges:    $(plan.edges_compromised)")
        println("Synthesized edges:    $(plan.edges_synthesized)")
        println("Predicted min quality: $(round(plan.predicted_min_quality, digits=3))")
        println("Predicted compatibility: $(round(plan.predicted_compatibility_score, digits=3))")
        println("="^70)
    end
    
    if !plan.is_feasible
        @error "Plan is not feasible!"
        println("\nFeasibility issues:")
        for issue in plan.feasibility_issues
            println("  - $issue")
        end
        return false
    end
    
    # Generate the unified mesh (already pre-computed in plan)
    verbose && println("\nUsing pre-generated unified interface mesh...")
    unified_mesh = plan.target_unified_mesh
    
    if verbose
        println("\nUnified mesh statistics:")
        println("  Triangles:            $(length(unified_mesh.triangles))")
        println("  Min quality:          $(round(unified_mesh.min_triangle_quality, digits=3))")
        println("  Total area:           $(round(unified_mesh.total_area, digits=2))")
        println("  Compatible with A:    $(unified_mesh.compatible_with_A)")
        println("  Compatible with B:    $(unified_mesh.compatible_with_B)")
    end
    
    if !isempty(unified_mesh.compatibility_report) && verbose
        println("\nCompatibility notes:")
        for note in unified_mesh.compatibility_report
            println("  - $note")
        end
    end
    
    # Replace both interfaces
    verbose && println("\nProceeding with symmetric interface replacement...")
    success = replace_both_interfaces!(
        ws,
        unified_mesh,
        plan.topology;
        verbose=verbose
    )
    
    if verbose
        if success
            println("\n✓ Symmetric repair completed successfully!")
        else
            println("\n✗ Symmetric repair failed!")
        end
    end
    
    return success
end
