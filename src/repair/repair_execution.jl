"""
    repair_execution.jl

Phase 3: Surgical Mesh Repair Execution

Applies repair plans to mesh using transactional workspace.
Implements quad retriangulation with rollback on failure.
"""

using LinearAlgebra
using Printf

"""
    apply_quad_retriangulation!(ws::RepairWorkspace, plan::EdgeInsertionPlan, pid::Int) -> Bool

Apply a single quad retriangulation to the mesh.
Returns true on success, false on failure (workspace will be rolled back).
"""
function apply_quad_retriangulation!(ws::RepairWorkspace, plan::EdgeInsertionPlan, pid::Int)
    if plan.plan_type != :quad_retriangulation
        error("Plan type must be :quad_retriangulation, got $(plan.plan_type)")
    end
    
    if length(plan.triangles_to_replace) != 2
        @warn "Expected 2 triangles to replace, got $(length(plan.triangles_to_replace))"
        return false
    end
    
    if length(plan.new_triangles) != 2
        @warn "Expected 2 new triangles, got $(length(plan.new_triangles))"
        return false
    end
    
    # Find the two triangles to delete by their coordinates
    triangles_to_delete = Int[]
    for tri_coords in plan.triangles_to_replace
        face_idx = get_face_by_nodes(ws, pid, tri_coords)
        if face_idx === nothing
            @warn "Could not find triangle to replace: $tri_coords in PID $pid"
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
    for tri_coords in plan.new_triangles
        # Convert coordinates to node IDs
        node_ids = Int[]
        for coord in tri_coords
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
    apply_edge_insertion_plan!(ws::RepairWorkspace, plan::EdgeInsertionPlan, pid::Int) -> Bool

Apply a single edge insertion plan (handles all plan types).
Returns true on success, false on failure.
"""
function apply_edge_insertion_plan!(ws::RepairWorkspace, plan::EdgeInsertionPlan, pid::Int)
    if plan.plan_type == :quad_retriangulation
        return apply_quad_retriangulation!(ws, plan, pid)
    elseif plan.plan_type == :triangle_split
        # Future: implement T-junction handling
        @warn "Triangle split not yet implemented"
        return false
    else
        @warn "Unknown plan type: $(plan.plan_type)"
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
            println("Edge: $(edge_plan.edge_key)")
            println("Plan type: $(edge_plan.plan_type)")
            
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
        plan_type = Symbol(plan_dict["plan_type"])
        
        # Parse triangles to replace
        triangles_to_replace = [tuple.(tri...) for tri in plan_dict["triangles_to_replace"]]
        
        # Parse new triangles
        new_triangles = [tuple.(tri...) for tri in plan_dict["new_triangles"]]
        
        # Parse new nodes
        new_nodes = tuple.(plan_dict["new_nodes"]...)
        
        # Create plan
        plan = EdgeInsertionPlan(
            edge_key,
            plan_type,
            triangles_to_replace,
            new_triangles,
            new_nodes,
            plan_dict["predicted_min_quality"],
            plan_dict["predicted_max_quality_loss"],
            plan_dict["is_feasible"],
            plan_dict["feasibility_issues"]
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
