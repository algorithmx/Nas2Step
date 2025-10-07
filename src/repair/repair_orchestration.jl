"""
    repair_orchestration.jl

Phase 5: Repair Orchestration

High-level functions that orchestrate the complete repair workflow:
topology → classification → strategy selection → execution → validation
"""

using Printf

"""
    RepairStrategy

Strategy for repairing an interface (placeholder structure).
"""
struct RepairStrategy
    interface_pair::Tuple{Int,Int}
    approach::Symbol  # :unidirectional, :bidirectional, :none, :infeasible, :error
    plans::Vector{RepairPlan}
    metadata::Dict{String,Any}
end

"""
    build_symmetric_repair_plan(
        topology::InterfaceTopology,
        sym_mismatches::Vector{SymmetricEdgeMismatch},
        constraints::BoundaryConstraints;
        thresholds::QualityThresholds = default_thresholds(),
        verbose::Bool = true
    ) -> SymmetricRepairPlan

Build a SymmetricRepairPlan from classified symmetric edge mismatches.

This function:
1. Assigns repair strategies to each symmetric mismatch using determine_repair_strategy()
2. Creates a placeholder UnifiedInterfaceMesh (Phase 5 would generate the actual mesh)
3. Computes statistics and feasibility
4. Returns a complete SymmetricRepairPlan

Note: This is a Phase 7 implementation that creates the plan structure.
Phase 5 (mesh generation) would fill in the actual unified mesh.
"""
function build_symmetric_repair_plan(
    topology::InterfaceTopology,
    sym_mismatches::Vector{SymmetricEdgeMismatch},
    constraints::BoundaryConstraints;
    thresholds::QualityThresholds = default_thresholds(),
    verbose::Bool = true
)::SymmetricRepairPlan
    
    if verbose
        println("      Building symmetric repair plan...")
    end
    
    # Assign strategies to all mismatches
    updated_mismatches = SymmetricEdgeMismatch[]
    edges_from_A = 0
    edges_from_B = 0
    edges_compromised = 0
    edges_skipped = 0
    min_quality = 1.0
    
    for sym in sym_mismatches
        # Determine strategy for this edge
        strategy, priority, reason = determine_repair_strategy(sym, constraints, thresholds=thresholds)
        
        # Create updated mismatch with strategy
        updated_sym = SymmetricEdgeMismatch(
            sym.edge_key,
            sym.classification_A_perspective,
            sym.classification_B_perspective,
            strategy,
            priority,
            reason
        )
        
        push!(updated_mismatches, updated_sym)
        
        # Update statistics
        if strategy == :use_A
            edges_from_A += 1
            if sym.classification_A_perspective !== nothing
                min_quality = min(min_quality, sym.classification_A_perspective.min_affected_triangle_quality)
            end
        elseif strategy == :use_B
            edges_from_B += 1
            if sym.classification_B_perspective !== nothing
                min_quality = min(min_quality, sym.classification_B_perspective.min_affected_triangle_quality)
            end
        elseif strategy == :compromise
            edges_compromised += 1
        elseif strategy == :skip
            edges_skipped += 1
        end
    end
    
    # Check feasibility
    total_repairable = edges_from_A + edges_from_B + edges_compromised
    is_feasible = total_repairable > 0
    feasibility_issues = String[]
    
    if edges_skipped == length(updated_mismatches)
        is_feasible = false
        push!(feasibility_issues, "All edges marked as skip - no repairs possible")
    end
    
    if !is_feasible
        push!(feasibility_issues, "Insufficient repairable edges ($total_repairable/$length(updated_mismatches))")
    end
    
    # NOTE: Phase 5 would generate the actual UnifiedInterfaceMesh here
    # For now, create a placeholder
    unified_mesh = UnifiedInterfaceMesh()
    
    # Compatibility score estimate (Phase 5 would compute actual score)
    compatibility_score = total_repairable / max(1, length(updated_mismatches))
    
    if verbose
        println("      Strategy statistics:")
        println("        Use A:       $edges_from_A")
        println("        Use B:       $edges_from_B")
        println("        Compromise:  $edges_compromised")
        println("        Skipped:     $edges_skipped")
        println("        Feasible:    $is_feasible")
    end
    
    return SymmetricRepairPlan(
        (topology.pidA, topology.pidB),
        updated_mismatches,
        unified_mesh,
        UnifiedMeshOperation[],  # Phase 5 would generate operations
        is_feasible,
        feasibility_issues,
        edges_from_A,
        edges_from_B,
        edges_compromised,
        0,  # edges_synthesized (Phase 5)
        min_quality,
        compatibility_score,
        topology,
        constraints
    )
end

"""
    select_repair_strategy_symmetric(topology, sym_mismatches, constraints; thresholds, verbose)

Select repair strategy for symmetric mismatches (simplified placeholder).
"""
function select_repair_strategy_symmetric(
    topology::InterfaceTopology,
    sym_mismatches::Vector{SymmetricEdgeMismatch},
    constraints::BoundaryConstraints;
    thresholds::QualityThresholds = default_thresholds(),
    verbose::Bool = true,
    use_symmetric::Bool = false
)
    if isempty(sym_mismatches)
        return RepairStrategy(
            (topology.pidA, topology.pidB),
            :none,
            RepairPlan[],
            Dict{String,Any}("reason" => "no_mismatches")
        )
    end
    
    # If symmetric mode requested, build SymmetricRepairPlan
    if use_symmetric
        if verbose
            println("      Using symmetric repair approach...")
        end
        
        sym_plan = build_symmetric_repair_plan(
            topology,
            sym_mismatches,
            constraints;
            thresholds=thresholds,
            verbose=verbose
        )
        
        # Return strategy wrapping the symmetric plan
        return RepairStrategy(
            (topology.pidA, topology.pidB),
            sym_plan.is_feasible ? :symmetric : :infeasible,
            RepairPlan[],  # Empty - using SymmetricRepairPlan instead
            Dict{String,Any}(
                "mode" => "symmetric",
                "symmetric_plan" => sym_plan,
                "edges_from_A" => sym_plan.edges_from_A,
                "edges_from_B" => sym_plan.edges_from_B,
                "edges_compromised" => sym_plan.edges_compromised
            )
        )
    end
    
    # Otherwise, use traditional bidirectional planning
    if verbose
        println("      Using traditional bidirectional repair approach...")
    end
    
    # Build classification from symmetric mismatches
    mismatches_A = [m.classification_A_perspective for m in sym_mismatches if m.classification_A_perspective !== nothing]
    mismatches_B = [m.classification_B_perspective for m in sym_mismatches if m.classification_B_perspective !== nothing]
    
    # Create simple classification for repair planning
    classification = InterfaceClassification(
        topology,
        mismatches_A,
        mismatches_B,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        length(mismatches_A) + length(mismatches_B),
        0,
        1.0
    )
    
    # Generate bidirectional plan
    plan_result = generate_repair_plan_bidirectional(
        topology,
        classification,
        constraints;
        thresholds=thresholds,
        verbose=verbose
    )
    
    plan, tried_both, direction_info = plan_result
    
    return RepairStrategy(
        (topology.pidA, topology.pidB),
        plan.is_feasible ? :bidirectional : :infeasible,
        [plan],
        direction_info
    )
end

"""
    RepairOrchestrationResult

Result of complete repair orchestration.
"""
struct RepairOrchestrationResult
    success::Bool
    interface_pair::Tuple{Int,Int}
    strategy::RepairStrategy
    repairs_attempted::Int
    repairs_succeeded::Int
    repairs_failed::Int
    validation::Union{ValidationResult,Nothing}
    elapsed_time::Float64
    error_message::Union{String,Nothing}
end

"""
    orchestrate_single_interface_repair(
        mesh_file::String,
        pid_a::Int,
        pid_b::Int;
        output_file::Union{String,Nothing} = nothing,
        thresholds::QualityThresholds = default_thresholds(),
        validate::Bool = true,
        verbose::Bool = true
    ) -> RepairOrchestrationResult

Orchestrate repair for a single interface between two PIDs.

Complete workflow:
1. Load mesh and build topology
2. Classify interface mismatches (symmetric)
3. Select repair strategy
4. Create repair workspace
5. Execute repairs
6. Validate results
7. Export repaired mesh (optional)

Returns comprehensive result with statistics and validation.
"""
function orchestrate_single_interface_repair(
    mesh_file::String,
    pid_a::Int,
    pid_b::Int;
    output_file::Union{String,Nothing} = nothing,
    thresholds::QualityThresholds = default_thresholds(),
    validate::Bool = true,
    verbose::Bool = true,
    use_symmetric::Bool = false
)
    start_time = time()
    
    if verbose
        println("\n" * "="^70)
        println("INTERFACE REPAIR ORCHESTRATION")
        println("="^70)
        println("Mesh:       $mesh_file")
        println("Interface:  PID $pid_a ↔ PID $pid_b")
        println("Output:     $(output_file === nothing ? "(validation only)" : output_file)")
        println("="^70)
    end
    
    try
        # Step 1: Load mesh
        if verbose
            println("\n[1/7] Loading mesh...")
        end
        mesh = read_nastran(mesh_file)
        if verbose
            println("      ✓ Loaded: $(length(mesh.nodes)) nodes, $(sum(length(f) for f in values(mesh.all_pid_surfaces))) faces")
        end
        
        # Step 2: Build topology
        if verbose
            println("\n[2/7] Building interface topology...")
        end
        topology = build_interface_topology(mesh, pid_a, pid_b)
        if verbose
            println("      ✓ Topology built: $(length(topology.edges_A)) edges in A, $(length(topology.edges_B)) edges in B")
        end
        
        # Step 3: Classify mismatches (symmetric)
        if verbose
            println("\n[3/7] Classifying interface mismatches (symmetric)...")
        end
        sym_mismatches = classify_interface_mismatches_symmetric(topology, tol=1e-4, verbose=verbose)
        
        total_mismatches = length(sym_mismatches)
        if verbose
            println("      ✓ Classified: $total_mismatches symmetric edge mismatches")
        end
        
        # Check if interface is already conforming
        if total_mismatches == 0
            if verbose
                println("\n✓ Interface is already conforming! No repairs needed.")
                println("="^70)
            end
            
            elapsed = time() - start_time
            return RepairOrchestrationResult(
                true,
                (pid_a, pid_b),
                RepairStrategy((pid_a, pid_b), :none, RepairPlan[], Dict{String,Any}("reason" => "already_conforming")),
                0, 0, 0,
                nothing,
                elapsed,
                nothing
            )
        end
        
        # Step 4: Select repair strategy
        if verbose
            println("\n[4/7] Selecting repair strategy...")
        end
        
        # Build constraints for strategy selection
        constraints = build_boundary_constraints(mesh, topology, pid_a, pid_b)
        
        # Use symmetric strategy selection
        strategy = select_repair_strategy_symmetric(
            topology,
            sym_mismatches,
            constraints;
            thresholds=thresholds,
            verbose=verbose,
            use_symmetric=use_symmetric
        )
        
        if verbose
            println("      ✓ Strategy selected: $(strategy.approach)")
        end
        
        # Check if strategy is feasible
        if strategy.approach == :infeasible
            if verbose
                println("\n✗ No feasible repair strategy found!")
                println("  Reasons:")
                for (i, reason) in enumerate(get(strategy.metadata, "reasons", String[]))
                    println("    $i. $reason")
                end
                println("="^70)
            end
            
            elapsed = time() - start_time
            return RepairOrchestrationResult(
                false,
                (pid_a, pid_b),
                strategy,
                0, 0, 0,
                nothing,
                elapsed,
                "No feasible repair strategy found"
            )
        end
        
        # Step 5: Create repair workspace
        if verbose
            println("\n[5/7] Creating repair workspace...")
        end
        ws = RepairWorkspace(mesh_file)
        if verbose
            println("      ✓ Workspace created")
        end
        
        # Step 6: Execute repairs
        if verbose
            println("\n[6/7] Executing repairs...")
        end
        
        # Get the primary repair plan from strategy
        primary_plan = isempty(strategy.plans) ? nothing : strategy.plans[1]
        
        if primary_plan === nothing
            if verbose
                println("      ✗ No repair plan available")
            end
            
            elapsed = time() - start_time
            return RepairOrchestrationResult(
                false,
                (pid_a, pid_b),
                strategy,
                0, 0, 0,
                nothing,
                elapsed,
                "No repair plan available"
            )
        end
        
        # Execute the repair plan
        # Check if this is a symmetric repair
        if strategy.approach == :symmetric && haskey(strategy.metadata, "symmetric_plan")
            # Execute symmetric repair
            sym_plan = strategy.metadata["symmetric_plan"]
            
            if verbose
                println("      Executing symmetric repair (Phase 6)...")
                println("      NOTE: Phase 5 (mesh generation) not yet implemented")
                println("            This will use placeholder unified mesh")
            end
            
            # Execute symmetric repair (Phase 6)
            execution_success = execute_symmetric_repair!(ws, sym_plan)
            
            repairs_attempted = sym_plan.edges_from_A + sym_plan.edges_from_B + sym_plan.edges_compromised
            repairs_succeeded = execution_success ? repairs_attempted : 0
            repairs_failed = repairs_attempted - repairs_succeeded
        else
            # Execute traditional unidirectional repair
            execution_success = apply_repair_plan!(ws, primary_plan)
            
            repairs_attempted = length(filter(p -> p.is_feasible, primary_plan.insertion_sequence))
            repairs_succeeded = ws.faces_added + ws.faces_deleted  # Approximate
            repairs_failed = repairs_attempted - repairs_succeeded
        end
        
        if verbose
            if execution_success
                println("      ✓ Repairs executed successfully")
            else
                println("      ⚠ Some repairs failed")
            end
            println("      Statistics:")
            println("        - Attempted: $repairs_attempted")
            println("        - Succeeded: $repairs_succeeded")
            println("        - Failed:    $repairs_failed")
        end
        
        # Step 7: Validate results
        validation_result = nothing
        if validate
            if verbose
                println("\n[7/7] Validating repairs...")
            end
            validation_result = validate_repairs(ws, (pid_a, pid_b), thresholds=thresholds)
            
            if verbose
                if validation_result.passed
                    println("      ✓ Validation passed")
                else
                    println("      ✗ Validation failed: $(length(validation_result.errors)) errors")
                end
            end
        else
            if verbose
                println("\n[7/7] Skipping validation (disabled)")
            end
        end
        
        # Export repaired mesh if output file specified
        if output_file !== nothing
            if verbose
                println("\nExporting repaired mesh to: $output_file")
            end
            modified_mesh = export_modified_mesh(ws)
            write_nastran(output_file, modified_mesh)
            if verbose
                println("✓ Mesh exported successfully")
            end
        end
        
        elapsed = time() - start_time
        
        if verbose
            println("\n" * "="^70)
            println("ORCHESTRATION COMPLETE")
            println("="^70)
            println("Time elapsed: $(round(elapsed, digits=2))s")
            println("Status: $(execution_success && (validation_result === nothing || validation_result.passed) ? "✓ SUCCESS" : "⚠ PARTIAL SUCCESS")")
            println("="^70)
        end
        
        overall_success = execution_success && (validation_result === nothing || validation_result.passed)
        
        return RepairOrchestrationResult(
            overall_success,
            (pid_a, pid_b),
            strategy,
            repairs_attempted,
            repairs_succeeded,
            repairs_failed,
            validation_result,
            elapsed,
            nothing
        )
        
    catch e
        elapsed = time() - start_time
        
        if verbose
            println("\n" * "="^70)
            println("✗ ORCHESTRATION FAILED")
            println("="^70)
            println("Error: $e")
            println("Time elapsed: $(round(elapsed, digits=2))s")
            println("="^70)
        end
        
        # Return failure result
        return RepairOrchestrationResult(
            false,
            (pid_a, pid_b),
            RepairStrategy((pid_a, pid_b), :error, RepairPlan[], Dict{String,Any}("error" => string(e))),
            0, 0, 0,
            nothing,
            elapsed,
            string(e)
        )
    end
end

"""
    orchestrate_multi_interface_repair(
        mesh_file::String,
        interface_pairs::Vector{Tuple{Int,Int}};
        output_file::Union{String,Nothing} = nothing,
        thresholds::QualityThresholds = default_thresholds(),
        validate::Bool = true,
        verbose::Bool = true
    ) -> Vector{RepairOrchestrationResult}

Orchestrate repairs for multiple interfaces in a single mesh.

Repairs are applied sequentially to a single workspace, allowing
cumulative modifications across interfaces.
"""
function orchestrate_multi_interface_repair(
    mesh_file::String,
    interface_pairs::Vector{Tuple{Int,Int}};
    output_file::Union{String,Nothing} = nothing,
    thresholds::QualityThresholds = default_thresholds(),
    validate::Bool = true,
    verbose::Bool = true
)
    start_time = time()
    results = RepairOrchestrationResult[]
    
    if verbose
        println("\n" * "="^70)
        println("MULTI-INTERFACE REPAIR ORCHESTRATION")
        println("="^70)
        println("Mesh:       $mesh_file")
        println("Interfaces: $(length(interface_pairs))")
        for (i, (pid_a, pid_b)) in enumerate(interface_pairs)
            println("  $i. PID $pid_a ↔ PID $pid_b")
        end
        println("Output:     $(output_file === nothing ? "(validation only)" : output_file)")
        println("="^70)
    end
    
    try
        # Load mesh once
        mesh = read_nastran(mesh_file)
        
        # Create workspace once
        ws = RepairWorkspace(mesh_file)
        
        # Process each interface
        for (interface_idx, (pid_a, pid_b)) in enumerate(interface_pairs)
            if verbose
                println("\n" * "#"^70)
                println("INTERFACE $interface_idx/$(length(interface_pairs)): PID $pid_a ↔ PID $pid_b")
                println("#"^70)
            end
            
            interface_start = time()
            
            # Build topology for this interface
            topology = build_interface_topology_from_workspace(ws, pid_a, pid_b)
            
            # Classify mismatches
            sym_mismatches = classify_interface_mismatches_symmetric(topology, tol=1e-4, verbose=verbose)
            
            if isempty(sym_mismatches)
                if verbose
                    println("✓ Interface already conforming, skipping.")
                end
                interface_elapsed = time() - interface_start
                
                push!(results, RepairOrchestrationResult(
                    true,
                    (pid_a, pid_b),
                    RepairStrategy((pid_a, pid_b), :none, RepairPlan[], Dict{String,Any}("reason" => "already_conforming")),
                    0, 0, 0,
                    nothing,
                    interface_elapsed,
                    nothing
                ))
                continue
            end
            
            # Build constraints
            constraints = build_boundary_constraints_from_workspace(ws, topology, pid_a, pid_b)
            
            # Select strategy
            strategy = select_repair_strategy_symmetric(
                topology,
                sym_mismatches,
                constraints;
                thresholds=thresholds,
                verbose=verbose
            )
            
            if strategy.approach == :infeasible
                if verbose
                    println("✗ No feasible strategy for this interface, skipping.")
                end
                interface_elapsed = time() - interface_start
                
                push!(results, RepairOrchestrationResult(
                    false,
                    (pid_a, pid_b),
                    strategy,
                    0, 0, 0,
                    nothing,
                    interface_elapsed,
                    "No feasible repair strategy"
                ))
                continue
            end
            
            # Execute repairs
            primary_plan = isempty(strategy.plans) ? nothing : strategy.plans[1]
            
            if primary_plan !== nothing
                execution_success = apply_repair_plan!(ws, primary_plan)
                
                repairs_attempted = length(filter(p -> p.is_feasible, primary_plan.insertion_sequence))
                repairs_succeeded = ws.faces_added + ws.faces_deleted
                repairs_failed = repairs_attempted - repairs_succeeded
                
                interface_elapsed = time() - interface_start
                
                push!(results, RepairOrchestrationResult(
                    execution_success,
                    (pid_a, pid_b),
                    strategy,
                    repairs_attempted,
                    repairs_succeeded,
                    repairs_failed,
                    nothing,  # Defer validation until all interfaces processed
                    interface_elapsed,
                    nothing
                ))
            end
        end
        
        # Final validation across all interfaces
        if validate && !isempty(results)
            if verbose
                println("\n" * "="^70)
                println("FINAL VALIDATION (all interfaces)")
                println("="^70)
            end
            
            # Run validation for each interface
            for (interface_idx, (pid_a, pid_b)) in enumerate(interface_pairs)
                validation_result = validate_repairs(ws, (pid_a, pid_b), thresholds=thresholds)
                
                # Update result with validation
                if interface_idx <= length(results)
                    old_result = results[interface_idx]
                    results[interface_idx] = RepairOrchestrationResult(
                        old_result.success && validation_result.passed,
                        old_result.interface_pair,
                        old_result.strategy,
                        old_result.repairs_attempted,
                        old_result.repairs_succeeded,
                        old_result.repairs_failed,
                        validation_result,
                        old_result.elapsed_time,
                        old_result.error_message
                    )
                end
            end
        end
        
        # Export final mesh
        if output_file !== nothing
            if verbose
                println("\nExporting final repaired mesh to: $output_file")
            end
            modified_mesh = export_modified_mesh(ws)
            write_nastran(output_file, modified_mesh)
            if verbose
                println("✓ Mesh exported successfully")
            end
        end
        
        elapsed = time() - start_time
        
        if verbose
            println("\n" * "="^70)
            println("MULTI-INTERFACE ORCHESTRATION COMPLETE")
            println("="^70)
            println("Total time:  $(round(elapsed, digits=2))s")
            println("Interfaces:  $(length(interface_pairs))")
            successful = count(r -> r.success, results)
            println("Successful:  $successful/$(length(results))")
            println("="^70)
        end
        
        return results
        
    catch e
        elapsed = time() - start_time
        
        if verbose
            println("\n" * "="^70)
            println("✗ MULTI-INTERFACE ORCHESTRATION FAILED")
            println("="^70)
            println("Error: $e")
            println("Time elapsed: $(round(elapsed, digits=2))s")
            println("="^70)
        end
        
        rethrow(e)
    end
end

"""
Helper: Build topology from workspace (after modifications).
"""
function build_interface_topology_from_workspace(ws::RepairWorkspace, pid_a::Int, pid_b::Int)
    # Extract faces from workspace
    faces_a = get(ws.working_faces, pid_a, Vector{Vector{Int}}[])
    faces_b = get(ws.working_faces, pid_b, Vector{Vector{Int}}[])
    
    # Convert to Triangle objects
    triangles_a = Triangle[]
    for face_nodes in faces_a
        if length(face_nodes) == 3
            coords = [ws.working_nodes[nid] for nid in face_nodes]
            tri = Triangle(coords[1], coords[2], coords[3])
            push!(triangles_a, tri)
        end
    end
    
    triangles_b = Triangle[]
    for face_nodes in faces_b
        if length(face_nodes) == 3
            coords = [ws.working_nodes[nid] for nid in face_nodes]
            tri = Triangle(coords[1], coords[2], coords[3])
            push!(triangles_b, tri)
        end
    end
    
    # Build topology using existing function (would need to be adapted)
    # For now, return a minimal topology
    # This would need proper implementation based on the actual topology builder
    
    return InterfaceTopology(
        pid_a, pid_b,
        triangles_a, triangles_b,
        EdgeKey[], EdgeKey[],  # Would need to extract edges
        Set{NTuple{3,Float64}}(),
        0, 0
    )
end

"""
Helper: Build constraints from workspace.
"""
function build_boundary_constraints_from_workspace(ws::RepairWorkspace, topology::InterfaceTopology, pid_a::Int, pid_b::Int)
    # Simplified constraint builder
    # Would need proper implementation
    
    return BoundaryConstraints(
        Set{EdgeKey}(),
        Set{EdgeKey}(),
        Set{NTuple{3,Float64}}(),
        Set{NTuple{3,Float64}}(),
        topology
    )
end
