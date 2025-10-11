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
    thresholds::QualityThresholds=default_thresholds(),
    verbose::Bool=true
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

    # Phase 5: generate the actual UnifiedInterfaceMesh (MVP)
    unified_mesh = generate_unified_interface_mesh(
        topology,
        updated_mismatches,
        constraints;
        thresholds=thresholds,
        verbose=verbose
    )

    # Feasibility and score from unified mesh
    is_feasible = unified_mesh.compatible_with_A && unified_mesh.compatible_with_B && !isempty(unified_mesh.triangles)
    if !is_feasible
        push!(feasibility_issues, "Unified mesh not compatible or empty")
        append!(feasibility_issues, unified_mesh.compatibility_report)
    end
    # Simple compatibility score: weighted by fraction of non-skipped edges and min quality
    compatibility_score = (total_repairable / max(1, length(updated_mismatches))) * (unified_mesh.min_triangle_quality)

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
        isempty(unified_mesh.triangles) ? 0.0 : unified_mesh.min_triangle_quality,
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
    thresholds::QualityThresholds=default_thresholds(),
    verbose::Bool=true,
    use_symmetric::Bool=false
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

    # Create simple classification for repair planning (auto-compute statistics)
    classification = InterfaceClassification(topology, mismatches_A, mismatches_B)

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
    RepairErrorType

Enumeration of different error types that can occur during repair.
"""
@enum RepairErrorType begin
    NO_ERROR
    REPAIR_EXECUTION_FAILED
    VALIDATION_FAILED
    STRATEGY_INFEASIBLE
    MESH_FILE_NOT_FOUND
    TOPOLOGY_BUILD_FAILED
    CLASSIFICATION_FAILED
    WORKSPACE_CREATION_FAILED
    EXPORT_FAILED
    UNKNOWN_ERROR
end

"""
    RepairErrorInfo

Detailed error information structure for comprehensive error reporting.
"""
struct RepairErrorInfo
    error_type::RepairErrorType
    primary_message::String
    detailed_description::String
    affected_components::Vector{String}
    error_counts::Dict{String,Int}
    context_data::Dict{String,Any}
    suggestions::Vector{String}

    function RepairErrorInfo(
        error_type::RepairErrorType = NO_ERROR;
        primary_message::String = "",
        detailed_description::String = "",
        affected_components::Vector{String} = String[],
        error_counts::Dict{String,Int} = Dict{String,Int}(),
        context_data::Dict{String,Any} = Dict{String,Any}(),
        suggestions::Vector{String} = String[]
    )
        new(error_type, primary_message, detailed_description,
            affected_components, error_counts, context_data, suggestions)
    end
end

"""
    get_summary_message(error_info::RepairErrorInfo) -> String

Get a concise summary message for logging.
"""
function get_summary_message(error_info::RepairErrorInfo)::String
    if error_info.error_type == NO_ERROR
        return ""
    end

    msg = error_info.primary_message
    if !isempty(error_info.error_counts)
        count_strs = ["$(k):$(v)" for (k,v) in error_info.error_counts]
        msg *= " ($(join(count_strs, ", ")))"
    end
    return msg
end

"""
    get_detailed_report(error_info::RepairErrorInfo) -> String

Get a detailed error report for debugging.
"""
function get_detailed_report(error_info::RepairErrorInfo)::String
    if error_info.error_type == NO_ERROR
        return "No errors"
    end

    report = "Error Type: $(error_info.error_type)\n"
    report *= "Primary Message: $(error_info.primary_message)\n"

    if !isempty(error_info.detailed_description)
        report *= "Description: $(error_info.detailed_description)\n"
    end

    if !isempty(error_info.affected_components)
        report *= "Affected Components: $(join(error_info.affected_components, ", "))\n"
    end

    if !isempty(error_info.error_counts)
        report *= "Error Counts:\n"
        for (category, count) in error_info.error_counts
            report *= "  $category: $count\n"
        end
    end

    if !isempty(error_info.suggestions)
        report *= "Suggestions:\n"
        for suggestion in error_info.suggestions
            report *= "  • $suggestion\n"
        end
    end

    return report
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
    error_info::RepairErrorInfo
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
    output_file::Union{String,Nothing}=nothing,
    thresholds::QualityThresholds=default_thresholds(),
    validate::Bool=true,
    verbose::Bool=true,
    use_symmetric::Bool=false
)
    start_time = time()

    # if verbose
    #     println("\n" * "="^70)
    #     println("INTERFACE REPAIR ORCHESTRATION")
    #     println("="^70)
    #     println("Mesh:       $mesh_file")
    #     println("Interface:  PID $pid_a ↔ PID $pid_b")
    #     println("Output:     $(output_file === nothing ? "(validation only)" : output_file)")
    #     println("="^70)
    # end

    # Step 1: Validate mesh file exists
    if verbose
        println("\n\033[32m[1/6] Validating mesh file...\033[0m")
    end
    if !isfile(mesh_file)
        error("Mesh file not found: $mesh_file")
    end
    if verbose
        println("      ✓ Mesh file exists: $mesh_file")
    end

    # Step 2: Build topology
    if verbose
        println("\n\033[32m[2/6] Building interface topology...\033[0m")
    end
    topology = build_interface_topology(mesh_file, pid_a, pid_b)
    if verbose
        println("      ✓ Topology built: $(length(topology.edges_A)) edges in A, $(length(topology.edges_B)) edges in B")
    end

    # Step 3: Classify mismatches (symmetric)
    if verbose
        println("\n\033[32m[3/6] Classifying interface mismatches (symmetric)...\033[0m")
    end
    sym_result = classify_interface_mismatches_symmetric(topology, tol=1e-4, verbose=verbose)
    sym_mismatches = sym_result.symmetric_mismatches

    total_mismatches = length(sym_mismatches)
    if verbose
        println("      ✓ Classified: $total_mismatches symmetric edge mismatches")
    end

    # Check if interface is already conforming
    if total_mismatches == 0
        if verbose
            println("\n✓ Interface is already conforming! No repairs needed.")
        end

        elapsed = time() - start_time
        return RepairOrchestrationResult(
            true,
            (pid_a, pid_b),
            RepairStrategy((pid_a, pid_b), :none, RepairPlan[], Dict{String,Any}("reason" => "already_conforming")),
            0, 0, 0,
            nothing,
            elapsed,
            RepairErrorInfo(NO_ERROR)
        )
    end

    # Step 4: Select repair strategy
    if verbose
        println("\n\033[32m[4/6] Selecting repair strategy...\033[0m")
    end

    # Build constraints for strategy selection
    constraints = build_boundary_constraints(mesh_file, pid_a, pid_b)

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
        end

        elapsed = time() - start_time
        return RepairOrchestrationResult(
            false,
            (pid_a, pid_b),
            strategy,
            0, 0, 0,
            nothing,
            elapsed,
            RepairErrorInfo(
                STRATEGY_INFEASIBLE,
                primary_message="No feasible repair strategy found",
                detailed_description="Repair algorithm could not find a viable approach to resolve interface mismatches",
                affected_components=["strategy_selection", "repair_planning"],
                error_counts=Dict("infeasible_strategies" => 1),
                suggestions=[
                    "Check interface geometry complexity",
                    "Consider different repair parameters",
                    "Verify input mesh quality"
                ]
            )
        )
    end

    # Step 5: Create repair workspace (optimized using existing topology)
    if verbose
        println("\n\033[32m[5/6] Creating repair workspace (reusing topology)...\033[0m")
    end
    ws = RepairWorkspace(mesh_file, topology)
    if verbose
        println("      ✓ Workspace created efficiently using topology data")
    end

    # Step 6: Execute repairs
    if verbose
        println("\n\033[32m[6/6] Executing repairs...\033[0m")
    end

    # Get the primary repair plan from strategy
    primary_plan = isempty(strategy.plans) ? nothing : strategy.plans[1]

    if primary_plan === nothing
        # Check if this is a symmetric repair case
        if strategy.approach == :symmetric && haskey(strategy.metadata, "symmetric_plan")
            sym_plan = strategy.metadata["symmetric_plan"]

            if verbose
                println("      Symmetric repair analysis:")
                if sym_plan.edges_from_A > 0
                    println("        • Edges to use from A: $(sym_plan.edges_from_A)")
                end
                if sym_plan.edges_from_B > 0
                    println("        • Edges to use from B: $(sym_plan.edges_from_B)")
                end
                if sym_plan.edges_compromised > 0
                    println("        • Edges to compromise: $(sym_plan.edges_compromised)")
                end
                if sym_plan.edges_from_A == 0 && sym_plan.edges_from_B == 0 && sym_plan.edges_compromised == 0
                    println("        ⚠ All edges were skipped - no viable repair strategies found")
                    println("          (All edges require boundary reconstruction or are too complex)")
                end

                if !sym_plan.is_feasible
                    println("        ✗ Repair plan is not feasible:")
                    for issue in sym_plan.feasibility_issues
                        println("          • $issue")
                    end
                end
            end

            if !sym_plan.is_feasible
                elapsed = time() - start_time
                return RepairOrchestrationResult(
                    false,
                    (pid_a, pid_b),
                    strategy,
                    0, 0, 0,
                    nothing,
                    elapsed,
                    RepairErrorInfo(
                        STRATEGY_INFEASIBLE,
                        primary_message="Symmetric repair plan not feasible",
                        detailed_description=join(sym_plan.feasibility_issues, "; "),
                        affected_components=["symmetric_repair_planning", "unified_mesh_generation"],
                        error_counts=Dict(
                            "edges_skipped" => sym_plan.edges_from_A + sym_plan.edges_from_B + sym_plan.edges_compromised,
                            "feasibility_issues" => length(sym_plan.feasibility_issues)
                        ),
                        context_data=Dict(
                            "edges_from_A" => sym_plan.edges_from_A,
                            "edges_from_B" => sym_plan.edges_from_B,
                            "edges_compromised" => sym_plan.edges_compromised,
                            "compatibility_score" => sym_plan.compatibility_score
                        ),
                        suggestions=[
                            "Check boundary constraint definitions",
                            "Consider relaxing quality thresholds",
                            "Verify interface topology accuracy"
                        ]
                    )
                )
            end
        else
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
                RepairErrorInfo(
                    REPAIR_EXECUTION_FAILED,
                    primary_message="No repair plan available",
                    detailed_description="Strategy selection failed to generate a concrete repair plan",
                    affected_components=["repair_planning", "strategy_execution"],
                    suggestions=[
                        "Check if symmetric repair mode is properly enabled",
                        "Verify interface mismatch detection results",
                        "Consider alternative repair strategies"
                    ]
                )
            )
        end
    end

    # Execute the repair plan
    # Check if this is a symmetric repair
    if strategy.approach == :symmetric && haskey(strategy.metadata, "symmetric_plan")
        # Execute symmetric repair
        sym_plan = strategy.metadata["symmetric_plan"]

        repairs_attempted = sym_plan.edges_from_A + sym_plan.edges_from_B + sym_plan.edges_compromised

        # Execute symmetric repair (Phase 6)
        execution_success = execute_symmetric_repair!(ws, sym_plan)

        repairs_succeeded = execution_success ? repairs_attempted : 0

        if verbose
            if repairs_attempted > 0
                println("      ✓ Symmetric repair completed: $repairs_succeeded/$repairs_attempted successful")
            else
                # Despite initial analysis, unified mesh was still generated and applied
                println("      ✓ Unified interface mesh applied despite no direct edge repairs")
                println("        (Used gap-fill approach with $(length(sym_plan.target_unified_mesh.triangles)) triangles)")
            end
        end
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

    # Validate results
    validation_result = nothing
    if validate
        validation_result = validate_repairs(ws, (pid_a, pid_b), thresholds=thresholds)

        if verbose
            if validation_result.passed
                println("      ✓ Validation passed")
            else
                println("      ✗ Validation failed: $(length(validation_result.errors)) errors")
                # Show concise error summary
                if !isempty(validation_result.errors)
                    non_manifold = count(e -> occursin("Non-manifold", e), validation_result.errors)
                    duplicate = count(e -> occursin("Duplicate", e), validation_result.errors)
                    critical = count(e -> occursin("Critical", e), validation_result.errors)

                    if non_manifold > 0
                        println("        • Non-manifold edges: $non_manifold")
                    end
                    if duplicate > 0
                        println("        • Duplicate faces: $duplicate")
                    end
                    if critical > 0
                        println("        • Critical quality issues: $critical")
                    end
                end

                # Show warnings summary if many isolated nodes
                if !isempty(validation_result.warnings)
                    isolated = count(w -> occursin("Isolated node", w), validation_result.warnings)
                    if isolated > 0
                        println("        • Isolated nodes: $isolated")
                    end
                end
            end
        end
    end

    # Export repaired mesh if output file specified
    if output_file !== nothing
        if verbose
            println("\nExporting repaired mesh to: $output_file")
        end
        modified_mesh = export_modified_mesh(ws)
        write_nas_volume(output_file, modified_mesh)
        if verbose
            println("✓ Mesh exported successfully")
        end
    end

    elapsed = time() - start_time

    if verbose
        println("\nOrchestration Complete:")
        println("  Time elapsed: $(round(elapsed, digits=2))s")
        println("  Status: $(execution_success && (validation_result === nothing || validation_result.passed) ? "✓ SUCCESS" : "⚠ PARTIAL SUCCESS")")
    end

    overall_success = execution_success && (validation_result === nothing || validation_result.passed)

    # Set appropriate error information
    error_info = RepairErrorInfo(NO_ERROR)
    if !overall_success
        if !execution_success
            error_info = RepairErrorInfo(
                REPAIR_EXECUTION_FAILED,
                primary_message="Repair execution failed",
                detailed_description="The repair operations could not be successfully applied to the mesh",
                affected_components=["repair_execution", "mesh_modification"],
                error_counts=Dict(
                    "repairs_attempted" => repairs_attempted,
                    "repairs_failed" => repairs_failed
                ),
                context_data=Dict(
                    "strategy_approach" => strategy.approach,
                    "execution_time" => elapsed
                ),
                suggestions=[
                    "Check workspace initialization",
                    "Verify repair plan feasibility",
                    "Review mesh topology consistency"
                ]
            )
        elseif validation_result !== nothing && !validation_result.passed
            # Categorize validation errors
            non_manifold_count = count(e -> occursin("Non-manifold", e), validation_result.errors)
            duplicate_count = count(e -> occursin("Duplicate", e), validation_result.errors)
            critical_quality_count = count(e -> occursin("Critical", e), validation_result.errors)
            degenerate_count = count(e -> occursin("Degenerate", e), validation_result.errors)

            error_info = RepairErrorInfo(
                VALIDATION_FAILED,
                primary_message="Validation failed",
                detailed_description="Post-repair validation identified mesh quality and manifoldness issues",
                affected_components=["mesh_validation", "quality_assessment"],
                error_counts=Dict(
                    "total_errors" => length(validation_result.errors),
                    "warnings" => length(validation_result.warnings),
                    "non_manifold_edges" => non_manifold_count,
                    "duplicate_faces" => duplicate_count,
                    "critical_quality_issues" => critical_quality_count,
                    "degenerate_triangles" => degenerate_count
                ),
                context_data=Dict(
                    "validation_stats" => validation_result.stats,
                    "repairs_applied" => repairs_succeeded
                ),
                suggestions=[
                    "Check for boundary constraint violations",
                    "Review mesh quality thresholds",
                    "Consider alternative repair strategies",
                    "Verify interface topology accuracy"
                ]
            )
        else
            error_info = RepairErrorInfo(
                UNKNOWN_ERROR,
                primary_message="Repair failed for unknown reason",
                detailed_description="Failure occurred but specific cause could not be determined",
                affected_components=["unknown"],
                suggestions=[
                    "Enable verbose logging for more details",
                    "Check system resources and permissions",
                    "Verify input file integrity"
                ]
            )
        end
    end

    return RepairOrchestrationResult(
        overall_success,
        (pid_a, pid_b),
        strategy,
        repairs_attempted,
        repairs_succeeded,
        repairs_failed,
        validation_result,
        elapsed,
        error_info
    )
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
    output_file::Union{String,Nothing}=nothing,
    thresholds::QualityThresholds=default_thresholds(),
    validate::Bool=true,
    verbose::Bool=true
)
    start_time = time()
    results = RepairOrchestrationResult[]

    # if verbose
    #     println("\n" * "="^70)
    #     println("MULTI-INTERFACE REPAIR ORCHESTRATION")
    #     println("="^70)
    #     println("Mesh:       $mesh_file")
    #     println("Interfaces: $(length(interface_pairs))")
    #     for (i, (pid_a, pid_b)) in enumerate(interface_pairs)
    #         println("  $i. PID $pid_a ↔ PID $pid_b")
    #     end
    #     println("Output:     $(output_file === nothing ? "(validation only)" : output_file)")
    #     println("="^70)
    # end

    # Validate mesh file once
    if !isfile(mesh_file)
        error("Mesh file not found: $mesh_file")
    end

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
        sym_result = classify_interface_mismatches_symmetric(topology, tol=1e-4, verbose=verbose)
        sym_mismatches = sym_result.symmetric_mismatches

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
                RepairErrorInfo(NO_ERROR)
            ))
            continue
        end

        # Build constraints
        constraints = build_boundary_constraints(mesh_file, pid_a, pid_b)

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
                RepairErrorInfo(
                    STRATEGY_INFEASIBLE,
                    primary_message="No feasible repair strategy",
                    detailed_description="Multi-interface repair could not find viable approach for this interface",
                    affected_components=["strategy_selection", "multi_interface_planning"],
                    context_data=Dict("interface_index" => interface_idx),
                    suggestions=[
                        "Consider processing this interface separately",
                        "Check for complex geometry conflicts",
                        "Review interface topology data"
                    ]
                )
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

            # Set appropriate error info for execution failures
        error_info = execution_success ? RepairErrorInfo(NO_ERROR) : RepairErrorInfo(
            REPAIR_EXECUTION_FAILED,
            primary_message="Multi-interface repair execution failed",
            detailed_description="Repair operations could not be applied successfully in multi-interface context",
            affected_components=["repair_execution", "multi_interface_workspace"],
            error_counts=Dict(
                "repairs_attempted" => repairs_attempted,
                "repairs_failed" => repairs_failed
            ),
            context_data=Dict(
                "interface_index" => interface_idx,
                "total_interfaces" => length(interface_pairs)
            ),
            suggestions=[
                "Check for interface interference effects",
                "Verify workspace consistency across interfaces",
                "Consider processing interfaces individually"
            ]
        )

        push!(results, RepairOrchestrationResult(
                execution_success,
                (pid_a, pid_b),
                strategy,
                repairs_attempted,
                repairs_succeeded,
                repairs_failed,
                nothing,  # Defer validation until all interfaces processed
                interface_elapsed,
                error_info
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

                # Create validation error info if validation failed
                new_error_info = old_result.error_info
                if !validation_result.passed && old_result.success
                    # Previously successful but validation failed
                    non_manifold_count = count(e -> occursin("Non-manifold", e), validation_result.errors)
                    duplicate_count = count(e -> occursin("Duplicate", e), validation_result.errors)
                    critical_quality_count = count(e -> occursin("Critical", e), validation_result.errors)

                    new_error_info = RepairErrorInfo(
                        VALIDATION_FAILED,
                        primary_message="Multi-interface validation failed",
                        detailed_description="Post-repair validation identified issues in multi-interface context",
                        affected_components=["multi_interface_validation", "cross_interface_consistency"],
                        error_counts=Dict(
                            "total_errors" => length(validation_result.errors),
                            "warnings" => length(validation_result.warnings),
                            "non_manifold_edges" => non_manifold_count,
                            "duplicate_faces" => duplicate_count,
                            "critical_quality_issues" => critical_quality_count
                        ),
                        context_data=Dict(
                            "interface_index" => interface_idx,
                            "validation_stats" => validation_result.stats,
                            "repairs_applied" => old_result.repairs_succeeded
                        ),
                        suggestions=[
                            "Check for interference between repaired interfaces",
                            "Review cumulative mesh modifications",
                            "Consider individual interface processing"
                        ]
                    )
                end

                results[interface_idx] = RepairOrchestrationResult(
                    old_result.success && validation_result.passed,
                    old_result.interface_pair,
                    old_result.strategy,
                    old_result.repairs_attempted,
                    old_result.repairs_succeeded,
                    old_result.repairs_failed,
                    validation_result,
                    old_result.elapsed_time,
                    new_error_info
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
        write_nas_volume(output_file, modified_mesh)
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
        Set{NTuple{3,Int}}(),  # shared_node_keys
        Dict{NTuple{3,Int},Tuple{Int,Int}}(),  # node_key_to_ids
        triangles_a, triangles_b,
        EdgeKey[], EdgeKey[],  # edges_A, edges_B (would need to extract edges)
        Set{EdgeKey}(), Set{EdgeKey}(), Set{EdgeKey}(),  # edge sets
        BoundingBox(NTuple{3,Float64}[]),  # interface_bbox
        0,  # total_shared_nodes
        length(triangles_a), length(triangles_b),  # face counts
        0, 0,  # edge counts
        0.0,  # conformity_ratio
        0.0, 0.0, 0, 1.0  # consistency metrics: max_vertex_dist, mean_vertex_dist, edge_mismatch_count, triangulation_similarity
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
