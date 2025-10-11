#!/usr/bin/env julia
"""
Simple Mesh Conformity Repair (using orchestration)

Repairs non-conforming interfaces in a Nastran mesh using the orchestration
pipeline in Nas2Step: topology → classification → strategy selection → execution → validation

Usage: julia repair_mesh_simple.jl input.nas output.nas [options]
"""

using Nas2Step
using Gmsh
using Printf

# Import error handling functions
import Nas2Step: get_summary_message, RepairErrorType

println("Mesh Conformity Repair (Orchestration)")

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
function main()
    if length(ARGS) < 2
        println("\nUsage: julia repair_mesh_simple.jl input.nas output.nas [options]")
        println("Example: julia repair_mesh_simple.jl examples/realistic/NC_Reduction_4.nas repaired.nas")
        println("\nOptions:")
        println("  --single-interface pid_a pid_b  Repair only specific interface")
        println("  --symmetric                     Use symmetric repair strategy")
        println("  --no-validation                Skip validation step")
        println("  --quiet                        Reduce verbosity")
        exit(1)
    end

    # Parse arguments
    INPUT_FILE = ARGS[1]
    OUTPUT_FILE = ARGS[2]

    # Parse optional arguments
    single_interface = nothing
    use_symmetric = false
    validate = true
    verbose = true
    i = 3
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--single-interface" && i + 2 <= length(ARGS)
            pid_a = parse(Int, ARGS[i+1])
            pid_b = parse(Int, ARGS[i+2])
            single_interface = (pid_a, pid_b)
            i += 3
        elseif arg == "--symmetric"
            use_symmetric = true
            i += 1
        elseif arg == "--no-validation"
            validate = false
            i += 1
        elseif arg == "--quiet"
            verbose = false
            i += 1
        else
            println("Unknown argument: $arg")
            exit(1)
        end
    end

    if !isfile(INPUT_FILE)
        error("Input file not found: $INPUT_FILE")
    end

    println("\nInput : $INPUT_FILE")
    println("Output: $OUTPUT_FILE")
    if single_interface !== nothing
        println("Interface: PID $(single_interface[1]) ↔ PID $(single_interface[2])")
    end

    # -----------------------------------------------------------------------------
    # Helper Functions
    # -----------------------------------------------------------------------------

    """Find all volume tag pairs that share a meaningful set of nodes (potential interfaces)."""
    function find_interface_pairs(nas_file::AbstractString; min_shared::Int=10)
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        try
            gmsh.open(nas_file)
            vols = [tag for (dim, tag) in gmsh.model.getEntities(3)]
            if length(vols) < 2
                return Tuple{Int,Int}[]
            end

            # Build a global node tag -> coordinate lookup once
            node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
            tag_to_coord = Dict{Int,NTuple{3,Float64}}()
            for (i, t) in enumerate(node_tags)
                base = (i - 1) * 3
                if base + 3 <= length(node_coords)
                    tag_to_coord[Int(t)] = (node_coords[base+1], node_coords[base+2], node_coords[base+3])
                end
            end

            # Nodes by volume
            vol_nodes = Dict{Int,Set{NTuple{3,Float64}}}()
            for v in vols
                etypes, _, enodes = gmsh.model.mesh.getElements(3, v)
                nodes_here = Set{NTuple{3,Float64}}()
                for (t, ev) in zip(etypes, enodes)
                    # 4 = tetra (Gmsh type)
                    if t != 4
                        continue
                    end
                    for i in 1:4:length(ev)
                        # all node ids of the tet
                        # We just record their coords as keys
                        for k in 0:3
                            nid = Int(ev[i+k])
                            coord = get(tag_to_coord, nid, nothing)
                            if coord === nothing
                                continue
                            end
                            push!(nodes_here, Nas2Step.coordinate_key_int(coord))
                        end
                    end
                end
                vol_nodes[v] = nodes_here
            end

            pairs = Tuple{Int,Int}[]
            for i in 1:length(vols)
                for j in (i+1):length(vols)
                    a, b = vols[i], vols[j]
                    if length(intersect(vol_nodes[a], vol_nodes[b])) >= min_shared
                        push!(pairs, (a, b))
                    end
                end
            end
            return pairs
        finally
            gmsh.finalize()
        end
    end

    # -----------------------------------------------------------------------------
    # Main Execution
    # -----------------------------------------------------------------------------

    if single_interface !== nothing
        # Single interface repair
        println("Processing single interface...")

        # Execute single interface repair
        result = Nas2Step.orchestrate_single_interface_repair(
            INPUT_FILE,
            single_interface[1],
            single_interface[2];
            output_file=OUTPUT_FILE,
            validate=validate,
            verbose=verbose,
            use_symmetric=use_symmetric
        )

        if result.success
            println("✓ Interface repair completed successfully")
            println("  Strategy: $(result.strategy.approach)")
            println("  Repairs: $(result.repairs_attempted) attempted, $(result.repairs_succeeded) succeeded")
            println("  Time: $(round(result.elapsed_time, digits=3))s")
            println("  Output: $OUTPUT_FILE")
        else
            # Use structured error information for better reporting
            error_summary = get_summary_message(result.error_info)
            println("✗ Interface repair failed: $error_summary")

            # Show additional details for specific error types
            if result.error_info.error_type == Nas2Step.VALIDATION_FAILED
                if haskey(result.error_info.error_counts, "non_manifold_edges") && result.error_info.error_counts["non_manifold_edges"] > 0
                    println("  • Non-manifold edges: $(result.error_info.error_counts["non_manifold_edges"])")
                end
                if haskey(result.error_info.error_counts, "duplicate_faces") && result.error_info.error_counts["duplicate_faces"] > 0
                    println("  • Duplicate faces: $(result.error_info.error_counts["duplicate_faces"])")
                end
                if haskey(result.error_info.error_counts, "critical_quality_issues") && result.error_info.error_counts["critical_quality_issues"] > 0
                    println("  • Critical quality issues: $(result.error_info.error_counts["critical_quality_issues"])")
                end
            elseif result.error_info.error_type == Nas2Step.STRATEGY_INFEASIBLE
                if !isempty(result.error_info.suggestions)
                    println("  Suggestions:")
                    for suggestion in result.error_info.suggestions[1:min(2, end)]  # Show first 2 suggestions
                        println("    • $suggestion")
                    end
                end
            end

            exit(1)
        end

    else
        # Multi-interface repair
        println("Analyzing interfaces...")
        pairs = find_interface_pairs(INPUT_FILE)
        if isempty(pairs)
            println("No interfaces detected. Nothing to do.")
            cp(INPUT_FILE, OUTPUT_FILE; force=true)
            exit(0)
        end

        for (i, (a, b)) in enumerate(pairs)
            println(@sprintf("  %2d) PID %d ↔ PID %d", i, a, b))
        end

        # Wrap multi-interface processing in a function to avoid scope issues
        function process_interfaces()
            successful_repairs = 0
            failed_repairs = 0
            session_start = time()

            for (idx, (pidA, pidB)) in enumerate(pairs)
                if verbose
                    println()
                    println("="^80)
                    println("\nInterface $(idx)/$(length(pairs)): PID $pidA ↔ PID $pidB")
                end

                # Execute repair for this interface
                result = Nas2Step.orchestrate_single_interface_repair(
                    INPUT_FILE,
                    pidA,
                    pidB;
                    output_file=OUTPUT_FILE,  # Export directly to final file
                    validate=validate,
                    verbose=verbose,
                    use_symmetric=use_symmetric
                )

                if result.success
                    successful_repairs += 1
                    if verbose
                        println("  ✓ Interface repaired successfully")
                    end
                else
                    failed_repairs += 1
                    if verbose
                        error_summary = get_summary_message(result.error_info)
                        println("  ✗ Interface repair failed: $error_summary")

                        # Show brief additional context based on error type
                        if result.error_info.error_type == Nas2Step.VALIDATION_FAILED
                            total_errors = get(result.error_info.error_counts, "total_errors", 0)
                            println("    ($total_errors validation errors)")
                        elseif result.error_info.error_type == Nas2Step.STRATEGY_INFEASIBLE
                            if haskey(result.error_info.context_data, "edges_skipped")
                                edges_skipped = result.error_info.context_data["edges_skipped"]
                                println("    ($edges_skipped edges skipped - no viable strategies)")
                            end
                        end
                    end
                end
            end

            session_duration = time() - session_start

            println("\nRepair Summary:")
            println("  Total interfaces: $(length(pairs))")
            println("  Successfully repaired: $successful_repairs")
            println("  Failed repairs: $failed_repairs")
            println("  Total time: $(round(session_duration, digits=2))s")
            println("  Output file: $OUTPUT_FILE")
        end

        process_interfaces()
    end

    println("\n✅ Repair script completed successfully!")
end

# Run the main function
main()