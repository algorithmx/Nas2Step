# Mesh quality verification functions for Nas2Step

using Statistics

include("./verify/utils.jl")

include("./verify/element_volumes.jl")

include("./verify/conformity.jl")

include("./verify/coordination.jl")

include("./verify/overlap.jl")

include("./verify/closure.jl")

include("./verify/export.jl")



"""
    test_node_swap_fix(filename; swap_pair=(1,2))

Test if swapping a pair of nodes in all elements would fix a convention mismatch.
This performs a virtual test without modifying the mesh file.

Returns a named tuple with:
- original_inverted: number of inverted elements before swap
- swapped_inverted: number of inverted elements after virtual swap
- improvement_ratio: fraction of elements that became non-inverted
- would_fix: boolean indicating if swap would fix >90% of inversions
- swap_pair: the node pair that was tested
"""
function test_node_swap_fix(filename; swap_pair::Tuple{Int,Int}=(1,2))
    # Check original
    original = check_element_volumes(filename)
    
    # Check with swapped nodes
    swapped = check_element_volumes(filename, swap_nodes=swap_pair)
    
    original_inv = original.inverted_count
    swapped_inv = swapped.inverted_count
    total = original.total_elements
    
    # Calculate improvement
    fixed_count = original_inv - swapped_inv
    improvement_ratio = total > 0 ? fixed_count / total : 0.0
    
    # Would this fix most issues?
    would_fix = (swapped_inv < 0.1 * total)  # Less than 10% still inverted
    
    return (
        original_inverted=original_inv,
        swapped_inverted=swapped_inv,
        fixed_count=fixed_count,
        improvement_ratio=improvement_ratio,
        would_fix=would_fix,
        swap_pair=swap_pair,
        original_ratio=original_inv/total,
        swapped_ratio=swapped_inv/total
    )
end


"""
    comprehensive_mesh_check(nas_file; output_json="mesh_quality_report.json", 
                            run_conversion=false, step_output=nothing, verbose=true)

Run all mesh quality checks and optionally test STEP conversion.
Exports unified JSON report including all quality metrics.

Arguments:
- nas_file: Input NAS mesh file
- output_json: Output JSON report file (default: "mesh_quality_report.json")
- run_conversion: Also run NAS to STEP conversion and include anomalies (default: false)
- step_output: STEP output path if run_conversion=true (default: auto-generate)
- verbose: Print detailed results (default: true)

Returns named tuple with verification results and JSON path.
"""
function comprehensive_mesh_check(nas_file::AbstractString; 
                                 output_json::AbstractString="mesh_quality_report.json",
                                 run_conversion::Bool=false,
                                 step_output::Union{Nothing,AbstractString}=nothing,
                                 verbose::Bool=true)
    
    if !isfile(nas_file)
        error("NAS file not found: $nas_file")
    end
    
    if verbose
        println("\n" * "="^70)
        println("COMPREHENSIVE MESH QUALITY CHECK")
        println("="^70)
        println("File: $nas_file")
        println("="^70)
    end
    
    # Run mesh verification
    verification = verify_nas_mesh(nas_file, verbose=verbose)
    
    anomaly_file = nothing
    
    # Optionally run conversion
    if run_conversion
        if verbose
            println("\n" * "="^70)
            println("TESTING NAS TO STEP CONVERSION")
            println("="^70)
        end
        
        step_path = step_output === nothing ? replace(nas_file, ".nas" => "_quality_test.step") : step_output
        anomaly_file = replace(step_path, ".step" => "_anomalies.json")
        
        try
            nas_to_step(nas_file, step_path=step_path, emit_anomaly_json=true, 
                       anomaly_json_path=anomaly_file)
            if verbose
                println("\n✓ Conversion successful: $step_path")
                if isfile(anomaly_file)
                    println("  Anomaly report: $anomaly_file")
                end
            end
        catch e
            if verbose
                println("\n⚠️  Conversion failed: $e")
            end
            anomaly_file = nothing
        end
    end
    
    # Export unified JSON report
    if verbose
        println("\n" * "="^70)
        println("EXPORTING UNIFIED REPORT")
        println("="^70)
    end
    
    export_mesh_quality_json(verification, output_json, include_anomalies_from=anomaly_file)
    
    if verbose
        println("\n" * "="^70)
        println("CHECK COMPLETE")
        println("="^70)
        println("Overall Status: ", 
                verification.overall_status == :ok ? "✓ PASSED" :
                verification.overall_status == :convention_mismatch ? "🔄 CONVENTION MISMATCH" :
                verification.overall_status == :warning ? "⚠️  WARNING" :
                "❌ FAILED")
        println("Full report: $output_json")
        println("="^70)
    end
    
    return (verification=verification, json_path=output_json, anomaly_path=anomaly_file)
end


"""
    verify_nas_mesh(filename; verbose=true)

Comprehensive verification of a NAS mesh file.

Checks for:
1. Inverted/negative volume elements (twists)
2. Element quality metrics
3. Surface closure
4. Vertex coordination numbers

Returns a named tuple with all check results and an overall status.
"""
function verify_nas_mesh(filename; verbose=true)
    if !isfile(filename)
        error("File $filename not found!")
    end
    
    if verbose
        println("="^70)
        println("Mesh Quality Verification")
        println("File: $filename")
        println("="^70)
        println()
    end
    
    # Check element volumes for inversions
    if verbose
        println("Checking for inverted elements...")
    end
    vol_check = check_element_volumes(filename)
    swap_test_result = nothing  # Will be populated if convention mismatch detected
    
    if verbose
        println("  Elements: $(vol_check.total_elements)")
        print_inverted_summary(vol_check)
        vol_check.volume_stats !== nothing && print_volume_stats(vol_check.volume_stats)
        
        if vol_check.status == :ok
            println("  ✓ OK")
        elseif vol_check.status == :warning
            println("  ⚠️  Degenerate elements detected")
        elseif vol_check.status == :convention_mismatch
            println("  🔄 CONVENTION MISMATCH ($(fmt_pct(vol_check.inverted_count, vol_check.total_elements))% inverted)")
            swap_test = test_node_swap_fix(filename, swap_pair=(1,2))
            swap_test_result = swap_test
            print_swap_test_result(swap_test)
            println("  Likely: RHR/LHR convention difference. Solution: swap node pairs or use abs(volume).")
        else
            println("  ⚠️  INVERTED ELEMENTS (mesh defects):")
            length(vol_check.inverted_elements) > 0 && print_inverted_details(vol_check.inverted_elements, 5)
        end
        println()
    end
    
    # Check element quality
    verbose && println("Checking element quality...")
    qual_check = check_element_quality(filename)
    
    if verbose
        poor_pct = fmt_pct(qual_check.poor_count, length(qual_check.qualities))
        println("  Quality: mean=$(round(qual_check.mean, digits=3)), " *
               "median=$(round(qual_check.median, digits=3)), " *
               "poor(<0.1)=$(qual_check.poor_count) ($(poor_pct)%)")
        
        status_msg = qual_check.status == :ok ? "✓ OK" :
                    qual_check.status == :warning ? "⚠️  Many poor elements" : "⚠️  Quality issues"
        println("  $status_msg")
        println()
    end
    
    # Check surface closure
    verbose && println("Checking surface closure...")
    closure_check = check_surface_closure(filename)
    
    if verbose
        msg = closure_check.is_closed ? "✓ Closed" : "⚠️  Open ($(closure_check.boundary_edge_count) boundary edges)"  
        println("  $msg")
        println()
    end
    
    # Check vertex coordination numbers
    verbose && println("Checking vertex coordination...")
    coord_check = check_vertex_coordination(filename)
    
    if verbose
        under_count = length(coord_check.undercoordinated)
        over_count = length(coord_check.overcoordinated)

        println("  Vertices: $(coord_check.total_vertices), mean=$(round(coord_check.mean_coordination, digits=1)), p95=$(round(coord_check.p95_coordination, digits=1)), p99=$(round(coord_check.p99_coordination, digits=1))")

        if coord_check.status == :ok
            println("  ✓ OK")
        else
            if under_count > 0
                println("  ⚠️  $(under_count) undercoordinated (<$(coord_check.min_coord))")
            end
            if over_count > 0
                println("  ⚠️  $(over_count) overcoordinated (>p$(Int(coord_check.overcoord_percentile))=$(coord_check.max_coord_threshold))")
            end
        end
        println()
    end
    
    # Determine overall status
    overall_status = :ok
    if vol_check.status == :error || qual_check.status == :error || closure_check.status == :error || coord_check.status == :error
        overall_status = :error
    elseif vol_check.status == :convention_mismatch
        overall_status = :convention_mismatch
    elseif qual_check.status == :warning || vol_check.status == :warning || coord_check.status == :warning
        overall_status = :warning
    end
    
    if verbose
        println("="^70)
        println("Overall Status: ", 
                overall_status == :ok ? "✓ PASSED" :
                overall_status == :convention_mismatch ? "🔄 CONVENTION MISMATCH" :
                overall_status == :warning ? "⚠️  WARNING" :
                "❌ FAILED")
        println("="^70)
    end
    
    return (
        volume_check=vol_check,
        quality_check=qual_check,
        closure_check=closure_check,
        coordination_check=coord_check,
        swap_test=swap_test_result,
        overall_status=overall_status
    )
end

