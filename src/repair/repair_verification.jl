#=
    repair_verification.jl

Phase 3: Surgical Mesh Repair Verification

Verifies mesh conformity after repairs and validates that adjacent interfaces remain intact.
=#

using Printf
using JSON

## Note: verify_interface_conformity is intentionally not defined here.
## Tests provide a stub implementation, and production integration will
## add a proper implementation within the package module to avoid method
## redefinition in the test environment.

"""
    compare_conformity(before::Dict, after::Dict) -> Dict

Compare conformity before and after repair.
Returns improvement statistics.
"""
function compare_conformity(before::Dict, after::Dict)
    improvement = Dict(
        "conformity_improvement" => after["overall_conformity"] - before["overall_conformity"],
        "mismatches_resolved" => before["total_mismatches"] - after["total_mismatches"],
        "shared_edges_added" => after["total_shared_edges"] - before["total_shared_edges"],
        "before" => before,
        "after" => after
    )
    
    return improvement
end

"""
    print_conformity_report(conformity::Dict, title::String="Conformity Report")

Print a formatted conformity report.
"""
function print_conformity_report(conformity::Dict, title::String="Conformity Report")
    println("\n" * "="^70)
    println(title)
    println("="^70)
    println("Interface: $(conformity["interface_pair"])")
    println()
    println("Topology:")
    println("  Shared nodes:  $(conformity["total_shared_nodes"])")
    println("  Edges in A:    $(conformity["total_edges_A"])")
    println("  Edges in B:    $(conformity["total_edges_B"])")
    println("  Shared edges:  $(conformity["total_shared_edges"])")
    println()
    println("Conformity:")
    println("  Ratio A:       $(round(conformity["conformity_ratio_A"] * 100, digits=2))%")
    println("  Ratio B:       $(round(conformity["conformity_ratio_B"] * 100, digits=2))%")
    println("  Overall:       $(round(conformity["overall_conformity"] * 100, digits=2))%")
    println()
    println("Mismatches:    $(conformity["total_mismatches"])")
    
    if conformity["total_mismatches"] > 0
        mismatch_types = conformity["mismatch_types"]
        println("  T-junctions:   $(mismatch_types["T_JUNCTION"])")
        println("  Diagonal:      $(mismatch_types["DIAGONAL"])")
        println("  Refinement:    $(mismatch_types["REFINEMENT"])")
        println("  Unknown:       $(mismatch_types["UNKNOWN"])")
    end
    
    println("="^70)
end

"""
    print_improvement_report(improvement::Dict)

Print a formatted before/after comparison report.
"""
function print_improvement_report(improvement::Dict)
    println("\n" * "="^70)
    println("Repair Improvement Report")
    println("="^70)
    
    before = improvement["before"]
    after = improvement["after"]
    
    println("\nBEFORE REPAIR:")
    println("  Overall conformity: $(round(before["overall_conformity"] * 100, digits=2))%")
    println("  Total mismatches:   $(before["total_mismatches"])")
    println("  Shared edges:       $(before["total_shared_edges"])")
    
    println("\nAFTER REPAIR:")
    println("  Overall conformity: $(round(after["overall_conformity"] * 100, digits=2))%")
    println("  Total mismatches:   $(after["total_mismatches"])")
    println("  Shared edges:       $(after["total_shared_edges"])")
    
    println("\nIMPROVEMENT:")
    conf_change = improvement["conformity_improvement"] * 100
    conf_symbol = conf_change >= 0 ? "+" : ""
    println("  Conformity change:  $(conf_symbol)$(round(conf_change, digits=2))%")
    
    mismatches_resolved = improvement["mismatches_resolved"]
    mismatch_symbol = mismatches_resolved >= 0 ? "-" : "+"
    println("  Mismatches resolved: $(mismatch_symbol)$(abs(mismatches_resolved))")
    
    edges_added = improvement["shared_edges_added"]
    edges_symbol = edges_added >= 0 ? "+" : ""
    println("  Shared edges added:  $(edges_symbol)$(edges_added)")
    
    if conf_change > 0
        println("\n✓ REPAIR SUCCESSFUL - Conformity improved!")
    elseif conf_change == 0
        println("\n⚠ REPAIR NEUTRAL - No conformity change")
    else
        println("\n✗ REPAIR DEGRADED - Conformity decreased!")
    end
    
    println("="^70)
end

"""
    verify_adjacent_interfaces(mesh, repaired_interface::Tuple{Int,Int}) -> Bool

Verify that adjacent interfaces to the repaired interface remain conforming.
Returns true if all adjacent interfaces are still valid, false otherwise.
"""
function verify_adjacent_interfaces(mesh, repaired_interface::Tuple{Int,Int})
    pid_a, pid_b = repaired_interface
    
    println("\n" * "="^70)
    println("Verifying Adjacent Interfaces")
    println("="^70)
    println("Repaired interface: $repaired_interface")
    println()
    
    # Get all PIDs
    all_pids = sort(collect(keys(mesh.all_pid_surfaces)))
    
    # Find adjacent interfaces (interfaces that share nodes with repaired PIDs)
    adjacent_interfaces = Set{Tuple{Int,Int}}()
    
    for pid_other in all_pids
        if pid_other != pid_a && pid_other != pid_b
            # Check if pid_other shares nodes with pid_a or pid_b
            for (node_id, coord) in mesh.nodes
                # Check if node is used by pid_other
                used_by_other = false
                for face in mesh.all_pid_surfaces[pid_other]
                    if node_id in face
                        used_by_other = true
                        break
                    end
                end
                
                if used_by_other
                    # Check if also used by pid_a
                    for face in mesh.all_pid_surfaces[pid_a]
                        if node_id in face
                            push!(adjacent_interfaces, minmax(pid_a, pid_other))
                            break
                        end
                    end
                    
                    # Check if also used by pid_b
                    for face in mesh.all_pid_surfaces[pid_b]
                        if node_id in face
                            push!(adjacent_interfaces, minmax(pid_b, pid_other))
                            break
                        end
                    end
                end
            end
        end
    end
    
    # Remove the repaired interface itself
    delete!(adjacent_interfaces, minmax(pid_a, pid_b))
    
    println("Found $(length(adjacent_interfaces)) adjacent interfaces:")
    for (i, interface) in enumerate(sort(collect(adjacent_interfaces)))
        println("  $i. PID $(interface[1]) ↔ PID $(interface[2])")
    end
    
    if isempty(adjacent_interfaces)
        println("\nNo adjacent interfaces to verify.")
        println("="^70)
        return true
    end
    
    println("\nChecking conformity of adjacent interfaces...")
    all_valid = true
    
    for interface in sort(collect(adjacent_interfaces))
        pid1, pid2 = interface
        conformity = verify_interface_conformity(mesh, pid1, pid2)
        conf_pct = round(conformity["overall_conformity"] * 100, digits=2)
        
        status = conformity["overall_conformity"] >= 0.90 ? "✓" : "✗"
        println("  $status PID $pid1 ↔ $pid2: $conf_pct% conformity")
        
        if conformity["overall_conformity"] < 0.90
            all_valid = false
        end
    end
    
    println()
    if all_valid
        println("✓ All adjacent interfaces remain valid")
    else
        println("✗ Some adjacent interfaces degraded!")
    end
    
    println("="^70)
    return all_valid
end

"""
    export_verification_report(improvement::Dict, output_file::String)

Export verification report to JSON file.
"""
function export_verification_report(improvement::Dict, output_file::String)
    open(output_file, "w") do f
        JSON.print(f, improvement, 2)
    end
    println("Verification report exported to: $output_file")
end
