#!/usr/bin/env julia
# Debug script for analyzing diagonal mismatch failures
# This script runs the edge classification with detailed diagnostics enabled

using Pkg
Pkg.activate(".")

using Nas2Step
using Gmsh: gmsh

"""
Run classification with debug output for a specific interface.
"""
function debug_interface(nas_file::String, pidA::Int, pidB::Int; 
                        tol::Real=1e-4, debug_samples::Int=5)
    
    println("="^70)
    println("DEBUG: Analyzing Interface PID $pidA ↔ PID $pidB")
    println("="^70)
    println("Input file: $nas_file")
    println("Tolerance: $tol")
    println()
    
    # Build topology
    println("Building interface topology...")
    topology = build_interface_topology(nas_file, pidA, pidB, tol=tol)
    
    println("\nTopology summary:")
    println("  Faces A: $(length(topology.faces_A))")
    println("  Faces B: $(length(topology.faces_B))")
    println("  Edges only in A: $(length(topology.edges_only_in_A))")
    println("  Edges only in B: $(length(topology.edges_only_in_B))")
    println("  Conformity: $(round(topology.conformity_ratio * 100, digits=2))%")
    
    # Run classification with debug enabled
    println("\n" * "="^70)
    println("Running classification with DEBUG enabled...")
    println("="^70)
    
    classification = classify_interface_mismatches(
        topology, 
        tol=tol, 
        debug=true,
        debug_samples=debug_samples
    )
    
    println("\n" * "="^70)
    println("SUMMARY")
    println("="^70)
    println("Total mismatches: $(length(classification.mismatches_A) + length(classification.mismatches_B))")
    println("  T-junctions: $(classification.t_junction_count)")
    println("  Diagonal: $(classification.diagonal_count)")
    println("  Refinement: $(classification.refinement_count)")
    println("  Unknown: $(classification.unknown_count)")
    println("\nFeasibility:")
    println("  Feasible: $(classification.total_feasible)")
    println("  Infeasible: $(classification.total_infeasible)")
    
    # Count diagonal mismatches with empty quads
    diagonal_mismatches = filter(m -> m.mismatch_type == Nas2Step.DIAGONAL, 
                                 vcat(classification.mismatches_A, classification.mismatches_B))
    empty_quad_count = count(m -> isempty(m.quad_vertices), diagonal_mismatches)
    
    println("\nDiagonal mismatch details:")
    println("  Total diagonal mismatches: $(length(diagonal_mismatches))")
    println("  With empty quads (failed to find): $empty_quad_count")
    println("  With valid quads: $(length(diagonal_mismatches) - empty_quad_count)")
    
    return classification
end

# Main execution
if length(ARGS) < 1
    println("Usage: julia debug_diagonal_mismatches.jl <nas_file> [pidA] [pidB] [num_samples]")
    println("\nExample:")
    println("  julia debug_diagonal_mismatches.jl examples/realistic/NC_Reduction_4.nas 1 3 5")
    println("\nIf pidA and pidB are not specified, will analyze the first interface found.")
    exit(1)
end

nas_file = ARGS[1]
pidA = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
pidB = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 3
num_samples = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 5

if !isfile(nas_file)
    println("ERROR: File not found: $nas_file")
    exit(1)
end

try
    classification = debug_interface(nas_file, pidA, pidB, debug_samples=num_samples)
    println("\n✓ Debug analysis complete!")
catch e
    println("\n✗ Error during analysis:")
    println(e)
    rethrow()
end
