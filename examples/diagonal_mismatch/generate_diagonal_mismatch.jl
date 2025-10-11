#!/usr/bin/env julia

# Diagonal Mismatch Example Generator
#
# This script creates two touching boxes with vertex-compatible interfaces
# where triangulation differs in a diagonal manner, then exports to NAS format.

using Pkg
Pkg.activate("../../")

using Nas2Step
using .Nas2Step: DIAGONAL

# Include the diagonal mismatch example functions
include("example_diagonal_mismatch.jl")

# Simple mesh structure for the example
struct SimpleMesh
    vertices::Vector{NTuple{3,Float64}}
    triangles::Vector{NTuple{3,Int}}
    pid::Int
end

"""
    create_mesh(vertices, triangles, pid)

Create a simple mesh structure.
"""
function create_mesh(vertices::Vector{NTuple{3,Float64}},
                    triangles::Vector{NTuple{3,Int}},
                    pid::Int)
    return SimpleMesh(vertices, triangles, pid)
end

"""
    export_volumes_to_nas(all_vertices, box_a_tets, box_b_tets, filename)

Export two box volume meshes to Nastran NAS format with different PIDs.
"""
function export_volumes_to_nas(all_vertices::Vector{NTuple{3,Float64}},
                              box_a_tets::Vector{NTuple{4,Int}},
                              box_b_tets::Vector{NTuple{4,Int}},
                              filename::String)
    println("\n" * "="^60)
    println("EXPORTING TO NAS FILE")
    println("="^60)

    # Create volume regions for each box
    volume_regions = [
        VolumeRegion("BoxA", all_vertices, box_a_tets),
        VolumeRegion("BoxB", all_vertices, box_b_tets)
    ]

    # Write to NAS file
    write_nas_volume(filename, volume_regions)

    println("Exported to NAS file: $filename")
    println("  Box A: PID 1 ($(length(box_a_tets)) tetrahedra)")
    println("  Box B: PID 2 ($(length(box_b_tets)) tetrahedra)")

    return filename
end

"""
    analyze_diagonal_mismatch_from_nas(nas_filename)

Analyze the diagonal mismatch from a NAS file using the repair system.
"""
function analyze_diagonal_mismatch_from_nas(nas_filename::String)
    println("\n" * "="^60)
    println("ANALYZING DIAGONAL MISMATCH FROM NAS FILE")
    println("="^60)

    # Build interface topology from the NAS file
    # PID 1 and 2 are the two boxes we created
    topology = build_interface_topology(nas_filename, 1, 2, verbose=true)

    # Classify mismatches to detect the diagonal issue
    classification = classify_interface_mismatches(topology, verbose=true)

    println("\nMismatch classification:")
    println("  Total mismatches: $(length(classification.mismatches_A) + length(classification.mismatches_B))")
    println("  T-junctions: $(classification.t_junction_count)")
    println("  Diagonal mismatches: $(classification.diagonal_count)")
    println("  Other types: $(classification.refinement_count + classification.quad_mismatch_count + classification.boundary_edge_count)")

    # Show details of diagonal mismatches
    diagonal_mismatches = filter(m -> m.mismatch_type == DIAGONAL,
                                vcat(classification.mismatches_A, classification.mismatches_B))

    if !isempty(diagonal_mismatches)
        println("\nDiagonal mismatch details:")
        for (i, diag) in enumerate(diagonal_mismatches)
            println("  Diagonal $i:")
            println("    Edge: $(diag.edge_key)")
            println("    Present in: $(diag.present_in)")
            println("    Quad vertices: $(diag.quad_vertices)")
            println("    Triangles to replace: $(diag.triangles_to_replace)")
            println("    Complexity: $(round(diag.complexity_score, digits=3))")
            println("    Repair feasible: $(diag.repair_feasible)")
        end
    end

    return topology, classification
end

"""
    main()

Main function to run the complete diagonal mismatch example.
"""
function main()
    println("DIAGONAL MISMATCH EXAMPLE GENERATOR")
    println("="^60)
    println("This script generates two boxes with diagonal triangulation mismatch")
    println("and exports them to a Nastran NAS file for repair testing.")
    println()

    # Step 1: Create the diagonal mismatch example
    println("Step 1: Creating diagonal mismatch example...")
    all_vertices, box_a_tets, box_b_tets, shared_vertices = create_diagonal_mismatch_example()

    # Step 2: Export to NAS file
    println("\nStep 2: Exporting to NAS format...")
    nas_filename = "diagonal_mismatch_example.nas"
    export_volumes_to_nas(all_vertices, box_a_tets, box_b_tets, nas_filename)

    # Step 3: Analyze the mismatch from the NAS file
    println("\nStep 3: Analyzing the mismatch from NAS file...")
    topology, classification = analyze_diagonal_mismatch_from_nas(nas_filename)

    # Step 4: Visualize the issue
    println("\nStep 4: Visualizing the mismatch...")
    visualize_diagonal_mismatch(topology, classification)

    # Step 5: Show repair instructions
    println("\n" * "="^60)
    println("GENERATION COMPLETE")
    println("="^60)
    println("✓ Created diagonal mismatch example")
    println("✓ Analyzed interface topology")
    println("✓ Detected $(classification.diagonal_count) diagonal mismatch(es)")
    println("✓ Exported to NAS file: $nas_filename")
    println()
    println("Next steps:")
    println("1. Run repair script:")
    println("   cd ../../")
    println("   julia repair_mesh.jl examples/diagonal_mismatch/$nas_filename")
    println()
    println("2. The repair script will:")
    println("   - Detect the diagonal mismatch")
    println("   - Choose optimal triangulation")
    println("   - Generate repaired mesh")
    println("   - Verify the repair")

    return nas_filename, box_a_tets, box_b_tets, topology, classification
end

# Run the example if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    nas_filename, box_a_tets, box_b_tets, topology, classification = main()
end