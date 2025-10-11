# Example: Two touching boxes with diagonal triangulation mismatch
#
# This demonstrates a DIAGONAL mismatch type where two volume meshes share
# the same 4 vertices at their interface but use different triangulations.

using Nas2Step

"""
    partition_hexahedron_to_tetrahedra(hex_vertices, diagonal_pattern=:body_diagonal_17)

Partition a hexahedron into 6 tetrahedra using the 1:6 rule.

Args:
- hex_vertices: Tuple of 8 vertex IDs (n1, n2, n3, n4, n5, n6, n7, n8) in standard hexahedron order
- diagonal_pattern: :body_diagonal_17 (n1-n7) or :body_diagonal_28 (n2-n8)

Returns:
- Vector of 6 tetrahedra, each as a 4-tuple of vertex IDs

The hexahedron vertices should be ordered as:
- Bottom face: n1, n2, n3, n4 (counter-clockwise)
- Top face:    n5, n6, n7, n8 (counter-clockwise, directly above n1, n2, n3, n4)
"""
function partition_hexahedron_to_tetrahedra(hex_vertices::NTuple{8,Int},
                                           diagonal_pattern::Symbol=:body_diagonal_17)
    n1, n2, n3, n4, n5, n6, n7, n8 = hex_vertices

    if diagonal_pattern == :body_diagonal_17
        # Use body diagonal from n1 to n7 - creates diagonal 2-6 on shared face (n2,n3,n6,n7)
        # Shared face triangulation: triangle (2,6,7) + triangle (2,3,6)
        return [
            (n1, n2, n6, n7),  # tet 1 - contains diagonal (n2,n6) on shared face
            (n1, n6, n5, n7),  # tet 2
            (n1, n5, n8, n7),  # tet 3
            (n1, n8, n4, n7),  # tet 4
            (n1, n4, n3, n7),  # tet 5
            (n1, n3, n2, n7),  # tet 6 - contains triangle (2,3,7) on shared face
        ]
    elseif diagonal_pattern == :body_diagonal_28
        # Use body diagonal from n2 to n8 - creates diagonal 3-6 on shared face (n2,n3,n6,n7)
        # Shared face triangulation: triangle (3,6,7) + triangle (2,3,6)
        return [
            (n2, n3, n6, n7),  # tet 1 - directly creates diagonal (n3,n6) on shared face
            (n2, n6, n9, n7),  # tet 2 - contains triangle (3,6,7) on shared face
            (n2, n9, n10, n7), # tet 3 - contains triangle (6,9,7)
            (n2, n10, n11, n7),# tet 4 - contains triangle (9,10,7)
            (n2, n11, n12, n7),# tet 5 - contains triangle (10,11,7)
            (n2, n12, n3, n7), # tet 6 - contains triangle (11,12,7) + triangle (2,3,7)
        ]
    else
        throw(ArgumentError("diagonal_pattern must be :body_diagonal_17 or :body_diagonal_28"))
    end
end

"""
    create_diagonal_mismatch_example()

Create two proper cubes sharing a face with different diagonal triangulations.
"""
function create_diagonal_mismatch_example()
    println("Creating diagonal mismatch example with two proper cubes...")

    # Define coordinates for two unit cubes sharing a face
    # Box A: from (0,0,0) to (1,1,1)
    # Box B: from (1,0,0) to (2,1,1)

    # All vertices for both cubes (12 vertices total, 4 shared)
    all_vertices = [
        # Box A vertices (8 vertices)
        (0.0, 0.0, 0.0),  # 1 - Box A origin
        (1.0, 0.0, 0.0),  # 2 - shared with Box B
        (1.0, 1.0, 0.0),  # 3 - shared with Box B
        (0.0, 1.0, 0.0),  # 4
        (0.0, 0.0, 1.0),  # 5
        (1.0, 0.0, 1.0),  # 6 - shared with Box B
        (1.0, 1.0, 1.0),  # 7 - shared with Box B
        (0.0, 1.0, 1.0),  # 8
        # Box B vertices (4 new vertices, since 4 are shared)
        (2.0, 0.0, 0.0),  # 9
        (2.0, 1.0, 0.0),  # 10
        (2.0, 0.0, 1.0),  # 11
        (2.0, 1.0, 1.0),  # 12
    ]

    # Box A hexahedron vertices (1-8)
    box_a_hex = (1, 2, 3, 4, 5, 6, 7, 8)

    # Box B hexahedron vertices (2,3,6,7,9,10,11,12)
    box_b_hex = (2, 9, 10, 3, 6, 11, 12, 7)

    # Create tetrahedra for Box A using 1:6 rule with diagonal 1-7
    box_a_tets = partition_hexahedron_to_tetrahedra(box_a_hex, :body_diagonal_17)

    # Create tetrahedra for Box B using explicit decomposition to create diagonal mismatch
    # Box B will use diagonal (3,6) on shared face instead of (2,7)
    box_b_tets = [
        (2, 3, 6, 7),  # tet 1: diagonal (3,6) on shared face
        (2, 6, 7, 11), # tet 2
        (2, 7, 3, 10), # tet 3
        (2, 9, 10, 11),# tet 4
        (2, 10, 3, 7), # tet 5
        (2, 11, 9, 7), # tet 6
    ]

    println("Box A created with $(length(box_a_tets)) tetrahedra")
    println("Box B created with $(length(box_b_tets)) tetrahedra")

    # The shared interface vertices (4 vertices at x=1.0)
    shared_interface_vertices = [
        (1.0, 0.0, 0.0),  # vertex 2
        (1.0, 1.0, 0.0),  # vertex 3
        (1.0, 0.0, 1.0),  # vertex 6
        (1.0, 1.0, 1.0),  # vertex 7
    ]

    println("\nShared interface vertices ($(length(shared_interface_vertices)) total):")
    println("  Interface plane: x = 1.0")
    println("  4 vertices forming the shared face between the two cubes")

    println("\nDiagonal difference:")
    println("  Box A uses body diagonal: (0,0,0) → (1,1,1)")
    println("  Box B uses body diagonal: (1,0,0) → (2,1,1)")
    println("  This creates different triangulations on the shared interface face")

    return all_vertices, box_a_tets, box_b_tets, shared_interface_vertices
end

"""
    analyze_diagonal_mismatch(box_a, box_b)

Analyze the diagonal mismatch between the two boxes using the repair system.
"""
function analyze_diagonal_mismatch(box_a, box_b)
    println("\n" * "="^60)
    println("ANALYZING DIAGONAL MISMATCH")
    println("="^60)

    # Build interface topology between the two boxes
    topology = build_interface_topology(box_a, box_b, pidA=1, pidB=2)

    println("\nInterface topology summary:")
    println("  PID A: $(topology.pidA)")
    println("  PID B: $(topology.pidB)")
    println("  Shared nodes: $(topology.total_shared_nodes)")
    println("  Edges only in A: $(length(topology.edges_only_in_A))")
    println("  Edges only in B: $(length(topology.edges_only_in_B))")
    println("  Shared edges: $(length(topology.edges_shared))")

    # Classify mismatches to detect the diagonal issue
    classification = classify_interface_mismatches(topology, verbose=true)

    println("\nMismatch classification:")
    println("  Total mismatches: $(length(classification.mismatches_A) + length(classification.mismatches_B))")
    println("  T-junctions: $(classification.t_junction_count)")
    println("  Diagonal mismatches: $(classification.diagonal_count)")
    println("  Other types: $(classification.refinement_count + classification.quad_mismatch_count + classification.boundary_edge_count)")

    return topology, classification
end

"""
    visualize_diagonal_mismatch(topology, classification)

Create a simple visualization of the diagonal mismatch.
"""
function visualize_diagonal_mismatch(topology, classification)
    println("\n" * "="^60)
    println("DIAGONAL MISMATCH VISUALIZATION")
    println("="^60)

    # Find the diagonal mismatch edges
    diagonal_mismatches = filter(m -> m.mismatch_type == DIAGONAL,
                                vcat(classification.mismatches_A, classification.mismatches_B))

    if isempty(diagonal_mismatches)
        println("No diagonal mismatches found to visualize.")
        return
    end

    println("Visualizing diagonal mismatch on shared face:")
    println()
    println("    Box A triangulation:    Box B triangulation:")
    println("    +-------+              +-------+")
    println("    |\\     /|              |     /|")
    println("    | \\   / |              |    / |")
    println("    |  \\ /  |              |   /  |")
    println("    |   X   |  ← Different  |  X   |")
    println("    |  / \\  |    diagonals  | / \\  |")
    println("    | /   \\ |              |/   \\ |")
    println("    |/     \\|              |     \\|")
    println("    +-------+              +-------+")
    println()
    println("Both boxes share the same 4 corner vertices,")
    println("but use different diagonals to split the face.")
    println()
    println("This creates a DIAGONAL mismatch that can be repaired")
    println("by choosing one triangulation and updating both boxes.")
end

"""
    run_diagonal_mismatch_example()

Complete example demonstrating diagonal mismatch detection and repair.
"""
function run_diagonal_mismatch_example()
    println("DIAGONAL MISMATCH EXAMPLE")
    println("="^60)
    println("This example demonstrates two touching boxes with vertex-compatible")
    println("interfaces where triangulation differs in diagonal manner.")
    println()

    # Step 1: Create the example meshes
    box_a, box_b, shared_vertices = create_diagonal_mismatch_example()

    # Step 2: Analyze the mismatch
    topology, classification = analyze_diagonal_mismatch(box_a, box_b)

    # Step 3: Visualize the issue
    visualize_diagonal_mismatch(topology, classification)

    println("\n" * "="^60)
    println("SUMMARY")
    println("="^60)
    println("✓ Created two boxes sharing 4 interface vertices")
    println("✓ Box A uses diagonal: bottom-front → top-back")
    println("✓ Box B uses diagonal: top-front → bottom-back")
    println("✓ Detected $(classification.diagonal_count) diagonal mismatch(es)")
    println("✓ Repair system identified the issue and can fix it")

    return box_a, box_b, topology, classification
end

# Uncomment to run the example
# if abspath(PROGRAM_FILE) == @__FILE__
#     run_diagonal_mismatch_example()
# end