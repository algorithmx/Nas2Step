"""
    test_classification.jl

Test suite for edge mismatch classification.
"""

using Test
using Nas2Step.NasMesh
using Nas2Step.Repair: build_interface_topology, classify_interface_mismatches, MismatchType

@testset "T-Junction Classification" begin
    # Create a mesh with a known T-junction
    #
    #   Surface A (PID 1): One large triangle
    #   (0,0,0) *-------* (2,0,0)
    #           | \      |
    #           |  \     |
    #           |   \    |
    #   (0,2,0) *----*---* (2,2,0)
    #
    #   Surface B (PID 2): Two smaller triangles, sharing a node at (1,0,0)
    #   (0,0,0) *-------* (1,0,0) *-------* (2,0,0)
    #           | \      |      /|
    #           |  \     |     / |
    #           |   \    |    /  |
    #   (0,2,0) *----*---*---*---* (2,2,0)

    nodes = Dict(
        1 => (0.0, 0.0, 0.0),
        2 => (2.0, 0.0, 0.0),
        3 => (0.0, 2.0, 0.0),
        4 => (1.0, 0.0, 0.0)
    )

    pid_surfaces = Dict(
        1 => [ [1, 2, 3] ],  # Surface A
        2 => [ [1, 4, 3], [4, 2, 3] ]   # Surface B
    )

    mesh = NasMesh(nodes, pid_surfaces)

    # Build topology and classify
    topology = build_interface_topology(mesh, 1, 2)
    classification = classify_interface_mismatches(topology)

    # Verification
    @test classification.t_junction_count == 1
    @test classification.diagonal_count == 0
    @test classification.refinement_count == 0
    @test classification.unknown_count == 0

    # Verify the details of the mismatch
    @test length(classification.mismatches_A) == 1
    @test length(classification.mismatches_B) == 0

    mismatch = classification.mismatches_A[1]
    @test mismatch.mismatch_type == MismatchType.T_JUNCTION
    @test mismatch.present_in == :B_only # The edge [1,3] is in B, but not A
    @test mismatch.should_be_in == 1 # A needs to be subdivided
    @test length(mismatch.hanging_nodes) == 1
    @test mismatch.hanging_nodes[1] == (1.0, 0.0, 0.0)
end


@testset "Diagonal Mismatch Classification" begin
    # Create a mesh with a known diagonal mismatch
    #
    #   Two surfaces covering the same square area, but with different internal diagonals.
    #
    #      (0,1) *-------* (1,1)
    #            | \     | 
    #            |  \    |
    #            |   \   |   Surface A: Diagonal from (0,0) to (1,1)
    #            |    \  |
    #            |     \ |
    #      (0,0) *-------* (1,0)
    # 
    #      (0,1) *-------* (1,1)
    #            |     / |
    #            |    /  |
    #            |   /   |   Surface B: Diagonal from (1,0) to (0,1)
    #            |  /    |
    #            | /     |
    #      (0,0) *-------* (1,0)

    nodes = Dict(
        1 => (0.0, 0.0, 0.0),
        2 => (1.0, 0.0, 0.0),
        3 => (1.0, 1.0, 0.0),
        4 => (0.0, 1.0, 0.0)
    )

    pid_surfaces = Dict(
        1 => [ [1, 2, 3], [1, 3, 4] ],  # Surface A
        2 => [ [1, 2, 4], [2, 3, 4] ]   # Surface B
    )

    mesh = NasMesh(nodes, pid_surfaces)

    # Build topology and classify
    topology = build_interface_topology(mesh, 1, 2)
    classification = classify_interface_mismatches(topology)

    # Verification
    @test classification.t_junction_count == 0
    @test classification.diagonal_count == 2 # One edge for each direction
    @test classification.refinement_count == 0
    @test classification.unknown_count == 0

    # Verify the details of the mismatch
    @test length(classification.mismatches_A) == 1
    @test length(classification.mismatches_B) == 1

    mismatch_A = classification.mismatches_A[1]
    @test mismatch_A.mismatch_type == MismatchType.DIAGONAL
    @test isempty(mismatch_A.hanging_nodes)

    mismatch_B = classification.mismatches_B[1]
    @test mismatch_B.mismatch_type == MismatchType.DIAGONAL
    @test isempty(mismatch_B.hanging_nodes)
end
