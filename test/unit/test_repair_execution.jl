"""
    test_repair_execution.jl

Unit tests for repair_execution.jl
Tests apply_quad_retriangulation!, apply_edge_insertion_plan!,
and apply_repair_plan! with success/failure and rollback scenarios.
"""

module RepairExecutionTests

using Test
using Main.Nas2StepTestUtils

# Import workspace API and types from Main (already included by runtests)
import Main: RepairWorkspace, begin_transaction!, commit_transaction!, rollback_transaction!
import Main: delete_face!, add_face!, add_node!, get_face_by_nodes, get_node_id_by_coords
import Main: export_modified_mesh, print_workspace_stats

# Alias package types used by included execution code
const EdgeKey = Main.Nas2Step.EdgeKey

# Minimal NasMesh placeholder to satisfy method signatures in execution file
mutable struct NasMesh
    all_pid_surfaces::Dict{Int, Vector{Vector{Int}}}
    nodes::Dict{Int, NTuple{3,Float64}}
end

# Define minimal local types to match repair_execution.jl expectations within this module scope.
mutable struct EdgeInsertionPlan
    edge_key::Main.Nas2Step.EdgeKey
    plan_type::Symbol
    triangles_to_replace::Vector{Vector{NTuple{3,Float64}}}  # each is length-3 coords
    new_triangles::Vector{Vector{NTuple{3,Float64}}}
    new_nodes::Vector{NTuple{3,Float64}}
    predicted_min_quality::Float64
    predicted_max_quality_loss::Float64
    is_feasible::Bool
    feasibility_issues::Vector{String}
end

mutable struct RepairPlan
    interface_pair::Tuple{Int,Int}
    repair_direction::Symbol  # :subdivide_A | :subdivide_B
    insertion_sequence::Vector{EdgeInsertionPlan}
    total_edges_to_insert::Int
    total_triangles_to_split::Int
    total_nodes_to_add::Int
    predicted_min_quality::Float64
    predicted_max_quality_loss::Float64
    is_feasible::Bool
    feasibility_issues::Vector{String}
end

# Bring the execution API into this module (methods will bind to types above)
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "..", "src", "repair", "repair_execution.jl"))

# Small helpers
coords_of(ws, ids::Vector{Int}) = [ws.working_nodes[i] for i in ids]

@testset "Repair Execution - Quad Retriangulation" begin
    @testset "Apply Quad Retriangulation - Success" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        # Use faces forming a quad: [1,2,3] and [1,3,4]
        triA = coords_of(ws, [1,2,3])
        triB = coords_of(ws, [1,3,4])
        new1 = coords_of(ws, [1,2,4])
        new2 = coords_of(ws, [2,3,4])
        plan = EdgeInsertionPlan(
            Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]),
            :quad_retriangulation,
            [triA, triB],
            [new1, new2],
            NTuple{3,Float64}[],
            0.0,
            0.0,
            true,
            String[]
        )

        faces_before = length(ws.working_faces[pid])
        mods_before = length(ws.modifications)
        ok = apply_quad_retriangulation!(ws, plan, pid)
        @test ok === true
        # Net face count should be unchanged (2 deleted, 2 added)
        @test length(ws.working_faces[pid]) == faces_before
        # Four modifications: 2 deletions + 2 additions
        @test length(ws.modifications) == mods_before + 4
    end

    @testset "Apply Quad Retriangulation - Early Invalid Plan" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        triA = coords_of(ws, [1,2,3])
        plan_bad_count = EdgeInsertionPlan(
            Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]),
            :quad_retriangulation,
            [triA],  # only one triangle to replace
            Vector{Vector{NTuple{3,Float64}}}(),
            NTuple{3,Float64}[], 0.0, 0.0, true, String[]
        )
        faces_before = deepcopy(ws.working_faces)
        ok = apply_quad_retriangulation!(ws, plan_bad_count, pid)
        @test ok === false
        @test ws.working_faces == faces_before  # No changes when counts invalid

        # Wrong plan type throws
        plan_wrong_type = EdgeInsertionPlan(
            Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]),
            :triangle_split,
            [triA, triA],
            [triA, triA],
            NTuple{3,Float64}[], 0.0, 0.0, true, String[]
        )
        @test_throws ErrorException apply_quad_retriangulation!(ws, plan_wrong_type, pid)
    end

    @testset "Apply Quad Retriangulation - New Triangle Node Missing" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        triA = coords_of(ws, [1,2,3])
        triB = coords_of(ws, [1,3,4])
        # One of the new triangle vertices uses a non-existent coordinate
        fake = (100.0, 100.0, 100.0)
        new_bad = [ws.working_nodes[1], ws.working_nodes[2], fake]
        plan = EdgeInsertionPlan(
            Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]),
            :quad_retriangulation,
            [triA, triB],
            [coords_of(ws, [2,3,4]), new_bad],
            NTuple{3,Float64}[], 0.0, 0.0, true, String[]
        )

        begin_transaction!(ws)
        faces_before = length(ws.working_faces[pid])
        ok = apply_quad_retriangulation!(ws, plan, pid)
        @test ok === false
        # Two deletions and one successful addition happened before failure => net -1
        @test length(ws.working_faces[pid]) == faces_before - 1
        rollback_transaction!(ws)  # Clean up to restore workspace
    end
end

@testset "Repair Execution - Edge Insertion" begin
    @testset "Apply Edge Insertion Plan dispatch" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        triA = coords_of(ws, [1,2,3])
        triB = coords_of(ws, [1,3,4])
        good = EdgeInsertionPlan(Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]), :quad_retriangulation,
                                 [triA, triB], [coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])
    bad_type = EdgeInsertionPlan(Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]), :triangle_split,
                                     [triA, triB], [coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])
    unknown = EdgeInsertionPlan(Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]), :foobar,
                                    [triA, triB], [coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])

        @test apply_edge_insertion_plan!(ws, good, pid) === true
        @test apply_edge_insertion_plan!(ws, bad_type, pid) === false
        @test apply_edge_insertion_plan!(ws, unknown, pid) === false
    end
end

@testset "Repair Execution - Full Repair Plan" begin
    @testset "Apply Repair Plan - All Succeed" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        # Two independent quads on PID 1: (1,2,3)+(1,3,4) and (5,6,7)+(5,7,8)
        p1 = EdgeInsertionPlan(
            Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]), :quad_retriangulation,
            [coords_of(ws,[1,2,3]), coords_of(ws,[1,3,4])],
            [coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])
        p2 = EdgeInsertionPlan(
            Main.Nas2Step.EdgeKey(ws.working_nodes[5], ws.working_nodes[7]), :quad_retriangulation,
            [coords_of(ws,[5,6,7]), coords_of(ws,[5,7,8])],
            [coords_of(ws,[5,6,8]), coords_of(ws,[6,7,8])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])

        plan = RepairPlan((pid, 2), :subdivide_A, [p1, p2], 2, 0, 0, 0.0, 0.0, true, String[])

        ok = apply_repair_plan!(ws, plan)
        @test ok === true
        @test ws.transaction_active == false
    end

    @testset "Apply Repair Plan - Some Fail (commit partial)" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        # First plan good, second plan uses a fake node coordinate to force failure
        good = EdgeInsertionPlan(
            Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]), :quad_retriangulation,
            [coords_of(ws,[1,2,3]), coords_of(ws,[1,3,4])],
            [coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])
        bad = EdgeInsertionPlan(
            Main.Nas2Step.EdgeKey(ws.working_nodes[5], ws.working_nodes[7]), :quad_retriangulation,
            [coords_of(ws,[5,6,7]), coords_of(ws,[5,7,8])],
            [[ws.working_nodes[5], ws.working_nodes[6], (9.9,9.9,9.9)], coords_of(ws,[6,7,8])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])

        plan = RepairPlan((pid, 2), :subdivide_A, [good, bad], 2, 0, 0, 0.0, 0.0, true, String[])
        ok = apply_repair_plan!(ws, plan)
        @test ok === false  # partial success committed
        @test ws.transaction_active == false
    end

    @testset "Apply Repair Plan - All Fail (rollback)" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        initial_state = Nas2StepTestUtils.capture_workspace_state(ws)
        # Both plans marked feasible but constructed to fail during execution
        bad1 = EdgeInsertionPlan(
            Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]), :quad_retriangulation,
            [coords_of(ws,[1,2,3]), coords_of(ws,[1,3,4])],
            [[ws.working_nodes[1], ws.working_nodes[2], (99.0,99.0,99.0)], coords_of(ws,[2,3,4])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])
        bad2 = EdgeInsertionPlan(
            Main.Nas2Step.EdgeKey(ws.working_nodes[5], ws.working_nodes[7]), :quad_retriangulation,
            [coords_of(ws,[5,6,7]), coords_of(ws,[5,7,8])],
            [[ws.working_nodes[5], ws.working_nodes[6], (77.0,77.0,77.0)], coords_of(ws,[6,7,8])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])

        plan = RepairPlan((pid, 2), :subdivide_A, [bad1, bad2], 2, 0, 0, 0.0, 0.0, true, String[])
        ok = apply_repair_plan!(ws, plan)
        @test ok === false
        # On all failures, transaction is rolled back and state restored
        @test Nas2StepTestUtils.workspace_state_equals(ws, initial_state)
    end

    @testset "Apply Repair Plan - No Feasible Plans" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        triA = coords_of(ws, [1,2,3])
        triB = coords_of(ws, [1,3,4])
    infeas = EdgeInsertionPlan(Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]), :quad_retriangulation,
                                   [triA, triB], [triA, triB], NTuple{3,Float64}[], 0.0, 0.0, false, ["not feasible"])  # is_feasible=false
        plan = RepairPlan((pid, 2), :subdivide_A, [infeas], 1, 0, 0, 0.0, 0.0, true, String[])
        faces_before = deepcopy(ws.working_faces)
        ok = apply_repair_plan!(ws, plan)
        @test ok === false
        # No transaction started, workspace unchanged
        @test ws.transaction_active == false
        @test ws.working_faces == faces_before
    end
end

@testset "Repair Execution - Rollback Testing" begin
    @testset "Exception During Repair triggers rollback and rethrow" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        initial_state = Nas2StepTestUtils.capture_workspace_state(ws)
        # Choose a target pid that DOES NOT exist to trigger KeyError inside execution
        fake_pid = 999
        triA = coords_of(ws, [1,2,3])
        triB = coords_of(ws, [1,3,4])
    p = EdgeInsertionPlan(Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]), :quad_retriangulation,
                              [triA, triB], [coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])
        plan = RepairPlan((fake_pid, 2), :subdivide_A, [p], 1, 0, 0, 0.0, 0.0, true, String[])

        @test_throws Any apply_repair_plan!(ws, plan)
        # After exception, transaction should have been rolled back and state unchanged
        @test Nas2StepTestUtils.workspace_state_equals(ws, initial_state)
    end
end

end # module


