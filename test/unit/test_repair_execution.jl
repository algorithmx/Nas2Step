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

# Import CoordinateKeys functions for consistent EdgeKey handling
import Nas2Step: create_edge_key_int
import Nas2Step.CoordinateKeys

# Alias package types used by included execution code
const EdgeKey = Main.Nas2Step.EdgeKey

# Minimal NasMesh placeholder to satisfy method signatures in execution file
mutable struct NasMesh
    all_pid_surfaces::Dict{Int, Vector{Vector{Int}}}
    nodes::Dict{Int, NTuple{3,Float64}}
end

"""
Local mirror of EdgeInsertionPlan matching src/repair/repair_planning.jl.
Only the fields used by repair_execution.jl are functionally relevant here.
"""
mutable struct EdgeInsertionPlan
    target_triangle::Int
    insert_edge::Main.Nas2Step.EdgeKey
    split_type::Symbol
    new_nodes::Vector{NTuple{3,Float64}}
    existing_nodes::Vector{Int}
    old_triangles::Vector{NTuple{9,Float64}}
    replacement_triangles::Vector{NTuple{9,Float64}}
    min_angle_before::Float64
    min_angle_after::Float64
    quality_acceptable::Bool
    depends_on::Vector{Int}
    violates_constraints::Bool
    constraint_violations::Vector{String}
    is_feasible::Bool
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

"""
Small helpers to extract coordinates and flatten triangles to NTuple{9,Float64}.
"""
coords_of(ws, ids::Vector{Int}) = [ws.working_nodes[i] for i in ids]
to_flat(tri::Vector{NTuple{3,Float64}})::NTuple{9,Float64} = (tri[1]..., tri[2]..., tri[3]...)

"""
Convenience builder for a quad retriangulation plan.
"""
function make_quad_plan(edge::Main.Nas2Step.EdgeKey,
                        old1::Vector{NTuple{3,Float64}}, old2::Vector{NTuple{3,Float64}},
                        new1::Vector{NTuple{3,Float64}}, new2::Vector{NTuple{3,Float64}})
    return EdgeInsertionPlan(
        0,                      # target_triangle (unused here)
        edge,                   # insert_edge
        :quad_retriangulation,  # split_type
        NTuple{3,Float64}[],    # new_nodes
        Int[],                  # existing_nodes
        [to_flat(old1), to_flat(old2)],
        [to_flat(new1), to_flat(new2)],
        0.0, 0.0, true,         # quality placeholders
        Int[],                  # depends_on
        false,                  # violates_constraints
        String[],               # constraint_violations
        true                    # is_feasible
    )
end

"""
Convenience builder for a triangle-split style plan (uses one old triangle).
"""
function make_split_plan(edge::Main.Nas2Step.EdgeKey,
                         old1::Vector{NTuple{3,Float64}},
                         newA::Vector{NTuple{3,Float64}}, newB::Vector{NTuple{3,Float64}};
                         split_type::Symbol = :triangle_split,
                         new_nodes::Vector{NTuple{3,Float64}} = NTuple{3,Float64}[])
    return EdgeInsertionPlan(
        0,
        edge,
        split_type,
        new_nodes,
        Int[],
        [to_flat(old1)],
        [to_flat(newA), to_flat(newB)],
        0.0, 0.0, true,
        Int[],
        false,
        String[],
        true
    )
end

@testset "Repair Execution - Quad Retriangulation" begin
    @testset "Apply Quad Retriangulation - Success" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        # Use faces forming a quad: [1,2,3] and [1,3,4]
        triA = coords_of(ws, [1,2,3])
        triB = coords_of(ws, [1,3,4])
        new1 = coords_of(ws, [1,2,4])
        new2 = coords_of(ws, [2,3,4])
        plan = make_quad_plan(create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]),
                              triA, triB, new1, new2)

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
            0,
            create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]),
            :quad_retriangulation,
            NTuple{3,Float64}[],
            Int[],
            [to_flat(triA)],               # only one triangle to replace
            NTuple{9,Float64}[],            # no new triangles
            0.0, 0.0, true,
            Int[], false, String[], true
        )
        faces_before = deepcopy(ws.working_faces)
        ok = apply_quad_retriangulation!(ws, plan_bad_count, pid)
        @test ok === false
        @test ws.working_faces == faces_before  # No changes when counts invalid

        # Wrong plan type throws
        plan_wrong_type = EdgeInsertionPlan(
            0,
            create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]),
            :triangle_split,
            NTuple{3,Float64}[],
            Int[],
            [to_flat(triA), to_flat(triA)],
            [to_flat(triA), to_flat(triA)],
            0.0, 0.0, true,
            Int[], false, String[], true
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
        plan = make_quad_plan(create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]),
                              triA, triB, coords_of(ws, [2,3,4]), new_bad)

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
        good = make_quad_plan(create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]),
                              triA, triB, coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4]))
    split = make_split_plan(create_edge_key_int(ws.working_nodes[1], ws.working_nodes[2]),
                triA, coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4]); split_type=:bisect)
        # Unknown type should return false
        unknown = EdgeInsertionPlan(0, create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]), :foobar,
                                    NTuple{3,Float64}[], Int[], [to_flat(triA), to_flat(triB)],
                                    [to_flat(coords_of(ws,[1,2,4])), to_flat(coords_of(ws,[2,3,4]))],
                                    0.0, 0.0, true, Int[], false, String[], true)

        @test apply_edge_insertion_plan!(ws, good, pid) === true
        @test apply_edge_insertion_plan!(ws, split, pid) === true
        @test apply_edge_insertion_plan!(ws, unknown, pid) === false
    end
end

@testset "Repair Execution - Full Repair Plan" begin
    @testset "Apply Repair Plan - All Succeed" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        # Two independent quads on PID 1: (1,2,3)+(1,3,4) and (5,6,7)+(5,7,8)
        p1 = make_quad_plan(create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]),
                            coords_of(ws,[1,2,3]), coords_of(ws,[1,3,4]),
                            coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4]))
        p2 = make_quad_plan(create_edge_key_int(ws.working_nodes[5], ws.working_nodes[7]),
                            coords_of(ws,[5,6,7]), coords_of(ws,[5,7,8]),
                            coords_of(ws,[5,6,8]), coords_of(ws,[6,7,8]))

        plan = RepairPlan((pid, 2), :subdivide_A, [p1, p2], 2, 0, 0, 0.0, 0.0, true, String[])

        ok = apply_repair_plan!(ws, plan)
        @test ok === true
        @test ws.transaction_active == false
    end

    @testset "Apply Repair Plan - Some Fail (commit partial)" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        pid = 1
        # First plan good, second plan uses a fake node coordinate to force failure
        good = make_quad_plan(create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]),
                              coords_of(ws,[1,2,3]), coords_of(ws,[1,3,4]),
                              coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4]))
        bad = make_quad_plan(create_edge_key_int(ws.working_nodes[5], ws.working_nodes[7]),
                             coords_of(ws,[5,6,7]), coords_of(ws,[5,7,8]),
                             [ws.working_nodes[5], ws.working_nodes[6], (9.9,9.9,9.9)], coords_of(ws,[6,7,8]))

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
        bad1 = make_quad_plan(create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]),
                              coords_of(ws,[1,2,3]), coords_of(ws,[1,3,4]),
                              [ws.working_nodes[1], ws.working_nodes[2], (99.0,99.0,99.0)], coords_of(ws,[2,3,4]))
        bad2 = make_quad_plan(create_edge_key_int(ws.working_nodes[5], ws.working_nodes[7]),
                              coords_of(ws,[5,6,7]), coords_of(ws,[5,7,8]),
                              [ws.working_nodes[5], ws.working_nodes[6], (77.0,77.0,77.0)], coords_of(ws,[6,7,8]))

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
    infeas = EdgeInsertionPlan(0, create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]), :quad_retriangulation,
                                   NTuple{3,Float64}[], Int[], [to_flat(triA), to_flat(triB)], [to_flat(triA), to_flat(triB)],
                                   0.0, 0.0, true, Int[], false, String[], false)  # is_feasible=false
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
    p = make_quad_plan(create_edge_key_int(ws.working_nodes[1], ws.working_nodes[3]),
                              triA, triB, coords_of(ws,[1,2,4]), coords_of(ws,[2,3,4]))
        plan = RepairPlan((fake_pid, 2), :subdivide_A, [p], 1, 0, 0, 0.0, 0.0, true, String[])

        @test_throws Any apply_repair_plan!(ws, plan)
        # After exception, transaction should have been rolled back and state unchanged
        @test Nas2StepTestUtils.workspace_state_equals(ws, initial_state)
    end
end

end # module


