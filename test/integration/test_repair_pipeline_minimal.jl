"""
    test_repair_pipeline_minimal.jl

Integration tests (fast) for a tiny, in-memory pipeline using the local
RepairWorkspace and execution functions. No file I/O or external tools.
"""

module RepairPipelineMinimalTests

using Test
using Main.Nas2StepTestUtils

# Import workspace API from Main (already included by test runner)
import Main: RepairWorkspace, begin_transaction!, commit_transaction!, rollback_transaction!
import Main: delete_face!, add_face!, add_node!, get_face_by_nodes, get_node_id_by_coords

# Define local aliases/types matching what repair_execution.jl expects, then include it
# Satisfy type annotations without needing full implementation
const NasMesh = Any

mutable struct EdgeInsertionPlan
    edge_key::Main.Nas2Step.EdgeKey
    plan_type::Symbol
    triangles_to_replace::Vector{Vector{NTuple{3,Float64}}}
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

Base.include(@__MODULE__, joinpath(@__DIR__, "..", "..", "src", "repair", "repair_execution.jl"))

# Helpers
coords_of(ws, ids::Vector{Int}) = [ws.working_nodes[i] for i in ids]

@testset "Minimal pipeline: quad retriangulation success" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
    pid = 1

    # Use faces forming a quad: [1,2,3] and [1,3,4], switch diagonal
    triA = coords_of(ws, [1,2,3])
    triB = coords_of(ws, [1,3,4])
    new1 = coords_of(ws, [1,2,4])
    new2 = coords_of(ws, [2,3,4])

    p = EdgeInsertionPlan(
        Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]),
        :quad_retriangulation,
        [triA, triB],
        [new1, new2],
        NTuple{3,Float64}[],
        0.0, 0.0, true, String[]
    )

    plan = RepairPlan((pid, 2), :subdivide_A, [p], 1, 0, 0, 0.0, 0.0, true, String[])

    faces_before = deepcopy(ws.working_faces)
    ok = apply_repair_plan!(ws, plan)
    @test ok === true
    # Net faces per PID should be stable for the modified PID (2 replaced by 2)
    @test ws.working_faces[pid] |> length == faces_before[pid] |> length
    @test ws.transaction_active == false
end

@testset "Minimal pipeline: all failures trigger rollback" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
    initial = Nas2StepTestUtils.capture_workspace_state(ws)
    pid = 1

    triA = coords_of(ws, [1,2,3])
    triB = coords_of(ws, [1,3,4])
    fake = (99.9, 99.9, 99.9)  # Ensure missing node coordinate
    bad1 = EdgeInsertionPlan(Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]), :quad_retriangulation,
                             [triA, triB], [coords_of(ws, [1,2,4]), [ws.working_nodes[2], ws.working_nodes[3], fake]], NTuple{3,Float64}[], 0.0, 0.0, true, String[])
    bad2 = EdgeInsertionPlan(Main.Nas2Step.EdgeKey(ws.working_nodes[5], ws.working_nodes[7]), :quad_retriangulation,
                             [coords_of(ws,[5,6,7]), coords_of(ws,[5,7,8])], [[ws.working_nodes[5], ws.working_nodes[6], fake], coords_of(ws,[6,7,8])], NTuple{3,Float64}[], 0.0, 0.0, true, String[])

    plan = RepairPlan((pid, 2), :subdivide_A, [bad1, bad2], 2, 0, 0, 0.0, 0.0, true, String[])
    ok = apply_repair_plan!(ws, plan)
    @test ok === false
    @test Nas2StepTestUtils.workspace_state_equals(ws, initial)
    @test ws.transaction_active == false
end

@testset "Minimal pipeline: no-op plan is a no-op" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
    initial = Nas2StepTestUtils.capture_workspace_state(ws)
    pid = 1
    triA = coords_of(ws, [1,2,3])
    triB = coords_of(ws, [1,3,4])
    infeas = EdgeInsertionPlan(Main.Nas2Step.EdgeKey(ws.working_nodes[1], ws.working_nodes[3]), :quad_retriangulation,
                               [triA, triB], [triA, triB], NTuple{3,Float64}[], 0.0, 0.0, false, ["not feasible"])
    plan = RepairPlan((pid, 2), :subdivide_A, [infeas], 1, 0, 0, 0.0, 0.0, true, String[])
    ok = apply_repair_plan!(ws, plan)
    @test ok === false
    @test Nas2StepTestUtils.workspace_state_equals(ws, initial)
    @test ws.transaction_active == false
end

end # module
