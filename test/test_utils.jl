module Nas2StepTestUtils

using Nas2Step
using Test

export ckey

"""
    ckey(p)

Round coordinates to match production code's 4-digit precision for coordinate keys.
This is a convenience alias for Nas2Step.coordinate_key_float to maintain test readability.
Use this helper in tests to avoid duplicate method definitions across files.
"""
ckey(p::NTuple{3,Float64}) = (round(p[1]; digits=4), round(p[2]; digits=4), round(p[3]; digits=4))

"""
    create_minimal_workspace()

Return a small synthetic workspace using the RepairWorkspace type included by the test runner.
"""
function create_minimal_workspace()
    # Simple cube-like set of 8 nodes
    working_nodes = Dict{Int,NTuple{3,Float64}}(
        1 => (0.0, 0.0, 0.0),
        2 => (1.0, 0.0, 0.0),
        3 => (1.0, 1.0, 0.0),
        4 => (0.0, 1.0, 0.0),
        5 => (0.0, 0.0, 1.0),
        6 => (1.0, 0.0, 1.0),
        7 => (1.0, 1.0, 1.0),
        8 => (0.0, 1.0, 1.0),
    )

    # A couple of triangular faces across two PIDs
    working_faces = Dict{Int,Vector{Vector{Int}}}(
        1 => [[1,2,3], [1,3,4], [5,6,7], [5,7,8]],
        2 => [[2,6,7], [2,7,3]],
    )

    return Nas2Step.RepairWorkspace(
        "synthetic_cube.nas",
        working_faces,
        working_nodes,
        8,
        Nas2Step.MeshModification[],
        0,
        false,
        0, 0, 0,
    )
end

"""
    capture_workspace_state(ws)

Snapshot a subset of workspace fields for equality checks.
"""
function capture_workspace_state(ws)
    return (
        working_faces = deepcopy(ws.working_faces),
        working_nodes = deepcopy(ws.working_nodes),
        max_node_id = ws.max_node_id,
        modifications = deepcopy(ws.modifications),
        faces_added = ws.faces_added,
        faces_deleted = ws.faces_deleted,
        nodes_added = ws.nodes_added,
        checkpoint_mod_count = ws.checkpoint_mod_count,
        transaction_active = ws.transaction_active,
    )
end

"""
    workspace_state_equals(ws, state)

Compare workspace against a captured state.
"""
function workspace_state_equals(ws, state)
    return (
        ws.working_faces == state.working_faces &&
        ws.working_nodes == state.working_nodes &&
        ws.max_node_id == state.max_node_id &&
        length(ws.modifications) == length(state.modifications) &&
        ws.faces_added == state.faces_added &&
        ws.faces_deleted == state.faces_deleted &&
        ws.nodes_added == state.nodes_added
    )
end

"""
    random_point_3d()

Generate a random 3D point for simple property tests.
"""
random_point_3d() = (rand(), rand(), rand())

end # module

