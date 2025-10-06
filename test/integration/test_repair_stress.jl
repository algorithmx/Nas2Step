"""
    test_repair_stress.jl

Stress/performance-flavored tests gated by ENV["STRESS"] == "1".
These are coarse and avoid strict timing assertions.
"""

module RepairStressTests

using Test
using Main.Nas2StepTestUtils

const RUN_STRESS = get(ENV, "STRESS", "0") == "1"

if RUN_STRESS

@testset "Deep rollback sequences" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
    begin_transaction!(ws)
    # Perform many small edits
    for i in 1:100
        nid = add_node!(ws, (i*0.001, i*0.001, i*0.001))
        # Add then delete a trivial face if enough nodes exist
        if length(ws.working_nodes) >= 3
            ids = collect(keys(ws.working_nodes))[1:3]
            add_face!(ws, first(keys(ws.working_faces)), ids)
            delete_face!(ws, first(keys(ws.working_faces)), length(ws.working_faces[first(keys(ws.working_faces))]))
        end
    end
    @test length(ws.modifications) >= 100
    rollback_transaction!(ws)
    @test ws.transaction_active == false
end

else
@testset "Stress tests gated off" begin
    @test true
end
end # RUN_STRESS

end # module
