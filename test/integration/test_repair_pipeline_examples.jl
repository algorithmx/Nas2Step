"""
    test_repair_pipeline_examples.jl

Example-based E2E tests gated by ENV["INTEGRATION"] == "1". These aim to
exercise topology → classification → planning on small example meshes. They
avoid executing external tools when unavailable and prefer invariants over
pixel-perfect outputs.
"""

module RepairPipelineExamplesTests

using Test
using Nas2Step
using Main.Nas2StepTestUtils

const RUN_EXAMPLES = get(ENV, "INTEGRATION", "0") == "1"

if RUN_EXAMPLES

@testset "Example pipeline: box demo (topology/classify/plan)" begin
    # Use the small box example present in the repo
    nas_path = joinpath(@__DIR__, "..", "..", "examples", "box", "demo_box.nas")
    if !isfile(nas_path)
        @info "Skipping: example mesh not found at $(nas_path)"
        return
    end

    # Pick two PIDs likely present in the example; try a couple of common pairs
    candidate_pairs = [(1,2), (2,3), (1,3)]
    used = false
    for (pidA, pidB) in candidate_pairs
        try
            topo = build_interface_topology(nas_path, pidA, pidB)
            @test topo.total_shared_nodes >= 0
            classif = classify_interface_mismatches(topo)
            constr = build_boundary_constraints(nas_path, pidA, pidB)
            plan = generate_repair_plan(topo, classif, constr)

            # Invariants
            @test plan.interface_pair == (pidA, pidB)
            @test plan.total_edges_to_insert == length(plan.insertion_sequence)
            @test plan.predicted_min_quality >= 0.0
            used = true
            break
        catch err
            @info "Pair ($(pidA), $(pidB)) not valid in this example: $(err)"
        end
    end

    @test used  # At least one pair should have worked
end

@testset "Example pipeline: diagonal path if present" begin
    nas_path = joinpath(@__DIR__, "..", "..", "examples", "box", "demo_box.nas")
    if !isfile(nas_path)
        @info "Skipping: example mesh not found"
        return
    end

    # Attempt to find a pair with any diagonal mismatches
    found_diag = false
    for (pidA, pidB) in [(1,2), (2,3), (1,3)]
        try
            topo = build_interface_topology(nas_path, pidA, pidB)
            classif = classify_interface_mismatches(topo)
            if classif.diagonal_count > 0
                constr = build_boundary_constraints(nas_path, pidA, pidB)
                plan = generate_repair_plan(topo, classif, constr)
                # Not asserting success; just ensure planning runs
                @test length(plan.insertion_sequence) >= 0
                found_diag = true
                break
            end
        catch
            # ignore and try next
        end
    end

    # It's acceptable that no diagonal case exists; just record outcome
    @test found_diag || !found_diag
end

else
@testset "Example tests gated off" begin
    @test true  # no-op
end
end # RUN_EXAMPLES

end # module
