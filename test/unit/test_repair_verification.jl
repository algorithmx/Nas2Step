"""
    test_repair_verification.jl

Unit tests for repair_verification.jl
Tests compare_conformity, print_conformity_report,
print_improvement_report, and verify_adjacent_interfaces.

Notes:
- Phase 3 mesh verification functions depend on a NasMesh type that is not
  defined in the package yet. We define a minimal NasMesh here so the module
  can be included, and we stub verify_interface_conformity for adjacency tests.
"""

using Test
using Nas2Step

# Define a minimal NasMesh type to satisfy method signatures in repair_verification.jl
mutable struct NasMesh
    all_pid_surfaces::Dict{Int, Vector{Vector{Int}}}
    nodes::Dict{Int, NTuple{3,Float64}}
end

# Bring Phase 3 verification API into scope for tests (not exported by the package)
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "..", "src", "repair", "repair_verification.jl"))

# Small helpers
function capture_stdout(f::Function)
    path = mktemp()[1]
    try
        open(path, "w") do io
            redirect_stdout(io) do
                f()
            end
        end
        return read(path, String)
    finally
        isfile(path) && rm(path; force=true)
    end
end

make_conf(; pair=(1,2), shared_nodes=0, edgesA=0, edgesB=0, shared_edges=0,
          ratioA=1.0, ratioB=1.0, overall=1.0,
          mismatches=0, tj=0, diag=0, refn=0, unk=0) = Dict(
    "interface_pair" => pair,
    "total_shared_nodes" => shared_nodes,
    "total_edges_A" => edgesA,
    "total_edges_B" => edgesB,
    "total_shared_edges" => shared_edges,
    "conformity_ratio_A" => ratioA,
    "conformity_ratio_B" => ratioB,
    "overall_conformity" => overall,
    "total_mismatches" => mismatches,
    "mismatch_types" => Dict(
        "T_JUNCTION" => tj,
        "DIAGONAL" => diag,
        "REFINEMENT" => refn,
        "UNKNOWN" => unk,
    ),
)

# Provide a stub for verify_interface_conformity so adjacency checks don't hit heavy topology.
# We'll define conformity per interface pair via this table.
const _adj_conf = Dict{Tuple{Int,Int},Float64}()
# If any method exists (e.g., from included source), remove it to avoid redefinition warnings.
if @isdefined(verify_interface_conformity)
    for m in methods(verify_interface_conformity)
        # Method signatures look like Tuple{typeof(verify_interface_conformity), NasMesh, Int64, Int64}
        sig = m.sig
        if sig isa Type{<:Tuple} && length(sig.parameters) == 4
            T = sig.parameters
            if T[2] === NasMesh && T[3] === Int && T[4] === Int
                Base.delete_method(m)
            end
        end
    end
end
function verify_interface_conformity(mesh::NasMesh, pid_a::Int, pid_b::Int)
    a, b = pid_a <= pid_b ? (pid_a, pid_b) : (pid_b, pid_a)
    overall = get(_adj_conf, (a, b), 1.0)
    # Return a minimal dictionary matching the expected schema
    return Dict(
        "interface_pair" => (a, b),
        "total_shared_nodes" => 0,
        "total_edges_A" => 0,
        "total_edges_B" => 0,
        "total_shared_edges" => 0,
        "conformity_ratio_A" => overall,
        "conformity_ratio_B" => overall,
        "overall_conformity" => overall,
        "total_mismatches" => 0,
        "mismatch_types" => Dict(
            "T_JUNCTION" => 0,
            "DIAGONAL" => 0,
            "REFINEMENT" => 0,
            "UNKNOWN" => 0,
        ),
    )
end

@testset "Repair Verification - Conformity Comparison" begin
    @testset "Compare Conformity - Improvement" begin
        before = make_conf(overall=0.5, shared_edges=4, mismatches=6)
        after  = make_conf(overall=0.8, shared_edges=7, mismatches=2)
        imp = compare_conformity(before, after)
        @test isapprox(imp["conformity_improvement"], 0.3; atol=1e-12)
        @test imp["mismatches_resolved"] == 4
        @test imp["shared_edges_added"] == 3
        # Smoke: contains before/after dicts
        @test imp["before"] === before
        @test imp["after"] === after
    end

    @testset "Compare Conformity - Degradation" begin
        before = make_conf(overall=0.9, shared_edges=9, mismatches=1)
        after  = make_conf(overall=0.6, shared_edges=6, mismatches=5)
        imp = compare_conformity(before, after)
        @test isapprox(imp["conformity_improvement"], -0.3; atol=1e-12)
        @test imp["mismatches_resolved"] == -4
        @test imp["shared_edges_added"] == -3
    end

    @testset "Compare Conformity - Unchanged" begin
        before = make_conf(overall=0.75, shared_edges=6, mismatches=3)
        after  = make_conf(overall=0.75, shared_edges=6, mismatches=3)
        imp = compare_conformity(before, after)
        @test isapprox(imp["conformity_improvement"], 0.0; atol=1e-12)
        @test imp["mismatches_resolved"] == 0
        @test imp["shared_edges_added"] == 0
    end
end

@testset "Repair Verification - Report Printing" begin
    @testset "Conformity Report - Fully Conforming" begin
        conf = make_conf(pair=(10, 20), shared_nodes=12, edgesA=10, edgesB=10, shared_edges=10,
                         ratioA=1.0, ratioB=1.0, overall=1.0, mismatches=0)
        out = capture_stdout() do
            print_conformity_report(conf, "Unit Test Conformity Report")
        end
        @test occursin("Unit Test Conformity Report", out)
        @test occursin("Interface: (10, 20)", out)
        @test occursin("Shared nodes:  12", out)
        @test occursin("Edges in A:    10", out)
        @test occursin("Edges in B:    10", out)
        @test occursin("Shared edges:  10", out)
        @test occursin("Overall:       100.0%", out)
        @test occursin("Mismatches:    0", out)
        # No per-type mismatch section printed when total_mismatches == 0
        @test !occursin("T-junctions:", out)
    end

    @testset "Conformity Report - With Mismatches" begin
        conf = make_conf(pair=(1, 2), shared_nodes=5, edgesA=8, edgesB=10, shared_edges=6,
                         ratioA=0.75, ratioB=0.6, overall=0.6, mismatches=3,
                         tj=1, diag=1, refn=1, unk=0)
        out = capture_stdout() do
            print_conformity_report(conf)
        end
        @test occursin("Interface: (1, 2)", out)
        @test occursin("Mismatches:    3", out)
        @test occursin("T-junctions:   1", out)
        @test occursin("Diagonal:      1", out)
        @test occursin("Refinement:    1", out)
        @test occursin("Unknown:       0", out)
    end

    @testset "Improvement Report - Status Markers" begin
        before = make_conf(overall=0.4, shared_edges=4, mismatches=8)
        after  = make_conf(overall=0.7, shared_edges=7, mismatches=3)
        out1 = capture_stdout() do
            print_improvement_report(compare_conformity(before, after))
        end
        @test occursin("✓ REPAIR SUCCESSFUL", out1)

        before2 = make_conf(overall=0.7, shared_edges=7, mismatches=3)
        after2  = make_conf(overall=0.7, shared_edges=7, mismatches=3)
        out2 = capture_stdout() do
            print_improvement_report(compare_conformity(before2, after2))
        end
        @test occursin("⚠ REPAIR NEUTRAL", out2)

        before3 = make_conf(overall=0.9, shared_edges=9, mismatches=1)
        after3  = make_conf(overall=0.8, shared_edges=8, mismatches=2)
        out3 = capture_stdout() do
            print_improvement_report(compare_conformity(before3, after3))
        end
        @test occursin("✗ REPAIR DEGRADED", out3)
    end
end

@testset "Repair Verification - Adjacent Interfaces" begin
    # Build a tiny synthetic mesh-like object
    function make_mesh(all_pid_surfaces, nodes)
        return NasMesh(all_pid_surfaces, nodes)
    end

    @testset "Verify Adjacent Interfaces - None" begin
        # Only two PIDs present; no other PIDs to be adjacent.
        nodes = Dict(1 => (0.0,0.0,0.0), 2 => (1.0,0.0,0.0), 3 => (0.0,1.0,0.0))
        faces = Dict(
            1 => [[1,2,3]],
            2 => [[2,3,1]],
        )
        mesh = make_mesh(faces, nodes)
        out = capture_stdout() do
            @test verify_adjacent_interfaces(mesh, (1,2)) === true
        end
        @test occursin("No adjacent interfaces to verify.", out)
    end

    @testset "Verify Adjacent Interfaces - Single" begin
        empty!(_adj_conf)
        # PIDs: 1, 2 (repaired), and 3 adjacent to 1 via shared node 1
        nodes = Dict(
            1 => (0.0,0.0,0.0), 2 => (1.0,0.0,0.0), 3 => (0.0,1.0,0.0),
            4 => (1.0,1.0,0.0), 5 => (2.0,0.0,0.0)
        )
        faces = Dict(
            1 => [[1,2,3]],
            2 => [[2,3,4]],
            3 => [[1,5,3]],  # shares node 1 with PID 1
        )
        mesh = make_mesh(faces, nodes)
        _adj_conf[(1,3)] = 0.95  # Above 0.90 threshold => valid
        out = capture_stdout() do
            ok = verify_adjacent_interfaces(mesh, (1,2))
            @test ok === true
        end
    # Should list adjacent interface(s) and report all valid
        @test occursin("✓ All adjacent interfaces remain valid", out)
    end

    @testset "Verify Adjacent Interfaces - Multiple (with degradation)" begin
        empty!(_adj_conf)
        # PIDs: 1,2 repaired; 3 adjacent to 1; 4 adjacent to both 1 and 2
        nodes = Dict(
            1 => (0.0,0.0,0.0), 2 => (1.0,0.0,0.0), 3 => (0.0,1.0,0.0),
            4 => (1.0,1.0,0.0), 5 => (2.0,0.0,0.0), 6 => (2.0,1.0,0.0)
        )
        faces = Dict(
            1 => [[1,2,3]],            # PID 1 uses nodes 1,2,3
            2 => [[2,3,4]],            # PID 2 uses nodes 2,3,4
            3 => [[1,5,3]],            # Adjacent to PID 1 via node 1 or 3
            4 => [[2,6,4]],            # Adjacent to PID 2 via node 2 or 4
        )
        mesh = make_mesh(faces, nodes)
        _adj_conf[(1,3)] = 0.92
        _adj_conf[(2,4)] = 0.80   # Below threshold -> should fail overall
        out = capture_stdout() do
            ok = verify_adjacent_interfaces(mesh, (1,2))
            @test ok === false
        end
    # Should list adjacent interface(s) and catch a degraded pair
        @test occursin("✗ Some adjacent interfaces degraded!", out)
    end
end


