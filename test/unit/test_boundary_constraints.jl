"""
    test_boundary_constraints.jl

Unit tests for boundary_constraints.jl
Tests build_boundary_constraints, check_constraint_violations,
corner node detection, external edge identification, and locked node/edge queries.
"""

using Test
using Nas2Step
using JSON
using .Nas2StepTestUtils: ckey
using Nas2Step: create_edge_key_int


"""
    make_constraints(; pidA=1, pidB=2, extA=EdgeKey[], extB=EdgeKey[], corners=NTuple{3,Float64}[], features=EdgeKey[])

Construct a BoundaryConstraints instance for unit testing without invoking Gmsh.
"""
function make_constraints(; pidA::Int=1, pidB::Int=2,
    extA::Vector{EdgeKey}=EdgeKey[], extB::Vector{EdgeKey}=EdgeKey[],
    corners::Vector{NTuple{3,Float64}}=NTuple{3,Float64}[], features::Vector{EdgeKey}=EdgeKey[])
    pidA_external = Set(extA)
    pidB_external = Set(extB)
    corner_nodes = Set(corners)
    feature_edges = Set(features)
    locked_nodes = copy(corner_nodes)
    locked_edges = union(pidA_external, pidB_external, feature_edges)
    return Nas2Step.BoundaryConstraints(
        pidA, pidB,
        pidA_external, pidB_external,
        corner_nodes,
        feature_edges,
        locked_nodes,
        locked_edges,
        length(pidA_external),
        length(pidB_external),
        length(corner_nodes),
        length(locked_edges),
    )
end

@testset "Boundary Constraints - Build Constraints" begin
    @testset "Build Boundary Constraints (demo_box.nas)" begin
        # Use the small example mesh with two regions (PID 1 and 2)
        nas_path = normpath(joinpath(@__DIR__, "..", "..", "examples", "box", "demo_box.nas"))
        @test isfile(nas_path)

        bc = Nas2Step.build_boundary_constraints(nas_path, 1, 2)

        # Basic type and pid checks
        @test bc isa Nas2Step.BoundaryConstraints
        @test bc.pidA == 1
        @test bc.pidB == 2

        # Invariants: counts match set sizes
        @test bc.total_external_edges_A == length(bc.pidA_external_edges)
        @test bc.total_external_edges_B == length(bc.pidB_external_edges)
        @test bc.total_corner_nodes == length(bc.corner_nodes)
        @test bc.total_locked_edges == length(bc.locked_edges)

        # Locked sets are derived from components
        @test bc.locked_nodes == bc.corner_nodes
        @test bc.locked_edges == union(bc.pidA_external_edges, bc.pidB_external_edges, bc.feature_edges)

        # With only two PIDs in this fixture, there should be no corner nodes (need 3+ PIDs)
        @test isempty(bc.corner_nodes)
    end
end

@testset "Boundary Constraints - Constraint Violations" begin
    @testset "Check Constraint Violations" begin
        # Synthetic geometry
        a = ckey((0.0, 0.0, 0.0))
        b = ckey((1.0, 0.0, 0.0))
        c = ckey((0.0, 1.0, 0.0))
        ek_ab = create_edge_key_int(a, b)
        ek_ac = create_edge_key_int(a, c)

        # Case 1: Edge locked
        bc_locked_edge = make_constraints(extA=[ek_ab])
        has_violation, reasons = Nas2Step.check_constraint_violations(ek_ab, bc_locked_edge)
        @test has_violation == true
        @test any(occursin("Edge is locked", r) for r in reasons)

        # Case 2: Endpoint locked (corner node)
        bc_locked_node = make_constraints(corners=[a])
        has_violation2, reasons2 = Nas2Step.check_constraint_violations(ek_ac, bc_locked_node)
        @test has_violation2 == true
        @test any(occursin("Endpoint node1 is locked", r) for r in reasons2)

        # Case 3: No violations
        bc_clear = make_constraints()
        has_violation3, reasons3 = Nas2Step.check_constraint_violations(ek_ab, bc_clear)
        @test has_violation3 == false
        @test isempty(reasons3)
    end
end

@testset "Boundary Constraints - Corner Nodes" begin
    @testset "Corner Node Detection/Locking" begin
        # Build manual constraints with two corner nodes
        p = ckey((0.0, 0.0, 0.0))
        q = ckey((1.0, 1.0, 1.0))
        bc = make_constraints(corners=[p, q])
        @test length(bc.corner_nodes) == 2
        @test bc.locked_nodes == bc.corner_nodes
        # Edges incident to these nodes should trigger endpoint violation
        ek = create_edge_key_int(p, q)
        hv, reasons = Nas2Step.check_constraint_violations(ek, bc)
        @test hv
        @test any(occursin("Endpoint node1 is locked", r) for r in reasons) || any(occursin("Endpoint node2 is locked", r) for r in reasons)
    end
end

@testset "Boundary Constraints - External Edges" begin
    @testset "External Edge Identification (locked union)" begin
        a = ckey((0.0, 0.0, 0.0))
        b = ckey((1.0, 0.0, 0.0))
        e = create_edge_key_int(a, b)
        # Mark as external on A and ensure it appears in locked_edges
        bc = make_constraints(extA=[e])
        @test e in bc.pidA_external_edges
        @test e in bc.locked_edges
    end
end

@testset "Boundary Constraints - Export JSON" begin
    @testset "Export and parse JSON report" begin
        # Small synthetic constraints
        a = ckey((0.0, 0.0, 0.0))
        b = ckey((1.0, 0.0, 0.0))
        e = create_edge_key_int(a, b)
        bc = make_constraints(pidA=10, pidB=20, extA=[e], corners=[a])

        tmpdir = mktempdir()
        outfile = joinpath(tmpdir, "constraints.json")
        path = Nas2Step.export_constraints_json(bc, outfile)
        @test path == outfile
        @test isfile(outfile)

        # Parse and validate a few fields
        data = JSON.parse(read(outfile, String))
        @test data["interface"]["pidA"] == 10
        @test data["interface"]["pidB"] == 20
        @test data["summary"]["corner_nodes"] == 1
        @test data["summary"]["external_edges_A"] == 1
        # Samples should be arrays
        @test haskey(data["details"], "corner_nodes_sample")
        @test data["details"]["corner_nodes_sample"] isa Vector
    end
end


