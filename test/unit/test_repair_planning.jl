"""
    test_repair_planning.jl

Unit tests for repair_planning.jl
Tests compute_edge_density, determine_dominant_side, plan_quad_retriangulation,
plan_triangle_split, generate_edge_insertion_plan, and generate_repair_plan.
"""

using Test
using Nas2Step
using .Nas2StepTestUtils: ckey

# --- Local helpers for tiny synthetic fixtures ---------------------------------

const P1 = (0.0, 0.0, 0.0)
const P2 = (1.0, 0.0, 0.0)
const P3 = (0.0, 1.0, 0.0)
const P4 = (1.0, 1.0, 0.0)


make_triangle(n1::Int, n2::Int, n3::Int, c1, c2, c3; eid::Int=1) = Triangle(n1, n2, n3, eid, c1, c2, c3)

function build_edges_map(tris::Vector{Triangle})
    m = Dict{EdgeKey, Vector{Int}}()
    for (idx, t) in enumerate(tris)
        k1, k2, k3 = ckey(t.coord1), ckey(t.coord2), ckey(t.coord3)
        for (a, b) in ((k1, k2), (k1, k3), (k2, k3))
            ek = EdgeKey(a, b)
            push!(get!(m, ek, Int[]), idx)
        end
    end
    m
end

 

@testset "Feasibility Assessment - assess_interface_feasibility" begin
    # Helper to craft synthetic EdgeInsertionPlan quickly
    make_plan(; qa::Bool=true, vc::Bool=false, reasons::Vector{String}=String[],
              min_b::Float64=0.5, min_a::Float64=0.5) = begin
        EdgeInsertionPlan(
            1,                                 # target_triangle
            EdgeKey(P1, P2),                   # insert_edge
            :quad_retriangulation,             # split_type
            NTuple{3,Float64}[],               # new_nodes
            Int[],                             # existing_nodes
            NTuple{9,Float64}[],               # old_triangles
            NTuple{9,Float64}[],               # replacement_triangles
            min_b,                             # min_angle_before (normalized)
            min_a,                             # min_angle_after (normalized)
            qa,                                # quality_acceptable
            Int[],                             # depends_on
            vc,                                # violates_constraints
            reasons,                           # constraint_violations
            qa && !vc                          # is_feasible (per-plan)
        )
    end

    thresholds = Nas2Step.default_thresholds()

    @testset "All feasible meets ratio" begin
        plans = [make_plan() for _ in 1:5]
    a = Nas2Step.assess_interface_feasibility(plans, 5, thresholds)
        @test a.is_feasible == true
        @test length(a.feasible_plans) == 5
        @test length(a.infeasible_plans) == 0
        @test isempty(a.issues)
        @test a.qstats.min_before ≈ 0.5
        @test a.qstats.min_after ≈ 0.5
    end

    @testset "Poor quality counted in issues" begin
        plans = [make_plan() for _ in 1:3]
        append!(plans, [make_plan(qa=false, min_b=0.1, min_a=0.1) for _ in 1:2])
    a = Nas2Step.assess_interface_feasibility(plans, 5, thresholds)
        @test length(a.feasible_plans) == 3
        @test length(a.infeasible_plans) == 2
    @test any(occursin("fail quality", issue) for issue in a.issues)
    end

    @testset "Constraint reasons aggregated" begin
        plans = [
            make_plan(vc=true, reasons=["Edge is locked (external or feature)"]),
            make_plan(vc=true, reasons=["Edge is locked (external or feature)", "Endpoint node1 is locked (corner node)"]),
            make_plan(vc=true, reasons=["Endpoint node2 is locked (corner node)"]),
            make_plan(),
        ]
    a = Nas2Step.assess_interface_feasibility(plans, 4, thresholds)
        @test a.is_feasible == false
        @test a.reason_counts["Edge is locked (external or feature)"] == 2
        @test a.reason_counts["Endpoint node1 is locked (corner node)"] == 1
        @test a.reason_counts["Endpoint node2 is locked (corner node)"] == 1
    @test any(occursin("blocked by constraints", issue) for issue in a.issues)
    end

    @testset "Ratio gating" begin
        plans = vcat([make_plan() for _ in 1:7], [make_plan(qa=false) for _ in 1:3])
        # With total_needed = 10 and default ratio 0.8, need >8 => 7 is not enough
    a = Nas2Step.assess_interface_feasibility(plans, 10, thresholds)
        @test a.is_feasible == false
    @test any(occursin("Feasible plans 7/10", issue) for issue in a.issues)
        # If we drop ratio to 0.5, it should pass
    a2 = Nas2Step.assess_interface_feasibility(plans, 10, thresholds; ratio=0.5)
        @test a2.is_feasible == true
    end

    @testset "Empty plans" begin
        a = Nas2Step.assess_interface_feasibility(EdgeInsertionPlan[], 0, thresholds)
        @test a.is_feasible == true
        @test length(a.feasible_plans) == 0
        @test length(a.infeasible_plans) == 0
        @test isempty(a.issues)
    end

    @testset "Ratio boundary is strict (>)" begin
    # Exactly at boundary should be false because gate uses '>' not '>='
        plans = vcat([make_plan() for _ in 1:8], [make_plan(qa=false) for _ in 1:2])
    a = Nas2Step.assess_interface_feasibility(plans, 10, thresholds)
        @test a.is_feasible == false
        # One more feasible crosses the boundary
        plans2 = vcat([make_plan() for _ in 1:9], [make_plan(qa=false)])
    a2 = Nas2Step.assess_interface_feasibility(plans2, 10, thresholds)
        @test a2.is_feasible == true
    end

    @testset "Zero needed but some feasible" begin
        # When total_needed=0, any feasible plans should pass ('>' 0)
        plans = [make_plan() for _ in 1:2]
    a = Nas2Step.assess_interface_feasibility(plans, 0, thresholds)
        @test a.is_feasible == true
        @test isempty(a.issues)
    end

    @testset "Combined issues with exact counts" begin
        # 3 constraint-violating plans and 2 poor-quality plans, total 7, feasible 2
        c1 = make_plan(vc=true, reasons=["Edge is locked (external or feature)"])
        c2 = make_plan(vc=true, reasons=["Endpoint node1 is locked (corner node)"])
        c3 = make_plan(vc=true, reasons=["Edge is locked (external or feature)", "Endpoint node2 is locked (corner node)"])
        q1 = make_plan(qa=false, min_b=0.05, min_a=0.05)
        q2 = make_plan(qa=false, min_b=0.02, min_a=0.02)
        ok1 = make_plan()
        ok2 = make_plan()
        plans = [c1, c2, c3, q1, q2, ok1, ok2]
    a = Nas2Step.assess_interface_feasibility(plans, 7, thresholds)
        # Expect not feasible, with all three issue categories present
        @test a.is_feasible == false
    @test any(occursin("Feasible plans 2/7", i) for i in a.issues)
    @test any(occursin("3 plans blocked by constraints", i) for i in a.issues)
    @test any(occursin("2 plans fail quality", i) for i in a.issues)
        # Reason counts aggregated correctly
        @test a.reason_counts["Edge is locked (external or feature)"] == 2
        @test a.reason_counts["Endpoint node1 is locked (corner node)"] == 1
        @test a.reason_counts["Endpoint node2 is locked (corner node)"] == 1
    end

    @testset "QStats and threshold normalization" begin
        # Use diverse min_before/after and a custom threshold to test normalization
        plans = [
            make_plan(min_b=0.6, min_a=0.55),
            make_plan(min_b=0.3, min_a=0.25),
            make_plan(min_b=0.8, min_a=0.7),
        ]
    th = QualityThresholds(12.0, 180.0, 10.0) # 12° => normalized 0.2
        a = Nas2Step.assess_interface_feasibility(plans, 3, th)
        @test isapprox(a.qstats.min_before, 0.3; atol=1e-12)
        @test isapprox(a.qstats.min_after, 0.25; atol=1e-12)
        @test isapprox(a.qstats.min_threshold_norm, 12.0/60.0; atol=1e-12)
    end
end

function bbox_from_tris(tris::Vector{Triangle})
    pts = NTuple{3,Float64}[]
    for t in tris
        push!(pts, t.coord1, t.coord2, t.coord3)
    end
    BoundingBox(pts)
end

function make_topology(fA::Vector{Triangle}, fB::Vector{Triangle}; pidA::Int=1, pidB::Int=2)
    edgesA = build_edges_map(fA)
    edgesB = build_edges_map(fB)
    shared = Set{NTuple{3,Float64}}(ckey.(vcat([t.coord1 for t in fA]..., [t.coord2 for t in fA]..., [t.coord3 for t in fA]...)))
    node_map = Dict{NTuple{3,Float64}, Tuple{Int,Int}}()
    bbox = bbox_from_tris(!isempty(fA) ? fA : fB)
    InterfaceTopology(
        pidA, pidB,
        shared,
        node_map,
        fA,
        fB,
        edgesA,
        edgesB,
        Set{EdgeKey}(),
        Set{EdgeKey}(),
        Set{EdgeKey}(),
        bbox,
        length(shared),
        length(fA),
        length(fB),
        length(edgesA),
        length(edgesB),
        0.0,
    )
end

make_constraints(pidA::Int, pidB::Int; extA=Set{EdgeKey}(), extB=Set{EdgeKey}(),
                 corners=Set{NTuple{3,Float64}}(), features=Set{EdgeKey}(),
                 locked_nodes=Set{NTuple{3,Float64}}(), locked_edges=Set{EdgeKey}()) =
    BoundaryConstraints(pidA, pidB, extA, extB, corners, features, locked_nodes, union(extA, extB, features, locked_edges), length(extA), length(extB), length(corners), length(union(extA, extB, features, locked_edges)))

 

@testset "Repair Planning - Edge Density" begin
    @testset "Compute Edge Density" begin
        # A: single right triangle (area 0.5), 3 unique edges => density 6.0
        tA = [make_triangle(1, 2, 3, P1, P2, P3; eid=1)]
        # B: square split into two triangles (area 1.0), 5 unique edges => density 5.0
        tB = [
            make_triangle(1, 2, 4, P1, P2, P4; eid=1),
            make_triangle(1, 4, 3, P1, P4, P3; eid=2),
        ]
        topo = make_topology(tA, tB)
        @test isapprox(Nas2Step.compute_edge_density(topo, :A), 6.0; atol=1e-12)
        @test isapprox(Nas2Step.compute_edge_density(topo, :B), 5.0; atol=1e-12)
    end
end

@testset "Repair Planning - Dominant Side" begin
    @testset "Determine Dominant Side by density" begin
        # Make A much denser than B: shrink A's area
        small = (0.2, 0.0, 0.0)
        tA = [make_triangle(1, 2, 3, P1, small, P3; eid=1)]  # area 0.1, 3 edges => density 30.0
        tB = [make_triangle(1, 2, 4, P1, P2, P4; eid=1), make_triangle(1, 4, 3, P1, P4, P3; eid=2)] # density 5.0
        topo = make_topology(tA, tB)
        cons = make_constraints(1, 2)
        @test Nas2Step.determine_dominant_side(topo, cons) == :A

        # Swap to make B denser
        tA2 = [make_triangle(1, 2, 4, P1, P2, P4; eid=1), make_triangle(1, 4, 3, P1, P4, P3; eid=2)]
        tB2 = [make_triangle(1, 2, 3, P1, small, P3; eid=1)]
        topo2 = make_topology(tA2, tB2)
        @test Nas2Step.determine_dominant_side(topo2, cons) == :B
    end

    @testset "Determine Dominant Side by constraints" begin
        # Equal-ish densities, prefer side with fewer external constraints
        tA = [make_triangle(1, 2, 3, P1, P2, P3; eid=1)]  # density 6.0
        tB = [make_triangle(1, 2, 3, P1, P2, P3; eid=1)]  # density 6.0
        topo = make_topology(tA, tB)
        # Make A have fewer locked (external) edges than B
        # Build some dummy edge keys to populate external sets
        eAB = collect(keys(build_edges_map(tA)))
        extA = Set{EdgeKey}([eAB[1]])                 # 1/3 ~ 0.333
        extB = Set{EdgeKey}([eAB[1], eAB[min(end,2)]]) # 2/3 ~ 0.666
        cons = make_constraints(1, 2; extA=extA, extB=extB)
        # A fewer constraints -> prefer modifying A -> dominant should be :B
        @test Nas2Step.determine_dominant_side(topo, cons) == :B
    end

    @testset "Determine Dominant Side tie-breaker by size" begin
        # Equal densities and constraints, tie-break by smaller side (fewer faces)
        tA = [make_triangle(1, 2, 3, P1, P2, P3; eid=1)]
        tB = [make_triangle(1, 2, 4, P1, P2, P4; eid=1), make_triangle(1, 4, 3, P1, P4, P3; eid=2)]
        topo = make_topology(tA, tB)
        cons = make_constraints(1, 2)
        # A is smaller -> prefer modifying A -> dominant should be :B
        @test Nas2Step.determine_dominant_side(topo, cons) == :B
    end
end

@testset "Repair Planning - Quad Retriangulation" begin
    @testset "Plan Quad Retriangulation" begin
        quad = [P1, P2, P3, P4]
        diag = EdgeKey(P1, P4)  # use the P1-P4 diagonal
        plan = Nas2Step.plan_quad_retriangulation(quad, diag)
        @test plan !== nothing
        @test plan.type == :quad_retriangulation
        @test isempty(plan.new_nodes)
        @test length(plan.new_triangles) == 2
        # Both triangles must include the diagonal endpoints
        for tri in plan.new_triangles
            @test tri[1:3] == P1
            @test tri[4:6] == P4
        end
    # Third vertices should be the two remaining corners (P2 and P3) in any order
    t1 = plan.new_triangles[1]
    t2 = plan.new_triangles[2]
    v1 = (t1[7], t1[8], t1[9])
    v2 = (t2[7], t2[8], t2[9])
    third_vertices = Set([v1, v2])
    @test third_vertices == Set([P2, P3])

        # Invalid input: not a quad
        @test Nas2Step.plan_quad_retriangulation([P1, P2, P3], diag) === nothing
    end
end

@testset "Repair Planning - Triangle Split" begin
    @testset "Plan Triangle Split" begin
        tri = make_triangle(1, 2, 3, P1, P2, P3; eid=10)
        edge = EdgeKey(P1, P2)
        plan = Nas2Step.plan_triangle_split(tri, edge)
        @test plan !== nothing
        @test plan.type == :bisect
        @test isempty(plan.new_nodes)
        @test length(plan.new_triangles) == 1
        t = plan.new_triangles[1]
        @test t[1:3] == P1
        @test t[4:6] == P2
        @test t[7:9] == P3

        # Edge not in triangle -> nothing
        edge2 = EdgeKey(P1, P4)
        @test Nas2Step.plan_triangle_split(tri, edge2) === nothing
    end
end

@testset "Repair Planning - Edge Insertion" begin
    @testset "Generate Edge Insertion Plan (diagonal success)" begin
        # Target side B has two triangles of a square (diagonal P2-P3); we want to insert P1-P4
        tA = [make_triangle(1, 2, 3, P1, P2, P3; eid=1)]
        tB = [
            make_triangle(2, 3, 4, P2, P3, P4; eid=1),
            make_triangle(2, 4, 1, P2, P4, P1; eid=2),
        ]
        topo = make_topology(tA, tB)
        diag_edge = EdgeKey(P1, P4)
        mismatch = EdgeMismatch(
            diag_edge,
            Nas2Step.DIAGONAL,
            :A_only,
            topo.pidB,
            NTuple{3,Float64}[],  # hanging
            Int[],                 # affected (not used in DIAGONAL path)
            [P1, P2, P3, P4],      # quad vertices
            [1, 2],                # triangles to replace in B
            0.4,                   # complexity
            1.0,                   # min_affected_triangle_quality (unused here)
            true                   # repair_feasible
        )
        cons = make_constraints(1, 2)
    plan = Nas2Step.generate_edge_insertion_plan(mismatch, topo, cons, Nas2Step.default_thresholds())
        @test plan isa EdgeInsertionPlan
        @test plan.split_type == :quad_retriangulation
        @test plan.target_triangle == 1
        @test plan.violates_constraints == false
        @test plan.quality_acceptable == true
        @test plan.is_feasible == true
        @test length(plan.old_triangles) == 2
        @test length(plan.replacement_triangles) == 2
    end

    @testset "Generate Edge Insertion Plan (diagonal missing info)" begin
        tB = [make_triangle(2, 3, 4, P2, P3, P4; eid=1), make_triangle(2, 4, 1, P2, P4, P1; eid=2)]
        topo = make_topology(Triangle[], tB)
        diag_edge = EdgeKey(P1, P4)
        mismatch = EdgeMismatch(diag_edge, Nas2Step.DIAGONAL, :A_only, topo.pidB,
                                NTuple{3,Float64}[], Int[], NTuple{3,Float64}[], Int[], 0.3, 1.0, false)
        cons = make_constraints(1, 2)
    plan = Nas2Step.generate_edge_insertion_plan(mismatch, topo, cons, Nas2Step.default_thresholds())
        @test plan isa EdgeInsertionPlan
        @test plan.split_type == :failed
        # Implementation sets violates_constraints=true in this branch
        @test plan.violates_constraints == true
        @test plan.is_feasible == false
    end

    @testset "Generate Edge Insertion Plan (constraint violation)" begin
        tB = [make_triangle(2, 3, 4, P2, P3, P4; eid=1), make_triangle(2, 4, 1, P2, P4, P1; eid=2)]
        topo = make_topology(Triangle[], tB)
        diag_edge = EdgeKey(P1, P4)
        mismatch = EdgeMismatch(diag_edge, Nas2Step.DIAGONAL, :A_only, topo.pidB,
                                NTuple{3,Float64}[], Int[], [P1, P2, P3, P4], [1,2], 0.3, 1.0, true)
    # Lock an endpoint node to guarantee a violation regardless of edge set behavior
    cons = make_constraints(1, 2; extA=Set{EdgeKey}(), extB=Set{EdgeKey}(), locked_nodes=Set([P1]))
    plan = Nas2Step.generate_edge_insertion_plan(mismatch, topo, cons, Nas2Step.default_thresholds())
        @test plan.violates_constraints == true
        @test plan.is_feasible == false
    end

    @testset "Generate Edge Insertion Plan (T-junction split)" begin
        # Target side A has the triangle; we need to insert edge P1-P2
        tri = make_triangle(1, 2, 3, P1, P2, P3; eid=10)
        topo = make_topology([tri], Triangle[])
        tedge = EdgeKey(P1, P2)
        mismatch = EdgeMismatch(tedge, Nas2Step.T_JUNCTION, :B_only, topo.pidA,
                                NTuple{3,Float64}[], [1], NTuple{3,Float64}[], Int[], 0.2, 1.0, true)
        cons = make_constraints(1, 2)
    plan = Nas2Step.generate_edge_insertion_plan(mismatch, topo, cons, Nas2Step.default_thresholds())
        @test plan.split_type == :bisect
        @test plan.violates_constraints == false
        @test plan.is_feasible == true
        @test length(plan.old_triangles) == 1
        @test length(plan.replacement_triangles) == 1

        # No affected triangles => failed plan
        mismatch2 = EdgeMismatch(tedge, Nas2Step.T_JUNCTION, :B_only, topo.pidA,
                                 NTuple{3,Float64}[], Int[], NTuple{3,Float64}[], Int[], 0.2, 1.0, false)
        plan2 = Nas2Step.generate_edge_insertion_plan(mismatch2, topo, cons, Nas2Step.default_thresholds())
        @test plan2 isa EdgeInsertionPlan
        @test plan2.split_type == :failed
        @test plan2.is_feasible == false
        @test !isempty(plan2.constraint_violations)
    end
end

@testset "Repair Planning - Full Repair Plan" begin
    @testset "Generate Repair Plan and direction selection" begin
        # Topology where A is much denser -> plan should subdivide B
        small = (0.2, 0.0, 0.0)
        tA = [make_triangle(1, 2, 3, P1, small, P3; eid=1)]  # very dense
        # B has the square with diagonal P2-P3; we want to add P1-P4
        tB = [make_triangle(2, 3, 4, P2, P3, P4; eid=1), make_triangle(2, 4, 1, P2, P4, P1; eid=2)]
        topo = make_topology(tA, tB)
        cons = make_constraints(1, 2)
        # Build classification with one DIAGONAL mismatch missing in B (edges only in A)
        diag = EdgeKey(P1, P4)
        mismatchB = EdgeMismatch(diag, Nas2Step.DIAGONAL, :A_only, topo.pidB,
                                  NTuple{3,Float64}[], Int[], [P1, P2, P3, P4], [1,2], 0.4, 1.0, true)
        cls = InterfaceClassification(topo, EdgeMismatch[], [mismatchB], 0, 1, 0, 0, 1, 0, 0.4)

    plan = Nas2Step.generate_repair_plan(topo, cls, cons; thresholds=Nas2Step.default_thresholds())
        @test plan.repair_direction == :subdivide_B
        @test length(plan.insertion_sequence) == 1
        @test plan.total_edges_to_insert == 1
        @test plan.is_feasible == true
        @test plan.predicted_min_quality >= 0.0
    end

    @testset "Generate Repair Plan direction override" begin
        # Equal densities; initial dominant may pick :A -> :subdivide_B
        # but provide mismatches only in A so it should flip to :subdivide_A
        tA = [make_triangle(1, 2, 3, P1, P2, P3; eid=1)]
        tB = [make_triangle(1, 2, 3, P1, P2, P3; eid=1)]
        topo = make_topology(tA, tB)
        cons = make_constraints(1, 2)
        # Mismatch missing in A (present only in B)
        m_edge = EdgeKey(P1, P2)
        mismatchA = EdgeMismatch(m_edge, Nas2Step.T_JUNCTION, :B_only, topo.pidA,
                                  NTuple{3,Float64}[], [1], NTuple{3,Float64}[], Int[], 0.2, 1.0, true)
        cls = InterfaceClassification(topo, [mismatchA], EdgeMismatch[], 1, 0, 0, 0, 1, 0, 0.2)

    plan = Nas2Step.generate_repair_plan(topo, cls, cons; thresholds=Nas2Step.default_thresholds())
        @test plan.repair_direction == :subdivide_A
        @test length(plan.insertion_sequence) == 1
    end
end


