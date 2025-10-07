"""
    test_edge_classification.jl

Unit tests for edge_classification.jl
Tests point_on_segment, find_hanging_nodes_on_edge, compute_triangle_quality,
find_quad_for_diagonal, and classify_edge_mismatch.
"""

using Test
using Nas2Step
using .Nas2StepTestUtils: ckey

const Z = 0.0

tri(n1, n2, n3, elem_id, c1, c2, c3) = Triangle(n1, n2, n3, elem_id, c1, c2, c3)

function make_topology(; pidA::Int=1, pidB::Int=2, facesA::Vector{Triangle}=Triangle[], facesB::Vector{Triangle}=Triangle[])
    # Build shared keys from provided triangles
    pts = NTuple{3,Float64}[]
    for t in facesA
        push!(pts, t.coord1, t.coord2, t.coord3)
    end
    for t in facesB
        push!(pts, t.coord1, t.coord2, t.coord3)
    end
    shared = Set(ckey.(unique(pts)))
    node_key_to_ids = Dict(k => (1, 1) for k in shared)
    edgesA = Dict{EdgeKey,Vector{Int}}()
    edgesB = Dict{EdgeKey,Vector{Int}}()
    bbox = BoundingBox(collect(shared))
    return InterfaceTopology(
        pidA, pidB,
        shared,
        node_key_to_ids,
        facesA,
        facesB,
        edgesA,
        edgesB,
        Set{EdgeKey}(),
        Set{EdgeKey}(),
        Set{EdgeKey}(),
        bbox,
        length(shared),
        length(facesA),
        length(facesB),
        0,
        0,
        1.0,
    )
end

@testset "Edge Classification - Point on Segment" begin
    a = (0.0, 0.0, Z)
    b = (10.0, 0.0, Z)

    @testset "Interior/Endpoint/Off-segment" begin
        # Interior point
        p_mid = (5.0, 0.0, Z)
        is_on, t, d2 = Nas2Step.point_on_segment(p_mid, a, b)
        @test is_on == true
        @test isapprox(t, 0.5, atol=1e-12)
        @test d2 <= 1e-16

        # Endpoint should not count as interior
    is_on_a, t_a, d2_a = Nas2Step.point_on_segment(a, a, b)
        @test is_on_a == false
        @test isapprox(t_a, 0.0; atol=1e-12)
        @test d2_a <= 1e-16

        # Off-segment point
        p_off = (5.0, 1.0, Z)
    is_on_off, t_off, d2_off = Nas2Step.point_on_segment(p_off, a, b)
        @test is_on_off == false
        @test isapprox(t_off, 0.5; atol=1e-12)  # projection
        @test d2_off > 0.9
    end

    @testset "Degenerate edge" begin
        a2 = (1.0, 1.0, Z)
        b2 = a2  # zero-length
        # Exactly at the point -> within tol
    is_on_deg, t_deg, d2_deg = Nas2Step.point_on_segment(a2, a2, b2)
        @test is_on_deg == true
        @test t_deg == 0.0
        @test d2_deg <= 1e-16
        # Far from the point -> not on segment
        p_far = (2.0, 2.0, Z)
    is_on_far, _, _ = Nas2Step.point_on_segment(p_far, a2, b2)
        @test is_on_far == false
    end
end

@testset "Edge Classification - Triangle Quality" begin
    @testset "Equilateral and Degenerate" begin
        c1 = (0.0, 0.0, Z); c2 = (1.0, 0.0, Z); c3 = (0.5, sqrt(3)/2, Z)
        t_eq = tri(1,2,3,1,c1,c2,c3)
    q_eq = Nas2Step.compute_triangle_quality(t_eq)
        @test isapprox(q_eq, 1.0; atol=1e-6)

        # Degenerate (collinear)
        d1 = (0.0, 0.0, Z); d2p = (1.0, 0.0, Z); d3 = (2.0, 0.0, Z)
        t_deg = tri(1,2,3,2,d1,d2p,d3)
    q_deg = Nas2Step.compute_triangle_quality(t_deg)
    @test q_deg < 1e-6
    end

    @testset "Right Isosceles (45-45-90)" begin
        r1 = (0.0, 0.0, Z); r2 = (1.0, 0.0, Z); r3 = (0.0, 1.0, Z)
        t_rt = tri(1,2,3,3,r1,r2,r3)
    q_rt = Nas2Step.compute_triangle_quality(t_rt)
        @test isapprox(q_rt, 45.0/60.0; atol=1e-6)  # 0.75
    end
end

@testset "Edge Classification - Hanging Nodes" begin
    # Edge from (0,0,0) to (10,0,0)
    a = (0.0, 0.0, Z); b = (10.0, 0.0, Z)
    edge = EdgeKey(ckey(a), ckey(b))
    nodes = Set([ckey(a), ckey(b), ckey((2.0,0.0,Z)), ckey((7.0,0.0,Z)), ckey((5.0,0.1,Z))])
    hanging = Nas2Step.find_hanging_nodes_on_edge(edge, nodes)
    @test length(hanging) == 2
    @test hanging[1] == ckey((2.0,0.0,Z))
    @test hanging[2] == ckey((7.0,0.0,Z))
end

@testset "Edge Classification - Quad Finding" begin
    # Square ABCD: A(0,0), B(1,0), C(1,1), D(0,1)
    A = (0.0,0.0,Z); B = (1.0,0.0,Z); C = (1.0,1.0,Z); D = (0.0,1.0,Z)
    # Source side (A) has diagonal BD, target side (B) has diagonal AC
    t1_source = tri(1,2,4,10,B,A,D)   # BAD (source has BD)
    t2_source = tri(2,3,4,11,B,C,D)   # BCD (source has BD)
    t1_target = tri(1,2,3,1,A,B,C)    # ABC (target has AC)
    t2_target = tri(1,3,4,2,A,C,D)    # ACD (target has AC)
    topo = make_topology(facesA=[t1_source,t2_source], facesB=[t1_target,t2_target])
    # Edge BD exists in A, missing in B
    ek_bd = EdgeKey(ckey(B), ckey(D))
    # present_in=:A (edge exists in A), target_pid=pidB (needs to be added to B)
    qverts, to_replace = Nas2Step.find_quad_for_diagonal(ek_bd, topo, :A, topo.pidB)
    @test length(qverts) == 4
    @test Set(qverts) == Set([ckey(A), ckey(B), ckey(C), ckey(D)])
    @test Set(to_replace) == Set([1,2])
end

@testset "Edge Classification - Mismatch Classification" begin
    # Shared: same square vertices
    A = (0.0,0.0,Z); B = (1.0,0.0,Z); C = (1.0,1.0,Z); D = (0.0,1.0,Z)
    # Target side faces include edge AB
    tri_abx = tri(1,2,3,1,A,B,C)
    topo_base = make_topology(facesB=[tri_abx])

    @testset "T-Junction (single hanging node)" begin
        # Edge AB, with a midpoint present in shared nodes
        ek_ab = EdgeKey(ckey(A), ckey(B))
        # Inject one hanging node into shared set
        push!(topo_base.shared_node_keys, ckey((0.5,0.0,Z)))
    m = Nas2Step.classify_edge_mismatch(ek_ab, topo_base, :A)
    @test m.mismatch_type == Nas2Step.T_JUNCTION
        @test m.present_in == :B_only
        @test m.should_be_in == topo_base.pidB
        @test length(m.hanging_nodes) == 1
        @test !isempty(m.affected_triangles)
        @test m.repair_feasible == true
    end

    @testset "Refinement (multiple hanging nodes)" begin
        # New topology to avoid shared mutation
        tri_abx2 = tri(1,2,3,1,A,B,C)
        topo_ref = make_topology(facesB=[tri_abx2])
        push!(topo_ref.shared_node_keys, ckey((0.25,0.0,Z)))
        push!(topo_ref.shared_node_keys, ckey((0.75,0.0,Z)))
        ek_ab = EdgeKey(ckey(A), ckey(B))
    m = Nas2Step.classify_edge_mismatch(ek_ab, topo_ref, :A)
    @test m.mismatch_type == Nas2Step.REFINEMENT
        @test length(m.hanging_nodes) == 2
        @test !isempty(m.affected_triangles)
        @test m.repair_feasible == true
    end

    @testset "Diagonal (no hanging nodes)" begin
        # Source side (A) has edge BD, target side (B) uses diagonal AC
        # Source mesh has BD diagonal
        t1_source = tri(1,2,4,10,B,A,D)   # BAD (source has BD)
        t2_source = tri(2,3,4,11,B,C,D)   # BCD (source has BD)
        # Target mesh has AC diagonal
        t1_target = tri(1,2,3,1,A,B,C)    # ABC (target has AC)
        t2_target = tri(1,3,4,2,A,C,D)    # ACD (target has AC)
        topo_diag = make_topology(facesA=[t1_source, t2_source], facesB=[t1_target, t2_target])
        # Edge BD exists in A, missing in B
        ek_bd = EdgeKey(ckey(B), ckey(D))
        # Classify as edge present in A, should be added to B
	m = Nas2Step.classify_edge_mismatch(ek_bd, topo_diag, :A)
	@test m.mismatch_type == Nas2Step.DIAGONAL
        @test length(m.hanging_nodes) == 0
        # For diagonal case, affected triangles list (sharing BD in B) is empty
        @test isempty(m.affected_triangles)
        # But the function should locate the quad and triangles to replace in B
        @test length(m.quad_vertices) == 4
        @test Set(m.triangles_to_replace) == Set([1,2])
        # Feasibility false due to no directly affected triangles on target side
        @test m.repair_feasible == false
    end
end


