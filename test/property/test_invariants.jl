"""
    test_invariants.jl

Property-based and invariant tests for core data structures and workspace.
Covers:
- EdgeKey commutativity and hashing
- Triangle quality bounds and degeneracy
- BoundingBox containment
- Workspace modification count consistency and rollback restores state
- InterfaceTopology simple invariant on shared edges bound
"""

using Test
using Nas2Step
using .Nas2StepTestUtils: random_point_3d, ckey, capture_workspace_state, workspace_state_equals, create_minimal_workspace

const Z = 0.0

# Local helper to build a simple topology without Gmsh using existing structs
tri_prop(n1, n2, n3, elem_id, c1, c2, c3) = Triangle(n1, n2, n3, elem_id, c1, c2, c3)

function make_simple_topology()
    # Square ABCD with two triangles per side
    A = (0.0,0.0,Z); B = (1.0,0.0,Z); C = (1.0,1.0,Z); D = (0.0,1.0,Z)
    tA1 = tri_prop(1,2,3, 1, A,B,C)
    tA2 = tri_prop(1,3,4, 2, A,C,D)
    # Mirror on B side
    tB1 = tri_prop(1,2,3, 1, A,B,C)
    tB2 = tri_prop(1,3,4, 2, A,C,D)

    pts = NTuple{3,Float64}[]
    for t in (tA1,tA2,tB1,tB2); push!(pts, t.coord1, t.coord2, t.coord3); end
    shared = Set(ckey.(unique(pts)))
    node_key_to_ids = Dict(k => (1,1) for k in shared)
    bbox = BoundingBox(collect(shared))

    # Minimal edge maps with shared edges equal to all unique triangle edges
    edgesA = Dict{EdgeKey,Vector{Int}}()
    edgesB = Dict{EdgeKey,Vector{Int}}()
    for (i, t) in enumerate((tA1,tA2))
        k1, k2, k3 = ckey(t.coord1), ckey(t.coord2), ckey(t.coord3)
        for (a,b) in ((k1,k2),(k2,k3),(k1,k3))
            ek = EdgeKey(a,b)
            push!(get!(edgesA, ek, Int[]), i)
        end
    end
    for (i, t) in enumerate((tB1,tB2))
        k1, k2, k3 = ckey(t.coord1), ckey(t.coord2), ckey(t.coord3)
        for (a,b) in ((k1,k2),(k2,k3),(k1,k3))
            ek = EdgeKey(a,b)
            push!(get!(edgesB, ek, Int[]), i)
        end
    end

    setA = Set(keys(edgesA)); setB = Set(keys(edgesB)); sharedE = intersect(setA, setB)

    return InterfaceTopology(
        1, 2,
        shared,
        node_key_to_ids,
        [tA1, tA2],
        [tB1, tB2],
        edgesA,
        edgesB,
        setdiff(setA, setB),
        setdiff(setB, setA),
        sharedE,
        bbox,
        length(shared),
        2, 2,
        length(setA), length(setB),
        length(union(setA, setB)) > 0 ? length(sharedE)/length(union(setA,setB)) : 1.0,
    )
end

@testset "Property: EdgeKey Commutativity" begin
    for _ in 1:50
        p1 = ckey(random_point_3d())
        p2 = ckey(random_point_3d())
        ek1 = EdgeKey(p1, p2)
        ek2 = EdgeKey(p2, p1)
        @test ek1 == ek2
        @test hash(ek1) == hash(ek2)
    end
end

@testset "Property: Triangle Quality Bounds" begin
    for _ in 1:50
        # Build non-degenerate random triangle by rejection sampling
        c1 = random_point_3d(); c2 = random_point_3d(); c3 = random_point_3d()
    t = tri_prop(1,2,3,1,c1,c2,c3)
        q = Nas2Step.compute_triangle_quality(t)
        @test 0.0 <= q <= 1.0
    end
    # Degenerate triangle should yield near-zero quality
    d1 = (0.0,0.0,Z); d2 = (1.0,0.0,Z); d3 = (2.0,0.0,Z)
    td = tri_prop(1,2,3,1,d1,d2,d3)
    qd = Nas2Step.compute_triangle_quality(td)
    @test qd <= 1e-6
end

@testset "Property: BoundingBox Contains Points" begin
    points = [random_point_3d() for _ in 1:100]
    bb = BoundingBox(points)
    for p in points
        @test p[1] >= bb.min_corner[1]
        @test p[1] <= bb.max_corner[1]
        @test p[2] >= bb.min_corner[2]
        @test p[2] <= bb.max_corner[2]
        @test p[3] >= bb.min_corner[3]
        @test p[3] <= bb.max_corner[3]
    end
    # Empty points: bbox is zeros
    bb_empty = BoundingBox(NTuple{3,Float64}[])
    @test bb_empty.min_corner == (0.0,0.0,0.0)
    @test bb_empty.max_corner == (0.0,0.0,0.0)
end

@testset "Property: Workspace Invariants" begin
    # Use in-memory minimal workspace; operations require Phase 3 include in runtests
    ws = create_minimal_workspace()
    # Begin/end transaction and ensure counters match modifications
    Main.begin_transaction!(ws)

    # Add a node
    nid = Main.add_node!(ws, (1.0, 2.0, 3.0))
    @test nid > 0
    # Add a face using existing nodes (pick first PID and its first two nodes plus new)
    pid = first(keys(ws.working_faces))
    existing_face = ws.working_faces[pid][1]
    Main.add_face!(ws, pid, [existing_face[1], existing_face[2], nid])
    # Delete a face
    Main.delete_face!(ws, pid, 1)

    expected = ws.faces_added + ws.faces_deleted + ws.nodes_added
    @test length(ws.modifications) == expected

    # Capture and rollback
    snap = capture_workspace_state(ws)
    Main.rollback_transaction!(ws)

    # After rollback, transaction becomes inactive and state matches checkpoint (initial)
    @test ws.transaction_active == false
    @test workspace_state_equals(ws, snap) == false  # snap was taken mid-transaction

    # Now check that rolling back restores to checkpoint
    # Start again and record state before modifications
    Main.begin_transaction!(ws)
    base = capture_workspace_state(ws)
    nid2 = Main.add_node!(ws, (9.0,9.0,9.0))
    Main.add_face!(ws, pid, [existing_face[1], existing_face[2], nid2])
    Main.delete_face!(ws, pid, 1)
    Main.rollback_transaction!(ws)
    @test workspace_state_equals(ws, base)
end

@testset "Property: Topology Shared Edge Bound" begin
    topo = make_simple_topology()
    shared = length(topo.edges_shared)
    edges_A = length(topo.edges_A)
    edges_B = length(topo.edges_B)
    @test shared <= min(edges_A, edges_B)
end

@testset "Property: Topology Face Incidence Sum" begin
    # Helper: given tetrahedra (as 4-tuples of node ids), compute face incidence map
    faces_of_tet(t) = (
        (t[1], t[2], t[3]),
        (t[1], t[2], t[4]),
        (t[1], t[3], t[4]),
        (t[2], t[3], t[4])
    )
    # Normalize a face as a sorted 3-tuple key
    face_key(f::NTuple{3,Int}) = Tuple(sort(collect(f)))
    function face_incidence(tets::Vector{NTuple{4,Int}})
        inc = Dict{NTuple{3,Int},Int}()
        for tet in tets
            for f in faces_of_tet(tet)
                fs = face_key(f)
                inc[fs] = get(inc, fs, 0) + 1
            end
        end
        return inc
    end

    # Case 1: two disjoint tetrahedra
    tets1 = NTuple{4,Int}[(1,2,3,4), (5,6,7,8)]
    inc1 = face_incidence(tets1)
    @test sum(values(inc1)) == 4 * length(tets1)

    # Case 2: two tetrahedra sharing one face (1,2,3)
    tets2 = NTuple{4,Int}[(1,2,3,4), (1,2,3,5)]
    inc2 = face_incidence(tets2)
    @test inc2[face_key((1,2,3))] == 2  # shared face counted twice
    @test sum(values(inc2)) == 4 * length(tets2)
end
