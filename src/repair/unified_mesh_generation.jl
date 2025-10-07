# Unified third-mesh generation (Phase 5 MVP)
#
# Minimal implementation that builds a unified interface mesh by copying
# triangles adjacent to selected edges from the chosen side (:use_A or :use_B),
# avoiding duplicates and producing basic node mappings and topology.

"""
    generate_unified_interface_mesh(topology, symmetric_mismatches, constraints; thresholds=default_thresholds(), verbose=true)

Build a unified interface mesh from per-edge symmetric strategies.

MVP behavior:
- For edges marked :use_A, copy the triangles in A that use that edge
- For edges marked :use_B, copy the triangles in B that use that edge
- Avoid duplicate triangles; skip conflicts conservatively
- Build basic node mappings, edge topology, quality metrics
- Mark compatibility true if we produced any triangles; otherwise false
"""
function generate_unified_interface_mesh(
    topology::InterfaceTopology,
    symmetric_mismatches::Vector{SymmetricEdgeMismatch},
    constraints::BoundaryConstraints;
    thresholds::QualityThresholds = default_thresholds(),
    verbose::Bool = true
)::UnifiedInterfaceMesh

    verbose && println("\n[Phase 5] Generating unified interface mesh (MVP)...")

    # Helper: coordinate key with 4-digit rounding (must match EdgeKey policy)
    ckey(p::NTuple{3,Float64}) = (round(p[1]; digits=4), round(p[2]; digits=4), round(p[3]; digits=4))

    # Accumulators
    unified_tris = Triangle[]
    provenance = Symbol[]

    # Duplicate filter: use sorted rounded coords as canonical signature
    function tri_signature(tri::Triangle)
        pts = [ckey(tri.coord1), ckey(tri.coord2), ckey(tri.coord3)]
        sort!(pts)
        return (pts[1], pts[2], pts[3])
    end
    seen = Set{Tuple{NTuple{3,Float64},NTuple{3,Float64},NTuple{3,Float64}}}()
    # Edge incidence map for conflict detection (rounded pairs)
    edge_incidence = Dict{Tuple{NTuple{3,Float64},NTuple{3,Float64}}, Int}()

    # Convenience: accessors
    faces_A = topology.faces_A
    faces_B = topology.faces_B
    edges_A = topology.edges_A
    edges_B = topology.edges_B

    # Sort mismatches by priority (desc), favor use_A/use_B, ignore :skip for MVP
    candidates = filter(m -> m.repair_strategy == :use_A || m.repair_strategy == :use_B, symmetric_mismatches)
    sort!(candidates, by = m -> m.repair_priority, rev = true)

    # Utility to push triangle if not duplicate or conflicting
    function try_add_triangle!(tri::Triangle, src::Symbol)
        sig = tri_signature(tri)
        if sig in seen
            return false
        end
        # Conservative conflict check: avoid creating non-manifold edges (>2)
        edges = ((ckey(tri.coord1), ckey(tri.coord2)), (ckey(tri.coord2), ckey(tri.coord3)), (ckey(tri.coord3), ckey(tri.coord1)))
        for (a, b) in edges
            e = a <= b ? (a, b) : (b, a)
            if get(edge_incidence, e, 0) >= 2
                return false
            end
        end
        push!(seen, sig)
        push!(unified_tris, tri)
        push!(provenance, src)
        # Update incidences
        for (a, b) in edges
            e = a <= b ? (a, b) : (b, a)
            edge_incidence[e] = get(edge_incidence, e, 0) + 1
        end
        return true
    end

    # Fetch triangles using an edge from a given side via edge maps
    function triangles_for_edge(edges_map::Dict{EdgeKey,Vector{Int}}, faces::Vector{Triangle}, ek::EdgeKey)
        idxs = get(edges_map, ek, Int[])
        return [faces[i] for i in idxs if 1 <= i <= length(faces)]
    end

    # Main accumulation
    for sym in candidates
        ek = sym.edge_key
        if sym.repair_strategy == :use_A
            tris = triangles_for_edge(edges_A, faces_A, ek)
            for t in tris
                try_add_triangle!(t, :from_A)
            end
        elseif sym.repair_strategy == :use_B
            tris = triangles_for_edge(edges_B, faces_B, ek)
            for t in tris
                try_add_triangle!(t, :from_B)
            end
        end
    end

    # If nothing was added, attempt a minimal shared-triangle fallback from A
    if isempty(unified_tris)
        verbose && println("  → No triangles from strategies; attempting shared-triangle fallback from A")
        for t in topology.faces_A[1:min(10, length(topology.faces_A))]
            try_add_triangle!(t, :from_A)
        end
    end
    if isempty(unified_tris)
        verbose && println("  → Fallback failed; returning incompatible empty mesh")
        return UnifiedInterfaceMesh(
            unified_tris,
            Dict{NTuple{3,Float64}, Union{Int, Nothing}}(),
            Dict{NTuple{3,Float64}, Union{Int, Nothing}}(),
            Dict{EdgeKey, Vector{Int}}(),
            Symbol[],
            0.0,
            0.0,
            false,
            false,
            ["Unified mesh is empty after fallback"]
        )
    end

    # Build node mappings: map coords to node ids where shared; else nothing
    node_mapping_A = Dict{NTuple{3,Float64}, Union{Int, Nothing}}()
    node_mapping_B = Dict{NTuple{3,Float64}, Union{Int, Nothing}}()
    for tri in unified_tris
        for c in (tri.coord1, tri.coord2, tri.coord3)
            if !haskey(node_mapping_A, c)
                key = ckey(c)
                if haskey(topology.node_key_to_ids, key)
                    nidA, nidB = topology.node_key_to_ids[key]
                    node_mapping_A[c] = nidA
                    node_mapping_B[c] = nidB
                else
                    node_mapping_A[c] = nothing
                    node_mapping_B[c] = nothing
                end
            end
        end
    end

    # Build edge topology
    function build_edge_topology(tris::Vector{Triangle})
        edges = Dict{EdgeKey, Vector{Int}}()
        for (i, t) in enumerate(tris)
            k1 = ckey(t.coord1); k2 = ckey(t.coord2); k3 = ckey(t.coord3)
            for (a, b) in ((k1, k2), (k2, k3), (k3, k1))
                ek = EdgeKey(a, b)
                push!(get!(edges, ek, Int[]), i)
            end
        end
        return edges
    end
    unified_edges = build_edge_topology(unified_tris)

    # ----------------------------------------------------------------------
    # Gap detection and fill (centroid fan) - simple, assumes near-planar loops
    # ----------------------------------------------------------------------
    # Build a map from rounded coord -> an original coord seen in unified_tris
    function build_ckey_to_coord_map(tris::Vector{Triangle})
        m = Dict{NTuple{3,Float64},NTuple{3,Float64}}()
        for t in tris
            for c in (t.coord1, t.coord2, t.coord3)
                k = ckey(c)
                if !haskey(m, k)
                    m[k] = c
                end
            end
        end
        return m
    end

    function detect_boundary_loops(edges::Dict{EdgeKey,Vector{Int}})
        # Boundary edges are incident to exactly 1 triangle
        boundary_edges = EdgeKey[]
        for (ek, tris) in edges
            if length(tris) == 1
                push!(boundary_edges, ek)
            end
        end
        # Build adjacency (rounded node -> neighbors)
        nbrs = Dict{NTuple{3,Float64}, Vector{NTuple{3,Float64}}}()
        for ek in boundary_edges
            a = ek.node1; b = ek.node2
            push!(get!(nbrs, a, NTuple{3,Float64}[]), b)
            push!(get!(nbrs, b, NTuple{3,Float64}[]), a)
        end
        # Traverse closed loops greedily
        visited = Set{Tuple{NTuple{3,Float64},NTuple{3,Float64}}}()
        loops = Vector{Vector{NTuple{3,Float64}}}()
        for start in keys(nbrs)
            if length(nbrs[start]) < 2
                continue
            end
            for nxt in nbrs[start]
                e = start <= nxt ? (start, nxt) : (nxt, start)
                if e in visited; continue; end
                path = NTuple{3,Float64}[start]
                prev = start
                curr = nxt
                push!(visited, e)
                while true
                    push!(path, curr)
                    neighbors = get(nbrs, curr, NTuple{3,Float64}[])
                    if isempty(neighbors); break; end
                    next_candidates = filter(x -> x != prev, neighbors)
                    if isempty(next_candidates); break; end
                    nxt2 = next_candidates[1]
                    ed = curr <= nxt2 ? (curr, nxt2) : (nxt2, curr)
                    if ed in visited
                        break
                    end
                    push!(visited, ed)
                    if nxt2 == start
                        push!(loops, path)
                        break
                    end
                    prev, curr = curr, nxt2
                end
            end
        end
        return loops
    end

    function triangulate_loop_centroid_fan(loop::Vector{NTuple{3,Float64}}, k2c::Dict{NTuple{3,Float64},NTuple{3,Float64}})
        if length(loop) < 3
            return Triangle[]
        end
        coords = [get(k2c, k, k) for k in loop]
        cx = mean(c[1] for c in coords)
        cy = mean(c[2] for c in coords)
        cz = mean(c[3] for c in coords)
        center = (cx, cy, cz)
        tris = Triangle[]
        for i in 2:(length(coords)-1)
            c1 = center; c2 = coords[i]; c3 = coords[i+1]
            push!(tris, Triangle(1,2,3, i, c1,c2,c3))
        end
        return tris
    end

    # Try to fill gaps using boundary loops
    if !isempty(unified_tris)
        k2c_map = build_ckey_to_coord_map(unified_tris)
        loops = detect_boundary_loops(unified_edges)
        added = 0
        for loop in loops
            new_tris = triangulate_loop_centroid_fan(loop, k2c_map)
            for nt in new_tris
                if try_add_triangle!(nt, :synthesized)
                    added += 1
                end
            end
        end
        if added > 0
            unified_edges = build_edge_topology(unified_tris)
            verbose && println("  → Gap fill added $added triangle(s)")
        end
    end

    # Metrics
    min_quality = 1.0
    total_area = 0.0
    for t in unified_tris
        # compute_triangle_quality is defined in edge_classification.jl
        q = compute_triangle_quality(t)
        min_quality = min(min_quality, q)
        total_area += t.area
    end

    # Compatibility: edge coverage vs each side
    function compare_edge_coverage(unified_edges::Dict{EdgeKey,Vector{Int}}, side_edges::Dict{EdgeKey,Vector{Int}})
        total = length(unified_edges)
        matched = 0
        for ek in keys(unified_edges)
            if haskey(side_edges, ek)
                matched += 1
            end
        end
        return total == 0 ? 0.0 : matched / total
    end
    coverage_A = compare_edge_coverage(unified_edges, topology.edges_A)
    coverage_B = compare_edge_coverage(unified_edges, topology.edges_B)
    coverage_threshold = 0.5
    compatible_A = coverage_A >= coverage_threshold
    compatible_B = coverage_B >= coverage_threshold
    report = String[]
    if !compatible_A
        push!(report, "Low edge coverage vs A: $(round(coverage_A*100, digits=1))% < $(Int(coverage_threshold*100))%")
    end
    if !compatible_B
        push!(report, "Low edge coverage vs B: $(round(coverage_B*100, digits=1))% < $(Int(coverage_threshold*100))%")
    end

    return UnifiedInterfaceMesh(
        unified_tris,
        node_mapping_A,
        node_mapping_B,
        unified_edges,
        provenance,
        min_quality,
        total_area,
        compatible_A,
        compatible_B,
        report
    )
end
