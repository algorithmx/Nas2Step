# Comprehensive interface analysis and JSON export

"""
    export_interface_mismatches_json(nas_file::AbstractString, output_file::AbstractString; min_shared::Int=10)

Analyze all volume pairs (PIDs) that share an interface in the given NAS file,
classify edge mismatches for each interface, compute relevant constraints and
context metrics, and export a complete structured JSON report.

The JSON contains per-interface topology summaries, densities, locked ratios,
suggested repair direction, and full mismatch entries (including constraint
violations and diagonal retriangulation context when available).
"""
function export_interface_mismatches_json(nas_file::AbstractString, output_file::AbstractString; min_shared::Int=10)
    isfile(nas_file) || error("NAS file not found: $(nas_file)")

    # Discover interface pairs (3D gmsh volume tag pairs sharing >= min_shared nodes)
    pairs = _find_interface_pairs(nas_file; min_shared=min_shared)

    interfaces = Vector{Dict{String,Any}}()

    # Global counters
    total_mismatch_count = 0
    total_by_type = Dict{String,Int}("T_JUNCTION"=>0, "DIAGONAL"=>0, "REFINEMENT"=>0, "UNKNOWN"=>0)

    for (pidA, pidB) in pairs
        topo = build_interface_topology(String(nas_file), pidA, pidB)
        cls  = classify_interface_mismatches(topo)
        cons = build_boundary_constraints(String(nas_file), pidA, pidB)

        # Edge densities and locked ratios (replicate determine_dominant_side without println)
        density_A = compute_interface_area(topo, :A) > 1e-10 ? length(topo.edges_A) / compute_interface_area(topo, :A) : 0.0
        density_B = compute_interface_area(topo, :B) > 1e-10 ? length(topo.edges_B) / compute_interface_area(topo, :B) : 0.0

        interface_edges_A = length(topo.edges_A)
        interface_edges_B = length(topo.edges_B)
        locked_ratio_A = min(1.0, cons.total_external_edges_A / max(1, interface_edges_A))
        locked_ratio_B = min(1.0, cons.total_external_edges_B / max(1, interface_edges_B))

        # Suggest repair direction using the same heuristics
        suggested = let
            if density_A > density_B * 1.2
                :subdivide_B
            elseif density_B > density_A * 1.2
                :subdivide_A
            elseif locked_ratio_A < locked_ratio_B * 0.8
                :subdivide_A  # modify A (target pattern from B)
            elseif locked_ratio_B < locked_ratio_A * 0.8
                :subdivide_B  # modify B (target pattern from A)
            else
                topo.total_faces_A < topo.total_faces_B ? :subdivide_A : :subdivide_B
            end
        end

        # Build a per-interface points table (node_id => coordinates)
        points_map = Dict{Int,NTuple{3,Float64}}()
        function _collect_points!(t::Triangle)
            points_map[t.node1] = t.coord1
            points_map[t.node2] = t.coord2
            points_map[t.node3] = t.coord3
        end
        for tri in topo.faces_A; _collect_points!(tri); end
        for tri in topo.faces_B; _collect_points!(tri); end

        # Rounding key for coordinate lookups (shared across helpers)
        ckey(p::NTuple{3,Float64}) = (round(p[1]; digits=4), round(p[2]; digits=4), round(p[3]; digits=4))

        # Helper to serialize a triangle via node IDs
        tri_dict = (t, idx) -> Dict(
            "index" => idx,
            "element_id" => t.element_id,
            "nodes" => [t.node1, t.node2, t.node3],
        )

        # Helper to map a set of vertex coordinates to node IDs on a given side ("A"|"B")
        function vertices_node_ids(vertices::Vector{NTuple{3,Float64}}, side::String)
            ids = Int[]
            for v in vertices
                pair = get(topo.node_key_to_ids, ckey(v), nothing)
                if pair === nothing
                    # If mapping is missing, skip this vertex (should be rare on interface)
                    continue
                end
                push!(ids, side == "A" ? pair[1] : pair[2])
            end
            return ids
        end

        # Helper to get edge node IDs on a given side ("A" or "B")
        function edge_node_ids(ek::EdgeKey, side::String)
            # EdgeKey nodes are coordinate tuples; look up shared mapping
            # node_key_to_ids maps coord_key => (nodeA_id, nodeB_id)
            nidpair1 = get(topo.node_key_to_ids, ek.node1, nothing)
            nidpair2 = get(topo.node_key_to_ids, ek.node2, nothing)
            if nidpair1 === nothing || nidpair2 === nothing
                return Int[]
            end
            if side == "A"
                return [nidpair1[1], nidpair2[1]]
            else
                return [nidpair1[2], nidpair2[2]]
            end
        end

        # Prepare per-mismatch serialization with constraint checks
        # Small helpers to reduce repetition
        faces_for_side = side -> (side == "A" ? topo.faces_A : topo.faces_B)
        target_side_for = m -> (m.should_be_in == topo.pidA ? "A" : "B")

        function serialize_triangles_from_indices(indices::Vector{Int}, faces::Vector{Triangle})
            out = Any[]
            for idx in indices
                if 1 <= idx <= length(faces)
                    tri = faces[idx]
                    push!(out, tri_dict(tri, idx))
                end
            end
            return out
        end

        # Reasons for emptiness tagging
        function diagonal_context_empty_reason(m, has_vio::Bool)
            # Empty if either quad vertices missing or no triangles selected
            if isempty(m.quad_vertices) || isempty(m.triangles_to_replace)
                if has_vio
                    return "BLOCKED_BY_CONSTRAINTS"
                end
                if isempty(m.triangles_to_replace)
                    return "NO_REPLACEABLE_TRIANGLES"
                end
                # Triangles chosen but no quad vertices extracted
                if length(m.triangles_to_replace) == 2 && isempty(m.quad_vertices)
                    return "QUAD_VERTICES_MISSING"
                end
                return "NO_TARGET_QUAD"
            end
            return nothing
        end

        function target_affected_empty_reason(m, affected::Vector{Any}, has_vio::Bool)
            if !isempty(affected)
                return nothing
            end
            if has_vio
                return "BLOCKED_BY_CONSTRAINTS"
            end
            if !m.repair_feasible
                return "REPAIR_INFEASIBLE"
            end
            if isempty(m.triangles_to_replace)
                return "NO_REPLACEABLE_TRIANGLES"
            end
            return "UNKNOWN"
        end

        function serialize_mismatch(m)
            has_vio, reasons = check_constraint_violations(m.edge_key, cons)
            tgt_side = target_side_for(m)
            tgt_faces = faces_for_side(tgt_side)
            affected = serialize_triangles_from_indices(m.affected_triangles, tgt_faces)
            replace_tris = serialize_triangles_from_indices(m.triangles_to_replace, tgt_faces)
            dc_reason = diagonal_context_empty_reason(m, has_vio)
            ta_reason = target_affected_empty_reason(m, affected, has_vio)
            d = Dict(
                "edge" => edge_node_ids(m.edge_key, tgt_side),
                "type" => string(m.mismatch_type),
                "present_in" => string(m.present_in),
                "should_be_in_pid" => m.should_be_in,
                "hanging_nodes" => [[n...] for n in m.hanging_nodes],
                "affected_triangles" => affected,
                "quad_vertices" => vertices_node_ids(m.quad_vertices, tgt_side),
                "triangles_to_replace" => replace_tris,
                "complexity_score" => round(m.complexity_score, digits=3),
                "min_quality" => round(m.min_affected_triangle_quality, digits=3),
                "repair_feasible" => m.repair_feasible,
                "violates_constraints" => has_vio,
                "constraint_violations" => reasons,
            )
            if dc_reason !== nothing
                d["diagonal_context_empty_reason"] = dc_reason
            end
            if ta_reason !== nothing
                d["target_affected_empty_reason"] = ta_reason
            end
            return d
        end

        mismA = [serialize_mismatch(m) for m in cls.mismatches_A]
        mismB = [serialize_mismatch(m) for m in cls.mismatches_B]

        # Build side-by-side contrast entries to place incompatible pieces together
        function triangles_with_edge(faces::Vector{Triangle}, ek::EdgeKey)
            out = Any[]
            for (idx, tri) in enumerate(faces)
                tri_keys = Set([ckey(tri.coord1), ckey(tri.coord2), ckey(tri.coord3)])
                if ek.node1 in tri_keys && ek.node2 in tri_keys
                    push!(out, tri_dict(tri, idx))
                end
            end
            return out
        end

        contrast_entries = Any[]
        # Contrast serialization helper
        function serialize_contrast_entry(m, target_side::String, source_side::String)
            has_vio, reasons = check_constraint_violations(m.edge_key, cons)
            source_faces = faces_for_side(source_side)
            target_faces = faces_for_side(target_side)
            source_tris = triangles_with_edge(source_faces, m.edge_key)
            affected = serialize_triangles_from_indices(m.affected_triangles, target_faces)
            replace_tris = serialize_triangles_from_indices(m.triangles_to_replace, target_faces)
            dc_reason = diagonal_context_empty_reason(m, has_vio)
            ta_reason = target_affected_empty_reason(m, affected, has_vio)
            d = Dict(
                "edge" => edge_node_ids(m.edge_key, target_side),
                "missing_in" => target_side,
                "source_side" => source_side,
                "target_side" => target_side,
                "classification" => string(m.mismatch_type),
                "complexity" => round(m.complexity_score, digits=3),
                "min_quality" => round(m.min_affected_triangle_quality, digits=3),
                "source_triangles" => source_tris,
                "target_affected" => affected,
                "diagonal_context" => Dict{String,Any}(
                    "quad_vertices" => vertices_node_ids(m.quad_vertices, target_side),
                    "triangles_to_replace" => replace_tris,
                ),
                "constraints_on_target" => Dict(
                    "violates_constraints" => has_vio,
                    "reasons" => reasons,
                )
            )
            if dc_reason !== nothing
                d["diagonal_context"]["empty_reason"] = dc_reason
            end
            if ta_reason !== nothing
                d["target_affected_empty_reason"] = ta_reason
            end
            return d
        end

        # Edges missing in A (present in B), and missing in B (present in A)
        for m in cls.mismatches_A
            push!(contrast_entries, serialize_contrast_entry(m, "A", "B"))
        end
        for m in cls.mismatches_B
            push!(contrast_entries, serialize_contrast_entry(m, "B", "A"))
        end

        # Update global counters
        total_mismatch_count += length(cls.mismatches_A) + length(cls.mismatches_B)
        total_by_type["T_JUNCTION"] += cls.t_junction_count
        total_by_type["DIAGONAL"]   += cls.diagonal_count
        total_by_type["REFINEMENT"] += cls.refinement_count
        total_by_type["UNKNOWN"]    += cls.unknown_count

        push!(interfaces, Dict(
            "pidA" => pidA,
            "pidB" => pidB,
            "conformity_ratio" => round(topo.conformity_ratio, digits=4),
            "topology" => Dict(
                "total_shared_nodes" => topo.total_shared_nodes,
                "faces_A" => topo.total_faces_A,
                "faces_B" => topo.total_faces_B,
                "edges_A" => topo.total_edges_A,
                "edges_B" => topo.total_edges_B,
                "edges_only_in_A" => length(topo.edges_only_in_A),
                "edges_only_in_B" => length(topo.edges_only_in_B),
                "edges_shared" => length(topo.edges_shared),
            ),
            "geometry" => Dict(
                "bbox_min" => [topo.interface_bbox.min_corner...],
                "bbox_max" => [topo.interface_bbox.max_corner...],
                "area_A" => compute_interface_area(topo, :A),
                "area_B" => compute_interface_area(topo, :B),
            ),
            "densities" => Dict("A"=>density_A, "B"=>density_B),
            "locked_ratios" => Dict("A"=>locked_ratio_A, "B"=>locked_ratio_B),
            "suggested_repair_direction" => string(suggested),
            "constraints_summary" => Dict(
                "external_edges_A" => cons.total_external_edges_A,
                "external_edges_B" => cons.total_external_edges_B,
                "corner_nodes" => cons.total_corner_nodes,
                "total_locked_edges" => cons.total_locked_edges,
            ),
            # points table mapping node_id -> coordinate
            "points" => Dict(string(k) => [v...] for (k,v) in points_map),
            "mismatches" => Dict(
                "missing_in_A" => mismA,
                "missing_in_B" => mismB,
            ),
            "mismatch_contrast" => contrast_entries
        ))
    end

    report = Dict(
        "input_file" => String(nas_file),
        "min_shared_nodes_for_pair" => min_shared,
        "interface_count" => length(pairs),
        "summary" => Dict(
            "total_mismatches" => total_mismatch_count,
            "by_type" => total_by_type,
        ),
        "interfaces" => interfaces,
    )

    open(output_file, "w") do io
        write_json(io, report, 0)
    end

    println("Interface mismatch report written to: $output_file")
    return output_file
end


"""
    _find_interface_pairs(nas_file::AbstractString; min_shared::Int=10) -> Vector{Tuple{Int,Int}}

Internal helper: return all gmsh volume tag pairs whose node coordinate key
sets intersect by at least `min_shared`.
"""
function _find_interface_pairs(nas_file::AbstractString; min_shared::Int=10)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.open(nas_file)

        vols = [tag for (dim, tag) in gmsh.model.getEntities(3)]
        if length(vols) < 2
            return Tuple{Int,Int}[]
        end

        # Build a global node tag -> coordinate lookup once
        tags, coords, _ = gmsh.model.mesh.getNodes()
        t2c = Dict{Int,NTuple{3,Float64}}()
        for (i, t) in enumerate(tags)
            j = (i-1)*3
            t2c[Int(t)] = (coords[j+1], coords[j+2], coords[j+3])
        end

        # Rounding key
        ckey(p::NTuple{3,Float64}) = (round(p[1]; digits=4), round(p[2]; digits=4), round(p[3]; digits=4))

        # Nodes by volume
        vol_nodes = Dict{Int,Set{NTuple{3,Float64}}}()
        for v in vols
            etypes, _, enodes = gmsh.model.mesh.getElements(3, v)
            nodes_here = Set{NTuple{3,Float64}}()
            for (et, ev) in zip(etypes, enodes)
                if et != 4; continue; end
                for i in 1:4:length(ev)
                    for k in 0:3
                        nid = Int(ev[i+k])
                        if haskey(t2c, nid)
                            push!(nodes_here, ckey(t2c[nid]))
                        end
                    end
                end
            end
            vol_nodes[v] = nodes_here
        end

        pairs = Tuple{Int,Int}[]
        for i in 1:length(vols)
            for j in (i+1):length(vols)
                a, b = vols[i], vols[j]
                if length(intersect(vol_nodes[a], vol_nodes[b])) >= min_shared
                    push!(pairs, (a, b))
                end
            end
        end
        return pairs
    finally
        gmsh.finalize()
    end
end
