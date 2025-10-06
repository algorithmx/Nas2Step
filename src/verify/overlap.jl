
"""
    check_region_overlap(nas_file; sample_fraction=0.05, max_checks_per_pair=1000, grid_cells=40, verbose=true)

Detect volumetric overlaps between different PIDs in a NAS tetrahedral mesh.

Method:
- Parse all CTETRA (etype=4) with their PIDs and node coordinates using Gmsh.
- Build a uniform 3D grid binning of tetrahedron AABBs per region (PID) to prune pair tests.
- For each region pair (A,B), sample up to `max_checks_per_pair` tets from each region (or `sample_fraction`) and test whether their centroids lie inside any tet of the other region using the binning.
- If centroid of a tet in A lies inside a tet in B (or vice versa), record an overlap example.

Return shape (simplified):
- If no overlaps detected: Dict("status" => "ok", "pairs_checked" => N)
- If overlaps detected: Dict("status" => "overlap", "overlaps" => [ { pidA, pidB, overlap_hits, sampleA, sampleB, examples=[...] } ... ])

Notes:
- This is probabilistic (sampling) for scalability; it reports likely overlaps with examples. Increase `sample_fraction` or `max_checks_per_pair` for thoroughness.
- Inside test uses barycentric coordinates; degenerate tets are skipped.
"""
struct TetInfo
    eid::Int
    nodes::NTuple{4,Int}
    centroid::NTuple{3,Float64}
    aabb_min::NTuple{3,Float64}
    aabb_max::NTuple{3,Float64}
end

function check_region_overlap(nas_file::AbstractString; sample_fraction::Real=0.05, max_checks_per_pair::Int=1000, grid_cells::Int=40, verbose::Bool=true)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.open(nas_file)

        # collect nodes
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        coords = Dict{Int,NTuple{3,Float64}}()
        for (i, tag) in enumerate(node_tags)
            idx = (i-1)*3
            coords[Int(tag)] = (node_coords[idx+1], node_coords[idx+2], node_coords[idx+3])
        end

    # collect all tetra of all PIDs by entity tag (PID proxy)
        # Map: pid => Vector{(eid, (n1,n2,n3,n4))}
        region_tets = Dict{Int,Vector{Tuple{Int,NTuple{4,Int}}}}()

        # Try to use 3D discrete entities to partition by PID-like tags
        for (dim, tag) in gmsh.model.getEntities(3)
            etypes, etags, enodes = gmsh.model.mesh.getElements(dim, tag)
            for (etype, etag_vec, enode_vec) in zip(etypes, etags, enodes)
                if etype != 4; continue; end
                nper = 4
                for i in 1:length(etag_vec)
                    eid = Int(etag_vec[i])
                    base = (i-1)*nper
                    n1 = Int(enode_vec[base+1]); n2 = Int(enode_vec[base+2]); n3 = Int(enode_vec[base+3]); n4 = Int(enode_vec[base+4])
                    # Use the 3D entity tag as PID proxy if physical group missing
                    pid = tag
                    push!(get!(region_tets, pid, Vector{Tuple{Int,NTuple{4,Int}}}()), (eid, (n1,n2,n3,n4)))
                end
            end
        end

        if isempty(region_tets)
            return Dict("status"=>"no_tets")
        end

        # Compute AABBs and centroids per tet
        function tet_info(eid::Int, nd::NTuple{4,Int})
            p = (coords[nd[1]], coords[nd[2]], coords[nd[3]], coords[nd[4]])
            cx = (p[1][1]+p[2][1]+p[3][1]+p[4][1])/4
            cy = (p[1][2]+p[2][2]+p[3][2]+p[4][2])/4
            cz = (p[1][3]+p[2][3]+p[3][3]+p[4][3])/4
            xs = (p[1][1],p[2][1],p[3][1],p[4][1]); ys=(p[1][2],p[2][2],p[3][2],p[4][2]); zs=(p[1][3],p[2][3],p[3][3],p[4][3])
            aabb_min = (minimum(xs), minimum(ys), minimum(zs))
            aabb_max = (maximum(xs), maximum(ys), maximum(zs))
            return TetInfo(eid, nd, (cx,cy,cz), aabb_min, aabb_max)
        end

        reg_infos = Dict{Int,Vector{TetInfo}}()
        for (pid, tets) in region_tets
            infos = TetInfo[]
            for (eid, nd) in tets
                push!(infos, tet_info(eid, nd))
            end
            reg_infos[pid] = infos
        end

        # Global bounds
        xs = Float64[]; ys = Float64[]; zs = Float64[]
        for infos in values(reg_infos)
            for t in infos
                push!(xs, t.aabb_min[1]); push!(xs, t.aabb_max[1])
                push!(ys, t.aabb_min[2]); push!(ys, t.aabb_max[2])
                push!(zs, t.aabb_min[3]); push!(zs, t.aabb_max[3])
            end
        end
        if isempty(xs)
            return Dict("status"=>"no_bounds")
        end
        minB = (minimum(xs), minimum(ys), minimum(zs))
        maxB = (maximum(xs), maximum(ys), maximum(zs))

        # Simple uniform grid indexer
        function grid_index(p::NTuple{3,Float64})
            gx = clamp(floor(Int, (p[1]-minB[1]) / max(1e-12, (maxB[1]-minB[1])) * grid_cells), 0, grid_cells-1)
            gy = clamp(floor(Int, (p[2]-minB[2]) / max(1e-12, (maxB[2]-minB[2])) * grid_cells), 0, grid_cells-1)
            gz = clamp(floor(Int, (p[3]-minB[3]) / max(1e-12, (maxB[3]-minB[3])) * grid_cells), 0, grid_cells-1)
            return (gx, gy, gz)
        end

        # Insert tets per PID into grid by AABB coverage
        function build_grid(infos::Vector{TetInfo})
            grid = Dict{Tuple{Int,Int,Int},Vector{Int}}()
            for (idx, ti) in enumerate(infos)
                # grid range over aabb
                gmin = grid_index(ti.aabb_min)
                gmax = grid_index(ti.aabb_max)
                for ix in gmin[1]:gmax[1], iy in gmin[2]:gmax[2], iz in gmin[3]:gmax[3]
                    push!(get!(grid, (ix,iy,iz), Int[]), idx)
                end
            end
            return grid
        end

        function point_in_tet(p::NTuple{3,Float64}, nd::NTuple{4,Int})
            a = coords[nd[1]]; b = coords[nd[2]]; c = coords[nd[3]]; d = coords[nd[4]]
            # barycentric using volumes
            function vol_s(u::NTuple{3,Float64}, v::NTuple{3,Float64}, w::NTuple{3,Float64}, x::NTuple{3,Float64})
                dv1 = (v[1]-u[1], v[2]-u[2], v[3]-u[3])
                dv2 = (w[1]-u[1], w[2]-u[2], w[3]-u[3])
                dv3 = (x[1]-u[1], x[2]-u[2], x[3]-u[3])
                det = dv1[1]*(dv2[2]*dv3[3] - dv2[3]*dv3[2]) - dv1[2]*(dv2[1]*dv3[3] - dv2[3]*dv3[1]) + dv1[3]*(dv2[1]*dv3[2] - dv2[2]*dv3[1])
                return det / 6.0
            end
            v0s = vol_s(a,b,c,d)
            if abs(v0s) < 1e-18
                return false
            end
            s1 = vol_s(p,b,c,d)
            s2 = vol_s(a,p,c,d)
            s3 = vol_s(a,b,p,d)
            s4 = vol_s(a,b,c,p)
            # same sign and sum equals v0 (within tolerance)
            same_sign = (s1>=0 && s2>=0 && s3>=0 && s4>=0) || (s1<=0 && s2<=0 && s3<=0 && s4<=0)
            if !same_sign
                return false
            end
            return abs((s1+s2+s3+s4) - v0s) <= 1e-9*abs(v0s)
        end

        # build grids per region
        reg_grids = Dict{Int,Any}()
        for (pid, infos) in reg_infos
            reg_grids[pid] = build_grid(infos)
        end

        # helper: candidate tet indices in B whose AABB grid cell matches centroid of candidate
        function candidate_indices(target_pid::Int, point::NTuple{3,Float64})
            g = reg_grids[target_pid]
            idx = grid_index(point)
            return get(g, idx, Int[])
        end

        pids = sort(collect(keys(reg_infos)))
    report = Dict("pairs"=>Any[], "status"=>"ok")

        # Deterministic uniform sampler to avoid Random dependency
        function uniform_sample_indices(n::Int, k::Int)
            k = clamp(k, 0, n)
            if k == 0
                return Int[]
            elseif k == n
                return collect(1:n)
            else
                # take approximately evenly spaced indices
                idxs = Int[]
                for i in 1:k
                    # position in [1,n]
                    pos = 1 + floor(Int, (i-1) * (n-1) / max(1, k-1))
                    push!(idxs, pos)
                end
                return idxs
            end
        end

        for i in 1:length(pids)-1
            for j in i+1:length(pids)
                pidA = pids[i]; pidB = pids[j]
                infosA = reg_infos[pidA]; infosB = reg_infos[pidB]
                # choose sample sizes
                sA = clamp(round(Int, length(infosA)*sample_fraction), 1, min(length(infosA), max_checks_per_pair))
                sB = clamp(round(Int, length(infosB)*sample_fraction), 1, min(length(infosB), max_checks_per_pair))

                idxA = uniform_sample_indices(length(infosA), sA)
                idxB = uniform_sample_indices(length(infosB), sB)

                overlaps = Int[]
                examples = Any[]

                # test A against B
                for ia in idxA
                    ti = infosA[ia]
                    cand = candidate_indices(pidB, ti.centroid)
                    for ib in cand
                        tj = infosB[ib]
                        # quick AABB check
                        if !(ti.aabb_min[1] <= tj.aabb_max[1] && ti.aabb_max[1] >= tj.aabb_min[1] &&
                             ti.aabb_min[2] <= tj.aabb_max[2] && ti.aabb_max[2] >= tj.aabb_min[2] &&
                             ti.aabb_min[3] <= tj.aabb_max[3] && ti.aabb_max[3] >= tj.aabb_min[3])
                            continue
                        end
                        if point_in_tet(ti.centroid, tj.nodes)
                            push!(overlaps, 1)
                            push!(examples, Dict(
                                "pidA"=>pidA, "tetA"=>ti.eid, "pidB"=>pidB, "tetB"=>tj.eid,
                                "centroid"=>[ti.centroid[1], ti.centroid[2], ti.centroid[3]]
                            ))
                            break
                        end
                    end
                end

                # test B against A
                for ib in idxB
                    tj = infosB[ib]
                    cand = candidate_indices(pidA, tj.centroid)
                    for ia2 in cand
                        ti = infosA[ia2]
                        if !(tj.aabb_min[1] <= ti.aabb_max[1] && tj.aabb_max[1] >= ti.aabb_min[1] &&
                             tj.aabb_min[2] <= ti.aabb_max[2] && tj.aabb_max[2] >= ti.aabb_min[2] &&
                             tj.aabb_min[3] <= ti.aabb_max[3] && tj.aabb_max[3] >= ti.aabb_min[3])
                            continue
                        end
                        if point_in_tet(tj.centroid, ti.nodes)
                            push!(overlaps, 1)
                            push!(examples, Dict(
                                "pidA"=>pidA, "tetA"=>ti.eid, "pidB"=>pidB, "tetB"=>tj.eid,
                                "centroid"=>[tj.centroid[1], tj.centroid[2], tj.centroid[3]]
                            ))
                            break
                        end
                    end
                end

                found = sum(overlaps)
                push!(report["pairs"], Dict(
                    "pidA"=>pidA, "pidB"=>pidB,
                    "sampleA"=>sA, "sampleB"=>sB,
                    "overlap_hits"=>found, "examples"=>examples
                ))
                if verbose
                    println("  Pair (PID=$(pidA), PID=$(pidB)): sampled $(sA)+$(sB), overlap hits=$(found)")
                end
            end
        end
        # Build simplified return: minimal on success, detailed only for overlapping pairs
        pairs = report["pairs"]
        overlaps_only = Any[]
        # Limit number of examples per pair to keep report compact
        local MAX_EXAMPLES = 10
        for p in pairs
            if get(p, "overlap_hits", 0) > 0
                ex = get(p, "examples", Any[])
                ex_trim = length(ex) > MAX_EXAMPLES ? ex[1:MAX_EXAMPLES] : ex
                push!(overlaps_only, Dict(
                    "pidA"=>p["pidA"], "pidB"=>p["pidB"],
                    "overlap_hits"=>p["overlap_hits"],
                    "sampleA"=>p["sampleA"], "sampleB"=>p["sampleB"],
                    "examples"=>ex_trim
                ))
            end
        end
        if isempty(overlaps_only)
            return Dict("status"=>"ok", "pairs_checked"=>length(pairs))
        else
            return Dict("status"=>"overlap", "overlaps"=>overlaps_only)
        end
    finally
        gmsh.finalize()
    end
end

"""
    export_region_overlap_json(nas_file; output_json="region_overlap_report.json", sample_fraction=0.05, max_checks_per_pair=1000, grid_cells=40)

Run `check_region_overlap` and write a JSON report with examples per overlapping pair.
"""
function export_region_overlap_json(nas_file::AbstractString; output_json::AbstractString="region_overlap_report.json", sample_fraction::Real=0.05, max_checks_per_pair::Int=1000, grid_cells::Int=40)
    report = check_region_overlap(nas_file; sample_fraction=sample_fraction, max_checks_per_pair=max_checks_per_pair, grid_cells=grid_cells, verbose=true)
    open(output_json, "w") do io
        write_json(io, report, 0)
    end
    println("Region overlap report exported to: $(output_json)")
    return output_json
end
