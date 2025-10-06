
"""
    check_interface_conformity(nas_file; tol=1e-4, max_examples_per_pair=20, verbose=true)

Check conformity of interfaces between PIDs on a tetrahedral NAS mesh.

Concept:
- Build boundary faces per PID (faces with odd incidence for that PID).
- Build coordinate-keyed nodes to allow cross-PID matching regardless of node IDs.
- For each PID pair (A,B) sharing node keys, build edge sets from boundary faces whose endpoints exist in both A and B.
- Compare edge sets: edges present on A but not on B (and vice versa) indicate non-conforming triangulations (splits/T-junctions).
- Additionally, detect hanging nodes: nodes from A that lie on B's edges (without being endpoints) and vice versa.

Return shape (simplified):
- If all pairs conform: Dict("status"=>"ok", "pairs_checked"=>N)
- If issues: Dict("status"=>"nonconforming", "problems"=>[ {pidA, pidB, missing_in_B_count, missing_in_B_examples, missing_in_A_count, missing_in_A_examples, hanging_A_on_B_count, hanging_A_on_B_examples, hanging_B_on_A_count, hanging_B_on_A_examples}... ])
"""
function check_interface_conformity(nas_file::AbstractString; tol::Real=1e-4, max_examples_per_pair::Int=20, verbose::Bool=true)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.open(nas_file)

        # Nodes and coordinates
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        coords = Dict{Int,NTuple{3,Float64}}()
        for (i, tag) in enumerate(node_tags)
            idx = (i-1)*3
            coords[Int(tag)] = (node_coords[idx+1], node_coords[idx+2], node_coords[idx+3])
        end

        # Helper: coordinate key with rounding
        function ckey(p::NTuple{3,Float64})
            return (round(p[1]; digits=4), round(p[2]; digits=4), round(p[3]; digits=4))
        end

        # Collect tets per PID (using 3D entity tag as PID proxy)
        region_tets = Dict{Int,Vector{Tuple{Int,NTuple{4,Int}}}}()
        for (dim, tag) in gmsh.model.getEntities(3)
            etypes, etags, enodes = gmsh.model.mesh.getElements(dim, tag)
            for (etype, etag_vec, enode_vec) in zip(etypes, etags, enodes)
                if etype != 4; continue; end
                for i in 1:length(etag_vec)
                    eid = Int(etag_vec[i])
                    base = (i-1)*4
                    nd = (Int(enode_vec[base+1]), Int(enode_vec[base+2]), Int(enode_vec[base+3]), Int(enode_vec[base+4]))
                    push!(get!(region_tets, tag, Vector{Tuple{Int,NTuple{4,Int}}}()), (eid, nd))
                end
            end
        end
        if isempty(region_tets)
            return Dict("status"=>"no_tets")
        end

        # Build face incidence: face (sorted 3-tuple of node IDs) -> Dict(pid=>count)
        face_inc = Dict{NTuple{3,Int},Dict{Int,Int}}()
        tet_faces = ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        for (pid, tets) in region_tets
            for (eid, nd) in tets
                n = (nd[1], nd[2], nd[3], nd[4])
                for f in tet_faces
                    tri = (n[f[1]], n[f[2]], n[f[3]])
                    tri_sorted_v = sort!(collect(tri))
                    d = get!(face_inc, (tri_sorted_v[1],tri_sorted_v[2],tri_sorted_v[3]), Dict{Int,Int}())
                    d[pid] = get(d, pid, 0) + 1
                end
            end
        end

        # Boundary faces per PID = faces with odd count under that PID
        boundary_faces_by_pid = Dict{Int,Vector{NTuple{3,Int}}}()
        for (face, pid_counts) in face_inc
            for (pid, cnt) in pid_counts
                if isodd(cnt)
                    push!(get!(boundary_faces_by_pid, pid, NTuple{3,Int}[]), face)
                end
            end
        end

        # Build node key maps per PID and global reverse map
        nodes_by_pid = Dict{Int,Vector{Int}}()
        pid_keyset = Dict{Int,Set{NTuple{3,Float64}}}()
        key_to_nodes_by_pid = Dict{Int,Dict{NTuple{3,Float64},Vector{Int}}}()
        for pid in keys(region_tets)
            nodes = Int[]
            for (_, nd) in region_tets[pid]
                append!(nodes, (nd[1], nd[2], nd[3], nd[4]))
            end
            nodes = unique(nodes)
            nodes_by_pid[pid] = nodes
            ks = Set{NTuple{3,Float64}}()
            kmap = Dict{NTuple{3,Float64},Vector{Int}}()
            for nid in nodes
                k = ckey(coords[nid])
                push!(ks, k)
                push!(get!(kmap, k, Int[]), nid)
            end
            pid_keyset[pid] = ks
            key_to_nodes_by_pid[pid] = kmap
        end

        # Helper: build edge set from faces for a PID restricted to shared node keys with another PID
        function edges_from_faces(pid::Int, other_pid::Int)
            faces = get(boundary_faces_by_pid, pid, NTuple{3,Int}[])
            shared_keys = intersect(pid_keyset[pid], pid_keyset[other_pid])
            if isempty(shared_keys)
                return Dict{Tuple{NTuple{3,Float64},NTuple{3,Float64}},Int}()
            end
            shared = Set(shared_keys)
            edge_set = Dict{Tuple{NTuple{3,Float64},NTuple{3,Float64}},Int}()
            for tri in faces
                k1 = ckey(coords[tri[1]])
                k2 = ckey(coords[tri[2]])
                k3 = ckey(coords[tri[3]])
                # Only keep edges whose endpoints are present in both PIDs
                for (a,b) in ((k1,k2),(k1,k3),(k2,k3))
                    if (a in shared) && (b in shared)
                        ek = a <= b ? (a,b) : (b,a)
                        edge_set[ek] = get(edge_set, ek, 0) + 1
                    end
                end
            end
            return edge_set
        end

        # Geometry helpers
        function dot(a,b); return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]; end
        function sub(a,b); return (a[1]-b[1], a[2]-b[2], a[3]-b[3]); end
        function norm2(a); return dot(a,a); end
        function point_on_segment_dist2(p,a,b)
            ab = sub(b,a)
            ap = sub(p,a)
            denom = norm2(ab)
            if denom <= (tol^2)
                return norm2(sub(p,a)), 0.0
            end
            t = clamp(dot(ap, ab)/denom, 0.0, 1.0)
            proj = (a[1]+t*ab[1], a[2]+t*ab[2], a[3]+t*ab[3])
            return norm2(sub(p, proj)), t
        end

        # Scan PID pairs
        pids = sort(collect(keys(region_tets)))
        problems = Any[]
        pairs_checked = 0
        for i in 1:length(pids)-1
            for j in i+1:length(pids)
                pidA = pids[i]; pidB = pids[j]
                # Quick adjacency test: share any node key
                if isempty(intersect(pid_keyset[pidA], pid_keyset[pidB]))
                    continue
                end
                pairs_checked += 1

                edgesA = edges_from_faces(pidA, pidB)
                edgesB = edges_from_faces(pidB, pidA)

                # Build set difference by keys
                setA = Set(keys(edgesA)); setB = Set(keys(edgesB))
                onlyA = setdiff(setA, setB)
                onlyB = setdiff(setB, setA)

                # Examples formatting helper
                function edge_example_list(es)
                    ex = Any[]
                    for ek in es
                        if length(ex) >= max_examples_per_pair; break; end
                        a = ek[1]; b = ek[2]
                        push!(ex, [ [a[1],a[2],a[3]], [b[1],b[2],b[3]] ])
                    end
                    return ex
                end

                # Hanging nodes: nodes from one PID lying on edges of the other
                function hanging_nodes_on_edges(pid_from::Int, edges_other::Set{Tuple{NTuple{3,Float64},NTuple{3,Float64}}})
                    # Candidate nodes restricted to shared key set (interface vicinity)
                    shared = intersect(pid_keyset[pid_from], pid_keyset[pid_from == pidA ? pidB : pidA])
                    # Build coordinate list of unique nodes
                    nodes_coords = NTuple{3,Float64}[]
                    for k in shared
                        push!(nodes_coords, k)
                    end
                    # Prepare edge geometries
                    egeom = [ (ek[1], ek[2]) for ek in edges_other ]
                    count = 0
                    examples = Any[]
                    for p in nodes_coords
                        # skip if p equals an endpoint to avoid double counting
                        for (a,b) in egeom
                            if p == a || p == b; continue; end
                            d2, t = point_on_segment_dist2(p, a, b)
                            if d2 <= (tol^2) && t > 1e-6 && t < 1-1e-6
                                count += 1
                                if length(examples) < max_examples_per_pair
                                    push!(examples, Dict("point"=>[p[1],p[2],p[3]], "edge"=>[[a[1],a[2],a[3]],[b[1],b[2],b[3]]]))
                                end
                                break
                            end
                        end
                        if length(examples) >= max_examples_per_pair
                            break
                        end
                    end
                    return count, examples
                end

                hA_count, hA_examples = hanging_nodes_on_edges(pidA, Set(keys(edgesB)))
                hB_count, hB_examples = hanging_nodes_on_edges(pidB, Set(keys(edgesA)))

                if !isempty(onlyA) || !isempty(onlyB) || hA_count > 0 || hB_count > 0
                    push!(problems, Dict(
                        "pidA"=>pidA, "pidB"=>pidB,
                        "missing_in_B_count"=>length(onlyA),
                        "missing_in_B_examples"=>edge_example_list(onlyA),
                        "missing_in_A_count"=>length(onlyB),
                        "missing_in_A_examples"=>edge_example_list(onlyB),
                        "hanging_A_on_B_count"=>hA_count,
                        "hanging_A_on_B_examples"=>hA_examples,
                        "hanging_B_on_A_count"=>hB_count,
                        "hanging_B_on_A_examples"=>hB_examples
                    ))
                    if verbose
                        println("  Pair (PID=$(pidA), PID=$(pidB)): non-conforming -> missingInB=$(length(onlyA)), missingInA=$(length(onlyB)), hangA=$(hA_count), hangB=$(hB_count)")
                    end
                else
                    if verbose
                        println("  Pair (PID=$(pidA), PID=$(pidB)): ✓ conforming")
                    end
                end
            end
        end

        if isempty(problems)
            return Dict("status"=>"ok", "pairs_checked"=>pairs_checked)
        else
            return Dict("status"=>"nonconforming", "problems"=>problems)
        end
    finally
        gmsh.finalize()
    end
end

"""
    export_interface_conformity_json(nas_file; output_json="interface_conformity_report.json", tol=1e-4)

Run `check_interface_conformity` and write a compact JSON report.
"""
function export_interface_conformity_json(nas_file::AbstractString; output_json::AbstractString="interface_conformity_report.json", tol::Real=1e-4)
    report = check_interface_conformity(nas_file; tol=tol, verbose=true)
    open(output_json, "w") do io
        write_json(io, report, 0)
    end
    println("Interface conformity report exported to: $(output_json)")
    return output_json
end

"""
    check_interface_conformity_json(nas_file; output_json="interface_conformity_report.json", tol=1e-4)

Back-compat alias for export_interface_conformity_json.
"""
function check_interface_conformity_json(nas_file::AbstractString; output_json::AbstractString="interface_conformity_report.json", tol::Real=1e-4)
    return export_interface_conformity_json(nas_file; output_json=output_json, tol=tol)
end
