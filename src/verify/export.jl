
"""
    export_mesh_quality_json(verification_result, output_file="mesh_quality.json"; 
                            include_anomalies_from=nothing)

Export comprehensive mesh quality verification results to JSON.
Includes volume checks, quality metrics, coordination numbers, and closure status.

Optional:
- include_anomalies_from: Path to *_anomalies.json file from nas_to_step conversion to merge
"""
function export_mesh_quality_json(verification_result, output_file="mesh_quality.json";
                                  include_anomalies_from::Union{Nothing,AbstractString}=nothing)
    # Build JSON structure
    report = Dict(
        "overall_status" => String(verification_result.overall_status),
        "volume_check" => Dict(
            "total_elements" => verification_result.volume_check.total_elements,
            "inverted_count" => verification_result.volume_check.inverted_count,
            "zero_volume_count" => verification_result.volume_check.zero_volume_count,
            "min_volume" => verification_result.volume_check.min_volume,
            "max_volume" => verification_result.volume_check.max_volume,
            "convention_mismatch" => verification_result.volume_check.convention_mismatch,
            "status" => String(verification_result.volume_check.status)
        ),
        "quality_check" => Dict(
            "min" => verification_result.quality_check.min,
            "max" => verification_result.quality_check.max,
            "mean" => verification_result.quality_check.mean,
            "median" => verification_result.quality_check.median,
            "poor_count" => verification_result.quality_check.poor_count,
            "bad_count" => verification_result.quality_check.bad_count,
            "status" => String(verification_result.quality_check.status)
        ),
        "coordination_check" => Dict(
            "total_vertices" => verification_result.coordination_check.total_vertices,
            "mean_coordination" => verification_result.coordination_check.mean_coordination,
            "median_coordination" => verification_result.coordination_check.median_coordination,
            "p90_coordination" => verification_result.coordination_check.p90_coordination,
            "p95_coordination" => verification_result.coordination_check.p95_coordination,
            "p99_coordination" => verification_result.coordination_check.p99_coordination,
            "min_coord_threshold" => verification_result.coordination_check.min_coord,
            "max_coord_threshold" => verification_result.coordination_check.max_coord_threshold,
            "overcoord_percentile" => verification_result.coordination_check.overcoord_percentile,
            "undercoordinated_count" => length(verification_result.coordination_check.undercoordinated),
            "overcoordinated_count" => length(verification_result.coordination_check.overcoordinated),
            "undercoordinated_vertices" => [Dict("vertex_id" => v[1], "coordination" => v[2]) 
                                            for v in verification_result.coordination_check.undercoordinated],
            "overcoordinated_vertices" => [Dict("vertex_id" => v[1], "coordination" => v[2]) 
                                           for v in verification_result.coordination_check.overcoordinated],
            "status" => String(verification_result.coordination_check.status)
        ),
        "closure_check" => Dict(
            "is_closed" => verification_result.closure_check.is_closed,
            "boundary_edge_count" => verification_result.closure_check.boundary_edge_count,
            "surface_count" => verification_result.closure_check.surface_count,
            "status" => String(verification_result.closure_check.status)
        )
    )
    
    # Add swap test if present
    if verification_result.swap_test !== nothing
        report["swap_test"] = Dict(
            "original_inverted" => verification_result.swap_test.original_inverted,
            "swapped_inverted" => verification_result.swap_test.swapped_inverted,
            "would_fix" => verification_result.swap_test.would_fix
        )
    end
    
    # Merge anomalies from conversion if provided
    if include_anomalies_from !== nothing && isfile(include_anomalies_from)
        try
            anomaly_content = read(include_anomalies_from, String)
            # Simple parse to extract anomalies array
            if occursin("\"anomalies\":", anomaly_content)
                report["conversion_anomalies_file"] = include_anomalies_from
                report["conversion_anomalies_note"] = "Manifoldness issues from NAS to STEP conversion (see separate file for details)"
            end
        catch e
            @warn "Could not read anomaly file: $include_anomalies_from"
        end
    end
    
    # Write to file (simple JSON without external deps)
    open(output_file, "w") do io
        write_json(io, report, 0)
    end
    
    println("Mesh quality report exported to: $output_file")
    return output_file
end



"""
    export_inverted_elements(inverted_elements, output_file="inverted_elements.txt")

Export detailed information about inverted elements to a text file.
Useful for debugging mesh inversions.
"""
function export_inverted_elements(inverted_elements, output_file="inverted_elements.txt")
    if length(inverted_elements) == 0
        println("No inverted elements to export.")
        return
    end
    
    open(output_file, "w") do io
        println(io, "="^80)
        println(io, "INVERTED ELEMENT REPORT")
        println(io, "Total inverted elements: $(length(inverted_elements))")
        println(io, "="^80)
        println(io)
        
        # Sort by volume (most negative first)
        sorted_elems = sort(inverted_elements, by = e -> e.volume)
        
        for (i, elem) in enumerate(sorted_elems)
            println(io, "[$i] Element ID: $(elem.element_id)")
            println(io, "    Volume: $(elem.volume)")
            println(io, "    Node IDs: $(elem.node_ids)")
            println(io, "    Centroid: $(elem.centroid)")
            println(io, "    Edge length: min=$(elem.min_edge), max=$(elem.max_edge)")
            println(io, "    Node coordinates:")
            for (j, coord) in enumerate(elem.coordinates)
                println(io, "      Node $(elem.node_ids[j]): $coord")
            end
            println(io)
        end
        
        println(io, "="^80)
        println(io, "END OF REPORT")
        println(io, "="^80)
    end
    
    println("Inverted element details exported to: $output_file")
end


"""
    export_interface_mismatch_surfaces(nas_file; out_dir="interface_mismatch_surfaces", tol=1e-4, verbose=true)

Create a NAS file per PID pair that shows only the mismatched boundary triangles on each side of the interface.
Each file contains two PSHELL regions: one for pidA mismatches and one for pidB mismatches.
Returns a Dict with status and the list of generated file paths.
"""
function export_interface_mismatch_surfaces(nas_file::AbstractString; out_dir::AbstractString="interface_mismatch_surfaces", tol::Real=1e-4, verbose::Bool=true, include_full::Bool=true)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.open(nas_file)

        isdir(out_dir) || mkpath(out_dir)

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
            return Dict("status"=>"no_tets", "files"=>String[])
        end

        # Build face incidence: face (sorted 3-tuple of node IDs) -> Dict(pid=>count)
        face_inc = Dict{NTuple{3,Int},Dict{Int,Int}}()
        tet_faces = ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        for (pid, tets) in region_tets
            for (_, nd) in tets
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

        # Node key sets per PID
        pid_keyset = Dict{Int,Set{NTuple{3,Float64}}}()
        for pid in keys(region_tets)
            ks = Set{NTuple{3,Float64}}()
            for (_, nd) in region_tets[pid]
                for k in (ckey(coords[nd[1]]), ckey(coords[nd[2]]), ckey(coords[nd[3]]), ckey(coords[nd[4]]))
                    push!(ks, k)
                end
            end
            pid_keyset[pid] = ks
        end

        # Helper: edges from faces for PID restricted to shared keys
        function edges_from_faces(pid::Int, other_pid::Int)
            faces = get(boundary_faces_by_pid, pid, NTuple{3,Int}[])
            shared_keys = intersect(pid_keyset[pid], pid_keyset[other_pid])
            shared = Set(shared_keys)
            edge_set = Dict{Tuple{NTuple{3,Float64},NTuple{3,Float64}},Int}()
            for tri in faces
                k1 = ckey(coords[tri[1]]); k2 = ckey(coords[tri[2]]); k3 = ckey(coords[tri[3]])
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
            ab = sub(b,a); ap = sub(p,a)
            denom = norm2(ab)
            if denom <= (tol^2)
                return norm2(sub(p,a)), 0.0
            end
            t = clamp(dot(ap, ab)/denom, 0.0, 1.0)
            proj = (a[1]+t*ab[1], a[2]+t*ab[2], a[3]+t*ab[3])
            return norm2(sub(p, proj)), t
        end

        # Hanging node detection returning key set
        function hanging_node_keys_on_edges(pid_from::Int, edges_other::Set{Tuple{NTuple{3,Float64},NTuple{3,Float64}}}, other_pid::Int)
            shared = intersect(pid_keyset[pid_from], pid_keyset[other_pid])
            nodes_coords = NTuple{3,Float64}[]
            for k in shared
                push!(nodes_coords, k)
            end
            egeom = [ (ek[1], ek[2]) for ek in edges_other ]
            hset = Set{NTuple{3,Float64}}()
            for p in nodes_coords
                for (a,b) in egeom
                    if p == a || p == b; continue; end
                    d2, t = point_on_segment_dist2(p, a, b)
                    if d2 <= (tol^2) && t > 1e-6 && t < 1-1e-6
                        push!(hset, p)
                        break
                    end
                end
            end
            return hset
        end

        # Face selection helpers
        function face_has_edge_from_set(face::NTuple{3,Int}, edgekeys::Set{Tuple{NTuple{3,Float64},NTuple{3,Float64}}})
            k1 = ckey(coords[face[1]]); k2 = ckey(coords[face[2]]); k3 = ckey(coords[face[3]])
            for (a,b) in ((k1,k2),(k1,k3),(k2,k3))
                ek = a <= b ? (a,b) : (b,a)
                if ek in edgekeys
                    return true
                end
            end
            return false
        end
        function face_has_node_from_set(face::NTuple{3,Int}, nodeset::Set{NTuple{3,Float64}})
            k1 = ckey(coords[face[1]]); k2 = ckey(coords[face[2]]); k3 = ckey(coords[face[3]])
            return (k1 in nodeset) || (k2 in nodeset) || (k3 in nodeset)
        end

        # Build files
        files = String[]
        pids = sort(collect(keys(region_tets)))
        for ia in 1:length(pids)-1
            for ib in ia+1:length(pids)
                pidA = pids[ia]; pidB = pids[ib]
                shared_keys = intersect(pid_keyset[pidA], pid_keyset[pidB])
                isempty(shared_keys) && continue

                edgesA = edges_from_faces(pidA, pidB)
                edgesB = edges_from_faces(pidB, pidA)
                setA = Set(keys(edgesA)); setB = Set(keys(edgesB))
                onlyA = setdiff(setA, setB)
                onlyB = setdiff(setB, setA)

                # Hanging node sets
                hA = hanging_node_keys_on_edges(pidA, setB, pidB)
                hB = hanging_node_keys_on_edges(pidB, setA, pidA)

                if isempty(onlyA) && isempty(onlyB) && isempty(hA) && isempty(hB)
                    if verbose
                        println("  Pair (PID=$(pidA), PID=$(pidB)): conforming, no surface export")
                    end
                    continue
                end

                # Select faces to include for each side
                facesA = get(boundary_faces_by_pid, pidA, NTuple{3,Int}[])
                facesB = get(boundary_faces_by_pid, pidB, NTuple{3,Int}[])
                selA = NTuple{3,Int}[]           # mismatch faces A
                selB = NTuple{3,Int}[]           # mismatch faces B
                selA_full = NTuple{3,Int}[]      # full interface faces A
                selB_full = NTuple{3,Int}[]      # full interface faces B
                onlyAset = Set(onlyA)
                onlyBset = Set(onlyB)
                shared_set = Set(shared_keys)

                # helper to check if all nodes in shared node set
                function face_all_nodes_in_set(face::NTuple{3,Int}, nodeset::Set{NTuple{3,Float64}})
                    k1 = ckey(coords[face[1]]); k2 = ckey(coords[face[2]]); k3 = ckey(coords[face[3]])
                    return (k1 in nodeset) && (k2 in nodeset) && (k3 in nodeset)
                end
                for f in facesA
                    # full interface faces for A
                    if include_full && face_all_nodes_in_set(f, shared_set)
                        push!(selA_full, f)
                    end
                    if face_has_edge_from_set(f, onlyAset) || (!isempty(hA) && face_has_node_from_set(f, hA))
                        push!(selA, f)
                    end
                end
                for f in facesB
                    if include_full && face_all_nodes_in_set(f, shared_set)
                        push!(selB_full, f)
                    end
                    if face_has_edge_from_set(f, onlyBset) || (!isempty(hB) && face_has_node_from_set(f, hB))
                        push!(selB, f)
                    end
                end

                if isempty(selA) && isempty(selB)
                    if verbose
                        println("  Pair (PID=$(pidA), PID=$(pidB)): no specific mismatch faces, skipping export")
                    end
                    continue
                end

                # Build SurfaceRegion inputs
                function make_region(name::String, faces::Vector{NTuple{3,Int}})
                    if isempty(faces)
                        return SurfaceRegion(name, NTuple{3,Float64}[], NTuple{3,Int}[], 1.0)
                    end
                    # Collect unique keys
                    klist = NTuple{3,Float64}[]
                    key2idx = Dict{NTuple{3,Float64},Int}()
                    tri_local = NTuple{3,Int}[]
                    for f in faces
                        ks = (ckey(coords[f[1]]), ckey(coords[f[2]]), ckey(coords[f[3]]))
                        idxs = Int[]
                        for kk in ks
                            id = get(key2idx, kk, 0)
                            if id == 0
                                push!(klist, kk)
                                id = length(klist)
                                key2idx[kk] = id
                            end
                            push!(idxs, id)
                        end
                        push!(tri_local, (idxs[1], idxs[2], idxs[3]))
                    end
                    return SurfaceRegion(name, klist, tri_local, 1.0)
                end

                regions = SurfaceRegion[]
                if include_full
                    push!(regions, make_region("PID$(pidA)_full", selA_full))
                    push!(regions, make_region("PID$(pidB)_full", selB_full))
                end
                push!(regions, make_region("PID$(pidA)_mismatch", selA))
                push!(regions, make_region("PID$(pidB)_mismatch", selB))

                out_path = joinpath(out_dir, "iface_pid$(pidA)_pid$(pidB)_mismatch.nas")
                write_nas_surface(out_path, regions; tol=1e-4)
                push!(files, out_path)
                if verbose
                    if include_full
                        println("  Exported interface+highlight: $(out_path)  (A full=$(length(selA_full)), B full=$(length(selB_full)), A mismatch=$(length(selA)), B mismatch=$(length(selB)))")
                    else
                        println("  Exported mismatch surfaces: $(out_path)  (A faces=$(length(selA)), B faces=$(length(selB)))")
                    end
                end
            end
        end

        status = isempty(files) ? "no_mismatches" : "ok"
        return Dict("status"=>status, "files"=>files)
    finally
        gmsh.finalize()
    end
end
