"""Helper: Create canonical face representation (sorted node IDs)"""
function _sort_face(n1::Int, n2::Int, n3::Int)
    # Avoid allocating an array; sort 3 ints via swaps
    a = n1
    b = n2
    c = n3
    if a > b
        a, b = b, a
    end
    if b > c
        b, c = c, b
    end
    if a > b
        a, b = b, a
    end
    return (a, b, c)
end


"""Helper: Get 4 triangular faces of a tetrahedron"""
function _get_tet_faces(n1::Int, n2::Int, n3::Int, n4::Int)
    return [
        _sort_face(n1, n2, n3),
        _sort_face(n1, n2, n4),
        _sort_face(n1, n3, n4),
        _sort_face(n2, n3, n4)
    ]
end


"""Helper: Compute face normal using cross product"""
function _compute_face_normal(p1, p2, p3)
    v1 = (p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3])
    v2 = (p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3])

    # Normal = v1 × v2
    nx = v1[2] * v2[3] - v1[3] * v2[2]
    ny = v1[3] * v2[1] - v1[1] * v2[3]
    nz = v1[1] * v2[2] - v1[2] * v2[1]

    return (nx, ny, nz)
end


"""Helper: Compute tetrahedron centroid"""
function _compute_tet_centroid(p1, p2, p3, p4)
    cx = (p1[1] + p2[1] + p3[1] + p4[1]) / 4.0
    cy = (p1[2] + p2[2] + p3[2] + p4[2]) / 4.0
    cz = (p1[3] + p2[3] + p3[3] + p4[3]) / 4.0
    return (cx, cy, cz)
end


"""Helper: Compute face centroid"""
function _compute_face_centroid(p1, p2, p3)
    cx = (p1[1] + p2[1] + p3[1]) / 3.0
    cy = (p1[2] + p2[2] + p3[2]) / 3.0
    cz = (p1[3] + p2[3] + p3[3]) / 3.0
    return (cx, cy, cz)
end


"""Helper: Dot product"""
function _dot(v1, v2)
    return v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]
end


"""
    has_surface_elements(nas_path::AbstractString) -> Bool

Check if a NAS file contains surface elements (CTRIA3/CQUAD4).
Returns true if surface elements are found, false otherwise.
"""
function has_surface_elements(nas_path::AbstractString)
    for line in eachline(nas_path)
        if startswith(line, "CTRIA3") || startswith(line, "CQUAD4")
            return true
        end
    end
    return false
end


"""
    extract_boundary_surfaces(input_nas::AbstractString, output_nas::AbstractString;
                             verbose::Bool=true, include_internal_interfaces::Bool=true)

Extract external boundary surfaces from a volume-only NAS file containing CTETRA elements,
and optionally the internal interfaces between regions (different PIDs).

This function:
1. Reads all CTETRA volume elements from the input NAS file
2. Identifies boundary faces (faces that appear only once = external boundaries)
3. Creates CTRIA3 surface elements with correct outward-facing orientation
     - for external faces (one adjacent tet): outward is to outside of the mesh
     - for internal interfaces (two adjacent tets with different PIDs): create two faces,
         one for each PID, oriented outward relative to that PID's tet
4. Writes a new NAS file with both surface and volume elements

Arguments:
- `input_nas`: Path to input NAS file (volume elements only)
- `output_nas`: Path to output NAS file (surfaces + volumes)
- `verbose`: Print progress messages (default: true)

Returns the output file path.

Note: When `include_internal_interfaces=true` (default), faces shared by two tets with
different PIDs are added once per region (duplicated in NAS with different PIDs), but
they will be merged into a single shared face in the resulting STEP by the converter.
"""
function extract_boundary_surfaces(input_nas::AbstractString, output_nas::AbstractString;
    verbose::Bool=true, include_internal_interfaces::Bool=true)
    verbose && println("Extracting boundary surfaces from volume mesh...")
    verbose && println("  Input:  $input_nas")
    verbose && println("  Output: $output_nas")

    # --- Efficient parsing helpers ---
    # These helpers avoid allocating temporary arrays from `split` for every
    # line in large NAS files. They use `SubString` slices and manual comma
    # index scanning to reduce garbage and speed up parsing. Each helper
    # returns `nothing` on parse failure so the code falls back to ignoring
    # malformed lines instead of throwing.
    # Reuse an index buffer to avoid reallocating per line
    # Try to parse a GRID line (returns (node_id, x,y,z) or nothing)
    function tryparse_grid_line(line::String, idxs::Vector{Int})
        # Expect fields comma-separated; avoid allocation of split when possible
        if !startswith(line, "GRID")
            return nothing
        end
        # Find commas positions quickly
        empty!(idxs)
        for (i, c) in enumerate(line)
            if c == ','
                push!(idxs, i)
            end
        end
        # Need at least 6 parts -> at least 5 commas
        if length(idxs) < 5
            return nothing
        end
        # Extract substrings without copying via SubString
        # parts: 1: GRID, 2: id, 3: CP (maybe empty), 4: X, 5: Y, 6: Z
        # compute ranges
        s = line
        # node id between first and second comma
        id_start = idxs[1] + 1
        id_end = idxs[2] - 1
        x_start = idxs[3] + 1
        x_end = idxs[4] - 1
        y_start = idxs[4] + 1
        y_end = idxs[5] - 1
        z_start = idxs[5] + 1
        # z goes until next comma or end
        nextz = findnext(',', s, z_start)
        z_end = nextz === nothing ? lastindex(s) : nextz - 1
        try
            node_id = parse(Int, SubString(s, id_start:id_end))
            x = parse(Float64, SubString(s, x_start:x_end))
            y = parse(Float64, SubString(s, y_start:y_end))
            z = parse(Float64, SubString(s, z_start:z_end))
            return (node_id, x, y, z)
        catch
            return nothing
        end
    end

    # Try to parse a CTETRA line (returns (eid,pid,n1,n2,n3,n4) or nothing)
    function tryparse_ctetra_line(line::String, idxs::Vector{Int})
        if !startswith(line, "CTETRA")
            return nothing
        end
        # gather comma positions; CTETRA lines typically have at least 7 parts
        empty!(idxs)
        for (i, c) in enumerate(line)
            if c == ','
                push!(idxs, i)
            end
        end
        if length(idxs) < 6
            return nothing
        end
        s = line
        # field ranges: 2: eid, 3: pid, 4:n1,5:n2,6:n3,7:n4
        ranges = Tuple{Int,Int}[]
        for k in 1:6
            st = idxs[k] + 1
            en = (k < length(idxs)) ? idxs[k+1] - 1 : lastindex(s)
            push!(ranges, (st, en))
        end
        try
            eid = parse(Int, SubString(s, ranges[1][1]:ranges[1][2]))
            pid = parse(Int, SubString(s, ranges[2][1]:ranges[2][2]))
            n1 = parse(Int, SubString(s, ranges[3][1]:ranges[3][2]))
            n2 = parse(Int, SubString(s, ranges[4][1]:ranges[4][2]))
            n3 = parse(Int, SubString(s, ranges[5][1]:ranges[5][2]))
            n4 = parse(Int, SubString(s, ranges[6][1]:ranges[6][2]))
            return (eid, pid, n1, n2, n3, n4)
        catch
            return nothing
        end
    end

    # Single-pass parse of the file: GRID and CTETRA in one sweep
    idxbuf = Int[]
    grids = Dict{Int,String}()
    grid_coords = Dict{Int,Tuple{Float64,Float64,Float64}}()
    face_count = Dict{Tuple{Int,Int,Int},Int}()
    face_to_elem = Dict{Tuple{Int,Int,Int},Vector{Tuple{Int,Int,Tuple{Int,Int,Int,Int}}}}()
    elements = String[]

    for line in eachline(input_nas)
        # Try GRID first
        if (res = tryparse_grid_line(line, idxbuf)) !== nothing
            node_id, x, y, z = res
            grids[node_id] = line
            grid_coords[node_id] = (x, y, z)
            continue
        end
        # Then CTETRA
        if (res2 = tryparse_ctetra_line(line, idxbuf)) !== nothing
            eid, pid, n1, n2, n3, n4 = res2
            push!(elements, line)
            faces = _get_tet_faces(n1, n2, n3, n4)
            for face in faces
                face_count[face] = get(face_count, face, 0) + 1
                if !haskey(face_to_elem, face)
                    face_to_elem[face] = Vector{Tuple{Int,Int,Tuple{Int,Int,Int,Int}}}()
                end
                push!(face_to_elem[face], (eid, pid, (n1, n2, n3, n4)))
            end
        end
    end
    verbose && println("  Found $(length(grids)) GRID nodes")

    verbose && println("  Found $(length(elements)) CTETRA elements")
    verbose && println("  Enumerated $(length(face_count)) unique faces")

    # Collect oriented faces to write: external boundaries, plus (optionally) internal interfaces per PID
    oriented_faces = Vector{NamedTuple{(:nodes,:pid),Tuple{NTuple{3,Int},Int}}}()

    ext_count = 0
    int_pairs = 0

    for (face, count) in face_count
        n1, n2, n3 = face
        # external face: appears once
        if count == 1
            adj = face_to_elem[face][1]
            eid, pid, tet_nodes = adj

            # coordinates
            if haskey(grid_coords, n1) && haskey(grid_coords, n2) && haskey(grid_coords, n3)
                p1 = grid_coords[n1]; p2 = grid_coords[n2]; p3 = grid_coords[n3]
                p1t = grid_coords[tet_nodes[1]]; p2t = grid_coords[tet_nodes[2]]; p3t = grid_coords[tet_nodes[3]]; p4t = grid_coords[tet_nodes[4]]
                tet_centroid = _compute_tet_centroid(p1t, p2t, p3t, p4t)
                face_normal = _compute_face_normal(p1, p2, p3)
                face_centroid = _compute_face_centroid(p1, p2, p3)
                outward = (face_centroid[1]-tet_centroid[1], face_centroid[2]-tet_centroid[2], face_centroid[3]-tet_centroid[3])
                nodes = _dot(face_normal, outward) > 0 ? (n1,n2,n3) : (n1,n3,n2)
                push!(oriented_faces, (nodes=nodes, pid=pid))
                ext_count += 1
            end
        elseif include_internal_interfaces
            # internal interface: appears with multiple adjacent tets; include pairs with different PIDs
            adjs = face_to_elem[face]
            # group by PID
            by_pid = Dict{Int,Vector{Tuple{Int,Int,NTuple{4,Int}}}}()
            for t in adjs
                pid = t[2]
                if !haskey(by_pid, pid); by_pid[pid] = Vector{Tuple{Int,Int,NTuple{4,Int}}}(); end
                push!(by_pid[pid], t)
            end
            if length(by_pid) >= 2
                # create one oriented face per distinct PID, oriented relative to one representative tet of that PID
                for (pid, tets) in by_pid
                    rep = tets[1]
                    tet_nodes = rep[3]
                    if haskey(grid_coords, n1) && haskey(grid_coords, n2) && haskey(grid_coords, n3)
                        p1 = grid_coords[n1]; p2 = grid_coords[n2]; p3 = grid_coords[n3]
                        p1t = grid_coords[tet_nodes[1]]; p2t = grid_coords[tet_nodes[2]]; p3t = grid_coords[tet_nodes[3]]; p4t = grid_coords[tet_nodes[4]]
                        tet_centroid = _compute_tet_centroid(p1t, p2t, p3t, p4t)
                        face_normal = _compute_face_normal(p1, p2, p3)
                        face_centroid = _compute_face_centroid(p1, p2, p3)
                        outward = (face_centroid[1]-tet_centroid[1], face_centroid[2]-tet_centroid[2], face_centroid[3]-tet_centroid[3])
                        nodes = _dot(face_normal, outward) > 0 ? (n1,n2,n3) : (n1,n3,n2)
                        push!(oriented_faces, (nodes=nodes, pid=pid))
                    end
                end
                int_pairs += 1
            end
        end
    end

    if isempty(oriented_faces)
        @warn "No faces found to write! The mesh might be malformed."
        return nothing
    end

    verbose && println("  Oriented $(length(oriented_faces)) faces (external=$(ext_count), internal-interfaces=$(include_internal_interfaces ? int_pairs : 0))")

    # CRITICAL: For non-conformal meshes, identify and duplicate interface faces
    # Interface faces have same coordinates but different node IDs (different PIDs)
    # We need to include them in BOTH regions to create closed volumes
    
    # Step 1: Group faces by coordinate-based key to find interfaces
    coord_to_faces = Dict{NTuple{3,NTuple{3,Float64}},Vector{NamedTuple{(:nodes,:pid),Tuple{NTuple{3,Int},Int}}}}()
    
    for face in oriented_faces
        n1, n2, n3 = face.nodes
        if haskey(grid_coords, n1) && haskey(grid_coords, n2) && haskey(grid_coords, n3)
            c1 = grid_coords[n1]
            c2 = grid_coords[n2]
            c3 = grid_coords[n3]
            
            # Round coordinates
            c1_round = (round(c1[1], digits=10), round(c1[2], digits=10), round(c1[3], digits=10))
            c2_round = (round(c2[1], digits=10), round(c2[2], digits=10), round(c2[3], digits=10))
            c3_round = (round(c3[1], digits=10), round(c3[2], digits=10), round(c3[3], digits=10))
            
            # Sort coordinates to create canonical key
            coords_sorted = sort([c1_round, c2_round, c3_round])
            coord_key = (coords_sorted[1], coords_sorted[2], coords_sorted[3])
            
            if !haskey(coord_to_faces, coord_key)
                coord_to_faces[coord_key] = Vector{NamedTuple{(:nodes,:pid),Tuple{NTuple{3,Int},Int}}}()
            end
            push!(coord_to_faces[coord_key], face)
        end
    end
    
    # Debug: Check coord_to_faces structure
    if verbose
        total_coord_keys = length(coord_to_faces)
        max_faces_at_coord = maximum([length(v) for v in values(coord_to_faces)]; init=0)
        println("  Debug: coord_to_faces has $(total_coord_keys) unique coordinate keys, max faces at one coord: $(max_faces_at_coord)")
        
        # Sample some entries
        sample_count = 0
        for (k, v) in coord_to_faces
            if length(v) > 1 && sample_count < 3
                pids_at_loc = unique([f.pid for f in v])
                println("    Sample: $(length(v)) faces at coord, PIDs: $(pids_at_loc)")
                sample_count += 1
            end
        end
    end
    
    # Step 2: Identify interface faces and create duplicates for each PID
    final_faces = Vector{NamedTuple{(:nodes,:pid),Tuple{NTuple{3,Int},Int}}}()
    interface_count = 0
    multi_pid_locations = 0
    
    for (coord_key, faces_at_coord) in coord_to_faces
        # Group by PID
        by_pid = Dict{Int,Vector{NamedTuple{(:nodes,:pid),Tuple{NTuple{3,Int},Int}}}}()
        for f in faces_at_coord
            if !haskey(by_pid, f.pid)
                by_pid[f.pid] = Vector{NamedTuple{(:nodes,:pid),Tuple{NTuple{3,Int},Int}}}()
            end
            push!(by_pid[f.pid], f)
        end
        
        if length(by_pid) >= 2
            # Interface between different PIDs!
            # Include one face per PID to ensure each region is closed
            for (pid, faces) in by_pid
                # Take first face for this PID at this coordinate
                push!(final_faces, faces[1])
            end
            interface_count += 1
            multi_pid_locations += length(by_pid)
        elseif length(by_pid) == 1
            # Single PID at this location
            pid = collect(keys(by_pid))[1]
            faces = by_pid[pid]
            # Just include one face (avoid duplicates within same PID)
            push!(final_faces, faces[1])
        end
    end
    
    if verbose && multi_pid_locations > 0
        println("  Debug: Found $(interface_count) coordinate locations with multiple PIDs ($(multi_pid_locations) total faces)")
    end
    
    oriented_faces = final_faces
    
    if interface_count > 0
        verbose && println("  Duplicated $(interface_count) interface faces to create closed volumes")
    end
    verbose && println("  Final face count: $(length(oriented_faces)) (includes $(interface_count) duplicated interfaces)")

    # Write output NAS file
    verbose && println("  Writing output file...")
    open(output_nas, "w") do io
        # Write header
        println(io, "CEND")
        println(io, "BEGIN BULK")
        println(io, "\$ Generated by Nas2Step.extract_boundary_surfaces()")
        println(io, "\$ Original file: $input_nas")
        println(io, "\$ Boundary faces extracted: $(length(oriented_faces))")
        println(io, "")

        # Write GRID nodes
        println(io, "\$ GRID Nodes")
        for (node_id, line) in sort(collect(grids); by=first)
            println(io, line)
        end
        println(io, "")

        # Write CTRIA3 boundary surface elements
        println(io, "\$ Boundary Surface Elements (CTRIA3)")
        surf_eid = 1
        for face in oriented_faces
            println(io, "CTRIA3,$(surf_eid),$(face.pid),$(face.nodes[1]),$(face.nodes[2]),$(face.nodes[3])")
            surf_eid += 1
        end
        println(io, "")

        # Write CTETRA volume elements
        println(io, "\$ Volume Elements (CTETRA)")
        for elem_line in elements
            println(io, elem_line)
        end

        println(io, "ENDDATA")
    end

    verbose && println("✓ Successfully extracted $(length(oriented_faces)) boundary surfaces")

    return output_nas
end
