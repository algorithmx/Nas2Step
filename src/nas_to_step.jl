"""
    nas_to_step(nas_path::AbstractString; step_path=nothing, classify_angle_deg=40.0,
                auto_extract_boundaries=true, emit_anomaly_json::Bool=true,
                anomaly_json_path::Union{Nothing,AbstractString}=nothing,
                anomaly_max_examples::Int=50)

Convert a NASTRAN .nas mesh to a STEP file using Gmsh.jl.

If the NAS file contains only volume elements (CTETRA) without surface elements,
and `auto_extract_boundaries=true`, boundary surfaces will be automatically extracted
before conversion.

Arguments:
- `nas_path`: Input NAS file path
- `step_path`: Output STEP file path (default: replaces .nas with .step)
- `classify_angle_deg`: Feature angle in degrees for surface classification (default: 40)
- `auto_extract_boundaries`: Automatically extract boundaries for volume-only meshes (default: true)

Returns the output STEP path.
"""
function nas_to_step(nas_path::AbstractString; step_path=nothing, classify_angle_deg=40.0,
    auto_extract_boundaries=true, emit_anomaly_json::Bool=true,
    anomaly_json_path::Union{Nothing,AbstractString}=nothing,
    anomaly_max_examples::Int=50)

    nas_path1 = nas_path
    # Check if NAS file has surface elements
    if auto_extract_boundaries # && !has_surface_elements(nas_path)
        # Create temporary file with boundary surfaces
        # Use include_internal_interfaces=true to include interface faces for non-conformal meshes
        # These duplicate interface faces will be merged during OCC geometry creation
        nas_path1 = replace(nas_path, ".nas" => "_with_surfaces_temp.nas")
        result_nas = extract_boundary_surfaces(nas_path, nas_path1; verbose=true, include_internal_interfaces=true)
        # Check if extraction succeeded
        if result_nas === nothing
            error("Failed to extract boundary surfaces from volume mesh. "
                  * "The mesh might be malformed or have no external boundaries. "
                  * "Try inspecting the NAS file to ensure it contains valid CTETRA elements.")
        end
    end

    # Compute default out paths
    out_step = isnothing(step_path) ? string(replace(nas_path1, r"\.[Nn][Aa][Ss]$" => ""), ".step") : String(step_path)
    out_json = nothing
    if emit_anomaly_json
        out_json = isnothing(anomaly_json_path) ? string(replace(out_step, r"\.[Ss][Tt][Ee][Pp]$" => ""), "_anomalies.json") : String(anomaly_json_path)
    end

    # Convert the enhanced file
    result = convert_nas_to_step(nas_path1; outpath=out_step, featureAngleDeg=classify_angle_deg,
        anomaly_json_path=out_json, anomaly_max_examples=anomaly_max_examples)

    # Clean up temporary file
    if auto_extract_boundaries && nas_path1 != nas_path
        rm(nas_path1; force=true)
    end

    return result
end



#############################
# Helper functions (top-level)
#############################

# Return a canonical undirected edge key
edge_key(a::Int, b::Int) = (a < b ? (a, b) : (b, a))

# Return a canonical undirected triangle key
function face_key(a::Int, b::Int, c::Int)
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

# Triangle area and centroid from node_map (NAS node ids -> coords)
function tri_area_and_centroid(node_map::Dict{Int,NTuple{3,Float64}}, n1::Int, n2::Int, n3::Int)
    p1 = node_map[n1]
    p2 = node_map[n2]
    p3 = node_map[n3]
    v1 = (p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3])
    v2 = (p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3])
    cx = v1[2] * v2[3] - v1[3] * v2[2]
    cy = v1[3] * v2[1] - v1[1] * v2[3]
    cz = v1[1] * v2[2] - v1[2] * v2[1]
    area = 0.5 * sqrt(cx * cx + cy * cy + cz * cz)
    cx0 = (p1[1] + p2[1] + p3[1]) / 3.0
    cy0 = (p1[2] + p2[2] + p3[2]) / 3.0
    cz0 = (p1[3] + p2[3] + p3[3]) / 3.0
    return area, (cx0, cy0, cz0)
end

# Discrete model setup and data extraction
function load_discrete_model!(inpath::AbstractString)
    gmsh.model.add("discrete")
    gmsh.open(inpath)
end

function extract_mesh_data()
    # Nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_map = Dict{Int,NTuple{3,Float64}}()
    for (i, tag) in enumerate(node_tags)
        idx = (i - 1) * 3
        node_map[Int(tag)] = (node_coords[idx+1], node_coords[idx+2], node_coords[idx+3])
    end

    # Names for 2D entities
    discrete_entity_names = Dict{Tuple{Int,Int},String}()
    for entity in gmsh.model.getEntities(2)
        dim, tag = entity
        discrete_entity_names[(dim, tag)] = gmsh.model.getEntityName(dim, tag)
    end

    # Triangles per-entity when possible
    entity_triangles = Dict{Tuple{Int,Int},Vector{Tuple{Int64,Int64,Int64}}}()
    for entity in gmsh.model.getEntities(2)
        dim, tag = entity
        elem_types_ent, _elem_tags_ent, node_tags_ent = gmsh.model.mesh.getElements(dim, tag)
        for (i, etype) in enumerate(elem_types_ent)
            if etype == 2 # 3-node triangle
                tags = node_tags_ent[i]
                tris = get!(entity_triangles, (dim, tag), Tuple{Int64,Int64,Int64}[])
                for j in 1:3:length(tags)
                    if j + 2 <= length(tags)
                        push!(tris, (tags[j], tags[j+1], tags[j+2]))
                    end
                end
            end
        end
    end

    # Fallback if no per-entity triangles found
    triangles = Tuple{Int64,Int64,Int64}[]
    if isempty(entity_triangles)
        elem_types, _elem_tags, node_tags_per_elem = gmsh.model.mesh.getElements()
        for (i, etype) in enumerate(elem_types)
            if etype == 2
                tags = node_tags_per_elem[i]
                for j in 1:3:length(tags)
                    if j + 2 <= length(tags)
                        push!(triangles, (tags[j], tags[j+1], tags[j+2]))
                    end
                end
            end
        end
    end

    println("  Extracted $(length(node_map)) mesh nodes and $(isempty(entity_triangles) ? length(triangles) : sum(length(v) for v in values(entity_triangles))) triangles")
    println("  Found $(length(entity_triangles)) distinct regions (entities)")

    return node_map, discrete_entity_names, entity_triangles, triangles
end

# OCC model construction helpers
function build_occ_points(node_map::Dict{Int,NTuple{3,Float64}})
    # Group nodes by rounded coordinate to merge duplicates
    coord_to_nodes = Dict{NTuple{3,Float64},Vector{Int}}()
    for (node_tag, coord) in node_map
        rounded = (round(coord[1], digits=10), round(coord[2], digits=10), round(coord[3], digits=10))
        push!(get!(coord_to_nodes, rounded, Int[]), Int(node_tag))
    end
    occ_node_map = Dict{Int,Int}()
    unique_coords = 0
    for (coord, node_list) in coord_to_nodes
        occ_tag = gmsh.model.occ.addPoint(coord[1], coord[2], coord[3])
        for nas_node in node_list
            occ_node_map[nas_node] = occ_tag
        end
        unique_coords += 1
    end
    duplicate_nodes = length(node_map) - unique_coords
    println("  Created $(unique_coords) OCC points (merged $(duplicate_nodes) duplicate nodes)")
    return occ_node_map, unique_coords, duplicate_nodes
end

function create_face_for_tri!(n1::Int, n2::Int, n3::Int,
    occ_node_map::Dict{Int,Int}, node_map::Dict{Int,NTuple{3,Float64}},
    edge_map::Dict{Tuple{Int,Int},Int}, face_map::Dict{Tuple{Int,Int,Int},Int},
    face_oriented_nodes::Dict{Tuple{Int,Int,Int},NTuple{3,Int}})

    area_eps = 1e-16
    occ_n1 = Int(occ_node_map[n1])
    occ_n2 = Int(occ_node_map[n2])
    occ_n3 = Int(occ_node_map[n3])
    # skip degenerate triangles (merged OCC points)
    if occ_n1 == occ_n2 || occ_n2 == occ_n3 || occ_n3 == occ_n1
        return nothing
    end
    k = face_key(occ_n1, occ_n2, occ_n3)
    if haskey(face_map, k)
        return face_map[k]
    end
    area, _ = tri_area_and_centroid(node_map, n1, n2, n3)
    if !(area > area_eps)
        return nothing
    end
    # ensure shared edges
    e12 = edge_key(occ_n1, occ_n2)
    e23 = edge_key(occ_n2, occ_n3)
    e31 = edge_key(occ_n3, occ_n1)
    if !haskey(edge_map, e12)
        edge_map[e12] = gmsh.model.occ.addLine(e12[1], e12[2])
    end
    if !haskey(edge_map, e23)
        edge_map[e23] = gmsh.model.occ.addLine(e23[1], e23[2])
    end
    if !haskey(edge_map, e31)
        edge_map[e31] = gmsh.model.occ.addLine(e31[1], e31[2])
    end
    # oriented curve loop (respecting triangle order)
    l1 = edge_map[e12]
    l2 = edge_map[e23]
    l3 = edge_map[e31]
    oriented_edges = Int[]
    push!(oriented_edges, (e12[1] == occ_n1 && e12[2] == occ_n2) ? l1 : -l1)
    push!(oriented_edges, (e23[1] == occ_n2 && e23[2] == occ_n3) ? l2 : -l2)
    push!(oriented_edges, (e31[1] == occ_n3 && e31[2] == occ_n1) ? l3 : -l3)
    loop = gmsh.model.occ.addCurveLoop(oriented_edges)
    surf = gmsh.model.occ.addPlaneSurface([loop])
    face_map[k] = surf
    face_oriented_nodes[k] = (occ_n1, occ_n2, occ_n3)
    return surf
end

function build_shared_geometry!(entity_triangles, triangles, node_map, occ_node_map)
    # caches for shared topology (using OCC node IDs)
    edge_map = Dict{Tuple{Int,Int},Int}()
    face_map = Dict{Tuple{Int,Int,Int},Int}()
    face_oriented_nodes = Dict{Tuple{Int,Int,Int},NTuple{3,Int}}()

    if isempty(entity_triangles)
        println("\033[38;5;208m  No entity information found; creating shared surfaces once for all triangles\033[0m")
        for tri in triangles
            try
                create_face_for_tri!(Int(tri[1]), Int(tri[2]), Int(tri[3]), occ_node_map, node_map, edge_map, face_map, face_oriented_nodes)
            catch
                # skip degenerate
            end
        end
    else
        println("\033[36m  Creating shared OCC edges and faces across all regions\033[0m")
        for (_, tris) in entity_triangles
            for tri in tris
                n1, n2, n3 = Int(tri[1]), Int(tri[2]), Int(tri[3])
                try
                    create_face_for_tri!(n1, n2, n3, occ_node_map, node_map, edge_map, face_map, face_oriented_nodes)
                catch
                    # skip degenerate
                end
            end
        end
    end

    return edge_map, face_map, face_oriented_nodes
end

function get_valid_surfaces()
    valid_surfaces = Set{Int}()
    for (dim, tag) in gmsh.model.occ.getEntities(2)
        push!(valid_surfaces, tag)
    end
    return valid_surfaces
end

function construct_regions_and_volumes(entity_triangles, face_map, face_oriented_nodes, occ_node_map, valid_surfaces, discrete_entity_names; anomaly_max_examples::Int=50)
    region_signed_surfs = Dict{Tuple{Int,Int},Vector{Int}}()
    region_volumes = Dict{Tuple{Int,Int},Vector{Int}}()

    # Build unsigned surface sets per region
    if !isempty(entity_triangles)
        for (entity, tris) in entity_triangles
            surf_set = Set{Int}()
            for tri in tris
                nas_nodes = (Int(tri[1]), Int(tri[2]), Int(tri[3]))
                if all(haskey(occ_node_map, n) for n in nas_nodes)
                    occ_n1 = Int(occ_node_map[nas_nodes[1]])
                    occ_n2 = Int(occ_node_map[nas_nodes[2]])
                    occ_n3 = Int(occ_node_map[nas_nodes[3]])
                    k = face_key(occ_n1, occ_n2, occ_n3)
                    if haskey(face_map, k)
                        surf = face_map[k]
                        if surf in valid_surfaces
                            push!(surf_set, surf)
                        end
                    end
                end
            end
            region_signed_surfs[entity] = collect(surf_set)
        end
    end

    # Precompute face->edge sets (with OCC node keys)
    face_edges = Dict{Int,NTuple{3,Tuple{Int,Int}}}()
    for (k, surf) in face_map
        on = face_oriented_nodes[k]
        e1 = edge_key(on[1], on[2])
        e2 = edge_key(on[2], on[3])
        e3 = edge_key(on[3], on[1])
        face_edges[surf] = (e1, e2, e3)
    end

    println("  Constructing closed shells per region with shared faces…")

    # Refresh current valid surfaces before loops
    valid_surfaces_current = get_valid_surfaces()

    anomalies = Vector{Any}()
    for (entity, surfs) in region_signed_surfs
        dim, tag = entity
        unsigned = unique(surfs)
        unsigned = filter(f -> f in valid_surfaces_current, unsigned)

        # Build edge adjacency
        edge_to_faces = Dict{Tuple{Int,Int},Vector{Int}}()
        for f in unsigned
            if !haskey(face_edges, f)
                continue
            end
            for e in face_edges[f]
                push!(get!(edge_to_faces, e, Int[]), f)
            end
        end

        # BFS components by shared edges
        comps = Vector{Vector{Int}}()
        visited = Set{Int}()
        for f in unsigned
            if f in visited
                continue
            end
            queue = [f]
            push!(visited, f)
            comp = Int[]
            while !isempty(queue)
                cur = popfirst!(queue)
                push!(comp, cur)
                for e in get(face_edges, cur, (edge_key(0, 0), edge_key(0, 0), edge_key(0, 0)))
                    for nb in get(edge_to_faces, e, Int[])
                        if !(nb in visited)
                            push!(visited, nb)
                            push!(queue, nb)
                        end
                    end
                end
            end
            push!(comps, comp)
        end

        # For each component, create shell and try to build volume
        for (comp_idx, comp) in enumerate(comps)
            comp_unsigned = Int[]
            invalid_surfs = Int[]
            for f in comp
                if f in valid_surfaces_current
                    push!(comp_unsigned, f)
                else
                    push!(invalid_surfs, f)
                end
            end
            if isempty(comp_unsigned)
                if !isempty(invalid_surfs)
                    println("    ⚠ Region PID=$(tag): component skipped - all $(length(invalid_surfs)) surfaces invalid")
                end
                continue
            end
            if !isempty(invalid_surfs)
                println("    ⚠ Region PID=$(tag): component has $(length(invalid_surfs)) invalid surfaces (kept $(length(comp_unsigned)) valid)")
            end

            # manifoldness check with evidence report
            closed, non2, boundary_edges, nonmanifold_edges, total_edges = analyze_component_manifoldness(comp, face_edges)
            if !closed
                # summarize a few examples
                max_show = 8
                if !isempty(boundary_edges)
                    println("      · Boundary edges (count=$(length(boundary_edges))): ")
                    for (i, (e, faces)) in enumerate(boundary_edges)
                        if i > max_show
                            println("        …")
                            break
                        end
                        println("        edge=$(e) incident_faces=$(faces)")
                    end
                end
                if !isempty(nonmanifold_edges)
                    println("      · Nonmanifold edges (count=$(length(nonmanifold_edges))): ")
                    for (i, (e, faces)) in enumerate(nonmanifold_edges)
                        if i > max_show
                            println("        …")
                            break
                        end
                        println("        edge=$(e) incident_faces=$(faces)")
                    end
                end
            end
            # Collect anomaly record (sample edges up to anomaly_max_examples)
            if !closed
                be_max = min(length(boundary_edges), anomaly_max_examples)
                nme_max = min(length(nonmanifold_edges), anomaly_max_examples)
                be_samples = boundary_edges[1:be_max]
                nme_samples = nonmanifold_edges[1:nme_max]
                entity_name = get(discrete_entity_names, (dim, tag), "")
                push!(anomalies, (
                    region_dim=dim,
                    region_tag=tag,
                    region_name=entity_name,
                    component_index=comp_idx,
                    faces_in_component=length(comp),
                    total_edges=total_edges,
                    boundary_edges_count=length(boundary_edges),
                    nonmanifold_edges_count=length(nonmanifold_edges),
                    boundary_edges=be_samples,
                    nonmanifold_edges=nme_samples,
                ))
            end
            try
                shell_tag = gmsh.model.occ.addSurfaceLoop(comp_unsigned)
                if closed
                    vol_tag = gmsh.model.occ.addVolume([shell_tag])
                    push!(get!(region_volumes, entity, Int[]), vol_tag)
                    println("    ✓ Region PID=$(tag): Created volume $(vol_tag) (component with $(length(comp_unsigned)) faces)")
                else
                    println("    ⚠ Region PID=$(tag): component is NOT closed (non-2 edge usages=$(non2)); keeping OPEN_SHELL ($(length(comp_unsigned)) faces)")
                end
            catch e
                @warn "Could not create shell/volume for region PID=$(tag) component: $e"
                println("    ⚠ Region PID=$(tag): component kept as OPEN_SHELL ($(length(comp_unsigned)) faces)")
            end
        end
    end

    return region_signed_surfs, region_volumes, anomalies
end

# Analyze whether a set of faces (comp) forms a 2-manifold closed shell by counting edge incidences.
# Returns:
# - closed::Bool: true if every edge is used exactly twice
# - non2::Int: number of edges with usage != 2
# - boundary_edges::Vector{Tuple{Tuple{Int,Int},Vector{Int}}}: edges with usage 1 and their incident faces
# - nonmanifold_edges::Vector{Tuple{Tuple{Int,Int},Vector{Int}}}: edges with usage >2 and their incident faces
function analyze_component_manifoldness(comp::Vector{Int}, face_edges::Dict{Int,NTuple{3,Tuple{Int,Int}}})
    edge_to_faces = Dict{Tuple{Int,Int},Vector{Int}}()
    for f in comp
        if !haskey(face_edges, f)
            continue
        end
        for e in face_edges[f]
            push!(get!(edge_to_faces, e, Int[]), f)
        end
    end
    usage = Dict{Tuple{Int,Int},Int}()
    for (e, faces) in edge_to_faces
        usage[e] = length(faces)
    end
    non2 = count(v -> v != 2, values(usage))
    closed = (non2 == 0)

    boundary_edges = Vector{Tuple{Tuple{Int,Int},Vector{Int}}}()
    nonmanifold_edges = Vector{Tuple{Tuple{Int,Int},Vector{Int}}}()
    for (e, faces) in edge_to_faces
        c = length(faces)
        if c == 1
            push!(boundary_edges, (e, faces))
        elseif c > 2
            push!(nonmanifold_edges, (e, faces))
        end
    end
    return closed, non2, boundary_edges, nonmanifold_edges, length(usage)
end

# Write anomalies as JSON to a file without external dependencies.
# The structure is: { "anomalies": [ { ... } ] }
function write_anomalies_json(path::AbstractString, anomalies)
    # local helpers
    json_escape(s::AbstractString) = replace(replace(String(s), '\\' => "\\\\"), '"' => "\\\"")
    function json_array_ints(v::AbstractVector{<:Integer})
        return "[" * join(v, ",") * "]"
    end
    function json_edge_obj(edge::Tuple{Int,Int}, faces::Vector{Int})
        return "{\"edge\":[" * string(edge[1]) * "," * string(edge[2]) * "],\"faces\":" * json_array_ints(faces) * "}"
    end
    open(path, "w") do io
        print(io, "{\"anomalies\":[")
        for (i, a) in enumerate(anomalies)
            if i > 1
                print(io, ",")
            end
            # unpack NamedTuple like fields
            rd = a.region_dim
            rt = a.region_tag
            rn = a.region_name
            ci = a.component_index
            fic = a.faces_in_component
            te = a.total_edges
            bec = a.boundary_edges_count
            nmc = a.nonmanifold_edges_count
            be = a.boundary_edges
            nme = a.nonmanifold_edges
            print(io, "{\"region_dim\":", rd,
                ",\"region_tag\":", rt,
                ",\"region_name\":\"", json_escape(rn), "\",")
            print(io, "\"component_index\":", ci,
                ",\"faces_in_component\":", fic,
                ",\"total_edges\":", te,
                ",\"boundary_edges_count\":", bec,
                ",\"nonmanifold_edges_count\":", nmc, ",")
            # boundary edges
            print(io, "\"boundary_edges\":[")
            for (j, item) in enumerate(be)
                if j > 1
                    print(io, ",")
                end
                print(io, json_edge_obj(item[1], item[2]))
            end
            print(io, "],")
            # nonmanifold edges
            print(io, "\"nonmanifold_edges\":[")
            for (j, item) in enumerate(nme)
                if j > 1
                    print(io, ",")
                end
                print(io, json_edge_obj(item[1], item[2]))
            end
            print(io, "]}")
        end
        print(io, "]}")
    end
end

function name_regions!(entity_triangles, region_volumes, region_signed_surfs, discrete_entity_names)
    if isempty(entity_triangles)
        return
    end
    region_id = 1
    for (entity, _) in entity_triangles
        dim, tag = entity
        entity_name = get(discrete_entity_names, entity, "")
        if isempty(entity_name)
            entity_name = "R$(region_id)"
        end
        if haskey(region_volumes, entity) && !isempty(region_volumes[entity])
            vols = region_volumes[entity]
            for (i, vt) in enumerate(vols)
                local_name = length(vols) == 1 ? entity_name : string(entity_name, "_", i)
                gmsh.model.setEntityName(3, vt, local_name)
            end
            phys_tag = gmsh.model.addPhysicalGroup(3, vols)
            gmsh.model.setPhysicalName(3, phys_tag, entity_name)
            println("    ✓ Named $(length(vols)) volume(s) as '$(entity_name)' (PID=$(tag))")
        else
            surfs = get(region_signed_surfs, entity, Int[])
            unsigned_surfs = unique(abs.(surfs))
            for s in unsigned_surfs
                gmsh.model.setEntityName(2, s, entity_name)
            end
            phys_tag = gmsh.model.addPhysicalGroup(2, unsigned_surfs)
            gmsh.model.setPhysicalName(2, phys_tag, entity_name)
            println("    ⚠ Named $(length(unsigned_surfs)) shared surfaces as '$(entity_name)' (PID=$(tag)) [no volume]")
        end
        region_id += 1
    end
end

"""
    convert_nas_to_step(inpath; outpath=nothing, featureAngleDeg=40, sewFaces=true, tol=1e-5)

Convert a NASTRAN mesh into a STEP file using a shared-topology OCC reconstruction.

This function orchestrates the workflow and delegates each step to small helpers
to keep logic clear and maintainable.
"""
function convert_nas_to_step(inpath::AbstractString; outpath::Union{Nothing,AbstractString}=nothing,
    featureAngleDeg::Real=40, sewFaces::Bool=true, tol::Real=1e-5,
    anomaly_json_path::Union{Nothing,AbstractString}=nothing, anomaly_max_examples::Int=50)

    # all helpers are now defined at top-level

    out = isnothing(outpath) ? string(replace(inpath, r"\.[Nn][Aa][Ss]$" => ""), ".step") : String(outpath)

    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)  # 0=silent, 1=errors, 2=warnings, 3=info
        gmsh.option.setNumber("Geometry.OCCSewFaces", sewFaces ? 1 : 0)
        gmsh.option.setNumber("Geometry.Tolerance", tol)

        # Step 1: Read mesh in "discrete" model and extract data
        load_discrete_model!(inpath)
        node_map, discrete_entity_names, entity_triangles, triangles = extract_mesh_data()

        # Step 2: Build OCC model (shared topology)
        gmsh.model.add("occ")
        occ_node_map, unique_coords, duplicate_nodes = build_occ_points(node_map)
        edge_map, face_map, face_oriented_nodes = build_shared_geometry!(entity_triangles, triangles, node_map, occ_node_map)
        gmsh.model.occ.synchronize()

        # Summary of created surfaces
        valid_surfaces = get_valid_surfaces()
        total_surfaces = length(valid_surfaces)
        println("  Created $(total_surfaces) unique OCC surfaces (duplicate nodes already merged)")

        # Step 3: Build regions, shells, and volumes (if possible)
        region_signed_surfs = Dict{Tuple{Int,Int},Vector{Int}}()
        region_volumes = Dict{Tuple{Int,Int},Vector{Int}}()
        anomalies = Any[]
        if !isempty(entity_triangles)
            region_signed_surfs, region_volumes, anomalies = construct_regions_and_volumes(
                entity_triangles, face_map, face_oriented_nodes, occ_node_map, valid_surfaces, discrete_entity_names;
                anomaly_max_examples=anomaly_max_examples)
            gmsh.model.occ.synchronize()
            if !isempty(region_volumes)
                total_vols = sum(length(v) for v in values(region_volumes))
                println("  Successfully created $(total_vols) volumetric components across $(length(region_volumes)) regions with shared boundaries")
            else
                println("  ⚠ Warning: Volumes not created - STEP will contain OPEN_SHELL surfaces (still shared)")
            end
            if anomaly_json_path !== nothing && !isempty(anomalies)
                try
                    write_anomalies_json(String(anomaly_json_path), anomalies)
                    println("  ⚠ Wrote manifoldness anomalies JSON to '$(anomaly_json_path)' (entries=$(length(anomalies)))")
                catch e
                    @warn "Failed to write anomalies JSON: $e"
                end
            end
        end

        # Step 4: Name entities
        name_regions!(entity_triangles, region_volumes, region_signed_surfs, discrete_entity_names)

        # Step 5: Write STEP file
        gmsh.write(out)
        return out
    finally
        gmsh.finalize()
    end
end
