"""
    check_nas_volumes(nas_path::AbstractString; angle_deg::Real=40.0,
                       keep_temp::Bool=false) -> Dict

Concise volume check for a NASTRAN .nas file containing CTETRA elements.

What it does:
- Splits CTETRA by PID into temporary NAS files (nodes are filtered per PID).
- For each PID file, opens it in Gmsh, attempts surface classification and
  geometry creation, and computes a simple bounding box and total tetra volume.
- Reports if a volume appears to have no boundary surfaces after geometry
  creation (a common cause of conversion failures).

Returns a Dict with keys:
- :pids :: Vector{Int}
- :reports :: Vector{Dict} (one per PID) with keys:
    :pid, :num_nodes, :num_elements, :bbox => (xmin,ymin,zmin,xmax,ymax,zmax),
    :total_volume, :has_boundary_error, :temp_file

Set keep_temp=true to retain the split NAS files for manual inspection.
"""
function check_nas_volumes(nas_path::AbstractString; angle_deg::Real=40.0,
    keep_temp::Bool=false)

    isfile(nas_path) || error("NAS file not found: $(nas_path)")

    # --- Minimal parse: collect GRID lines and CTETRA grouped by PID ---
    grids = Dict{Int,String}()
    pid_elements = Dict{Int,Vector{String}}()
    pid_nodes = Dict{Int,Set{Int}}()

    for line in eachline(nas_path)
        if startswith(line, "GRID")
            parts = split(line, ',')
            if length(parts) >= 2
                try
                    nid = parse(Int, parts[2])
                    grids[nid] = line
                catch
                end
            end
        elseif startswith(line, "CTETRA")
            parts = split(line, ',')
            if length(parts) >= 7
                try
                    pid = parse(Int, parts[3])
                    n1 = parse(Int, parts[4])
                    n2 = parse(Int, parts[5])
                    n3 = parse(Int, parts[6])
                    n4 = parse(Int, parts[7])
                    push!(get!(pid_elements, pid, String[]), line)
                    nodeset = get!(pid_nodes, pid, Set{Int}())
                    push!(nodeset, n1); push!(nodeset, n2); push!(nodeset, n3); push!(nodeset, n4)
                catch
                end
            end
        end
    end

    pids = sort(collect(keys(pid_elements)))
    if isempty(pids)
        return Dict(:pids => Int[], :reports => Vector{Dict}())
    end

    # Write temp split files
    tmpdir = mktempdir()
    temp_files = Dict{Int,String}()
    for pid in pids
        path = joinpath(tmpdir, "volume_$(pid).nas")
        open(path, "w") do io
            println(io, "CEND")
            println(io, "BEGIN BULK")
            println(io, "\$ Extracted Volume PID=$(pid)")
            println(io, "")
            for nid in sort(collect(pid_nodes[pid]))
                if haskey(grids, nid)
                    println(io, grids[nid])
                end
            end
            println(io, "")
            for elem_line in pid_elements[pid]
                println(io, elem_line)
            end
            println(io, "ENDDATA")
        end
        temp_files[pid] = path
    end

    # --- Helper: tet volume ---
    calculate_tetrahedron_volume(p1, p2, p3, p4) = begin
        v1 = (p2[1]-p1[1], p2[2]-p1[2], p2[3]-p1[3])
        v2 = (p3[1]-p1[1], p3[2]-p1[2], p3[3]-p1[3])
        v3 = (p4[1]-p1[1], p4[2]-p1[2], p4[3]-p1[3])
        det = v1[1]*(v2[2]*v3[3]-v2[3]*v3[2]) -
              v1[2]*(v2[1]*v3[3]-v2[3]*v3[1]) +
              v1[3]*(v2[1]*v3[2]-v2[2]*v3[1])
        abs(det) / 6.0
    end

    # --- Inspect each PID file with Gmsh ---
    reports = Dict[]
    for pid in pids
        f = temp_files[pid]
        info = Dict{Symbol,Any}(:pid => pid, :temp_file => f)
        try
            Gmsh.gmsh.initialize()
            try
                Gmsh.gmsh.option.setNumber("General.Terminal", 0)
                Gmsh.gmsh.model.add("inspect_vol")
                Gmsh.gmsh.open(f)

                node_tags, node_coords, _ = Gmsh.gmsh.model.mesh.getNodes()
                info[:num_nodes] = length(node_tags)
                etypes, etags, _ = Gmsh.gmsh.model.mesh.getElements()
                info[:num_elements] = sum(length.(etags))

                # bbox
                if !isempty(node_tags)
                    xs = [node_coords[3i-2] for i in 1:length(node_tags)]
                    ys = [node_coords[3i-1] for i in 1:length(node_tags)]
                    zs = [node_coords[3i]   for i in 1:length(node_tags)]
                    info[:bbox] = (minimum(xs), minimum(ys), minimum(zs),
                                   maximum(xs), maximum(ys), maximum(zs))
                end

                # classify and create geometry to detect boundary problems
                has_boundary_error = false
                try
                    Gmsh.gmsh.model.mesh.classifySurfaces(angle_deg*pi/180, true, true, true)
                catch
                    # ignore classification errors; we'll still try geometry
                end
                try
                    Gmsh.gmsh.model.mesh.createGeometry()
                    vols = Gmsh.gmsh.model.getEntities(3)
                    if !isempty(vols)
                        for v in vols
                            try
                                bnd = Gmsh.gmsh.model.getBoundary([v], false, false, false)
                                if isempty(bnd)
                                    has_boundary_error = true
                                    break
                                end
                            catch
                                has_boundary_error = true
                                break
                            end
                        end
                    end
                catch e
                    msg = sprint(showerror, e)
                    if occursin("no surfaces on their boundary", msg) || occursin("Discrete volume", msg)
                        has_boundary_error = true
                    end
                end
                info[:has_boundary_error] = has_boundary_error

                # volume sum from tets
                total_vol = 0.0
                if !isempty(node_tags)
                    # map coords
                    coord_map = Dict{Int,NTuple{3,Float64}}()
                    for (i, t) in enumerate(node_tags)
                        coord_map[Int(t)] = (node_coords[3i-2], node_coords[3i-1], node_coords[3i])
                    end
                    etypes, _etags, nodelists = Gmsh.gmsh.model.mesh.getElements()
                    for (i, et) in enumerate(etypes)
                        if et == 4 # 4-node tetra
                            tags = nodelists[i]
                            for j in 1:4:length(tags)
                                if j+3 <= length(tags)
                                    n1 = Int(tags[j]); n2 = Int(tags[j+1]); n3 = Int(tags[j+2]); n4 = Int(tags[j+3])
                                    if haskey(coord_map, n1) && haskey(coord_map, n2) && haskey(coord_map, n3) && haskey(coord_map, n4)
                                        total_vol += calculate_tetrahedron_volume(coord_map[n1], coord_map[n2], coord_map[n3], coord_map[n4])
                                    end
                                end
                            end
                        end
                    end
                end
                info[:total_volume] = total_vol
            finally
                Gmsh.gmsh.finalize()
            end
        catch e
            info[:error] = string(e)
        end
        push!(reports, info)
    end

    # Cleanup temp unless retained
    if !keep_temp
        try
            for f in values(temp_files)
                isfile(f) && rm(f; force=true)
            end
            isdir(tmpdir) && rm(tmpdir; force=true, recursive=true)
        catch
        end
    end

    return Dict(:pids => pids, :reports => reports)
end

"""
    test_closed_shells(nas_path::AbstractString; step_path::Union{Nothing,AbstractString}=nothing,
                        verbose::Bool=true) -> Dict

Converts a NAS file to STEP using `nas_to_step`, then verifies that CLOSED_SHELL
entities are present and that Gmsh detects volume entities. This consolidates the
"closed shell" test previously provided as a standalone script.

Inputs:
- nas_path: Path to input .nas file
- step_path: Optional output path for the STEP file; defaults to `<nas>_closed_shell.step`
- verbose: Whether to print a human-friendly report

Returns a Dict with keys:
- :input_nas, :output_step, :file_size_bytes, :file_size_mb
- :step_counts => Dict(:CLOSED_SHELL, :OPEN_SHELL, :MANIFOLD_SOLID_BREP)
- :gmsh_counts => Dict(:volumes, :surfaces, :curves, :points)
- :success_closed_shell :: Bool
- :success_volumes_detected :: Bool
- :success :: Bool (both of the above)
"""
function test_closed_shells(nas_path::AbstractString; step_path::Union{Nothing,AbstractString}=nothing,
    verbose::Bool=true)

    isfile(nas_path) || error("NAS file not found: $(nas_path)")

    # Helper to count substring occurrences in a large text
    count_substr(text::AbstractString, sub::AbstractString) = begin
        n = 0
        i = firstindex(text)
        while true
            r = findnext(sub, text, i)
            r === nothing && break
            n += 1
            i = last(r) + 1
        end
        n
    end

    # Determine output path
    out_step = isnothing(step_path) ? replace(nas_path, ".nas" => "_closed_shell.step") : String(step_path)

    if verbose
        println("="^70)
        println("Testing CLOSED_SHELL Construction in STEP Conversion")
        println("="^70)
        println("\n[1] Input NAS file: $nas_path")
        if has_surface_elements(nas_path)
            println("    ✓ File has surface elements (CTRIA3/CQUAD4)")
        else
            println("    ⚠ Volume-only file (will auto-extract boundaries)")
        end
        println("\n[2] Converting NAS to STEP with closed shell construction...")
    end

    step_file = nas_to_step(nas_path; step_path=out_step)

    # Read STEP file to count entity labels
    step_text = read(step_file, String)
    closed_shell_count = count_substr(step_text, "CLOSED_SHELL")
    open_shell_count   = count_substr(step_text, "OPEN_SHELL")
    manifold_solid_count = count_substr(step_text, "MANIFOLD_SOLID_BREP")

    if verbose
        println("\n[3] Verifying STEP file structure...")
        println("\n" * "="^70)
        println("STEP Structure Analysis")
        println("="^70)
        println("  CLOSED_SHELL entities:        $(closed_shell_count)")
        println("  OPEN_SHELL entities:          $(open_shell_count)")
        println("  MANIFOLD_SOLID_BREP entities: $(manifold_solid_count)")
        if closed_shell_count > 0
            println("\n  ✓ SUCCESS: STEP file contains CLOSED_SHELL entities!")
            println("            This means regions are now represented as closed volumes.")
        else
            println("\n  ✗ FAILURE: STEP file still only contains OPEN_SHELL entities.")
            println("            Closed shell construction may have failed.")
        end
        println("\n[4] Loading STEP file in Gmsh for verification...")
    end

    info = verify_step(step_file)
    counts = info[:counts]

    if verbose
        println("\n" * "="^70)
        println("Gmsh Entity Verification")
        println("="^70)
        println("  Volumes:  $(counts[:volumes])")
        println("  Surfaces: $(counts[:surfaces])")
        println("  Curves:   $(counts[:curves])")
        println("  Points:   $(counts[:points])")
        if counts[:volumes] > 0
            println("\n  ✓ SUCCESS: Gmsh detected $(counts[:volumes]) volume entities!")
            println("            Region information is preserved as volumetric structures.")
        else
            println("\n  ⚠ WARNING: No volume entities detected by Gmsh.")
            println("            The STEP file may only contain surface geometry.")
        end

        println("\n" * "="^70)
        println("Test Complete!")
        println("="^70)
        sz = stat(step_file).size
        println("\nOutput STEP file: $step_file")
        println("File size: $(round(sz / 1024 / 1024, digits=2)) MB")

        if closed_shell_count > 0 && counts[:volumes] > 0
            println("\n🎉 All tests passed! The STEP file now contains proper volumetric regions.")
        else
            println("\n⚠ Some tests failed. Check the output above for details.")
            println("\nNote: If closed shell construction failed, this likely means:")
            println("  • The boundary surface mesh is not perfectly closed (has gaps)")
            println("  • Surface normals are inconsistent")
            println("  • Duplicate vertices/edges prevent proper shell formation")
        end
    end

    sz_bytes = stat(step_file).size
    result = Dict(
        :input_nas => nas_path,
        :output_step => step_file,
        :file_size_bytes => sz_bytes,
        :file_size_mb => sz_bytes / 1024 / 1024,
        :step_counts => Dict(
            :CLOSED_SHELL => closed_shell_count,
            :OPEN_SHELL => open_shell_count,
            :MANIFOLD_SOLID_BREP => manifold_solid_count,
        ),
        :gmsh_counts => counts,
        :success_closed_shell => closed_shell_count > 0,
        :success_volumes_detected => counts[:volumes] > 0,
    )
    result[:success] = result[:success_closed_shell] && result[:success_volumes_detected]
    return result
end
