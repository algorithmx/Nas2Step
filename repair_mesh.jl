#!/usr/bin/env julia
"""
Mesh Conformity Repair (straightforward version)

Repairs non-conforming interfaces in a Nastran mesh using the planning
APIs in Nas2Step and applies quad retriangulation flips directly.

Usage: julia repair_mesh.jl input.nas output.nas
"""

using Nas2Step

using Gmsh
using Printf

println("="^70)
println("Mesh Conformity Repair")
println("="^70)

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
if length(ARGS) < 2
    println("\nUsage: julia repair_mesh.jl input.nas output.nas")
    println("Example: julia repair_mesh.jl examples/realistic/NC_Reduction_4.nas repaired.nas")
    exit(1)
end

const INPUT_FILE = ARGS[1]
const OUTPUT_FILE = ARGS[2]

if !isfile(INPUT_FILE)
    error("Input file not found: $INPUT_FILE")
end

println("\nInput : $INPUT_FILE")
println("Output: $OUTPUT_FILE")
println("-"^70)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

"""Round a coordinate triple for robust dict keys."""
ckey(p::NTuple{3,Float64}) = (round(p[1]; digits=4), round(p[2]; digits=4), round(p[3]; digits=4))

"""Build a lookup from rounded coords to node id for a NAS file (single gmsh open)."""
function build_node_lookup(nas_file::AbstractString)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.open(nas_file)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        lut = Dict{NTuple{3,Float64},Int}()
        for (i, tag) in enumerate(tags)
            j = (i-1)*3
            p = (coords[j+1], coords[j+2], coords[j+3])
            lut[ckey(p)] = Int(tag)
        end
        return lut
    finally
        gmsh.finalize()
    end
end

"""Find all volume tag pairs that share a meaningful set of nodes (potential interfaces)."""
function find_interface_pairs(nas_file::AbstractString; min_shared::Int=10)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.open(nas_file)
        vols = [tag for (dim, tag) in gmsh.model.getEntities(3)]
        if length(vols) < 2
            return Tuple{Int,Int}[]
        end

        # Build a global node tag -> coordinate lookup once
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        tag_to_coord = Dict{Int,NTuple{3,Float64}}()
        for (i, t) in enumerate(node_tags)
            base = (i-1)*3
            if base + 3 <= length(node_coords)
                tag_to_coord[Int(t)] = (node_coords[base+1], node_coords[base+2], node_coords[base+3])
            end
        end

        # Nodes by volume
        vol_nodes = Dict{Int,Set{NTuple{3,Float64}}}()
        for v in vols
            etypes, _, enodes = gmsh.model.mesh.getElements(3, v)
            nodes_here = Set{NTuple{3,Float64}}()
            for (t, ev) in zip(etypes, enodes)
                # 4 = tetra (Gmsh type)
                if t != 4; continue; end
                for i in 1:4:length(ev)
                    # all node ids of the tet
                    # We just record their coords as keys
                    for k in 0:3
                        nid = Int(ev[i+k])
                        coord = get(tag_to_coord, nid, nothing)
                        if coord === nothing
                            continue
                        end
                        push!(nodes_here, ckey(coord))
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

"""Parse a CTRIA3 line; return (eid, pid, n1, n2, n3) or nothing."""
function parse_ctria3(line::AbstractString)
    s = strip(line)
    startswith(s, "CTRIA3") || return nothing
    tok = split(s)
    length(tok) < 6 && return nothing
    try
        return (
            eid = parse(Int, tok[2]),
            pid = parse(Int, tok[3]),
            n1 = parse(Int, tok[4]),
            n2 = parse(Int, tok[5]),
            n3 = parse(Int, tok[6])
        )
    catch
        return nothing
    end
end

"""Batch-apply quad retriangulation flips into a new NAS file."""
function apply_quad_flips_batch!(nas_file::AbstractString,
                                output_file::AbstractString,
                                target_pid::Int,
                                flips::Vector{NamedTuple{(:old_tris,:new_tris),Tuple{Vector{NTuple{3,Int}},Vector{NTuple{3,Int}}}}})
    lines = readlines(nas_file)

    # Build a set of triangles to delete (unordered nodes)
    dels = Set{NTuple{3,Int}}()
    for fl in flips
        for tri in fl.old_tris
            vec = [tri[1], tri[2], tri[3]]
            sort!(vec)
            nsorted = (vec[1], vec[2], vec[3])
            push!(dels, nsorted)
        end
    end

    kept = String[]
    deleted = 0
    for ln in lines
        info = parse_ctria3(ln)
        if info === nothing
            push!(kept, ln)
            continue
        end
        nodes_sorted = sort((info.n1, info.n2, info.n3))
        if info.pid == target_pid && (nodes_sorted in dels)
            deleted += 1
            continue
        end
        push!(kept, ln)
    end

    enddata_idx = findfirst(l -> startswith(l, "\$ENDDATA") || startswith(l, "ENDDATA"), kept)
    body = enddata_idx === nothing ? kept : kept[1:enddata_idx-1]

    # Find max CTRIA3 element id among kept lines
    max_eid = 0
    for ln in body
        info = parse_ctria3(ln)
        if info !== nothing
            max_eid = max(max_eid, info.eid)
        end
    end

    # Append new triangles
    next_eid = max_eid + 1
    for fl in flips
        for tri in fl.new_tris
            push!(body, @sprintf("CTRIA3  %8d%8d%8d%8d%8d", next_eid, target_pid, tri[1], tri[2], tri[3]))
            next_eid += 1
        end
    end

    # Reconstruct final text
    out = enddata_idx === nothing ? body : vcat(body, kept[enddata_idx:end])
    open(output_file, "w") do io
        for ln in out
            println(io, ln)
        end
    end

    return deleted
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

println("Analyzing interfaces...")
pairs = find_interface_pairs(INPUT_FILE)
if isempty(pairs)
    println("No interfaces detected. Nothing to do.")
    cp(INPUT_FILE, OUTPUT_FILE; force=true)
    exit(0)
end

for (i, (a, b)) in enumerate(pairs)
    println(@sprintf("  %2d) PID %d ↔ PID %d", i, a, b))
end

current_file = INPUT_FILE
total_applied = 0

for (idx, (pidA, pidB)) in enumerate(pairs)
    # Ensure we refer to and update the global bindings in this soft-scope loop
    global current_file
    global total_applied
    # Build topology first to allow early skip for perfect interfaces
    topo = build_interface_topology(current_file, pidA, pidB)
    if isapprox(topo.conformity_ratio, 1.0; atol=1e-12)
        # For 100% conformity, print nothing and skip work for this pair
        continue
    end

    println("\n" * "="^70)
    println("Interface $(idx)/$(length(pairs)): PID $pidA ↔ PID $pidB")
    println("-"^70)
    println(@sprintf("  Conformity before: %.2f%%", topo.conformity_ratio*100))
    cls = classify_interface_mismatches(topo)
    cons = build_boundary_constraints(current_file, pidA, pidB)
    plan = generate_repair_plan(topo, cls, cons)

    feas = count(p -> p.is_feasible, plan.insertion_sequence)
    println("  Feasible plans: $feas/$(length(plan.insertion_sequence))")

    if feas == 0
        println("  → No feasible repairs for this interface. Skipping.")
        continue
    end

    # We only apply quad retriangulation flips in this script
    target_pid = plan.repair_direction === :subdivide_A ? pidA : pidB

    # Build node lookup once per current file
    lut = build_node_lookup(current_file)
    function tri_coords_to_ids(tri_flat::NTuple{9,Float64})
        p1 = ckey((tri_flat[1], tri_flat[2], tri_flat[3]))
        p2 = ckey((tri_flat[4], tri_flat[5], tri_flat[6]))
        p3 = ckey((tri_flat[7], tri_flat[8], tri_flat[9]))
        n1 = get(lut, p1, 0); n2 = get(lut, p2, 0); n3 = get(lut, p3, 0)
        return n1 > 0 && n2 > 0 && n3 > 0 ? (n1, n2, n3) : nothing
    end

    flips = NamedTuple{(:old_tris,:new_tris),Tuple{Vector{NTuple{3,Int}},Vector{NTuple{3,Int}}}}[]
    applied_here = 0

    for p in plan.insertion_sequence
        if !p.is_feasible || p.split_type != :quad_retriangulation
            continue
        end

        # Expect 2 old and 2 replacement triangles for a quad flip
        if length(p.old_triangles) < 2 || length(p.replacement_triangles) < 2
            continue
        end

        old1 = tri_coords_to_ids(p.old_triangles[1])
        old2 = tri_coords_to_ids(p.old_triangles[2])
        new1 = tri_coords_to_ids(p.replacement_triangles[1])
        new2 = tri_coords_to_ids(p.replacement_triangles[2])

        if any(x -> x === nothing, (old1, old2, new1, new2))
            @warn "  Skipping a plan due to missing node IDs in lookup"
            continue
        end

        push!(flips, (old_tris = [old1::NTuple{3,Int}, old2::NTuple{3,Int}],
                       new_tris = [new1::NTuple{3,Int}, new2::NTuple{3,Int}]))
        applied_here += 1
    end

    if isempty(flips)
        println("  → No applicable quad flips in this plan. Skipping.")
        continue
    end

    tmp_out = @sprintf("temp_repair_%03d.nas", idx)
    deleted = apply_quad_flips_batch!(current_file, tmp_out, target_pid, flips)
    if deleted == 0
        println("  → No matching triangles were found to delete. Skipping output swap.")
        rm(tmp_out; force=true)
        continue
    end

    current_file = tmp_out
    total_applied += applied_here

    # Quick post-check for this interface on the updated file
    topo_after = build_interface_topology(current_file, pidA, pidB)
    println(@sprintf("  Conformity after : %.2f%%", topo_after.conformity_ratio*100))
end

println("\n" * "="^70)
if total_applied > 0
    cp(current_file, OUTPUT_FILE; force=true)
    # Clean up temps
    for f in readdir(".")
        startswith(f, "temp_repair_") && endswith(f, ".nas") && rm(f; force=true)
    end
    println("Applied $total_applied quad flip(s).")
    println("Wrote repaired mesh to: $OUTPUT_FILE")
else
    println("No repairs applied. Copying input to output unchanged.")
    cp(INPUT_FILE, OUTPUT_FILE; force=true)
end

