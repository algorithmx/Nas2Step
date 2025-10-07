#!/usr/bin/env julia
"""
Mesh Conformity Repair (symmetric unified mesh)

Repairs non-conforming interfaces in a Nastran mesh using the new symmetric
pipeline in Nas2Step:
    topology → symmetric classification → strategy selection → unified mesh → replacement

Usage: julia repair_mesh.jl input.nas output.nas
"""

using Nas2Step

using Gmsh
using Printf
using Dates

# Include JSON reporting infrastructure
include("src/repair/json_reporting_types.jl")
include("src/repair/json_export.jl")

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

"""
Write the current workspace's interface surfaces to a NASTRAN surface file (.nas).
One surface region per PID with CTRIA3 elements. Thickness is set to 1.0 by default.
"""
function write_workspace_surfaces!(ws::Nas2Step.RepairWorkspace, out_path::AbstractString)
    regions = Nas2Step.SurfaceRegion[]
    # Build one region per PID
    for (pid, faces) in ws.working_faces
        # Collect unique points used by this PID's faces
        point_index = Dict{NTuple{3,Float64}, Int}()
        points = NTuple{3,Float64}[]
        triangles = NTuple{3,Int}[]
        for face_nodes in faces
            length(face_nodes) == 3 || continue
            local_ids = Int[]
            for nid in face_nodes
                coord = ws.working_nodes[nid]
                idx = get(point_index, coord, 0)
                if idx == 0
                    push!(points, coord)
                    idx = length(points)
                    point_index[coord] = idx
                end
                push!(local_ids, idx)
            end
            push!(triangles, (local_ids[1], local_ids[2], local_ids[3]))
        end
        isempty(triangles) && continue
        push!(regions, Nas2Step.SurfaceRegion("PID $(pid)", points, triangles, 1.0))
    end
    Nas2Step.write_nas_surface(out_path, regions)
    return out_path
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

# Legacy quad-flip utilities removed in favor of unified symmetric repair

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

"""
Process all detected interfaces in a single workspace with reduced verbosity.
"""

# Create one workspace for the original input mesh
ws = Nas2Step.RepairWorkspace(INPUT_FILE)
total_interfaces_repaired = 0

# Initialize session-level reporting
session_start = time()
interface_reports = InterfaceRepairReport[]

for (idx, (pidA, pidB)) in enumerate(pairs)
    interface_start = time()
    
    # Build topology from the original file (pre-repair) for accurate classification
    topo = build_interface_topology(INPUT_FILE, pidA, pidB)
    if isapprox(topo.conformity_ratio, 1.0; atol=1e-12)
        continue
    end

    println("\n" * "="^50)
    println("Interface $(idx)/$(length(pairs)): PID $pidA ↔ PID $pidB")
    println(@sprintf("  Conformity before: %.2f%%", topo.conformity_ratio*100))

    # Symmetric classification (lower verbosity)
    sym_result = classify_interface_mismatches_symmetric(topo, tol=1e-4, verbose=false)
    sym_mismatches = sym_result.symmetric_mismatches
    isempty(sym_mismatches) && continue

    # Boundary constraints from the original file
    cons = build_boundary_constraints(INPUT_FILE, pidA, pidB)

    # Build symmetric repair plan (reduced verbosity)
    sym_plan = Nas2Step.build_symmetric_repair_plan(
        topo,
        sym_mismatches,
        cons;
        thresholds=default_thresholds(),
        verbose=false
    )

    if !sym_plan.is_feasible
        println("  → Plan not feasible. Skipping.")
        
        # Still collect report for failed interface
        push!(interface_reports, InterfaceRepairReport(
            (pidA, pidB),
            topo.conformity_ratio,
            length(topo.edges_A),
            length(topo.edges_B),
            length(topo.faces_A),
            length(topo.faces_B);
            success = false,
            error_message = "Repair plan not feasible",
            duration_seconds = time() - interface_start,
            symmetric_mismatches = [serialize_symmetric_mismatch(sm) for sm in sym_mismatches],
            edges_only_A = sym_result.edges_only_in_A,
            edges_only_B = sym_result.edges_only_in_B,
            edges_both = sym_result.edges_in_both,
            agreement_rate = sym_result.agreement_rate
        ))
        
        continue
    end

    # Execute symmetric replacement into the shared workspace (reduced verbosity)
    exec_ok = Nas2Step.execute_symmetric_repair!(ws, sym_plan; verbose=false)
    
    # Collect post-repair statistics
    nm_stats = collect_nonmanifold_stats_from_workspace(ws, pidA, pidB)
    conformity_after_val = compute_conformity_after(ws, pidA, pidB)
    
    # Strategy counts
    use_A_count = count(sm -> sm.repair_strategy == :use_A, sym_mismatches)
    use_B_count = count(sm -> sm.repair_strategy == :use_B, sym_mismatches)
    compromise_count = count(sm -> sm.repair_strategy == :compromise, sym_mismatches)
    skip_count = count(sm -> sm.repair_strategy == :skip, sym_mismatches)
    
    # Build interface report
    push!(interface_reports, InterfaceRepairReport(
        (pidA, pidB),
        topo.conformity_ratio,
        length(topo.edges_A),
        length(topo.edges_B),
        length(topo.faces_A),
        length(topo.faces_B);
        symmetric_mismatches = [serialize_symmetric_mismatch(sm) for sm in sym_mismatches],
        edges_only_A = sym_result.edges_only_in_A,
        edges_only_B = sym_result.edges_only_in_B,
        edges_both = sym_result.edges_in_both,
        agreement_rate = sym_result.agreement_rate,
        edges_use_A = use_A_count,
        edges_use_B = use_B_count,
        edges_compromise = compromise_count,
        edges_skipped = skip_count,
        unified_triangles = length(sym_plan.target_unified_mesh.triangles),
        unified_min_quality = sym_plan.target_unified_mesh.min_triangle_quality,
        unified_total_area = sym_plan.target_unified_mesh.total_area,
        unified_compatible_A = sym_plan.target_unified_mesh.compatible_with_A,
        unified_compatible_B = sym_plan.target_unified_mesh.compatible_with_B,
        nm_total_edges = nm_stats[1],
        nm_boundary_edges = nm_stats[2],
        nm_manifold_edges = nm_stats[3],
        nm_nonmanifold_edges = nm_stats[4],
        nm_interface_count = nm_stats[5],
        nm_within_A_count = nm_stats[6],
        nm_within_B_count = nm_stats[7],
        nm_max_incidence = nm_stats[8],
        nm_incidence_distribution = nm_stats[9],
        bv_loops_A = 0,  # Would need boundary verification integration
        bv_loops_B = 0,
        bv_loops_unified = 0,
        bv_matched_A = 0,
        bv_matched_B = 0,
        bv_unmatched_A = 0,
        bv_unmatched_B = 0,
        bv_normals_consistent = false,
        bv_normals_flipped = 0,
        bv_normals_degenerate = 0,
        bv_normal_alignment_score = 0.0,
        success = exec_ok,
        modifications_count = length(sym_plan.target_unified_mesh.triangles),  # Approximation
        conformity_after = exec_ok ? conformity_after_val : nothing,
        error_message = exec_ok ? nothing : "Execution failed",
        duration_seconds = time() - interface_start,
        compatibility_issues = sym_plan.target_unified_mesh.compatibility_report,
        boundary_issues = String[]
    ))
    
    if exec_ok
        global total_interfaces_repaired
        total_interfaces_repaired += 1
    else
        println("  → Execution failed for this interface, continuing to next.")
    end
end

println("\n" * "="^70)
if total_interfaces_repaired > 0
    # Export surfaces for updated workspace (surface-only export)
    write_workspace_surfaces!(ws, OUTPUT_FILE)
    println("Repaired $total_interfaces_repaired interface(s). (surface-only export)")
    println("Wrote repaired mesh to: $OUTPUT_FILE")
else
    println("No repairs applied. Copying input to output unchanged.")
    cp(INPUT_FILE, OUTPUT_FILE; force=true)
end

# Build and export session report
session_duration = time() - session_start
total_modifications = sum(r.modifications_count for r in interface_reports)
interfaces_repaired_count = count(r -> r.success, interface_reports)
interfaces_failed = length(interface_reports) - interfaces_repaired_count
success_rate = length(interface_reports) > 0 ? interfaces_repaired_count / length(interface_reports) : 0.0

session_report = RepairSessionReport(
    INPUT_FILE,
    OUTPUT_FILE,
    string(Dates.now()),
    length(pairs),
    length(interface_reports),
    interfaces_repaired_count,
    interfaces_failed,
    success_rate,
    total_modifications,
    session_duration,
    interface_reports
)

# Export JSON report
json_output = replace(OUTPUT_FILE, ".nas" => "_repair_report.json")
export_session_report_json(session_report, json_output)
println("="^70)
