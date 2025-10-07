# json_export.jl
# JSON serialization and export functions for symmetric repair reporting
# Part of the symmetric repair architecture

using JSON
using Statistics
using Dates

include("json_reporting_types.jl")
include("symmetric_classification.jl")
include("unified_mesh_generation.jl")
include("boundary_verification.jl")
include("repair_workspace.jl")

"""
    serialize_symmetric_mismatch(sym::SymmetricEdgeMismatch)::Dict{String,Any}

Serialize a symmetric edge mismatch to JSON-compatible dictionary.
"""
function serialize_symmetric_mismatch(sym::SymmetricEdgeMismatch)::Dict{String,Any}
    d = Dict{String,Any}(
        "edge" => Dict(
            "node1" => [sym.edge_key.node1...],
            "node2" => [sym.edge_key.node2...]
        ),
        "presence" => Dict(
            "in_A" => sym.present_in_A,
            "in_B" => sym.present_in_B,
            "in_both" => sym.present_in_both
        ),
        "agreement" => Dict(
            "on_type" => sym.agree_on_type,
            "on_feasibility" => sym.agree_on_feasibility
        ),
        "strategy" => string(sym.repair_strategy),
        "priority" => round(sym.repair_priority, digits=3),
        "reason" => sym.resolution_reason
    )
    
    # Add A perspective if present
    if sym.classification_A_perspective !== nothing
        m = sym.classification_A_perspective
        d["perspective_A"] = Dict(
            "type" => string(m.mismatch_type),
            "feasible" => m.repair_feasible,
            "complexity" => round(m.complexity_score, digits=3),
            "min_quality" => round(m.min_affected_triangle_quality, digits=3),
            "hanging_nodes_count" => length(m.hanging_nodes),
            "triangles_to_replace_count" => length(m.triangles_to_replace),
            "affected_triangles_count" => length(m.affected_triangles)
        )
    end
    
    # Add B perspective if present  
    if sym.classification_B_perspective !== nothing
        m = sym.classification_B_perspective
        d["perspective_B"] = Dict(
            "type" => string(m.mismatch_type),
            "feasible" => m.repair_feasible,
            "complexity" => round(m.complexity_score, digits=3),
            "min_quality" => round(m.min_affected_triangle_quality, digits=3),
            "hanging_nodes_count" => length(m.hanging_nodes),
            "triangles_to_replace_count" => length(m.triangles_to_replace),
            "affected_triangles_count" => length(m.affected_triangles)
        )
    end
    
    return d
end

"""
    serialize_unified_mesh(mesh::UnifiedInterfaceMesh)::Dict{String,Any}

Serialize a unified interface mesh to JSON-compatible dictionary.
"""
function serialize_unified_mesh(mesh::UnifiedInterfaceMesh)::Dict{String,Any}
    # Compute mean quality
    qualities = [compute_triangle_quality(t) for t in mesh.triangles]
    mean_quality = isempty(qualities) ? 0.0 : mean(qualities)
    
    # Count provenance
    from_A_count = count(==(from_A), mesh.triangle_provenance)
    from_B_count = count(==(from_B), mesh.triangle_provenance)
    synthesized_count = count(==(synthesized), mesh.triangle_provenance)
    from_compromise_count = count(==(from_compromise), mesh.triangle_provenance)
    
    return Dict{String,Any}(
        "triangles" => length(mesh.triangles),
        "quality" => Dict(
            "min" => round(mesh.min_triangle_quality, digits=4),
            "mean" => round(mean_quality, digits=4)
        ),
        "area" => Dict(
            "total" => round(mesh.total_area, digits=4)
        ),
        "compatibility" => Dict(
            "with_A" => mesh.compatible_with_A,
            "with_B" => mesh.compatible_with_B,
            "issues" => mesh.compatibility_report
        ),
        "provenance" => Dict(
            "from_A" => from_A_count,
            "from_B" => from_B_count,
            "synthesized" => synthesized_count,
            "from_compromise" => from_compromise_count
        ),
        "topology" => Dict(
            "edges" => length(mesh.edges),
            "nodes_mapped_A" => count(v -> v !== nothing, values(mesh.node_mapping_A)),
            "nodes_mapped_B" => count(v -> v !== nothing, values(mesh.node_mapping_B))
        )
    )
end

"""
    serialize_nonmanifold_stats(...)::Dict{String,Any}

Serialize non-manifold statistics to JSON-compatible dictionary.
"""
function serialize_nonmanifold_stats(
    total::Int,
    boundary::Int,
    manifold::Int,
    nonmanifold::Int,
    nm_interface::Int,
    nm_within_A::Int,
    nm_within_B::Int,
    max_inc::Int,
    incidence_hist::Dict{Int,Int}
)::Dict{String,Any}
    # Convert incidence histogram keys to strings for JSON compatibility
    incidence_hist_str = Dict(string(k) => v for (k, v) in incidence_hist)
    
    # Compute health assessment
    health = if nonmanifold == 0
        "excellent"
    elseif nm_interface > nonmanifold * 0.8
        "good_interface_expected"
    elseif nm_within_A > 0 || nm_within_B > 0
        "concerning_internal_issues"
    else
        "moderate"
    end
    
    return Dict{String,Any}(
        "total_edges" => total,
        "boundary_edges" => boundary,
        "manifold_edges" => manifold,
        "nonmanifold_edges" => nonmanifold,
        "nonmanifold_percentage" => round(100 * nonmanifold / max(1, total), digits=2),
        "breakdown" => Dict(
            "interface" => nm_interface,
            "within_A" => nm_within_A,
            "within_B" => nm_within_B
        ),
        "severity" => Dict(
            "max_incidence" => max_inc,
            "distribution" => incidence_hist_str
        ),
        "health" => health
    )
end

"""
    serialize_boundary_verification(report::BoundaryVerificationReport)::Dict{String,Any}

Serialize boundary verification report to JSON-compatible dictionary.
"""
function serialize_boundary_verification(report::BoundaryVerificationReport)::Dict{String,Any}
    return Dict{String,Any}(
        "loops" => Dict(
            "A" => report.loops_A,
            "B" => report.loops_B,
            "unified" => report.loops_unified,
            "matched_A" => report.matched_A,
            "matched_B" => report.matched_B,
            "unmatched_A" => report.unmatched_A,
            "unmatched_B" => report.unmatched_B
        ),
        "normals" => Dict(
            "consistent" => report.normals_consistent,
            "flipped" => report.normals_flipped,
            "degenerate" => report.normals_degenerate,
            "alignment_score" => round(report.normal_alignment_score, digits=3)
        ),
        "boundary_equivalent" => Dict(
            "A" => report.boundary_equivalent_A,
            "B" => report.boundary_equivalent_B
        ),
        "issues" => report.issues
    )
end

"""
    serialize_interface_report(report::InterfaceRepairReport)::Dict{String,Any}

Serialize a complete interface repair report to JSON-compatible dictionary.
"""
function serialize_interface_report(report::InterfaceRepairReport)::Dict{String,Any}
    # Compute conformity improvement
    conformity_improvement = if report.conformity_after !== nothing
        report.conformity_after - report.conformity_before
    else
        nothing
    end
    
    return Dict{String,Any}(
        "pair" => [report.interface_pair[1], report.interface_pair[2]],
        "conformity" => Dict(
            "before" => round(report.conformity_before, digits=4),
            "after" => report.conformity_after !== nothing ? round(report.conformity_after, digits=4) : nothing,
            "improvement" => conformity_improvement !== nothing ? round(conformity_improvement, digits=4) : nothing
        ),
        "topology" => Dict(
            "edges_A" => report.edges_A_before,
            "edges_B" => report.edges_B_before,
            "faces_A" => report.faces_A_before,
            "faces_B" => report.faces_B_before
        ),
        "classification" => Dict(
            "unique_edges" => report.edges_only_A + report.edges_only_B + report.edges_both,
            "only_A" => report.edges_only_A,
            "only_B" => report.edges_only_B,
            "both" => report.edges_both,
            "agreement_rate" => round(report.agreement_rate, digits=3),
            "symmetric_mismatches" => report.symmetric_mismatches
        ),
        "strategies" => Dict(
            "use_A" => report.edges_use_A,
            "use_B" => report.edges_use_B,
            "compromise" => report.edges_compromise,
            "skipped" => report.edges_skipped
        ),
        "unified_mesh" => Dict(
            "triangles" => report.unified_triangles,
            "quality" => Dict(
                "min" => round(report.unified_min_quality, digits=4),
            ),
            "area" => Dict(
                "total" => round(report.unified_total_area, digits=4)
            ),
            "compatibility" => Dict(
                "with_A" => report.unified_compatible_A,
                "with_B" => report.unified_compatible_B,
                "issues" => report.compatibility_issues
            )
        ),
        "non_manifold" => serialize_nonmanifold_stats(
            report.nm_total_edges,
            report.nm_boundary_edges,
            report.nm_manifold_edges,
            report.nm_nonmanifold_edges,
            report.nm_interface_count,
            report.nm_within_A_count,
            report.nm_within_B_count,
            report.nm_max_incidence,
            report.nm_incidence_distribution
        ),
        "boundary_verification" => Dict(
            "loops" => Dict(
                "A" => report.bv_loops_A,
                "B" => report.bv_loops_B,
                "unified" => report.bv_loops_unified,
                "matched_A" => report.bv_matched_A,
                "matched_B" => report.bv_matched_B,
                "unmatched_A" => report.bv_unmatched_A,
                "unmatched_B" => report.bv_unmatched_B
            ),
            "normals" => Dict(
                "consistent" => report.bv_normals_consistent,
                "flipped" => report.bv_normals_flipped,
                "degenerate" => report.bv_normals_degenerate,
                "alignment_score" => round(report.bv_normal_alignment_score, digits=3)
            ),
            "issues" => report.boundary_issues
        ),
        "execution" => Dict(
            "success" => report.success,
            "modifications" => report.modifications_count,
            "duration_seconds" => round(report.duration_seconds, digits=2),
            "error" => report.error_message
        )
    )
end

"""
    serialize_session_report(report::RepairSessionReport)::Dict{String,Any}

Serialize a complete repair session report to JSON-compatible dictionary.
"""
function serialize_session_report(report::RepairSessionReport)::Dict{String,Any}
    return Dict{String,Any}(
        "session" => Dict(
            "input_file" => report.input_file,
            "output_file" => report.output_file,
            "timestamp" => report.timestamp,
            "duration_seconds" => round(report.total_duration_seconds, digits=2)
        ),
        "summary" => Dict(
            "interfaces_detected" => report.interfaces_detected,
            "interfaces_attempted" => report.interfaces_attempted,
            "interfaces_repaired" => report.interfaces_repaired,
            "interfaces_failed" => report.interfaces_failed,
            "success_rate" => round(report.success_rate, digits=3),
            "total_modifications" => report.total_modifications
        ),
        "interfaces" => [serialize_interface_report(iface) for iface in report.interface_reports]
    )
end

"""
    export_interface_report_json(report::InterfaceRepairReport, output_file::String)

Export a single interface repair report to JSON file.
"""
function export_interface_report_json(report::InterfaceRepairReport, output_file::String)
    data = serialize_interface_report(report)
    open(output_file, "w") do io
        JSON.print(io, data, 2)
    end
    println("✓ Interface report exported: $output_file")
end

"""
    export_session_report_json(report::RepairSessionReport, output_file::String)

Export a complete repair session report to JSON file.
"""
function export_session_report_json(report::RepairSessionReport, output_file::String)
    data = serialize_session_report(report)
    open(output_file, "w") do io
        JSON.print(io, data, 2)
    end
    println("✓ Session report exported: $output_file")
    println("  Interfaces: $(report.interfaces_attempted) attempted, $(report.interfaces_repaired) repaired")
    println("  Success rate: $(round(report.success_rate * 100, digits=1))%")
    println("  Total modifications: $(report.total_modifications)")
end

# ============================================================================
# Helper Functions for Data Collection
# ============================================================================

"""
    collect_nonmanifold_stats_from_workspace(ws::RepairWorkspace, pidA::Int, pidB::Int) -> Tuple

Extract non-manifold statistics from workspace for JSON reporting.
Returns tuple of:
(total, boundary, manifold, nonmanifold, nm_interface, nm_within_A, nm_within_B, max_inc, incidence_hist)
"""
function collect_nonmanifold_stats_from_workspace(
    ws,  # RepairWorkspace (duck-typed to avoid module scoping issues)
    pidA::Int,
    pidB::Int
)
    
    edge_pid_incidence = Dict{Tuple{Int,Int}, Dict{Int,Int}}()
    
    for pid in (pidA, pidB)
        if !haskey(ws.working_faces, pid)
            # Return empty stats if PID not found
            return (0, 0, 0, 0, 0, 0, 0, 0, Dict{Int,Int}())
        end
        
        faces = ws.working_faces[pid]
        for face_nodes in faces
            # Create edges (sorted pairs)
            edges = [
                (min(face_nodes[1], face_nodes[2]), max(face_nodes[1], face_nodes[2])),
                (min(face_nodes[2], face_nodes[3]), max(face_nodes[2], face_nodes[3])),
                (min(face_nodes[3], face_nodes[1]), max(face_nodes[3], face_nodes[1]))
            ]
            
            for edge in edges
                d = get!(edge_pid_incidence, edge, Dict{Int,Int}())
                d[pid] = get(d, pid, 0) + 1
            end
        end
    end
    
    # Analyze edge topology
    total_edges = length(edge_pid_incidence)
    boundary_edges = 0
    manifold_edges = 0
    nonmanifold_edges = 0
    
    nm_within_A = 0
    nm_within_B = 0
    nm_interface = 0
    max_incidence = 0
    
    incidence_histogram = Dict{Int,Int}()
    
    for (edge, perpid) in edge_pid_incidence
        a_cnt = get(perpid, pidA, 0)
        b_cnt = get(perpid, pidB, 0)
        total = a_cnt + b_cnt
        
        if total == 1
            boundary_edges += 1
        elseif total == 2
            manifold_edges += 1
        else  # total > 2: non-manifold
            nonmanifold_edges += 1
            incidence_histogram[total] = get(incidence_histogram, total, 0) + 1
            
            # Categorize by location
            if a_cnt > 0 && b_cnt > 0
                nm_interface += 1
            elseif a_cnt > 0
                nm_within_A += 1
            else
                nm_within_B += 1
            end
            
            max_incidence = max(max_incidence, total)
        end
    end
    
    return (
        total_edges,
        boundary_edges,
        manifold_edges,
        nonmanifold_edges,
        nm_interface,
        nm_within_A,
        nm_within_B,
        max_incidence,
        incidence_histogram
    )
end

"""
    compute_conformity_after(workspace::RepairWorkspace, pidA::Int, pidB::Int) -> Float64

Compute interface conformity after repair by rebuilding topology.
Returns conformity ratio (shared_edges / total_interface_edges).
"""
function compute_conformity_after(
    workspace,  # RepairWorkspace (duck-typed to avoid module scoping issues)
    pidA::Int,
    pidB::Int
)
    
    # Get faces from both PIDs
    if !haskey(workspace.working_faces, pidA) || !haskey(workspace.working_faces, pidB)
        return 0.0
    end
    
    faces_a = workspace.working_faces[pidA]
    faces_b = workspace.working_faces[pidB]
    
    # Build edge sets using coordinates
    edges_a = Set{NTuple{2,NTuple{3,Float64}}}()
    edges_b = Set{NTuple{2,NTuple{3,Float64}}}()
    
    for face_nodes in faces_a
        if length(face_nodes) == 3
            coords = [workspace.working_nodes[nid] for nid in face_nodes]
            for i in 1:3
                j = mod1(i+1, 3)
                c1, c2 = coords[i], coords[j]
                edge = c1 < c2 ? (c1, c2) : (c2, c1)
                push!(edges_a, edge)
            end
        end
    end
    
    for face_nodes in faces_b
        if length(face_nodes) == 3
            coords = [workspace.working_nodes[nid] for nid in face_nodes]
            for i in 1:3
                j = mod1(i+1, 3)
                c1, c2 = coords[i], coords[j]
                edge = c1 < c2 ? (c1, c2) : (c2, c1)
                push!(edges_b, edge)
            end
        end
    end
    
    # Compute conformity
    shared_edges = length(edges_a ∩ edges_b)
    total_edges_a = length(edges_a)
    total_edges_b = length(edges_b)
    
    if total_edges_a == 0 && total_edges_b == 0
        return 0.0
    end
    
    # Conformity ratio: shared / total unique edges
    total_unique = length(edges_a ∪ edges_b)
    return shared_edges / max(1, total_unique)
end
