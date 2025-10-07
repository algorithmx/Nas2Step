# json_reporting_types.jl
# Data structures for comprehensive JSON repair reporting
# Part of the symmetric repair architecture

using Statistics

"""
    InterfaceRepairReport

Complete report for a single interface repair operation.
Contains pre/post state, classification results, strategies, quality metrics, and diagnostics.
"""
struct InterfaceRepairReport
    # Interface identification
    interface_pair::Tuple{Int,Int}
    
    # Pre-repair topology state
    conformity_before::Float64
    edges_A_before::Int
    edges_B_before::Int
    faces_A_before::Int
    faces_B_before::Int
    
    # Symmetric classification results
    symmetric_mismatches::Vector{Dict{String,Any}}  # Serialized SymmetricEdgeMismatch
    edges_only_A::Int
    edges_only_B::Int
    edges_both::Int
    agreement_rate::Float64
    
    # Repair strategy decisions
    edges_use_A::Int
    edges_use_B::Int
    edges_compromise::Int
    edges_skipped::Int
    
    # Unified mesh results
    unified_triangles::Int
    unified_min_quality::Float64
    unified_total_area::Float64
    unified_compatible_A::Bool
    unified_compatible_B::Bool
    
    # Non-manifold analysis
    nm_total_edges::Int
    nm_boundary_edges::Int
    nm_manifold_edges::Int
    nm_nonmanifold_edges::Int
    nm_interface_count::Int
    nm_within_A_count::Int
    nm_within_B_count::Int
    nm_max_incidence::Int
    nm_incidence_distribution::Dict{Int,Int}
    
    # Boundary verification
    bv_loops_A::Int
    bv_loops_B::Int
    bv_loops_unified::Int
    bv_matched_A::Int
    bv_matched_B::Int
    bv_unmatched_A::Int
    bv_unmatched_B::Int
    bv_normals_consistent::Bool
    bv_normals_flipped::Int
    bv_normals_degenerate::Int
    bv_normal_alignment_score::Float64
    
    # Execution results
    success::Bool
    modifications_count::Int
    conformity_after::Union{Float64,Nothing}
    error_message::Union{String,Nothing}
    duration_seconds::Float64
    
    # Detailed diagnostics (optional, can be large)
    compatibility_issues::Vector{String}
    boundary_issues::Vector{String}
end

"""
    RepairSessionReport

Top-level report for an entire repair session covering all interfaces.
Contains summary statistics and per-interface detailed reports.
"""
struct RepairSessionReport
    # Input/output files
    input_file::String
    output_file::String
    timestamp::String
    
    # Summary statistics
    interfaces_detected::Int
    interfaces_attempted::Int
    interfaces_repaired::Int
    interfaces_failed::Int
    success_rate::Float64
    
    # Aggregated metrics
    total_modifications::Int
    total_duration_seconds::Float64
    
    # Per-interface detailed reports
    interface_reports::Vector{InterfaceRepairReport}
end

# Constructor helpers for safer initialization

"""
    InterfaceRepairReport constructor with safe defaults
"""
function InterfaceRepairReport(
    interface_pair::Tuple{Int,Int},
    conformity_before::Float64,
    edges_A_before::Int,
    edges_B_before::Int,
    faces_A_before::Int,
    faces_B_before::Int;
    symmetric_mismatches::Vector{Dict{String,Any}} = Dict{String,Any}[],
    edges_only_A::Int = 0,
    edges_only_B::Int = 0,
    edges_both::Int = 0,
    agreement_rate::Float64 = 0.0,
    edges_use_A::Int = 0,
    edges_use_B::Int = 0,
    edges_compromise::Int = 0,
    edges_skipped::Int = 0,
    unified_triangles::Int = 0,
    unified_min_quality::Float64 = 0.0,
    unified_total_area::Float64 = 0.0,
    unified_compatible_A::Bool = false,
    unified_compatible_B::Bool = false,
    nm_total_edges::Int = 0,
    nm_boundary_edges::Int = 0,
    nm_manifold_edges::Int = 0,
    nm_nonmanifold_edges::Int = 0,
    nm_interface_count::Int = 0,
    nm_within_A_count::Int = 0,
    nm_within_B_count::Int = 0,
    nm_max_incidence::Int = 0,
    nm_incidence_distribution::Dict{Int,Int} = Dict{Int,Int}(),
    bv_loops_A::Int = 0,
    bv_loops_B::Int = 0,
    bv_loops_unified::Int = 0,
    bv_matched_A::Int = 0,
    bv_matched_B::Int = 0,
    bv_unmatched_A::Int = 0,
    bv_unmatched_B::Int = 0,
    bv_normals_consistent::Bool = false,
    bv_normals_flipped::Int = 0,
    bv_normals_degenerate::Int = 0,
    bv_normal_alignment_score::Float64 = 0.0,
    success::Bool = false,
    modifications_count::Int = 0,
    conformity_after::Union{Float64,Nothing} = nothing,
    error_message::Union{String,Nothing} = nothing,
    duration_seconds::Float64 = 0.0,
    compatibility_issues::Vector{String} = String[],
    boundary_issues::Vector{String} = String[]
)
    return InterfaceRepairReport(
        interface_pair,
        conformity_before,
        edges_A_before,
        edges_B_before,
        faces_A_before,
        faces_B_before,
        symmetric_mismatches,
        edges_only_A,
        edges_only_B,
        edges_both,
        agreement_rate,
        edges_use_A,
        edges_use_B,
        edges_compromise,
        edges_skipped,
        unified_triangles,
        unified_min_quality,
        unified_total_area,
        unified_compatible_A,
        unified_compatible_B,
        nm_total_edges,
        nm_boundary_edges,
        nm_manifold_edges,
        nm_nonmanifold_edges,
        nm_interface_count,
        nm_within_A_count,
        nm_within_B_count,
        nm_max_incidence,
        nm_incidence_distribution,
        bv_loops_A,
        bv_loops_B,
        bv_loops_unified,
        bv_matched_A,
        bv_matched_B,
        bv_unmatched_A,
        bv_unmatched_B,
        bv_normals_consistent,
        bv_normals_flipped,
        bv_normals_degenerate,
        bv_normal_alignment_score,
        success,
        modifications_count,
        conformity_after,
        error_message,
        duration_seconds,
        compatibility_issues,
        boundary_issues
    )
end
