# Interface Conformity Check
# Examines the fundamental basis of interface conformity: whether adjacent regions
# share the same set of vertices at their common interface

using Printf

# Note: geometric_utilities.jl is included by the parent module

"""
Conformity level classification for an interface.
"""
@enum ConformityLevel begin
    PERFECTLY_CONFORMING    # All vertices match exactly
    VERTEX_CONFORMING       # Vertices match, but edge triangulation differs
    PARTIALLY_CONFORMING    # Most vertices match, some discrepancies
    NON_CONFORMING          # Significant vertex mismatch
    DISCONNECTED            # No shared vertices at all
end

"""
Result of vertex-level conformity check for an interface.
"""
struct VertexConformityReport
    pidA::Int
    pidB::Int
    
    # Vertex analysis
    vertices_A::Set{NTuple{3,Float64}}       # Vertices on A's boundary
    vertices_B::Set{NTuple{3,Float64}}       # Vertices on B's boundary
    shared_vertices::Set{NTuple{3,Float64}}  # Vertices in both
    
    vertices_only_A::Set{NTuple{3,Float64}}  # In A but not B
    vertices_only_B::Set{NTuple{3,Float64}}  # In B but not A
    
    # Conformity metrics
    conformity_level::ConformityLevel
    vertex_match_ratio::Float64              # shared / total unique
    vertex_coverage_A::Float64               # shared / |A|
    vertex_coverage_B::Float64               # shared / |B|
    
    # Spatial distribution of mismatches
    mismatch_clusters::Vector{Vector{NTuple{3,Float64}}}  # Groups of nearby mismatched vertices
    max_mismatch_distance::Float64           # Furthest distance between any two mismatched vertices
    
    # Diagnostics
    is_acceptable::Bool                       # Whether this interface can be repaired
    issues::Vector{String}                    # Description of problems found
end

# Note: extract_boundary_vertices and find_matching_vertex are now in geometric_utilities.jl

"""
    compute_vertex_correspondence(vertices_A, vertices_B; tol=1e-4)

Compute which vertices are shared and which are unique to each side.
Returns (shared, only_A, only_B).
"""
function compute_vertex_correspondence(vertices_A::Set{NTuple{3,Float64}}, 
                                      vertices_B::Set{NTuple{3,Float64}}; 
                                      tol::Real=1e-4)
    shared = Set{NTuple{3,Float64}}()
    only_A = Set{NTuple{3,Float64}}()
    only_B = Set{NTuple{3,Float64}}(vertices_B)
    
    for v_a in vertices_A
        match = find_matching_vertex(v_a, only_B, tol=tol)
        if match !== nothing
            push!(shared, v_a)  # Use A's coordinate as canonical
            delete!(only_B, match)
        else
            push!(only_A, v_a)
        end
    end
    
    return (shared, only_A, only_B)
end

# ============================================================================
# Spatial clustering of mismatches
# ============================================================================

"""
    cluster_nearby_vertices(vertices; max_distance=10.0)

Group vertices into spatial clusters based on proximity.
"""
function cluster_nearby_vertices(vertices::Set{NTuple{3,Float64}}; 
                                max_distance::Float64=10.0)::Vector{Vector{NTuple{3,Float64}}}
    if isempty(vertices)
        return Vector{NTuple{3,Float64}}[]
    end
    
    vertex_list = collect(vertices)
    n = length(vertex_list)
    visited = falses(n)
    clusters = Vector{Vector{NTuple{3,Float64}}}()
    
    max_dist2 = max_distance * max_distance
    
    for i in 1:n
        if visited[i]
            continue
        end
        
        # Start new cluster
        cluster = NTuple{3,Float64}[]
        queue = [i]
        visited[i] = true
        
        while !isempty(queue)
            idx = popfirst!(queue)
            push!(cluster, vertex_list[idx])
            
            # Find nearby vertices
            for j in 1:n
                if visited[j]
                    continue
                end
                
                dx = vertex_list[idx][1] - vertex_list[j][1]
                dy = vertex_list[idx][2] - vertex_list[j][2]
                dz = vertex_list[idx][3] - vertex_list[j][3]
                dist2 = dx*dx + dy*dy + dz*dz
                
                if dist2 <= max_dist2
                    visited[j] = true
                    push!(queue, j)
                end
            end
        end
        
        push!(clusters, cluster)
    end
    
    # Sort clusters by size (largest first)
    sort!(clusters, by=length, rev=true)
    
    return clusters
end

"""
    compute_max_distance(vertices)

Compute the maximum distance between any two vertices.
"""
function compute_max_distance(vertices::Set{NTuple{3,Float64}})::Float64
    if length(vertices) <= 1
        return 0.0
    end
    
    vertex_list = collect(vertices)
    max_dist = 0.0
    
    for i in 1:length(vertex_list)
        for j in (i+1):length(vertex_list)
            dx = vertex_list[i][1] - vertex_list[j][1]
            dy = vertex_list[i][2] - vertex_list[j][2]
            dz = vertex_list[i][3] - vertex_list[j][3]
            dist = sqrt(dx*dx + dy*dy + dz*dz)
            max_dist = max(max_dist, dist)
        end
    end
    
    return max_dist
end

# ============================================================================
# Conformity classification
# ============================================================================

"""
    classify_conformity_level(vertex_match_ratio, coverage_A, coverage_B)

Classify the conformity level based on vertex matching statistics.
"""
function classify_conformity_level(vertex_match_ratio::Float64, 
                                  coverage_A::Float64, 
                                  coverage_B::Float64)::ConformityLevel
    # Perfect conformity: all vertices match
    if vertex_match_ratio >= 0.9999 && coverage_A >= 0.9999 && coverage_B >= 0.9999
        return PERFECTLY_CONFORMING
    end
    
    # Vertex conforming: vertices match well, but may have different triangulations
    if vertex_match_ratio >= 0.95 && coverage_A >= 0.90 && coverage_B >= 0.90
        return VERTEX_CONFORMING
    end
    
    # Partially conforming: most vertices match
    if vertex_match_ratio >= 0.70 && coverage_A >= 0.60 && coverage_B >= 0.60
        return PARTIALLY_CONFORMING
    end
    
    # Non-conforming: significant mismatch
    if vertex_match_ratio >= 0.10
        return NON_CONFORMING
    end
    
    # Disconnected: almost no shared vertices
    return DISCONNECTED
end

# ============================================================================
# Main conformity check function
# ============================================================================

"""
    check_vertex_conformity(topology::InterfaceTopology; tol=1e-4, cluster_distance=10.0)

Perform a comprehensive vertex-level conformity check on an interface.
"""
function check_vertex_conformity(topology::InterfaceTopology; 
                                tol::Real=1e-4,
                                cluster_distance::Float64=10.0)::VertexConformityReport
    
    # Extract vertices from boundary faces
    vertices_A = extract_boundary_vertices(topology.faces_A)
    vertices_B = extract_boundary_vertices(topology.faces_B)
    
    # Compute correspondence
    shared, only_A, only_B = compute_vertex_correspondence(vertices_A, vertices_B, tol=tol)
    
    # Compute metrics
    total_unique = length(vertices_A) + length(vertices_B) - length(shared)
    vertex_match_ratio = total_unique > 0 ? length(shared) / total_unique : 0.0
    
    coverage_A = length(vertices_A) > 0 ? length(shared) / length(vertices_A) : 0.0
    coverage_B = length(vertices_B) > 0 ? length(shared) / length(vertices_B) : 0.0
    
    # Classify conformity level
    conformity_level = classify_conformity_level(vertex_match_ratio, coverage_A, coverage_B)
    
    # Cluster mismatched vertices
    all_mismatches = union(only_A, only_B)
    mismatch_clusters = cluster_nearby_vertices(all_mismatches, max_distance=cluster_distance)
    max_mismatch_distance = compute_max_distance(all_mismatches)
    
    # Assess acceptability and generate issues
    issues = String[]
    is_acceptable = true
    
    if conformity_level == DISCONNECTED
        push!(issues, "CRITICAL: Interface is disconnected - no shared vertices found")
        is_acceptable = false
    elseif conformity_level == NON_CONFORMING
        push!(issues, "SEVERE: Interface is non-conforming - only $(round(vertex_match_ratio*100, digits=1))% vertices match")
        is_acceptable = false
    elseif conformity_level == PARTIALLY_CONFORMING
        push!(issues, "WARNING: Interface is partially conforming - $(round(vertex_match_ratio*100, digits=1))% vertices match")
        if coverage_A < 0.70
            push!(issues, "Side A coverage is low: $(round(coverage_A*100, digits=1))%")
        end
        if coverage_B < 0.70
            push!(issues, "Side B coverage is low: $(round(coverage_B*100, digits=1))%")
        end
    end
    
    if length(only_A) > 0
        push!(issues, "$(length(only_A)) vertices unique to side A (PID $(topology.pidA))")
    end
    
    if length(only_B) > 0
        push!(issues, "$(length(only_B)) vertices unique to side B (PID $(topology.pidB))")
    end
    
    if length(mismatch_clusters) > 0
        push!(issues, "Mismatches form $(length(mismatch_clusters)) spatial cluster(s)")
        if length(mismatch_clusters) == 1
            push!(issues, "Largest cluster: $(length(mismatch_clusters[1])) vertices")
        else
            push!(issues, "Largest cluster: $(length(mismatch_clusters[1])) vertices, smallest: $(length(mismatch_clusters[end])) vertices")
        end
    end
    
    return VertexConformityReport(
        topology.pidA,
        topology.pidB,
        vertices_A,
        vertices_B,
        shared,
        only_A,
        only_B,
        conformity_level,
        vertex_match_ratio,
        coverage_A,
        coverage_B,
        mismatch_clusters,
        max_mismatch_distance,
        is_acceptable,
        issues
    )
end

# ============================================================================
# Reporting functions
# ============================================================================

"""
    print_conformity_report(report::VertexConformityReport; verbose=false)

Print a human-readable conformity report.
"""
function print_conformity_report(report::VertexConformityReport; verbose::Bool=false)
    println("\n" * "="^70)
    println("Vertex Conformity Report: PID $(report.pidA) ↔ PID $(report.pidB)")
    println("="^70)
    
    # Overall status
    status_symbol = report.is_acceptable ? "✓" : "✗"
    println("\n$status_symbol Overall Status: $(report.conformity_level)")
    
    # Vertex counts
    println("\nVertex Counts:")
    println("  Side A (PID $(report.pidA)): $(length(report.vertices_A)) vertices")
    println("  Side B (PID $(report.pidB)): $(length(report.vertices_B)) vertices")
    println("  Shared: $(length(report.shared_vertices)) vertices")
    println("  Unique to A: $(length(report.vertices_only_A)) vertices")
    println("  Unique to B: $(length(report.vertices_only_B)) vertices")
    
    # Metrics
    println("\nConformity Metrics:")
    println("  Vertex match ratio: $(round(report.vertex_match_ratio * 100, digits=2))%")
    println("  Coverage A: $(round(report.vertex_coverage_A * 100, digits=2))%")
    println("  Coverage B: $(round(report.vertex_coverage_B * 100, digits=2))%")
    
    # Mismatch distribution
    if !isempty(report.vertices_only_A) || !isempty(report.vertices_only_B)
        println("\nMismatch Distribution:")
        println("  Number of clusters: $(length(report.mismatch_clusters))")
        println("  Max mismatch span: $(round(report.max_mismatch_distance, digits=2))")
        
        if verbose && !isempty(report.mismatch_clusters)
            println("\n  Cluster details:")
            for (i, cluster) in enumerate(report.mismatch_clusters[1:min(5, end)])
                println("    Cluster $i: $(length(cluster)) vertices")
            end
            if length(report.mismatch_clusters) > 5
                println("    ... and $(length(report.mismatch_clusters) - 5) more")
            end
        end
    end
    
    # Issues
    if !isempty(report.issues)
        println("\nIssues Detected:")
        for issue in report.issues
            println("  • $issue")
        end
    end
    
    # Recommendations
    println("\nRecommendations:")
    if report.conformity_level == DISCONNECTED
        println("  ⚠ Interface is fundamentally broken - requires complete remeshing")
    elseif report.conformity_level == NON_CONFORMING
        println("  ⚠ Interface has severe vertex mismatch - surgical repair not recommended")
        println("  → Consider constrainted remeshing of the interface region")
    elseif report.conformity_level == PARTIALLY_CONFORMING
        println("  ⚠ Interface has significant vertex mismatch")
        println("  → Review mismatch clusters to understand the source")
        println("  → Edge-level repairs may have limited success")
    elseif report.conformity_level == VERTEX_CONFORMING
        println("  ✓ Vertices are well-aligned")
        println("  → Edge mismatches are likely due to different triangulations")
        println("  → Surgical repair (edge insertion, quad flips) should be effective")
    else  # PERFECTLY_CONFORMING
        println("  ✓ Perfect vertex conformity")
        println("  → Any edge mismatches can be resolved with simple retriangulation")
    end
    
    println("="^70)
end

"""
    export_conformity_report_json(report::VertexConformityReport, output_file::String)

Export conformity report to JSON format.
"""
function export_conformity_report_json(report::VertexConformityReport, output_file::String)
    report_dict = Dict(
        "interface" => Dict(
            "pidA" => report.pidA,
            "pidB" => report.pidB
        ),
        "conformity_level" => string(report.conformity_level),
        "is_acceptable" => report.is_acceptable,
        "vertex_counts" => Dict(
            "side_A" => length(report.vertices_A),
            "side_B" => length(report.vertices_B),
            "shared" => length(report.shared_vertices),
            "unique_to_A" => length(report.vertices_only_A),
            "unique_to_B" => length(report.vertices_only_B)
        ),
        "metrics" => Dict(
            "vertex_match_ratio" => round(report.vertex_match_ratio, digits=4),
            "coverage_A" => round(report.vertex_coverage_A, digits=4),
            "coverage_B" => round(report.vertex_coverage_B, digits=4),
            "max_mismatch_distance" => round(report.max_mismatch_distance, digits=4)
        ),
        "mismatch_clusters" => Dict(
            "count" => length(report.mismatch_clusters),
            "sizes" => [length(c) for c in report.mismatch_clusters]
        ),
        "issues" => report.issues
    )
    
    open(output_file, "w") do io
        write_json(io, report_dict, 0)
    end
    
    println("Conformity report exported to: $output_file")
end
