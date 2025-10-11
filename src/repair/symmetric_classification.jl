# Symmetric Classification Implementation
# Phase 3: Bidirectional edge classification that preserves both perspectives
#
# This module implements symmetric classification where each edge is classified
# from BOTH A→B and B→A perspectives, and both results are preserved for
# local repair strategy selection.

using Statistics

# ============================================================================
# Core symmetric classification function
# ============================================================================

"""
    classify_interface_mismatches_symmetric(topology; tol=1e-4, verbose=true)

Classify all edge mismatches from BOTH A→B and B→A perspectives.

Unlike the original `classify_interface_mismatches_bidirectional()` which merges
results, this function PRESERVES both classifications for each edge, enabling
local per-edge repair strategy decisions.

# Arguments
- `topology::InterfaceTopology`: Interface topology between two PIDs
- `tol::Real=1e-4`: Geometric tolerance for comparisons
- `verbose::Bool=true`: Print progress information

# Returns
`SymmetricClassificationResult` containing:
- `symmetric_mismatches`: Vector of SymmetricEdgeMismatch with both perspectives
- `classification_AB`: A→B perspective classification (reference)
- `classification_BA`: B→A perspective classification (reference)
- Agreement statistics and edge distribution

# Algorithm
1. Classify from A→B perspective (A as source, B as target)
2. Create swapped topology for B→A perspective
3. Classify from B→A perspective (B as source, A as target)
4. Build symmetric mismatches preserving BOTH classifications
5. Compute agreement and distribution statistics

# Example
```julia
topology = build_interface_topology(mesh_file, 4, 5)
result = classify_interface_mismatches_symmetric(topology)

# Access symmetric mismatches
for sym in result.symmetric_mismatches
    if sym.present_in_both
        println("Edge in both, strategies differ") if !sym.agree_on_type
    end
end

# Check agreement
println("Agreement rate: \$(result.agreement_rate * 100)%")
```

# Note
This is the foundation for symmetric repair. The actual repair strategy selection
(Phase 4) will use these symmetric mismatches to make local decisions.
"""
function classify_interface_mismatches_symmetric(
    topology::InterfaceTopology;
    tol::Real=1e-4,
    verbose::Bool=true
)::SymmetricClassificationResult

    # if verbose
    #     println("\n" * "="^70)
    #     println("SYMMETRIC EDGE CLASSIFICATION")
    #     println("="^70)
    #     println("Interface: PID $(topology.pidA) ↔ PID $(topology.pidB)")
    #     println("Strategy: Classify from BOTH perspectives and preserve results")
    #     println()
    # end

    # ========================================================================
    # Step 1: Classify from A→B perspective
    # ========================================================================
    # Suppress individual classification printing - we'll handle it in the merged output
    classification_AB = classify_interface_mismatches(topology, tol=tol, debug=false, verbose=false)

    # ========================================================================
    # Step 2: Create swapped topology for B→A perspective
    # ========================================================================
    topology_swapped = InterfaceTopology(
        topology.pidB,                 # Swap: B becomes pidA
        topology.pidA,                 # Swap: A becomes pidB
        topology.shared_node_keys,
        topology.node_key_to_ids,
        topology.faces_B,              # Swap faces
        topology.faces_A,
        topology.edges_B,              # Swap edges
        topology.edges_A,
        topology.edges_only_in_B,      # Swap edge sets
        topology.edges_only_in_A,
        topology.edges_shared,
        topology.interface_bbox,
        topology.total_shared_nodes,
        topology.total_faces_B,        # Swap face counts
        topology.total_faces_A,
        topology.total_edges_B,        # Swap edge counts
        topology.total_edges_A,
        topology.conformity_ratio,
        0.0, 0.0, 0, 1.0            # consistency metrics: max_vertex_dist, mean_vertex_dist, edge_mismatch_count, triangulation_similarity
    )

    # ========================================================================
    # Step 3: Classify from B→A perspective
    # ========================================================================
    # Suppress individual classification printing - we'll handle it in the merged output
    classification_BA = classify_interface_mismatches(topology_swapped, tol=tol, debug=false, verbose=false)

    # ========================================================================
    # Step 4: Build symmetric mismatch map
    # ========================================================================
    if verbose
        total_AB = length(classification_AB.mismatches_A) + length(classification_AB.mismatches_B)
        total_BA = length(classification_BA.mismatches_A) + length(classification_BA.mismatches_B)
        println("Analyzing interface mismatches...")
        println("  Found $(max(total_AB, total_BA)) total mismatches")
    end

    # Map edge keys to symmetric mismatches
    symmetric_map = Dict{EdgeKey,SymmetricEdgeMismatch}()

    # Process edges present in A (missing in B from A's perspective)
    for mismatch_AB in classification_AB.mismatches_B
        edge_key = mismatch_AB.edge_key

        # Find corresponding B→A classification
        # Since this edge is in A, from B's perspective it's "missing in A" 
        # which would be in classification_BA.mismatches_A (but swapped back)
        mismatch_BA = find_corresponding_mismatch_in_swapped(
            classification_BA.mismatches_A,
            edge_key
        )

        # For now, create symmetric mismatch without strategy (Phase 4)
        symmetric_map[edge_key] = SymmetricEdgeMismatch(
            edge_key,
            mismatch_AB,   # A's perspective
            mismatch_BA,   # B's perspective (may be nothing)
            :pending,      # Strategy selection in Phase 4
            0.5,           # Default priority
            "Strategy not yet determined"
        )
    end

    # Process edges present in B (missing in A from B's perspective)
    for mismatch_BA in classification_BA.mismatches_B
        edge_key = mismatch_BA.edge_key

        # Skip if already processed
        if haskey(symmetric_map, edge_key)
            continue
        end

        # Find corresponding A→B classification
        mismatch_AB = find_corresponding_mismatch_in_swapped(
            classification_AB.mismatches_A,
            edge_key
        )

        symmetric_map[edge_key] = SymmetricEdgeMismatch(
            edge_key,
            mismatch_AB,   # A's perspective (may be nothing)
            mismatch_BA,   # B's perspective
            :pending,      # Strategy selection in Phase 4
            0.5,           # Default priority
            "Strategy not yet determined"
        )
    end

    symmetric_mismatches = collect(values(symmetric_map))

    if verbose
        println("  Created $(length(symmetric_mismatches)) symmetric mismatches")
    end

    # ========================================================================
    # Step 5: Compute statistics
    # ========================================================================
    edges_only_in_A = count(m -> m.present_in_A && !m.present_in_B, symmetric_mismatches)
    edges_only_in_B = count(m -> m.present_in_B && !m.present_in_A, symmetric_mismatches)
    edges_in_both = count(m -> m.present_in_both, symmetric_mismatches)

    # Agreement statistics (only for edges present in both)
    edges_with_both_classifications = filter(m -> m.present_in_both, symmetric_mismatches)

    if !isempty(edges_with_both_classifications)
        agree_on_type_count = count(m -> m.agree_on_type, edges_with_both_classifications)
        agree_on_feasibility_count = count(m -> m.agree_on_feasibility, edges_with_both_classifications)
        agreement_rate = (agree_on_type_count + agree_on_feasibility_count) / (2.0 * length(edges_with_both_classifications))
    else
        agreement_rate = 1.0  # No conflicts if no edges in both
    end

    # Build comparison metrics
    comparison_metrics = Dict{String,Any}(
        "perspective_AB" => Dict(
            "edges_only_in_A" => length(classification_AB.mismatches_B),
            "edges_only_in_B" => length(classification_AB.mismatches_A),
            "total" => length(classification_AB.mismatches_A) + length(classification_AB.mismatches_B)
        ),
        "perspective_BA" => Dict(
            "edges_only_in_B" => length(classification_BA.mismatches_B),
            "edges_only_in_A" => length(classification_BA.mismatches_A),
            "total" => length(classification_BA.mismatches_A) + length(classification_BA.mismatches_B)
        ),
        "symmetric" => Dict(
            "total_unique_edges" => length(symmetric_mismatches),
            "edges_only_in_A" => edges_only_in_A,
            "edges_only_in_B" => edges_only_in_B,
            "edges_in_both" => edges_in_both
        ),
        "agreement" => Dict(
            "edges_with_both_classifications" => length(edges_with_both_classifications),
            "agreement_rate" => agreement_rate
        )
    )

    # Compute interface-local non-manifold edge stats (pre-repair)
    interface_manifold_stats = compute_interface_nonmanifold_stats(topology)

    if verbose
        # If no mismatches found (already conforming), show minimal report
        if isempty(symmetric_mismatches)
            # Only show non-manifold info if there are any
            if !isempty(interface_manifold_stats) && interface_manifold_stats["non_manifold_count"] > 0
                nm_total = interface_manifold_stats["non_manifold_count"]
                nm_pct = round(interface_manifold_stats["percent_non_manifold"] * 100, digits=2)
                max_inc = interface_manifold_stats["max_incidence"]
                println("Non-manifold edges: $nm_total ($nm_pct%, max incidence: $max_inc)")
            end
        else
            # Show simplified report when there are mismatches
            # Special case: if all edges are in both and 100% agreement, show minimal info
            if edges_only_in_A == 0 && edges_only_in_B == 0 && edges_in_both > 0 && agreement_rate >= 0.99
                println("Interface Analysis:")
                println("  Total edges with different triangulation: $edges_in_both")
                println("  Agreement rate: 100.0%")
            else
                # Show detailed breakdown for complex cases
                println("Interface Analysis:")
                println("  Total unique edges: $(length(symmetric_mismatches))")
                if edges_only_in_A > 0 || edges_only_in_B > 0 || edges_in_both > 0
                    println("  Distribution:")
                    if edges_only_in_A > 0
                        println("    • Only in A: $edges_only_in_A")
                    end
                    if edges_only_in_B > 0
                        println("    • Only in B: $edges_only_in_B")
                    end
                    if edges_in_both > 0
                        println("    • Different triangulation: $edges_in_both")
                    end
                end

                if edges_in_both > 0 && agreement_rate < 0.99
                    println("  Agreement rate: $(round(agreement_rate * 100, digits=1))%")
                end
            end

            # Non-manifold overview (interface scope)
            if !isempty(interface_manifold_stats)
                nm_total = interface_manifold_stats["non_manifold_count"]
                if nm_total > 0
                    nm_pct = round(interface_manifold_stats["percent_non_manifold"] * 100, digits=2)
                    max_inc = interface_manifold_stats["max_incidence"]
                    dist = interface_manifold_stats["overload_distribution"]

                    # Simplified reporting for common case: all non-manifold edges are in both
                    if get(dist, "A_only", 0) == 0 && get(dist, "B_only", 0) == 0 && get(dist, "both", 0) == nm_total
                        println("  Non-manifold edges: $nm_total ($nm_pct%, max incidence: $max_inc)")
                    else
                        # Detailed reporting for complex cases
                        println("  Non-manifold edges: $nm_total ($nm_pct%)")
                        println("    • Max incidence: $max_inc")
                        println("    • Distribution: A-only=$(get(dist, "A_only", 0)), B-only=$(get(dist, "B_only", 0)), both=$(get(dist, "both", 0))")
                    end
                end
            end

            # Add unified classification results merging both perspectives
            println("\nUNIFIED CLASSIFICATION RESULTS")
            print_unified_classification_results(classification_AB, classification_BA)
        end
    end

    return SymmetricClassificationResult(
        symmetric_mismatches,
        classification_AB,
        classification_BA,
        agreement_rate,
        length(symmetric_mismatches),
        edges_only_in_A,
        edges_only_in_B,
        edges_in_both,
        comparison_metrics,
        interface_manifold_stats
    )
end

# ============================================================================
# Helper functions for merged classification printing
# ============================================================================

"""
    print_classification_counts_perspective(classification, indent="")

Print classification counts for a perspective, showing only non-zero categories.
Uses the same format as the original classification printing.
"""
function print_classification_counts_perspective(classification::InterfaceClassification, indent::String="")
    all_mismatches = vcat(classification.mismatches_A, classification.mismatches_B)

    if length(all_mismatches) == 0
        println("$(indent)No mismatches found")
        return
    end

    println("$(indent)Classification complete:")

    # Count by type (using the classification's computed counts)
    if classification.t_junction_count > 0
        println("$(indent)  T-junctions: $(classification.t_junction_count)")
    end
    if classification.diagonal_count > 0
        println("$(indent)  Diagonal mismatches: $(classification.diagonal_count)")
    end
    if classification.refinement_count > 0
        println("$(indent)  Refinement differences: $(classification.refinement_count)")
    end
    if classification.quad_mismatch_count > 0
        println("$(indent)  Quad mismatches: $(classification.quad_mismatch_count)")
    end
    if classification.boundary_edge_count > 0
        println("$(indent)  Boundary edges: $(classification.boundary_edge_count)")
    end
    if classification.non_manifold_count > 0
        println("$(indent)  Non-manifold: $(classification.non_manifold_count)")
    end
    if classification.unshared_endpoint_count > 0
        println("$(indent)  Unshared endpoints: $(classification.unshared_endpoint_count)")
    end
    if classification.degenerate_edge_count > 0
        println("$(indent)  Degenerate edges: $(classification.degenerate_edge_count)")
    end
    if classification.source_edge_absent_count > 0
        println("$(indent)  Source edge absent: $(classification.source_edge_absent_count)")
    end
    if classification.quad_not_found_in_source_count > 0
        println("$(indent)  Quad not found in source: $(classification.quad_not_found_in_source_count)")
    end
    if classification.target_uses_finer_triangulation_count > 0
        println("$(indent)  Target uses finer triangulation: $(classification.target_uses_finer_triangulation_count)")
    end
    if classification.unknown_count > 0
        println("$(indent)  Unknown: $(classification.unknown_count)")
    end

    # Always show repair feasibility and complexity when there are mismatches
    println("$(indent)  Feasible for repair: $(classification.total_feasible) / $(length(all_mismatches))")
    println("$(indent)  Average complexity: $(round(classification.average_complexity, digits=2))")
end

"""
    print_classification_counts_perspective_swapped(classification, indent="")

Print classification counts for a swapped perspective (B→A), converting back to original A/B terminology.
"""
function print_classification_counts_perspective_swapped(classification::InterfaceClassification, indent::String="")
    all_mismatches = vcat(classification.mismatches_A, classification.mismatches_B)

    if length(all_mismatches) == 0
        println("$(indent)No mismatches found")
        return
    end

    println("$(indent)Classification complete:")

    # For swapped perspective, mismatches_A corresponds to edges missing in original A
    # and mismatches_B corresponds to edges missing in original B
    # But the type counts remain the same

    if classification.t_junction_count > 0
        println("$(indent)  T-junctions: $(classification.t_junction_count)")
    end
    if classification.diagonal_count > 0
        println("$(indent)  Diagonal mismatches: $(classification.diagonal_count)")
    end
    if classification.refinement_count > 0
        println("$(indent)  Refinement differences: $(classification.refinement_count)")
    end
    if classification.quad_mismatch_count > 0
        println("$(indent)  Quad mismatches: $(classification.quad_mismatch_count)")
    end
    if classification.boundary_edge_count > 0
        println("$(indent)  Boundary edges: $(classification.boundary_edge_count)")
    end
    if classification.non_manifold_count > 0
        println("$(indent)  Non-manifold: $(classification.non_manifold_count)")
    end
    if classification.unshared_endpoint_count > 0
        println("$(indent)  Unshared endpoints: $(classification.unshared_endpoint_count)")
    end
    if classification.degenerate_edge_count > 0
        println("$(indent)  Degenerate edges: $(classification.degenerate_edge_count)")
    end
    if classification.source_edge_absent_count > 0
        println("$(indent)  Source edge absent: $(classification.source_edge_absent_count)")
    end
    if classification.quad_not_found_in_source_count > 0
        println("$(indent)  Quad not found in source: $(classification.quad_not_found_in_source_count)")
    end
    if classification.target_uses_finer_triangulation_count > 0
        println("$(indent)  Target uses finer triangulation: $(classification.target_uses_finer_triangulation_count)")
    end
    if classification.unknown_count > 0
        println("$(indent)  Unknown: $(classification.unknown_count)")
    end

    # Always show repair feasibility and complexity when there are mismatches
    println("$(indent)  Feasible for repair: $(classification.total_feasible) / $(length(all_mismatches))")
    println("$(indent)  Average complexity: $(round(classification.average_complexity, digits=2))")
end

# ============================================================================
# Helper functions
# ============================================================================

"""
    print_unified_classification_results(classification_AB, classification_BA, indent="")

Print unified classification results merging both A→B and B→A perspectives
into a single compact report, avoiding redundancy.
"""
function print_unified_classification_results(
    classification_AB::InterfaceClassification,
    classification_BA::InterfaceClassification,
    indent::String=""
)
    ab_total = length(classification_AB.mismatches_A) + length(classification_AB.mismatches_B)
    ba_total = length(classification_BA.mismatches_A) + length(classification_BA.mismatches_B)

    if ab_total == 0 && ba_total == 0
        println("$(indent)No mismatches found in either perspective")
        return
    end

    # Show total from both perspectives (should be the same)
    total_mismatches = max(ab_total, ba_total)
    println("$(indent)Classification complete:")
    println("$(indent)  Total mismatches: $total_mismatches")

    # Count by type across both perspectives (use the higher count if they differ)
    t_junctions = max(classification_AB.t_junction_count, classification_BA.t_junction_count)
    diagonal = max(classification_AB.diagonal_count, classification_BA.diagonal_count)
    refinement = max(classification_AB.refinement_count, classification_BA.refinement_count)
    quad = max(classification_AB.quad_mismatch_count, classification_BA.quad_mismatch_count)
    boundary = max(classification_AB.boundary_edge_count, classification_BA.boundary_edge_count)
    non_manifold = max(classification_AB.non_manifold_count, classification_BA.non_manifold_count)
    unshared = max(classification_AB.unshared_endpoint_count, classification_BA.unshared_endpoint_count)
    degenerate = max(classification_AB.degenerate_edge_count, classification_BA.degenerate_edge_count)
    source_absent = max(classification_AB.source_edge_absent_count, classification_BA.source_edge_absent_count)
    quad_not_found = max(classification_AB.quad_not_found_in_source_count, classification_BA.quad_not_found_in_source_count)
    finer = max(classification_AB.target_uses_finer_triangulation_count, classification_BA.target_uses_finer_triangulation_count)
    unknown = max(classification_AB.unknown_count, classification_BA.unknown_count)

    # Only show non-zero categories
    if t_junctions > 0
        println("$(indent)  T-junctions: $t_junctions")
    end
    if diagonal > 0
        println("$(indent)  Diagonal mismatches: $diagonal")
    end
    if refinement > 0
        println("$(indent)  Refinement differences: $refinement")
    end
    if quad > 0
        println("$(indent)  Quad mismatches: $quad")
    end
    if boundary > 0
        println("$(indent)  Boundary edges: $boundary")
    end
    if non_manifold > 0
        println("$(indent)  Non-manifold: $non_manifold")
    end
    if unshared > 0
        println("$(indent)  Unshared endpoints: $unshared")
    end
    if degenerate > 0
        println("$(indent)  Degenerate edges: $degenerate")
    end
    if source_absent > 0
        println("$(indent)  Source edge absent: $source_absent")
    end
    if quad_not_found > 0
        println("$(indent)  Quad not found in source: $quad_not_found")
    end
    if finer > 0
        println("$(indent)  Target uses finer triangulation: $finer")
    end
    if unknown > 0
        println("$(indent)  Unknown: $unknown")
    end

    # Show repair feasibility and complexity (use averages from both perspectives)
    avg_feasible = round((classification_AB.total_feasible + classification_BA.total_feasible) / 2, digits=0)
    avg_complexity = round((classification_AB.average_complexity + classification_BA.average_complexity) / 2, digits=2)
    println("$(indent)  Feasible for repair: $avg_feasible / $total_mismatches")
    println("$(indent)  Average complexity: $avg_complexity")
end

"""
    find_corresponding_mismatch_in_swapped(mismatches, edge_key)

Find a mismatch with the given edge key in a list of mismatches from swapped perspective.

Since the topology was swapped, we need to match edges carefully by coordinates.
Returns the matching EdgeMismatch or nothing if not found.
"""
function find_corresponding_mismatch_in_swapped(
    mismatches::Vector{EdgeMismatch},
    edge_key::EdgeKey
)::Union{EdgeMismatch,Nothing}

    for mismatch in mismatches
        # Check if edge keys match (EdgeKey already handles canonical ordering)
        if mismatch.edge_key == edge_key
            return mismatch
        end
    end

    return nothing
end

"""
    update_symmetric_mismatches_with_strategies!(symmetric_mismatches, strategy_fn)

Update a vector of symmetric mismatches with repair strategies.

This is a placeholder for Phase 4 integration. The strategy function will
determine the optimal repair approach for each edge based on both perspectives.

# Arguments
- `symmetric_mismatches`: Vector of SymmetricEdgeMismatch to update
- `strategy_fn`: Function (sym_mismatch, constraints) -> (strategy, priority, reason)

# Example
```julia
result = classify_interface_mismatches_symmetric(topology)
update_symmetric_mismatches_with_strategies!(
    result.symmetric_mismatches,
    (sym, constraints) -> determine_repair_strategy(sym, constraints)
)
```
"""
function update_symmetric_mismatches_with_strategies!(
    symmetric_mismatches::Vector{SymmetricEdgeMismatch},
    strategy_fn::Function
)
    # This will be implemented in Phase 4
    # For now, just a placeholder
    @warn "Strategy selection not yet implemented (Phase 4)"
end

# ============================================================================
# Export utilities
# ============================================================================

"""
    export_symmetric_classification_json(result, output_file)

Export symmetric classification result to JSON for analysis.

Includes:
- All symmetric mismatches with both perspectives
- Agreement statistics
- Comparison metrics between perspectives
- Edge distribution analysis

# Example
```julia
result = classify_interface_mismatches_symmetric(topology)
export_symmetric_classification_json(result, "symmetric_classification.json")
```
"""
function export_symmetric_classification_json(
    result::SymmetricClassificationResult,
    output_file::String
)

    # Convert symmetric mismatches to JSON-friendly format
    function sym_mismatch_to_dict(sym::SymmetricEdgeMismatch)
        d = Dict(
            "edge" => [[sym.edge_key.node1...], [sym.edge_key.node2...]],
            "present_in_A" => sym.present_in_A,
            "present_in_B" => sym.present_in_B,
            "present_in_both" => sym.present_in_both,
            "agree_on_type" => sym.agree_on_type,
            "agree_on_feasibility" => sym.agree_on_feasibility,
            "repair_strategy" => string(sym.repair_strategy),
            "repair_priority" => round(sym.repair_priority, digits=3),
            "resolution_reason" => sym.resolution_reason
        )

        # Add A perspective if present
        if sym.classification_A_perspective !== nothing
            m_A = sym.classification_A_perspective
            d["classification_A"] = Dict(
                "mismatch_type" => string(m_A.mismatch_type),
                "repair_feasible" => m_A.repair_feasible,
                "complexity_score" => round(m_A.complexity_score, digits=3),
                "min_quality" => round(m_A.min_affected_triangle_quality, digits=3)
            )
        end

        # Add B perspective if present
        if sym.classification_B_perspective !== nothing
            m_B = sym.classification_B_perspective
            d["classification_B"] = Dict(
                "mismatch_type" => string(m_B.mismatch_type),
                "repair_feasible" => m_B.repair_feasible,
                "complexity_score" => round(m_B.complexity_score, digits=3),
                "min_quality" => round(m_B.min_affected_triangle_quality, digits=3)
            )
        end

        return d
    end

    report = Dict(
        "interface" => Dict(
            "pidA" => result.classification_AB.topology.pidA,
            "pidB" => result.classification_AB.topology.pidB
        ),
        "summary" => Dict(
            "total_unique_edges" => result.total_unique_edges,
            "edges_only_in_A" => result.edges_only_in_A,
            "edges_only_in_B" => result.edges_only_in_B,
            "edges_in_both" => result.edges_in_both,
            "agreement_rate" => round(result.agreement_rate, digits=4)
        ),
        "comparison_metrics" => result.comparison_metrics,
        "interface_manifold_stats" => result.interface_manifold_stats,
        "symmetric_mismatches" => [
            sym_mismatch_to_dict(sym)
            for sym in result.symmetric_mismatches[1:min(200, length(result.symmetric_mismatches))]
        ]
    )

    open(output_file, "w") do io
        write_json(io, report, 0)
    end

    println("Symmetric classification exported to: $output_file")
    return output_file
end

"""
    print_symmetric_classification_summary(result)

Print a human-readable summary of symmetric classification results.

# Example
```julia
result = classify_interface_mismatches_symmetric(topology)
print_symmetric_classification_summary(result)
```
"""
function print_symmetric_classification_summary(result::SymmetricClassificationResult)
    println("\n" * "="^70)
    println("SYMMETRIC CLASSIFICATION SUMMARY")
    println("="^70)

    println("\nEdge Distribution:")
    println("  Total unique edges: $(result.total_unique_edges)")
    println("  • Only in A: $(result.edges_only_in_A)")
    println("  • Only in B: $(result.edges_only_in_B)")
    println("  • In both: $(result.edges_in_both)")

    if result.edges_in_both > 0
        println("\nAgreement Analysis:")
        println("  Agreement rate: $(round(result.agreement_rate * 100, digits=1))%")

        agree_type = count(m -> m.agree_on_type, result.symmetric_mismatches)
        agree_feas = count(m -> m.agree_on_feasibility, result.symmetric_mismatches)

        println("  • Agree on mismatch type: $agree_type")
        println("  • Agree on feasibility: $agree_feas")
    end

    println("\nStrategy Distribution:")
    counts = count_by_strategy(result.symmetric_mismatches)
    println("  • Use A: $(counts[:use_A])")
    println("  • Use B: $(counts[:use_B])")
    println("  • Compromise: $(counts[:compromise])")
    println("  • Skip: $(counts[:skip])")
    println("  • Pending: $(get(counts, :pending, 0))")

    # Also echo interface-local non-manifold statistics if present
    if haskey(result.interface_manifold_stats, "non_manifold_count")
        println("\nInterface Non-manifold Summary:")
        nm = result.interface_manifold_stats
        println("  • Non-manifold edges: $(nm["non_manifold_count"]) of $(nm["total_unique_edges"]) (max incidence=$(nm["max_incidence"]))")
        dist = nm["overload_distribution"]
        println("  • Distribution: A-only=$(get(dist, "A_only", 0)), B-only=$(get(dist, "B_only", 0)), both=$(get(dist, "both", 0))")
    end

    println("="^70)
end

# ============================================================================
# Interface-local non-manifold statistics
# ============================================================================

"""
    compute_interface_nonmanifold_stats(topology::InterfaceTopology) -> Dict{String,Any}

Compute non-manifold edge statistics limited to the interface region before any repair.
We build an edge incidence map for both sides (A and B) using edge maps in the topology.
An edge is considered non-manifold if its total incidence across both sides exceeds 2.

Returns a dictionary with counts, incidence histogram, distribution, and examples.
"""
function compute_interface_nonmanifold_stats(topology::InterfaceTopology)::Dict{String,Any}
    # Build combined edge incidence from edge->face indices maps
    incidence = Dict{EdgeKey,Int}()
    sideA = Set{EdgeKey}()
    sideB = Set{EdgeKey}()

    for (ek, faces) in topology.edges_A
        incidence[ek] = get(incidence, ek, 0) + length(faces)
        push!(sideA, ek)
    end
    for (ek, faces) in topology.edges_B
        incidence[ek] = get(incidence, ek, 0) + length(faces)
        push!(sideB, ek)
    end

    if isempty(incidence)
        return Dict{String,Any}()
    end

    # Identify non-manifold edges (>2 incidences)
    nonmanifold = [(ek, c) for (ek, c) in incidence if c > 2]
    nn = length(nonmanifold)
    total_unique = length(incidence)
    max_incidence = nn > 0 ? maximum(c for (_, c) in nonmanifold) : 0

    # Histogram of incidences (3,4,5,...)
    hist = Dict{Int,Int}()
    for (_, c) in nonmanifold
        hist[c] = get(hist, c, 0) + 1
    end

    # Where do they occur relative to sides?
    dist = Dict{String,Int}("A_only" => 0, "B_only" => 0, "both" => 0)
    for (ek, _) in nonmanifold
        inA = ek in sideA
        inB = ek in sideB
        if inA && inB
            dist["both"] += 1
        elseif inA
            dist["A_only"] += 1
        elseif inB
            dist["B_only"] += 1
        end
    end

    # Compact examples for printing
    examples = String[]
    for (ek, c) in nonmanifold[1:min(10, nn)]
        push!(examples, "[($(ek.node1[1]),$(ek.node1[2]),$(ek.node1[3]))—($(ek.node2[1]),$(ek.node2[2]),$(ek.node2[3]))]×$(c)")
    end

    return Dict{String,Any}(
        "total_unique_edges" => total_unique,
        "non_manifold_count" => nn,
        "percent_non_manifold" => total_unique > 0 ? nn / total_unique : 0.0,
        "max_incidence" => max_incidence,
        "incidence_histogram" => hist,
        "overload_distribution" => dist,
        "top_examples" => examples
    )
end
