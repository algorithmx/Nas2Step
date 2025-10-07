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
    
    if verbose
        println("\n" * "="^70)
        println("SYMMETRIC EDGE CLASSIFICATION")
        println("="^70)
        println("Interface: PID $(topology.pidA) ↔ PID $(topology.pidB)")
        println("Strategy: Classify from BOTH perspectives and preserve results")
        println()
    end
    
    # ========================================================================
    # Step 1: Classify from A→B perspective
    # ========================================================================
    if verbose
        println("-"^70)
        println("PERSPECTIVE 1: A (PID=$(topology.pidA)) → B (PID=$(topology.pidB))")
        println("-"^70)
    end
    
    classification_AB = classify_interface_mismatches(topology, tol=tol, debug=false)
    
    if verbose
        println("  Edges only in A (missing in B): $(length(classification_AB.mismatches_B))")
        println("  Edges only in B (missing in A): $(length(classification_AB.mismatches_A))")
        println("  Total classified: $(length(classification_AB.mismatches_A) + length(classification_AB.mismatches_B))")
    end
    
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
        topology.conformity_ratio
    )
    
    # ========================================================================
    # Step 3: Classify from B→A perspective
    # ========================================================================
    if verbose
        println("\n" * "-"^70)
        println("PERSPECTIVE 2: B (PID=$(topology.pidB)) → A (PID=$(topology.pidA))")
        println("-"^70)
    end
    
    classification_BA = classify_interface_mismatches(topology_swapped, tol=tol, debug=false)
    
    if verbose
        println("  Edges only in B (missing in A): $(length(classification_BA.mismatches_B))")
        println("  Edges only in A (missing in B): $(length(classification_BA.mismatches_A))")
        println("  Total classified: $(length(classification_BA.mismatches_A) + length(classification_BA.mismatches_B))")
    end
    
    # ========================================================================
    # Step 4: Build symmetric mismatch map
    # ========================================================================
    if verbose
        println("\n" * "-"^70)
        println("BUILDING SYMMETRIC MISMATCHES")
        println("-"^70)
    end
    
    # Map edge keys to symmetric mismatches
    symmetric_map = Dict{EdgeKey, SymmetricEdgeMismatch}()
    
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
    
    if verbose
        println("\n" * "="^70)
        println("SYMMETRIC CLASSIFICATION SUMMARY")
        println("="^70)
        println("Total unique edges: $(length(symmetric_mismatches))")
        println("  • Only in A: $edges_only_in_A")
        println("  • Only in B: $edges_only_in_B")
        println("  • In both (different triangulation): $edges_in_both")
        
        if edges_in_both > 0
            println("\nAgreement analysis (edges in both):")
            println("  • Edges with both classifications: $(length(edges_with_both_classifications))")
            println("  • Agreement rate: $(round(agreement_rate * 100, digits=1))%")
        end
        
        println("="^70)
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
        comparison_metrics
    )
end

# ============================================================================
# Helper functions
# ============================================================================

"""
    find_corresponding_mismatch_in_swapped(mismatches, edge_key)

Find a mismatch with the given edge key in a list of mismatches from swapped perspective.

Since the topology was swapped, we need to match edges carefully by coordinates.
Returns the matching EdgeMismatch or nothing if not found.
"""
function find_corresponding_mismatch_in_swapped(
    mismatches::Vector{EdgeMismatch},
    edge_key::EdgeKey
)::Union{EdgeMismatch, Nothing}
    
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
    
    println("="^70)
end
