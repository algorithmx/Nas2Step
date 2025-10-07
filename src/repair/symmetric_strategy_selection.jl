# Symmetric Strategy Selection Implementation
# Phase 4: Local per-edge repair strategy selection
#
# This module implements strategy selection for symmetric edge mismatches.
# Each edge independently chooses the optimal repair approach based on
# classifications from BOTH A→B and B→A perspectives.

# ============================================================================
# Core strategy selection function
# ============================================================================

"""
    determine_repair_strategy(sym_mismatch, constraints; thresholds=default_thresholds())

Determine the optimal repair strategy for a symmetric edge mismatch.

Analyzes both A→B and B→A perspectives to select the best local approach:
- `:use_A`: Use A's triangulation for this edge in unified mesh
- `:use_B`: Use B's triangulation for this edge in unified mesh
- `:compromise`: Synthesize a compromise triangulation (not yet implemented)
- `:skip`: Cannot repair this edge (infeasible or irreconcilable)

# Decision Tree

## Case 1: Edge only in A
- If A's classification is feasible and acceptable quality → `:use_A`
- Otherwise → `:skip`

## Case 2: Edge only in B
- If B's classification is feasible and acceptable quality → `:use_B`
- Otherwise → `:skip`

## Case 3: Edge in both (different triangulation)
- Compare quality predictions from both perspectives
- If one significantly better (>20% improvement) → use better
- If quality similar, prefer feasible over infeasible
- If both feasible/infeasible, use A by default (tiebreaker)

# Arguments
- `sym_mismatch::SymmetricEdgeMismatch`: Edge with both perspectives
- `constraints::BoundaryConstraints`: Constraints from both sides
- `thresholds::QualityThresholds`: Quality requirements (optional)

# Returns
`(strategy::Symbol, priority::Float64, reason::String)`

# Example
```julia
strategy, priority, reason = determine_repair_strategy(sym, constraints)

if strategy == :use_A
    println("Using A's triangulation: \$reason")
elseif strategy == :skip
    println("Skipping edge: \$reason")
end
```
"""
function determine_repair_strategy(
    sym_mismatch::SymmetricEdgeMismatch,
    constraints::BoundaryConstraints;
    thresholds::QualityThresholds = default_thresholds()
)::Tuple{Symbol, Float64, String}
    
    m_A = sym_mismatch.classification_A_perspective
    m_B = sym_mismatch.classification_B_perspective
    
    # ========================================================================
    # Case 1: Edge only in A (present in A, missing in B)
    # ========================================================================
    if sym_mismatch.present_in_A && !sym_mismatch.present_in_B
        if m_A === nothing
            return (:skip, 0.0, "No A classification available (unexpected)")
        end
        
        # Check if A's edge is repairable
        if !m_A.repair_feasible
            return (:skip, 0.0, "A's classification marked as infeasible: $(m_A.mismatch_type)")
        end
        
        # Check mismatch type
        if m_A.mismatch_type in [BOUNDARY_EDGE, NON_MANIFOLD, UNSHARED_ENDPOINT]
            return (:skip, 0.0, "A has irreparable mismatch: $(m_A.mismatch_type)")
        end
        
        # Check quality
        quality_A = m_A.min_affected_triangle_quality
        min_quality_threshold = thresholds.min_angle / 60.0  # Normalized
        
        if quality_A < min_quality_threshold
            return (:skip, 0.1, "A's quality too poor ($(round(quality_A*60, digits=1))° < $(thresholds.min_angle)°)")
        end
        
        # Check constraints
        has_violation, violations = check_constraint_violations(sym_mismatch.edge_key, constraints)
        if has_violation
            return (:skip, 0.1, "Constraint violations: $(join(violations, ", "))")
        end
        
        # Use A - compute priority based on mismatch type
        priority = compute_priority_for_mismatch_type(m_A.mismatch_type, quality_A)
        reason = "Edge from A ($(m_A.mismatch_type)) with quality $(round(quality_A*60, digits=1))°"
        
        return (:use_A, priority, reason)
    end
    
    # ========================================================================
    # Case 2: Edge only in B (present in B, missing in A)
    # ========================================================================
    if sym_mismatch.present_in_B && !sym_mismatch.present_in_A
        if m_B === nothing
            return (:skip, 0.0, "No B classification available (unexpected)")
        end
        
        # Check if B's edge is repairable
        if !m_B.repair_feasible
            return (:skip, 0.0, "B's classification marked as infeasible: $(m_B.mismatch_type)")
        end
        
        # Check mismatch type
        if m_B.mismatch_type in [BOUNDARY_EDGE, NON_MANIFOLD, UNSHARED_ENDPOINT]
            return (:skip, 0.0, "B has irreparable mismatch: $(m_B.mismatch_type)")
        end
        
        # Check quality
        quality_B = m_B.min_affected_triangle_quality
        min_quality_threshold = thresholds.min_angle / 60.0
        
        if quality_B < min_quality_threshold
            return (:skip, 0.1, "B's quality too poor ($(round(quality_B*60, digits=1))° < $(thresholds.min_angle)°)")
        end
        
        # Check constraints
        has_violation, violations = check_constraint_violations(sym_mismatch.edge_key, constraints)
        if has_violation
            return (:skip, 0.1, "Constraint violations: $(join(violations, ", "))")
        end
        
        # Use B - compute priority
        priority = compute_priority_for_mismatch_type(m_B.mismatch_type, quality_B)
        reason = "Edge from B ($(m_B.mismatch_type)) with quality $(round(quality_B*60, digits=1))°"
        
        return (:use_B, priority, reason)
    end
    
    # ========================================================================
    # Case 3: Edge in both (different triangulation)
    # ========================================================================
    if sym_mismatch.present_in_both
        if m_A === nothing || m_B === nothing
            return (:skip, 0.0, "Missing classification for edge in both (unexpected)")
        end
        
        quality_A = m_A.min_affected_triangle_quality
        quality_B = m_B.min_affected_triangle_quality
        feasible_A = m_A.repair_feasible
        feasible_B = m_B.repair_feasible
        
        # Check if both are infeasible
        if !feasible_A && !feasible_B
            return (:skip, 0.0, "Both perspectives marked as infeasible")
        end
        
        # Check irreparable types
        irreparable_A = m_A.mismatch_type in [BOUNDARY_EDGE, NON_MANIFOLD, UNSHARED_ENDPOINT]
        irreparable_B = m_B.mismatch_type in [BOUNDARY_EDGE, NON_MANIFOLD, UNSHARED_ENDPOINT]
        
        if irreparable_A && irreparable_B
            return (:skip, 0.0, "Both have irreparable types: A=$(m_A.mismatch_type), B=$(m_B.mismatch_type)")
        end
        
        if irreparable_A && !irreparable_B
            return (:use_B, 0.6, "A is irreparable ($(m_A.mismatch_type)), using B")
        end
        
        if irreparable_B && !irreparable_A
            return (:use_A, 0.6, "B is irreparable ($(m_B.mismatch_type)), using A")
        end
        
        # Compare quality (both are repairable types)
        quality_ratio_A_to_B = quality_A / (quality_B + 1e-10)
        
        # Significant quality difference (>20%)
        if quality_ratio_A_to_B > 1.2
            priority = compute_priority_for_mismatch_type(m_A.mismatch_type, quality_A)
            return (:use_A, priority, "A has significantly better quality ($(round(quality_A*60, digits=1))° vs $(round(quality_B*60, digits=1))°)")
        elseif quality_ratio_A_to_B < 0.833  # 1/1.2 ≈ 0.833
            priority = compute_priority_for_mismatch_type(m_B.mismatch_type, quality_B)
            return (:use_B, priority, "B has significantly better quality ($(round(quality_B*60, digits=1))° vs $(round(quality_A*60, digits=1))°)")
        end
        
        # Quality similar - check feasibility
        if feasible_A && !feasible_B
            priority = compute_priority_for_mismatch_type(m_A.mismatch_type, quality_A)
            return (:use_A, priority, "Quality similar, A is feasible, B is not")
        elseif feasible_B && !feasible_A
            priority = compute_priority_for_mismatch_type(m_B.mismatch_type, quality_B)
            return (:use_B, priority, "Quality similar, B is feasible, A is not")
        end
        
        # Both feasible, quality similar - check mismatch type priority
        type_priority_A = mismatch_type_priority(m_A.mismatch_type)
        type_priority_B = mismatch_type_priority(m_B.mismatch_type)
        
        if type_priority_A > type_priority_B
            priority = compute_priority_for_mismatch_type(m_A.mismatch_type, quality_A)
            return (:use_A, priority, "Both feasible, A has higher type priority ($(m_A.mismatch_type) > $(m_B.mismatch_type))")
        elseif type_priority_B > type_priority_A
            priority = compute_priority_for_mismatch_type(m_B.mismatch_type, quality_B)
            return (:use_B, priority, "Both feasible, B has higher type priority ($(m_B.mismatch_type) > $(m_A.mismatch_type))")
        end
        
        # Complete tie - use A as tiebreaker
        priority = compute_priority_for_mismatch_type(m_A.mismatch_type, quality_A)
        return (:use_A, priority, "Complete tie, using A as tiebreaker")
    end
    
    # Should never reach here
    return (:skip, 0.0, "Unexpected case: edge not in A, B, or both")
end

# ============================================================================
# Helper functions for priority computation
# ============================================================================

"""
    mismatch_type_priority(mtype::MismatchType)::Float64

Return priority value for a mismatch type. Higher = more important to repair.

Priority ordering:
- T_JUNCTION: 1.0 (highest - critical for conformity)
- DIAGONAL: 0.8 (high - improves consistency)
- REFINEMENT: 0.6 (medium - improves resolution matching)
- QUAD_MISMATCH: 0.4 (lower - structural issue)
- Others: 0.2 (lowest)
"""
function mismatch_type_priority(mtype::MismatchType)::Float64
    if mtype == T_JUNCTION
        return 1.0
    elseif mtype == DIAGONAL
        return 0.8
    elseif mtype == REFINEMENT
        return 0.6
    elseif mtype == QUAD_MISMATCH
        return 0.4
    else
        return 0.2
    end
end

"""
    compute_priority_for_mismatch_type(mtype::MismatchType, quality::Float64)::Float64

Compute repair priority combining mismatch type and quality.

Priority = base_priority * quality_factor

Where:
- base_priority is from mismatch_type_priority()
- quality_factor is normalized quality (0.0 to 1.0)

This ensures high-priority mismatches with good quality get repaired first.
"""
function compute_priority_for_mismatch_type(mtype::MismatchType, quality::Float64)::Float64
    base_priority = mismatch_type_priority(mtype)
    quality_factor = clamp(quality, 0.0, 1.0)
    
    # Combined priority with slight quality weighting
    return base_priority * (0.7 + 0.3 * quality_factor)
end

# ============================================================================
# Batch strategy selection
# ============================================================================

"""
    apply_strategy_selection!(result::SymmetricClassificationResult, constraints::BoundaryConstraints; thresholds=default_thresholds(), verbose=true)

Apply strategy selection to all symmetric mismatches in a classification result.

Updates each `SymmetricEdgeMismatch` with:
- `repair_strategy`: Selected strategy (`:use_A`, `:use_B`, `:skip`)
- `repair_priority`: Priority for repair ordering (0.0 to 1.0)
- `resolution_reason`: Explanation of why strategy was chosen

# Note
Since `SymmetricEdgeMismatch` is immutable, this function creates NEW instances
with updated strategies and replaces them in the result.

# Arguments
- `result::SymmetricClassificationResult`: Classification result to update
- `constraints::BoundaryConstraints`: Constraints from both sides
- `thresholds::QualityThresholds`: Quality requirements (optional)
- `verbose::Bool`: Print progress information (optional)

# Returns
New `SymmetricClassificationResult` with updated strategies

# Example
```julia
result = classify_interface_mismatches_symmetric(topology)
result_with_strategies = apply_strategy_selection!(result, constraints)

# Now check strategies
for sym in result_with_strategies.symmetric_mismatches
    println("Edge: \$(sym.edge_key) → Strategy: \$(sym.repair_strategy)")
end
```
"""
function apply_strategy_selection!(
    result::SymmetricClassificationResult,
    constraints::BoundaryConstraints;
    thresholds::QualityThresholds = default_thresholds(),
    verbose::Bool = true
)::SymmetricClassificationResult
    
    if verbose
        println("\n" * "="^70)
        println("STRATEGY SELECTION")
        println("="^70)
        println("Processing $(length(result.symmetric_mismatches)) edges...")
    end
    
    # Create new symmetric mismatches with updated strategies
    updated_mismatches = SymmetricEdgeMismatch[]
    
    strategy_counts = Dict{Symbol,Int}(
        :use_A => 0,
        :use_B => 0,
        :compromise => 0,
        :skip => 0
    )
    
    for (idx, sym) in enumerate(result.symmetric_mismatches)
        # Determine strategy
        strategy, priority, reason = determine_repair_strategy(sym, constraints, thresholds=thresholds)
        
        # Create new SymmetricEdgeMismatch with updated strategy
        updated_sym = SymmetricEdgeMismatch(
            sym.edge_key,
            sym.classification_A_perspective,
            sym.classification_B_perspective,
            strategy,
            priority,
            reason
        )
        
        push!(updated_mismatches, updated_sym)
        strategy_counts[strategy] = get(strategy_counts, strategy, 0) + 1
        
        # Progress reporting
        if verbose && idx % 50 == 0
            println("  Progress: $idx/$(length(result.symmetric_mismatches))")
        end
    end
    
    if verbose
        println("\n" * "-"^70)
        println("STRATEGY SELECTION SUMMARY")
        println("-"^70)
        println("Total edges: $(length(updated_mismatches))")
        println("  • Use A: $(strategy_counts[:use_A])")
        println("  • Use B: $(strategy_counts[:use_B])")
        println("  • Compromise: $(strategy_counts[:compromise])")
        println("  • Skip: $(strategy_counts[:skip])")
        
        feasible = strategy_counts[:use_A] + strategy_counts[:use_B] + strategy_counts[:compromise]
        if length(updated_mismatches) > 0
            feasibility_rate = feasible / length(updated_mismatches) * 100
            println("\nFeasibility: $feasible/$(length(updated_mismatches)) ($(round(feasibility_rate, digits=1))%)")
        end
        println("="^70)
    end
    
    # Create new result with updated mismatches
    return SymmetricClassificationResult(
        updated_mismatches,
        result.classification_AB,
        result.classification_BA,
        result.agreement_rate,
        result.total_unique_edges,
        result.edges_only_in_A,
        result.edges_only_in_B,
        result.edges_in_both,
        result.comparison_metrics
    )
end

# ============================================================================
# Strategy analysis utilities
# ============================================================================

"""
    analyze_strategy_distribution(result::SymmetricClassificationResult)

Analyze the distribution of repair strategies and provide insights.

Returns a comprehensive report including:
- Strategy counts by type
- Quality statistics per strategy
- Feasibility breakdown
- Priority distribution

# Example
```julia
result = apply_strategy_selection!(result, constraints)
report = analyze_strategy_distribution(result)

println("Use A: \$(report[:use_A][:count])")
println("Avg priority: \$(report[:use_A][:avg_priority])")
```
"""
function analyze_strategy_distribution(result::SymmetricClassificationResult)::Dict{Symbol,Any}
    
    report = Dict{Symbol,Any}()
    
    for strategy in [:use_A, :use_B, :compromise, :skip]
        edges_with_strategy = filter(m -> m.repair_strategy == strategy, result.symmetric_mismatches)
        
        if isempty(edges_with_strategy)
            report[strategy] = Dict(
                :count => 0,
                :avg_priority => 0.0,
                :avg_quality => 0.0,
                :mismatch_types => Dict{String,Int}()
            )
            continue
        end
        
        # Compute statistics
        priorities = [m.repair_priority for m in edges_with_strategy]
        avg_priority = mean(priorities)
        
        # Quality statistics
        qualities = Float64[]
        for m in edges_with_strategy
            if strategy == :use_A && m.classification_A_perspective !== nothing
                push!(qualities, m.classification_A_perspective.min_affected_triangle_quality)
            elseif strategy == :use_B && m.classification_B_perspective !== nothing
                push!(qualities, m.classification_B_perspective.min_affected_triangle_quality)
            end
        end
        avg_quality = isempty(qualities) ? 0.0 : mean(qualities)
        
        # Mismatch type distribution
        type_counts = Dict{String,Int}()
        for m in edges_with_strategy
            if strategy == :use_A && m.classification_A_perspective !== nothing
                tname = string(m.classification_A_perspective.mismatch_type)
                type_counts[tname] = get(type_counts, tname, 0) + 1
            elseif strategy == :use_B && m.classification_B_perspective !== nothing
                tname = string(m.classification_B_perspective.mismatch_type)
                type_counts[tname] = get(type_counts, tname, 0) + 1
            end
        end
        
        report[strategy] = Dict(
            :count => length(edges_with_strategy),
            :avg_priority => avg_priority,
            :avg_quality => avg_quality,
            :mismatch_types => type_counts
        )
    end
    
    return report
end

"""
    print_strategy_analysis(result::SymmetricClassificationResult)

Print a detailed analysis of strategy distribution.

# Example
```julia
result = apply_strategy_selection!(result, constraints)
print_strategy_analysis(result)
```
"""
function print_strategy_analysis(result::SymmetricClassificationResult)
    report = analyze_strategy_distribution(result)
    
    println("\n" * "="^70)
    println("STRATEGY DISTRIBUTION ANALYSIS")
    println("="^70)
    
    for strategy in [:use_A, :use_B, :compromise, :skip]
        stats = report[strategy]
        println("\n$(uppercase(string(strategy))):")
        println("  Count: $(stats[:count])")
        println("  Avg Priority: $(round(stats[:avg_priority], digits=3))")
        println("  Avg Quality: $(round(stats[:avg_quality] * 60, digits=1))°")
        
        if !isempty(stats[:mismatch_types])
            println("  Mismatch Types:")
            for (tname, count) in sort(collect(stats[:mismatch_types]), by=x->-x[2])
                println("    • $tname: $count")
            end
        end
    end
    
    println("\n" * "="^70)
end

"""
    get_high_priority_edges(result::SymmetricClassificationResult, min_priority::Float64=0.7)

Get edges with priority above threshold, sorted by descending priority.

# Example
```julia
high_priority = get_high_priority_edges(result, 0.8)
println("Found \$(length(high_priority)) high-priority edges")
```
"""
function get_high_priority_edges(
    result::SymmetricClassificationResult,
    min_priority::Float64 = 0.7
)::Vector{SymmetricEdgeMismatch}
    
    high_priority = filter(
        m -> m.repair_priority >= min_priority && m.repair_strategy != :skip,
        result.symmetric_mismatches
    )
    
    return sort(high_priority, by = m -> m.repair_priority, rev = true)
end
