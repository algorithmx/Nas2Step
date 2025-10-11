# Repair strategy generation for interface conformity
# Phase 2: Plan surgical mesh repairs

using .CoordinateKeys: convert_to_float

# ============================================================================
# Repair plan data structures
# ============================================================================

"""
Plan for inserting a single edge into a triangle.
"""
struct EdgeInsertionPlan
    target_triangle::Int              # Triangle index that will be split
    insert_edge::EdgeKey              # Edge to insert
    
    # Operation type
    split_type::Symbol                # :bisect | :trisect | :quadrisect | :quad_retriangulation
    
    # Nodes to create/reuse
    new_nodes::Vector{NTuple{3,Float64}}  # Coordinates of nodes to create
    existing_nodes::Vector{Int}           # Node IDs to reuse from original mesh
    
    # Triangles for quad retriangulation
    old_triangles::Vector{NTuple{9,Float64}}      # Original triangles to delete
    replacement_triangles::Vector{NTuple{9,Float64}}  # New triangles to add
    
    # Quality validation
    min_angle_before::Float64
    min_angle_after::Float64
    quality_acceptable::Bool
    
    # Dependencies
    depends_on::Vector{Int}  # Other insertion indices that must happen first
    
    # Constraint check
    violates_constraints::Bool
    constraint_violations::Vector{String}
    
    # Overall feasibility (computed from quality_acceptable && !violates_constraints)
    is_feasible::Bool
end

"""
Complete repair plan for an interface.
"""
struct RepairPlan
    interface_pair::Tuple{Int,Int}
    repair_direction::Symbol  # :subdivide_A | :subdivide_B | :mutual
    
    # Ordered list of edge insertions
    insertion_sequence::Vector{EdgeInsertionPlan}
    
    # Summary statistics
    total_edges_to_insert::Int
    total_triangles_to_split::Int
    total_nodes_to_add::Int
    
    # Quality predictions
    predicted_min_quality::Float64
    predicted_max_quality_loss::Float64
    
    # Feasibility
    is_feasible::Bool
    feasibility_issues::Vector{String}
    
    # Metadata
    topology::InterfaceTopology
    classification::InterfaceClassification
    constraints::BoundaryConstraints
end

"""
Quality thresholds for repair validation.
"""
struct QualityThresholds
    min_angle::Float64      # degrees, e.g., 15.0
    max_angle::Float64      # degrees, e.g., 150.0
    max_aspect_ratio::Float64  # e.g., 10.0
end

function default_thresholds()
    return QualityThresholds(10.0, 170.0, 20.0)
end

# ============================================================================
# Feasibility assessment helpers
# ============================================================================

const PLAN_FEASIBLE_RATIO_DEFAULT = 0.8

# ============================================================================
# Helper functions for feasibility assessment
# ============================================================================

"""
Compute quantiles (p25, median, p75) from a vector of values.
Handles empty vectors gracefully.
"""
function _compute_quantiles(x::Vector{Float64})
    if isempty(x)
        return (p25=0.0, p50=0.0, p75=0.0)
    end
    
    s = sort(x)
    n = length(s)
    
    # Compute indices for 25th, 50th (median), and 75th percentiles
    idx25 = max(1, min(n, Int(floor(0.25*(n-1))) + 1))
    idx50a = ((n + 1) >>> 1)  # floor((n+1)/2)
    idx50b = isodd(n) ? idx50a : idx50a + (idx50a < n ? 1 : 0)
    med = isodd(n) ? s[idx50a] : 0.5*(s[idx50a] + s[idx50b])
    idx75 = max(1, min(n, Int(floor(0.75*(n-1))) + 1))
    
    return (p25=s[idx25], p50=med, p75=s[idx75])
end

"""
Categorize plans into feasible and infeasible groups, with detailed breakdown.

Returns:
- n_feasible: Plans that pass both quality and constraint checks
- n_infeasible: Plans that fail at least one check
- n_only_quality: Plans that fail only quality (constraints OK)
- n_only_constraints: Plans that fail only constraints (quality OK) - these are actionable!
- n_both: Plans that fail both quality and constraints
"""
function _categorize_plans(plans::Vector{EdgeInsertionPlan})
    n_feasible = count(p -> p.quality_acceptable && !p.violates_constraints, plans)
    n_infeasible = length(plans) - n_feasible
    
    # Detailed breakdown to understand what's blocking repairs
    n_only_quality = count(p -> (!p.quality_acceptable) && !p.violates_constraints, plans)
    n_only_constraints = count(p -> p.quality_acceptable && p.violates_constraints, plans)
    n_both = count(p -> (!p.quality_acceptable) && p.violates_constraints, plans)
    
    return (
        n_feasible = n_feasible,
        n_infeasible = n_infeasible,
        n_only_quality = n_only_quality,
        n_only_constraints = n_only_constraints,
        n_both = n_both,
        n_constraints = n_only_constraints + n_both,
        n_poor_quality = n_only_quality + n_both
    )
end

"""
Collect all constraint violations from plans and track which plans hit each reason.

Returns:
- reason_counts: How many times each constraint violation occurred
- plans_per_reason: Which plan indices hit each constraint
"""
function _collect_constraint_violations(plans::Vector{EdgeInsertionPlan})
    reason_counts = Dict{String,Int}()
    plans_per_reason = Dict{String,Vector{Int}}()
    
    for (idx, p) in enumerate(plans)
        if p.violates_constraints
            for reason in p.constraint_violations
                # Count occurrences
                reason_counts[reason] = get(reason_counts, reason, 0) + 1
                
                # Track which plans have this violation
                if !haskey(plans_per_reason, reason)
                    plans_per_reason[reason] = Int[]
                end
                push!(plans_per_reason[reason], idx)
            end
        end
    end
    
    return reason_counts, plans_per_reason
end

"""
Analyze the impact of each constraint violation.

For each constraint, determine:
- How many plans it affects
- What percentage of plans it blocks
- How many plans would become feasible if ONLY this constraint were removed
- Whether it's truly blocking progress (actionable) or a symptom of quality issues
"""
function _analyze_constraint_impact(
    plans::Vector{EdgeInsertionPlan},
    reason_counts::Dict{String,Int},
    plans_per_reason::Dict{String,Vector{Int}},
    n_plans::Int,
    n_infeasible::Int
)
    constraint_impact = Dict{String,Any}()
    
    for (reason, count) in reason_counts
        # Calculate percentages to understand scope of impact
        pct_of_infeasible = n_infeasible > 0 ? (count / n_infeasible) * 100 : 0.0
        pct_of_all = n_plans > 0 ? (count / n_plans) * 100 : 0.0
        
        # Critical question: How many plans would become feasible if we fixed ONLY this constraint?
        # These are "actionable" constraints worth investigating
        plans_with_this_reason = Set(get(plans_per_reason, reason, Int[]))
        would_be_feasible = 0
        
        for idx in plans_with_this_reason
            p = plans[idx]
            # A plan would become feasible if:
            # 1. Quality is already acceptable
            # 2. This is the ONLY constraint violation
            if p.quality_acceptable && length(p.constraint_violations) == 1
                would_be_feasible += 1
            end
        end
        
        constraint_impact[reason] = (
            count = count,
            pct_of_infeasible = pct_of_infeasible,
            pct_of_all = pct_of_all,
            would_be_feasible = would_be_feasible,
            is_blocking = would_be_feasible > 0  # True if fixing this would unlock plans
        )
    end
    
    return constraint_impact
end

"""
Analyze constraint violations with detailed geometric breakdown.

Breaks down constraint blocking by:
- Geometric reasons (min angle, max angle, aspect ratio)
- Topological issues (no affected triangles, invalid indices)
- Boundary constraints (locked edges/nodes)
"""
function _analyze_constraint_breakdown(
    plans::Vector{EdgeInsertionPlan},
    thresholds::QualityThresholds
)
    breakdown = Dict{String,Any}(
        "topological" => Dict{String,Int}(),
        "boundary" => Dict{String,Int}(),
        "geometric" => Dict{String,Any}(
            "min_angle_violations" => 0,
            "max_angle_violations" => 0,
            "aspect_ratio_violations" => 0,
            "quality_stats" => Dict{String,Vector{Float64}}()
        ),
        "by_mismatch_type" => Dict{Symbol,Int}()
    )
    
    # Track quality metrics for constrained plans
    min_angles_constrained = Float64[]
    max_angles_constrained = Float64[]
    
    for plan in plans
        if plan.violates_constraints
            # Categorize each violation
            for violation in plan.constraint_violations
                if contains(violation, "No affected triangles")
                    breakdown["topological"]["no_affected_triangles"] = 
                        get(breakdown["topological"], "no_affected_triangles", 0) + 1
                elseif contains(violation, "Invalid triangle index")
                    breakdown["topological"]["invalid_triangle_index"] = 
                        get(breakdown["topological"], "invalid_triangle_index", 0) + 1
                elseif contains(violation, "Cannot plan")
                    breakdown["topological"]["cannot_plan_operation"] = 
                        get(breakdown["topological"], "cannot_plan_operation", 0) + 1
                elseif contains(violation, "locked") || contains(violation, "external") || contains(violation, "feature")
                    breakdown["boundary"]["locked_edges_or_nodes"] = 
                        get(breakdown["boundary"], "locked_edges_or_nodes", 0) + 1
                elseif contains(violation, "corner")
                    breakdown["boundary"]["corner_nodes"] = 
                        get(breakdown["boundary"], "corner_nodes", 0) + 1
                end
            end
            
            # Track quality metrics for plans with constraint violations
            push!(min_angles_constrained, plan.min_angle_after)
            
            # Detect implicit geometric constraint violations from quality thresholds
            # (even if not explicitly listed in constraint_violations)
            min_angle_deg = plan.min_angle_after * 60.0  # Convert from normalized to degrees
            if min_angle_deg < thresholds.min_angle
                breakdown["geometric"]["min_angle_violations"] += 1
            end
        end
    end
    
    # Store quality statistics
    breakdown["geometric"]["quality_stats"]["min_angles_constrained"] = min_angles_constrained
    
    return breakdown
end

"""
Compute quality statistics (before/after repair) for all plans.

Returns triangle quality metrics including min, percentiles, and changes.
"""
function _compute_quality_stats(plans::Vector{EdgeInsertionPlan}, thresholds::QualityThresholds)
    q_before = [p.min_angle_before for p in plans]
    q_after = [p.min_angle_after for p in plans]
    deltas = [p.min_angle_after - p.min_angle_before for p in plans]
    
    quantiles_after = _compute_quantiles(q_after)
    quantiles_delta = _compute_quantiles(deltas)
    
    return (
        min_before = isempty(q_before) ? 0.0 : minimum(q_before),
        min_after = isempty(q_after) ? 0.0 : minimum(q_after),
        min_threshold_norm = thresholds.min_angle / 60.0,
        p25_after = quantiles_after.p25,
        median_after = quantiles_after.p50,
        p75_after = quantiles_after.p75,
        delta_min = isempty(deltas) ? 0.0 : minimum(deltas),
        delta_max = isempty(deltas) ? 0.0 : maximum(deltas),
        delta_med = quantiles_delta.p50
    )
end

"""
Generate human-readable issue descriptions based on feasibility assessment.
"""
function _generate_issue_descriptions(
    n_plans::Int,
    n_feasible::Int,
    total_needed::Int,
    required_threshold::Float64,
    required_min::Int,
    is_feasible::Bool,
    categories
)
    issues = String[]
    
    # Check for coverage gap
    if n_plans < total_needed
        gap = total_needed - n_plans
        push!(issues, "Coverage gap: generated $(n_plans)/$(total_needed) plans (short by $gap)")
    end
    
    # Check if feasibility threshold is met
    if !is_feasible
        # Detailed phrasing for analysis
        push!(issues, "Feasible plans $(n_feasible)/$(total_needed) below required > $(round(required_threshold, digits=2)) (need at least $(required_min))")
        # Legacy phrasing for compatibility
        push!(issues, "Only $(n_feasible)/$(total_needed) plans are feasible")
    end
    
    # Report constraint violations with breakdown
    if categories.n_constraints > 0
        push!(issues, "$(categories.n_constraints) plans blocked by constraints ($(categories.n_only_constraints) only-constraints, $(categories.n_both) both)")
        push!(issues, "$(categories.n_constraints) plans violate constraints")  # Legacy
    end
    
    # Report quality failures with breakdown
    if categories.n_poor_quality > 0
        push!(issues, "$(categories.n_poor_quality) plans fail quality ($(categories.n_only_quality) only-quality, $(categories.n_both) both)")
        push!(issues, "$(categories.n_poor_quality) plans would create poor quality triangles")  # Legacy
    end
    
    return issues
end

# ============================================================================
# Main feasibility assessment function
# ============================================================================

"""
    assess_interface_feasibility(plans, total_needed, thresholds; ratio=0.8)

Compute comprehensive feasibility assessment from a set of EdgeInsertionPlans.

Returns a named tuple with:
- feasible_plans / infeasible_plans: Filtered plan lists
- is_feasible: Overall decision (true if enough plans are feasible)
- issues: Human-readable issue descriptions
- reason_counts: Count of each constraint violation type
- constraint_impact: Detailed impact analysis for each constraint
- constraint_breakdown: Categorized breakdown (topological, boundary, geometric)
- qstats: Quality statistics (before/after, percentiles, deltas)
- counts: Detailed counters (feasible, infeasible, categories)
- ratios: Coverage and feasibility ratios
- split_types: Count of plans by operation type
- top_reasons: Top constraint violations by frequency
"""
function assess_interface_feasibility(
    plans::Vector{EdgeInsertionPlan},
    total_needed::Int,
    thresholds::QualityThresholds;
    ratio::Float64 = PLAN_FEASIBLE_RATIO_DEFAULT,
)
    # Separate plans into feasible and infeasible groups
    feasible_plans = filter(p -> p.quality_acceptable && !p.violates_constraints, plans)
    infeasible_plans = filter(p -> !p.quality_acceptable || p.violates_constraints, plans)
    
    n_plans = length(plans)
    n_feasible = length(feasible_plans)
    n_infeasible = length(infeasible_plans)
    
    # Perform detailed constraint breakdown analysis
    constraint_breakdown = _analyze_constraint_breakdown(plans, thresholds)
    
    # Categorize plans to understand failure modes
    categories = _categorize_plans(plans)
    
    # Determine if we have enough feasible plans to meet threshold
    # Special case: if nothing is needed (total_needed == 0), automatically feasible
    required_threshold = total_needed * ratio
    required_min = total_needed == 0 ? 0 : (floor(Int, required_threshold) + 1)
    is_feasible = total_needed == 0 ? true : (n_feasible > required_threshold)
    
    # Analyze constraint violations in detail
    reason_counts, plans_per_reason = _collect_constraint_violations(plans)
    constraint_impact = _analyze_constraint_impact(
        plans, reason_counts, plans_per_reason, n_plans, n_infeasible
    )
    
    # Get top reasons sorted by frequency (up to 5)
    top_reasons = sort(collect(reason_counts); by = x -> -x[2])[1:min(5, length(reason_counts))]
    
    # Compute quality statistics
    qstats = _compute_quality_stats(plans, thresholds)
    
    # Count split types to understand operation distribution
    split_types = Dict{Symbol,Int}()
    for p in plans
        split_types[p.split_type] = get(split_types, p.split_type, 0) + 1
    end
    
    # Generate human-readable issue descriptions
    issues = _generate_issue_descriptions(
        n_plans, n_feasible, total_needed, required_threshold, 
        required_min, is_feasible, categories
    )
    
    # Return comprehensive assessment
    return (
        feasible_plans = feasible_plans,
        infeasible_plans = infeasible_plans,
        is_feasible = is_feasible,
        issues = issues,
        reason_counts = reason_counts,
        constraint_impact = constraint_impact,
        constraint_breakdown = constraint_breakdown,  # NEW: Detailed breakdown by category
        qstats = qstats,
        counts = (
            total_needed = total_needed,
            generated = n_plans,
            feasible = n_feasible,
            infeasible = n_infeasible,
            only_quality = categories.n_only_quality,
            only_constraints = categories.n_only_constraints,
            both = categories.n_both,
            constraints = categories.n_constraints,
            poor_quality = categories.n_poor_quality,
            required_min = required_min,
        ),
        ratios = (
            coverage = total_needed == 0 ? 1.0 : n_plans / total_needed,
            feasible = total_needed == 0 ? 1.0 : n_feasible / total_needed,
            required = ratio,
        ),
        split_types = split_types,
        top_reasons = top_reasons,
    )
end

# ============================================================================
# Step 2.1: Determine dominant side
# ============================================================================

"""
    compute_edge_density(topology, use_pid)

Compute edge density (edges per unit area) for one side of interface.
"""
function compute_edge_density(topology::InterfaceTopology, use_pid::Symbol)::Float64
    edges = use_pid == :A ? topology.edges_A : topology.edges_B
    area = compute_interface_area(topology, use_pid)
    
    if area < 1e-10
        return 0.0
    end
    
    return length(edges) / area
end

"""
    determine_dominant_side(topology, constraints)

Determine which side should provide the target edge pattern.
Returns :A or :B.
"""
function determine_dominant_side(topology::InterfaceTopology, 
                                constraints::BoundaryConstraints)::Symbol
    
    # Heuristic 1: Edge density
    density_A = compute_edge_density(topology, :A)
    density_B = compute_edge_density(topology, :B)
    
    println("  Edge densities: A=$(round(density_A, digits=4)), B=$(round(density_B, digits=4))")
    
    if density_A > density_B * 1.2
        println("  → A is finer ($(round(density_A/density_B, digits=2))x), should subdivide B to match A")
        return :A
    elseif density_B > density_A * 1.2
        println("  → B is finer ($(round(density_B/density_A, digits=2))x), should subdivide A to match B")
        return :B
    end
    
    # Heuristic 2: External connectivity (fewer locked edges = easier to modify)
    # Calculate the ratio of external edges to the interface edges
    interface_edges_A = length(topology.edges_A)
    interface_edges_B = length(topology.edges_B)
    locked_ratio_A = min(1.0, length(constraints.pidA_external_edges) / max(1, interface_edges_A))
    locked_ratio_B = min(1.0, length(constraints.pidB_external_edges) / max(1, interface_edges_B))
    
    println("  Locked edge ratios: A=$(round(locked_ratio_A*100, digits=1))%, B=$(round(locked_ratio_B*100, digits=1))%")
    
    # Prefer modifying the side with fewer external constraints
    if locked_ratio_A < locked_ratio_B * 0.8
        println("  → A has fewer constraints, prefer modifying A")
        return :B  # B is target pattern, modify A
    elseif locked_ratio_B < locked_ratio_A * 0.8
        println("  → B has fewer constraints, prefer modifying B")
        return :A  # A is target pattern, modify B
    end
    
    # Tie-breaker: prefer smaller side (fewer triangles to modify)
    if topology.total_faces_A < topology.total_faces_B
        println("  → A is smaller, prefer modifying A")
        return :B
    else
        println("  → B is smaller, prefer modifying B")
        return :A
    end
end

# ============================================================================
# Step 2.2: Generate edge insertion plans
# ============================================================================

"""
    plan_quad_retriangulation(quad_vertices, diagonal_edge)

Plan to retriangulate a quad using the specified diagonal edge.
Returns plan with 2 new triangles.
"""
function plan_quad_retriangulation(quad_vertices::Vector{NTuple{3,Float64}}, diagonal_edge::EdgeKey)
    if length(quad_vertices) != 4
        return nothing
    end
    
    # The diagonal connects 2 opposite corners
    # EdgeKey stores integer coordinates, but quad_vertices are float coordinates
    # Convert EdgeKey integer coordinates back to float for comparison
    corner1 = convert_to_float(diagonal_edge.node1)
    corner2 = convert_to_float(diagonal_edge.node2)

    # Find the other 2 vertices (not on the diagonal)
    # NOTE: EdgeKey coordinates are rounded to 4 digits, but quad_vertices may not be.
    # We need to check both exact match and rounded match.
    round_coord(c) = (round(c[1], digits=4), round(c[2], digits=4), round(c[3], digits=4))

    other_vertices = filter(quad_vertices) do v
        # Not a corner if it matches neither corner1 nor corner2
        not_corner1 = (v != corner1) && (round_coord(v) != corner1)
        not_corner2 = (v != corner2) && (round_coord(v) != corner2)
        not_corner1 && not_corner2
    end
    
    if length(other_vertices) != 2
        return nothing
    end
    
    vertex3 = other_vertices[1]
    vertex4 = other_vertices[2]
    
    # Create 2 new triangles using the desired diagonal
    # Triangle 1: (corner1, corner2, vertex3)
    # Triangle 2: (corner1, corner2, vertex4)
    new_triangles = [
        (corner1[1], corner1[2], corner1[3],
         corner2[1], corner2[2], corner2[3],
         vertex3[1], vertex3[2], vertex3[3]),
        (corner1[1], corner1[2], corner1[3],
         corner2[1], corner2[2], corner2[3],
         vertex4[1], vertex4[2], vertex4[3])
    ]
    
    return (
        type = :quad_retriangulation,
        new_nodes = NTuple{3,Float64}[],  # No new nodes needed
        new_triangles = new_triangles
    )
end

"""
    plan_triangle_split(triangle, edge_to_insert)

Plan how to split a triangle to insert an edge.
Returns split plan with new triangles.
"""
function plan_triangle_split(triangle::Triangle, edge_to_insert::EdgeKey)
    # Check if edge endpoints are in the triangle
    (has_edge, opposite_vertex) = TriangleHasEdge(triangle, edge_to_insert)

    if !has_edge
        # Edge doesn't fully belong to this triangle
        return nothing
    end
    @assert opposite_vertex !== nothing
    
    # Create two new triangles by bisecting
    # Triangle 1: (edge.node1, edge.node2, opposite)
    # This is essentially just reconfirming the edge exists
    opposite_vertex = opposite_vertex==1 ? triangle.coord1 : (opposite_vertex==2 ? triangle.coord2 : triangle.coord3)

    # Convert EdgeKey integer coordinates back to float for consistency
    edge_start = convert_to_float(edge_to_insert.node1)
    edge_end = convert_to_float(edge_to_insert.node2)

    new_triangles = [
        (edge_start[1], edge_start[2], edge_start[3],
         edge_end[1], edge_end[2], edge_end[3],
         opposite_vertex[1], opposite_vertex[2], opposite_vertex[3])
    ]
    
    return (
        type = :bisect,
        new_nodes = NTuple{3,Float64}[],  # No new nodes needed for edge that exists
        new_triangles = new_triangles
    )
end

"""
    generate_edge_insertion_plan(mismatch, topology, constraints, thresholds)

Generate detailed plan for inserting one missing edge.
Always returns an EdgeInsertionPlan (may be marked as :failed if planning fails).
"""
function generate_edge_insertion_plan(
    mismatch::EdgeMismatch,
    topology::InterfaceTopology,
    constraints::BoundaryConstraints,
    thresholds::QualityThresholds
)::EdgeInsertionPlan

    # Determine target side (where edge needs to be inserted)
    target_faces = mismatch.should_be_in == topology.pidA ? topology.faces_A : topology.faces_B
    
    # Check constraint violations
    has_violation, violations = check_constraint_violations(mismatch.edge_key, constraints)
    
    # Handle non-repairable mismatch types with clear error messages
    if mismatch.mismatch_type == QUAD_MISMATCH
        return EdgeInsertionPlan(
            1, mismatch.edge_key, :failed,
            NTuple{3,Float64}[], Int[], NTuple{3,Float64}[], NTuple{9,Float64}[],
            0.0, 0.0, false, Int[], true,
            ["Quad topology mismatch: same 4 vertices but different boundary edges - cannot repair with retriangulation"],
            false
        )
    elseif mismatch.mismatch_type == BOUNDARY_EDGE
        return EdgeInsertionPlan(
            1, mismatch.edge_key, :failed,
            NTuple{3,Float64}[], Int[], NTuple{3,Float64}[], NTuple{9,Float64}[],
            0.0, 0.0, false, Int[], true,
            ["Boundary edge: appears in only 1 triangle - not a true interface edge, cannot repair"],
            false
        )
    elseif mismatch.mismatch_type == NON_MANIFOLD
        return EdgeInsertionPlan(
            1, mismatch.edge_key, :failed,
            NTuple{3,Float64}[], Int[], NTuple{3,Float64}[], NTuple{9,Float64}[],
            0.0, 0.0, false, Int[], true,
            ["Non-manifold edge: shared by more than 2 triangles - topology is broken, cannot repair"],
            false
        )
    elseif mismatch.mismatch_type == UNSHARED_ENDPOINT
        return EdgeInsertionPlan(
            1, mismatch.edge_key, :failed,
            NTuple{3,Float64}[], Int[], NTuple{3,Float64}[], NTuple{9,Float64}[],
            0.0, 0.0, false, Int[], true,
            ["Unshared endpoint: one or both edge endpoints not in shared vertex set - fundamental topology issue"],
            false
        )
    elseif mismatch.mismatch_type == DIAGONAL
        return handle_diagonal_mismatch(mismatch, has_violation, violations, target_faces, thresholds)
    else
        return handle_Tjunction_other_mismatch(mismatch, has_violation, violations, target_faces, thresholds)
    end
end

function handle_Tjunction_other_mismatch(
        mismatch::EdgeMismatch, 
        has_violation::Bool, 
        violations::Vector{String},
        target_faces::Vector{Triangle}, 
        thresholds::QualityThresholds
    )::EdgeInsertionPlan

    # Helper function to create a failed plan with consistent structure
    function create_failed_plan(target_idx::Int, reason::String, min_angle::Float64=0.0)::EdgeInsertionPlan
        return EdgeInsertionPlan(
            target_idx,
            mismatch.edge_key,
            :failed,
            NTuple{3,Float64}[],
            Int[],
            NTuple{9,Float64}[],
            NTuple{9,Float64}[],
            min_angle,
            0.0,
            false,
            Int[],
            true,
            [reason],
            false
        )
    end

    # For T-junctions and other types, we expect exactly one affected triangle
    if isempty(mismatch.affected_triangles)
        return create_failed_plan(1, "No affected triangles found for edge insertion")
    end
    
    # Take first affected triangle
    target_tri_idx = mismatch.affected_triangles[1]
    if target_tri_idx < 1 || target_tri_idx > length(target_faces)
        return create_failed_plan(target_tri_idx, "Invalid triangle index: $target_tri_idx")
    end
    
    target_triangle = target_faces[target_tri_idx]
    
    # Compute current quality
    min_angle_before = compute_triangle_quality(target_triangle)
    
    # Plan the split
    split_plan = plan_triangle_split(target_triangle, mismatch.edge_key)
    
    if split_plan === nothing
        return create_failed_plan(target_tri_idx, "Cannot plan split for this triangle", min_angle_before)
    end
    
    # Quality after is essentially the same
    min_angle_after = min_angle_before
    quality_acceptable = min_angle_after >= thresholds.min_angle / 60.0
    is_feasible = quality_acceptable && !has_violation
    
    # Extract old triangle coordinates
    old_tri = target_triangle
    old_tri_flat = (old_tri.coord1[1], old_tri.coord1[2], old_tri.coord1[3],
                   old_tri.coord2[1], old_tri.coord2[2], old_tri.coord2[3],
                   old_tri.coord3[1], old_tri.coord3[2], old_tri.coord3[3])
    
    # Ensure we always return EdgeInsertionPlan (never nothing)
    return EdgeInsertionPlan(
        target_tri_idx,
        mismatch.edge_key,
        split_plan.type,
        split_plan.new_nodes,
        Int[],  # Will be filled during execution
        [old_tri_flat],  # old_triangles
        split_plan.new_triangles,  # replacement_triangles
        min_angle_before,
        min_angle_after,
        quality_acceptable,
        Int[],
        has_violation,
        violations,
        is_feasible
    )
end


function handle_diagonal_mismatch(
        mismatch::EdgeMismatch, 
        has_violation::Bool, 
        violations::Vector{String},
        target_faces::Vector{Triangle}, 
        thresholds::QualityThresholds
    )::EdgeInsertionPlan

    # Helper function to create a failed plan with consistent structure
    function create_failed_plan(target_idx::Int, reason::String, min_angle::Float64=0.0)::EdgeInsertionPlan
        return EdgeInsertionPlan(
            target_idx,
            mismatch.edge_key,
            :failed,
            NTuple{3,Float64}[],
            Int[],
            NTuple{9,Float64}[],
            NTuple{9,Float64}[],
            min_angle,
            0.0,
            false,
            Int[],
            true,
            [reason],
            false
        )
    end

    # For diagonal mismatches, use quad retriangulation
    if isempty(mismatch.quad_vertices) || isempty(mismatch.triangles_to_replace)
        return create_failed_plan(1, "Quad vertices or triangles to replace not found")
    end
    
    # Compute quality of triangles being replaced
    min_angle_before = 1.0
    for tri_idx in mismatch.triangles_to_replace
        if tri_idx > 0 && tri_idx <= length(target_faces)
            quality = compute_triangle_quality(target_faces[tri_idx])
            min_angle_before = min(min_angle_before, quality)
        end
    end
    
    # Plan quad retriangulation
    retri_plan = plan_quad_retriangulation(mismatch.quad_vertices, mismatch.edge_key)
    
    if retri_plan === nothing
        first_tri_idx = isempty(mismatch.triangles_to_replace) ? 1 : mismatch.triangles_to_replace[1]
        return create_failed_plan(first_tri_idx, "Cannot plan quad retriangulation", min_angle_before)
    end
    
    # Quality after should be similar (same vertices, different diagonal)
    min_angle_after = min_angle_before  # Approximation
    quality_acceptable = min_angle_after >= thresholds.min_angle / 60.0
    is_feasible = quality_acceptable && !has_violation
    
    # Extract old triangle coordinates
    old_triangles = NTuple{9,Float64}[]
    for tri_idx in mismatch.triangles_to_replace
        if tri_idx > 0 && tri_idx <= length(target_faces)
            tri = target_faces[tri_idx]
            tri_flat = (tri.coord1[1], tri.coord1[2], tri.coord1[3],
                        tri.coord2[1], tri.coord2[2], tri.coord2[3],
                        tri.coord3[1], tri.coord3[2], tri.coord3[3])
            push!(old_triangles, tri_flat)
        end
    end
    
    # Ensure we always return EdgeInsertionPlan (never nothing)
    return EdgeInsertionPlan(
        mismatch.triangles_to_replace[1],  # First triangle to replace
        mismatch.edge_key,
        retri_plan.type,
        retri_plan.new_nodes,
        Int[],  # Will be filled during execution
        old_triangles,  # old_triangles
        retri_plan.new_triangles,  # replacement_triangles
        min_angle_before,
        min_angle_after,
        quality_acceptable,
        Int[],
        has_violation,
        violations,
        is_feasible
    )
end

# ============================================================================
# Step 2.3: Generate and validate complete repair plan
# ============================================================================

"""
    generate_repair_plan(topology, classification, constraints; thresholds)

Generate complete repair plan for an interface.
"""
function generate_repair_plan(
    topology::InterfaceTopology,
    classification::InterfaceClassification,
    constraints::BoundaryConstraints;
    thresholds::QualityThresholds = default_thresholds()
)::RepairPlan
    
    println("\nGenerating repair plan for PID=$(topology.pidA) ↔ PID=$(topology.pidB)...")
    
    # Step 2.1: Determine dominant side
    dominant = determine_dominant_side(topology, constraints)
    
    # Determine repair direction
    repair_direction = if dominant == :A
        :subdivide_B  # B needs edges from A
    else
        :subdivide_A  # A needs edges from B
    end
    
    println("  Repair direction: $repair_direction")
    
    # Step 2.2: Generate insertion plans
    # Note: mismatches_A contains edges missing in A (present only in B)
    #       mismatches_B contains edges missing in B (present only in A)
    # If we want to subdivide_A, we need to add edges that A is missing (i.e., mismatches_A)
    # If we want to subdivide_B, we need to add edges that B is missing (i.e., mismatches_B)
    mismatches_to_fix = if repair_direction == :subdivide_A
        classification.mismatches_A
    else
        classification.mismatches_B
    end
    
    # Override if the selected direction has no mismatches but the other does
    if isempty(mismatches_to_fix)
        if repair_direction == :subdivide_A && !isempty(classification.mismatches_B)
            # We wanted to subdivide A but A has no missing edges. B has edges A doesn't have.
            # Switch to subdividing B to add A's edges to B
            println("  Note: No edges missing in A, switching to subdivide B with edges from A")
            repair_direction = :subdivide_B
            mismatches_to_fix = classification.mismatches_B
        elseif repair_direction == :subdivide_B && !isempty(classification.mismatches_A)
            # We wanted to subdivide B but B has no missing edges. A has edges B doesn't have.
            # Switch to subdividing A to add B's edges to A
            println("  Note: No edges missing in B, switching to subdivide A with edges from B")
            repair_direction = :subdivide_A
            mismatches_to_fix = classification.mismatches_A
        end
    end
    
    println("  Generating insertion plans for $(length(mismatches_to_fix)) edges...")
    
    insertion_plans = EdgeInsertionPlan[]
    for (idx, mismatch) in enumerate(mismatches_to_fix)
        plan = generate_edge_insertion_plan(mismatch, topology, constraints, thresholds)
        push!(insertion_plans, plan)
        
        if idx % 50 == 0
            println("    Progress: $idx/$(length(mismatches_to_fix))")
        end
    end
    
    # Count plans that were marked as :failed during planning
    failed_count = count(p -> p.split_type == :failed, insertion_plans)
    if failed_count > 0
        println("  Generated $(length(insertion_plans)) insertion plans ($(failed_count) marked as :failed)")
    else
        println("  Generated $(length(insertion_plans)) insertion plans")
    end
    
    # Compute statistics
    total_edges = length(insertion_plans)
    total_triangles = length(insertion_plans) > 0 ? length(unique(p.target_triangle for p in insertion_plans)) : 0
    total_nodes = length(insertion_plans) > 0 ? sum(length(p.new_nodes) for p in insertion_plans) : 0
    
    # Quality predictions
    qualities = [p.min_angle_after for p in insertion_plans if p.quality_acceptable]
    predicted_min_quality = isempty(qualities) ? 0.0 : minimum(qualities)
    
    quality_losses = [p.min_angle_before - p.min_angle_after for p in insertion_plans]
    predicted_max_loss = isempty(quality_losses) ? 0.0 : maximum(quality_losses)
    
    # Feasibility assessment (modular)
    assessment = assess_interface_feasibility(insertion_plans, length(mismatches_to_fix), thresholds)

    # Streamlined feasibility output
    println("\n  Feasibility: $(assessment.counts.feasible)/$(assessment.counts.total_needed) feasible (need ≥$(assessment.counts.required_min)) → $(assessment.is_feasible ? "✓ PASS" : "✗ FAIL")")
    
    if !assessment.is_feasible
        println("    Breakdown: $(assessment.counts.only_quality) quality-only, $(assessment.counts.only_constraints) constraints-only, $(assessment.counts.both) both")
        
        # Detailed constraint blocking analysis - show ALL root causes
        if assessment.counts.constraints > 0
            bd = assessment.constraint_breakdown
            
            println("    Blocking criteria:")
            
            # Show topological blockers
            if !isempty(bd["topological"])
                for (issue, count) in sort(collect(bd["topological"]); by=x->-x[2])
                    pct = round(count/assessment.counts.total_needed*100, digits=1)
                    println("      • Topological: $(replace(issue, "_" => " ")) - $(count)/$(assessment.counts.total_needed) plans ($(pct)%)")
                end
            end
            
            # Show boundary constraint blockers
            if !isempty(bd["boundary"])
                for (issue, count) in sort(collect(bd["boundary"]); by=x->-x[2])
                    pct = round(count/assessment.counts.total_needed*100, digits=1)
                    println("      • Boundary: $(replace(issue, "_" => " ")) - $(count)/$(assessment.counts.total_needed) plans ($(pct)%)")
                end
            end
            
            # Show geometric/quality blockers with details
            geom = bd["geometric"]
            if geom["min_angle_violations"] > 0
                angles = geom["quality_stats"]["min_angles_constrained"]
                pct = round(geom["min_angle_violations"]/assessment.counts.total_needed*100, digits=1)
                
                if !isempty(angles)
                    angles_deg = [a * 60.0 for a in angles]
                    sorted_angles = sort(angles_deg)
                    min_deg = round(minimum(angles_deg), digits=1)
                    med_deg = round(sorted_angles[(length(sorted_angles)+1)÷ 2], digits=1)
                    max_deg = round(maximum(angles_deg), digits=1)
                    
                    println("      • Geometric: min_angle < $(thresholds.min_angle)° - $(geom["min_angle_violations"])/$(assessment.counts.total_needed) plans ($(pct)%)")
                    println("        Quality distribution: min=$(min_deg)°, median=$(med_deg)°, max=$(max_deg)°")
                else
                    println("      • Geometric: min_angle < $(thresholds.min_angle)° - $(geom["min_angle_violations"])/$(assessment.counts.total_needed) plans ($(pct)%)")
                end
            end
        end
    end
    
    return RepairPlan(
        (topology.pidA, topology.pidB),
        repair_direction,
        insertion_plans,
        total_edges,
        total_triangles,
        total_nodes,
        predicted_min_quality,
        predicted_max_loss,
        assessment.is_feasible,
        assessment.issues,
        topology,
        classification,
        constraints
    )
end

"""
    export_repair_plan_json(plan, output_file)

Export repair plan to JSON for inspection.
"""
function export_repair_plan_json(plan::RepairPlan, output_file::String)
    
    function insertion_to_dict(p::EdgeInsertionPlan)
        # Convert old_triangles from NTuple{9,Float64} to arrays of 3 points
        old_triangles_array = []
        for tri_flat in p.old_triangles
            tri = [
                [tri_flat[1], tri_flat[2], tri_flat[3]],
                [tri_flat[4], tri_flat[5], tri_flat[6]],
                [tri_flat[7], tri_flat[8], tri_flat[9]]
            ]
            push!(old_triangles_array, tri)
        end
        
        # Convert replacement_triangles from NTuple{9,Float64} to arrays of 3 points
        new_triangles_array = []
        for tri_flat in p.replacement_triangles
            tri = [
                [tri_flat[1], tri_flat[2], tri_flat[3]],
                [tri_flat[4], tri_flat[5], tri_flat[6]],
                [tri_flat[7], tri_flat[8], tri_flat[9]]
            ]
            push!(new_triangles_array, tri)
        end
        
        return Dict(
            "target_triangle" => p.target_triangle,
            "edge" => [[p.insert_edge.node1...], [p.insert_edge.node2...]],
            "plan_type" => string(p.split_type),
            "split_type" => string(p.split_type),
            "new_nodes" => [[n...] for n in p.new_nodes],
            "old_triangles" => old_triangles_array,
            "triangles_to_replace" => old_triangles_array,  # Alias for compatibility
            "replacement_triangles" => new_triangles_array,
            "new_triangles" => new_triangles_array,  # Alias for compatibility
            "num_replacement_triangles" => length(p.replacement_triangles),
            "min_angle_before" => round(p.min_angle_before, digits=3),
            "min_angle_after" => round(p.min_angle_after, digits=3),
            "quality_acceptable" => p.quality_acceptable,
            "violates_constraints" => p.violates_constraints,
            "constraint_violations" => p.constraint_violations,
            "is_feasible" => p.is_feasible
        )
    end
    
    # Sample plans for export (limit size)
    feasible_plans = filter(p -> p.quality_acceptable && !p.violates_constraints, plan.insertion_sequence)
    infeasible_plans = filter(p -> !p.quality_acceptable || p.violates_constraints, plan.insertion_sequence)
    
    report = Dict(
        "interface" => Dict(
            "pidA" => plan.interface_pair[1],
            "pidB" => plan.interface_pair[2]
        ),
        "strategy" => Dict(
            "repair_direction" => string(plan.repair_direction),
            "total_edges_to_insert" => plan.total_edges_to_insert,
            "total_triangles_to_split" => plan.total_triangles_to_split,
            "total_nodes_to_add" => plan.total_nodes_to_add
        ),
        "quality_prediction" => Dict(
            "predicted_min_quality" => round(plan.predicted_min_quality, digits=3),
            "predicted_max_quality_loss" => round(plan.predicted_max_quality_loss, digits=3)
        ),
        "feasibility" => Dict(
            "is_feasible" => plan.is_feasible,
            "feasible_plans" => length(feasible_plans),
            "infeasible_plans" => length(infeasible_plans),
            "issues" => plan.feasibility_issues
        ),
        "insertion_sequence" => [insertion_to_dict(p) for p in plan.insertion_sequence],
        "sample_plans" => Dict(
            "feasible" => [insertion_to_dict(p) for p in feasible_plans[1:min(20, length(feasible_plans))]],
            "infeasible" => [insertion_to_dict(p) for p in infeasible_plans[1:min(20, length(infeasible_plans))]]
        )
    )
    
    open(output_file, "w") do io
        write_json(io, report, 0)
    end
    
    println("Repair plan exported to: $output_file")
    return output_file
end

# ============================================================================
# Bidirectional repair planning (strategic improvement)
# ============================================================================

"""
    generate_repair_plan_bidirectional(topology, classification, constraints; thresholds, verbose=true)

Generate repair plans for BOTH possible source-target orientations and select the better one.

This addresses the strategic limitation where the source-target direction is fixed upfront
based on heuristics. By trying both orientations, we:
  1. Avoid premature commitment to a poor direction
  2. Discover repairs that would be missed by heuristic-only selection
  3. Choose the direction with more feasible repairs

Returns:
  - best_plan: The RepairPlan with the most feasible repairs
  - tried_both: Boolean indicating if both directions were actually attempted
  - direction_info: Dictionary with comparison metrics
"""
function generate_repair_plan_bidirectional(
    topology::InterfaceTopology,
    classification::InterfaceClassification,
    constraints::BoundaryConstraints;
    thresholds::QualityThresholds = default_thresholds(),
    verbose::Bool = true
)::Tuple{RepairPlan, Bool, Dict{String,Any}}
    
    if verbose
        println("\n" * "="^70)
        println("BIDIRECTIONAL REPAIR PLANNING")
        println("="^70)
        println("Strategy: Try both possible source ↔ target orientations")
        println("          and select the direction with more feasible repairs")
        println()
    end
    
    # Check if there are mismatches in both directions
    has_A_mismatches = !isempty(classification.mismatches_A)
    has_B_mismatches = !isempty(classification.mismatches_B)
    
    if verbose
        println("Mismatch distribution:")
        println("  • Edges missing in A (present in B): $(length(classification.mismatches_A))")
        println("  • Edges missing in B (present in A): $(length(classification.mismatches_B))")
    end
    
    # If only one direction has mismatches, skip bidirectional and use unidirectional
    if !has_A_mismatches && !has_B_mismatches
        if verbose
            println("\n→ No edge mismatches detected. Skipping repair planning.")
            println("="^70)
        end
        # Return a dummy feasible plan
        dummy_plan = generate_repair_plan(topology, classification, constraints, thresholds=thresholds)
        return (dummy_plan, false, Dict("reason" => "no_mismatches"))
    elseif !has_A_mismatches || !has_B_mismatches
        if verbose
            println("\n→ Only one direction has mismatches. Using unidirectional planning.")
            println("="^70)
        end
        unidirectional_plan = generate_repair_plan(topology, classification, constraints, thresholds=thresholds)
        direction_info = Dict(
            "reason" => "unidirectional_only",
            "chosen_direction" => string(unidirectional_plan.repair_direction)
        )
        return (unidirectional_plan, false, direction_info)
    end
    
    # Both directions have mismatches - try both!
    if verbose
        println("\n→ Both directions have mismatches. Testing both orientations...\n")
    end
    
    # -------------------------------------------------------------------------
    # Direction 1: subdivide_A (insert edges from B into A)
    # -------------------------------------------------------------------------
    if verbose
        println("-"^70)
        println("DIRECTION 1: Subdivide A (modify PID=$(topology.pidA))")
        println("-"^70)
    end
    
    # Temporarily suppress verbose output for individual plan generation
    plan_A = generate_repair_plan(topology, classification, constraints, thresholds=thresholds)
    
    feasible_A = count(p -> p.is_feasible, plan_A.insertion_sequence)
    total_A = length(plan_A.insertion_sequence)
    
    if verbose
        println("  Result: $feasible_A/$total_A feasible plans")
        if plan_A.is_feasible
            println("  Status: ✓ FEASIBLE")
        else
            println("  Status: ✗ INFEASIBLE")
        end
    end
    
    # -------------------------------------------------------------------------
    # Direction 2: subdivide_B (insert edges from A into B)
    # -------------------------------------------------------------------------
    if verbose
        println("\n" * "-"^70)
        println("DIRECTION 2: Subdivide B (modify PID=$(topology.pidB))")
        println("-"^70)
    end
    
    # Force the opposite direction by temporarily swapping mismatches
    # Create a new classification with swapped mismatch lists (auto-compute statistics)
    cls_swapped = InterfaceClassification(
        topology,
        classification.mismatches_B,  # Swap: now A's list contains B's mismatches
        classification.mismatches_A   # Swap: now B's list contains A's mismatches
    )
    
    plan_B = generate_repair_plan(topology, cls_swapped, constraints, thresholds=thresholds)
    
    feasible_B = count(p -> p.is_feasible, plan_B.insertion_sequence)
    total_B = length(plan_B.insertion_sequence)
    
    if verbose
        println("  Result: $feasible_B/$total_B feasible plans")
        if plan_B.is_feasible
            println("  Status: ✓ FEASIBLE")
        else
            println("  Status: ✗ INFEASIBLE")
        end
    end
    
    # -------------------------------------------------------------------------
    # Selection: Choose the direction with more feasible repairs
    # -------------------------------------------------------------------------
    if verbose
        println("\n" * "="^70)
        println("DIRECTION SELECTION")
        println("="^70)
    end
    
    direction_info = Dict(
        "tried_both" => true,
        "direction_A" => Dict(
            "repair_direction" => string(plan_A.repair_direction),
            "feasible" => feasible_A,
            "total" => total_A,
            "feasibility_ratio" => total_A > 0 ? feasible_A / total_A : 0.0,
            "is_feasible" => plan_A.is_feasible
        ),
        "direction_B" => Dict(
            "repair_direction" => string(plan_B.repair_direction),
            "feasible" => feasible_B,
            "total" => total_B,
            "feasibility_ratio" => total_B > 0 ? feasible_B / total_B : 0.0,
            "is_feasible" => plan_B.is_feasible
        )
    )
    
    # Selection criteria (in order of priority):
    # 1. Prefer the direction that meets the feasibility threshold (is_feasible == true)
    # 2. If both or neither meet threshold, prefer the one with more feasible plans
    # 3. If tied, prefer the one with fewer total repairs needed (less invasive)
    
    best_plan = nothing
    chosen_direction = ""
    selection_reason = ""
    
    if plan_A.is_feasible && !plan_B.is_feasible
        best_plan = plan_A
        chosen_direction = "Direction 1 (subdivide A)"
        selection_reason = "Only Direction 1 meets feasibility threshold"
    elseif plan_B.is_feasible && !plan_A.is_feasible
        best_plan = plan_B
        chosen_direction = "Direction 2 (subdivide B)"
        selection_reason = "Only Direction 2 meets feasibility threshold"
    elseif plan_A.is_feasible && plan_B.is_feasible
        # Both feasible - prefer the one with higher feasibility ratio
        ratio_A = total_A > 0 ? feasible_A / total_A : 0.0
        ratio_B = total_B > 0 ? feasible_B / total_B : 0.0
        
        if ratio_A > ratio_B
            best_plan = plan_A
            chosen_direction = "Direction 1 (subdivide A)"
            selection_reason = "Higher feasibility ratio ($(round(ratio_A*100, digits=1))% vs $(round(ratio_B*100, digits=1))%)"
        elseif ratio_B > ratio_A
            best_plan = plan_B
            chosen_direction = "Direction 2 (subdivide B)"
            selection_reason = "Higher feasibility ratio ($(round(ratio_B*100, digits=1))% vs $(round(ratio_A*100, digits=1))%)"
        else
            # Tied - prefer fewer total repairs (less invasive)
            if total_A <= total_B
                best_plan = plan_A
                chosen_direction = "Direction 1 (subdivide A)"
                selection_reason = "Equal feasibility ratio, fewer repairs needed ($total_A vs $total_B)"
            else
                best_plan = plan_B
                chosen_direction = "Direction 2 (subdivide B)"
                selection_reason = "Equal feasibility ratio, fewer repairs needed ($total_B vs $total_A)"
            end
        end
    else
        # Neither meets threshold - choose the one with more feasible plans
        if feasible_A > feasible_B
            best_plan = plan_A
            chosen_direction = "Direction 1 (subdivide A)"
            selection_reason = "More feasible plans ($feasible_A vs $feasible_B), though neither meets threshold"
        elseif feasible_B > feasible_A
            best_plan = plan_B
            chosen_direction = "Direction 2 (subdivide B)"
            selection_reason = "More feasible plans ($feasible_B vs $feasible_A), though neither meets threshold"
        else
            # Tied at 0 or equal - prefer fewer total repairs
            if total_A <= total_B
                best_plan = plan_A
                chosen_direction = "Direction 1 (subdivide A)"
                selection_reason = "Equal feasible plans, fewer total repairs ($total_A vs $total_B)"
            else
                best_plan = plan_B
                chosen_direction = "Direction 2 (subdivide B)"
                selection_reason = "Equal feasible plans, fewer total repairs ($total_B vs $total_A)"
            end
        end
    end
    
    direction_info["chosen_direction"] = chosen_direction
    direction_info["selection_reason"] = selection_reason
    
    if verbose
        println("\nComparison:")
        println("  Direction 1 (subdivide A): $feasible_A/$total_A feasible")
        println("  Direction 2 (subdivide B): $feasible_B/$total_B feasible")
        println("\n✓ Selected: $chosen_direction")
        println("  Reason: $selection_reason")
        println("="^70)
    end
    
    return (best_plan, true, direction_info)
end
