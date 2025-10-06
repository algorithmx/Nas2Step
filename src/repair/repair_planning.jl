# Repair strategy generation for interface conformity
# Phase 2: Plan surgical mesh repairs

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

"""
    assess_interface_feasibility(plans, total_needed, thresholds; ratio=0.8)

Compute per-interface feasibility from a set of EdgeInsertionPlans.
Returns a named tuple with:
- feasible_plans / infeasible_plans
- is_feasible: overall decision (treats total_needed==0 as feasible)
- issues: list of human-readable issues
- reason_counts: constraint violation reasons and counts
- qstats: light quality stats (min-before/after, normalized threshold, delta stats)
- counts: rich counters (coverage, category breakdown, required minimum)
- ratios: coverage and feasible ratios vs. needs
- split_types: count of plans by split_type
- top_reasons: sorted top-3 constraint reasons
"""
function assess_interface_feasibility(
    plans::Vector{EdgeInsertionPlan},
    total_needed::Int,
    thresholds::QualityThresholds;
    ratio::Float64 = PLAN_FEASIBLE_RATIO_DEFAULT,
)
    # Helper for simple quantiles without extra deps
    function _quantiles(x::Vector{Float64})
        if isempty(x)
            return (p25=0.0, p50=0.0, p75=0.0)
        end
        s = sort(x)
        n = length(s)
        idx25 = max(1, min(n, Int(floor(0.25*(n-1))) + 1))
        idx50a = ((n + 1) >>> 1)  # floor((n+1)/2)
        idx50b = isodd(n) ? idx50a : idx50a + (idx50a < n ? 1 : 0)
        med = isodd(n) ? s[idx50a] : 0.5*(s[idx50a] + s[idx50b])
        idx75 = max(1, min(n, Int(floor(0.75*(n-1))) + 1))
        return (p25=s[idx25], p50=med, p75=s[idx75])
    end

    feasible_plans = filter(p -> p.quality_acceptable && !p.violates_constraints, plans)
    infeasible_plans = filter(p -> !p.quality_acceptable || p.violates_constraints, plans)

    n_plans = length(plans)
    n_feasible = length(feasible_plans)
    n_infeasible = length(infeasible_plans)

    # Category breakdown
    n_only_quality = count(p -> (!p.quality_acceptable) && !p.violates_constraints, plans)
    n_only_constraints = count(p -> p.quality_acceptable && p.violates_constraints, plans)
    n_both = count(p -> (!p.quality_acceptable) && p.violates_constraints, plans)
    n_constraints = n_only_constraints + n_both
    n_poor_quality = n_only_quality + n_both

    # Overall feasibility gate
    # Treat zero-needed as feasible (nothing to do)
    required_threshold = total_needed * ratio
    required_min = total_needed == 0 ? 0 : (floor(Int, required_threshold) + 1)  # due to strict >
    is_feasible = total_needed == 0 ? true : (n_feasible > required_threshold)

    issues = String[]
    if n_plans < total_needed
        push!(issues, "Coverage gap: generated $(n_plans)/$(total_needed) plans (short by $(total_needed - n_plans))")
    end
    if !is_feasible
        # New detailed phrasing
        push!(issues, "Feasible plans $(n_feasible)/$(total_needed) below required > $(round(required_threshold, digits=2)) (need at least $(required_min))")
        # Legacy phrasing (compat)
        push!(issues, "Only $(n_feasible)/$(total_needed) plans are feasible")
    end
    if n_constraints > 0
        push!(issues, "$(n_constraints) plans blocked by constraints ($(n_only_constraints) only-constraints, $(n_both) both)")
        # Legacy wording for compatibility
        push!(issues, "$(n_constraints) plans violate constraints")
    end
    if n_poor_quality > 0
        push!(issues, "$(n_poor_quality) plans fail quality ($(n_only_quality) only-quality, $(n_both) both)")
        # Legacy wording for compatibility
        push!(issues, "$(n_poor_quality) plans would create poor quality triangles")
    end

    # Reason breakdown for constraints
    reason_counts = Dict{String,Int}()
    for p in plans
        if p.violates_constraints
            for r in p.constraint_violations
                reason_counts[r] = get(reason_counts, r, 0) + 1
            end
        end
    end
    # Top 3 reasons
    top_reasons = sort(collect(reason_counts); by = x -> -x[2])[1:min(3, length(reason_counts))]

    # Quality stats
    q_before = [p.min_angle_before for p in plans]
    q_after  = [p.min_angle_after for p in plans]
    deltas = [p.min_angle_after - p.min_angle_before for p in plans]
    qstats = (
        min_before = isempty(q_before) ? 0.0 : minimum(q_before),
        min_after = isempty(q_after) ? 0.0 : minimum(q_after),
        min_threshold_norm = thresholds.min_angle / 60.0,
        p25_after = _quantiles(q_after).p25,
        median_after = _quantiles(q_after).p50,
        p75_after = _quantiles(q_after).p75,
        delta_min = isempty(deltas) ? 0.0 : minimum(deltas),
        delta_max = isempty(deltas) ? 0.0 : maximum(deltas),
        delta_med = _quantiles(deltas).p50,
    )

    # Split type counts
    split_types = Dict{Symbol,Int}()
    for p in plans
        split_types[p.split_type] = get(split_types, p.split_type, 0) + 1
    end

    return (
        feasible_plans = feasible_plans,
        infeasible_plans = infeasible_plans,
        is_feasible = is_feasible,
        issues = issues,
        reason_counts = reason_counts,
        qstats = qstats,
        counts = (
            total_needed = total_needed,
            generated = n_plans,
            feasible = n_feasible,
            infeasible = n_infeasible,
            only_quality = n_only_quality,
            only_constraints = n_only_constraints,
            both = n_both,
            constraints = n_constraints,
            poor_quality = n_poor_quality,
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
    corner1 = diagonal_edge.node1
    corner2 = diagonal_edge.node2
    
    # Find the other 2 vertices (not on the diagonal)
    other_vertices = filter(v -> v != corner1 && v != corner2, quad_vertices)
    
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
    tri_keys = Set([
        (round(triangle.coord1[1], digits=4), round(triangle.coord1[2], digits=4), round(triangle.coord1[3], digits=4)),
        (round(triangle.coord2[1], digits=4), round(triangle.coord2[2], digits=4), round(triangle.coord2[3], digits=4)),
        (round(triangle.coord3[1], digits=4), round(triangle.coord3[2], digits=4), round(triangle.coord3[3], digits=4))
    ])
    
    has_node1 = edge_to_insert.node1 in tri_keys
    has_node2 = edge_to_insert.node2 in tri_keys
    
    if !has_node1 || !has_node2
        # Edge doesn't fully belong to this triangle
        return nothing
    end
    
    # Find which vertex is NOT on the edge
    opposite_vertex = nothing
    if edge_to_insert.node1 ∉ [tri_keys...] || edge_to_insert.node2 ∉ [tri_keys...]
        # This shouldn't happen if we filtered correctly
        return nothing
    end
    
    tri_coords = [triangle.coord1, triangle.coord2, triangle.coord3]
    edge_coords = Set([edge_to_insert.node1, edge_to_insert.node2])
    
    opposite_idx = findfirst(c -> begin
        ck = (round(c[1], digits=4), round(c[2], digits=4), round(c[3], digits=4))
        ck ∉ edge_coords
    end, tri_coords)
    
    if opposite_idx === nothing
        return nothing
    end
    
    opposite_vertex = tri_coords[opposite_idx]
    
    # Create two new triangles by bisecting
    # Triangle 1: (edge.node1, edge.node2, opposite)
    # This is essentially just reconfirming the edge exists
    
    new_triangles = [
        (edge_to_insert.node1[1], edge_to_insert.node1[2], edge_to_insert.node1[3],
         edge_to_insert.node2[1], edge_to_insert.node2[2], edge_to_insert.node2[3],
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
"""
function generate_edge_insertion_plan(
    mismatch::EdgeMismatch,
    topology::InterfaceTopology,
    constraints::BoundaryConstraints,
    thresholds::QualityThresholds
)::Union{EdgeInsertionPlan, Nothing}
    
    # Determine target side (where edge needs to be inserted)
    target_faces = mismatch.should_be_in == topology.pidA ? topology.faces_A : topology.faces_B
    
    # Check constraint violations
    has_violation, violations = check_constraint_violations(mismatch.edge_key, constraints)
    
    # Handle DIAGONAL mismatches differently
    if mismatch.mismatch_type == DIAGONAL
        # For diagonal mismatches, use quad retriangulation
        if isempty(mismatch.quad_vertices) || isempty(mismatch.triangles_to_replace)
            return EdgeInsertionPlan(
                1,  # Placeholder
                mismatch.edge_key,
                :failed,
                NTuple{3,Float64}[],
                Int[],
                NTuple{9,Float64}[],  # old_triangles
                NTuple{9,Float64}[],  # replacement_triangles
                0.0,
                0.0,
                false,
                Int[],
                true,
                ["Quad vertices or triangles to replace not found"],
                false  # is_feasible
            )
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
            return EdgeInsertionPlan(
                mismatch.triangles_to_replace[1],
                mismatch.edge_key,
                :failed,
                NTuple{3,Float64}[],
                Int[],
                NTuple{9,Float64}[],  # old_triangles
                NTuple{9,Float64}[],  # replacement_triangles
                min_angle_before,
                0.0,
                false,
                Int[],
                true,
                ["Cannot plan quad retriangulation"],
                false  # is_feasible
            )
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
    
    # Handle T-junction and other types
    if isempty(mismatch.affected_triangles)
        return nothing
    end
    
    # Take first affected triangle
    target_tri_idx = mismatch.affected_triangles[1]
    if target_tri_idx < 1 || target_tri_idx > length(target_faces)
        return nothing
    end
    
    target_triangle = target_faces[target_tri_idx]
    
    # Compute current quality
    min_angle_before = compute_triangle_quality(target_triangle)
    
    # Plan the split
    split_plan = plan_triangle_split(target_triangle, mismatch.edge_key)
    
    if split_plan === nothing
        return EdgeInsertionPlan(
            target_tri_idx,
            mismatch.edge_key,
            :failed,
            NTuple{3,Float64}[],
            Int[],
            NTuple{9,Float64}[],  # old_triangles
            NTuple{9,Float64}[],  # replacement_triangles
            min_angle_before,
            0.0,
            false,
            Int[],
            true,
            ["Cannot plan split for this triangle"],
            false  # is_feasible
        )
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
    failed_count = 0
    for (idx, mismatch) in enumerate(mismatches_to_fix)
        plan = generate_edge_insertion_plan(mismatch, topology, constraints, thresholds)
        if plan !== nothing
            push!(insertion_plans, plan)
        else
            failed_count += 1
        end
        
        if idx % 50 == 0
            println("    Progress: $idx/$(length(mismatches_to_fix))")
        end
    end
    
    println("  Generated $(length(insertion_plans)) insertion plans ($failed_count failed)")
    
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

    # Contextual prints
    println("\n  Validation thresholds:")
    println("    min_angle ≥ $(round(thresholds.min_angle, digits=2))°, max_angle ≤ $(round(thresholds.max_angle, digits=2))°, max_aspect ≤ $(round(thresholds.max_aspect_ratio, digits=2))")
    if !isempty(insertion_plans)
        println("  Quality (normalized, 1.0 = 60°):")
        println("    min before: $(round(assessment.qstats.min_before, digits=3)), min after: $(round(assessment.qstats.min_after, digits=3))")
    end
    if !isempty(assessment.reason_counts)
        println("  Constraint violations breakdown:")
        for (reason, cnt) in sort(collect(assessment.reason_counts); by = x -> -x[2])
            println("    $(reason): $(cnt)")
        end
    end
    
    println("\n  Feasibility assessment:")
    println("    Generated: $(assessment.counts.generated) / Needed: $(assessment.counts.total_needed)  (coverage=$(round(assessment.ratios.coverage*100, digits=1))%)")
    println("    Feasible: $(assessment.counts.feasible) / Needed: $(assessment.counts.total_needed)  (feasible ratio=$(round(assessment.ratios.feasible*100, digits=1))%, required > $(round(assessment.counts.total_needed * assessment.ratios.required, digits=2)))")
    println("    Infeasible: $(assessment.counts.infeasible)  [only-quality=$(assessment.counts.only_quality), only-constraints=$(assessment.counts.only_constraints), both=$(assessment.counts.both)]")
    println("    Overall feasible: $(assessment.is_feasible)")
    
    if !isempty(assessment.issues)
        println("    Issues:")
        for issue in assessment.issues
            println("      - $issue")
        end
    end
    if !isempty(assessment.top_reasons)
        println("    Top constraint reasons:")
        for (reason, cnt) in assessment.top_reasons
            println("      - $(reason): $(cnt)")
        end
    end
    if !isempty(insertion_plans)
        println("    Split types:")
        for (stype, cnt) in sort(collect(assessment.split_types); by=x->-x[2])
            println("      - $(stype): $(cnt)")
        end
        println("    Quality after percentiles (normalized): p25=$(round(assessment.qstats.p25_after, digits=3)), p50=$(round(assessment.qstats.median_after, digits=3)), p75=$(round(assessment.qstats.p75_after, digits=3))")
        println("    Quality delta (after-before): min=$(round(assessment.qstats.delta_min, digits=3)), med=$(round(assessment.qstats.delta_med, digits=3)), max=$(round(assessment.qstats.delta_max, digits=3))")
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
