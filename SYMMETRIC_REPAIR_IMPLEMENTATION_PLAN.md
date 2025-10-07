# Symmetric Bidirectional Mesh Repair - Implementation Plan

## Executive Summary

This document outlines the implementation plan for transforming the mesh interface repair system from a **unidirectional source→target** approach to a **symmetric bidirectional** approach that generates a **unified third interface mesh** compatible with both sides.

### Current Architecture Limitations

1. **Early commitment**: Source/target sides are determined upfront by heuristics (edge density, constraint analysis)
2. **Asymmetric classification**: Each edge mismatch is classified from only one perspective
3. **One-sided repair**: Only the "target" side is modified to match the "source"
4. **Incomplete coverage**: Some mismatches may only be detectable from one direction

### Proposed Architecture Benefits

1. **Deferred decision**: No upfront source/target distinction
2. **Complete classification**: Every edge is classified from BOTH A→B and B→A perspectives
3. **Symmetric repair**: Generate a **third mesh** that satisfies constraints from both sides
4. **Bilateral replacement**: Replace BOTH original interface meshes with the unified result
5. **Local flexibility**: Each mismatch can independently choose the best repair strategy

---

## Phase 1: Current State Analysis

### What Already Exists (✓)

The codebase already has some bidirectional infrastructure:

1. **`classify_interface_mismatches_bidirectional()`** in `edge_classification.jl` (lines 1129-1346)
   - Classifies from both A→B and B→A perspectives
   - Creates swapped topology for reverse classification
   - Merges unique mismatches from both perspectives
   - Returns comparison metrics

2. **`generate_repair_plan_bidirectional()`** in `repair_planning.jl` (lines 1124-1267)
   - Generates repair plans in both directions
   - Selects the direction with more feasible repairs
   - **LIMITATION**: Still commits to ONE direction for repair

### What Needs to Change (✗)

1. **Classification results are merged but not both preserved**
   - Currently: merges to eliminate duplicates
   - Needed: preserve BOTH classifications for each edge for comparison

2. **Repair planning still picks one direction**
   - Currently: selects A→B OR B→A based on feasibility
   - Needed: use information from BOTH directions to build third mesh

3. **Repair execution modifies one side**
   - Currently: `apply_repair_plan!()` modifies only target_pid
   - Needed: create new third mesh and replace both sides

4. **Data structures assume unidirectional repair**
   - `EdgeMismatch`: stores single `present_in` and `should_be_in`
   - `RepairPlan`: has single `repair_direction`
   - Needed: new structures for symmetric repair

---

## Phase 2: New Data Structures

### 2.1 SymmetricEdgeMismatch

Stores classification results from BOTH perspectives for a single edge:

```julia
"""
Complete bidirectional classification for a single edge mismatch.
Stores the results of classifying the edge from both A→B and B→A perspectives.
"""
struct SymmetricEdgeMismatch
    edge_key::EdgeKey
    
    # Classification from A's perspective (A as source, B as target)
    # "This edge exists in A, is missing in B, how should B be modified?"
    classification_A_perspective::Union{EdgeMismatch, Nothing}
    
    # Classification from B's perspective (B as source, A as target)
    # "This edge exists in B, is missing in A, how should A be modified?"
    classification_B_perspective::Union{EdgeMismatch, Nothing}
    
    # Synthesis: which perspective(s) detected this edge?
    present_in_A::Bool
    present_in_B::Bool
    present_in_both::Bool  # True for shared edges
    
    # Synthesis: agreement analysis
    agree_on_type::Bool  # Do both perspectives assign same mismatch type?
    agree_on_feasibility::Bool  # Do both perspectives agree on repairability?
    
    # Recommended repair strategy (determined by local analysis)
    repair_strategy::Symbol  # :use_A | :use_B | :compromise | :skip
    repair_priority::Float64  # 0.0 (low) to 1.0 (high)
    
    # Conflict resolution notes
    resolution_reason::String
end
```

### 2.2 UnifiedInterfaceMesh

Represents the unified third mesh that will replace both original interfaces:

```julia
"""
Unified interface mesh generated from symmetric repair.
This mesh is compatible with BOTH sides and will replace both original interfaces.
"""
struct UnifiedInterfaceMesh
    # Triangulation
    triangles::Vector{Triangle}  # New unified triangulation
    
    # Node mapping to original meshes
    # Maps new node coords → (A_node_id, B_node_id)
    # Some nodes may be new (not in A or B originally)
    node_mapping_A::Dict{NTuple{3,Float64}, Union{Int, Nothing}}
    node_mapping_B::Dict{NTuple{3,Float64}, Union{Int, Nothing}}
    
    # Edge topology
    edges::Dict{EdgeKey, Vector{Int}}  # EdgeKey → triangle indices
    
    # Provenance tracking (for debugging)
    triangle_provenance::Vector{Symbol}  # :from_A | :from_B | :synthesized
    
    # Quality metrics
    min_triangle_quality::Float64
    total_area::Float64
    
    # Compatibility verification
    compatible_with_A::Bool
    compatible_with_B::Bool
    compatibility_report::Vector{String}
end
```

### 2.3 SymmetricRepairPlan

Replaces the unidirectional `RepairPlan`:

```julia
"""
Complete symmetric repair plan that generates a third unified mesh.
"""
struct SymmetricRepairPlan
    interface_pair::Tuple{Int,Int}  # (pidA, pidB)
    
    # Bidirectional classification results
    symmetric_mismatches::Vector{SymmetricEdgeMismatch}
    
    # The unified mesh to be generated
    target_unified_mesh::UnifiedInterfaceMesh
    
    # Operations to perform
    operations::Vector{UnifiedMeshOperation}  # Abstract operations
    
    # Feasibility
    is_feasible::Bool
    feasibility_issues::Vector{String}
    
    # Statistics
    edges_from_A::Int
    edges_from_B::Int
    edges_compromised::Int
    edges_synthesized::Int
    
    # Quality predictions
    predicted_min_quality::Float64
    predicted_compatibility_score::Float64
    
    # Original topology and constraints (for reference)
    topology::InterfaceTopology
    constraints::BoundaryConstraints
end
```

### 2.4 UnifiedMeshOperation

Abstract operation for building the unified mesh:

```julia
"""
Abstract operation for constructing the unified interface mesh.
Each operation describes a local triangulation decision.
"""
struct UnifiedMeshOperation
    operation_type::Symbol  # :copy_from_A | :copy_from_B | :retriangulate | :synthesize
    
    # Affected region (in terms of nodes)
    boundary_nodes::Vector{NTuple{3,Float64}}
    
    # Source triangles (if copying/adapting from A or B)
    source_triangles_A::Vector{Triangle}
    source_triangles_B::Vector{Triangle}
    
    # Result triangles for this operation
    result_triangles::Vector{NTuple{9,Float64}}  # Flattened triangle coords
    
    # Quality and feasibility
    min_quality::Float64
    is_feasible::Bool
    feasibility_notes::String
end
```

---

## Phase 3: Bidirectional Edge Classification Enhancement

### 3.1 Modify `classify_interface_mismatches_bidirectional()`

**Location**: `src/repair/edge_classification.jl`, lines 1129-1346

**Current behavior**: 
- Classifies from both perspectives
- Merges results to eliminate duplicates
- Returns merged classification

**New behavior**:
- Classify from both perspectives (keep this)
- **DO NOT merge** - preserve both classifications
- Create `SymmetricEdgeMismatch` for each unique edge
- Perform agreement analysis

**Implementation steps**:

```julia
function classify_interface_mismatches_symmetric(
    topology::InterfaceTopology;
    tol::Real=1e-4,
    verbose::Bool=true
)::Vector{SymmetricEdgeMismatch}
    
    # Step 1: Classify from A→B perspective
    classification_AB = classify_interface_mismatches(topology, tol=tol, debug=false)
    
    # Step 2: Create swapped topology for B→A perspective
    topology_swapped = create_swapped_topology(topology)
    
    # Step 3: Classify from B→A perspective
    classification_BA = classify_interface_mismatches(topology_swapped, tol=tol, debug=false)
    
    # Step 4: Build symmetric mismatch map
    # Key insight: Each unique edge needs BOTH classifications
    symmetric_map = Dict{EdgeKey, SymmetricEdgeMismatch}()
    
    # Process edges present in A
    for mismatch_AB in classification_AB.mismatches_B  # Edges in A, missing in B
        edge_key = mismatch_AB.edge_key
        
        # Find corresponding B→A classification (if exists)
        mismatch_BA = find_corresponding_mismatch(classification_BA.mismatches_A, edge_key)
        
        symmetric = SymmetricEdgeMismatch(
            edge_key,
            mismatch_AB,  # A's perspective
            mismatch_BA,  # B's perspective (may be nothing)
            true,   # present_in_A
            false,  # present_in_B
            false,  # present_in_both
            # ... compute agreement and strategy
        )
        
        symmetric_map[edge_key] = symmetric
    end
    
    # Process edges present in B (similar logic)
    for mismatch_BA in classification_BA.mismatches_B  # Edges in B, missing in A
        edge_key = mismatch_BA.edge_key
        if haskey(symmetric_map, edge_key)
            continue  # Already processed
        end
        
        mismatch_AB = find_corresponding_mismatch(classification_AB.mismatches_A, edge_key)
        
        symmetric = SymmetricEdgeMismatch(
            edge_key,
            mismatch_AB,  # A's perspective (may be nothing)
            mismatch_BA,  # B's perspective
            false,  # present_in_A
            true,   # present_in_B
            false,  # present_in_both
            # ... compute agreement and strategy
        )
        
        symmetric_map[edge_key] = symmetric
    end
    
    return collect(values(symmetric_map))
end
```

### 3.2 Agreement Analysis Logic

For each `SymmetricEdgeMismatch`, determine:

1. **Type agreement**: Do both perspectives assign the same `MismatchType`?
2. **Feasibility agreement**: Do both perspectives agree on `repair_feasible`?
3. **Quality agreement**: How do quality predictions compare?

This information guides the repair strategy selection.

---

## Phase 4: Local Repair Strategy Selection

### 4.1 Strategy Selection Algorithm

For each edge mismatch, select the best local repair approach:

```julia
"""
Determine the optimal repair strategy for a symmetric edge mismatch.

Decision tree:
1. If edge present in both (shared edge with different triangulation):
   - Compare quality predictions from both perspectives
   - Choose the triangulation with better quality
   - Strategy: :use_A or :use_B

2. If edge only in A:
   - Check B→A classification feasibility
   - If feasible, use A's triangulation: :use_A
   - If infeasible, check if A→B suggests alternative: :compromise or :skip

3. If edge only in B:
   - Check A→B classification feasibility
   - If feasible, use B's triangulation: :use_B
   - If infeasible, check if B→A suggests alternative: :compromise or :skip

4. Conflict resolution:
   - If both perspectives disagree on type/feasibility
   - Compute compromise solution or skip edge
"""
function determine_repair_strategy(
    sym_mismatch::SymmetricEdgeMismatch,
    constraints::BoundaryConstraints
)::Tuple{Symbol, Float64, String}
    
    # Unpack classifications
    m_A = sym_mismatch.classification_A_perspective
    m_B = sym_mismatch.classification_B_perspective
    
    # Case 1: Edge in both (shared edge with different triangulation)
    if sym_mismatch.present_in_both
        if m_A !== nothing && m_B !== nothing
            # Compare quality
            if m_A.min_affected_triangle_quality > m_B.min_affected_triangle_quality * 1.2
                return (:use_A, 0.8, "A has significantly better triangle quality")
            elseif m_B.min_affected_triangle_quality > m_A.min_affected_triangle_quality * 1.2
                return (:use_B, 0.8, "B has significantly better triangle quality")
            else
                # Quality similar - check feasibility
                if m_A.repair_feasible && !m_B.repair_feasible
                    return (:use_A, 0.6, "A is feasible, B is not")
                elseif m_B.repair_feasible && !m_A.repair_feasible
                    return (:use_B, 0.6, "B is feasible, A is not")
                else
                    return (:use_A, 0.5, "Tiebreaker: use A by default")
                end
            end
        end
    end
    
    # Case 2: Edge only in A
    if sym_mismatch.present_in_A && !sym_mismatch.present_in_B
        if m_A !== nothing && m_A.repair_feasible
            return (:use_A, 0.7, "Edge from A is feasible")
        else
            return (:skip, 0.0, "Edge from A is infeasible")
        end
    end
    
    # Case 3: Edge only in B
    if sym_mismatch.present_in_B && !sym_mismatch.present_in_A
        if m_B !== nothing && m_B.repair_feasible
            return (:use_B, 0.7, "Edge from B is feasible")
        else
            return (:skip, 0.0, "Edge from B is infeasible")
        end
    end
    
    # Fallback
    return (:skip, 0.0, "Cannot determine strategy")
end
```

### 4.2 Priority Assignment

Assign repair priority based on:
- Mismatch type (T_JUNCTION > DIAGONAL > REFINEMENT)
- Feasibility agreement
- Quality impact
- Constraint violations

High priority repairs are applied first.

---

## Phase 5: Third Mesh Generation Algorithm

### 5.1 Overview

Build the unified interface mesh by:
1. Starting with shared vertices (nodes present in both A and B)
2. For each edge mismatch, apply the selected repair strategy
3. Build a consistent triangulation covering the entire interface
4. Verify compatibility with both A and B

### 5.2 Algorithm Pseudocode

```julia
function generate_unified_interface_mesh(
    topology::InterfaceTopology,
    symmetric_mismatches::Vector{SymmetricEdgeMismatch},
    constraints::BoundaryConstraints
)::UnifiedInterfaceMesh
    
    # Step 1: Initialize with shared vertices
    unified_nodes = copy(topology.shared_node_keys)
    
    # Step 2: Sort mismatches by priority
    sorted_mismatches = sort(symmetric_mismatches, by=m->m.repair_priority, rev=true)
    
    # Step 3: Build region map (which triangles cover which boundary regions)
    # This helps detect overlaps and gaps
    region_map = initialize_region_map(topology)
    
    # Step 4: Process each mismatch according to strategy
    unified_triangles = Triangle[]
    operations = UnifiedMeshOperation[]
    
    for sym_mismatch in sorted_mismatches
        if sym_mismatch.repair_strategy == :skip
            continue
        end
        
        # Get the triangulation from selected source
        source_triangles = if sym_mismatch.repair_strategy == :use_A
            get_triangles_for_edge(topology.faces_A, sym_mismatch.edge_key)
        elseif sym_mismatch.repair_strategy == :use_B
            get_triangles_for_edge(topology.faces_B, sym_mismatch.edge_key)
        elseif sym_mismatch.repair_strategy == :compromise
            synthesize_compromise_triangulation(sym_mismatch, topology)
        end
        
        # Check for conflicts with existing unified triangulation
        if has_conflict(source_triangles, region_map)
            # Resolve conflict (may need to retriangulate)
            resolved_triangles = resolve_triangulation_conflict(
                source_triangles, region_map, unified_triangles
            )
            append!(unified_triangles, resolved_triangles)
        else
            # No conflict, add directly
            append!(unified_triangles, source_triangles)
            update_region_map!(region_map, source_triangles)
        end
        
        # Record operation
        push!(operations, create_operation_record(sym_mismatch, source_triangles))
    end
    
    # Step 5: Fill gaps (if any boundary regions remain uncovered)
    gaps = detect_gaps(region_map, topology)
    for gap in gaps
        gap_triangles = fill_gap(gap, unified_triangles, topology)
        append!(unified_triangles, gap_triangles)
    end
    
    # Step 6: Build node mappings
    node_mapping_A, node_mapping_B = build_node_mappings(
        unified_triangles, topology
    )
    
    # Step 7: Build edge topology
    unified_edges = build_edge_topology(unified_triangles)
    
    # Step 8: Compute quality metrics
    min_quality = compute_min_triangle_quality(unified_triangles)
    total_area = sum(tri.area for tri in unified_triangles)
    
    # Step 9: Verify compatibility
    compatible_A, issues_A = verify_compatibility(unified_triangles, topology.faces_A)
    compatible_B, issues_B = verify_compatibility(unified_triangles, topology.faces_B)
    
    return UnifiedInterfaceMesh(
        unified_triangles,
        node_mapping_A,
        node_mapping_B,
        unified_edges,
        determine_provenance(unified_triangles, topology),
        min_quality,
        total_area,
        compatible_A,
        compatible_B,
        vcat(issues_A, issues_B)
    )
end
```

### 5.3 Key Challenges

1. **Overlap detection**: Multiple repair strategies may produce overlapping triangles
2. **Gap filling**: Some regions may not be covered by any repair strategy
3. **Boundary consistency**: Unified mesh must respect both A's and B's boundaries
4. **Node merging**: Same geometric location may have different node IDs in A and B

---

## Phase 6: Mesh Replacement Operations

### 6.1 Replace Both Interface Meshes

After generating the unified mesh, replace BOTH original interface meshes:

```julia
"""
Replace both A's and B's interface faces with the unified mesh.

This is a transactional operation:
1. Begin transaction in RepairWorkspace
2. Delete all interface faces from PID A
3. Delete all interface faces from PID B
4. Insert unified mesh triangles into both PIDs (with appropriate node mappings)
5. Commit transaction if successful, rollback if failure
"""
function replace_both_interfaces!(
    workspace::RepairWorkspace,
    unified_mesh::UnifiedInterfaceMesh,
    topology::InterfaceTopology
)::Bool
    
    begin_transaction!(workspace)
    
    try
        # Step 1: Delete old interface faces from A
        for tri_idx in 1:length(topology.faces_A)
            delete_interface_face!(workspace, topology.pidA, tri_idx)
        end
        
        # Step 2: Delete old interface faces from B
        for tri_idx in 1:length(topology.faces_B)
            delete_interface_face!(workspace, topology.pidB, tri_idx)
        end
        
        # Step 3: Insert unified triangles into A
        for (tri_idx, tri) in enumerate(unified_mesh.triangles)
            # Map unified coordinates to A's node IDs
            node_ids_A = map_nodes_to_pid(tri, unified_mesh.node_mapping_A, topology.pidA)
            add_face!(workspace, topology.pidA, node_ids_A)
        end
        
        # Step 4: Insert unified triangles into B
        for (tri_idx, tri) in enumerate(unified_mesh.triangles)
            # Map unified coordinates to B's node IDs
            node_ids_B = map_nodes_to_pid(tri, unified_mesh.node_mapping_B, topology.pidB)
            add_face!(workspace, topology.pidB, node_ids_B)
        end
        
        # Step 5: Verify mesh integrity
        if !verify_mesh_integrity(workspace, topology.pidA, topology.pidB)
            @error "Mesh integrity check failed after replacement"
            rollback_transaction!(workspace)
            return false
        end
        
        # Success!
        commit_transaction!(workspace)
        return true
        
    catch e
        @error "Exception during interface replacement: $e"
        rollback_transaction!(workspace)
        return false
    end
end
```

### 6.2 Node Mapping Strategy

Critical issue: The same geometric point may have different node IDs in A and B.

**Solution**:
1. Create new node IDs in the workspace if needed
2. Use coordinate-based lookup to find existing node IDs
3. Merge nodes that are within tolerance
4. Update all references consistently

---

## Phase 7: Update Repair Execution Pipeline

### 7.1 New Execution Flow

**Old flow** (`repair_execution.jl`):
```
RepairPlan → determine target_pid → apply repairs to target_pid → commit
```

**New flow**:
```
SymmetricRepairPlan → generate UnifiedInterfaceMesh → 
replace interfaces in BOTH PIDs → verify → commit
```

### 7.2 Modified `repair_execution.jl`

Add new function:

```julia
"""
Execute symmetric repair plan by generating and installing unified interface mesh.
"""
function execute_symmetric_repair!(
    workspace::RepairWorkspace,
    plan::SymmetricRepairPlan
)::Bool
    
    println("\\n" * "="^70)
    println("Executing Symmetric Repair Plan")
    println("="^70)
    println("Interface: $(plan.interface_pair)")
    println("Edges from A: $(plan.edges_from_A)")
    println("Edges from B: $(plan.edges_from_B)")
    println("Compromised edges: $(plan.edges_compromised)")
    println("Predicted min quality: $(round(plan.predicted_min_quality, digits=3))")
    println("="^70)
    
    if !plan.is_feasible
        @error "Plan is not feasible!"
        println("Issues:")
        for issue in plan.feasibility_issues
            println("  - $issue")
        end
        return false
    end
    
    # Generate the unified mesh
    println("\\nGenerating unified interface mesh...")
    unified_mesh = plan.target_unified_mesh
    
    println("Unified mesh statistics:")
    println("  Triangles: $(length(unified_mesh.triangles))")
    println("  Min quality: $(round(unified_mesh.min_triangle_quality, digits=3))")
    println("  Compatible with A: $(unified_mesh.compatible_with_A)")
    println("  Compatible with B: $(unified_mesh.compatible_with_B)")
    
    # Replace both interfaces
    println("\\nReplacing both interface meshes...")
    success = replace_both_interfaces!(
        workspace,
        unified_mesh,
        plan.topology
    )
    
    if success
        println("\\n✓ Symmetric repair completed successfully!")
    else
        println("\\n✗ Symmetric repair failed!")
    end
    
    return success
end
```

---

## Phase 8: Testing Strategy

### 8.1 Unit Tests

Create comprehensive tests for each component:

1. **Symmetric classification tests**:
   - Edge present in A only
   - Edge present in B only
   - Edge present in both (different triangulation)
   - Agreement/disagreement cases

2. **Strategy selection tests**:
   - Quality-based selection
   - Feasibility-based selection
   - Conflict resolution
   - Priority ordering

3. **Mesh generation tests**:
   - Simple case (no conflicts)
   - Overlapping triangulations
   - Gap filling
   - Boundary preservation

4. **Replacement tests**:
   - Node mapping correctness
   - Topology preservation
   - Transaction rollback on failure

### 8.2 Integration Tests

End-to-end tests with realistic meshes:

1. **Simple T-junction case**: A has edge, B has hanging node
2. **Diagonal mismatch case**: Same quad, different diagonal
3. **Mixed case**: Multiple mismatch types
4. **Constraint-heavy case**: Many locked edges
5. **Large-scale case**: Thousands of mismatches

### 8.3 Comparison with Old Approach

For each test case:
- Run old unidirectional repair
- Run new symmetric repair
- Compare results:
  - Number of edges matched
  - Quality metrics
  - Conformity ratio
  - Processing time

---

## Phase 9: Migration Path

### 9.1 Backward Compatibility

Keep the old unidirectional approach available:

```julia
# Old API (preserved)
generate_repair_plan(topology, classification, constraints)

# New API (symmetric)
generate_symmetric_repair_plan(topology, symmetric_classification, constraints)
```

### 9.2 Gradual Rollout

1. **Phase 9.1**: Implement symmetric classification only
2. **Phase 9.2**: Implement strategy selection and mesh generation
3. **Phase 9.3**: Implement replacement operations
4. **Phase 9.4**: Test extensively with both approaches
5. **Phase 9.5**: Switch default to symmetric approach
6. **Phase 9.6**: Deprecate old approach after validation period

---

## Implementation Priority and Effort Estimates

| Phase | Component | Complexity | Est. Time | Priority |
|-------|-----------|------------|-----------|----------|
| 1 | Analysis | Low | 1 day | P0 |
| 2 | Data structures | Medium | 2 days | P0 |
| 3 | Symmetric classification | Medium | 3 days | P0 |
| 4 | Strategy selection | High | 4 days | P1 |
| 5 | Mesh generation | **Very High** | 7 days | P1 |
| 6 | Mesh replacement | High | 3 days | P1 |
| 7 | Execution pipeline | Medium | 2 days | P2 |
| 8 | Testing | High | 5 days | P2 |
| 9 | Migration | Low | 2 days | P3 |

**Total estimated effort**: ~29 days

**Critical path**: Phases 2 → 3 → 5 → 6

**Riskiest component**: Phase 5 (Mesh generation) due to complexity of:
- Overlap detection and resolution
- Gap filling algorithm
- Maintaining boundary consistency
- Ensuring manifoldness

---

## Key Technical Challenges

### Challenge 1: Triangulation Conflicts

**Problem**: When using edges from both A and B, their triangulations may overlap or conflict.

**Example**:
```
A's triangulation: [(v1,v2,v3), (v2,v3,v4)]
B's triangulation: [(v1,v2,v5), (v2,v4,v5)]  # Same region, different split
```

**Solution approaches**:
1. **Priority-based**: Higher-priority mismatches overwrite lower-priority ones
2. **Quality-based**: Keep triangulation with better quality
3. **Retriangulation**: Detect conflict, retriangulate the conflicting region
4. **Hybrid**: Use constraints from both to find compromise triangulation

### Challenge 2: Gap Filling

**Problem**: After applying all repair strategies, some boundary regions may remain uncovered.

**Solution**:
1. Detect boundary gaps using edge traversal
2. Identify boundary loops
3. Use Delaunay triangulation or advancing front method to fill gaps
4. Ensure quality thresholds are met

### Challenge 3: Node Identity

**Problem**: Same geometric location has different node IDs in A and B.

**Solution**:
- Use coordinate-based mapping with tolerance
- Create unified node ID space in the workspace
- Maintain bidirectional mapping: coord ↔ A_node_id ↔ B_node_id

### Challenge 4: Boundary Constraints

**Problem**: Unified mesh must respect locked edges from BOTH sides.

**Solution**:
- Merge constraint sets from A and B
- If A and B have conflicting constraints, flag as irreconcilable
- Prioritize stronger constraints (e.g., feature edges > regular boundaries)

---

## Success Criteria

The implementation will be considered successful if:

1. **Completeness**: System can handle all mismatch types from both perspectives
2. **Quality**: Unified mesh maintains triangle quality ≥ input meshes
3. **Conformity**: Achieves higher conformity ratio than unidirectional approach
4. **Compatibility**: Unified mesh integrates seamlessly with both A and B
5. **Performance**: Execution time ≤ 2× unidirectional approach
6. **Robustness**: Handles conflicts and gaps gracefully

---

## Next Steps

1. **Review and approve** this implementation plan
2. **Set up feature branch** for symmetric repair development
3. **Begin Phase 1**: Analyze current bidirectional infrastructure
4. **Prototype Phase 2**: Design and implement new data structures
5. **Implement Phase 3**: Enhance classification to preserve both perspectives
6. **Proceed incrementally** through remaining phases with testing at each step

---

## Appendix A: Key Files to Modify

| File | Current Lines | Modifications Needed |
|------|---------------|---------------------|
| `edge_classification.jl` | 1427 | Add symmetric classification functions (~300 lines) |
| `repair_planning.jl` | 1267 | Add symmetric planning functions (~500 lines) |
| `repair_execution.jl` | 331 | Add symmetric execution (~200 lines) |
| `interface_topology.jl` | 554 | Minor additions for mesh merging (~50 lines) |
| `repair_workspace.jl` | (to check) | Add node remapping functions (~100 lines) |

**Estimated total new code**: ~1150 lines

---

## Appendix B: Terminology Reference

- **Source mesh**: The mesh containing an edge
- **Target mesh**: The mesh missing the edge (where repair is needed)
- **Symmetric classification**: Classifying each edge from BOTH A→B and B→A perspectives
- **Unified mesh**: The third interface mesh compatible with both A and B
- **Bidirectional repair**: Repair approach that considers both perspectives equally
- **Local repair strategy**: Per-edge decision on how to handle the mismatch
- **Third mesh**: The result of symmetric repair that replaces both original interfaces

---

## Appendix C: Advantages Over Current Approach

| Aspect | Current (Unidirectional) | Proposed (Symmetric) |
|--------|-------------------------|---------------------|
| Source/target decision | Upfront heuristic | Deferred per-edge |
| Classification coverage | One perspective | Both perspectives |
| Repair flexibility | Fixed direction | Local flexibility |
| Result placement | Replaces one side | Replaces both sides |
| Quality optimization | One-sided | Global optimum |
| Constraint satisfaction | Target side only | Both sides equally |
| Mesh compatibility | One-way | Bidirectional |

The symmetric approach provides **comprehensive coverage**, **local flexibility**, and **true bidirectional compatibility**.
