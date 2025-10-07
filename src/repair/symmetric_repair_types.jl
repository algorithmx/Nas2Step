# Symmetric Repair Data Structures
# Phase 2: Core types for bidirectional symmetric mesh repair
#
# This module defines the data structures needed for symmetric mesh repair,
# where a unified third interface mesh is generated that is compatible with
# both sides (A and B) and replaces both original interfaces.

# ============================================================================
# SymmetricEdgeMismatch - Bidirectional edge classification
# ============================================================================

"""
    SymmetricEdgeMismatch

Complete bidirectional classification for a single edge mismatch.

Stores the results of classifying the edge from both A→B and B→A perspectives,
enabling local per-edge decisions about repair strategy without upfront
commitment to a source/target direction.

# Fields
- `edge_key::EdgeKey`: The edge being classified
- `classification_A_perspective::Union{EdgeMismatch, Nothing}`: Classification when A is source, B is target
- `classification_B_perspective::Union{EdgeMismatch, Nothing}`: Classification when B is source, A is target
- `present_in_A::Bool`: Edge exists in mesh A
- `present_in_B::Bool`: Edge exists in mesh B
- `present_in_both::Bool`: Edge exists in both (with different triangulation)
- `agree_on_type::Bool`: Both perspectives assign same MismatchType
- `agree_on_feasibility::Bool`: Both perspectives agree on repairability
- `repair_strategy::Symbol`: Selected strategy (:use_A | :use_B | :compromise | :skip)
- `repair_priority::Float64`: Priority 0.0 (low) to 1.0 (high)
- `resolution_reason::String`: Explanation of strategy selection

# Classification Perspective Meanings
- **A perspective** (A→B): "This edge exists in A, is missing in B. How should B be modified?"
- **B perspective** (B→A): "This edge exists in B, is missing in A. How should A be modified?"

# Strategy Semantics
- `:use_A`: Use A's triangulation for this edge in unified mesh
- `:use_B`: Use B's triangulation for this edge in unified mesh
- `:compromise`: Synthesize a compromise triangulation from both
- `:skip`: Cannot repair this edge (mark for manual review)

# Example
```julia
# Edge present only in A, feasible to add to B
sym = SymmetricEdgeMismatch(
    edge_key,
    mismatch_from_A,  # A's perspective
    nothing,          # B doesn't have this edge
    true,   # in A
    false,  # not in B
    false,  # not in both
    false,  # can't agree (only one perspective)
    true,   # A perspective says feasible
    :use_A, # Use A's triangulation
    0.7,    # Medium-high priority
    "Edge from A is feasible and has good quality"
)
```
"""
struct SymmetricEdgeMismatch
    edge_key::EdgeKey
    
    # Classification from A's perspective (A as source, B as target)
    # "This edge exists in A, is missing in B, how should B be modified?"
    classification_A_perspective::Union{EdgeMismatch, Nothing}
    
    # Classification from B's perspective (B as source, A as target)
    # "This edge exists in B, is missing in A, how should A be modified?"
    classification_B_perspective::Union{EdgeMismatch, Nothing}
    
    # Presence analysis
    present_in_A::Bool
    present_in_B::Bool
    present_in_both::Bool  # True for shared edges with different triangulation
    
    # Agreement analysis
    agree_on_type::Bool          # Do both perspectives assign same mismatch type?
    agree_on_feasibility::Bool   # Do both perspectives agree on repairability?
    
    # Repair decision (computed by strategy selection algorithm)
    repair_strategy::Symbol      # :use_A | :use_B | :compromise | :skip
    repair_priority::Float64     # 0.0 (low) to 1.0 (high)
    
    # Explanation
    resolution_reason::String
end

"""
    SymmetricEdgeMismatch(edge_key, class_A, class_B, strategy, priority, reason)

Convenience constructor that automatically computes presence and agreement fields.
"""
function SymmetricEdgeMismatch(
    edge_key::EdgeKey,
    classification_A::Union{EdgeMismatch, Nothing},
    classification_B::Union{EdgeMismatch, Nothing},
    repair_strategy::Symbol,
    repair_priority::Float64,
    resolution_reason::String
)
    # Determine presence
    present_in_A = classification_A !== nothing
    present_in_B = classification_B !== nothing
    present_in_both = present_in_A && present_in_B
    
    # Compute agreement
    agree_on_type = false
    agree_on_feasibility = false
    
    if classification_A !== nothing && classification_B !== nothing
        # Both classifications exist - can compare
        agree_on_type = classification_A.mismatch_type == classification_B.mismatch_type
        agree_on_feasibility = classification_A.repair_feasible == classification_B.repair_feasible
    end
    
    return SymmetricEdgeMismatch(
        edge_key,
        classification_A,
        classification_B,
        present_in_A,
        present_in_B,
        present_in_both,
        agree_on_type,
        agree_on_feasibility,
        repair_strategy,
        repair_priority,
        resolution_reason
    )
end

# ============================================================================
# UnifiedInterfaceMesh - The third mesh that replaces both interfaces
# ============================================================================

"""
    UnifiedInterfaceMesh

Unified interface mesh generated from symmetric repair.

This mesh is compatible with BOTH sides and will replace both original interface
meshes. It represents the "third mesh" that satisfies constraints from both A and B.

# Fields
- `triangles::Vector{Triangle}`: The unified triangulation
- `node_mapping_A::Dict{NTuple{3,Float64}, Union{Int, Nothing}}`: Maps coords → A's node IDs
- `node_mapping_B::Dict{NTuple{3,Float64}, Union{Int, Nothing}}`: Maps coords → B's node IDs
- `edges::Dict{EdgeKey, Vector{Int}}`: Edge topology (EdgeKey → triangle indices)
- `triangle_provenance::Vector{Symbol}`: Origin of each triangle (:from_A | :from_B | :synthesized)
- `min_triangle_quality::Float64`: Worst triangle quality in mesh
- `total_area::Float64`: Total surface area
- `compatible_with_A::Bool`: Can replace A's interface without issues
- `compatible_with_B::Bool`: Can replace B's interface without issues
- `compatibility_report::Vector{String}`: Issues found during compatibility check

# Node Mapping
The node mappings handle the fact that the same geometric location may have
different node IDs in meshes A and B. For each coordinate in the unified mesh:
- `node_mapping_A[coord]`: The node ID in mesh A (or `nothing` if new node)
- `node_mapping_B[coord]`: The node ID in mesh B (or `nothing` if new node)

# Triangle Provenance
Tracks where each triangle came from:
- `:from_A`: Copied/adapted from mesh A's triangulation
- `:from_B`: Copied/adapted from mesh B's triangulation
- `:synthesized`: Created by conflict resolution or gap filling

# Compatibility
The mesh is considered compatible if:
1. All boundary edges from original meshes are preserved
2. All shared vertices are included
3. No topological violations (manifoldness, etc.)
4. Quality thresholds are met

# Example
```julia
unified = UnifiedInterfaceMesh(
    triangles,           # 150 triangles
    node_mapping_A,      # Maps to A's node IDs
    node_mapping_B,      # Maps to B's node IDs
    edges,               # Edge topology
    provenance,          # [:from_A, :from_A, :from_B, ...]
    0.35,                # min quality (acceptable)
    42.7,                # total area
    true,                # compatible with A
    true,                # compatible with B
    String[]             # no compatibility issues
)
```
"""
struct UnifiedInterfaceMesh
    # Triangulation
    triangles::Vector{Triangle}
    
    # Node mapping to original meshes
    # Maps unified node coordinates → original node IDs
    # Nothing indicates a new node not present in original mesh
    node_mapping_A::Dict{NTuple{3,Float64}, Union{Int, Nothing}}
    node_mapping_B::Dict{NTuple{3,Float64}, Union{Int, Nothing}}
    
    # Edge topology
    edges::Dict{EdgeKey, Vector{Int}}  # EdgeKey → triangle indices in unified mesh
    
    # Provenance tracking (for debugging and analysis)
    triangle_provenance::Vector{Symbol}  # :from_A | :from_B | :synthesized
    
    # Quality metrics
    min_triangle_quality::Float64
    total_area::Float64
    
    # Compatibility verification
    compatible_with_A::Bool
    compatible_with_B::Bool
    compatibility_report::Vector{String}
end

"""
    UnifiedInterfaceMesh()

Create an empty UnifiedInterfaceMesh (useful for incremental construction).
"""
function UnifiedInterfaceMesh()
    return UnifiedInterfaceMesh(
        Triangle[],
        Dict{NTuple{3,Float64}, Union{Int, Nothing}}(),
        Dict{NTuple{3,Float64}, Union{Int, Nothing}}(),
        Dict{EdgeKey, Vector{Int}}(),
        Symbol[],
        1.0,  # Perfect quality initially
        0.0,  # Zero area initially
        false,  # Not yet verified
        false,  # Not yet verified
        String[]
    )
end

# ============================================================================
# UnifiedMeshOperation - Atomic operation for building unified mesh
# ============================================================================

"""
    UnifiedMeshOperation

Atomic operation for constructing the unified interface mesh.

Each operation describes a local triangulation decision for a region of the
interface. Operations are applied in priority order to build the unified mesh.

# Fields
- `operation_type::Symbol`: Type of operation (:copy_from_A | :copy_from_B | :retriangulate | :synthesize)
- `boundary_nodes::Vector{NTuple{3,Float64}}`: Nodes defining the affected region
- `source_triangles_A::Vector{Triangle}`: Triangles from A covering this region
- `source_triangles_B::Vector{Triangle}`: Triangles from B covering this region
- `result_triangles::Vector{NTuple{9,Float64}}`: Triangles to add to unified mesh (flattened coords)
- `min_quality::Float64`: Worst triangle quality in result
- `is_feasible::Bool`: Whether operation can be performed
- `feasibility_notes::String`: Explanation of feasibility status

# Operation Types
- `:copy_from_A`: Copy triangulation directly from mesh A
- `:copy_from_B`: Copy triangulation directly from mesh B
- `:retriangulate`: Retriangulate a quad with different diagonal
- `:synthesize`: Create new triangulation (gap filling or conflict resolution)

# Triangle Format
Result triangles are stored as flattened NTuple{9,Float64}:
```
(v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z)
```

This format is used for serialization and matches the repair plan format.

# Example
```julia
op = UnifiedMeshOperation(
    :copy_from_A,                    # Use A's triangulation
    boundary_nodes,                  # Nodes forming region boundary
    [tri1_A, tri2_A],               # A's triangles
    Triangle[],                      # B has no triangles here
    [(0,0,0, 1,0,0, 0,1,0), ...],  # Result triangles
    0.42,                            # min quality
    true,                            # feasible
    "Using A's triangulation with good quality"
)
```
"""
struct UnifiedMeshOperation
    operation_type::Symbol  # :copy_from_A | :copy_from_B | :retriangulate | :synthesize
    
    # Affected region (in terms of boundary nodes)
    boundary_nodes::Vector{NTuple{3,Float64}}
    
    # Source triangles (if copying/adapting from A or B)
    source_triangles_A::Vector{Triangle}
    source_triangles_B::Vector{Triangle}
    
    # Result triangles for this operation (flattened coordinates)
    result_triangles::Vector{NTuple{9,Float64}}
    
    # Quality assessment
    min_quality::Float64
    is_feasible::Bool
    feasibility_notes::String
end

"""
    UnifiedMeshOperation(type, result_tris, quality, feasible, notes)

Convenience constructor for operations without source triangle references.
"""
function UnifiedMeshOperation(
    operation_type::Symbol,
    result_triangles::Vector{NTuple{9,Float64}},
    min_quality::Float64,
    is_feasible::Bool,
    feasibility_notes::String
)
    return UnifiedMeshOperation(
        operation_type,
        NTuple{3,Float64}[],  # No boundary nodes specified
        Triangle[],            # No source A triangles
        Triangle[],            # No source B triangles
        result_triangles,
        min_quality,
        is_feasible,
        feasibility_notes
    )
end

# ============================================================================
# SymmetricRepairPlan - Complete symmetric repair plan
# ============================================================================

"""
    SymmetricRepairPlan

Complete symmetric repair plan that generates a unified third interface mesh.

Unlike the unidirectional `RepairPlan` which modifies only one side, this plan
generates a unified mesh that replaces BOTH original interface meshes.

# Fields
- `interface_pair::Tuple{Int,Int}`: (pidA, pidB)
- `symmetric_mismatches::Vector{SymmetricEdgeMismatch}`: All classified edges
- `target_unified_mesh::UnifiedInterfaceMesh`: The unified mesh to be installed
- `operations::Vector{UnifiedMeshOperation}`: Sequence of operations to build mesh
- `is_feasible::Bool`: Overall feasibility
- `feasibility_issues::Vector{String}`: Issues preventing feasibility
- `edges_from_A::Int`: Count of edges using A's triangulation
- `edges_from_B::Int`: Count of edges using B's triangulation
- `edges_compromised::Int`: Count of edges needing compromise triangulation
- `edges_synthesized::Int`: Count of edges needing synthesis (gap filling)
- `predicted_min_quality::Float64`: Expected worst triangle quality
- `predicted_compatibility_score::Float64`: 0.0 to 1.0, measures compatibility
- `topology::InterfaceTopology`: Original topology (reference)
- `constraints::BoundaryConstraints`: Constraints from both sides

# Statistics
The edge counts break down how the unified mesh is constructed:
- `edges_from_A`: Regions where A's triangulation is superior
- `edges_from_B`: Regions where B's triangulation is superior
- `edges_compromised`: Regions requiring compromise between A and B
- `edges_synthesized`: New triangulation created for gaps or conflicts

# Feasibility
A plan is feasible if:
1. All critical edges can be repaired
2. Quality thresholds are met
3. No irreconcilable conflicts exist
4. Unified mesh is compatible with both A and B

# Example
```julia
plan = SymmetricRepairPlan(
    (4, 5),              # Interface between PID 4 and 5
    sym_mismatches,      # 200 symmetric edge mismatches
    unified_mesh,        # The resulting third mesh
    operations,          # 150 operations to build it
    true,                # feasible
    String[],            # no issues
    85,                  # 85 edges from A
    110,                 # 110 edges from B
    5,                   # 5 compromise edges
    0,                   # no synthesized edges
    0.35,                # min quality acceptable
    0.94,                # high compatibility
    topology,
    constraints
)
```
"""
struct SymmetricRepairPlan
    interface_pair::Tuple{Int,Int}  # (pidA, pidB)
    
    # Bidirectional classification results
    symmetric_mismatches::Vector{SymmetricEdgeMismatch}
    
    # The unified mesh to be generated
    target_unified_mesh::UnifiedInterfaceMesh
    
    # Operations to perform
    operations::Vector{UnifiedMeshOperation}
    
    # Feasibility
    is_feasible::Bool
    feasibility_issues::Vector{String}
    
    # Statistics (breakdown by strategy)
    edges_from_A::Int          # Edges using A's triangulation
    edges_from_B::Int          # Edges using B's triangulation
    edges_compromised::Int     # Edges needing compromise
    edges_synthesized::Int     # Edges needing synthesis (gap filling)
    
    # Quality predictions
    predicted_min_quality::Float64
    predicted_compatibility_score::Float64  # 0.0 to 1.0
    
    # Reference data
    topology::InterfaceTopology
    constraints::BoundaryConstraints
end

"""
    SymmetricRepairPlan(interface_pair, mismatches, topology, constraints)

Create a SymmetricRepairPlan with placeholder unified mesh (to be filled later).
"""
function SymmetricRepairPlan(
    interface_pair::Tuple{Int,Int},
    symmetric_mismatches::Vector{SymmetricEdgeMismatch},
    topology::InterfaceTopology,
    constraints::BoundaryConstraints
)
    return SymmetricRepairPlan(
        interface_pair,
        symmetric_mismatches,
        UnifiedInterfaceMesh(),  # Empty placeholder
        UnifiedMeshOperation[],  # No operations yet
        false,                   # Not feasible until mesh is generated
        ["Unified mesh not yet generated"],
        0, 0, 0, 0,             # No edge statistics yet
        0.0,                    # No quality prediction yet
        0.0,                    # No compatibility score yet
        topology,
        constraints
    )
end

# ============================================================================
# SymmetricClassificationResult - Result of bidirectional classification
# ============================================================================

"""
    SymmetricClassificationResult

Complete result of bidirectional symmetric edge classification.

This structure wraps the classification results and provides comparison metrics
between the two perspectives (A→B and B→A).

# Fields
- `symmetric_mismatches::Vector{SymmetricEdgeMismatch}`: All classified edges
- `classification_AB::InterfaceClassification`: A→B perspective results
- `classification_BA::InterfaceClassification`: B→A perspective results
- `agreement_rate::Float64`: Fraction of edges where both perspectives agree
- `total_unique_edges::Int`: Total unique edges across both perspectives
- `edges_only_in_A::Int`: Edges present only in A
- `edges_only_in_B::Int`: Edges present only in B
- `edges_in_both::Int`: Edges present in both (different triangulation)
- `comparison_metrics::Dict{String,Any}`: Detailed comparison statistics

# Example
```julia
result = SymmetricClassificationResult(
    symmetric_mismatches,  # 200 edges
    classification_AB,      # A→B classification
    classification_BA,      # B→A classification
    0.85,                   # 85% agreement rate
    200,                    # 200 unique edges
    75,                     # 75 only in A
    100,                    # 100 only in B
    25,                     # 25 in both
    comparison_metrics      # Detailed stats
)
```
"""
struct SymmetricClassificationResult
    # Symmetric mismatches (primary result)
    symmetric_mismatches::Vector{SymmetricEdgeMismatch}
    
    # Original directional classifications (for reference)
    classification_AB::InterfaceClassification  # A→B perspective
    classification_BA::InterfaceClassification  # B→A perspective
    
    # Agreement analysis
    agreement_rate::Float64  # 0.0 to 1.0
    
    # Edge distribution
    total_unique_edges::Int
    edges_only_in_A::Int
    edges_only_in_B::Int
    edges_in_both::Int
    
    # Detailed comparison metrics
    comparison_metrics::Dict{String,Any}
end

# ============================================================================
# Utility functions for working with symmetric types
# ============================================================================

"""
    count_by_strategy(mismatches::Vector{SymmetricEdgeMismatch})

Count mismatches by repair strategy.

Returns a dictionary mapping strategy symbols to counts:
- `:use_A`: Number using A's triangulation
- `:use_B`: Number using B's triangulation
- `:compromise`: Number needing compromise
- `:skip`: Number being skipped

# Example
```julia
counts = count_by_strategy(sym_mismatches)
# Dict(:use_A => 85, :use_B => 110, :compromise => 5, :skip => 0)
```
"""
function count_by_strategy(mismatches::Vector{SymmetricEdgeMismatch})::Dict{Symbol,Int}
    counts = Dict{Symbol,Int}(
        :use_A => 0,
        :use_B => 0,
        :compromise => 0,
        :skip => 0
    )
    
    for m in mismatches
        counts[m.repair_strategy] = get(counts, m.repair_strategy, 0) + 1
    end
    
    return counts
end

"""
    compute_agreement_statistics(mismatches::Vector{SymmetricEdgeMismatch})

Compute detailed agreement statistics between A and B perspectives.

Returns a named tuple with:
- `total`: Total number of mismatches
- `both_present`: Edges present in both A and B
- `agree_on_type`: Number agreeing on mismatch type
- `agree_on_feasibility`: Number agreeing on feasibility
- `agreement_rate`: Overall agreement rate (0.0 to 1.0)

# Example
```julia
stats = compute_agreement_statistics(sym_mismatches)
# (total=200, both_present=50, agree_on_type=42, agree_on_feasibility=45, agreement_rate=0.85)
```
"""
function compute_agreement_statistics(mismatches::Vector{SymmetricEdgeMismatch})
    total = length(mismatches)
    both_present = count(m -> m.present_in_both, mismatches)
    agree_on_type = count(m -> m.agree_on_type, mismatches)
    agree_on_feasibility = count(m -> m.agree_on_feasibility, mismatches)
    
    # Agreement rate considers only edges present in both
    agreement_rate = if both_present > 0
        (agree_on_type + agree_on_feasibility) / (2.0 * both_present)
    else
        1.0  # No conflicts if no shared edges
    end
    
    return (
        total = total,
        both_present = both_present,
        agree_on_type = agree_on_type,
        agree_on_feasibility = agree_on_feasibility,
        agreement_rate = agreement_rate
    )
end

"""
    flatten_triangle(tri::Triangle)::NTuple{9,Float64}

Convert a Triangle to flattened coordinate tuple for serialization.

# Example
```julia
tri = Triangle(...)
flat = flatten_triangle(tri)
# (v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z)
```
"""
function flatten_triangle(tri::Triangle)::NTuple{9,Float64}
    return (
        tri.coord1[1], tri.coord1[2], tri.coord1[3],
        tri.coord2[1], tri.coord2[2], tri.coord2[3],
        tri.coord3[1], tri.coord3[2], tri.coord3[3]
    )
end

"""
    unflatten_triangle(flat::NTuple{9,Float64}, elem_id::Int)::Triangle

Convert flattened coordinate tuple back to Triangle (requires creating temporary node IDs).

Note: This creates a Triangle with dummy node IDs (1, 2, 3). Use only for geometry/quality
calculations, not for mesh operations requiring actual node IDs.

# Example
```julia
flat = (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
tri = unflatten_triangle(flat, 1)
```
"""
function unflatten_triangle(flat::NTuple{9,Float64}, elem_id::Int)::Triangle
    c1 = (flat[1], flat[2], flat[3])
    c2 = (flat[4], flat[5], flat[6])
    c3 = (flat[7], flat[8], flat[9])
    
    # Create Triangle with dummy node IDs (1, 2, 3)
    # This is sufficient for geometric/quality calculations
    return Triangle(1, 2, 3, elem_id, c1, c2, c3)
end

"""
    extract_boundary_nodes(triangles::Vector{Triangle})::Set{NTuple{3,Float64}}

Extract all unique node coordinates from a set of triangles.

# Example
```julia
nodes = extract_boundary_nodes(triangles)
# Set of all unique (x, y, z) coordinates
```
"""
function extract_boundary_nodes(triangles::Vector{Triangle})::Set{NTuple{3,Float64}}
    nodes = Set{NTuple{3,Float64}}()
    for tri in triangles
        push!(nodes, tri.coord1, tri.coord2, tri.coord3)
    end
    return nodes
end
