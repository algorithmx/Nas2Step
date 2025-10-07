# Analysis of UNKNOWN Classification and Refinement Proposals

## Current State: UNKNOWN is Over-Populated

### Observed Statistics (NC_Reduction_4.nas):
- Interface 1/14: **83 UNKNOWN**, 0 other types
- Interface 3/14: **6 UNKNOWN**, 0 other types  
- Interface 4/14: **3 UNKNOWN**, 0 other types
- Interface 7/14: **257 UNKNOWN**, 10 DIAGONAL
- Interface 8/14: **78 UNKNOWN**, 4 DIAGONAL
- Interface 9/14: **10 UNKNOWN**, 1 DIAGONAL
- Interface 10/14: **35 UNKNOWN**, 1 DIAGONAL
- Interface 11/14: **82 UNKNOWN**, 5 DIAGONAL
- Interface 13/14: **1 UNKNOWN**, 0 other types

**Total**: ~556 UNKNOWN mismatches vs ~21 DIAGONAL, 0 T-JUNCTION, 0 REFINEMENT, 0 QUAD_MISMATCH

## Current Classification Logic Flow

```
Edge mismatch found
  ↓
Check hanging nodes (nodes lying ON the edge)
  ├─ 1 hanging node → T_JUNCTION
  ├─ >1 hanging nodes → REFINEMENT  
  └─ 0 hanging nodes (clean edge) →
       ↓
     Check if both endpoints are shared vertices
       ├─ Both shared → Try find_quad_for_diagonal()
       │    ├─ Found quad + matching boundaries → DIAGONAL
       │    ├─ Found quad + mismatched boundaries → QUAD_MISMATCH
       │    └─ No quad found → **UNKNOWN** (Path 1)
       └─ Endpoints not shared → **UNKNOWN** (Path 2)
```

## Why UNKNOWN is Over-Populated

### Path 1: "Quad finding failed" (Lines 700-710)
**Triggers when**:
- No hanging nodes (clean edge)
- Both endpoints are shared vertices
- BUT `find_quad_for_diagonal()` returns empty `([], [])`

**Why this happens**:
1. Edge exists in SOURCE but only in 1 triangle (not 2) → not a diagonal
2. Edge is on the boundary of source mesh
3. Target mesh has no triangles using those 4 vertices
4. Source mesh has >2 or <2 triangles with the edge

**What these edges actually are**:
- **BOUNDARY_EDGE**: Edge on the boundary of one mesh, not truly shared
- **ISOLATED_EDGE**: Edge exists but not part of quad structure
- **NON_CONFORMING_EDGE**: Edge exists but geometric structure doesn't allow matching

### Path 2: "Vertices not shared" (Lines 711-723)
**Triggers when**:
- No hanging nodes
- One or both endpoints are NOT in shared vertex set

**What this means**:
- Fundamental topology issue
- Edge endpoints don't actually exist on both sides
- This should be caught earlier in vertex conformity check!

## Logical Completeness Assessment

### Well-Covered Cases ✅:
1. **T_JUNCTION**: Single hanging node → subdivide edge
2. **REFINEMENT**: Multiple hanging nodes → hierarchical refinement
3. **DIAGONAL**: Valid quad with matching boundaries → retriangulate
4. **QUAD_MISMATCH**: Same vertices, different quad topology

### Under-Classified Cases ❌:
1. **BOUNDARY_EDGE**: Edge on mesh boundary, not internal interface
2. **PARTIAL_EDGE**: Edge partially exists (one endpoint shared, one not)
3. **NON_PLANAR**: 4 vertices exist but not coplanar
4. **DEGENERATE**: Edge with zero or near-zero length
5. **MULTIPLICITY**: Edge appears more than twice in source (non-manifold)

## Proposed Refinements

### New Mismatch Types

```julia
@enum MismatchType begin
    T_JUNCTION          # One region has intermediate node, other doesn't
    DIAGONAL            # Same VALID quad, different triangulation
    REFINEMENT          # Hierarchical mesh refinement difference
    QUAD_MISMATCH       # Same 4 vertices but different quad boundary topology
    BOUNDARY_EDGE       # Edge on mesh boundary, not true interface edge
    NON_MANIFOLD        # Edge shared by != 2 triangles (topology issue)
    UNSHARED_ENDPOINT   # One or both endpoints not in shared vertex set
    UNKNOWN             # Cannot classify (catch-all)
end
```

### Refined Classification Logic

```julia
# After checking hanging nodes = 0...

# NEW: Check if edge is on boundary (appears in only 1 triangle in source)
source_triangles_with_edge = find_triangles_with_edge(edge, source_faces)
if length(source_triangles_with_edge) == 1
    mismatch_type = BOUNDARY_EDGE
    complexity = 0.9
    # Cannot repair boundary edges - they're not true interface edges
    return ...
end

# NEW: Check for non-manifold edges (>2 triangles share this edge)
if length(source_triangles_with_edge) > 2
    mismatch_type = NON_MANIFOLD
    complexity = 0.95
    # Topology is broken - cannot repair
    return ...
end

# Check if both endpoints are shared
endpoint1_shared = is_vertex_in_set(edge.node1, shared_vertices)
endpoint2_shared = is_vertex_in_set(edge.node2, shared_vertices)

# NEW: Distinguish unshared endpoints
if !endpoint1_shared || !endpoint2_shared
    mismatch_type = UNSHARED_ENDPOINT
    complexity = 0.95
    # Log which endpoint is not shared for diagnostics
    return ...
end

# Proceed with quad finding (existing logic)...
```

### Enhanced Diagnostic Information

For UNKNOWN cases that remain, capture diagnostic data:

```julia
struct UnknownDiagnostics
    source_triangle_count::Int
    target_triangle_count::Int
    has_affected_triangles::Bool
    endpoints_shared::(Bool, Bool)
    quad_vertices_found::Bool
    reason::String
end
```

## Implementation Priority

### High Priority (Immediate Impact):
1. **BOUNDARY_EDGE detection**: Check `length(source_triangles_with_edge) == 1`
   - Expected to capture 50-70% of current UNKNOWN cases
   - Clear diagnostic: "Edge is on mesh boundary, not an interface edge"

2. **UNSHARED_ENDPOINT detection**: Explicit check instead of falling through to UNKNOWN
   - Captures vertex-level topology issues
   - Should have been caught in Phase 0, investigate why it wasn't

### Medium Priority (Refine Diagnostics):
3. **NON_MANIFOLD detection**: Check `length(source_triangles_with_edge) > 2`
   - Rare but indicates serious topology problems
   - Cannot be repaired automatically

4. **Enhanced UNKNOWN diagnostics**: Capture why quad finding failed
   - Log source/target triangle counts
   - Log what geometric check failed

### Low Priority (Edge Cases):
5. **NON_PLANAR detection**: Check if 4 vertices are coplanar
6. **DEGENERATE detection**: Check edge length

## Expected Outcomes

### Before:
- UNKNOWN: ~556 (95% of mismatches)
- Classified: ~21 (5%)

### After (Estimated):
- **BOUNDARY_EDGE**: ~350-400 (60-70%)
- **DIAGONAL**: ~21 (same)
- **UNSHARED_ENDPOINT**: ~50-100 (10-15%)
- **NON_MANIFOLD**: ~10-20 (2-3%)
- **UNKNOWN**: ~50-100 (10-15% - irreducible)

## Testing Strategy

1. Add debug logging to track which UNKNOWN path is taken
2. Implement BOUNDARY_EDGE detection first
3. Re-run NC_Reduction_4.nas and compare statistics
4. Use debug mode on a few UNKNOWN cases to understand geometry
5. Iterate on remaining UNKNOWN cases

## Code Changes Required

### Files to Modify:
1. **`src/repair/edge_classification.jl`**:
   - Add BOUNDARY_EDGE, NON_MANIFOLD, UNSHARED_ENDPOINT to enum
   - Add early checks before quad finding
   - Update statistics tracking

2. **`src/repair/repair_planning.jl`**:
   - Handle new mismatch types (all non-repairable)
   - Add clear error messages for each type

3. **Statistics and reporting**:
   - Update InterfaceClassification struct
   - Update JSON export

## Success Criteria

1. **UNKNOWN reduced to <20%** of total mismatches
2. **Each mismatch type has clear meaning** and diagnostic message
3. **Users understand why repairs cannot proceed** (not "unknown" reason)
4. **No false classifications** (validate with known test cases)
