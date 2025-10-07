# Enhanced Mismatch Classification Plan

## Problem Analysis

### Current Issue
The current `DIAGONAL` classification is **too permissive** and fails to verify that both sides form valid, consistently-oriented quads. Specifically:

1. **Current logic** (lines 515-520 in `edge_classification.jl`):
   - Finds quad vertices in SOURCE mesh (where diagonal edge exists)
   - Checks if target mesh has triangles using those 4 vertices
   - If found → classifies as `DIAGONAL`
   
2. **What's missing**:
   - **No verification** that the 4 vertices form a proper PLANAR quad on both sides
   - **No check** that all 4 boundary edges of the quad exist on BOTH sides
   - **No validation** of quad orientation/winding consistency
   - **Result**: Classifies non-conforming cases as DIAGONAL, leading to "no affected triangles" errors

### Example of Misclassification

**Side A**: Vertices {1, 2, 3, 4} with edges: 1-2, 2-3, 3-4, 4-1, **1-3** (diagonal)
- Forms proper quad with triangles sharing diagonal 1-3

**Side B**: Vertices {1, 2, 3, 4} with edges: 1-3, 3-2, 2-4, 4-1, **3-4**
- Has same 4 vertices BUT quad is "twisted" (different boundary edges)
- Edge 3-4 exists but is NOT the opposite diagonal
- This should **NOT** be classified as DIAGONAL

## Root Cause

The "No affected triangles found for edge insertion" error occurs because:
1. Edge is classified as DIAGONAL
2. Repair planning tries quad retriangulation
3. But the quad structure doesn't actually match between sides
4. Target side doesn't have proper triangulation to retrieve

## Proposed Solution

### New Mismatch Types

Add new enum values to `MismatchType`:

```julia
@enum MismatchType begin
    T_JUNCTION          # One region has intermediate node, other doesn't
    DIAGONAL            # Same VALID quad, different triangulation (strict check)
    REFINEMENT          # Hierarchical mesh refinement difference
    QUAD_MISMATCH       # Same 4 vertices but different quad topology
    BOUNDARY_MISMATCH   # Edge on boundary, not truly interface edge
    UNKNOWN             # Cannot classify
end
```

### Enhanced Validation Functions

#### 1. **Verify Quad Validity** (`verify_quad_structure`)
Check that 4 vertices form a valid planar quad with correct boundary edges:
- Extract boundary edges from 2 triangles forming the quad
- Verify exactly 4 boundary edges (the quad perimeter)
- Check planarity (4 points should be coplanar within tolerance)

#### 2. **Verify Quad Consistency** (`verify_quad_consistency_across_sides`)
For DIAGONAL classification, verify BOTH sides have matching quad structure:
- Both sides must have the SAME 4 boundary edges
- Both sides must triangulate the quad (2 triangles using those 4 vertices)
- Only the diagonal should differ between sides

#### 3. **Classify Quad Mismatch** 
When 4 vertices are shared but quad structures don't match:
- Extract boundary edges from both sides
- Compare: if different → `QUAD_MISMATCH`
- If one side has extra internal edges → complexity issue

### Implementation Steps

#### Step 1: Add `QUAD_MISMATCH` type
- Update `MismatchType` enum
- Update classification statistics tracking
- Update reporting

#### Step 2: Implement quad validation functions
```julia
function extract_quad_boundary_edges(tri1, tri2, diagonal; tol=1e-4)
    # Extract the 4 edges forming quad perimeter (excluding diagonal)
    # Returns: Vector{EdgeKey} of boundary edges
end

function verify_quad_planarity(vertices; tol=1e-4)
    # Check if 4 vertices are coplanar
    # Returns: (is_planar::Bool, max_deviation::Float64)
end

function get_quad_boundary_edges_both_sides(edge, topology, present_in, quad_vertices; tol=1e-4)
    # Extract boundary edges from both source and target
    # Returns: (source_boundary_edges, target_boundary_edges)
end

function verify_quad_boundary_match(source_edges, target_edges; tol=1e-4)
    # Check if both sides have same 4 boundary edges
    # Returns: (match::Bool, missing_edges, extra_edges)
end
```

#### Step 3: Enhance `find_quad_for_diagonal`
Add validation that returns richer diagnostics:
```julia
function find_quad_for_diagonal(...)
    # ... existing logic ...
    
    # NEW: Verify quad validity
    is_valid, reason = verify_quad_structure(quad_vertices, source_triangles)
    if !is_valid
        return ([], [], :invalid_quad, reason)
    end
    
    # NEW: Verify consistency across sides
    match, details = verify_quad_consistency_across_sides(...)
    if !match
        return (quad_vertices, [], :quad_mismatch, details)
    end
    
    return (quad_vertices, target_triangles, :valid, "")
end
```

#### Step 4: Update classification logic
In `classify_edge_mismatch`, handle new return format:
```julia
quad_vertices, triangles_to_replace, status, details = find_quad_for_diagonal(...)

if status == :valid
    mismatch_type = DIAGONAL
    complexity = 0.4
elseif status == :quad_mismatch
    mismatch_type = QUAD_MISMATCH
    complexity = 0.8
    # Store diagnostic details
else  # :invalid_quad or other
    mismatch_type = UNKNOWN
    complexity = 0.9
end
```

#### Step 5: Update repair planning
Skip repair for `QUAD_MISMATCH` type (cannot be fixed by simple retriangulation):
```julia
if mismatch.mismatch_type == QUAD_MISMATCH
    return create_failed_plan(1, "Quad topology mismatch: boundary edges don't match")
end
```

## Expected Outcomes

1. **Accurate classification**: Only true diagonal mismatches classified as `DIAGONAL`
2. **Clear diagnostics**: `QUAD_MISMATCH` clearly indicates the problem
3. **Reduced failures**: Fewer "no affected triangles" errors
4. **Better reporting**: Users understand why repairs can't be automated

## Testing Strategy

1. Create unit tests with synthetic quad mismatch cases
2. Re-run NC_Reduction_4.nas to see reclassification
3. Verify that true diagonal cases still work
4. Check that QUAD_MISMATCH cases are properly identified

## Files to Modify

1. `src/repair/edge_classification.jl` - Core classification logic
2. `src/repair/repair_planning.jl` - Handle new mismatch type
3. `test/unit/test_edge_classification.jl` - Add test cases
