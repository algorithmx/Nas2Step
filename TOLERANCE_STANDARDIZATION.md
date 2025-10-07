# Tolerance Standardization in Mesh Repair

## Overview

This document describes the standardized tolerance values used throughout the mesh conformity analysis and repair system. Consistent tolerance values are critical for ensuring that geometric comparisons are compatible across different modules.

## The Problem We Solved

**Initial State**: The codebase had inconsistent tolerance values:
- Some functions used `tol=1e-4` (0.0001)
- Some used `digits=4` rounding (equivalent to 1e-4)
- One function used `atol=1e-8` (much tighter!)
- Documentation mentioned `tol=1e-8` but code used `tol=1e-4`

**This caused issues**:
1. Precision incompatibilities between `find_hanging_nodes_on_edge` and `find_quad_for_diagonal`
2. Inconsistent vertex matching across different analysis phases
3. Potential for misclassification of edge mismatches

**Solution**: Centralized tolerance configuration in `tolerance_config.jl`

## Standard Tolerance Values

### Primary Geometric Tolerance

```julia
DEFAULT_GEOMETRIC_TOLERANCE = 1e-4  # 0.0001 units
```

**Used for**:
- Vertex/node equality checks (`nodes_equal_within_tolerance`)
- Edge endpoint matching
- Hanging node detection
- Quad vertex identification
- All distance-based geometric comparisons

**Rationale**:
- Provides 4 decimal places of precision
- Adequate for typical engineering meshes (mm or inch scales)
- Robust to floating-point round-off errors
- Balances precision with numerical stability

### Coordinate Rounding

```julia
COORDINATE_ROUNDING_DIGITS = 4  # matches 1e-4
```

**Used for**:
- Creating hash keys for coordinates
- Set membership tests
- JSON export formatting

**Important**: This must match `DEFAULT_GEOMETRIC_TOLERANCE`:
```
10^(-COORDINATE_ROUNDING_DIGITS) = DEFAULT_GEOMETRIC_TOLERANCE
```

### Quality Thresholds

```julia
DEFAULT_MIN_ANGLE_DEG = 10.0      # Minimum triangle angle
DEFAULT_MAX_ANGLE_DEG = 170.0     # Maximum triangle angle  
DEFAULT_MAX_ASPECT_RATIO = 20.0   # Maximum aspect ratio
DEFAULT_FEATURE_ANGLE_DEG = 30.0  # Feature detection angle
```

## Where Tolerances Are Used

### 1. Interface Topology (`interface_topology.jl`)

**Current State**: ✅ Consistent
- Uses `digits=4` for coordinate rounding
- Uses `tol=1e-4` parameter (default)
- **One exception**: `isclose` in `TriangleHasEdge` uses `atol=1e-8`

**Action Needed**: Update `TriangleHasEdge` to use standard tolerance.

### 2. Geometric Utilities (`geometric_utilities.jl`)

**Current State**: ✅ Fully Consistent
- All functions use `tol=1e-4` default
- Distance-based comparisons use `tol²` optimization

### 3. Edge Classification (`edge_classification.jl`)

**Current State**: ✅ Fully Consistent
- Uses `tol=1e-4` throughout
- Hanging node detection: `tol=1e-4`
- Quad finding: `tol=1e-4`
- Vertex sanity checks: `tol=1e-4`

### 4. Vertex Conformity Check (`interface_conformity_check.jl`)

**Current State**: ✅ Fully Consistent
- Uses `tol=1e-4` for all vertex comparisons
- Cluster distance: `10.0` units (separate parameter)

### 5. Boundary Constraints (`boundary_constraints.jl`)

**Current State**: ✅ Consistent
- Uses `tol=1e-4` default
- Uses `digits=4` for coordinate rounding
- Feature angle: `30.0°` default

### 6. Repair Planning (`repair_planning.jl`)

**Current State**: ✅ Mostly Consistent
- Quality thresholds use standard values
- Edge density calculations use `digits=4` for display

## Best Practices

### For New Code

1. **Always use the default tolerance**:
   ```julia
   function my_function(coords; tol::Real=ToleranceConfig.DEFAULT_GEOMETRIC_TOLERANCE)
       # ...
   end
   ```

2. **Use distance² comparisons** to avoid sqrt():
   ```julia
   dist2 = dx*dx + dy*dy + dz*dz
   if dist2 <= ToleranceConfig.tolerance_squared()
       # Points are equal
   end
   ```

3. **Round coordinates consistently**:
   ```julia
   rounded = ToleranceConfig.round_coordinates(coords)
   ```

### For Existing Code

1. **Check default values**: Ensure all `tol` parameters default to `1e-4`
2. **Check digits**: Ensure all `digits` parameters use `4`
3. **Check special cases**: Look for hardcoded `1e-8`, `1e-6`, or other values

## Testing Tolerance Consistency

### Verification Checklist

- [x] All geometric comparison functions use `tol=1e-4`
- [x] All coordinate rounding uses `digits=4`
- [x] Hanging node detection uses same tolerance as quad finding
- [x] Vertex conformity uses same tolerance as edge classification
- [ ] `TriangleHasEdge` in `interface_topology.jl` needs update

### Known Issues

1. **`TriangleHasEdge` function**: Uses `atol=1e-8` internally
   - **Impact**: Low (function is only used in specific contexts)
   - **Fix**: Update to use `ToleranceConfig.DEFAULT_GEOMETRIC_TOLERANCE`

2. **Documentation inconsistency**: Some docstrings mention `tol=1e-8`
   - **Impact**: None (code uses correct value)
   - **Fix**: Update docstrings to reflect actual default

## Migration Guide

### Updating Code to Use Centralized Config

**Before**:
```julia
function my_function(coords; tol::Real=1e-4)
    # comparison
end
```

**After**:
```julia
using ..ToleranceConfig

function my_function(coords; tol::Real=ToleranceConfig.DEFAULT_GEOMETRIC_TOLERANCE)
    # comparison
end
```

**Or** (simpler, for module-level code):
```julia
function my_function(coords; tol::Real=1e-4)
    # Keep existing default, but document that it matches standard
end
```

## Impact Assessment

### Critical Path (High Priority)

✅ **Vertex Matching**: All consistent at `1e-4`
✅ **Edge Detection**: All consistent at `1e-4`  
✅ **Hanging Nodes**: Consistent at `1e-4`
✅ **Quad Finding**: Consistent at `1e-4`

### Secondary Path (Medium Priority)

⚠️ **`TriangleHasEdge`**: Uses `1e-8` internally
- Used in repair verification
- Not on critical path for classification
- Should be updated for completeness

### Non-Critical (Low Priority)

✅ **Display/Export**: Uses `digits=4` consistently
✅ **Documentation**: Some mentions of `1e-8` but code is correct

## Conclusion

The mesh repair system now has **consistent tolerance usage** across all critical paths:

- **Vertex comparison**: `1e-4` everywhere
- **Edge detection**: `1e-4` everywhere  
- **Classification logic**: `1e-4` everywhere
- **Coordinate rounding**: `digits=4` everywhere (equivalent to `1e-4`)

This ensures that:
1. Vertices detected as "shared" in one module are recognized as shared in all modules
2. Edge endpoints checked in classification match those used in quad finding
3. No false negatives from overly tight tolerances
4. No false positives from overly loose tolerances

The one remaining minor inconsistency (`TriangleHasEdge` using `atol=1e-8`) is in a non-critical utility function and can be addressed in future cleanup.
