# Coordinate Comparison Guide

**Date**: 2025-10-07  
**Issue**: EdgeKey uses rounded coordinates, Triangle uses unrounded  
**Status**: ✅ RESOLVED - All inappropriate comparisons fixed

## The Core Issue

### Coordinate Storage Differences

In the Nas2Step codebase, coordinates are stored in **two different ways**:

| Structure | Coordinates | Rounding |
|-----------|------------|----------|
| **EdgeKey** | `node1`, `node2` | **Rounded to 4 digits** |
| **Triangle** | `coord1`, `coord2`, `coord3` | **Unrounded (original precision)** |

### Where Rounding Happens

**File**: `src/repair/interface_topology.jl`, line 195

```julia
function ckey(p::NTuple{3,Float64})
    return (round(p[1]; digits=4), round(p[2]; digits=4), round(p[3]; digits=4))
end
```

When building the interface topology:
- **EdgeKey coordinates** are created using `ckey()` which rounds to 4 digits
- **Triangle coordinates** are stored directly from Gmsh without rounding

### Why This Matters

When comparing coordinates from EdgeKey vs Triangle, they may differ slightly:

```julia
edge.node1  # (215.8, 202.646, 43.2351)     ← ROUNDED
tri.coord1  # (215.800003, 202.646011, 43.235065)  ← UNROUNDED
```

Using exact equality (`==`) will fail:
```julia
tri.coord1 == edge.node1  # FALSE! (even though they're the "same" point)
```

## Correct Comparison Methods

### Method 1: Tolerance-Based Comparison (Recommended)

**Use**: When comparing coordinates in geometric algorithms

```julia
function are_nodes_equal(node1::NTuple{3,Float64}, 
                        node2::NTuple{3,Float64}; 
                        tol::Real=1e-4)::Bool
    dx = node1[1] - node2[1]
    dy = node1[2] - node2[2]
    dz = node1[3] - node2[3]
    dist2 = dx*dx + dy*dy + dz*dz
    return dist2 <= tol*tol
end

# Usage
if are_nodes_equal(tri.coord1, edge.node1, tol=1e-4)
    # Coordinates match within tolerance
end
```

**When to use**:
- ✅ Geometric comparisons (e.g., checking if a triangle contains a point)
- ✅ Finding triangles with an edge
- ✅ Checking if vertices match

### Method 2: Rounding-Based Comparison

**Use**: When you need to match the EdgeKey rounding logic exactly

```julia
# Helper: round coordinate to 4 digits (matching EdgeKey rounding)
round_coord(c) = (round(c[1], digits=4), round(c[2], digits=4), round(c[3], digits=4))

# Check both exact match AND rounded match
is_same = (node == edge.node1) || (round_coord(node) == edge.node1)
```

**When to use**:
- ✅ Filtering vertices that are NOT edge endpoints
- ✅ Checking if a coordinate would create the same EdgeKey
- ✅ When you need consistency with topology construction logic

### Method 3: Combined Approach (Most Robust)

**Use**: When you want to handle both exact and near-exact matches

```julia
function is_endpoint(node::NTuple{3,Float64}, edge::EdgeKey)::Bool
    # Helper: round coordinate to 4 digits (matching EdgeKey rounding)
    round_coord(c) = (round(c[1], digits=4), round(c[2], digits=4), round(c[3], digits=4))
    
    # Check both exact match and rounded match
    return (node == edge.node1) || (node == edge.node2) ||
           (round_coord(node) == edge.node1) || (round_coord(node) == edge.node2)
end
```

**When to use**:
- ✅ Critical comparisons where you can't afford false negatives
- ✅ Endpoint detection in hanging node algorithms
- ✅ When coordinates might be from different sources

## Incorrect Comparison (BUG!)

### ❌ DON'T: Use Exact Equality

```julia
# WRONG! Will fail due to rounding differences
if node == edge.node1
    # This will miss matches like:
    # node =  (215.800003, 202.646011, 43.235065)
    # edge.node1 = (215.8, 202.646, 43.2351)
end

# WRONG! Will miss vertices that should match
other_vertices = filter(v -> v != corner1 && v != corner2, quad_vertices)
```

## Fixed Bugs in This Codebase

### Bug #1: False T-Junction Detection (FIXED ✅)

**File**: `src/repair/edge_classification.jl`, line 162-173

**Problem**: Used `node == edge.node1` to check if node is endpoint
- Endpoints with slight coordinate differences were NOT recognized as endpoints
- These "non-endpoints" were treated as hanging nodes
- Led to **32 false T-junction classifications** in test mesh

**Fix**: Added rounding-based comparison
```julia
# Helper: round coordinate to 4 digits (matching EdgeKey rounding)
round_coord(c) = (round(c[1], digits=4), round(c[2], digits=4), round(c[3], digits=4))

is_endpoint = (node == edge.node1) || (node == edge.node2) ||
             (round_coord(node) == edge.node1) || (round_coord(node) == edge.node2)
```

**Result**: T-junction count dropped from 32 to 0 (all were false positives!)

### Bug #2: Quad Vertex Filtering (FIXED ✅)

**File**: `src/repair/repair_planning.jl`, line 557-566

**Problem**: Used `v != corner1 && v != corner2` to filter quad vertices
- Compared unrounded quad vertices with rounded EdgeKey corners
- Could incorrectly include corner vertices in "other_vertices"

**Fix**: Added rounding-based comparison
```julia
round_coord(c) = (round(c[1], digits=4), round(c[2], digits=4), round(c[3], digits=4))

other_vertices = filter(quad_vertices) do v
    not_corner1 = (v != corner1) && (round_coord(v) != corner1)
    not_corner2 = (v != corner2) && (round_coord(v) != corner2)
    not_corner1 && not_corner2
end
```

## Correct Usage Examples

### ✅ Example 1: Checking if Triangle Contains Edge

```julia
function find_triangles_with_edge(edge::EdgeKey, faces::Vector{Triangle}; tol::Real=1e-4)
    triangles_with_edge = Int[]
    
    for (idx, tri) in enumerate(faces)
        # CORRECT: Uses tolerance-based comparison
        has_node1 = triangle_has_node(tri, edge.node1, tol=tol)
        has_node2 = triangle_has_node(tri, edge.node2, tol=tol)
        
        if has_node1 && has_node2
            push!(triangles_with_edge, idx)
        end
    end
    
    return triangles_with_edge
end

# triangle_has_node() internally uses are_nodes_equal() with tolerance ✅
```

### ✅ Example 2: Finding Hanging Nodes

```julia
function find_hanging_nodes_on_edge(edge::EdgeKey, nodes::Set{NTuple{3,Float64}}; tol::Real=1e-4)
    hanging = NTuple{3,Float64}[]
    
    for node in nodes
        # CORRECT: Uses both exact and rounded comparison
        round_coord(c) = (round(c[1], digits=4), round(c[2], digits=4), round(c[3], digits=4))
        
        is_endpoint = (node == edge.node1) || (node == edge.node2) ||
                     (round_coord(node) == edge.node1) || (round_coord(node) == edge.node2)
        
        if is_endpoint
            continue  # Skip endpoints
        end
        
        # Check if node is on the edge interior
        is_on, t, _ = point_on_segment(node, edge.node1, edge.node2, tol=tol)
        if is_on
            push!(hanging, node)
        end
    end
    
    return hanging
end
```

### ✅ Example 3: Matching Edge Nodes

```julia
function edges_match(edges1::Vector{EdgeKey}, edges2::Vector{EdgeKey}; tol::Real=1e-4)
    for e1 in edges1
        found = false
        for e2 in edges2
            # CORRECT: Uses tolerance-based comparison for EdgeKey nodes
            if (are_nodes_equal(e1.node1, e2.node1, tol=tol) && 
                are_nodes_equal(e1.node2, e2.node2, tol=tol)) ||
               (are_nodes_equal(e1.node1, e2.node2, tol=tol) && 
                are_nodes_equal(e1.node2, e2.node1, tol=tol))
                found = true
                break
            end
        end
        if !found
            return false
        end
    end
    return true
end
```

## Design Recommendation

### Why EdgeKey Uses Rounding

The rounding in EdgeKey serves an important purpose:
1. **Hash consistency**: Nearby points hash to the same value
2. **Topology stability**: Slight coordinate differences don't create separate edges
3. **Memory efficiency**: Fewer unique edge keys

### Why Triangle Doesn't Round

Triangles preserve original precision because:
1. **Geometric accuracy**: Need precise coordinates for quality calculations
2. **Mesh fidelity**: Want to preserve original mesh geometry
3. **Export accuracy**: Output mesh should match input precision

### Best Practice

**Always use tolerance-based or rounding-aware comparisons** when comparing:
- EdgeKey coordinates with Triangle coordinates
- Any coordinates from different sources
- When filtering/searching for vertices

**Never use exact equality (`==`)** unless:
- Comparing coordinates from the same source
- You know both are rounded the same way
- You're writing tests with controlled data

## Testing

### Verification Test

```julia
@testset "Coordinate Comparison with Rounding" begin
    # EdgeKey with rounded coordinates
    edge = EdgeKey((215.8, 202.646, 43.2351), (215.8, 224.9031, 37.3982))
    
    # Triangle with unrounded coordinates (as from Gmsh)
    tri = Triangle(1, 2, 3, 1,
                  (215.800003, 202.646011, 43.235065),
                  (215.800003, 224.903076, 37.398232),
                  (216.0, 213.0, 50.0))
    
    # ❌ Exact equality fails
    @test tri.coord1 != edge.node1
    @test tri.coord2 != edge.node2
    
    # ✅ Tolerance-based comparison succeeds
    @test are_nodes_equal(tri.coord1, edge.node1, tol=1e-4)
    @test are_nodes_equal(tri.coord2, edge.node2, tol=1e-4)
    
    # ✅ Rounding-based comparison succeeds
    round_coord(c) = (round(c[1], digits=4), round(c[2], digits=4), round(c[3], digits=4))
    @test round_coord(tri.coord1) == edge.node1
    @test round_coord(tri.coord2) == edge.node2
end
```

## Summary

### Key Takeaways

1. **EdgeKey coordinates are rounded to 4 digits**, Triangle coordinates are not
2. **Never use exact equality (`==`)** to compare EdgeKey vs Triangle coordinates
3. **Always use tolerance-based comparison** for geometric operations
4. **Use rounding-based comparison** when you need to match EdgeKey logic exactly
5. **Combine both approaches** for maximum robustness

### Files with Correct Comparisons ✅

- `src/repair/edge_classification.jl` - Fixed hanging node detection
- `src/repair/repair_planning.jl` - Fixed quad vertex filtering
- `src/repair/geometric_utilities.jl` - Provides `are_nodes_equal()`

### Common Functions

- `are_nodes_equal(node1, node2; tol=1e-4)` - Tolerance-based comparison
- `triangle_has_node(tri, node; tol=1e-4)` - Check if triangle contains node
- `find_triangles_with_edge(edge, faces; tol=1e-4)` - Find triangles with edge

---

**Remember**: When in doubt, use tolerance-based comparison!
