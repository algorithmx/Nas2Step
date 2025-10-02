# Mesh Quality Checks in Nas2Step

This document describes the mesh quality verification capabilities added to the Nas2Step module.

## Overview

The module now includes comprehensive mesh quality checking with **unified JSON output**. All quality metrics are exported to a single JSON file for easy integration with other tools.

## Features

### 1. Vertex Coordination Check
Examines the coordination number (number of connected elements) for each vertex.

**Anomalies detected:**
- **Undercoordinated vertices**: < 3 connections (potential mesh defects)
- **Overcoordinated vertices**: > 90th percentile (default) - adaptive threshold based on mesh distribution

### 2. Volume Quality Checks
- Inverted elements (negative volumes)
- Zero-volume (degenerate) elements
- Volume statistics (min, max, mean, median)
- Convention mismatch detection (RHR vs LHR)

### 3. Element Quality Metrics
- Radius ratio (gamma metric)
- Poor quality element detection
- Quality distribution statistics

### 4. Surface Closure
- Checks if mesh forms closed manifolds
- Boundary edge detection

### 5. Conversion Anomalies (Optional)
- Manifoldness issues from NAS to STEP conversion
- Non-manifold edges
- Boundary edges in regions

## Usage

### Quick Check - Individual Functions

```julia
using Nas2Step

# Check vertex coordination only (90th percentile for overcoordination)
coord = check_vertex_coordination("mesh.nas", min_coord=3, overcoord_percentile=95)

# Check element volumes
vols = check_element_volumes("mesh.nas")

# Check element quality
qual = check_element_quality("mesh.nas")
```

### Comprehensive Check - Unified Report

```julia
using Nas2Step

# Run ALL checks and export unified JSON
result = comprehensive_mesh_check(
    "mesh.nas",
    output_json="quality_report.json",
    run_conversion=false,  # Set true to include conversion anomalies
    verbose=true
)

# Access results
println("Status: ", result.verification.overall_status)
println("Report: ", result.json_path)
```

### Advanced - With Conversion Testing

```julia
# Run all checks PLUS test NAS→STEP conversion
result = comprehensive_mesh_check(
    "mesh.nas",
    output_json="full_quality_report.json",
    run_conversion=true,
    step_output="test_output.step",
    verbose=true
)

# This includes manifoldness anomalies from the conversion
if result.anomaly_path !== nothing
    println("Conversion anomalies: ", result.anomaly_path)
end
```

## JSON Output Structure

The unified JSON report contains:

```json
{
  "overall_status": "ok|warning|error|convention_mismatch",
  
  "volume_check": {
    "total_elements": 1000,
    "inverted_count": 0,
    "zero_volume_count": 0,
    "min_volume": 0.001,
    "max_volume": 0.01,
    "convention_mismatch": false,
    "status": "ok"
  },
  
  "quality_check": {
    "min": 0.4,
    "max": 0.99,
    "mean": 0.75,
    "median": 0.76,
    "poor_count": 5,
    "bad_count": 0,
    "status": "ok"
  },
  
  "coordination_check": {
    "total_vertices": 500,
    "mean_coordination": 12.5,
    "median_coordination": 12.0,
    "p90_coordination": 24.0,
    "min_coord_threshold": 3,
    "max_coord_threshold": 24,
    "overcoord_percentile": 90,
    "undercoordinated_count": 2,
    "overcoordinated_count": 50,
    "undercoordinated_vertices": [
      {"vertex_id": 42, "coordination": 2}
    ],
    "overcoordinated_vertices": [
      {"vertex_id": 123, "coordination": 30}
    ],
    "status": "warning"
  },
  
  "closure_check": {
    "is_closed": true,
    "boundary_edge_count": 0,
    "surface_count": 0,
    "status": "ok"
  },
  
  "conversion_anomalies_file": "mesh_anomalies.json",
  "conversion_anomalies_note": "Manifoldness issues from conversion..."
}
```

## Exported Functions

All functions are exported from the `Nas2Step` module:

- `check_vertex_coordination(filename; min_coord=3, overcoord_percentile=95)`
- `check_element_volumes(filename; swap_nodes=nothing)`
- `check_element_quality(filename; metric="gamma")`
- `check_surface_closure(filename)`
- `verify_nas_mesh(filename; verbose=true)` - Runs all checks
- `export_mesh_quality_json(result, output_file; include_anomalies_from=nothing)`
- `comprehensive_mesh_check(nas_file; output_json, run_conversion, verbose)` - **Recommended**

## Integration Notes

### With Previous Functionality
The module already had `write_anomalies_json()` for conversion manifoldness issues. The new system:
- **Keeps separate anomaly JSON** from `nas_to_step()` conversion
- **References it** in the unified quality report
- **Merges all other checks** into one JSON file

### JSON Output Philosophy
- **One main report**: All verification metrics in `mesh_quality_report.json`
- **Optional reference**: Points to separate `*_anomalies.json` if conversion was tested
- **No external dependencies**: Uses built-in JSON writer (no JSON3.jl required)

## Example Workflow

```bash
# Run the demo
cd examples
julia --project=.. check_coordination_demo.jl

# Examine the unified report
cat box/unified_quality_report.json
```

## Interpreting Results

### Status Levels
- `ok`: All checks passed
- `warning`: Minor issues detected (e.g., overcoordinated vertices)
- `error`: Serious defects (e.g., undercoordinated vertices, inverted elements)
- `convention_mismatch`: Node ordering convention difference (not necessarily a defect)

### Coordination Numbers
- **Typical range**: 4-30 for tetrahedral meshes (depends on mesh density)
- **< 3**: Serious defect (isolated vertices, hanging edges)
- **> P90**: Automatically flagged based on distribution (top 10% outliers)
- **Very high (> 40)**: May indicate mesh quality issues or refinement boundaries
- **Adaptive threshold**: Uses 90th percentile so it adapts to each mesh

### Recommended Actions
1. **Undercoordinated vertices** → Inspect and fix mesh topology
2. **Inverted elements** → Check node ordering convention
3. **Poor quality elements** → Consider remeshing
4. **Non-closed surfaces** → Fix boundary extraction or meshing
