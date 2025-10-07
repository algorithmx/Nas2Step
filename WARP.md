# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Overview

Nas2Step is a **Julia package** that converts NASTRAN (.nas) mesh files to STEP CAD files using Gmsh.jl. The package handles both surface-only and volume-only meshes, automatically extracting boundary surfaces from tetrahedral volume meshes when needed.

## Essential Commands

### Basic Usage
```bash
# Run conversion on a NAS file
julia main.jl input.nas

# Run conversion with verification 
julia -e "using Nas2Step; main(ARGS[1]; is_verify=true)" input.nas
```

### Development and Testing
```bash
# Start Julia REPL with project
julia --project=.

# Run example demonstrations
julia examples/box/demo_box.jl
julia examples/torus/demo_torus.jl
julia examples/umbilic/demo_umbilic_torus.jl

# Install dependencies
julia --project=. -e "using Pkg; Pkg.instantiate()"

# Update dependencies
julia --project=. -e "using Pkg; Pkg.update()"
```

### Package Development
```bash
# Activate development mode
julia --project=. -e "using Pkg; Pkg.develop(PackageSpec(path=\".\"))"

# Run with module in REPL
julia --project=. -e "using Nas2Step"
```

## Architecture Overview

### Core Conversion Pipeline
The main conversion process follows these stages:
1. **Boundary extraction** (`extract.jl`) - Extracts external faces from volume meshes
2. **NAS to STEP conversion** (`nas_to_step.jl`) - Core conversion using shared-topology OCC reconstruction
3. **Mesh analysis** (`analyze.jl`) - Volume checking and region validation
4. **Verification** (`verify.jl`) - Quality checks for element volumes and mesh integrity

### Key Components

#### Mesh Processing (`extract.jl`)
- `extract_boundary_surfaces()` - Identifies external boundaries from CTETRA elements
- Uses face adjacency analysis to distinguish between external and internal interfaces
- Handles non-conformal meshes by including internal interface faces

#### Conversion Engine (`nas_to_step.jl`)
- `convert_nas_to_step()` - Main conversion function using Gmsh OCC backend
- Builds shared topology with edge/face deduplication
- Creates volumetric shells when possible, falls back to open surfaces
- Generates manifoldness anomaly reports in JSON format

#### Data Structures (`write_nas.jl`)
- `SurfaceRegion` - Surface mesh with thickness (PSHELL elements)
- `VolumeRegion` - Volume mesh (CTETRA elements)
- `write_nas_surface()` / `write_nas_volume()` - Generate free-field NASTRAN format

#### Quality Analysis (`analyze.jl`, `verify.jl`)
- `check_nas_volumes()` - Splits meshes by PID and validates boundary surfaces
- `check_element_volumes()` - Detects inverted elements (negative volumes)
- `verify_step()` - Inspects output STEP file geometry

### Geometric Processing Details

The conversion uses Gmsh's OCC (OpenCASCADE) backend for robust geometric operations:
- **Node deduplication**: Merges duplicate nodes within tolerance (1e-10)
- **Shared topology**: Creates single OCC edges/faces shared by adjacent regions
- **Manifold analysis**: Checks edge usage (boundary=1, manifold=2, non-manifold>2)
- **Feature classification**: Uses angle-based surface classification (default 40°)

### Error Handling and Diagnostics

The package includes comprehensive diagnostics:
- **Volume validation**: Checks if regions form closed manifold shells
- **Boundary analysis**: Reports boundary edges and non-manifold conditions
- **Element quality**: Detects inverted or degenerate tetrahedra
- **Anomaly reporting**: Exports manifoldness issues to JSON for debugging

## Working with Dependencies

### Required Packages
- **Gmsh.jl**: Core meshing and geometry operations
- **Statistics.jl**: Statistical analysis of mesh quality

### Version Constraints
- Gmsh.jl: 0.3.1
- Statistics.jl: 1.11.1

The project uses strict version pinning for reproducibility.

## File Organization

```
src/
├── Nas2Step.jl       # Main module with exports
├── nas_to_step.jl    # Core conversion pipeline
├── extract.jl        # Boundary surface extraction
├── write_nas.jl      # NASTRAN file generation utilities
├── analyze.jl        # Volume and mesh analysis
├── verify.jl         # Quality verification functions
└── inspect.jl        # STEP file inspection utilities

examples/
├── box/              # Multi-region conformal hex-to-tet example
├── torus/            # Parametric torus surface meshing
└── umbilic/          # Complex umbilic torus geometry
```

## Common Development Patterns

### Adding New Mesh Types
1. Extend `VolumeRegion` or `SurfaceRegion` structs in `write_nas.jl`
2. Update export functions to handle new element types
3. Add boundary extraction logic in `extract.jl` if needed

### Debugging Conversion Issues
1. Use `check_nas_volumes()` to validate input mesh by PID
2. Enable anomaly JSON output to identify manifoldness problems
3. Use `verify_step()` to inspect output geometry
4. Check element volumes with `check_element_volumes()` for inversions

### Performance Optimization
- The extraction phase uses efficient substring parsing to avoid allocation overhead
- Face adjacency uses canonical edge/face keys for fast lookup
- Gmsh operations are batched where possible to minimize OCC synchronization

### Testing New Features
Use the example generators as test cases:
- `demo_box.jl` tests multi-region conformal interfaces
- `demo_torus.jl` tests curved surface boundary extraction
- The examples generate reproducible test meshes with known geometric properties