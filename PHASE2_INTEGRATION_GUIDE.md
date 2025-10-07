# Phase 2 Integration Guide

This guide explains how to integrate the new symmetric repair data structures into the Nas2Step module.

---

## Step 1: Update Module File

Add the new symmetric repair types file to `src/Nas2Step.jl`:

```julia
# src/Nas2Step.jl

module Nas2Step

# ... existing includes ...

# Repair subsystem
include("repair/tolerance_config.jl")
include("repair/geometric_utilities.jl")
include("repair/interface_topology.jl")
include("repair/edge_classification.jl")
include("repair/boundary_constraints.jl")

# NEW: Add symmetric repair types
include("repair/symmetric_repair_types.jl")

include("repair/interface_analysis_export.jl")
include("repair/interface_conformity_check.jl")
include("repair/repair_workspace.jl")
include("repair/repair_planning.jl")
include("repair/repair_execution.jl")
include("repair/repair_verification.jl")

# ... rest of module ...

# Export symmetric repair types (optional - can be accessed via Nas2Step.SymmetricEdgeMismatch)
export SymmetricEdgeMismatch
export UnifiedInterfaceMesh
export UnifiedMeshOperation
export SymmetricRepairPlan
export SymmetricClassificationResult

# Export utility functions
export count_by_strategy
export compute_agreement_statistics
export flatten_triangle
export unflatten_triangle
export extract_boundary_nodes

end # module
```

**Important**: The include order matters! Place `symmetric_repair_types.jl` **after**:
- `interface_topology.jl` (provides `EdgeKey`, `Triangle`, `InterfaceTopology`)
- `edge_classification.jl` (provides `EdgeMismatch`, `MismatchType`)
- `boundary_constraints.jl` (provides `BoundaryConstraints`)

---

## Step 2: Verify Dependencies

The symmetric repair types depend on these existing types:

### From `interface_topology.jl`:
- `EdgeKey` - Used in `SymmetricEdgeMismatch`
- `Triangle` - Used in `UnifiedInterfaceMesh` and `UnifiedMeshOperation`
- `InterfaceTopology` - Used in `SymmetricRepairPlan`

### From `edge_classification.jl`:
- `EdgeMismatch` - Wrapped by `SymmetricEdgeMismatch`
- `InterfaceClassification` - Used in `SymmetricClassificationResult`
- `MismatchType` - Referenced for type agreement checks

### From `boundary_constraints.jl`:
- `BoundaryConstraints` - Used in `SymmetricRepairPlan`

All these types are already defined in the existing codebase, so no changes are needed.

---

## Step 3: Test Integration

Create a simple integration test to verify the module loads correctly:

```julia
# test/test_integration_symmetric_types.jl

using Test
using Nas2Step

@testset "Symmetric Types Integration" begin
    # Test that types are accessible
    @test isdefined(Nas2Step, :SymmetricEdgeMismatch)
    @test isdefined(Nas2Step, :UnifiedInterfaceMesh)
    @test isdefined(Nas2Step, :UnifiedMeshOperation)
    @test isdefined(Nas2Step, :SymmetricRepairPlan)
    @test isdefined(Nas2Step, :SymmetricClassificationResult)
    
    # Test that utility functions are accessible
    @test isdefined(Nas2Step, :count_by_strategy)
    @test isdefined(Nas2Step, :compute_agreement_statistics)
    @test isdefined(Nas2Step, :flatten_triangle)
    @test isdefined(Nas2Step, :unflatten_triangle)
    @test isdefined(Nas2Step, :extract_boundary_nodes)
    
    println("✓ All symmetric types successfully integrated!")
end
```

Run with:
```bash
julia --project=. -e 'using Pkg; Pkg.test("Nas2Step")'
```

---

## Step 4: Usage in External Code

Once integrated, users can access the types in two ways:

### Option 1: Direct import (if exported)
```julia
using Nas2Step

# Types are directly available
edge = EdgeKey((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
sym = SymmetricEdgeMismatch(edge, nothing, nothing, :skip, 0.0, "test")
unified = UnifiedInterfaceMesh()
```

### Option 2: Qualified access (always available)
```julia
using Nas2Step

# Access via module prefix
edge = Nas2Step.EdgeKey((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
sym = Nas2Step.SymmetricEdgeMismatch(edge, nothing, nothing, :skip, 0.0, "test")
unified = Nas2Step.UnifiedInterfaceMesh()
```

---

## Step 5: Update Documentation

Add entries to the module documentation explaining the new types:

### In README.md or docs/:

```markdown
## Symmetric Mesh Repair

The symmetric mesh repair system generates a unified interface mesh compatible with both sides:

- **SymmetricEdgeMismatch**: Stores bidirectional edge classifications
- **UnifiedInterfaceMesh**: Represents the third mesh that replaces both interfaces
- **UnifiedMeshOperation**: Atomic operation for building the unified mesh
- **SymmetricRepairPlan**: Complete repair plan for symmetric approach
- **SymmetricClassificationResult**: Classification results with comparison metrics

See `PHASE2_COMPLETION_SUMMARY.md` for detailed usage examples.
```

---

## Step 6: Compatibility with Existing Code

The new types are designed to **coexist** with the existing unidirectional repair system:

### Existing unidirectional approach (still works):
```julia
# Old API - still functional
topology = build_interface_topology(mesh_file, pidA, pidB)
classification = classify_interface_mismatches(topology)
plan = generate_repair_plan(topology, classification, constraints)
```

### New symmetric approach (when ready):
```julia
# New API - uses symmetric types
topology = build_interface_topology(mesh_file, pidA, pidB)
sym_classification = classify_interface_mismatches_symmetric(topology)  # Phase 3
sym_plan = generate_symmetric_repair_plan(topology, sym_classification, constraints)  # Phase 4-5
```

Both approaches can coexist during the transition period.

---

## Step 7: Update Exports List

Consider which symbols should be exported. Recommendation:

### Always export (core functionality):
- `SymmetricEdgeMismatch`
- `UnifiedInterfaceMesh`
- `SymmetricRepairPlan`
- `SymmetricClassificationResult`

### Maybe export (utilities):
- `count_by_strategy`
- `compute_agreement_statistics`
- `flatten_triangle` / `unflatten_triangle`

### Don't export (internal):
- `UnifiedMeshOperation` (used internally by mesh generation)
- `extract_boundary_nodes` (utility, not primary API)

This keeps the public API clean while allowing advanced users to access everything via `Nas2Step.symbol`.

---

## Step 8: Precompilation Check

After integration, verify the module precompiles correctly:

```bash
julia --project=. -e 'using Nas2Step'
```

If precompilation fails, check:
1. All dependencies are properly loaded
2. Include order is correct
3. No circular dependencies exist
4. All types are fully defined before use

---

## Troubleshooting

### Issue: "UndefVarError: EdgeKey not defined"

**Cause**: `symmetric_repair_types.jl` included before `interface_topology.jl`

**Solution**: Reorder includes to ensure dependencies come first

### Issue: "Method overwrite warning"

**Cause**: Function name collision with existing code

**Solution**: Rename utility functions with unique prefix (e.g., `symmetric_count_by_strategy`)

### Issue: "Cannot construct SymmetricEdgeMismatch with EdgeMismatch"

**Cause**: Type mismatch between expected and actual EdgeMismatch type

**Solution**: Verify `edge_classification.jl` exports correct EdgeMismatch type

---

## Verification Checklist

After integration, verify:

- [ ] Module loads without errors: `using Nas2Step`
- [ ] Types are accessible: `Nas2Step.SymmetricEdgeMismatch`
- [ ] Constructors work: `UnifiedInterfaceMesh()`
- [ ] Utility functions work: `count_by_strategy([])`
- [ ] Documentation is updated
- [ ] Tests pass: `Pkg.test("Nas2Step")`
- [ ] No precompilation warnings
- [ ] Existing functionality unaffected

---

## Next Phase Integration

When implementing Phase 3 (symmetric classification), you'll add:

```julia
# In src/Nas2Step.jl, after symmetric_repair_types.jl
include("repair/symmetric_classification.jl")

# Export new functions
export classify_interface_mismatches_symmetric
export determine_repair_strategy
```

Each subsequent phase will follow a similar pattern.

---

## Quick Integration Command Sequence

```bash
# 1. Ensure you're in the project directory
cd /Users/dabajabaza/Documents/Workspace/Nas2Step

# 2. Add the include line to src/Nas2Step.jl (manual edit)

# 3. Test that it loads
julia --project=. -e 'using Nas2Step; println("✓ Module loaded successfully")'

# 4. Run unit tests
julia --project=. test/test_symmetric_repair_types.jl

# 5. Run full test suite
julia --project=. -e 'using Pkg; Pkg.test("Nas2Step")'
```

---

## Summary

Phase 2 integration is **straightforward** because:

1. ✅ No external dependencies added
2. ✅ No changes to existing types
3. ✅ Clean separation from existing code
4. ✅ Backward compatible

Simply add the include statement, export desired symbols, and the new types are ready to use!

**Status**: Ready for integration
