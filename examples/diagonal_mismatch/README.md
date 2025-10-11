# Diagonal Mismatch Example

This example demonstrates two touching boxes with vertex-compatible interfaces where triangulation differs in a diagonal manner.

## Overview

- **Box A (PID 1)**: Uses diagonal from bottom-front to top-back
- **Box B (PID 2)**: Uses opposite diagonal from top-front to bottom-back
- **Shared Interface**: Both boxes share exactly 4 vertices but use different triangulations

## Files

- `generate_diagonal_mismatch.jl` - Main script to generate the example and create NAS file
- `README.md` - This documentation file

## Running the Example

1. Generate the NAS file:
```bash
julia generate_diagonal_mismatch.jl
```

2. Repair the generated NAS file:
```bash
cd ../../
julia repair_mesh.jl examples/diagonal_mismatch/diagonal_mismatch_example.nas
```

## Expected Output

The example will:
1. Create two boxes with different diagonal triangulations
2. Generate a Nastran NAS file with both meshes
3. Detect and classify the diagonal mismatch
4. Show visualization of the issue
5. The repair script will then fix the mismatch by choosing one triangulation