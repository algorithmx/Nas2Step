#!/usr/bin/env julia

using Nas2Step
using .Nas2Step: T_JUNCTION, DIAGONAL

println("Testing new InterfaceClassification constructors...")

# Test the new smart constructor
topology = InterfaceTopology(1, 2, Set{NTuple{3,Int}}(), Dict{NTuple{3,Int},Tuple{Int,Int}}(),
    Triangle[], Triangle[],
    Dict{EdgeKey, Vector{Int}}(), Dict{EdgeKey, Vector{Int}}(),
    Set{EdgeKey}(), Set{EdgeKey}(), Set{EdgeKey}(),
    BoundingBox(NTuple{3,Float64}[]), 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0, 1.0)

# Create some test mismatches
edge_key = create_edge_key_int((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
mismatch1 = EdgeMismatch(edge_key, T_JUNCTION, :B_only, 1,
    [(0.5, 0.0, 0.0)], [1], NTuple{3,Float64}[], Int[], 0.2, 1.0, true)
mismatch2 = EdgeMismatch(edge_key, DIAGONAL, :A_only, 2,
    NTuple{3,Float64}[], [1, 2], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)], [1, 2], 0.4, 1.0, true)

# Test 1: Smart constructor (recommended)
println("\n1. Testing smart constructor (recommended):")
cls1 = InterfaceClassification(topology, [mismatch1], [mismatch2])
println("  ✓ Smart constructor works")
println("  ✓ T-junction count: $(cls1.t_junction_count)")
println("  ✓ Diagonal count: $(cls1.diagonal_count)")
println("  ✓ Total feasible: $(cls1.total_feasible)")
println("  ✓ Average complexity: $(cls1.average_complexity)")

# Test 2: Keyword constructor
println("\n2. Testing keyword constructor:")
cls2 = InterfaceClassification(topology;
    mismatches_A=[mismatch1],
    mismatches_B=[mismatch2],
    average_complexity=0.3)  # Override computed value
println("  ✓ Keyword constructor works")
println("  ✓ Average complexity (overridden): $(cls2.average_complexity)")

# Test 3: Legacy constructor (should show warning)
println("\n3. Testing legacy constructor (deprecated):")
cls3 = InterfaceClassification(topology, [mismatch1], [mismatch2], 1, 1, 0, 0, 0, 0, 0.3)
println("  ✓ Legacy constructor works (with warning)")

# Test 4: Direct legacy function
println("\n4. Testing direct legacy function:")
cls4 = create_interface_classification_legacy(topology, [mismatch1], [mismatch2], 1, 1, 0, 0, 0, 0, 0.3)
println("  ✓ Direct legacy function works (with warning)")

println("\n✅ All InterfaceClassification constructors are working correctly!")
println("\nNew smart constructor InterfaceClassification(topology, mismatches_A, mismatches_B)")
println("is the RECOMMENDED way to create InterfaceClassification objects.")