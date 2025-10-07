#!/usr/bin/env julia
# Diagnostic script to verify T-junction classification

using Pkg
Pkg.activate(".")

using Nas2Step
using Gmsh: gmsh

"""
Diagnose T-junction classifications to ensure they're correct.
"""
function diagnose_tjunctions(nas_file::String, pidA::Int, pidB::Int; limit=10)
    println("="^70)
    println("T-Junction Classification Diagnostic")
    println("="^70)
    println("File: $nas_file")
    println("Interface: PID $pidA ↔ PID $pidB")
    println()
    
    # Build topology
    println("Building interface topology...")
    topology = build_interface_topology(nas_file, pidA, pidB, tol=1e-4)
    
    println("Interface statistics:")
    println("  Vertices shared: $(length(topology.shared_node_keys))")
    println("  Triangles in A: $(length(topology.faces_A))")
    println("  Triangles in B: $(length(topology.faces_B))")
    println("  Edges only in A: $(length(topology.edges_only_in_A))")
    println("  Edges only in B: $(length(topology.edges_only_in_B))")
    println()
    
    # Analyze T-junctions
    println("Analyzing T-junction classifications...")
    println()
    
    tjunction_count = 0
    examples_shown = 0
    
    for (i, edge) in enumerate(topology.edges_only_in_A)
        if examples_shown >= limit
            break
        end
        
        # Classify the edge
        mismatch = Nas2Step.classify_edge_mismatch(edge, topology, :A, tol=1e-4, debug=false)
        
        if mismatch.mismatch_type == Nas2Step.T_JUNCTION
            tjunction_count += 1
            
            if examples_shown < limit
                examples_shown += 1
                
                println("="^70)
                println("T-JUNCTION EXAMPLE #$examples_shown")
                println("="^70)
                println("Edge endpoints:")
                println("  Node1: $(edge.node1)")
                println("  Node2: $(edge.node2)")
                println()
                
                # Compute edge length
                dx = edge.node2[1] - edge.node1[1]
                dy = edge.node2[2] - edge.node1[2]
                dz = edge.node2[3] - edge.node1[3]
                edge_length = sqrt(dx*dx + dy*dy + dz*dz)
                println("Edge length: $edge_length")
                println()
                
                # Get target nodes
                target_faces = topology.faces_B
                target_nodes = Set{NTuple{3,Float64}}()
                for tri in target_faces
                    push!(target_nodes, tri.coord1, tri.coord2, tri.coord3)
                end
                
                println("Total target vertices: $(length(target_nodes))")
                println()
                
                # Find hanging nodes manually with detailed output
                hanging = NTuple{3,Float64}[]
                
                println("Checking each target vertex:")
                check_count = 0
                for node in target_nodes
                    # Check if node is endpoint (exact match)
                    is_endpoint_exact = (node == edge.node1) || (node == edge.node2)
                    
                    # Check if node is endpoint (tolerance match)
                    dist1 = sqrt((node[1]-edge.node1[1])^2 + (node[2]-edge.node1[2])^2 + (node[3]-edge.node1[3])^2)
                    dist2 = sqrt((node[1]-edge.node2[1])^2 + (node[2]-edge.node2[2])^2 + (node[3]-edge.node2[3])^2)
                    is_endpoint_tol = (dist1 < 1e-4) || (dist2 < 1e-4)
                    
                    # Check if on segment
                    is_on, t, dist_from_line = Nas2Step.point_on_segment(node, edge.node1, edge.node2, tol=1e-4)
                    
                    if is_on || is_endpoint_exact || is_endpoint_tol
                        check_count += 1
                        if check_count <= 5  # Only print first 5
                            status = if is_on
                                "ON INTERIOR"
                            elseif is_endpoint_exact
                                "ENDPOINT (exact)"
                            elseif is_endpoint_tol
                                "ENDPOINT (tolerance)"
                            else
                                "UNKNOWN"
                            end
                            
                            println("  Node: $node")
                            println("    Status: $status")
                            println("    Dist to node1: $(round(dist1, digits=8))")
                            println("    Dist to node2: $(round(dist2, digits=8))")
                            println("    Parameter t: $(round(t, digits=6))")
                            println("    Dist from line²: $(round(dist_from_line, digits=12))")
                            println()
                        end
                        
                        if is_on && !is_endpoint_exact
                            push!(hanging, node)
                        end
                    end
                end
                
                if check_count > 5
                    println("  ... and $(check_count - 5) more vertices checked")
                    println()
                end
                
                println("Hanging nodes found: $(length(hanging))")
                if !isempty(hanging)
                    println("Hanging node coordinates:")
                    for (i, h) in enumerate(hanging)
                        println("  $i: $h")
                    end
                    println()
                end
                
                # Verify this is actually a T-junction
                println("Verification:")
                if length(hanging) == 1
                    println("  ✓ CORRECT: Exactly 1 hanging node → Valid T-junction")
                elseif length(hanging) == 0
                    println("  ✗ ERROR: No hanging nodes found → Should NOT be T-junction!")
                    println("  → This suggests endpoint detection is failing")
                else
                    println("  ✗ ERROR: $(length(hanging)) hanging nodes → Should be REFINEMENT, not T-junction!")
                end
                println()
            end
        end
    end
    
    println("="^70)
    println("SUMMARY")
    println("="^70)
    println("Total edges in A: $(length(topology.edges_only_in_A))")
    println("Classified as T-junctions: $tjunction_count")
    println("Examples shown: $examples_shown")
    
    # Also check edges from B
    println()
    println("Checking edges from B side...")
    tjunction_count_b = 0
    for edge in topology.edges_only_in_B
        mismatch = Nas2Step.classify_edge_mismatch(edge, topology, :B, tol=1e-4, debug=false)
        if mismatch.mismatch_type == Nas2Step.T_JUNCTION
            tjunction_count_b += 1
        end
    end
    println("Edges in B classified as T-junctions: $tjunction_count_b")
    println()
    println("Total T-junctions reported: $(tjunction_count + tjunction_count_b)")
end

# Main execution
if length(ARGS) < 3
    println("""
Usage: julia diagnose_tjunctions.jl <nas_file> <pidA> <pidB> [limit]

Example:
  julia --project diagnose_tjunctions.jl examples/realistic/NC_Reduction_4.nas 1 3 5
""")
    exit(1)
end

nas_file = ARGS[1]
pidA = parse(Int, ARGS[2])
pidB = parse(Int, ARGS[3])
limit = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 10

try
    diagnose_tjunctions(nas_file, pidA, pidB, limit=limit)
catch e
    println("\n✗ Error during diagnosis:")
    println(e)
    rethrow()
end
