#!/usr/bin/env julia
# Diagnostic script to understand "Quad vertices not triangulated in target" issue

using Pkg
Pkg.activate(".")

using Nas2Step
using Gmsh: gmsh

"""
Diagnose why quad vertices from source are not found in target triangulation.
"""
function diagnose_quad_issue(nas_file::String, pidA::Int, pidB::Int; debug_limit=5)
    println("="^70)
    println("Diagnosing Quad Vertex Mismatch Issue")
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
    
    # Find edges classified as "Quad vertices not triangulated in target"
    println("Analyzing edges in A that are missing in B...")
    
    count_quad_not_triangulated = 0
    examples_shown = 0
    
    for (i, edge) in enumerate(topology.edges_only_in_A)
        if examples_shown >= debug_limit
            break
        end
        
        # Classify the edge
        mismatch = Nas2Step.classify_edge_mismatch(edge, topology, :A, tol=1e-4, debug=false)
        
        if mismatch.mismatch_type == Nas2Step.TARGET_USES_FINER_TRIANGULATION
            count_quad_not_triangulated += 1
            
            if examples_shown < debug_limit
                examples_shown += 1
                
                println("\n" * "="^70)
                println("EXAMPLE #$examples_shown: Edge with quad vertices not in target")
                println("="^70)
                println("Edge endpoints:")
                println("  Node1: $(edge.node1)")
                println("  Node2: $(edge.node2)")
                println()
                
                # Manually find the quad vertices in source
                source_faces = topology.faces_A
                source_tris_with_edge = Nas2Step.find_triangles_with_edge(edge, source_faces, tol=1e-4)
                
                println("Source triangles containing edge: $(length(source_tris_with_edge))")
                
                if length(source_tris_with_edge) == 2
                    tri1 = source_faces[source_tris_with_edge[1]]
                    tri2 = source_faces[source_tris_with_edge[2]]
                    
                    println("\nTriangle 1:")
                    println("  V1: $(tri1.coord1)")
                    println("  V2: $(tri1.coord2)")
                    println("  V3: $(tri1.coord3)")
                    
                    println("\nTriangle 2:")
                    println("  V1: $(tri2.coord1)")
                    println("  V2: $(tri2.coord2)")
                    println("  V3: $(tri2.coord3)")
                    
                    # Extract quad vertices
                    quad_verts = Nas2Step.extract_quad_vertices(tri1, tri2, tol=1e-4)
                    
                    println("\nExtracted quad vertices ($(length(quad_verts)) vertices):")
                    for (i, v) in enumerate(quad_verts)
                        println("  V$i: $v")
                    end
                    
                    # Check if these vertices are in the shared set
                    println("\nChecking if quad vertices are in shared vertex set:")
                    for (i, v) in enumerate(quad_verts)
                        is_shared = Nas2Step.is_vertex_in_set(v, topology.shared_node_keys, tol=1e-4)
                        println("  V$i: $(is_shared ? "✓ SHARED" : "✗ NOT SHARED")")
                    end
                    
                    # Now check target triangles
                    target_faces = topology.faces_B
                    println("\nSearching for triangles in target using these vertices...")
                    println("Total target triangles: $(length(target_faces))")
                    
                    # Find triangles using quad vertices
                    target_tris_using_quad = Nas2Step.find_triangles_using_vertices(quad_verts, target_faces, tol=1e-4)
                    println("Triangles found: $(length(target_tris_using_quad))")
                    
                    if isempty(target_tris_using_quad)
                        println("\n⚠ NO TRIANGLES FOUND IN TARGET USING QUAD VERTICES")
                        println("\nLet's check manually what triangles touch the quad vertices...")
                        
                        for (i, qv) in enumerate(quad_verts)
                            touching_tris = []
                            for (idx, tri) in enumerate(target_faces)
                                if Nas2Step.triangle_has_node(tri, qv, tol=1e-4)
                                    push!(touching_tris, idx)
                                end
                            end
                            println("  V$i: Found in $(length(touching_tris)) target triangles")
                        end
                        
                        # Check if any target triangle uses ALL 3 of its vertices from the quad set
                        println("\nChecking if any target triangle uses 3 vertices from quad set...")
                        for (idx, tri) in enumerate(target_faces)
                            tri_verts = [tri.coord1, tri.coord2, tri.coord3]
                            
                            function is_quad_vertex(v)
                                for qv in quad_verts
                                    if Nas2Step.are_nodes_equal(v, qv, tol=1e-4)
                                        return true
                                    end
                                end
                                return false
                            end
                            
                            matches = [is_quad_vertex(v) for v in tri_verts]
                            if all(matches)
                                println("  Triangle $idx: ALL 3 vertices match quad set!")
                                println("    V1: $(tri.coord1) - match: $(matches[1])")
                                println("    V2: $(tri.coord2) - match: $(matches[2])")
                                println("    V3: $(tri.coord3) - match: $(matches[3])")
                            end
                        end
                    else
                        println("\n✓ Found $(length(target_tris_using_quad)) triangles in target")
                        for idx in target_tris_using_quad[1:min(3, end)]
                            tri = target_faces[idx]
                            println("  Triangle $idx:")
                            println("    V1: $(tri.coord1)")
                            println("    V2: $(tri.coord2)")
                            println("    V3: $(tri.coord3)")
                        end
                    end
                end
            end
        end
    end
    
    println("\n" * "="^70)
    println("SUMMARY")
    println("="^70)
    println("Total edges in A missing from B: $(length(topology.edges_only_in_A))")
    println("Edges classified as 'Target uses finer triangulation': $count_quad_not_triangulated")
    println("Examples shown: $examples_shown")
end

# Main execution
if length(ARGS) < 3
    println("""
Usage: julia diagnose_quad_mismatch.jl <nas_file> <pidA> <pidB>

Example:
  julia --project diagnose_quad_mismatch.jl examples/realistic/NC_Reduction_4.nas 1 3
""")
    exit(1)
end

nas_file = ARGS[1]
pidA = parse(Int, ARGS[2])
pidB = parse(Int, ARGS[3])

try
    diagnose_quad_issue(nas_file, pidA, pidB)
catch e
    println("\n✗ Error during diagnosis:")
    println(e)
    rethrow()
end
