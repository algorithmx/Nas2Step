#!/usr/bin/env julia
using Pkg
Pkg.activate(".")

using Nas2Step

println("="^70)
println("TOLERANCE INVESTIGATION")
println("="^70)

# Load mesh
mesh_file = length(ARGS) >= 1 ? ARGS[1] : "./examples/realistic/NC_Reduction_4.nas"
println("\nLoading mesh: $mesh_file")

# Build topology directly from file
topo = build_interface_topology(mesh_file, 1, 3, tol=1e-4)

println("\nTopology:")
println("  Faces A: $(length(topo.faces_A))")
println("  Faces B: $(length(topo.faces_B))")
println("  Shared vertices: $(length(topo.shared_node_keys))")

# Get one diagonal mismatch for analysis
println("\n" * "="^70)
println("ANALYZING FIRST DIAGONAL MISMATCH")
println("="^70)

# Classify to find diagonal mismatches
classification = classify_interface_mismatches(topo, tol=1e-4, debug=false)

diagonal_mismatches = filter(m -> m.mismatch_type == Nas2Step.DIAGONAL, 
                              vcat(classification.mismatches_A, classification.mismatches_B))

if isempty(diagonal_mismatches)
    println("\n❌ No diagonal mismatches found!")
    exit(0)
end

mismatch = diagonal_mismatches[1]
println("\nDiagonal mismatch:")
println("  Edge: $(mismatch.edge_key.node1) → $(mismatch.edge_key.node2)")
println("  Present in: $(mismatch.present_in)")
println("  Should be in PID: $(mismatch.should_be_in)")
println("  Quad vertices found: $(length(mismatch.quad_vertices))")
println("  Target triangles: $(mismatch.triangles_to_replace)")

# Analyze why quad vertices aren't found in target
if isempty(mismatch.quad_vertices)
    println("\n❌ No quad vertices - edge might not be a proper diagonal")
else
    println("\nQuad vertices from SOURCE mesh:")
    for (i, v) in enumerate(mismatch.quad_vertices)
        println("  V$i: $v")
    end
    
    # Determine target faces
    target_faces = mismatch.should_be_in == topo.pidA ? topo.faces_A : topo.faces_B
    println("\nSearching TARGET mesh ($(length(target_faces)) faces)...")
    
    # For each quad vertex, find closest match in target
    println("\nVertex matching analysis:")
    for (i, qv) in enumerate(mismatch.quad_vertices)
        closest_dist = Inf
        closest_vertex = nothing
        matches_within_tol = 0
        
        for tri in target_faces
            for tv in [tri.coord1, tri.coord2, tri.coord3]
                dx = tv[1] - qv[1]
                dy = tv[2] - qv[2]
                dz = tv[3] - qv[3]
                dist = sqrt(dx*dx + dy*dy + dz*dz)
                
                if dist < closest_dist
                    closest_dist = dist
                    closest_vertex = tv
                end
                
                if dist <= 1e-4
                    matches_within_tol += 1
                end
            end
        end
        
        status = matches_within_tol > 0 ? "✓" : "✗"
        println("  $status V$i: closest_dist=$(closest_dist), matches=$(matches_within_tol)")
        if closest_dist > 1e-6
            println("      Quad: $qv")
            println("      Closest: $closest_vertex")
            println("      Diff: $((closest_vertex[1]-qv[1], closest_vertex[2]-qv[2], closest_vertex[3]-qv[3]))")
        end
    end
    
    # Check if ANY triangles use ALL 4 quad vertices
    println("\nChecking for triangles using these 4 vertices...")
    
    function is_quad_vertex(v, quad_verts, tol=1e-4)
        for qv in quad_verts
            dx = v[1] - qv[1]
            dy = v[2] - qv[2]
            dz = v[3] - qv[3]
            if dx*dx + dy*dy + dz*dz <= tol*tol
                return true
            end
        end
        return false
    end
    
    matching_triangles = Int[]
    for (idx, tri) in enumerate(target_faces)
        tri_verts = [tri.coord1, tri.coord2, tri.coord3]
        if all(v -> is_quad_vertex(v, mismatch.quad_vertices), tri_verts)
            push!(matching_triangles, idx)
        end
    end
    
    println("  Found $(length(matching_triangles)) triangles using only these 4 vertices")
    if length(matching_triangles) > 0
        println("  Triangle indices: $matching_triangles")
    end
end

println("\n" * "="^70)
println("DIAGNOSIS")
println("="^70)

if !isempty(mismatch.quad_vertices) && isempty(mismatch.triangles_to_replace)
    println("""
    The quad vertices were successfully extracted from the SOURCE mesh,
    but no triangles in the TARGET mesh use exactly those 4 vertices.
    
    This suggests:
    1. The vertices aren't perfectly shared (precision mismatch), OR
    2. The target mesh has additional vertices/triangulation, OR
    3. The tolerance (1e-4) is too strict for this mesh
    
    Recommendation: Try increasing tolerance or investigate vertex alignment.
    """)
end

println("\n✓ Analysis complete!")
