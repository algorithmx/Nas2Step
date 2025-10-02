#!/usr/bin/env julia

"""
Generate an Umbilic Torus volumetric mesh using gmsh.

Approach:
1. Generate boundary point cloud from parametric equations
2. Use gmsh API to create surface mesh
3. Use gmsh API to generate volumetric mesh
4. Export to NAS format

The Umbilic Torus is a (3,1) torus knot defined by:
    x = sin(u) * (7 + cos(u/3 - 2v) + 2cos(u/3 + v))
    y = cos(u) * (7 + cos(u/3 - 2v) + 2cos(u/3 + v))
    z = sin(u/3 - 2v) + 2sin(u/3 + v)
where u ∈ [0, 2π] and v ∈ [0, 2π]
"""

using Gmsh: gmsh

"""
    umbilic_torus_parametric(u, v)

Calculate a point on the Umbilic Torus surface.

# Arguments
- `u`: First parameter (0 to 2π) - NOT periodic!
- `v`: Second parameter (0 to 2π) - periodic

# Returns
- Tuple (x, y, z) coordinates
"""
function umbilic_torus_parametric(u, v)
    r = 7 + cos(u/3 - 2*v) + 2*cos(u/3 + v)
    x = sin(u) * r
    y = cos(u) * r
    z = sin(u/3 - 2*v) + 2*sin(u/3 + v)
    return (x, y, z)
end


"""
    generate_boundary_points(u_res=60, v_res=40)

Generate boundary point cloud for Umbilic Torus from parametric equations.

# Arguments
- `u_res`: Resolution in u direction (default: 60)
- `v_res`: Resolution in v direction (default: 40)

# Returns
- Vector of (x,y,z) tuples representing the boundary point cloud
"""
function generate_boundary_points(u_res=60, v_res=40)
    # Generate parameter space
    # u ∈ [0, 2π] traces the complete surface (not periodic!)
    # v ∈ [0, 2π] is periodic
    u_vals = range(0, stop=2π, length=u_res+1)[1:end-1]  # Exclude endpoint
    v_vals = range(0, stop=2π, length=v_res+1)[1:end-1]  # Exclude endpoint
    
    # Generate boundary point cloud
    points = []
    for u in u_vals
        for v in v_vals
            pt = umbilic_torus_parametric(u, v)
            push!(points, pt)
        end
    end
    
    println("  Generated $(length(points)) boundary points")
    return points
end


"""
    demo_umbilic_torus(; u_res=60, v_res=40, mesh_size=0.5, output="demo_umbilic_torus.nas")

Generate Umbilic Torus volumetric mesh using gmsh.

Steps:
1. Generate boundary point cloud from parametric equations
2. Add points to gmsh
3. Use gmsh to create surface mesh (Delaunay or other algorithm)
4. Use gmsh to fill interior with tetrahedral mesh
5. Export to NAS format

# Arguments
- `u_res`: Resolution in u direction (default: 60)
- `v_res`: Resolution in v direction (default: 40)
- `mesh_size`: Target mesh element size for gmsh (default: 0.5)
- `output`: Output NAS filename (default: "demo_umbilic_torus.nas")
"""
function demo_umbilic_torus(; u_res=60, v_res=40, mesh_size=0.5, output="demo_umbilic_torus.nas")
    println("=== Umbilic Torus Volumetric Mesh Generation ===")
    println("  Resolution: u=$u_res, v=$v_res")
    println("  Target mesh size: $mesh_size")
    println()
    
    # Step 1: Generate boundary point cloud from parametric equations
    println("Step 1: Generating boundary point cloud...")
    points = generate_boundary_points(u_res, v_res)
    
    # Initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)  # 0=silent, 1=errors, 2=warnings, 3=info
    gmsh.model.add("umbilic_torus")
    
    try
        # Step 2: Add points to gmsh
        println()
        println("Step 2: Adding points to gmsh...")
        point_tags = []
        for (idx, pt) in enumerate(points)
            tag = gmsh.model.geo.addPoint(pt[1], pt[2], pt[3], mesh_size, idx)
            push!(point_tags, tag)
        end
        println("  Added $(length(point_tags)) points")
        
        gmsh.model.geo.synchronize()
        
        # Step 3: Create structured connectivity for surface
        # Since we have a structured grid, create quads with proper wrapping
        println()
        println("Step 3: Creating surface connectivity...")
        
        # Helper function to get point index from (i,j)
        # Only j (v-direction) wraps periodically, not i (u-direction)
        function point_index(i, j)
            return (i * v_res + (j % v_res)) + 1
        end
        
        # Find optimal cyclic shift for u-boundary connection
        println("  Finding optimal shift for u-boundary...")
        best_shift = 0
        best_dist = Inf
        for shift in 0:(v_res-1)
            total_dist = 0.0
            for j in 0:(v_res-1)
                p_last = points[point_index(u_res-1, j)]
                p_first = points[point_index(0, (j + shift) % v_res)]
                total_dist += sqrt(sum((p_last .- p_first).^2))
            end
            if total_dist < best_dist
                best_dist = total_dist
                best_shift = shift
            end
        end
        println("    Optimal shift: $best_shift, avg gap: $(round(best_dist/v_res, digits=3))")
        
        # Create lines for structured quad mesh
        lines = Dict{Tuple{Int,Int},Int}()
        line_tag = 1
        
        function get_or_create_line(p1, p2)
            key = p1 < p2 ? (p1, p2) : (p2, p1)
            if !haskey(lines, key)
                gmsh.model.geo.addLine(key[1], key[2], line_tag)
                lines[key] = line_tag
                line_tag += 1
            end
            # Return with correct orientation
            return lines[key] * (p1 < p2 ? 1 : -1)
        end
        
        # Create surface patches as triangles (split quads)
        surface_tags = []
        surf_tag = 1
        
        for i in 0:(u_res-1)
            for j in 0:(v_res-1)
                # Get the 4 corners of the quad
                if i < u_res - 1
                    # Regular quads
                    p1 = point_index(i, j)
                    p2 = point_index(i+1, j)
                    p3 = point_index(i+1, j+1)
                    p4 = point_index(i, j+1)
                else
                    # Boundary bridge quads with optimal shift
                    p1 = point_index(i, j)
                    p2 = point_index(0, (j + best_shift) % v_res)
                    p3 = point_index(0, (j + 1 + best_shift) % v_res)
                    p4 = point_index(i, j+1)
                end
                
                # Triangle 1: p1-p2-p3
                l1 = get_or_create_line(p1, p2)
                l2 = get_or_create_line(p2, p3)
                l3 = get_or_create_line(p3, p1)
                
                loop1 = gmsh.model.geo.addCurveLoop([l1, l2, l3], surf_tag)
                s1 = gmsh.model.geo.addPlaneSurface([loop1], surf_tag)
                push!(surface_tags, s1)
                surf_tag += 1
                
                # Triangle 2: p1-p3-p4
                l4 = get_or_create_line(p1, p3)
                l5 = get_or_create_line(p3, p4)
                l6 = get_or_create_line(p4, p1)
                
                loop2 = gmsh.model.geo.addCurveLoop([l4, l5, l6], surf_tag)
                s2 = gmsh.model.geo.addPlaneSurface([loop2], surf_tag)
                push!(surface_tags, s2)
                surf_tag += 1
            end
        end
        
        println("  Created $(length(surface_tags)) surface triangles")
        
        gmsh.model.geo.synchronize()
        
        # Step 4: Create volume from closed surface
        println()
        println("Step 4: Creating volume from closed surface...")
        surface_loop = gmsh.model.geo.addSurfaceLoop(surface_tags, 1)
        volume = gmsh.model.geo.addVolume([surface_loop], 1)
        
        gmsh.model.geo.synchronize()
        
        # Set mesh size
        println("  Setting mesh size parameters...")
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
        
        # Step 5: Generate volumetric mesh
        println()
        println("Step 5: Generating volumetric tetrahedral mesh...")
        gmsh.model.mesh.generate(3)
        
        # Optimize mesh quality
        println("  Optimizing mesh quality...")
        gmsh.model.mesh.optimize("Netgen")
        
        # Step 6: Write to NAS format
        println()
        println("Step 6: Writing to NAS format...")
        gmsh.write(output)
        
        # Get mesh statistics
        node_tags, _, _ = gmsh.model.mesh.getNodes()
        n_nodes = length(node_tags)
        
        elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3, -1)
        n_tets = 0
        for (etype, etags) in zip(elem_types, elem_tags)
            if etype == 4  # Tetrahedra
                n_tets += length(etags)
            end
        end
        
        println()
        println("=== Mesh Statistics ===")
        println("  Total nodes: $n_nodes")
        println("  Total tetrahedra: $n_tets")
        println("  Output file: $output")
        
    finally
        gmsh.finalize()
    end
    
    println()
    println("=== Complete ===")
    return output
end


# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    demo_umbilic_torus()
end
