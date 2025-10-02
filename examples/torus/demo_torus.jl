#!/usr/bin/env julia

"""
Generate an ordinary torus volumetric mesh using gmsh.

The ordinary torus (donut shape) is defined by:
    x = (R + r*cos(v)) * cos(u)
    y = (R + r*cos(v)) * sin(u)
    z = r * sin(v)
where:
    R = major radius (distance from center to tube center)
    r = minor radius (tube radius)
    u ∈ [0, 2π] - angle around the major circle (NOT periodic)
    v ∈ [0, 2π] - angle around the tube (periodic)

This version uses gmsh's built-in 3D meshing capabilities to generate
a volumetric tetrahedral mesh from the surface definition.
"""

using Gmsh: gmsh

"""
    torus_parametric(u, v, R, r)

Calculate a point on the torus surface.

# Arguments
- `u`: Angle around major circle (0 to 2π)
- `v`: Angle around tube (0 to 2π) - periodic
- `R`: Major radius
- `r`: Minor radius

# Returns
- Tuple (x, y, z) coordinates
"""
function torus_parametric(u, v, R, r)
    x = (R + r * cos(v)) * cos(u)
    y = (R + r * cos(v)) * sin(u)
    z = r * sin(v)
    return (x, y, z)
end


"""
    generate_torus_boundary_points(u_res=60, v_res=40, R=2.0, r=0.8)

Generate boundary point cloud for a torus from parametric equations.

# Arguments
- `u_res`: Resolution in u direction (around major circle)
- `v_res`: Resolution in v direction (around tube)
- `R`: Major radius (default: 2.0)
- `r`: Minor radius (default: 0.8)

# Returns
- Vector of (x,y,z) tuples representing the boundary point cloud
"""
function generate_torus_boundary_points(u_res=60, v_res=40, R=2.0, r=0.8)
    # Generate parameter space
    # u ∈ [0, 2π] - wraps around (periodic)
    # v ∈ [0, 2π] - wraps around (periodic)
    u_vals = range(0, stop=2π, length=u_res+1)[1:end-1]  # Exclude endpoint
    v_vals = range(0, stop=2π, length=v_res+1)[1:end-1]  # Exclude endpoint
    
    # Generate boundary point cloud
    points = []
    for u in u_vals
        for v in v_vals
            pt = torus_parametric(u, v, R, r)
            push!(points, pt)
        end
    end
    
    println("  Generated $(length(points)) boundary points")
    return points
end


"""
    demo_torus(; u_res=60, v_res=40, R=2.0, r=0.8, mesh_size=0.2, output="demo_torus.nas")

Generate ordinary torus volumetric mesh using gmsh.

Steps:
1. Generate boundary point cloud from parametric equations
2. Add points to gmsh
3. Create structured surface mesh with proper periodic wrapping
4. Use gmsh to fill interior with tetrahedral mesh
5. Export to NAS format

# Arguments
- `u_res`: Resolution in u direction (default: 60)
- `v_res`: Resolution in v direction (default: 40)
- `R`: Major radius (default: 2.0)
- `r`: Minor radius (default: 0.8)
- `mesh_size`: Target mesh element size for gmsh (default: 0.2)
- `output`: Output NAS filename (default: "demo_torus.nas")
"""
function demo_torus(; u_res=60, v_res=40, R=2.0, r=0.8, mesh_size=0.2, output="demo_torus.nas")
    println("=== Ordinary Torus Volumetric Mesh Generation ===")
    println("  Parameters: R=$R, r=$r")
    println("  Resolution: u=$u_res, v=$v_res")
    println("  Target mesh size: $mesh_size")
    println()
    
    # Step 1: Generate boundary point cloud from parametric equations
    println("Step 1: Generating boundary point cloud...")
    points = generate_torus_boundary_points(u_res, v_res, R, r)
    
    # Initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("torus")
    
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
        
        # Step 3: Create structured surface connectivity with curve loops
        # Both u and v are periodic for an ordinary torus
        println()
        println("Step 3: Creating surface connectivity...")
        
        # Helper function to get point index from (i,j)
        # Both directions wrap periodically
        function point_index(i, j)
            return ((i % u_res) * v_res + (j % v_res)) + 1
        end
        
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
        # Using addSurfaceFilling instead of addPlaneSurface for curved surfaces
        surface_tags = []
        surf_tag = 1
        
        for i in 0:(u_res-1)
            for j in 0:(v_res-1)
                # Get the 4 corners of the quad (all periodic)
                p1 = point_index(i, j)
                p2 = point_index(i+1, j)
                p3 = point_index(i+1, j+1)
                p4 = point_index(i, j+1)
                
                # Triangle 1: p1-p2-p3
                l1 = get_or_create_line(p1, p2)
                l2 = get_or_create_line(p2, p3)
                l3 = get_or_create_line(p3, p1)
                
                loop1 = gmsh.model.geo.addCurveLoop([l1, l2, l3], surf_tag)
                s1 = gmsh.model.geo.addSurfaceFilling([loop1], surf_tag)
                push!(surface_tags, s1)
                surf_tag += 1
                
                # Triangle 2: p1-p3-p4
                l4 = get_or_create_line(p1, p3)
                l5 = get_or_create_line(p3, p4)
                l6 = get_or_create_line(p4, p1)
                
                loop2 = gmsh.model.geo.addCurveLoop([l4, l5, l6], surf_tag)
                s2 = gmsh.model.geo.addSurfaceFilling([loop2], surf_tag)
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
        println("  Torus parameters: R=$R, r=$r")
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
    demo_torus()
end
