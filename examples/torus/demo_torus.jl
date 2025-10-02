#!/usr/bin/env julia

using Gmsh: gmsh
using Nas2Step

# Torus parametric equations: x = (R + r*cos(v)) * cos(u), y = (R + r*cos(v)) * sin(u), z = r * sin(v)
torus_parametric(u, v, R, r) = ((R + r * cos(v)) * cos(u), (R + r * cos(v)) * sin(u), r * sin(v))


function demo_torus(; u_res=24, v_res=12, R=2.0, r=0.8, mesh_size=0.5, output="demo_torus.nas")
    println("Generating torus mesh: R=$R, r=$r, resolution=($u_res, $v_res)")
    
    # Generate boundary points
    u_vals = range(0, 2π, length=u_res+1)[1:end-1]
    v_vals = range(0, 2π, length=v_res+1)[1:end-1]
    boundary_points = [torus_parametric(u, v, R, r) for u in u_vals for v in v_vals]
    
    output_path = "$(@__DIR__)/$output"
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)  # 0=silent, 1=errors, 2=warnings, 3=info
    gmsh.model.add("torus")
    
    try
        # Add boundary points to gmsh
        for (idx, pt) in enumerate(boundary_points)
            gmsh.model.geo.addPoint(pt[1], pt[2], pt[3], mesh_size, idx)
        end
        gmsh.model.geo.synchronize()
        
        # Create structured surface mesh (periodic in both directions)
        point_index(i, j) = ((i % u_res) * v_res + (j % v_res)) + 1
        
        lines = Dict{Tuple{Int,Int},Int}()
        line_tag = 1
        get_or_create_line(p1, p2) = begin
            key = minmax(p1, p2)
            if !haskey(lines, key)
                gmsh.model.geo.addLine(key[1], key[2], line_tag)
                lines[key] = line_tag
                line_tag += 1
            end
            lines[key] * (p1 < p2 ? 1 : -1)
        end
        
        # Create surface triangles from structured grid
        surface_tags = []
        surf_tag = 1
        for i in 0:(u_res-1), j in 0:(v_res-1)
            p1, p2, p3, p4 = point_index(i, j), point_index(i+1, j), point_index(i+1, j+1), point_index(i, j+1)
            
            # Triangle 1: p1-p2-p3
            loop1 = gmsh.model.geo.addCurveLoop([get_or_create_line(p1, p2), get_or_create_line(p2, p3), get_or_create_line(p3, p1)], surf_tag)
            push!(surface_tags, gmsh.model.geo.addSurfaceFilling([loop1], surf_tag))
            surf_tag += 1
            
            # Triangle 2: p1-p3-p4
            loop2 = gmsh.model.geo.addCurveLoop([get_or_create_line(p1, p3), get_or_create_line(p3, p4), get_or_create_line(p4, p1)], surf_tag)
            push!(surface_tags, gmsh.model.geo.addSurfaceFilling([loop2], surf_tag))
            surf_tag += 1
        end
        gmsh.model.geo.synchronize()
        
        # Create volume and generate tetrahedral mesh
        gmsh.model.geo.addVolume([gmsh.model.geo.addSurfaceLoop(surface_tags, 1)], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")
        
        # Extract mesh and write to NAS
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        points = [(node_coords[3i-2], node_coords[3i-1], node_coords[3i]) for i in 1:length(node_tags)]
        
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(3, -1)
        tets = [(enode_tags[4i-3], enode_tags[4i-2], enode_tags[4i-1], enode_tags[4i]) 
                for (etype, etags, enode_tags) in zip(elem_types, elem_tags, elem_node_tags) if etype == 4 
                for i in 1:length(etags)]
        
        Nas2Step.write_nas_volume(output_path, [Nas2Step.VolumeRegion("Torus", points, tets)])
        
        println("Wrote $(length(points)) nodes, $(length(tets)) tetrahedra to $output_path")
    finally
        gmsh.finalize()
    end
    return output_path
end


# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    demo_torus()
end
