

"""
    check_surface_closure(filename)

Check if the surface mesh is properly closed (no boundary edges).

Returns a named tuple with:
- is_closed: boolean indicating if surface is closed
- boundary_edge_count: number of boundary edges found
- surface_count: number of surface entities
- status: :ok if closed, :error if open
"""
function check_surface_closure(filename)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    
    result = (is_closed=false, boundary_edge_count=-1, surface_count=0, status=:error)
    
    try
        gmsh.open(filename)
        
        # Create edges and check for boundary
        gmsh.model.mesh.createEdges()
        
        # Get all surfaces
        surfaces = gmsh.model.getEntities(2)
        
        if length(surfaces) > 0
            boundaries = gmsh.model.getBoundary(surfaces, false, false, true)
            
            is_closed = length(boundaries) == 0
            status = is_closed ? :ok : :error
            
            result = (
                is_closed=is_closed,
                boundary_edge_count=length(boundaries),
                surface_count=length(surfaces),
                status=status
            )
        else
            # No explicit surfaces - might be volume-only mesh
            result = (is_closed=true, boundary_edge_count=0, 
                     surface_count=0, status=:ok)
        end
        
    finally
        gmsh.finalize()
    end
    
    return result
end
