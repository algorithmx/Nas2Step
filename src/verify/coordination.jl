
"""
    check_vertex_coordination(filename; min_coord=3, overcoord_percentile=95)

Check vertex coordination numbers (number of elements connected to each vertex).
Identify anomalous vertices with too few or too many connections.

Arguments:
- `min_coord`: minimum acceptable coordination number (default: 3)
- `overcoord_percentile`: percentile threshold for overcoordination (default: 90)
  Vertices with coordination > this percentile are flagged as overcoordinated

Returns a named tuple with:
- total_vertices: total number of unique vertices
- coordination_numbers: Dict mapping vertex_id => coordination_number
- undercoordinated: array of vertices with coordination < min_coord
- overcoordinated: array of vertices with coordination > percentile threshold
- coord_distribution: Dict mapping coordination_number => count
- mean_coordination: average coordination number
- median_coordination: median coordination number
- p90_coordination: 90th percentile coordination
- max_coord_threshold: computed threshold for overcoordination
- status: :ok, :warning, or :error
"""
function check_vertex_coordination(filename; min_coord::Int=3, overcoord_percentile::Real=95)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    
    result = (total_vertices=0, coordination_numbers=Dict{Int,Int}(),
             undercoordinated=[], overcoordinated=[],
             coord_distribution=Dict{Int,Int}(),
             mean_coordination=NaN, median_coordination=NaN,
             status=:error)
    
    try
        gmsh.open(filename)
        
        # Get all tetrahedral elements
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(3, -1)
        
        if length(elem_types) == 0
            return result
        end
        
        # Build vertex coordination map: vertex_id => count of connected elements
        vertex_coord = Dict{Int,Int}()
        
        for (etype, etags, node_tags) in zip(elem_types, elem_tags, elem_node_tags)
            if etype == 4  # Tetrahedra
                nodes_per_tet = 4
                num_tets = length(etags)
                
                # Count element connections for each vertex
                for i in 1:num_tets
                    base = (i-1) * nodes_per_tet
                    for j in 1:nodes_per_tet
                        node_id = Int(node_tags[base + j])
                        vertex_coord[node_id] = get(vertex_coord, node_id, 0) + 1
                    end
                end
            end
        end
        
        # Analyze coordination numbers
        total_vertices = length(vertex_coord)
        
        if total_vertices == 0
            return result
        end
        
        # Find undercoordinated vertices
        undercoordinated = Tuple{Int,Int}[]  # (vertex_id, coord_num)
        overcoordinated = Tuple{Int,Int}[]
        
        for (vid, coord) in vertex_coord
            if coord < min_coord
                push!(undercoordinated, (vid, coord))
            end
        end
        
        # Sort undercoordinated by coordination number (most problematic first)
        sort!(undercoordinated, by = x -> x[2])
        
        # Build coordination distribution histogram
        coord_distribution = Dict{Int,Int}()
        for coord in values(vertex_coord)
            coord_distribution[coord] = get(coord_distribution, coord, 0) + 1
        end
        
        # Calculate statistics
        coord_values = collect(values(vertex_coord))
        mean_coord = mean(coord_values)
        median_coord = median(coord_values)
        
        # Calculate percentile threshold for overcoordination
        # quantile is already available from Statistics module imported at top
        max_coord_threshold = quantile(coord_values, overcoord_percentile / 100.0)
        p90_coord = quantile(coord_values, 0.90)
        p95_coord = quantile(coord_values, 0.95)
        p99_coord = quantile(coord_values, 0.99)
        
        # Find overcoordinated vertices based on percentile
        for (vid, coord) in vertex_coord
            if coord > max_coord_threshold && coord >= min_coord
                push!(overcoordinated, (vid, coord))
            end
        end
        
        # Sort by coordination number (most problematic first)
        sort!(overcoordinated, by = x -> x[2], rev=true)
        
        # Determine status
        status = if !isempty(undercoordinated)
            :error  # Undercoordinated vertices are serious
        elseif !isempty(overcoordinated)
            :warning  # Overcoordinated might be acceptable
        else
            :ok
        end
        
        result = (
            total_vertices=total_vertices,
            coordination_numbers=vertex_coord,
            undercoordinated=undercoordinated,
            overcoordinated=overcoordinated,
            coord_distribution=coord_distribution,
            mean_coordination=mean_coord,
            median_coordination=median_coord,
            p90_coordination=p90_coord,
            p95_coordination=p95_coord,
            p99_coordination=p99_coord,
            min_coord=min_coord,
            max_coord_threshold=ceil(Int, max_coord_threshold),
            overcoord_percentile=overcoord_percentile,
            status=status
        )
        
    finally
        gmsh.finalize()
    end
    
    return result
end
