
"""
    check_element_volumes(filename; swap_nodes=nothing)

Check for inverted elements (negative volumes) which indicate mesh twists or inversions.

Optional keyword:
- swap_nodes::Union{Nothing,Tuple{Int,Int}}: if provided (e.g., (1,2)), compute volumes as if
  the given pair of local node positions were swapped for every element. This is useful to test
  node ordering convention mismatches without modifying the mesh.

Returns a named tuple with:
- total_elements: total number of elements
- inverted_count: number of inverted elements
- volumes: array of all element volumes
- inverted_elements: array of inverted element info (element_id, volume, node_ids)
- volume_stats: volume statistics (min, max, mean, median)
- min_volume: most negative volume (worst inversion)
- max_volume: largest positive volume
- zero_volume_count: number of degenerate (zero volume) elements
- convention_mismatch: true if all/most elements inverted (convention issue, not mesh defect)
- status: :ok, :warning, :error, or :convention_mismatch
"""
function check_element_volumes(filename; swap_nodes::Union{Nothing,Tuple{Int,Int}}=nothing)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    
    result = (total_elements=0, inverted_count=0, volumes=Float64[], 
             inverted_elements=[], volume_stats=nothing, min_volume=0.0, 
             max_volume=0.0, zero_volume_count=0, convention_mismatch=false, 
             status=:error)
    
    try
        gmsh.open(filename)
        
        # Get all tetrahedral elements
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(3, -1)
        
        if length(elem_types) == 0
            return result
        end
        
        total_tets = 0
        inverted_count = 0
        zero_volume_count = 0
        volumes = Float64[]
        inverted_elements = []
        
        for (etype, etags, node_tags) in zip(elem_types, elem_tags, elem_node_tags)
            if etype == 4  # Tetrahedra
                total_tets = length(etags)
                
                # Get node coordinates
                all_nodes, coords, _ = gmsh.model.mesh.getNodes()
                node_coords = Dict{Int,Tuple{Float64,Float64,Float64}}()
                for i in 1:length(all_nodes)
                    idx = Int(all_nodes[i])
                    node_coords[idx] = (coords[3*(i-1)+1], coords[3*(i-1)+2], coords[3*(i-1)+3])
                end
                
                # Check each tetrahedron
                nodes_per_tet = 4
                for i in 1:total_tets
                    element_id = Int(etags[i])
                    
                    # Determine local node ordering, with optional swap test
                    idxs = [1, 2, 3, 4]
                    if swap_nodes !== nothing
                        a, b = swap_nodes
                        if 1 <= a <= 4 && 1 <= b <= 4
                            idxs[a], idxs[b] = idxs[b], idxs[a]
                        end
                    end
                    
                    # Get the 4 nodes of this tet (potentially swapped order)
                    base = (i-1)*nodes_per_tet
                    n1 = Int(node_tags[base + idxs[1]])
                    n2 = Int(node_tags[base + idxs[2]])
                    n3 = Int(node_tags[base + idxs[3]])
                    n4 = Int(node_tags[base + idxs[4]])
                    
                    # Get coordinates
                    p1 = node_coords[n1]
                    p2 = node_coords[n2]
                    p3 = node_coords[n3]
                    p4 = node_coords[n4]
                    
                    # Calculate signed volume: V = 1/6 * det([p2-p1, p3-p1, p4-p1])
                    v1 = (p2[1]-p1[1], p2[2]-p1[2], p2[3]-p1[3])
                    v2 = (p3[1]-p1[1], p3[2]-p1[2], p3[3]-p1[3])
                    v3 = (p4[1]-p1[1], p4[2]-p1[2], p4[3]-p1[3])
                    
                    # Determinant
                    det = v1[1]*(v2[2]*v3[3] - v2[3]*v3[2]) -
                          v1[2]*(v2[1]*v3[3] - v2[3]*v3[1]) +
                          v1[3]*(v2[1]*v3[2] - v2[2]*v3[1])
                    
                    volume = det / 6.0
                    push!(volumes, volume)
                    
                    # Track inverted elements
                    if volume < 0
                        inverted_count += 1
                        
                        # Calculate centroid for spatial information
                        centroid = ((p1[1]+p2[1]+p3[1]+p4[1])/4, 
                                   (p1[2]+p2[2]+p3[2]+p4[2])/4,
                                   (p1[3]+p2[3]+p3[3]+p4[3])/4)
                        
                        # Calculate edge lengths for size information
                        edge_lengths = [
                            sqrt((p2[1]-p1[1])^2 + (p2[2]-p1[2])^2 + (p2[3]-p1[3])^2),
                            sqrt((p3[1]-p1[1])^2 + (p3[2]-p1[2])^2 + (p3[3]-p1[3])^2),
                            sqrt((p4[1]-p1[1])^2 + (p4[2]-p1[2])^2 + (p4[3]-p1[3])^2),
                            sqrt((p3[1]-p2[1])^2 + (p3[2]-p2[2])^2 + (p3[3]-p2[3])^2),
                            sqrt((p4[1]-p2[1])^2 + (p4[2]-p2[2])^2 + (p4[3]-p2[3])^2),
                            sqrt((p4[1]-p3[1])^2 + (p4[2]-p3[2])^2 + (p4[3]-p3[3])^2)
                        ]
                        
                        inverted_elem_info = (
                            element_id=element_id,
                            volume=volume,
                            node_ids=(n1, n2, n3, n4),
                            centroid=centroid,
                            min_edge=minimum(edge_lengths),
                            max_edge=maximum(edge_lengths),
                            coordinates=(p1, p2, p3, p4)
                        )
                        
                        push!(inverted_elements, inverted_elem_info)
                    elseif abs(volume) < 1e-12  # Effectively zero volume
                        zero_volume_count += 1
                    end
                end
            end
        end
        
        # Calculate volume statistics
        volume_stats = if length(volumes) > 0
            (
                min=minimum(volumes),
                max=maximum(volumes), 
                mean=mean(volumes),
                median=median(volumes),
                std=std(volumes)
            )
        else
            nothing
        end
        
        min_volume = length(volumes) > 0 ? minimum(volumes) : 0.0
        max_volume = length(volumes) > 0 ? maximum(volumes) : 0.0
        
        # Check for convention mismatch: if all or nearly all elements are inverted
        # This likely indicates node ordering convention difference, not actual mesh defects
        inversion_ratio = total_tets > 0 ? inverted_count / total_tets : 0.0
        convention_mismatch = inversion_ratio > 0.95  # >95% inverted suggests convention issue
        
        # Determine status
        status = if convention_mismatch
            :convention_mismatch  # All/most elements inverted - convention issue
        elseif inverted_count > 0
            :error  # Some elements inverted - actual mesh defect
        elseif zero_volume_count > 0
            :warning  # Degenerate elements
        else
            :ok
        end
        
        result = (total_elements=total_tets, inverted_count=inverted_count, 
                 volumes=volumes, inverted_elements=inverted_elements,
                 volume_stats=volume_stats, min_volume=min_volume,
                 max_volume=max_volume, zero_volume_count=zero_volume_count,
                 convention_mismatch=convention_mismatch, status=status)
        
    finally
        gmsh.finalize()
    end
    
    return result
end


"""
    check_element_quality(filename; metric="gamma")

Check element quality metrics.

Available metrics:
- "gamma" - radius ratio (0=bad, 1=perfect)
- "eta" - normalized eta quality measure

Returns a named tuple with quality statistics and distribution.
"""
function check_element_quality(filename; metric="gamma")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    
    result = (min=NaN, max=NaN, mean=NaN, median=NaN, 
             poor_count=0, bad_count=0, inverted_count=0,
             qualities=Float64[], status=:error)
    
    try
        gmsh.open(filename)
        
        # Get all tetrahedral elements
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(3, -1)
        
        qualities = Float64[]
        
        for (etype, etags) in zip(elem_types, elem_tags)
            if etype == 4  # Tetrahedra
                for etag in etags
                    quality = gmsh.model.mesh.getElementQualities([Int(etag)], metric)[1]
                    push!(qualities, quality)
                end
            end
        end
        
        if length(qualities) > 0
            # Analyze quality distribution
            poor_count = count(q -> q < 0.1, qualities)
            bad_count = count(q -> q < 0.01, qualities)
            inverted_count = count(q -> q < 0, qualities)
            
            status = if inverted_count > 0
                :error
            elseif bad_count > length(qualities) * 0.05
                :warning
            else
                :ok
            end
            
            result = (
                min=minimum(qualities),
                max=maximum(qualities),
                mean=mean(qualities),
                median=median(qualities),
                poor_count=poor_count,
                bad_count=bad_count,
                inverted_count=inverted_count,
                qualities=qualities,
                status=status
            )
        end
        
    finally
        gmsh.finalize()
    end
    
    return result
end
