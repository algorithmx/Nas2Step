# Mesh quality verification functions for Nas2Step

using Statistics

# ============================================================================
# Helper functions for concise reporting
# ============================================================================

function fmt_pct(n, total)
    return round(100 * n / max(total, 1), digits=1)
end

function fmt_stat(x; digits=4)
    return round(x, sigdigits=digits)
end

function print_volume_stats(vs)
    println("  Volume: [$(fmt_stat(vs.min)), $(fmt_stat(vs.max))], " *
           "mean=$(fmt_stat(vs.mean)), std=$(fmt_stat(vs.std))")
end

function print_inverted_summary(vol_check)
    inv_pct = fmt_pct(vol_check.inverted_count, vol_check.total_elements)
    println("  $(vol_check.inverted_count)/$(vol_check.total_elements) inverted ($(inv_pct)%), " *
           "$(vol_check.zero_volume_count) zero-volume")
end

function print_inverted_details(inverted_elements, max_show=10)
    n_show = min(max_show, length(inverted_elements))
    sorted = sort(inverted_elements, by = e -> e.volume)
    
    for i in 1:n_show
        elem = sorted[i]
        println("  [$(i)] Elem $(elem.element_id): V=$(fmt_stat(elem.volume, digits=5)), " *
               "nodes=$(elem.node_ids), center=$(round.(elem.centroid, digits=2))")
    end
    
    if length(inverted_elements) > n_show
        println("  ... $(length(inverted_elements) - n_show) more")
    end
end

function print_swap_test_result(swap_test)
    println("  🔬 Swap test (1⇄2): $(swap_test.original_inverted) → $(swap_test.swapped_inverted) " *
           "[$(swap_test.fixed_count) fixed, $(fmt_pct(swap_test.improvement_ratio, 1))%]")
    if swap_test.would_fix
        println("  ✅ Would fix! This confirms convention mismatch.")
    else
        println("  ⚠️  Wouldn't fully fix. May have real defects.")
    end
end

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


"""
    test_node_swap_fix(filename; swap_pair=(1,2))

Test if swapping a pair of nodes in all elements would fix a convention mismatch.
This performs a virtual test without modifying the mesh file.

Returns a named tuple with:
- original_inverted: number of inverted elements before swap
- swapped_inverted: number of inverted elements after virtual swap
- improvement_ratio: fraction of elements that became non-inverted
- would_fix: boolean indicating if swap would fix >90% of inversions
- swap_pair: the node pair that was tested
"""
function test_node_swap_fix(filename; swap_pair::Tuple{Int,Int}=(1,2))
    # Check original
    original = check_element_volumes(filename)
    
    # Check with swapped nodes
    swapped = check_element_volumes(filename, swap_nodes=swap_pair)
    
    original_inv = original.inverted_count
    swapped_inv = swapped.inverted_count
    total = original.total_elements
    
    # Calculate improvement
    fixed_count = original_inv - swapped_inv
    improvement_ratio = total > 0 ? fixed_count / total : 0.0
    
    # Would this fix most issues?
    would_fix = (swapped_inv < 0.1 * total)  # Less than 10% still inverted
    
    return (
        original_inverted=original_inv,
        swapped_inverted=swapped_inv,
        fixed_count=fixed_count,
        improvement_ratio=improvement_ratio,
        would_fix=would_fix,
        swap_pair=swap_pair,
        original_ratio=original_inv/total,
        swapped_ratio=swapped_inv/total
    )
end


"""
    export_inverted_elements(inverted_elements, output_file="inverted_elements.txt")

Export detailed information about inverted elements to a text file.
Useful for debugging mesh inversions.
"""
function export_inverted_elements(inverted_elements, output_file="inverted_elements.txt")
    if length(inverted_elements) == 0
        println("No inverted elements to export.")
        return
    end
    
    open(output_file, "w") do io
        println(io, "="^80)
        println(io, "INVERTED ELEMENT REPORT")
        println(io, "Total inverted elements: $(length(inverted_elements))")
        println(io, "="^80)
        println(io)
        
        # Sort by volume (most negative first)
        sorted_elems = sort(inverted_elements, by = e -> e.volume)
        
        for (i, elem) in enumerate(sorted_elems)
            println(io, "[$i] Element ID: $(elem.element_id)")
            println(io, "    Volume: $(elem.volume)")
            println(io, "    Node IDs: $(elem.node_ids)")
            println(io, "    Centroid: $(elem.centroid)")
            println(io, "    Edge length: min=$(elem.min_edge), max=$(elem.max_edge)")
            println(io, "    Node coordinates:")
            for (j, coord) in enumerate(elem.coordinates)
                println(io, "      Node $(elem.node_ids[j]): $coord")
            end
            println(io)
        end
        
        println(io, "="^80)
        println(io, "END OF REPORT")
        println(io, "="^80)
    end
    
    println("Inverted element details exported to: $output_file")
end


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


"""
    verify_nas_mesh(filename; verbose=true)

Comprehensive verification of a NAS mesh file.

Checks for:
1. Inverted/negative volume elements (twists)
2. Element quality metrics
3. Surface closure

Returns a named tuple with all check results and an overall status.
"""
function verify_nas_mesh(filename; verbose=true)
    if !isfile(filename)
        error("File $filename not found!")
    end
    
    if verbose
        println("="^70)
        println("Mesh Quality Verification")
        println("File: $filename")
        println("="^70)
        println()
    end
    
    # Check element volumes for inversions
    if verbose
        println("Checking for inverted elements...")
    end
    vol_check = check_element_volumes(filename)
    swap_test_result = nothing  # Will be populated if convention mismatch detected
    
    if verbose
        println("  Elements: $(vol_check.total_elements)")
        print_inverted_summary(vol_check)
        vol_check.volume_stats !== nothing && print_volume_stats(vol_check.volume_stats)
        
        if vol_check.status == :ok
            println("  ✓ OK")
        elseif vol_check.status == :warning
            println("  ⚠️  Degenerate elements detected")
        elseif vol_check.status == :convention_mismatch
            println("  🔄 CONVENTION MISMATCH ($(fmt_pct(vol_check.inverted_count, vol_check.total_elements))% inverted)")
            swap_test = test_node_swap_fix(filename, swap_pair=(1,2))
            swap_test_result = swap_test
            print_swap_test_result(swap_test)
            println("  Likely: RHR/LHR convention difference. Solution: swap node pairs or use abs(volume).")
        else
            println("  ⚠️  INVERTED ELEMENTS (mesh defects):")
            length(vol_check.inverted_elements) > 0 && print_inverted_details(vol_check.inverted_elements, 5)
        end
        println()
    end
    
    # Check element quality
    verbose && println("Checking element quality...")
    qual_check = check_element_quality(filename)
    
    if verbose
        poor_pct = fmt_pct(qual_check.poor_count, length(qual_check.qualities))
        println("  Quality: mean=$(round(qual_check.mean, digits=3)), " *
               "median=$(round(qual_check.median, digits=3)), " *
               "poor(<0.1)=$(qual_check.poor_count) ($(poor_pct)%)")
        
        status_msg = qual_check.status == :ok ? "✓ OK" :
                    qual_check.status == :warning ? "⚠️  Many poor elements" : "⚠️  Quality issues"
        println("  $status_msg")
        println()
    end
    
    # Check surface closure
    verbose && println("Checking surface closure...")
    closure_check = check_surface_closure(filename)
    
    if verbose
        msg = closure_check.is_closed ? "✓ Closed" : "⚠️  Open ($(closure_check.boundary_edge_count) boundary edges)"  
        println("  $msg")
        println()
    end
    
    # Determine overall status
    overall_status = :ok
    if vol_check.status == :error || qual_check.status == :error || closure_check.status == :error
        overall_status = :error
    elseif vol_check.status == :convention_mismatch
        overall_status = :convention_mismatch
    elseif qual_check.status == :warning || vol_check.status == :warning
        overall_status = :warning
    end
    
    if verbose
        println("="^70)
        println("Overall Status: ", 
                overall_status == :ok ? "✓ PASSED" :
                overall_status == :convention_mismatch ? "🔄 CONVENTION MISMATCH" :
                overall_status == :warning ? "⚠️  WARNING" :
                "❌ FAILED")
        println("="^70)
    end
    
    return (
        volume_check=vol_check,
        quality_check=qual_check,
        closure_check=closure_check,
        swap_test=swap_test_result,
        overall_status=overall_status
    )
end
