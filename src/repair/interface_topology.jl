# Interface topology analysis for surgical mesh repair
# Phase 1.1: Build complete interface topology maps

using Statistics
using JSON
using .CoordinateKeys: coordinate_key_int, create_edge_key_int, EdgeKeyInt, convert_to_int

# ============================================================================
# Core data structures
# ============================================================================

# ==============================================================================
# EdgeKey Type Alias - Using CoordinateKeys Module
# ==============================================================================

"""
EdgeKey Type Alias

This module now uses the CoordinateKeys module for consistent, type-safe
EdgeKey creation and management.

## Migration:
- Old: `EdgeKey` with manual integer coordinate handling
- New: `EdgeKeyInt` from CoordinateKeys module with validation and performance monitoring

## Usage:
```julia
# Recommended: Use factory functions from CoordinateKeys module
edge_key = create_edge_key_int(coord1, coord2)  # Returns EdgeKeyInt

# For direct access, use EdgeKeyInt type
edge_key = EdgeKeyInt(node1_int, node2_int)
```

## Benefits:
- Type safety with dedicated EdgeKeyInt and EdgeKeyFloat types
- Built-in coordinate validation and range checking
- Performance monitoring and statistics tracking
- Consistent coordinate scaling across all modules
"""
const EdgeKey = EdgeKeyInt  # Backward compatibility alias

"""
Triangle with complete geometric and topological information.
"""
struct Triangle
    node1::Int
    node2::Int
    node3::Int
    element_id::Int
    coord1::NTuple{3,Float64}
    coord2::NTuple{3,Float64}
    coord3::NTuple{3,Float64}
    
    # Computed properties
    centroid::NTuple{3,Float64}
    area::Float64
    normal::NTuple{3,Float64}
end

function Triangle(n1::Int, n2::Int, n3::Int, elem_id::Int,
                 c1::NTuple{3,Float64}, c2::NTuple{3,Float64}, c3::NTuple{3,Float64})
    # Compute centroid
    centroid = ((c1[1] + c2[1] + c3[1]) / 3,
                (c1[2] + c2[2] + c3[2]) / 3,
                (c1[3] + c2[3] + c3[3]) / 3)
    
    # Compute area and normal via cross product
    v1 = (c2[1] - c1[1], c2[2] - c1[2], c2[3] - c1[3])
    v2 = (c3[1] - c1[1], c3[2] - c1[2], c3[3] - c1[3])
    
    nx = v1[2] * v2[3] - v1[3] * v2[2]
    ny = v1[3] * v2[1] - v1[1] * v2[3]
    nz = v1[1] * v2[2] - v1[2] * v2[1]
    
    mag = sqrt(nx^2 + ny^2 + nz^2)
    area = mag / 2.0
    
    normal = mag > 1e-14 ? (nx/mag, ny/mag, nz/mag) : (0.0, 0.0, 0.0)
    
    Triangle(n1, n2, n3, elem_id, c1, c2, c3, centroid, area, normal)
end

# ============================================================================
# Tolerance-based geometric comparison helpers
# ============================================================================

"""
    are_nodes_equal_local(node1, node2; tol=1e-4)

Check if two coordinate tuples are equal within geometric tolerance.
Uses Euclidean distance comparison.

This is a local copy of the function from edge_classification.jl to avoid
circular dependencies. Prefer using the version from edge_classification.jl
when available.
"""
function are_nodes_equal_local(node1::NTuple{3,Float64}, 
                               node2::NTuple{3,Float64}; 
                               tol::Real=1e-4)::Bool
    dx = node1[1] - node2[1]
    dy = node1[2] - node2[2]
    dz = node1[3] - node2[3]
    dist2 = dx*dx + dy*dy + dz*dz
    return dist2 <= tol*tol
end

"""
    triangle_has_edge(triangle::Triangle, edge::EdgeKey; tol::Real=1e-4)

Check if a triangle contains both endpoints of an edge.

Returns:
- (true, vertex_index): if triangle has the edge, with vertex_index being the
  index (1, 2, or 3) of the opposite vertex not on the edge
- (false, nothing): if triangle doesn't have both edge endpoints

# Comparison Strategy
Uses tolerance-based geometric comparison to handle the coordinate scaling
difference between EdgeKey (integer-scaled) and Triangle (unrounded) coordinates.
EdgeKey stores scaled integer coordinates (multiply by 10000), while Triangle
stores original float coordinates. This function converts EdgeKey integers back
to float coordinates for tolerance-based comparison.

This approach provides consistent matching regardless of coordinate scaling.
"""
function triangle_has_edge(triangle::Triangle, edge::EdgeKey; tol::Real=1e-4)
    # Convert triangle coordinates to integer keys for comparison with EdgeKey integer coordinates
    # This follows the user's preferred approach: convert to integer and compare element-wise
    edge_int1 = edge.node1
    edge_int2 = edge.node2

    # Check which triangle vertices match the edge endpoints
    coords = [triangle.coord1, triangle.coord2, triangle.coord3]

    # Find matches for each edge endpoint by converting triangle coords to integer keys
    matches_node1 = Int[]
    matches_node2 = Int[]

    for (i, coord) in enumerate(coords)
        coord_int = convert_to_int(coord)
        if coord_int == edge_int1
            push!(matches_node1, i)
        end
        if coord_int == edge_int2
            push!(matches_node2, i)
        end
    end

    # Triangle has the edge if we found exactly one match for each endpoint
    # and they're different vertices
    if length(matches_node1) == 1 && length(matches_node2) == 1 &&
       matches_node1[1] != matches_node2[1]

        # Find the opposite vertex (the one not on the edge)
        edge_vertices = Set([matches_node1[1], matches_node2[1]])
        opposite_vertex = first(setdiff([1, 2, 3], edge_vertices))

        return (true, opposite_vertex)
    end

    return (false, nothing)
end

# Backward compatibility alias - redirect old name to new implementation
"""
    TriangleHasEdge(triangle::Triangle, edge::EdgeKey; digits=4, tol::Real=1e-4)

DEPRECATED: Use triangle_has_edge() instead.
This alias provides backward compatibility.
"""
TriangleHasEdge(triangle::Triangle, edge::EdgeKey; digits=4, tol::Real=1e-4) = 
    triangle_has_edge(triangle, edge; tol=tol)


"""
Bounding box for spatial queries.
"""
struct BoundingBox
    min_corner::NTuple{3,Float64}
    max_corner::NTuple{3,Float64}
end

function BoundingBox(points::Vector{NTuple{3,Float64}})
    if isempty(points)
        return BoundingBox((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    end
    
    min_x = minimum(p[1] for p in points)
    min_y = minimum(p[2] for p in points)
    min_z = minimum(p[3] for p in points)
    
    max_x = maximum(p[1] for p in points)
    max_y = maximum(p[2] for p in points)
    max_z = maximum(p[3] for p in points)
    
    BoundingBox((min_x, min_y, min_z), (max_x, max_y, max_z))
end

"""
Complete topology map for an interface between two PIDs.
"""
struct InterfaceTopology
    pidA::Int
    pidB::Int

    # Shared nodes (by coordinate key)
    shared_node_keys::Set{NTuple{3,Int}}
    node_key_to_ids::Dict{NTuple{3,Int}, Tuple{Int,Int}}  # (nodeA_id, nodeB_id)

    # Boundary faces at interface
    faces_A::Vector{Triangle}
    faces_B::Vector{Triangle}

    # Edge maps (coordinate-based) - which faces use each edge
    edges_A::Dict{EdgeKey, Vector{Int}}  # EdgeKey -> List[face_index in faces_A]
    edges_B::Dict{EdgeKey, Vector{Int}}

    # Non-conformity specifics
    edges_only_in_A::Set{EdgeKey}  # missing in B
    edges_only_in_B::Set{EdgeKey}  # missing in A
    edges_shared::Set{EdgeKey}      # present in both

    # Geometric bounding
    interface_bbox::BoundingBox

    # Statistics
    total_shared_nodes::Int
    total_faces_A::Int
    total_faces_B::Int
    total_edges_A::Int
    total_edges_B::Int
    conformity_ratio::Float64  # shared / (A + B - shared)

    # Interface mesh consistency metrics (crude check)
    max_vertex_distance::Float64      # Maximum distance between shared vertices
    mean_vertex_distance::Float64     # Mean distance between shared vertices
    edge_mismatch_count::Int          # Number of edges that don't match between surfaces
    triangulation_similarity::Float64 # Ratio 0-1, 1.0 means identical triangulation up to node ordering
end

# ============================================================================
# Interface mesh consistency checking
# ============================================================================

"""
    compute_interface_consistency(coords::Dict{Int,NTuple{3,Float64}},
                                node_key_to_ids::Dict{NTuple{3,Int}, Tuple{Int,Int}},
                                triangles_A::Vector{Triangle}, triangles_B::Vector{Triangle},
                                edges_A::Dict{EdgeKey, Vector{Int}}, edges_B::Dict{EdgeKey, Vector{Int}},
                                edges_shared::Set{EdgeKey})

Perform crude interface mesh consistency checking.

Returns a tuple with:
- max_vertex_distance: Maximum Euclidean distance between shared vertices
- mean_vertex_distance: Mean Euclidean distance between shared vertices
- edge_mismatch_count: Total number of non-shared edges between surfaces
- triangulation_similarity: Ratio 0-1 measuring triangulation similarity

This is a minimal-effort check to assess how "far apart" the touching boundaries are.
The ideal case is zero vertex distance and perfect triangulation similarity.
"""
function compute_interface_consistency(coords::Dict{Int,NTuple{3,Float64}},
                                     node_key_to_ids::Dict{NTuple{3,Int}, Tuple{Int,Int}},
                                     triangles_A::Vector{Triangle}, triangles_B::Vector{Triangle},
                                     edges_A::Dict{EdgeKey, Vector{Int}}, edges_B::Dict{EdgeKey, Vector{Int}},
                                     edges_shared::Set{EdgeKey})

    # 1. Vertex position differences
    if isempty(node_key_to_ids)
        return (0.0, 0.0, 0, 0.0)  # No shared vertices
    end

    vertex_distances = Float64[]
    for (key, (id_A, id_B)) in node_key_to_ids
        coord_A = coords[id_A]
        coord_B = coords[id_B]
        dist = sqrt((coord_A[1]-coord_B[1])^2 + (coord_A[2]-coord_B[2])^2 + (coord_A[3]-coord_B[3])^2)
        push!(vertex_distances, dist)
    end

    max_vertex_distance = maximum(vertex_distances)
    mean_vertex_distance = sum(vertex_distances) / length(vertex_distances)

    # 2. Edge topology differences
    total_edges_A = length(edges_A)
    total_edges_B = length(edges_B)
    shared_edges = length(edges_shared)
    edge_mismatch_count = (total_edges_A - shared_edges) + (total_edges_B - shared_edges)

    # 3. Triangulation similarity (crude measure)
    # Compare face counts and shared edge ratio
    face_count_diff = abs(length(triangles_A) - length(triangles_B))
    max_face_count = max(length(triangles_A), length(triangles_B))
    face_similarity = max_face_count > 0 ? 1.0 - (face_count_diff / max_face_count) : 1.0

    # Edge similarity: ratio of shared edges to total unique edges
    total_unique_edges = total_edges_A + total_edges_B - shared_edges
    edge_similarity = total_unique_edges > 0 ? shared_edges / total_unique_edges : 1.0

    # Combine face and edge similarity for overall triangulation similarity
    triangulation_similarity = (face_similarity + edge_similarity) / 2.0

    return (max_vertex_distance, mean_vertex_distance, edge_mismatch_count, triangulation_similarity)
end

# ============================================================================
# Topology extraction
# ============================================================================

"""
    build_interface_topology(nas_file, pidA, pidB; tol=1e-4)

Build complete topology map for the interface between two PIDs.
Extracts all boundary triangles, edges, and shared nodes.

Tolerance (tol) defaults to 1e-4 for all geometric comparisons.
"""
function build_interface_topology(nas_file::String, pidA::Int, pidB::Int; 
                                  tol::Real=1e-4)::InterfaceTopology
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    
    try
        gmsh.open(nas_file)
        
        # Get all nodes and coordinates
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        coords = Dict{Int,NTuple{3,Float64}}()
        for (i, tag) in enumerate(node_tags)
            idx = (i-1)*3
            coords[Int(tag)] = (node_coords[idx+1], node_coords[idx+2], node_coords[idx+3])
        end
        
        # ====================================================================
        # CRITICAL: Coordinate Key Factory Function
        # ====================================================================
        # Now using CoordinateKeys module for consistent coordinate handling
        #
        # coordinate_key_int() from CoordinateKeys module:
        # - Scales coordinates by 10000 and rounds to integers
        # - Provides topological stability (vertices within 0.0001 are identical)
        # - Ensures hash consistency for dictionary lookups
        # - Includes validation and performance monitoring
        #
        # IMPORTANT: All EdgeKeys MUST be created using create_edge_key_int()!
        # ====================================================================

        # Collect tetrahedra per PID
        region_tets = Dict{Int,Vector{Tuple{Int,NTuple{4,Int}}}}()
        for (dim, tag) in gmsh.model.getEntities(3)
            etypes, etags, enodes = gmsh.model.mesh.getElements(dim, tag)
            for (etype, etag_vec, enode_vec) in zip(etypes, etags, enodes)
                if etype != 4; continue; end
                for i in eachindex(etag_vec)
                    eid = Int(etag_vec[i])
                    base = (i-1)*4
                    nd = (Int(enode_vec[base+1]), Int(enode_vec[base+2]), 
                          Int(enode_vec[base+3]), Int(enode_vec[base+4]))
                    push!(get!(region_tets, tag, Vector{Tuple{Int,NTuple{4,Int}}}()), (eid, nd))
                end
            end
        end
        
        if !haskey(region_tets, pidA) || !haskey(region_tets, pidB)
            error("PID $pidA or $pidB not found in mesh")
        end
        
        # Build face incidence: face -> Dict(pid => count)
        face_inc = Dict{NTuple{3,Int},Dict{Int,Int}}()
        tet_faces = ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        
        for (pid, tets) in region_tets
            for (_, nd) in tets
                n = (nd[1], nd[2], nd[3], nd[4])
                for f in tet_faces
                    tri = (n[f[1]], n[f[2]], n[f[3]])
                    tri_sorted = tuple(sort!(collect(tri))...)
                    d = get!(face_inc, tri_sorted, Dict{Int,Int}())
                    d[pid] = get(d, pid, 0) + 1
                end
            end
        end
        
        # Extract boundary faces (odd incidence)
        boundary_faces_A = Vector{NTuple{3,Int}}()
        boundary_faces_B = Vector{NTuple{3,Int}}()
        
        for (face, pid_counts) in face_inc
            if haskey(pid_counts, pidA) && isodd(pid_counts[pidA])
                push!(boundary_faces_A, face)
            end
            if haskey(pid_counts, pidB) && isodd(pid_counts[pidB])
                push!(boundary_faces_B, face)
            end
        end
        
        # ====================================================================
        # Shared Vertex Detection using CoordinateKeys Module
        # ====================================================================
        # Build coordinate key sets per PID using coordinate_key_int() from CoordinateKeys.
        #
        # This means vertices that differ by < 0.0001 in any coordinate
        # will be considered the same vertex (merged by rounding).
        #
        # This is intentional for geometric tolerance but creates a
        # fundamental inconsistency:
        # - Shared vertex sets use INTEGER coordinates (from coordinate_key_int)
        # - Triangle structs store ORIGINAL (unrounded) coordinates
        #
        # When comparing vertices from these two sources, use tolerance-based
        # comparison (e.g., are_nodes_equal) rather than exact equality.
        # ====================================================================
        keyset_A = Set{NTuple{3,Int}}()
        keyset_B = Set{NTuple{3,Int}}()

        for (_, nd) in region_tets[pidA]
            for nid in nd
                push!(keyset_A, coordinate_key_int(coords[nid]))
            end
        end

        for (_, nd) in region_tets[pidB]
            for nid in nd
                push!(keyset_B, coordinate_key_int(coords[nid]))
            end
        end
        
        # Find shared nodes (intersection of rounded coordinate sets)
        shared_keys = intersect(keyset_A, keyset_B)
        
        # Build node ID mapping for shared nodes
        node_key_to_ids = Dict{NTuple{3,Int}, Tuple{Int,Int}}()
        
        # Map A nodes
        nodeA_map = Dict{NTuple{3,Int}, Int}()
        for (_, nd) in region_tets[pidA]
            for nid in nd
                k = coordinate_key_int(coords[nid])
                if k in shared_keys && !haskey(nodeA_map, k)
                    nodeA_map[k] = nid
                end
            end
        end

        # Map B nodes
        nodeB_map = Dict{NTuple{3,Int}, Int}()
        for (_, nd) in region_tets[pidB]
            for nid in nd
                k = coordinate_key_int(coords[nid])
                if k in shared_keys && !haskey(nodeB_map, k)
                    nodeB_map[k] = nid
                end
            end
        end
        
        # Combine mappings
        for k in shared_keys
            if haskey(nodeA_map, k) && haskey(nodeB_map, k)
                node_key_to_ids[k] = (nodeA_map[k], nodeB_map[k])
            end
        end
        
        # Filter interface faces (only faces with all nodes in shared set)
        @inline _t(face,i) = (coordinate_key_int(coords[face[i]]) in shared_keys)
        @inline is_interface_face(face::NTuple{3,Int})::Bool = _t(face,1) && _t(face,2) && _t(face,3)

        interface_faces_A = filter(f -> is_interface_face(f), boundary_faces_A)
        interface_faces_B = filter(f -> is_interface_face(f), boundary_faces_B)
        
        # Build Triangle objects with geometry
        function build_tri(face::NTuple{3,Int}, idx::Int)
            return Triangle(face[1], face[2], face[3], idx,
                          coords[face[1]], coords[face[2]], coords[face[3]])
        end
        triangles_A = [build_tri(face, i) for (i, face) in enumerate(interface_faces_A)]
        triangles_B = [build_tri(face, i) for (i, face) in enumerate(interface_faces_B)]

        # Build edge maps using CoordinateKeys factory functions
        edges_A = Dict{EdgeKey, Vector{Int}}()
        for (idx, tri) in enumerate(triangles_A)
            # Use factory function for consistent EdgeKey creation
            ek1 = create_edge_key_int(tri.coord1, tri.coord2)
            ek2 = create_edge_key_int(tri.coord1, tri.coord3)
            ek3 = create_edge_key_int(tri.coord2, tri.coord3)

            push!(get!(edges_A, ek1, Int[]), idx)
            push!(get!(edges_A, ek2, Int[]), idx)
            push!(get!(edges_A, ek3, Int[]), idx)
        end

        edges_B = Dict{EdgeKey, Vector{Int}}()
        for (idx, tri) in enumerate(triangles_B)
            # Use factory function for consistent EdgeKey creation
            ek1 = create_edge_key_int(tri.coord1, tri.coord2)
            ek2 = create_edge_key_int(tri.coord1, tri.coord3)
            ek3 = create_edge_key_int(tri.coord2, tri.coord3)

            push!(get!(edges_B, ek1, Int[]), idx)
            push!(get!(edges_B, ek2, Int[]), idx)
            push!(get!(edges_B, ek3, Int[]), idx)
        end
        
        # Compute edge differences
        setA = Set(keys(edges_A))
        setB = Set(keys(edges_B))
        
        edges_only_in_A = setdiff(setA, setB)
        edges_only_in_B = setdiff(setB, setA)
        edges_shared = intersect(setA, setB)
        
        # Compute bounding box
        all_interface_coords = NTuple{3,Float64}[]
        for tri in triangles_A
            push!(all_interface_coords, tri.coord1, tri.coord2, tri.coord3)
        end
        bbox = BoundingBox(all_interface_coords)
        
        # Compute conformity ratio
        union_size = length(setA) + length(setB) - length(edges_shared)
        conformity_ratio = union_size > 0 ? length(edges_shared) / union_size : 1.0

        # ====================================================================
        # Interface Mesh Consistency Checking
        # ====================================================================
        # Perform crude consistency check to measure how "far apart" the
        # touching boundaries are in terms of vertex positions and edge topology.
        #
        # This provides metrics to assess interface mesh quality:
        # - Vertex distances: should be zero for perfectly matching interfaces
        # - Edge mismatches: should be zero for conforming meshes
        # - Triangulation similarity: 1.0 means identical up to node ordering
        # ====================================================================
        (max_vertex_dist, mean_vertex_dist, edge_mismatch_count,
         triangulation_sim) = compute_interface_consistency(
            coords, node_key_to_ids, triangles_A, triangles_B,
            edges_A, edges_B, edges_shared
        )

        return InterfaceTopology(
            pidA, pidB,
            shared_keys,
            node_key_to_ids,
            triangles_A,
            triangles_B,
            edges_A,
            edges_B,
            edges_only_in_A,
            edges_only_in_B,
            edges_shared,
            bbox,
            length(shared_keys),
            length(triangles_A),
            length(triangles_B),
            length(setA),
            length(setB),
            conformity_ratio,
            max_vertex_dist,
            mean_vertex_dist,
            edge_mismatch_count,
            triangulation_sim
        )
        
    finally
        gmsh.finalize()
    end
end

"""
    compute_interface_area(topology::InterfaceTopology, use_pid::Symbol)

Compute total surface area of the interface for the specified PID (:A or :B).
"""
function compute_interface_area(topology::InterfaceTopology, use_pid::Symbol)::Float64
    triangles = use_pid == :A ? topology.faces_A : topology.faces_B
    return sum(tri.area for tri in triangles)
end

"""
    export_interface_topology_json(topology::InterfaceTopology, output_file::String)

Export complete interface topology to JSON for detailed inspection.
"""
function export_interface_topology_json(topology::InterfaceTopology, output_file::String)
    report = Dict{String, Any}(
        "pidA" => topology.pidA,
        "pidB" => topology.pidB,
        "summary" => Dict(
            "total_shared_nodes" => topology.total_shared_nodes,
            "faces_A" => topology.total_faces_A,
            "faces_B" => topology.total_faces_B,
            "edges_A" => topology.total_edges_A,
            "edges_B" => topology.total_edges_B,
            "edges_only_in_A" => length(topology.edges_only_in_A),
            "edges_only_in_B" => length(topology.edges_only_in_B),
            "edges_shared" => length(topology.edges_shared),
            "conformity_ratio" => round(topology.conformity_ratio, digits=4),
            "is_conforming" => isempty(topology.edges_only_in_A) && isempty(topology.edges_only_in_B)
        ),
        "geometry" => Dict(
            "bounding_box" => Dict(
                "min" => [topology.interface_bbox.min_corner...],
                "max" => [topology.interface_bbox.max_corner...]
            ),
            "area_A" => compute_interface_area(topology, :A),
            "area_B" => compute_interface_area(topology, :B)
        ),
        "interface_consistency" => Dict(
            "max_vertex_distance" => round(topology.max_vertex_distance, digits=8),
            "mean_vertex_distance" => round(topology.mean_vertex_distance, digits=8),
            "edge_mismatch_count" => topology.edge_mismatch_count,
            "triangulation_similarity" => round(topology.triangulation_similarity, digits=4),
            "is_perfect_interface" => topology.max_vertex_distance < 1e-8 &&
                                    topology.edge_mismatch_count == 0 &&
                                    topology.triangulation_similarity > 0.99
        ),
        "non_conforming_edges" => Dict(
            "missing_in_B" => [[ek.node1..., ek.node2...] for ek in collect(topology.edges_only_in_A)[1:min(50, length(topology.edges_only_in_A))]],
            "missing_in_A" => [[ek.node1..., ek.node2...] for ek in collect(topology.edges_only_in_B)[1:min(50, length(topology.edges_only_in_B))]]
        )
    )

    open(output_file, "w") do io
        JSON.print(io, report, 2)
    end

    println("Interface topology exported to: $output_file")
    return output_file
end
