# Interface topology analysis for surgical mesh repair
# Phase 1.1: Build complete interface topology maps

using Statistics
using JSON

# ============================================================================
# Core data structures
# ============================================================================

"""
Coordinate-based edge key for geometric matching.
Nodes are ordered (smaller coordinate tuple first) for consistent hashing.
"""
struct EdgeKey
    node1::NTuple{3,Float64}
    node2::NTuple{3,Float64}
    
    function EdgeKey(n1::NTuple{3,Float64}, n2::NTuple{3,Float64})
        # Canonical ordering for consistent comparison
        if n1 <= n2
            new(n1, n2)
        else
            new(n2, n1)
        end
    end
end

Base.hash(ek::EdgeKey, h::UInt) = hash((ek.node1, ek.node2), h)
Base.:(==)(a::EdgeKey, b::EdgeKey) = (a.node1 == b.node1) && (a.node2 == b.node2)

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
    shared_node_keys::Set{NTuple{3,Float64}}
    node_key_to_ids::Dict{NTuple{3,Float64}, Tuple{Int,Int}}  # (nodeA_id, nodeB_id)
    
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
end

# ============================================================================
# Topology extraction
# ============================================================================

"""
    build_interface_topology(nas_file, pidA, pidB; tol=1e-8)

Build complete topology map for the interface between two PIDs.
Extracts all boundary triangles, edges, and shared nodes.
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
        
        # Helper: coordinate key with rounding
        function ckey(p::NTuple{3,Float64})
            return (round(p[1]; digits=4), round(p[2]; digits=4), round(p[3]; digits=4))
        end
        
        # Collect tetrahedra per PID
        region_tets = Dict{Int,Vector{Tuple{Int,NTuple{4,Int}}}}()
        for (dim, tag) in gmsh.model.getEntities(3)
            etypes, etags, enodes = gmsh.model.mesh.getElements(dim, tag)
            for (etype, etag_vec, enode_vec) in zip(etypes, etags, enodes)
                if etype != 4; continue; end
                for i in 1:length(etag_vec)
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
        
        # Build coordinate key sets per PID
        keyset_A = Set{NTuple{3,Float64}}()
        keyset_B = Set{NTuple{3,Float64}}()
        
        for (_, nd) in region_tets[pidA]
            for nid in nd
                push!(keyset_A, ckey(coords[nid]))
            end
        end
        
        for (_, nd) in region_tets[pidB]
            for nid in nd
                push!(keyset_B, ckey(coords[nid]))
            end
        end
        
        # Find shared nodes
        shared_keys = intersect(keyset_A, keyset_B)
        
        # Build node ID mapping for shared nodes
        node_key_to_ids = Dict{NTuple{3,Float64}, Tuple{Int,Int}}()
        
        # Map A nodes
        nodeA_map = Dict{NTuple{3,Float64}, Int}()
        for (_, nd) in region_tets[pidA]
            for nid in nd
                k = ckey(coords[nid])
                if k in shared_keys && !haskey(nodeA_map, k)
                    nodeA_map[k] = nid
                end
            end
        end
        
        # Map B nodes
        nodeB_map = Dict{NTuple{3,Float64}, Int}()
        for (_, nd) in region_tets[pidB]
            for nid in nd
                k = ckey(coords[nid])
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
        function is_interface_face(face::NTuple{3,Int}, shared::Set{NTuple{3,Float64}})
            k1 = ckey(coords[face[1]])
            k2 = ckey(coords[face[2]])
            k3 = ckey(coords[face[3]])
            return (k1 in shared) && (k2 in shared) && (k3 in shared)
        end
        
        interface_faces_A = filter(f -> is_interface_face(f, shared_keys), boundary_faces_A)
        interface_faces_B = filter(f -> is_interface_face(f, shared_keys), boundary_faces_B)
        
        # Build Triangle objects with geometry
        triangles_A = Vector{Triangle}()
        for (idx, face) in enumerate(interface_faces_A)
            n1, n2, n3 = face
            tri = Triangle(n1, n2, n3, idx, coords[n1], coords[n2], coords[n3])
            push!(triangles_A, tri)
        end
        
        triangles_B = Vector{Triangle}()
        for (idx, face) in enumerate(interface_faces_B)
            n1, n2, n3 = face
            tri = Triangle(n1, n2, n3, idx, coords[n1], coords[n2], coords[n3])
            push!(triangles_B, tri)
        end
        
        # Build edge maps
        edges_A = Dict{EdgeKey, Vector{Int}}()
        for (idx, tri) in enumerate(triangles_A)
            k1 = ckey(tri.coord1)
            k2 = ckey(tri.coord2)
            k3 = ckey(tri.coord3)
            
            for (a, b) in ((k1, k2), (k1, k3), (k2, k3))
                ek = EdgeKey(a, b)
                push!(get!(edges_A, ek, Int[]), idx)
            end
        end
        
        edges_B = Dict{EdgeKey, Vector{Int}}()
        for (idx, tri) in enumerate(triangles_B)
            k1 = ckey(tri.coord1)
            k2 = ckey(tri.coord2)
            k3 = ckey(tri.coord3)
            
            for (a, b) in ((k1, k2), (k1, k3), (k2, k3))
                ek = EdgeKey(a, b)
                push!(get!(edges_B, ek, Int[]), idx)
            end
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
            conformity_ratio
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
