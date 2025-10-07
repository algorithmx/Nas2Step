# Geometric Utilities for Interface Analysis
# Shared functions used by multiple interface analysis modules

# ============================================================================
# Node/Vertex comparison utilities
# ============================================================================

"""
    nodes_equal_within_tolerance(node1, node2; tol=1e-4)

Check if two nodes are equal within a given tolerance using Euclidean distance.
"""
function nodes_equal_within_tolerance(node1::NTuple{3,Float64}, 
                                     node2::NTuple{3,Float64}; 
                                     tol::Real=1e-4)::Bool
    dx = node1[1] - node2[1]
    dy = node1[2] - node2[2]
    dz = node1[3] - node2[3]
    dist2 = dx*dx + dy*dy + dz*dz
    return dist2 <= tol*tol
end

"""
    find_matching_vertex(vertex, vertex_set; tol=1e-4)

Find a vertex in vertex_set that matches the given vertex within tolerance.
Returns the matching vertex or nothing.
"""
function find_matching_vertex(vertex::NTuple{3,Float64}, 
                             vertex_set::Set{NTuple{3,Float64}}; 
                             tol::Real=1e-4)::Union{NTuple{3,Float64}, Nothing}
    tol2 = tol * tol
    for v in vertex_set
        dx = vertex[1] - v[1]
        dy = vertex[2] - v[2]
        dz = vertex[3] - v[3]
        dist2 = dx*dx + dy*dy + dz*dz
        if dist2 <= tol2
            return v
        end
    end
    return nothing
end

"""
    is_vertex_in_set(vertex, vertex_set; tol=1e-4)

Check if a vertex exists in a set within the given tolerance.
"""
function is_vertex_in_set(vertex::NTuple{3,Float64}, 
                         vertex_set::Set{NTuple{3,Float64}}; 
                         tol::Real=1e-4)::Bool
    return find_matching_vertex(vertex, vertex_set, tol=tol) !== nothing
end

# ============================================================================
# Vertex extraction and analysis
# ============================================================================

"""
    extract_boundary_vertices(faces::Vector{Triangle})

Extract unique vertices from a set of boundary faces.
"""
function extract_boundary_vertices(faces::Vector{Triangle})::Set{NTuple{3,Float64}}
    vertices = Set{NTuple{3,Float64}}()
    
    for face in faces
        push!(vertices, face.coord1)
        push!(vertices, face.coord2)
        push!(vertices, face.coord3)
    end
    
    return vertices
end

"""
    compute_shared_vertices(vertices_A, vertices_B; tol=1e-4)

Compute which vertices are shared between two vertex sets.
Returns the set of shared vertices (using coordinates from set A as canonical).
"""
function compute_shared_vertices(vertices_A::Set{NTuple{3,Float64}}, 
                                vertices_B::Set{NTuple{3,Float64}}; 
                                tol::Real=1e-4)::Set{NTuple{3,Float64}}
    shared = Set{NTuple{3,Float64}}()
    
    for v_a in vertices_A
        if is_vertex_in_set(v_a, vertices_B, tol=tol)
            push!(shared, v_a)
        end
    end
    
    return shared
end
