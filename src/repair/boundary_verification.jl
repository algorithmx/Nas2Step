"""
    boundary_verification.jl

Boundary consistency verification for unified mesh generation.

Critical for ensuring repaired meshes are suitable for downstream remeshing:
- Boundary loops must be topologically equivalent to originals
- Normal orientations must be consistent
- No flipped or degenerate triangles
"""

using Statistics

# ============================================================================
# Boundary Loop Extraction
# ============================================================================

"""
    extract_boundary_loops(faces::Vector{Triangle}) -> Vector{Vector{NTuple{3,Float64}}}

Extract closed boundary loops from a mesh by finding edges with incidence=1.

A boundary edge is incident to exactly one triangle. This function:
1. Identifies all boundary edges
2. Builds an adjacency graph of boundary nodes
3. Traverses to find closed loops

Returns a vector of loops, where each loop is a vector of rounded node coordinates.

# Example
```julia
loops = extract_boundary_loops(mesh_faces)
# Returns: [[node1, node2, node3, node1], [node4, node5, node6, node4], ...]
```

# Notes
- Uses 4-digit rounding for coordinate keys (matches EdgeKey policy)
- Handles disconnected boundary components (multiple loops)
- Filters out degenerate loops (< 3 nodes)
"""
function extract_boundary_loops(faces::Vector{Triangle})
    # Helper: coordinate rounding (must match EdgeKey policy)
    ckey(p::NTuple{3,Float64}) = (round(p[1]; digits=4), round(p[2]; digits=4), round(p[3]; digits=4))
    
    isempty(faces) && return Vector{Vector{NTuple{3,Float64}}}()
    
    # Step 1: Build edge incidence map
    edge_count = Dict{Tuple{NTuple{3,Float64},NTuple{3,Float64}}, Int}()
    for tri in faces
        k1, k2, k3 = ckey(tri.coord1), ckey(tri.coord2), ckey(tri.coord3)
        for (a, b) in ((k1,k2), (k2,k3), (k3,k1))
            e = a <= b ? (a, b) : (b, a)  # Canonical form
            edge_count[e] = get(edge_count, e, 0) + 1
        end
    end
    
    # Step 2: Find boundary edges (incidence = 1)
    boundary_edges = Tuple{NTuple{3,Float64},NTuple{3,Float64}}[]
    for (e, count) in edge_count
        if count == 1
            push!(boundary_edges, e)
        end
    end
    
    isempty(boundary_edges) && return Vector{Vector{NTuple{3,Float64}}}()
    
    # Step 3: Build adjacency graph (node -> neighbors)
    adj = Dict{NTuple{3,Float64}, Set{NTuple{3,Float64}}}()
    for (a, b) in boundary_edges
        push!(get!(adj, a, Set{NTuple{3,Float64}}()), b)
        push!(get!(adj, b, Set{NTuple{3,Float64}}()), a)
    end
    
    # Step 4: Traverse loops using DFS
    visited_edges = Set{Tuple{NTuple{3,Float64},NTuple{3,Float64}}}()
    loops = Vector{Vector{NTuple{3,Float64}}}()
    
    for start_node in keys(adj)
        neighbors = collect(adj[start_node])
        isempty(neighbors) && continue
        
        # Try each neighbor as potential loop start
        for first_neighbor in neighbors
            edge = start_node <= first_neighbor ? (start_node, first_neighbor) : (first_neighbor, start_node)
            edge in visited_edges && continue
            
            # Trace loop from start_node → first_neighbor
            loop = [start_node]
            prev = start_node
            curr = first_neighbor
            push!(visited_edges, edge)
            
            max_iterations = length(boundary_edges) + 1  # Safety limit
            iterations = 0
            
            while curr != start_node && iterations < max_iterations
                push!(loop, curr)
                
                # Find next unvisited neighbor (not prev)
                next_candidates = setdiff(get(adj, curr, Set{NTuple{3,Float64}}()), [prev])
                isempty(next_candidates) && break
                
                next_node = first(next_candidates)
                edge = curr <= next_node ? (curr, next_node) : (next_node, curr)
                edge in visited_edges && break
                push!(visited_edges, edge)
                
                prev, curr = curr, next_node
                iterations += 1
            end
            
            # Valid closed loop?
            if curr == start_node && length(loop) >= 3
                push!(loops, loop)
            end
        end
    end
    
    return loops
end

# ============================================================================
# Loop Normalization and Comparison
# ============================================================================

"""
    normalize_loop(loop::Vector{NTuple{3,Float64}}) -> Vector{NTuple{3,Float64}}

Normalize a loop to canonical form for comparison.

Canonical form:
1. Start at lexicographically smallest node
2. Traverse in canonical direction (toward smaller next neighbor)

This allows two loops representing the same boundary to compare as equal
even if they start at different nodes or traverse in opposite directions.

# Example
```julia
loop1 = [n2, n3, n1]  # Same boundary as loop2
loop2 = [n1, n2, n3]
normalize_loop(loop1) == normalize_loop(loop2)  # true
```
"""
function normalize_loop(loop::Vector{NTuple{3,Float64}})
    isempty(loop) && return loop
    length(loop) == 1 && return loop
    
    # Find lexicographically minimum node
    min_idx = argmin(loop)
    
    # Rotate to start at minimum
    normalized = vcat(loop[min_idx:end], loop[1:min_idx-1])
    
    # Determine canonical direction
    # Compare next and previous nodes from start
    if length(normalized) >= 2
        next_node = normalized[2]
        prev_node = normalized[end]
        
        # If previous < next, reverse to maintain canonical direction
        if prev_node < next_node
            # Reverse but keep first element in place
            normalized = vcat([normalized[1]], reverse(normalized[2:end]))
        end
    end
    
    return normalized
end

"""
    compare_loop_sets(loops_A, loops_B; tol=1e-6) -> (matched::Int, unmatched_A::Int, unmatched_B::Int)

Compare two sets of boundary loops for topological equivalence.

Loops are normalized before comparison, so they match even if:
- Started at different nodes
- Traversed in opposite directions
- Specified in different order

Returns:
- `matched`: Number of loops that appear in both sets
- `unmatched_A`: Loops in A but not in B
- `unmatched_B`: Loops in B but not in A

# Example
```julia
loops_unified = extract_boundary_loops(unified_mesh)
loops_original = extract_boundary_loops(original_mesh)
(matched, missing, extra) = compare_loop_sets(loops_unified, loops_original)

if missing > 0
    @warn "Unified mesh is missing \$missing boundary loops from original"
end
```
"""
function compare_loop_sets(
    loops_A::Vector{Vector{NTuple{3,Float64}}},
    loops_B::Vector{Vector{NTuple{3,Float64}}};
    tol::Real = 1e-6
)
    # Normalize all loops for comparison
    normalized_A = [normalize_loop(l) for l in loops_A]
    normalized_B = [normalize_loop(l) for l in loops_B]
    
    # Use Set for fast membership testing
    set_A = Set(normalized_A)
    set_B = Set(normalized_B)
    
    # Count matches and mismatches
    matched = length(intersect(set_A, set_B))
    unmatched_A = length(setdiff(set_A, set_B))
    unmatched_B = length(setdiff(set_B, set_A))
    
    return (matched, unmatched_A, unmatched_B)
end

# ============================================================================
# Normal Orientation Verification
# ============================================================================

"""
    compute_triangle_normal(tri::Triangle) -> NTuple{3,Float64}

Compute unit normal vector for a triangle using right-hand rule.

Given vertices v1, v2, v3 in counter-clockwise order when viewed from
the normal direction, computes: normal = (v2-v1) × (v3-v1)

Returns (0,0,0) for degenerate triangles (zero area).

# Example
```julia
normal = compute_triangle_normal(tri)
# Returns: (nx, ny, nz) with ||(nx,ny,nz)|| = 1.0
```
"""
function compute_triangle_normal(tri::Triangle)
    c1, c2, c3 = tri.coord1, tri.coord2, tri.coord3
    
    # Edge vectors from c1
    v1 = (c2[1] - c1[1], c2[2] - c1[2], c2[3] - c1[3])
    v2 = (c3[1] - c1[1], c3[2] - c1[2], c3[3] - c1[3])
    
    # Cross product: v1 × v2
    nx = v1[2] * v2[3] - v1[3] * v2[2]
    ny = v1[3] * v2[1] - v1[1] * v2[3]
    nz = v1[1] * v2[2] - v1[2] * v2[1]
    
    # Normalize to unit length
    mag = sqrt(nx^2 + ny^2 + nz^2)
    if mag < 1e-12
        return (0.0, 0.0, 0.0)  # Degenerate triangle
    end
    
    return (nx/mag, ny/mag, nz/mag)
end

"""
    verify_normal_consistency(faces::Vector{Triangle}; verbose=false) 
        -> (consistent::Bool, flipped::Int, degenerate::Int)

Check if all triangles have consistent normal orientation.

Uses first non-degenerate triangle as reference orientation.
A triangle is "flipped" if its normal points opposite to the reference
(dot product < 0).

This is critical for remeshing algorithms which expect consistent orientation.

# Returns
- `consistent`: true if no triangles are flipped
- `flipped`: count of triangles with opposite orientation
- `degenerate`: count of triangles with zero area

# Example
```julia
(ok, flipped, degen) = verify_normal_consistency(mesh_faces, verbose=true)
if !ok
    @warn "Mesh has \$flipped flipped triangles and \$degen degenerate triangles"
end
```
"""
function verify_normal_consistency(faces::Vector{Triangle}; verbose::Bool=false)
    isempty(faces) && return (true, 0, 0)
    
    # Compute all normals
    normals = [compute_triangle_normal(tri) for tri in faces]
    
    # Count degenerate triangles
    degenerate = count(n -> sum(abs.(n)) < 1e-6, normals)
    
    # Find reference normal (first non-degenerate)
    ref_normal = nothing
    for n in normals
        if sum(abs.(n)) > 1e-6
            ref_normal = n
            break
        end
    end
    
    # All degenerate?
    ref_normal === nothing && return (true, 0, degenerate)
    
    # Count flipped normals (dot product < 0)
    flipped = 0
    for n in normals
        if sum(abs.(n)) > 1e-6  # Skip degenerate
            dot = n[1]*ref_normal[1] + n[2]*ref_normal[2] + n[3]*ref_normal[3]
            if dot < -0.1  # Use threshold to avoid numerical noise
                flipped += 1
            end
        end
    end
    
    consistent = (flipped == 0)
    
    if verbose && !consistent
        println("  ⚠ Normal inconsistency detected:")
        println("    Flipped triangles:    $flipped / $(length(faces))")
        println("    Degenerate triangles: $degenerate / $(length(faces))")
        if flipped > 0
            flip_pct = 100.0 * flipped / length(faces)
            println("    Flip percentage:      $(round(flip_pct, digits=2))%")
        end
    end
    
    return (consistent, flipped, degenerate)
end

"""
    analyze_normal_distribution(faces::Vector{Triangle}) -> Dict{String, Float64}

Analyze the distribution of normal directions in a mesh.

Returns statistics about how normals are distributed, useful for
debugging orientation issues.

# Returns
Dictionary with:
- "mean_x", "mean_y", "mean_z": Average normal direction
- "std_x", "std_y", "std_z": Standard deviation of each component
- "alignment_score": How aligned normals are (0=random, 1=perfectly aligned)

# Example
```julia
stats = analyze_normal_distribution(faces)
println("Alignment: \$(stats[\"alignment_score\"])")  # Should be close to 1.0
```
"""
function analyze_normal_distribution(faces::Vector{Triangle})
    isempty(faces) && return Dict{String,Float64}()
    
    normals = [compute_triangle_normal(tri) for tri in faces]
    
    # Filter out degenerate
    valid_normals = filter(n -> sum(abs.(n)) > 1e-6, normals)
    isempty(valid_normals) && return Dict{String,Float64}()
    
    # Compute statistics
    nx_vals = [n[1] for n in valid_normals]
    ny_vals = [n[2] for n in valid_normals]
    nz_vals = [n[3] for n in valid_normals]
    
    mean_n = (mean(nx_vals), mean(ny_vals), mean(nz_vals))
    std_n = (Statistics.std(nx_vals), Statistics.std(ny_vals), Statistics.std(nz_vals))
    
    # Alignment score: how well do normals align with mean direction?
    # Compute average dot product with mean normal
    mean_mag = sqrt(mean_n[1]^2 + mean_n[2]^2 + mean_n[3]^2)
    if mean_mag > 1e-6
        mean_unit = (mean_n[1]/mean_mag, mean_n[2]/mean_mag, mean_n[3]/mean_mag)
        dots = [n[1]*mean_unit[1] + n[2]*mean_unit[2] + n[3]*mean_unit[3] for n in valid_normals]
        alignment = mean(dots)  # 1.0 = perfect alignment, 0.0 = random, -1.0 = opposite
    else
        alignment = 0.0
    end
    
    return Dict(
        "mean_x" => mean_n[1],
        "mean_y" => mean_n[2],
        "mean_z" => mean_n[3],
        "std_x" => std_n[1],
        "std_y" => std_n[2],
        "std_z" => std_n[3],
        "alignment_score" => alignment,
        "valid_normals" => Float64(length(valid_normals)),
        "degenerate" => Float64(length(normals) - length(valid_normals))
    )
end

# ============================================================================
# Comprehensive Boundary Verification
# ============================================================================

"""
    BoundaryVerificationReport

Structured report from boundary consistency verification.
"""
struct BoundaryVerificationReport
    # Loop comparison
    loops_A::Int
    loops_B::Int
    loops_unified::Int
    matched_A::Int
    matched_B::Int
    unmatched_A::Int
    unmatched_B::Int
    
    # Normal consistency
    normals_consistent::Bool
    normals_flipped::Int
    normals_degenerate::Int
    normal_alignment_score::Float64
    
    # Overall assessment
    boundary_equivalent_A::Bool
    boundary_equivalent_B::Bool
    issues::Vector{String}
end

"""
    verify_boundary_consistency(
        unified_mesh::Vector{Triangle},
        mesh_A::Vector{Triangle},
        mesh_B::Vector{Triangle};
        verbose::Bool = false
    ) -> BoundaryVerificationReport

Comprehensive boundary consistency verification.

Checks:
1. Boundary loops match between unified and original meshes
2. Normal orientations are consistent
3. No topological defects

Returns detailed report suitable for compatibility assessment.

# Example
```julia
report = verify_boundary_consistency(unified_tris, faces_A, faces_B, verbose=true)

if !report.boundary_equivalent_A || !report.boundary_equivalent_B
    @error "Boundary mismatch detected!"
    for issue in report.issues
        println("  - \$issue")
    end
end
```
"""
function verify_boundary_consistency(
    unified_mesh::Vector{Triangle},
    mesh_A::Vector{Triangle},
    mesh_B::Vector{Triangle};
    verbose::Bool = false
)
    verbose && println("  → Performing boundary consistency verification...")
    
    issues = String[]
    
    # Extract boundary loops
    loops_A = extract_boundary_loops(mesh_A)
    loops_B = extract_boundary_loops(mesh_B)
    loops_unified = extract_boundary_loops(unified_mesh)
    
    verbose && println("    Boundary loops: A=$(length(loops_A)), B=$(length(loops_B)), unified=$(length(loops_unified))")
    
    # Compare loops
    (matched_A, unmatched_A_in_unified, unmatched_unified_vs_A) = compare_loop_sets(loops_unified, loops_A)
    (matched_B, unmatched_B_in_unified, unmatched_unified_vs_B) = compare_loop_sets(loops_unified, loops_B)
    
    # Boundary equivalence: all original loops must be in unified
    boundary_equivalent_A = (unmatched_A_in_unified == 0) && (matched_A == length(loops_A))
    boundary_equivalent_B = (unmatched_B_in_unified == 0) && (matched_B == length(loops_B))
    
    if !boundary_equivalent_A
        push!(issues, "Boundary mismatch vs A: $unmatched_A_in_unified of $(length(loops_A)) loops not in unified mesh")
    end
    
    if !boundary_equivalent_B
        push!(issues, "Boundary mismatch vs B: $unmatched_B_in_unified of $(length(loops_B)) loops not in unified mesh")
    end
    
    if unmatched_unified_vs_A > 0 || unmatched_unified_vs_B > 0
        push!(issues, "Unified mesh has extra boundary loops not in originals (A: $unmatched_unified_vs_A, B: $unmatched_unified_vs_B)")
    end
    
    verbose && println("    Loop matching: A=$(matched_A)/$(length(loops_A)), B=$(matched_B)/$(length(loops_B))")
    
    # Normal consistency check
    (normal_ok, flipped, degenerate) = verify_normal_consistency(unified_mesh, verbose=verbose)
    
    if !normal_ok
        push!(issues, "Normal inconsistency: $flipped flipped triangles, $degenerate degenerate")
    end
    
    if degenerate > 0
        push!(issues, "Warning: $degenerate degenerate triangles detected")
    end
    
    # Normal alignment analysis
    normal_stats = analyze_normal_distribution(unified_mesh)
    alignment_score = get(normal_stats, "alignment_score", 0.0)
    
    if alignment_score < 0.8
        push!(issues, "Poor normal alignment: $(round(alignment_score, digits=3)) < 0.8 (normals are scattered)")
    end
    
    verbose && println("    Normals: $(normal_ok ? "consistent" : "inconsistent ($flipped flipped)"), alignment=$(round(alignment_score, digits=3))")
    
    return BoundaryVerificationReport(
        length(loops_A),
        length(loops_B),
        length(loops_unified),
        matched_A,
        matched_B,
        unmatched_A_in_unified,
        unmatched_B_in_unified,
        normal_ok,
        flipped,
        degenerate,
        alignment_score,
        boundary_equivalent_A,
        boundary_equivalent_B,
        issues
    )
end
