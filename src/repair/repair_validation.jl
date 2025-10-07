"""
    repair_validation.jl

Phase 5: Post-Repair Validation

Validates mesh quality, manifoldness, and interface conformity after repairs are applied.
"""

using LinearAlgebra
using Printf

"""
    ValidationResult

Result of mesh validation check.
"""
struct ValidationResult
    passed::Bool
    errors::Vector{String}
    warnings::Vector{String}
    stats::Dict{String,Any}
end

"""
    validate_mesh_quality(ws::RepairWorkspace; thresholds::QualityThresholds) -> ValidationResult

Check mesh quality after repairs.
Validates:
- Minimum angle constraints
- Maximum angle constraints
- Aspect ratio constraints
- Degenerate triangles
"""
function validate_mesh_quality(ws::RepairWorkspace; thresholds::QualityThresholds = default_thresholds())
    errors = String[]
    warnings = String[]
    
    # Statistics
    min_angle_global = Inf
    max_angle_global = 0.0
    max_aspect_ratio = 0.0
    degenerate_count = 0
    poor_quality_count = 0
    
    total_faces = 0
    
    for (pid, faces) in ws.working_faces
        for face_nodes in faces
            total_faces += 1
            
            if length(face_nodes) != 3
                push!(errors, "Non-triangular face in PID $pid with $(length(face_nodes)) nodes")
                continue
            end
            
            # Get node coordinates
            coords = [ws.working_nodes[nid] for nid in face_nodes]
            
            # Check for degenerate triangle (area ≈ 0)
            v1 = [coords[2][i] - coords[1][i] for i in 1:3]
            v2 = [coords[3][i] - coords[1][i] for i in 1:3]
            cross_prod = [
                v1[2]*v2[3] - v1[3]*v2[2],
                v1[3]*v2[1] - v1[1]*v2[3],
                v1[1]*v2[2] - v1[2]*v2[1]
            ]
            area = 0.5 * sqrt(sum(cross_prod.^2))
            
            if area < 1e-10
                push!(errors, "Degenerate triangle in PID $pid (area ≈ 0)")
                degenerate_count += 1
                continue
            end
            
            # Compute angles
            edges = [
                sqrt(sum((coords[2][i] - coords[1][i])^2 for i in 1:3)),
                sqrt(sum((coords[3][i] - coords[2][i])^2 for i in 1:3)),
                sqrt(sum((coords[1][i] - coords[3][i])^2 for i in 1:3))
            ]
            
            # Use law of cosines for angles
            angles_deg = Float64[]
            for i in 1:3
                j = mod1(i+1, 3)
                k = mod1(i+2, 3)
                cos_angle = (edges[i]^2 + edges[k]^2 - edges[j]^2) / (2 * edges[i] * edges[k])
                cos_angle = clamp(cos_angle, -1.0, 1.0)
                angle_rad = acos(cos_angle)
                angle_deg = rad2deg(angle_rad)
                push!(angles_deg, angle_deg)
            end
            
            min_angle = minimum(angles_deg)
            max_angle = maximum(angles_deg)
            
            min_angle_global = min(min_angle_global, min_angle)
            max_angle_global = max(max_angle_global, max_angle)
            
            # Check against thresholds
            if min_angle < thresholds.min_angle
                poor_quality_count += 1
                if min_angle < thresholds.min_angle * 0.5
                    push!(errors, "Critical min angle violation in PID $pid: $(round(min_angle, digits=2))° < $(thresholds.min_angle)°")
                else
                    push!(warnings, "Min angle below threshold in PID $pid: $(round(min_angle, digits=2))°")
                end
            end
            
            if max_angle > thresholds.max_angle
                if max_angle > thresholds.max_angle * 1.1
                    push!(errors, "Critical max angle violation in PID $pid: $(round(max_angle, digits=2))° > $(thresholds.max_angle)°")
                else
                    push!(warnings, "Max angle above threshold in PID $pid: $(round(max_angle, digits=2))°")
                end
            end
            
            # Aspect ratio (longest edge / shortest edge)
            aspect_ratio = maximum(edges) / minimum(edges)
            max_aspect_ratio = max(max_aspect_ratio, aspect_ratio)
            
            if aspect_ratio > thresholds.max_aspect_ratio
                push!(warnings, "High aspect ratio in PID $pid: $(round(aspect_ratio, digits=2))")
            end
        end
    end
    
    stats = Dict{String,Any}(
        "total_faces" => total_faces,
        "min_angle_global" => round(min_angle_global, digits=2),
        "max_angle_global" => round(max_angle_global, digits=2),
        "max_aspect_ratio" => round(max_aspect_ratio, digits=2),
        "degenerate_count" => degenerate_count,
        "poor_quality_count" => poor_quality_count
    )
    
    passed = isempty(errors)
    
    return ValidationResult(passed, errors, warnings, stats)
end

"""
    validate_manifoldness(ws::RepairWorkspace) -> ValidationResult

Check mesh manifoldness after repairs.
Validates:
- Each edge is shared by exactly 1 or 2 faces
- No duplicate faces
- No isolated nodes
"""
function validate_manifoldness(ws::RepairWorkspace)
    errors = String[]
    warnings = String[]
    
    # Build edge incidence map
    edge_incidence = Dict{NTuple{2,Int}, Vector{Tuple{Int,Int}}}()  # edge => [(pid, face_idx), ...]
    
    for (pid, faces) in ws.working_faces
        for (face_idx, face_nodes) in enumerate(faces)
            if length(face_nodes) != 3
                continue
            end
            
            # Add each edge
            for i in 1:3
                j = mod1(i+1, 3)
                n1, n2 = face_nodes[i], face_nodes[j]
                edge = n1 < n2 ? (n1, n2) : (n2, n1)  # Canonical form
                
                if !haskey(edge_incidence, edge)
                    edge_incidence[edge] = []
                end
                push!(edge_incidence[edge], (pid, face_idx))
            end
        end
    end
    
    # Check edge incidence
    boundary_edge_count = 0
    non_manifold_edge_count = 0
    
    for (edge, incidents) in edge_incidence
        n_inc = length(incidents)
        
        if n_inc == 1
            boundary_edge_count += 1
        elseif n_inc > 2
            non_manifold_edge_count += 1
            pids = unique([pid for (pid, _) in incidents])
            push!(errors, "Non-manifold edge $edge shared by $(n_inc) faces across PIDs: $pids")
        end
    end
    
    # Check for duplicate faces
    face_set = Set{NTuple{3,Int}}()
    duplicate_count = 0
    
    for (pid, faces) in ws.working_faces
        for face_nodes in faces
            if length(face_nodes) != 3
                continue
            end
            
            # Canonical form (sorted)
            face_canonical = tuple(sort(face_nodes)...)
            
            if face_canonical ∈ face_set
                duplicate_count += 1
                push!(errors, "Duplicate face $face_canonical in PID $pid")
            else
                push!(face_set, face_canonical)
            end
        end
    end
    
    # Check for isolated nodes
    used_nodes = Set{Int}()
    for faces in values(ws.working_faces)
        for face_nodes in faces
            for nid in face_nodes
                push!(used_nodes, nid)
            end
        end
    end
    
    isolated_count = 0
    for node_id in keys(ws.working_nodes)
        if node_id ∉ used_nodes
            isolated_count += 1
            push!(warnings, "Isolated node $node_id (not used by any face)")
        end
    end
    
    stats = Dict{String,Any}(
        "total_edges" => length(edge_incidence),
        "boundary_edges" => boundary_edge_count,
        "non_manifold_edges" => non_manifold_edge_count,
        "duplicate_faces" => duplicate_count,
        "isolated_nodes" => isolated_count
    )
    
    passed = isempty(errors)
    
    return ValidationResult(passed, errors, warnings, stats)
end

"""
    validate_interface_conformity(ws::RepairWorkspace, interface_pair::Tuple{Int,Int}) -> ValidationResult

Check interface conformity after repairs.
Validates that the interface between two PIDs is now conforming.
"""
function validate_interface_conformity(ws::RepairWorkspace, interface_pair::Tuple{Int,Int})
    errors = String[]
    warnings = String[]
    
    pid_a, pid_b = interface_pair
    
    # Get faces from both PIDs
    if !haskey(ws.working_faces, pid_a) || !haskey(ws.working_faces, pid_b)
        push!(errors, "One or both PIDs not found in workspace")
        return ValidationResult(false, errors, warnings, Dict{String,Any}())
    end
    
    faces_a = ws.working_faces[pid_a]
    faces_b = ws.working_faces[pid_b]
    
    # Build edge sets
    edges_a = Set{NTuple{2,NTuple{3,Float64}}}()
    edges_b = Set{NTuple{2,NTuple{3,Float64}}}()
    
    for face_nodes in faces_a
        if length(face_nodes) == 3
            coords = [ws.working_nodes[nid] for nid in face_nodes]
            for i in 1:3
                j = mod1(i+1, 3)
                c1, c2 = coords[i], coords[j]
                edge = c1 < c2 ? (c1, c2) : (c2, c1)
                push!(edges_a, edge)
            end
        end
    end
    
    for face_nodes in faces_b
        if length(face_nodes) == 3
            coords = [ws.working_nodes[nid] for nid in face_nodes]
            for i in 1:3
                j = mod1(i+1, 3)
                c1, c2 = coords[i], coords[j]
                edge = c1 < c2 ? (c1, c2) : (c2, c1)
                push!(edges_b, edge)
            end
        end
    end
    
    # Find interface edges (shared between both PIDs)
    interface_edges_a = edges_a ∩ edges_b
    interface_edges_b = edges_b ∩ edges_a
    
    # Check for edge mismatches
    edges_only_a = setdiff(edges_a, edges_b)
    edges_only_b = setdiff(edges_b, edges_a)
    
    # Count mismatches on interface (edges that should be shared but aren't)
    # This is a simplified check - a proper check would need the full topology
    mismatch_count = 0
    
    if !isempty(edges_only_a) || !isempty(edges_only_b)
        push!(warnings, "Interface may still have edge mismatches: $(length(edges_only_a)) unique to A, $(length(edges_only_b)) unique to B")
        mismatch_count = length(edges_only_a) + length(edges_only_b)
    end
    
    stats = Dict{String,Any}(
        "edges_a" => length(edges_a),
        "edges_b" => length(edges_b),
        "shared_edges" => length(interface_edges_a),
        "unique_to_a" => length(edges_only_a),
        "unique_to_b" => length(edges_only_b),
        "potential_mismatches" => mismatch_count
    )
    
    passed = mismatch_count == 0
    
    return ValidationResult(passed, errors, warnings, stats)
end

"""
    validate_repairs(ws::RepairWorkspace, interface_pair::Tuple{Int,Int}; thresholds::QualityThresholds) -> ValidationResult

Comprehensive validation of repairs.
Combines all validation checks.
"""
function validate_repairs(
    ws::RepairWorkspace,
    interface_pair::Tuple{Int,Int};
    thresholds::QualityThresholds = default_thresholds()
)
    println("\n" * "="^70)
    println("Repair Validation")
    println("="^70)
    
    all_errors = String[]
    all_warnings = String[]
    all_stats = Dict{String,Any}()
    
    # 1. Quality validation
    println("\n1. Mesh Quality Check...")
    quality_result = validate_mesh_quality(ws, thresholds=thresholds)
    append!(all_errors, quality_result.errors)
    append!(all_warnings, quality_result.warnings)
    all_stats["quality"] = quality_result.stats
    
    if quality_result.passed
        println("   ✓ PASSED")
    else
        println("   ✗ FAILED ($(length(quality_result.errors)) errors)")
    end
    println("     Min angle: $(quality_result.stats["min_angle_global"])°")
    println("     Max angle: $(quality_result.stats["max_angle_global"])°")
    println("     Max aspect ratio: $(quality_result.stats["max_aspect_ratio"])")
    
    # 2. Manifoldness validation
    println("\n2. Manifoldness Check...")
    manifold_result = validate_manifoldness(ws)
    append!(all_errors, manifold_result.errors)
    append!(all_warnings, manifold_result.warnings)
    all_stats["manifoldness"] = manifold_result.stats
    
    if manifold_result.passed
        println("   ✓ PASSED")
    else
        println("   ✗ FAILED ($(length(manifold_result.errors)) errors)")
    end
    println("     Total edges: $(manifold_result.stats["total_edges"])")
    println("     Boundary edges: $(manifold_result.stats["boundary_edges"])")
    println("     Non-manifold edges: $(manifold_result.stats["non_manifold_edges"])")
    
    # 3. Interface conformity validation
    println("\n3. Interface Conformity Check...")
    conformity_result = validate_interface_conformity(ws, interface_pair)
    append!(all_errors, conformity_result.errors)
    append!(all_warnings, conformity_result.warnings)
    all_stats["conformity"] = conformity_result.stats
    
    if conformity_result.passed
        println("   ✓ PASSED")
    else
        println("   ✗ FAILED ($(length(conformity_result.errors)) errors)")
    end
    println("     Shared edges: $(conformity_result.stats["shared_edges"])")
    println("     Potential mismatches: $(conformity_result.stats["potential_mismatches"])")
    
    # Overall result
    println("\n" * "="^70)
    overall_passed = quality_result.passed && manifold_result.passed && conformity_result.passed
    
    if overall_passed
        println("✓ ALL VALIDATIONS PASSED")
    else
        println("✗ VALIDATION FAILED")
        println("\nErrors: $(length(all_errors))")
        for err in all_errors[1:min(10, length(all_errors))]
            println("  • $err")
        end
        if length(all_errors) > 10
            println("  ... and $(length(all_errors)-10) more")
        end
    end
    
    if !isempty(all_warnings)
        println("\nWarnings: $(length(all_warnings))")
        for warn in all_warnings[1:min(5, length(all_warnings))]
            println("  • $warn")
        end
        if length(all_warnings) > 5
            println("  ... and $(length(all_warnings)-5) more")
        end
    end
    
    println("="^70)
    
    return ValidationResult(overall_passed, all_errors, all_warnings, all_stats)
end
