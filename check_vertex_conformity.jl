#!/usr/bin/env julia
# Check vertex conformity for mesh interfaces
# This examines whether interfaces share the same vertices (fundamental requirement)

using Pkg
Pkg.activate(".")

using Nas2Step
using Gmsh: gmsh

# Load the conformity check module
include("src/repair/interface_conformity_check.jl")

"""
Check vertex conformity for a specific interface or all interfaces.
"""
function check_conformity(nas_file::String; 
                         pidA::Union{Int,Nothing}=nothing,
                         pidB::Union{Int,Nothing}=nothing,
                         tol::Real=1e-4,
                         verbose::Bool=false,
                         export_json::Bool=false)
    
    println("="^70)
    println("Vertex Conformity Check")
    println("="^70)
    println("File: $nas_file")
    println("Tolerance: $tol")
    println()
    
    # Check if we're examining a specific interface or all
    if pidA !== nothing && pidB !== nothing
        # Single interface
        println("Building topology for PID $pidA ↔ PID $pidB...")
        topology = build_interface_topology(nas_file, pidA, pidB, tol=tol)
        
        report = check_vertex_conformity(topology, tol=tol)
        print_conformity_report(report, verbose=verbose)
        
        if export_json
            output_file = "conformity_$(pidA)_$(pidB).json"
            export_conformity_report_json(report, output_file)
        end
        
        return [report]
    else
        # All interfaces
        println("Finding all interfaces...")
        
        # Find interface pairs (simplified version)
        interface_pairs = Tuple{Int,Int}[]
        
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        
        try
            gmsh.open(nas_file)
            vols = [tag for (dim, tag) in gmsh.model.getEntities(3)]
            
            if length(vols) < 2
                println("Found only $(length(vols)) volume(s). Need at least 2 for interfaces.")
                return VertexConformityReport[]
            end
            
            println("Found $(length(vols)) volumes")
            
            # Build node lookups
            node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
            tag_to_coord = Dict{Int,NTuple{3,Float64}}()
            for (i, t) in enumerate(node_tags)
                base = (i-1)*3
                if base + 3 <= length(node_coords)
                    tag_to_coord[Int(t)] = (node_coords[base+1], node_coords[base+2], node_coords[base+3])
                end
            end
            
            # Get nodes per volume
            ckey(p) = (round(p[1], digits=4), round(p[2], digits=4), round(p[3], digits=4))
            
            vol_nodes = Dict{Int,Set{NTuple{3,Float64}}}()
            for v in vols
                etypes, _, enodes = gmsh.model.mesh.getElements(3, v)
                nodes_here = Set{NTuple{3,Float64}}()
                for (t, ev) in zip(etypes, enodes)
                    if t != 4; continue; end  # 4 = tetra
                    for i in 1:4:length(ev)
                        for k in 0:3
                            nid = Int(ev[i+k])
                            coord = get(tag_to_coord, nid, nothing)
                            if coord !== nothing
                                push!(nodes_here, ckey(coord))
                            end
                        end
                    end
                end
                vol_nodes[v] = nodes_here
            end
            
            # Find pairs with shared nodes
            for i in 1:length(vols)
                for j in (i+1):length(vols)
                    a, b = vols[i], vols[j]
                    if length(intersect(vol_nodes[a], vol_nodes[b])) >= 10
                        push!(interface_pairs, (a, b))
                    end
                end
            end
            
            println("Found $(length(interface_pairs)) interface(s) with shared nodes")
            println()
            
        finally
            gmsh.finalize()
        end
        
        if isempty(interface_pairs)
            println("No interfaces found.")
            return VertexConformityReport[]
        end
        
        # Check each interface
        reports = VertexConformityReport[]
        for (i, (a, b)) in enumerate(interface_pairs)
            println("\n[$i/$(length(interface_pairs))] Checking PID $a ↔ PID $b...")
            
            topology = build_interface_topology(nas_file, a, b, tol=tol)
            report = check_vertex_conformity(topology, tol=tol)
            
            # Print summary
            status = report.is_acceptable ? "✓" : "✗"
            println("  $status $(report.conformity_level)")
            println("  Vertex match: $(round(report.vertex_match_ratio * 100, digits=1))%")
            println("  Vertices: $(length(report.shared_vertices)) shared, " *
                   "$(length(report.vertices_only_A)) unique to A, " *
                   "$(length(report.vertices_only_B)) unique to B")
            
            if !report.is_acceptable
                println("  Issues: $(length(report.issues))")
            end
            
            push!(reports, report)
            
            if export_json
                output_file = "conformity_$(a)_$(b).json"
                export_conformity_report_json(report, output_file)
            end
        end
        
        # Summary across all interfaces
        println("\n" * "="^70)
        println("SUMMARY ACROSS ALL INTERFACES")
        println("="^70)
        
        level_counts = Dict{ConformityLevel,Int}()
        for level in instances(ConformityLevel)
            level_counts[level] = count(r -> r.conformity_level == level, reports)
        end
        
        println("\nConformity Level Distribution:")
        for level in instances(ConformityLevel)
            count = level_counts[level]
            if count > 0
                pct = round(100 * count / length(reports), digits=1)
                println("  $(level): $count ($pct%)")
            end
        end
        
        acceptable_count = count(r -> r.is_acceptable, reports)
        println("\nRepairability:")
        println("  Acceptable: $acceptable_count / $(length(reports)) " *
               "($(round(100 * acceptable_count / length(reports), digits=1))%)")
        
        return reports
    end
end

# Main execution
if length(ARGS) < 1
    println("""
Usage: julia check_vertex_conformity.jl <nas_file> [pidA] [pidB] [options]

Arguments:
  nas_file    Path to NAS file
  pidA        (Optional) First PID - if specified, must also specify pidB
  pidB        (Optional) Second PID - analyzes single interface
  
Options (must come after PIDs if specified):
  --verbose   Show detailed cluster information
  --json      Export reports to JSON files
  --tol=X     Set tolerance (default: 1e-4)

Examples:
  # Check all interfaces
  julia check_vertex_conformity.jl mesh.nas
  
  # Check specific interface
  julia check_vertex_conformity.jl mesh.nas 1 3
  
  # Check with options
  julia check_vertex_conformity.jl mesh.nas 1 3 --verbose --json
  julia check_vertex_conformity.jl mesh.nas --tol=1e-3 --json
""")
    exit(1)
end

nas_file = ARGS[1]

if !isfile(nas_file)
    println("ERROR: File not found: $nas_file")
    exit(1)
end

# Parse arguments
pidA = nothing
pidB = nothing
verbose = false
export_json = false
tol = 1e-4

arg_idx = 2
# Check for PIDs (must be first two non-option args)
if length(ARGS) >= 2 && !startswith(ARGS[2], "--")
    pidA = parse(Int, ARGS[2])
    if length(ARGS) >= 3 && !startswith(ARGS[3], "--")
        pidB = parse(Int, ARGS[3])
        arg_idx = 4
    else
        println("ERROR: If pidA is specified, pidB must also be specified")
        exit(1)
    end
end

# Parse options
for i in arg_idx:length(ARGS)
    arg = ARGS[i]
    if arg == "--verbose"
        global verbose = true
    elseif arg == "--json"
        global export_json = true
    elseif startswith(arg, "--tol=")
        global tol = parse(Float64, arg[7:end])
    else
        println("WARNING: Unknown option: $arg")
    end
end

try
    reports = check_conformity(nas_file, pidA=pidA, pidB=pidB, tol=tol, 
                              verbose=verbose, export_json=export_json)
    
    println("\n✓ Conformity check complete!")
    
    if !isempty(reports)
        problem_count = count(r -> !r.is_acceptable, reports)
        if problem_count > 0
            println("\n⚠ Warning: $problem_count interface(s) have vertex conformity issues")
        end
    end
catch e
    println("\n✗ Error during conformity check:")
    println(e)
    rethrow()
end
