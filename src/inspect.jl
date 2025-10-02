module __Nas2StepInspect
using Gmsh

function _bbox_of_entities(ents::Vector{Tuple{Int32,Int32}})
    if isempty(ents)
        return nothing
    end
    # gmsh.model.getBoundingBox(dim, tag) returns (xmin,ymin,zmin,xmax,ymax,zmax)
    firstBox = gmsh.model.getBoundingBox(ents[1]...)
    xmin, ymin, zmin, xmax, ymax, zmax = firstBox
    for e in Iterators.drop(ents, 1)
        bx = gmsh.model.getBoundingBox(e...)
        xmin = min(xmin, bx[1])
        ymin = min(ymin, bx[2])
        zmin = min(zmin, bx[3])
        xmax = max(xmax, bx[4])
        ymax = max(ymax, bx[5])
        zmax = max(zmax, bx[6])
    end
    return (xmin, ymin, zmin, xmax, ymax, zmax)
end

function _summarize()
    vols = gmsh.model.getEntities(3)
    surfs = gmsh.model.getEntities(2)
    curves = gmsh.model.getEntities(1)
    pts = gmsh.model.getEntities(0)
    bbox3 = _bbox_of_entities(vols)
    bbox2 = _bbox_of_entities(surfs)
    return Dict(
        :counts => Dict(
            :volumes => length(vols),
            :surfaces => length(surfs),
            :curves => length(curves),
            :points => length(pts),
        ),
        :bbox => Dict(
            :volumes => bbox3,
            :surfaces => bbox2,
        ),
        :names => Dict(
            :volumes => [(t, gmsh.model.getEntityName(3, t)) for (_, t) in vols],
            :surfaces => [(t, gmsh.model.getEntityName(2, t)) for (_, t) in surfs],
        ),
    )
end
end # module __Nas2StepInspect


using .__Nas2StepInspect


"""
    load_step_summary(path::AbstractString) -> Dict

Open a STEP file in Gmsh and return a summary containing counts of entities,
bounding boxes for volumes/surfaces, and entity names.
"""
function load_step_summary(path::AbstractString)
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)  # Suppress terminal output
        gmsh.model.add("inspect_step")
        # Preferred: OCC import
        try
            gmsh.model.occ.importShapes(path)
            gmsh.model.occ.synchronize()
        catch
            # Fallback: generic open
            gmsh.clear()
            gmsh.model.add("inspect_step_open")
            gmsh.open(path)
            try
                gmsh.model.occ.synchronize()
            catch
            end
        end
        # If still empty, try alternate path
        if isempty(gmsh.model.getEntities())
            gmsh.clear()
            gmsh.model.add("inspect_step_retry")
            gmsh.open(path)
            try
                gmsh.model.occ.synchronize()
            catch
            end
        end
        return __Nas2StepInspect._summarize()
    finally
        gmsh.finalize()
    end
end


"""
    verify_step(path::AbstractString; verbose::Bool=false)

Print a concise summary of a STEP file and return the underlying Dict.
Set `verbose=true` for detailed entity names.
"""
function verify_step(path::AbstractString; verbose::Bool=false)
    info = load_step_summary(path)
    counts = info[:counts]
    bbox = info[:bbox]

    println("STEP file: $(basename(path))")
    println("  Entities: $(counts[:surfaces]) surfaces, $(counts[:curves]) curves, $(counts[:points]) points")

    if counts[:volumes] > 0
        println("  Volumes:  $(counts[:volumes])")
    end

    # Print bounding box in compact format
    if bbox[:surfaces] !== nothing
        b = bbox[:surfaces]
        println("  BBox: [$(round(b[1],digits=2)), $(round(b[2],digits=2)), $(round(b[3],digits=2))] → "
                *
                "[$(round(b[4],digits=2)), $(round(b[5],digits=2)), $(round(b[6],digits=2))]")
    elseif bbox[:volumes] !== nothing
        b = bbox[:volumes]
        println("  BBox: [$(round(b[1],digits=2)), $(round(b[2],digits=2)), $(round(b[3],digits=2))] → "
                *
                "[$(round(b[4],digits=2)), $(round(b[5],digits=2)), $(round(b[6],digits=2))]")
    end

    # Only show names in verbose mode
    if verbose
        names = info[:names]
        if !isempty(names[:volumes])
            println("  Volume names: ", [n for (t, n) in names[:volumes]])
        end
        if !isempty(names[:surfaces])
            println("  Surface names: ", [n for (t, n) in names[:surfaces]])
        end
    end

    return info
end
