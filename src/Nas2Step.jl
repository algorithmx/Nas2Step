module Nas2Step

using Gmsh
using Statistics

export SurfaceRegion, write_nas_surface
export VolumeRegion, write_nas_volume
export load_step_summary, verify_step, convert_nas_to_step, nas_to_step
export extract_boundary_surfaces, has_surface_elements
export check_nas_volumes, test_closed_shells
export verify_nas_mesh, check_element_volumes, check_element_quality, check_surface_closure, export_inverted_elements, test_node_swap_fix, check_vertex_coordination, export_mesh_quality_json, comprehensive_mesh_check

include("write_nas.jl")

include("inspect.jl")

include("extract.jl")

include("nas_to_step.jl")

include("analyze.jl")

include("verify.jl")


end # module Nas2Step
