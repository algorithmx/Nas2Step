module Nas2Step

using Gmsh
using Statistics

export SurfaceRegion, write_nas_surface
export VolumeRegion, write_nas_volume
export load_step_summary, verify_step, convert_nas_to_step, nas_to_step
export extract_boundary_surfaces, has_surface_elements
export check_nas_volumes, test_closed_shells
export verify_nas_mesh, check_element_volumes, check_element_quality, check_surface_closure, export_inverted_elements, test_node_swap_fix, check_vertex_coordination, export_mesh_quality_json, comprehensive_mesh_check
export check_region_overlap, export_region_overlap_json
export check_interface_conformity, export_interface_conformity_json, check_interface_conformity_json, export_interface_mismatch_surfaces

# Phase 1: Enhanced Interface Analysis (Surgical Repair)
export build_interface_topology, InterfaceTopology, EdgeKey, Triangle, BoundingBox
export classify_interface_mismatches, InterfaceClassification, EdgeMismatch, MismatchType
export classify_interface_mismatches_bidirectional  # Strategic improvement: classify from both perspectives
export build_boundary_constraints, BoundaryConstraints, check_constraint_violations
export export_interface_topology_json, export_classification_json, export_constraints_json

# Phase 1.0: Vertex Conformity Check (Fundamental Sanity Check)
export check_vertex_conformity, VertexConformityReport, ConformityLevel
export print_conformity_report, export_conformity_report_json

# Phase 2: Repair Strategy Generation
export generate_repair_plan, RepairPlan, EdgeInsertionPlan, QualityThresholds
export generate_repair_plan_bidirectional  # Strategic improvement: try both directions
export determine_dominant_side, export_repair_plan_json, default_thresholds
export export_interface_mismatches_json

# Phase 3: Surgical Mesh Repair Execution (disabled pending integration)
# export RepairWorkspace, create_checkpoint!, begin_transaction!, commit_transaction!, rollback_transaction!
# export delete_face!, add_face!, add_node!, get_face_by_nodes, get_node_id_by_coords
# export export_modified_mesh, print_workspace_stats
# export apply_quad_retriangulation!, apply_edge_insertion_plan!, apply_repair_plan!
# export execute_repairs_from_json, load_repair_plan_from_json
# export verify_interface_conformity, compare_conformity, print_conformity_report
# export print_improvement_report, verify_adjacent_interfaces, export_verification_report

include("write_nas.jl")

include("inspect.jl")

include("extract.jl")

include("nas_to_step.jl")

include("analyze.jl")

include("verify.jl")

include("repair/interface_topology.jl")
include("repair/geometric_utilities.jl")
include("repair/interface_conformity_check.jl")
include("repair/edge_classification.jl")
include("repair/boundary_constraints.jl")
include("repair/repair_planning.jl")
include("repair/repair_workspace.jl")
include("repair/repair_execution.jl")
include("repair/repair_verification.jl")
include("repair/interface_analysis_export.jl")

end # module Nas2Step
