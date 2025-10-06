using Test
using Nas2Step

const TEST_DIR = @__DIR__

# Try to import Logging stdlib for log suppression; fallback if unavailable
const _LOGGING_AVAILABLE = let ok = true
    try
        @eval import Logging
    catch
        ok = false
        # Avoid using @warn here to not require Logging
        println("[tests] Logging stdlib not available; test logs won't be suppressed.")
    end
    ok
end


# Bring Phase 3 workspace API into scope for tests (not exported by the package)
Base.include(@__MODULE__, joinpath(TEST_DIR, "..", "src", "repair", "repair_workspace.jl"))

# Minimal test helpers (kept tiny and local)
Base.include(@__MODULE__, joinpath(TEST_DIR, "test_utils.jl"))
using .Nas2StepTestUtils

@testset "Nas2Step unit tests" begin
    if _LOGGING_AVAILABLE
        Logging.with_logger(Logging.NullLogger()) do
            Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_interface_topology.jl"))
            Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_repair_workspace.jl"))
            Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_edge_classification.jl"))
            Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_boundary_constraints.jl"))
            Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_repair_planning.jl"))
            Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_repair_execution.jl"))
            Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_repair_verification.jl"))
        end
    else
        # Run without log suppression
        Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_interface_topology.jl"))
        Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_repair_workspace.jl"))
        Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_edge_classification.jl"))
        Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_boundary_constraints.jl"))
        Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_repair_planning.jl"))
        Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_repair_execution.jl"))
        Base.include(@__MODULE__, joinpath(TEST_DIR, "unit", "test_repair_verification.jl"))
    end
end

@testset "Nas2Step integration tests" begin
    if _LOGGING_AVAILABLE
        Logging.with_logger(Logging.NullLogger()) do
            # Always-on minimal pipeline
            Base.include(@__MODULE__, joinpath(TEST_DIR, "integration", "test_repair_pipeline_minimal.jl"))

            if get(ENV, "INTEGRATION", "0") == "1"
                Base.include(@__MODULE__, joinpath(TEST_DIR, "integration", "test_repair_pipeline_examples.jl"))
            end

            if get(ENV, "STRESS", "0") == "1"
                Base.include(@__MODULE__, joinpath(TEST_DIR, "integration", "test_repair_stress.jl"))
            end
        end
    else
        Base.include(@__MODULE__, joinpath(TEST_DIR, "integration", "test_repair_pipeline_minimal.jl"))
        if get(ENV, "INTEGRATION", "0") == "1"
            Base.include(@__MODULE__, joinpath(TEST_DIR, "integration", "test_repair_pipeline_examples.jl"))
        end
        if get(ENV, "STRESS", "0") == "1"
            Base.include(@__MODULE__, joinpath(TEST_DIR, "integration", "test_repair_stress.jl"))
        end
    end
end

@testset "Nas2Step property tests" begin
    if _LOGGING_AVAILABLE
        Logging.with_logger(Logging.NullLogger()) do
            Base.include(@__MODULE__, joinpath(TEST_DIR, "property", "test_invariants.jl"))
        end
    else
        Base.include(@__MODULE__, joinpath(TEST_DIR, "property", "test_invariants.jl"))
    end
end
