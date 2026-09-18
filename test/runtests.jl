# test/runtests.jl – HydroPowerDynamics test suite
using Test

@testset "HydroPowerDynamics" begin
    include("test_connectors.jl")
    include("test_reservoir.jl")
    include("test_trollheim_reservoir_validation.jl")
    include("test_rigid_pipe.jl")
    include("test_standard_signals.jl")
    include("test_turbine_lookup.jl")
    include("test_utils.jl")
    include("test_francis_step.jl")
    include("test_openhpl_surge_tank.jl")
end
