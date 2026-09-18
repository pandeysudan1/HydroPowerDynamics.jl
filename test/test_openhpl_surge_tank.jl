using Test
using ModelingToolkit
using OrdinaryDiffEq
using HydroPowerDynamics

@testset "OpenHPLSurgeTank" begin
    # Hydrostatic equilibrium: bottom pressure exactly balances gravity.
    @mtkmodel HydrostaticBoundary begin
        @parameters begin
            rho = 1000.0
            g = 9.81
            p_atm = 101_325.0
            h0 = 50.0
        end
        @components begin
            port = HydraulicPort()
        end
        @equations begin
            port.p ~ p_atm + rho * g * h0
        end
    end

    @mtkbuild sys begin
        boundary = HydrostaticBoundary()
        tank = OpenHPLSurgeTank(
            H = 87.0,
            L = 87.0,
            D_tank = 3.4,
            h_0 = 50.0,
            Vdot_0 = 0.0,
            f = 0.02,
        )
        @equations begin
            connect(boundary.port, tank.port)
        end
    end

    prob = ODEProblem(sys, [], (0.0, 2.0))
    sol = solve(prob, Rodas5P(); reltol = 1e-8, abstol = 1e-9)

    @test SciMLBase.successful_retcode(sol)
    @test abs(sol[sys.tank.h][end] - 50.0) < 1e-5
    @test abs(sol[sys.tank.Vdot][end]) < 1e-5
end
