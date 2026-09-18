using Test
using ModelingToolkit
using OrdinaryDiffEq
using HydroPowerDynamics

@testset "TurbineLookup hill-chart model" begin
    const rho = 1000.0
    const g = 9.81
    const H = 371.0
    const D = 2.5
    const gate0 = 1.0
    const n11_target = 50.0
    const Q_expected = 37.0
    const eta_expected = 0.97
    const n_rpm = n11_target * sqrt(H) / D
    const omega = 2pi * n_rpm / 60
    const P_expected = rho * g * Q_expected * H * eta_expected

    @mtkmodel HydraulicBoundaryZ begin
        @parameters begin
            p = 101_325.0
            z = 0.0
        end
        @components port = HydraulicPortZ()
        @equations begin
            port.p ~ p
            port.z ~ z
        end
    end

    @mtkmodel ShaftSpeedBoundary begin
        @parameters begin
            omega0 = omega
        end
        @components port = RotationalPort()
        @equations begin
            port.omega ~ omega0
            port.phi ~ omega0 * t
        end
    end

    @mtkmodel GateBoundary begin
        @parameters begin
            gate0 = 1.0
        end
        @components port = SignalOutPort()
        @equations port.u ~ gate0
    end

    @mtkcompile sys begin
        upper = HydraulicBoundaryZ(p = 101_325.0 + rho*g*H, z = 0.0)
        turbine = TurbineLookup(rho=rho, g=g, D_runner=D)
        lower = HydraulicBoundaryZ(p = 101_325.0, z = 0.0)
        shaft = ShaftSpeedBoundary(omega0 = omega)
        gate = GateBoundary(gate0 = gate0)

        @equations begin
            connect(upper.port, turbine.port_a)
            connect(turbine.port_b, lower.port)
            connect(turbine.shaft, shaft.port)
            connect(gate.port, turbine.gate)
        end
    end

    report = model_structure_report(sys)
    println("TurbineLookup structural report: ", report)
    @test report.equations == report.unknowns
    @test report.balanced
    @test report.parameters > 0

    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas5P(); reltol=1e-10, abstol=1e-12)

    @test isapprox(sol[sys.turbine.H][end], H; rtol=1e-9)
    @test isapprox(sol[sys.turbine.n11][end], n11_target; rtol=1e-8)
    @test isapprox(sol[sys.turbine.Q][end], Q_expected; rtol=1e-8)
    @test isapprox(sol[sys.turbine.eta][end], eta_expected; rtol=1e-10)
    @test isapprox(sol[sys.turbine.P_mech][end], P_expected; rtol=1e-8)

    # Mid-cell interpolation sanity check:
    # gate=0.8 lies halfway between gate_2=0.6 and gate_3=1.0 at n11=50.
    q11_mid_expected = 0.5 * (0.190 + 0.30735108593749)
    eta_mid_expected = 0.5 * (0.93 + 0.97)
    @test isapprox(
        bilinear3x3_clamped(
            0.8, 50.0,
            0.2,0.6,1.0, 40.0,50.0,60.0,
            0.060,0.062,0.060,
            0.185,0.190,0.184,
            0.300,0.30735108593749,0.298
        ),
        q11_mid_expected;
        atol=1e-12
    )
    @test isapprox(
        bilinear3x3_clamped(
            0.8, 50.0,
            0.2,0.6,1.0, 40.0,50.0,60.0,
            0.72,0.76,0.70,
            0.89,0.93,0.88,
            0.94,0.97,0.92
        ),
        eta_mid_expected;
        atol=1e-12
    )
end
