using Test
using ModelingToolkit
using OrdinaryDiffEq
using HydroPowerDynamics

@testset "Trollheim reservoir analytical validation" begin
    # Trollheim benchmark values used in quick_examples/sim_test/
    const rho = 1000.0
    const g = 9.81
    const p_atm = 101_325.0
    const H_TROLLHEIM = 371.0

    @testset "Fixed reservoir hydrostatic pressure" begin
        @mtkmodel ReservoirProbe begin
            @variables begin
                x(t) = 0.0
                p_seen(t)
                z_seen(t)
            end
            @components begin
                port = HydraulicPortZ()
            end
            @equations begin
                port.dm ~ 0.0
                p_seen ~ port.p
                z_seen ~ port.z
                D(x) ~ 0.0
            end
        end

        @mtkcompile sys begin
            res = ReservoirBoundaryZ(
                h = H_TROLLHEIM,
                z_out = 0.0,
                rho = rho,
                g = g,
                p_atm = p_atm,
            )
            probe = ReservoirProbe()
            @equations begin
                connect(res.port, probe.port)
            end
        end

        prob = ODEProblem(sys, [sys.probe.x => 0.0], (0.0, 1.0))
        sol = solve(prob, Rodas5P(); reltol = 1e-10, abstol = 1e-12)

        p_expected = p_atm + rho * g * H_TROLLHEIM
        @test isapprox(sol[sys.probe.p_seen][end], p_expected; rtol = 1e-10, atol = 1e-6)
        @test isapprox(sol[sys.probe.z_seen][end], 0.0; atol = 1e-12)
        @test isapprox(sol[sys.res.h_abs][end], H_TROLLHEIM; atol = 1e-12)
    end

    @testset "Dynamic reservoir exact storage solution" begin
        # For alpha = 0:
        # A = h*W, V = L*W*h, m = rho*L*W*h
        # With constant outflow dm_out:
        # h(t) = h0 - dm_out/(rho*L*W)*t
        const L = 500.0
        const W = 100.0
        const h0 = H_TROLLHEIM
        const dm_out = 1000.0
        const tf = 10.0

        @mtkmodel ConstantMassFlowSinkZ begin
            @parameters begin
                dm = dm_out
            end
            @components begin
                port = HydraulicPortZ()
            end
            @equations begin
                port.dm ~ dm
            end
        end

        @mtkcompile sys begin
            res = DynamicReservoir(
                h_0 = h0,
                z_out = 0.0,
                L = L,
                W = W,
                alpha = 0.0,
                rho = rho,
                g = g,
                p_atm = p_atm,
            )
            sink = ConstantMassFlowSinkZ()
            @equations begin
                connect(res.port, sink.port)
            end
        end

        prob = ODEProblem(sys, [sys.res.h => h0], (0.0, tf))
        sol = solve(prob, Rodas5P(); reltol = 1e-10, abstol = 1e-12)

        h_expected = h0 - dm_out / (rho * L * W) * tf
        p_expected = p_atm + rho * g * h_expected
        q_expected = dm_out / rho

        @test isapprox(sol[sys.res.h][end], h_expected; rtol = 1e-9, atol = 1e-9)
        @test isapprox(sol[sys.res.h_abs][end], h_expected; rtol = 1e-9, atol = 1e-9)
        @test isapprox(sol[sys.res.port.p][end], p_expected; rtol = 1e-9, atol = 1e-5)
        @test isapprox(sol[sys.res.Q_out][end], q_expected; rtol = 1e-10, atol = 1e-10)
    end

    @testset "Trollheim hydraulic power reference" begin
        # Benchmark rated values from the existing Trollheim example.
        Q_rated = 37.0
        eta = 0.97
        P_hyd = rho * g * Q_rated * H_TROLLHEIM
        P_mech = P_hyd * eta

        @test isapprox(P_hyd / 1e6, 134.661; atol = 1e-3)
        @test isapprox(P_mech / 1e6, 130.6220139; atol = 1e-6)
    end
end
