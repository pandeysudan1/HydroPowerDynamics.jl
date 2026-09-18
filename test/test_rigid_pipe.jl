using Test
using ModelingToolkit
using OrdinaryDiffEq
using HydroPowerDynamics

@testset "Trollheim rigid-pipe analytical validation" begin
    const rho = 1000.0
    const g = 9.81
    const mu = 1.0e-3
    const roughness = 1.5e-5
    const L = 4496.5
    const Dpipe = 6.3
    const Qref = 37.0
    const Aref = pi * Dpipe^2 / 4
    const vref = Qref / Aref
    const dmref = rho * Qref
    const Reref = dmref * Dpipe / (mu * Aref)
    const fref = darcy_factor(Reref, Dpipe, roughness)
    const dpfref = fref * (L / Dpipe) * (rho / 2) * vref * abs(vref)
    const href = dpfref / (rho * g)

    @test isapprox(Aref, 31.1724531052; atol=1e-10)
    @test isapprox(vref, 1.18694540577; atol=1e-10)
    @test isapprox(Reref, 7.477756056e6; rtol=1e-10)
    @test isapprox(fref, 0.00869835224; rtol=1e-8)
    @test isapprox(dpfref, 4373.23213; atol=1e-3)
    @test isapprox(href, 0.445793286; atol=1e-6)

    @mtkmodel PressureElevationBoundary begin
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

    @mtkcompile steady_sys begin
        upstream = PressureElevationBoundary(p = 101_325.0 + dpfref, z = 0.0)
        pipe = RigidPipe(L=L, D_pipe=Dpipe, rho=rho, g=g, mu=mu, roughness=roughness)
        downstream = PressureElevationBoundary(p = 101_325.0, z = 0.0)
        @equations begin
            connect(upstream.port, pipe.port_a)
            connect(pipe.port_b, downstream.port)
        end
    end

    prob = ODEProblem(steady_sys, [steady_sys.pipe.dm => dmref], (0.0, 2.0))
    sol = solve(prob, Rodas5P(); reltol=1e-9, abstol=1e-11)

    @test isapprox(sol[steady_sys.pipe.Q][end], Qref; rtol=1e-6)
    @test isapprox(sol[steady_sys.pipe.Re][end], Reref; rtol=1e-6)
    @test isapprox(sol[steady_sys.pipe.f_D][end], fref; rtol=1e-6)
    @test isapprox(sol[steady_sys.pipe.dp_f][end], dpfref; rtol=1e-6)
    @test isapprox(sol[steady_sys.pipe.h_f][end], href; rtol=1e-6)

    # Zero-flow acceleration check: friction is zero initially, so the
    # momentum equation has a direct analytical derivative.
    const dp_step = 10_000.0
    const dmdt_expected = Aref / L * dp_step

    @mtkcompile accel_sys begin
        upstream = PressureElevationBoundary(p = 101_325.0 + dp_step, z = 0.0)
        pipe = RigidPipe(L=L, D_pipe=Dpipe, rho=rho, g=g, mu=mu, roughness=roughness)
        downstream = PressureElevationBoundary(p = 101_325.0, z = 0.0)
        @equations begin
            connect(upstream.port, pipe.port_a)
            connect(pipe.port_b, downstream.port)
        end
    end

    prob2 = ODEProblem(accel_sys, [accel_sys.pipe.dm => 0.0], (0.0, 1e-3))
    sol2 = solve(prob2, Rodas5P(); reltol=1e-10, abstol=1e-12, saveat=[0.0,1e-3])
    dmdt_numeric = (sol2[accel_sys.pipe.dm][2] - sol2[accel_sys.pipe.dm][1]) / 1e-3
    @test isapprox(dmdt_numeric, dmdt_expected; rtol=5e-3)
end
