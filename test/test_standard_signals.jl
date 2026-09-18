using Test
using ModelingToolkit
using OrdinaryDiffEq
using HydroPowerDynamics

@testset "Standard input-signal validation" begin
    const rho = 1000.0
    const g = 9.81
    const p_atm = 101_325.0

    # ------------------------------------------------------------------
    # 1. Dynamic reservoir: exact step and ramp outflow response
    # For alpha = 0: m = rho*L*W*h, hence
    #   dh/dt = -Qout/(L*W).
    # Step Qout = Qs => h = h0 - Qs*t/(L*W).
    # Ramp Qout = kr*t => h = h0 - kr*t^2/(2*L*W).
    # ------------------------------------------------------------------
    @mtkmodel TimeVaryingMassFlowSinkZ begin
        @parameters begin
            mode = 0.0
            Q_step = 1.0
            k_ramp = 0.2
            t0 = 1.0
        end
        @components port = HydraulicPortZ()
        @equations begin
            port.dm ~ rho * ifelse(mode < 0.5,
                ifelse(t < t0, 0.0, Q_step),
                ifelse(t < t0, 0.0, k_ramp * (t - t0)))
        end
    end

    const Lr = 500.0
    const Wr = 100.0
    const h0 = 50.0

    @mtkcompile step_res_sys begin
        res = DynamicReservoir(h_0=h0, z_out=0.0, L=Lr, W=Wr, alpha=0.0, rho=rho, g=g, p_atm=p_atm)
        sink = TimeVaryingMassFlowSinkZ(mode=0.0, Q_step=1.0, t0=1.0)
        @equations connect(res.port, sink.port)
    end
    sol_step = solve(ODEProblem(step_res_sys, [step_res_sys.res.h => h0], (0.0, 6.0)),
                     Rodas5P(); tstops=[1.0], reltol=1e-10, abstol=1e-12, saveat=0.1)
    h_step_expected = h0 - 1.0 * (6.0 - 1.0) / (Lr * Wr)
    @test isapprox(sol_step[step_res_sys.res.h][end], h_step_expected; atol=1e-9, rtol=1e-9)

    @mtkcompile ramp_res_sys begin
        res = DynamicReservoir(h_0=h0, z_out=0.0, L=Lr, W=Wr, alpha=0.0, rho=rho, g=g, p_atm=p_atm)
        sink = TimeVaryingMassFlowSinkZ(mode=1.0, k_ramp=0.2, t0=1.0)
        @equations connect(res.port, sink.port)
    end
    sol_ramp = solve(ODEProblem(ramp_res_sys, [ramp_res_sys.res.h => h0], (0.0, 6.0)),
                     Rodas5P(); tstops=[1.0], reltol=1e-10, abstol=1e-12, saveat=0.1)
    h_ramp_expected = h0 - 0.2 * (6.0 - 1.0)^2 / (2 * Lr * Wr)
    @test isapprox(sol_ramp[ramp_res_sys.res.h][end], h_ramp_expected; atol=1e-9, rtol=1e-9)

    # ------------------------------------------------------------------
    # 2. RigidPipe: initial response to pressure step and pressure ramp
    # At dm = 0, Darcy loss is zero, so
    #   d(dm)/dt = A/L * Delta p.
    # ------------------------------------------------------------------
    @mtkmodel PressureSignalBoundaryZ begin
        @parameters begin
            p0 = p_atm
            dp = 10_000.0
            ramp_rate = 20_000.0
            mode = 0.0
        end
        @components port = HydraulicPortZ()
        @equations begin
            port.p ~ p0 + ifelse(mode < 0.5, dp, ramp_rate * t)
            port.z ~ 0.0
        end
    end

    const Lp = 4496.5
    const Dp = 6.3
    const Ap = pi * Dp^2 / 4

    @mtkcompile pipe_step_sys begin
        up = PressureSignalBoundaryZ(mode=0.0, dp=10_000.0)
        pipe = RigidPipe(L=Lp, D_pipe=Dp, rho=rho, g=g, mu=1e-3, roughness=1.5e-5)
        dn = PressureSignalBoundaryZ(mode=0.0, dp=0.0)
        @equations begin
            connect(up.port, pipe.port_a)
            connect(pipe.port_b, dn.port)
        end
    end
    dt = 1e-4
    sol_pipe_step = solve(ODEProblem(pipe_step_sys, [pipe_step_sys.pipe.dm => 0.0], (0.0, dt)),
                          Rodas5P(); reltol=1e-11, abstol=1e-13, saveat=[0.0, dt])
    dmdt_num = (sol_pipe_step[pipe_step_sys.pipe.dm][2] - sol_pipe_step[pipe_step_sys.pipe.dm][1]) / dt
    dmdt_exact = Ap / Lp * 10_000.0
    @test isapprox(dmdt_num, dmdt_exact; rtol=2e-3)

    # For a ramp Delta p = r*t, initially d(dm)/dt = 0 and
    # d^2(dm)/dt^2 = A/L*r. Thus dm ≈ 0.5*A/L*r*t^2.
    @mtkcompile pipe_ramp_sys begin
        up = PressureSignalBoundaryZ(mode=1.0, ramp_rate=20_000.0)
        pipe = RigidPipe(L=Lp, D_pipe=Dp, rho=rho, g=g, mu=1e-3, roughness=1.5e-5)
        dn = PressureSignalBoundaryZ(mode=0.0, dp=0.0)
        @equations begin
            connect(up.port, pipe.port_a)
            connect(pipe.port_b, dn.port)
        end
    end
    dt2 = 1e-3
    sol_pipe_ramp = solve(ODEProblem(pipe_ramp_sys, [pipe_ramp_sys.pipe.dm => 0.0], (0.0, dt2)),
                          Rodas5P(); reltol=1e-11, abstol=1e-13, saveat=[0.0, dt2])
    dm_ramp_exact = 0.5 * (Ap / Lp) * 20_000.0 * dt2^2
    @test isapprox(sol_pipe_ramp[pipe_ramp_sys.pipe.dm][end], dm_ramp_exact; rtol=5e-3, atol=1e-12)

    # ------------------------------------------------------------------
    # 3. OpenHPLSurgeTank: hydrostatic equilibrium + pressure-step onset
    # Vertical shaft H=L. At equilibrium p = p_atm + rho*g*h0.
    # A pressure step Delta p gives initial dVdot/dt = Delta p/(rho*l).
    # ------------------------------------------------------------------
    @mtkmodel SurgePressureBoundary begin
        @parameters begin
            h0 = 50.0
            dp = 0.0
        end
        @components port = HydraulicPort()
        @equations port.p ~ p_atm + rho*g*h0 + dp
    end

    const Hs = 87.0
    const Ls = 87.0
    const Ds = 3.4
    const hs0 = 50.0

    @mtkcompile surge_step_sys begin
        b = SurgePressureBoundary(h0=hs0, dp=1000.0)
        tank = OpenHPLSurgeTank(H=Hs, L=Ls, D_tank=Ds, h_0=hs0, Vdot_0=0.0, rho=rho, g=g, p_atm=p_atm, f=0.02)
        @equations connect(b.port, tank.port)
    end
    dts = 1e-4
    sol_surge = solve(ODEProblem(surge_step_sys, [surge_step_sys.tank.h => hs0, surge_step_sys.tank.Vdot => 0.0], (0.0, dts)),
                      Rodas5P(); reltol=1e-11, abstol=1e-13, saveat=[0.0, dts])
    l0 = hs0 / (Hs/Ls)
    dVdt_exact = 1000.0 / (rho * l0)
    dVdt_num = (sol_surge[surge_step_sys.tank.Vdot][2] - sol_surge[surge_step_sys.tank.Vdot][1]) / dts
    @test isapprox(dVdt_num, dVdt_exact; rtol=5e-3)
end
