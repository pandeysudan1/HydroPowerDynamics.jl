# Trollheim AGC comparison using HydroPowerDynamics.jl
# Same 75 -> 90 MW electrical disturbance; only hydraulic topology differs.

using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq
using Printf

const PBASE = 150e6
const P0 = 75e6
const P1 = 90e6
const TSTEP = 5.0
const F0 = 50.0
const POLES = 16
const RPM0 = 120 * F0 / POLES
const W0 = 2π * RPM0 / 60
const H_GROSS = 371.0
const RHO = 1000.0
const G = 9.81
const P_ATM = 101325.0
const P_UP = P_ATM + RHO * G * H_GROSS
const Q_RATED = 37.0
const D_RUNNER = 2.5
const ETA_MAX = 0.97
const C_ETA = 0.25
const TAU_RATED = 0.90
const KQ = Q_RATED / (TAU_RATED * D_RUNNER^2 * sqrt(H_GROSS))
const TAU0 = 0.5384797177203611
const Q0 = 22.13749950628151
const DM0 = RHO * Q0
const H_INERTIA = 1.83
const JTOTAL = 2 * H_INERTIA * PBASE / W0^2
const R_PU = 0.50
const TGOV = 0.30
const KI = 0.10
const L_HEADRACE = 500.0
const D_HEADRACE = 6.0
const L_PENSTOCK = 500.0
const D_PENSTOCK = 4.0
const ROUGHNESS = 1.5e-5
const ATANK = π / 4 * 3.4^2
const ZSURGE0 = 69.97
const D_RISER = 3.4
const L_RISER = 87.0
const E_RISER = 1e-2

function pipe_operating_guesses(L, Dpipe)
    A = π * Dpipe^2 / 4
    mu = 1e-3
    Re = DM0 * Dpipe / (mu * A)
    fD = 1.0 / (2.0 * log10(ROUGHNESS / (3.7 * Dpipe) + 5.7 / Re^0.9))^2
    dpf = fD * (L / Dpipe) * (DM0 / (RHO * A))^2 * RHO / 2
    (; Re, fD, dpf)
end

const HRG = pipe_operating_guesses(L_HEADRACE, D_HEADRACE)
const PNG = pipe_operating_guesses(L_PENSTOCK, D_PENSTOCK)
const ETA0 = ETA_MAX * (1 - C_ETA * (Q0 / Q_RATED - 1)^2)
const TAU_SHAFT0 = P0 / W0

function build_trollheim(; with_surge::Bool)
    @named upper = Reservoir(H=H_GROSS, rho=RHO, g=G, p_atm=P_ATM)
    @named tail = Reservoir(H=0.0, rho=RHO, g=G, p_atm=P_ATM)
    @named headrace = Penstock(L=L_HEADRACE, D_pipe=D_HEADRACE, rho=RHO, B=2e9, e=ROUGHNESS)
    @named penstock = Penstock(L=L_PENSTOCK, D_pipe=D_PENSTOCK, rho=RHO, B=2e9, e=ROUGHNESS)
    @named turbine = FrancisTurbineAffinity(D=D_RUNNER, K_q=KQ, eta_max=ETA_MAX, c_eta=C_ETA,
                                            Q_rated=Q_RATED, rho=RHO, g=G)
    @named rotor = RotorInertia(J=JTOTAL, b_v=0.0, omega_0=W0)
    @named generator = VariableLoadGenerator(omega_rated=W0, omega_s=W0, D_d=0.0, eta_gen=1.0)
    @named speed = RotationalSpeedSensor()
    @named governor = SimpleGovernorAGC(R=R_PU, T_g=TGOV, K_i=KI, omega_ref=W0,
                                        gate_bias=TAU0, gate_min=0.05, gate_max=1.0)

    systems = Any[upper, tail, headrace, penstock, turbine, rotor, generator, speed, governor]
    eqs = Equation[
        connect(upper.port, headrace.port_a),
        connect(headrace.port_b, penstock.port_a),
        connect(penstock.port_b, turbine.port_in),
        connect(turbine.port_out, tail.port),
        connect(governor.gate_out, turbine.opening),
        connect(turbine.shaft, rotor.flange_turbine),
        connect(rotor.flange_generator, generator.flange, speed.flange),
        connect(speed.w, governor.speed_in),
        generator.P_load.u ~ ifelse(t < TSTEP, P0, P1),
    ]

    surge = nothing
    if with_surge
        @named surge_comp = SurgeTank(A_t=ATANK, Z_0=ZSURGE0, D_riser=D_RISER, L_riser=L_RISER,
                                      e_riser=E_RISER, rho=RHO, g=G, p_atm=P_ATM)
        surge = surge_comp
        push!(systems, surge)
        deleteat!(eqs, 2)
        insert!(eqs, 2, connect(headrace.port_b, penstock.port_a, surge.port))
    end

    raw = ODESystem(eqs, t; name=with_surge ? :TrollheimWithSurge : :TrollheimNoSurge,
                    systems=systems)
    sys = structural_simplify(raw)
    (; sys, raw, headrace, penstock, turbine, rotor, generator, governor, surge)
end

function make_guesses(m; with_surge::Bool)
    g = Dict(
        m.headrace.dm => DM0,
        m.headrace.p_avg => P_UP,
        m.headrace.Re => HRG.Re,
        m.headrace.f_D => HRG.fD,
        m.headrace.dp_f => HRG.dpf,
        m.penstock.dm => DM0,
        m.penstock.p_avg => P_UP,
        m.penstock.Re => PNG.Re,
        m.penstock.f_D => PNG.fD,
        m.penstock.dp_f => PNG.dpf,
        m.turbine.H => H_GROSS,
        m.turbine.Q => Q0,
        m.turbine.eta => ETA0,
        m.turbine.P_mech => P0,
        m.turbine.dm => DM0,
        m.turbine.tau_shaft => TAU_SHAFT0,
        m.rotor.omega => W0,
        m.rotor.theta => 0.0,
        m.governor.xi => 0.0,
        m.governor.gate => TAU0,
    )
    with_surge && (g[m.surge.Z] = ZSURGE0)
    g
end

function simulate_case(; with_surge::Bool)
    m = build_trollheim(with_surge=with_surge)

    # Important: these are operating-point guesses, not hard initial constraints.
    # Let MTK solve a consistent nonlinear initialization point.
    u0 = Dict{Any,Any}()
    guesses = make_guesses(m; with_surge=with_surge)

    println("  building initialization problem with operating-point guesses only")
    try
        iprob = ModelingToolkit.InitializationProblem(m.sys, -300.0, u0; guesses=guesses)
        @printf("  initialization unknowns=%d\n", length(iprob.u0))
    catch err
        println("  initialization diagnostic construction failed: ", sprint(showerror, err))
    end

    prob = ODEProblem(m.sys, u0, (-300.0, 65.0); guesses=guesses)
    try
        if isdefined(ModelingToolkit, :analyze_initialization_jacobian)
            diag = ModelingToolkit.analyze_initialization_jacobian(prob)
            println("  initialization Jacobian diagnostic: ", diag)
        end
    catch err
        println("  initialization Jacobian diagnostic failed: ", sprint(showerror, err))
    end

    sol = solve(prob, Rodas5P(); tstops=[TSTEP], abstol=1e-8, reltol=1e-8, saveat=0.01)

    t_end = isempty(sol.t) ? -Inf : sol.t[end]
    @printf("  solver retcode=%s; saved points=%d; t_end=%.6f s\n", string(sol.retcode), length(sol.t), t_end)
    if !isempty(sol.t)
        f_end = sol[m.rotor.omega][end] / W0 * F0
        gate_end = sol[m.governor.tau_o][end]
        @printf("  terminal diagnostic: f=%.6f Hz; gate=%.6f pu\n", f_end, gate_end)
    end
    t_end >= 64.999 || error("Simulation terminated before 65 s (retcode=$(sol.retcode), t_end=$t_end)")

    (; m, sol)
end

function extract(run)
    m, sol = run.m, run.sol
    tv = sol.t
    f = sol[m.rotor.omega] ./ W0 .* F0
    pm = sol[m.turbine.P_mech] ./ 1e6
    pe = sol[m.generator.P_elec] ./ 1e6
    gate = sol[m.governor.tau_o]
    q = sol[m.turbine.Q]
    pre = findall(tt -> 0 <= tt < TSTEP, tv)
    post = findall(>=(TSTEP), tv)
    isempty(pre) && error("No saved pre-disturbance samples in [0,5) s")
    isempty(post) && error("No saved post-disturbance samples at/after 5 s")
    i = post[argmin(f[post])]
    iae = 0.0
    for k in 2:length(tv)
        tv[k] >= TSTEP || continue
        iae += 0.5 * (abs(f[k] - F0) + abs(f[k-1] - F0)) * (tv[k] - tv[k-1])
    end
    settle = NaN
    for j in post
        if all(abs.(f[j:end] .- F0) .<= 0.05)
            settle = tv[j]
            break
        end
    end
    base = (t=tv, f=f, pm=pm, pe=pe, gate=gate, q=q, nadir=f[i], t_nadir=tv[i],
            final_f=f[end], final_pm=pm[end], final_pe=pe[end], iae=iae, settling=settle,
            pre_f_range=maximum(f[pre]) - minimum(f[pre]),
            pre_pm_range=maximum(pm[pre]) - minimum(pm[pre]))
    m.surge === nothing ? base : merge(base, (Z=sol[m.surge.Z], Qs=sol[m.surge.port.dm] ./ RHO))
end

function write_csv(path, r; surge=false)
    open(path, "w") do io
        println(io, surge ? "time_s,frequency_Hz,mechanical_power_MW,electrical_power_MW,gate_pu,turbine_flow_m3s,surge_level_m,surge_flow_m3s" :
                            "time_s,frequency_Hz,mechanical_power_MW,electrical_power_MW,gate_pu,turbine_flow_m3s")
        for i in eachindex(r.t)
            r.t[i] >= 0 || continue
            if surge
                @printf(io, "%.6f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                        r.t[i], r.f[i], r.pm[i], r.pe[i], r.gate[i], r.q[i], r.Z[i], r.Qs[i])
            else
                @printf(io, "%.6f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                        r.t[i], r.f[i], r.pm[i], r.pe[i], r.gate[i], r.q[i])
            end
        end
    end
end

println("=== HydroPowerDynamics.jl Trollheim AGC surge comparison ===")
@printf("Electrical load %.1f -> %.1f MW at t=%.1f s; H=%.2f s; R=%.2f; Ki=%.2f 1/s\n",
        P0/1e6, P1/1e6, TSTEP, H_INERTIA, R_PU, KI)
println("Solving WITHOUT surge tank...")
r0 = extract(simulate_case(with_surge=false))
println("Solving WITH surge tank...")
r1 = extract(simulate_case(with_surge=true))

@printf("%-25s %14s %14s\n", "Metric", "No surge", "With surge")
@printf("%-25s %14.6f %14.6f\n", "Pre f range [Hz]", r0.pre_f_range, r1.pre_f_range)
@printf("%-25s %14.6f %14.6f\n", "Pre Pm range [MW]", r0.pre_pm_range, r1.pre_pm_range)
@printf("%-25s %14.6f %14.6f\n", "Nadir [Hz]", r0.nadir, r1.nadir)
@printf("%-25s %14.3f %14.3f\n", "Nadir time [s]", r0.t_nadir, r1.t_nadir)
@printf("%-25s %14.6f %14.6f\n", "f(65) [Hz]", r0.final_f, r1.final_f)
@printf("%-25s %14.6f %14.6f\n", "Pm(65) [MW]", r0.final_pm, r1.final_pm)
@printf("%-25s %14.6f %14.6f\n", "Pe(65) [MW]", r0.final_pe, r1.final_pe)
@printf("%-25s %14.6f %14.6f\n", "IAE [Hz s]", r0.iae, r1.iae)
@printf("%-25s %14.3f %14.3f\n", "Settling [s]", r0.settling, r1.settling)
@printf("Nadir delta (with - without) = %.6f Hz\n", r1.nadir - r0.nadir)

for ts in (0.0, 4.9, 5.0, 5.1, 6.0, 7.0, 10.0, 20.0, 40.0, 65.0)
    i0 = argmin(abs.(r0.t .- ts))
    i1 = argmin(abs.(r1.t .- ts))
    @printf("SAMPLE,%.1f,no_surge,%.8f,%.6f,%.6f,%.6f\n", ts, r0.f[i0], r0.pm[i0], r0.gate[i0], r0.q[i0])
    @printf("SAMPLE,%.1f,with_surge,%.8f,%.6f,%.6f,%.6f,%.6f,%.6f\n", ts, r1.f[i1], r1.pm[i1], r1.gate[i1], r1.q[i1], r1.Z[i1], r1.Qs[i1])
end

outdir = joinpath(@__DIR__, "results")
mkpath(outdir)
write_csv(joinpath(outdir, "without_surgetank.csv"), r0)
write_csv(joinpath(outdir, "with_surgetank.csv"), r1; surge=true)
open(joinpath(outdir, "metrics.csv"), "w") do io
    println(io, "metric,without_surgetank,with_surgetank")
    @printf(io, "nadir_Hz,%.10f,%.10f\n", r0.nadir, r1.nadir)
    @printf(io, "nadir_time_s,%.10f,%.10f\n", r0.t_nadir, r1.t_nadir)
    @printf(io, "final_frequency_Hz,%.10f,%.10f\n", r0.final_f, r1.final_f)
    @printf(io, "final_mechanical_power_MW,%.10f,%.10f\n", r0.final_pm, r1.final_pm)
    @printf(io, "IAE_Hz_s,%.10f,%.10f\n", r0.iae, r1.iae)
    @printf(io, "settling_s,%.10f,%.10f\n", r0.settling, r1.settling)
end
println("RESULT_DIR=", outdir)
