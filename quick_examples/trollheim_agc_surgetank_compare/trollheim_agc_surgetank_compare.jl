# =============================================================================
# Trollheim AGC comparison: with vs without surge tank
# HydroPowerDynamics.jl / ModelingToolkit.jl
#
# Same plant, same governor, same 75 -> 90 MW load step at t = 5 s.
# The ONLY topology change is the surge-tank branch at the headrace/penstock
# junction. A preconditioning interval (-300...0 s) is used so both systems
# reach their nonlinear hydraulic/control operating point before the visible test.
# =============================================================================

using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq
using Plots
using Printf

const PBASE      = 150e6
const P0         = 75e6
const P1         = 90e6
const TSTEP      = 5.0
const F0         = 50.0
const POLES      = 16                    # 375 rpm at 50 Hz
const RPM0       = 120F0 / POLES
const W0         = 2π * RPM0 / 60
const H_GROSS    = 371.0
const RHO        = 1000.0
const G          = 9.81
const Q_RATED    = 37.0
const D_RUNNER   = 2.5
const ETA_MAX    = 0.97
const C_ETA      = 0.25
const TAU_RATED  = 0.90
const KQ         = Q_RATED / (TAU_RATED * D_RUNNER^2 * sqrt(H_GROSS))

# Gate/flow estimate giving ~75 MW at nominal head before pipe losses.
const TAU0       = 0.5384797177203611
const Q0         = 22.13749950628151
const DM0        = RHO * Q0

# Match the earlier Trollheim AGC benchmark approximately: H ~= 1.83 s.
const H_INERTIA  = 1.83
const JTOTAL     = 2H_INERTIA * PBASE / W0^2

# Conservative hydro governor/AGC tuning carried from the OpenHPL study.
# GGOV1Governor uses R_droop in physical units [rad s / W].
const R_PU       = 0.50
const R_DROOP    = R_PU * W0 / PBASE
const TP         = 0.05
const TGOV       = 0.30
const KI         = 0.10 / PBASE

# Same hydraulic geometry in both cases.
const L_HEADRACE = 500.0
const D_HEADRACE = 6.0
const L_PENSTOCK = 500.0
const D_PENSTOCK = 4.0
const ROUGHNESS  = 1.5e-5

# Surge branch geometry.
const ATANK      = π/4 * 3.4^2
const ZSURGE0    = 69.97
const D_RISER    = 3.4
const L_RISER    = 87.0
const E_RISER    = 1.0e-2

"""Build one Trollheim plant. `with_surge=true` adds only the surge branch."""
function build_trollheim(; with_surge::Bool)
    @named upper = Reservoir(H=H_GROSS, rho=RHO, g=G, p_atm=101325.0)
    @named tail  = Reservoir(H=0.0,     rho=RHO, g=G, p_atm=101325.0)

    @named headrace = Penstock(L=L_HEADRACE, D_pipe=D_HEADRACE,
                               rho=RHO, B=2.0e9, e=ROUGHNESS)
    @named penstock = Penstock(L=L_PENSTOCK, D_pipe=D_PENSTOCK,
                               rho=RHO, B=2.0e9, e=ROUGHNESS)

    @named turbine = FrancisTurbineAffinity(D=D_RUNNER, K_q=KQ,
                                            eta_max=ETA_MAX, c_eta=C_ETA,
                                            Q_rated=Q_RATED, rho=RHO, g=G)
    @named rotor = RotorInertia(J=JTOTAL, b_v=0.0, omega_0=W0)
    @named generator = SimpleGenerator(P_rated=P0, omega_rated=W0,
                                       omega_s=W0, D_d=0.0, eta_gen=1.0)
    @named speed = RotationalSpeedSensor()
    @named pshaft = MechanicalPowerSensor()
    @named governor = GGOV1Governor(R_droop=R_DROOP, T_p=TP, T_gov=TGOV,
                                    omega_ref=W0, P_set=P0,
                                    tau_min=0.05, tau_max=1.0,
                                    K_i=KI, tau_0=TAU0)

    systems = Any[upper, tail, headrace, penstock, turbine, rotor,
                  generator, speed, pshaft, governor]

    eqs = Equation[
        connect(upper.port, headrace.port_a),
        connect(headrace.port_b, penstock.port_a),
        connect(penstock.port_b, turbine.port_in),
        connect(turbine.port_out, tail.port),
        connect(governor.gate_out, turbine.opening),
        connect(turbine.shaft, pshaft.flange_a),
        connect(pshaft.flange_b, rotor.flange_turbine),
        connect(rotor.flange_generator, generator.flange, speed.flange),
        connect(speed.w, governor.speed_in),
        connect(pshaft.P, governor.power_in),
    ]

    surge = nothing
    if with_surge
        @named surge_comp = SurgeTank(A_t=ATANK, Z_0=ZSURGE0,
                                      D_riser=D_RISER, L_riser=L_RISER,
                                      e_riser=E_RISER, rho=RHO, g=G,
                                      p_atm=101325.0)
        surge = surge_comp
        push!(systems, surge)
        # Replace the ordinary two-port junction by a three-way acausal node.
        deleteat!(eqs, 2)
        insert!(eqs, 2, connect(headrace.port_b, penstock.port_a, surge.port))
    end

    model_name = with_surge ? :TrollheimWithSurge : :TrollheimNoSurge
    raw = ODESystem(eqs, t; name=model_name, systems=systems)
    sys = structural_simplify(raw)
    return (; sys, raw, upper, tail, headrace, penstock, turbine, rotor,
            generator, speed, pshaft, governor, surge)
end

function simulate_case(; with_surge::Bool)
    m = build_trollheim(with_surge=with_surge)

    # Starts are deliberately treated as guesses; the hidden -300...0 s interval
    # is the operating-point conditioner.
    u0 = Dict(
        m.headrace.dm       => DM0,
        m.penstock.dm       => DM0,
        m.rotor.omega       => W0,
        m.governor.P_meas   => P0,
        m.governor.x_gov    => 0.0,
        m.governor.x_int    => TAU0,
    )
    if with_surge
        u0[m.surge.Z] = ZSURGE0
    end

    p0 = Dict(m.generator.P_rated => P0)
    prob = ODEProblem(m.sys, merge(u0, p0), (-300.0, 65.0))

    load_step = DiscreteCallback(
        (u, tt, integrator) -> tt == TSTEP,
        integrator -> (integrator.ps[m.generator.P_rated] = P1; nothing),
        save_positions=(true, true),
    )

    sol = solve(prob, Rodas5P(); callback=load_step, tstops=[TSTEP],
                abstol=1e-8, reltol=1e-8, saveat=0.01)
    return m, sol
end

function metrics(m, sol)
    tv    = sol.t
    omega = sol[m.rotor.omega]
    freq  = omega ./ W0 .* F0
    pm    = sol[m.turbine.P_mech] ./ 1e6
    gate  = sol[m.governor.tau_o]
    q     = sol[m.turbine.Q]

    pre  = findall(tt -> 0.0 <= tt < TSTEP, tv)
    post = findall(>=(TSTEP), tv)
    i_nadir = post[argmin(freq[post])]

    # Integral of absolute frequency error after the disturbance.
    e = abs.(freq .- F0)
    iae = 0.0
    for k in 2:length(tv)
        if tv[k] >= TSTEP
            dt = tv[k] - tv[k-1]
            iae += 0.5 * (e[k] + e[k-1]) * dt
        end
    end

    # First instant after which frequency remains within ±0.05 Hz.
    settle = NaN
    for i in post
        if all(abs.(freq[i:end] .- F0) .<= 0.05)
            settle = tv[i]
            break
        end
    end

    out = (
        t=tv, f=freq, pm=pm, gate=gate, q=q,
        nadir=freq[i_nadir], t_nadir=tv[i_nadir],
        final_f=freq[end], final_pm=pm[end], iae=iae, settling=settle,
        pre_f_range=maximum(freq[pre]) - minimum(freq[pre]),
        pre_pm_range=maximum(pm[pre]) - minimum(pm[pre]),
    )
    if m.surge !== nothing
        return merge(out, (Z=sol[m.surge.Z], Qs=sol[m.surge.port.dm] ./ RHO))
    end
    return out
end

println("\n=== HydroPowerDynamics.jl: Trollheim AGC / surge-tank comparison ===")
@printf "Load: %.1f -> %.1f MW at t=%.1f s\n" P0/1e6 P1/1e6 TSTEP
@printf "Initial gate estimate = %.6f pu; initial Q estimate = %.3f m3/s\n" TAU0 Q0
@printf "Equivalent inertia H = %.2f s; J = %.3e kg m2\n\n" H_INERTIA JTOTAL

println("Solving WITHOUT surge tank...")
m0, sol0 = simulate_case(with_surge=false)
r0 = metrics(m0, sol0)
println("Solving WITH surge tank...")
m1, sol1 = simulate_case(with_surge=true)
r1 = metrics(m1, sol1)

@printf "%-25s %14s %14s\n" "Metric" "No surge" "With surge"
@printf "%-25s %14.6f %14.6f\n" "Pre f range [Hz]" r0.pre_f_range r1.pre_f_range
@printf "%-25s %14.6f %14.6f\n" "Pre Pm range [MW]" r0.pre_pm_range r1.pre_pm_range
@printf "%-25s %14.5f %14.5f\n" "Nadir [Hz]" r0.nadir r1.nadir
@printf "%-25s %14.3f %14.3f\n" "Nadir time [s]" r0.t_nadir r1.t_nadir
@printf "%-25s %14.5f %14.5f\n" "f(65) [Hz]" r0.final_f r1.final_f
@printf "%-25s %14.4f %14.4f\n" "Pm(65) [MW]" r0.final_pm r1.final_pm
@printf "%-25s %14.5f %14.5f\n" "IAE [Hz s]" r0.iae r1.iae
@printf "%-25s %14.3f %14.3f\n" "Settling [s]" r0.settling r1.settling
@printf "\nNadir improvement with surge tank = %.5f Hz\n" (r1.nadir-r0.nadir)

mkpath(joinpath(@__DIR__, "plots"))
mask0 = r0.t .>= 0
mask1 = r1.t .>= 0

p_f = plot(r0.t[mask0], r0.f[mask0], ls=:dash, lw=2.2,
           label="Without surge tank", xlabel="Time [s]", ylabel="Frequency [Hz]",
           title="Trollheim AGC — frequency")
plot!(p_f, r1.t[mask1], r1.f[mask1], ls=:solid, lw=2.2, label="With surge tank")
vline!(p_f, [TSTEP], ls=:dot, label="75→90 MW step")
hline!(p_f, [F0], ls=:dot, label=false)

p_p = plot(r0.t[mask0], r0.pm[mask0], ls=:dash, lw=2.2,
           label="Without surge tank", xlabel="Time [s]", ylabel="Mechanical power [MW]",
           title="Trollheim AGC — turbine power")
plot!(p_p, r1.t[mask1], r1.pm[mask1], ls=:solid, lw=2.2, label="With surge tank")
hline!(p_p, [P1/1e6], ls=:dot, label="90 MW load")
vline!(p_p, [TSTEP], ls=:dot, label=false)

p_g = plot(r0.t[mask0], r0.gate[mask0], ls=:dash, lw=2.2,
           label="Without surge tank", xlabel="Time [s]", ylabel="Gate [pu]",
           title="Trollheim AGC — guide vane")
plot!(p_g, r1.t[mask1], r1.gate[mask1], ls=:solid, lw=2.2, label="With surge tank")
vline!(p_g, [TSTEP], ls=:dot, label=false)

p_q = plot(r0.t[mask0], r0.q[mask0], ls=:dash, lw=2.2,
           label="Without surge tank", xlabel="Time [s]", ylabel="Turbine flow [m3/s]",
           title="Trollheim AGC — turbine flow")
plot!(p_q, r1.t[mask1], r1.q[mask1], ls=:solid, lw=2.2, label="With surge tank")
vline!(p_q, [TSTEP], ls=:dot, label=false)

savefig(p_f, joinpath(@__DIR__, "plots", "frequency_compare.png"))
savefig(p_p, joinpath(@__DIR__, "plots", "power_compare.png"))
savefig(p_g, joinpath(@__DIR__, "plots", "gate_compare.png"))
savefig(p_q, joinpath(@__DIR__, "plots", "flow_compare.png"))

p_st = plot(r1.t[mask1], r1.Z[mask1], lw=2.2, xlabel="Time [s]",
            ylabel="Surge level [m]", label="Z", title="Surge-tank state")
vline!(p_st, [TSTEP], ls=:dot, label=false)
savefig(p_st, joinpath(@__DIR__, "plots", "surge_level.png"))

p_sf = plot(r1.t[mask1], r1.Qs[mask1], lw=2.2, xlabel="Time [s]",
            ylabel="Surge flow [m3/s]", label="Qs", title="Surge-tank branch flow")
hline!(p_sf, [0.0], ls=:dot, label=false)
vline!(p_sf, [TSTEP], ls=:dot, label=false)
savefig(p_sf, joinpath(@__DIR__, "plots", "surge_flow.png"))

println("\nPlots written to: ", joinpath(@__DIR__, "plots"))
