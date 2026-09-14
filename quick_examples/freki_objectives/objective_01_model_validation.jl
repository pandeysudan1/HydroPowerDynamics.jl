# FREKI Objective 1: validate an HPD plant model from normal-operation data.
#
# The nonlinear Trollheim-inspired plant generates synthetic measurements under
# small load variations. Governor droop R and servo time constant T_g are then
# identified from the measured frequency and guide-vane motion. The identified
# parameters are replayed on the nonlinear HydroPowerDynamics.jl plant and
# checked on a held-out validation interval.

using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq
using Statistics
using Random
using LinearAlgebra
using Printf
using CSV
using DataFrames
using Plots

const PBASE = 150e6
const P0 = 75e6
const F0 = 50.0
const POLES = 16
const W0 = 2π * (120F0 / POLES) / 60
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
const GATE0 = 0.5384797177203611
const Q0 = 22.13749950628151
const DM0 = RHO * Q0
const H_INERTIA = 1.83
const JTOTAL = 2 * H_INERTIA * PBASE / W0^2
const KI = 0.10
const R_TRUE = 0.50
const TG_TRUE = 0.30
const L_HEADRACE = 500.0
const D_HEADRACE = 6.0
const L_PENSTOCK = 500.0
const D_PENSTOCK = 4.0
const ROUGHNESS = 1.5e-5
const ETA0 = ETA_MAX * (1 - C_ETA * (Q0 / Q_RATED - 1)^2)
const T_END = 65.0
const T_CAL_END = 35.0
const DT = 0.05

# Small deterministic variations represent ordinary plant operation, not a
# dedicated prequalification step test. At t=0 the perturbation is zero.
function normal_load(tt)
    tt < 0 && return P0
    return P0 +
        0.55e6 * sin(2π * tt / 17.0) +
        0.30e6 * sin(2π * tt / 7.5) +
        0.15e6 * sin(2π * tt / 31.0)
end

function pipe_guess(L, Dpipe)
    A = π * Dpipe^2 / 4
    μ = 1e-3
    Re = DM0 * Dpipe / (μ * A)
    fD = 1 / (2log10(ROUGHNESS / (3.7Dpipe) + 5.7 / Re^0.9))^2
    dpf = fD * (L / Dpipe) * (DM0 / (RHO * A))^2 * RHO / 2
    (; Re, fD, dpf)
end

const HRG = pipe_guess(L_HEADRACE, D_HEADRACE)
const PNG = pipe_guess(L_PENSTOCK, D_PENSTOCK)

function build_model(; R=R_TRUE, Tg=TG_TRUE)
    @named upper = Reservoir(H=H_GROSS, rho=RHO, g=G, p_atm=P_ATM)
    @named tail = Reservoir(H=0.0, rho=RHO, g=G, p_atm=P_ATM)
    @named headrace = Penstock(L=L_HEADRACE, D_pipe=D_HEADRACE, rho=RHO, B=2e9, e=ROUGHNESS)
    @named penstock = Penstock(L=L_PENSTOCK, D_pipe=D_PENSTOCK, rho=RHO, B=2e9, e=ROUGHNESS)
    @named turbine = FrancisTurbineAffinity(D=D_RUNNER, K_q=KQ, eta_max=ETA_MAX,
        c_eta=C_ETA, Q_rated=Q_RATED, rho=RHO, g=G)
    @named rotor = RotorInertia(J=JTOTAL, b_v=0.0, omega_0=W0)
    @named generator = VariableLoadGenerator(omega_rated=W0, omega_s=W0, D_d=0.0, eta_gen=1.0)
    @named speed = RotationalSpeedSensor()
    @named governor = SimpleGovernorAGC(R=R, T_g=Tg, K_i=KI, omega_ref=W0,
        gate_bias=GATE0, gate_min=0.05, gate_max=1.0)

    eqs = [
        connect(upper.port, headrace.port_a),
        connect(headrace.port_b, penstock.port_a),
        connect(penstock.port_b, turbine.port_in),
        connect(turbine.port_out, tail.port),
        connect(governor.gate_out, turbine.opening),
        connect(turbine.shaft, rotor.flange_turbine),
        connect(rotor.flange_generator, generator.flange, speed.flange),
        connect(speed.w, governor.speed_in),
        generator.P_load.u ~ ifelse(t < 0.0, P0,
            P0 + 0.55e6*sin(2π*t/17.0) + 0.30e6*sin(2π*t/7.5) + 0.15e6*sin(2π*t/31.0)),
    ]

    raw = ODESystem(eqs, t; name=:FREKIObjective1,
        systems=[upper, tail, headrace, penstock, turbine, rotor, generator, speed, governor])
    sys = structural_simplify(raw)
    (; sys, headrace, penstock, turbine, rotor, generator, governor)
end

function simulate(; R=R_TRUE, Tg=TG_TRUE)
    m = build_model(R=R, Tg=Tg)
    u0 = Dict(
        m.headrace.dm => DM0,
        m.headrace.p_avg => P_UP,
        m.penstock.dm => DM0,
        m.penstock.p_avg => P_UP,
        m.rotor.omega => W0,
        m.rotor.theta => 0.0,
        m.governor.xi => 0.0,
        m.governor.gate => GATE0,
    )
    guesses = Dict(
        m.headrace.Re => HRG.Re,
        m.headrace.f_D => HRG.fD,
        m.headrace.dp_f => HRG.dpf,
        m.penstock.Re => PNG.Re,
        m.penstock.f_D => PNG.fD,
        m.penstock.dp_f => PNG.dpf,
        m.turbine.H => H_GROSS,
        m.turbine.Q => Q0,
        m.turbine.eta => ETA0,
        m.turbine.P_mech => P0,
        m.turbine.dm => DM0,
        m.turbine.tau_shaft => P0 / W0,
    )
    prob = ODEProblem(m.sys, u0, (-300.0, T_END); guesses)
    sol = solve(prob, Rodas5P(); abstol=1e-8, reltol=1e-8, saveat=DT)
    sol.t[end] >= T_END - 1e-6 || error("Simulation stopped early: $(sol.retcode)")
    keep = findall(>=(0.0), sol.t)
    tv = sol.t[keep]
    (; t=tv,
       f=sol[m.rotor.omega][keep] ./ W0 .* F0,
       pm=sol[m.turbine.P_mech][keep] ./ 1e6,
       pe=sol[m.generator.P_elec][keep] ./ 1e6,
       gate=sol[m.governor.tau_o][keep],
       q=sol[m.turbine.Q][keep])
end

function moving_average(x, halfwindow)
    y = similar(x)
    for i in eachindex(x)
        lo = max(firstindex(x), i-halfwindow)
        hi = min(lastindex(x), i+halfwindow)
        y[i] = mean(@view x[lo:hi])
    end
    y
end

function cumulative_trapezoid(tvec, x)
    z = zeros(length(x))
    for i in 2:length(x)
        z[i] = z[i-1] + 0.5*(x[i] + x[i-1])*(tvec[i]-tvec[i-1])
    end
    z
end

function estimate_governor(tvec, fmeas, gatemeas)
    # Smooth only for derivative estimation. The raw measurements remain in the
    # reported validation metrics.
    gs = moving_average(gatemeas, 6)
    dg = zeros(length(gs))
    for i in 2:length(gs)-1
        dg[i] = (gs[i+1] - gs[i-1]) / (tvec[i+1] - tvec[i-1])
    end
    dg[1] = dg[2]
    dg[end] = dg[end-1]

    ef = (F0 .- fmeas) ./ F0
    xi = cumulative_trapezoid(tvec, ef)
    x1 = GATE0 .+ KI .* xi .- gs
    x2 = ef

    idx = findall(i -> 2.0 <= tvec[i] <= T_CAL_END, eachindex(tvec))
    X = hcat(x1[idx], x2[idx])
    β = X \ dg[idx]
    a, b = β
    Tg_hat = 1 / a
    R_hat = a / b
    (; R_hat, Tg_hat, beta=β, gate_smooth=gs, dgate=dg, ef, xi, calibration_idx=idx)
end

rmse(a, b) = sqrt(mean((a .- b).^2))
function fit_percent(y, yhat)
    den = norm(y .- mean(y))
    den < eps() && return NaN
    100 * (1 - norm(y .- yhat) / den)
end

println("=== FREKI Objective 1: normal-operation model validation ===")
@printf("Truth: R=%.4f pu/pu, Tg=%.4f s; calibration 0-%.1f s; validation %.1f-%.1f s\n",
    R_TRUE, TG_TRUE, T_CAL_END, T_CAL_END, T_END)
println("Solving nonlinear HPD truth model...")
truth = simulate()

Random.seed!(240914)
σf = 0.0020       # Hz
σp = 0.020        # MW
σg = 0.00020      # pu
σq = 0.0030       # m3/s
f_meas = truth.f .+ σf .* randn(length(truth.f))
pe_meas = truth.pe .+ σp .* randn(length(truth.pe))
gate_meas = truth.gate .+ σg .* randn(length(truth.gate))
q_meas = truth.q .+ σq .* randn(length(truth.q))

est = estimate_governor(truth.t, f_meas, gate_meas)
@printf("Identified: R=%.6f pu/pu (error %.2f%%), Tg=%.6f s (error %.2f%%)\n",
    est.R_hat, 100*(est.R_hat/R_TRUE-1), est.Tg_hat, 100*(est.Tg_hat/TG_TRUE-1))

(0.2 < est.R_hat < 1.0) || error("Identified R outside physical screening range: $(est.R_hat)")
(0.05 < est.Tg_hat < 1.0) || error("Identified Tg outside physical screening range: $(est.Tg_hat)")

println("Replaying nonlinear HPD model with identified parameters...")
identified = simulate(R=est.R_hat, Tg=est.Tg_hat)
length(identified.t) == length(truth.t) || error("Truth and replay grids differ")

ival = findall(>(T_CAL_END), truth.t)
metrics = Dict(
    "RMSE_frequency_Hz" => rmse(f_meas[ival], identified.f[ival]),
    "RMSE_electrical_power_MW" => rmse(pe_meas[ival], identified.pe[ival]),
    "RMSE_gate_pu" => rmse(gate_meas[ival], identified.gate[ival]),
    "RMSE_flow_m3s" => rmse(q_meas[ival], identified.q[ival]),
    "FIT_frequency_pct" => fit_percent(f_meas[ival], identified.f[ival]),
    "FIT_electrical_power_pct" => fit_percent(pe_meas[ival], identified.pe[ival]),
)

@printf("Validation RMSE: f=%.6f Hz, Pe=%.6f MW, gate=%.8f pu, Q=%.6f m3/s\n",
    metrics["RMSE_frequency_Hz"], metrics["RMSE_electrical_power_MW"],
    metrics["RMSE_gate_pu"], metrics["RMSE_flow_m3s"])
@printf("Validation FIT: f=%.2f%%, Pe=%.2f%%\n",
    metrics["FIT_frequency_pct"], metrics["FIT_electrical_power_pct"])

outdir = joinpath(@__DIR__, "results")
plotdir = joinpath(@__DIR__, "plots")
mkpath(outdir); mkpath(plotdir)

load_MW = normal_load.(truth.t) ./ 1e6
res_f = f_meas .- identified.f
res_pe = pe_meas .- identified.pe

CSV.write(joinpath(outdir, "objective_01_timeseries.csv"), DataFrame(
    time_s=truth.t,
    load_MW=load_MW,
    measured_frequency_Hz=f_meas,
    model_frequency_Hz=identified.f,
    measured_electrical_power_MW=pe_meas,
    model_electrical_power_MW=identified.pe,
    measured_gate_pu=gate_meas,
    model_gate_pu=identified.gate,
    measured_flow_m3s=q_meas,
    model_flow_m3s=identified.q,
    frequency_residual_Hz=res_f,
    power_residual_MW=res_pe,
))

summary = DataFrame(
    quantity=["R", "Tg", "RMSE_frequency", "RMSE_electrical_power", "RMSE_gate", "RMSE_flow", "FIT_frequency", "FIT_electrical_power"],
    reference=[R_TRUE, TG_TRUE, NaN, NaN, NaN, NaN, NaN, NaN],
    identified_or_metric=[est.R_hat, est.Tg_hat,
        metrics["RMSE_frequency_Hz"], metrics["RMSE_electrical_power_MW"],
        metrics["RMSE_gate_pu"], metrics["RMSE_flow_m3s"],
        metrics["FIT_frequency_pct"], metrics["FIT_electrical_power_pct"]],
    unit=["pu/pu", "s", "Hz", "MW", "pu", "m3/s", "%", "%"],
)
CSV.write(joinpath(outdir, "objective_01_summary.csv"), summary)

p1 = plot(truth.t, load_MW, xlabel="Time [s]", ylabel="Load [MW]",
    title="Normal-operation excitation", legend=false)
savefig(p1, joinpath(plotdir, "01_normal_operation_input.png"))

p2 = plot(truth.t, f_meas, label="synthetic measurement", xlabel="Time [s]", ylabel="Frequency [Hz]",
    title="Frequency: measured vs identified HPD model")
plot!(p2, truth.t, identified.f, label="identified HPD model")
vline!(p2, [T_CAL_END], label="validation starts")
savefig(p2, joinpath(plotdir, "02_frequency_validation.png"))

p3 = plot(truth.t, pe_meas, label="synthetic measurement", xlabel="Time [s]", ylabel="Electrical power [MW]",
    title="Power: measured vs identified HPD model")
plot!(p3, truth.t, identified.pe, label="identified HPD model")
vline!(p3, [T_CAL_END], label="validation starts")
savefig(p3, joinpath(plotdir, "03_power_validation.png"))

p4 = plot(truth.t, gate_meas, label="measured gate", xlabel="Time [s]", ylabel="Gate [pu]",
    title="Guide-vane validation")
plot!(p4, truth.t, identified.gate, label="identified model")
savefig(p4, joinpath(plotdir, "04_gate_validation.png"))

p5 = plot(truth.t, res_f, label="frequency residual [Hz]", xlabel="Time [s]", ylabel="Residual",
    title="Held-out model residuals")
plot!(p5, truth.t, res_pe ./ 10, label="power residual / 10 [MW]")
vline!(p5, [T_CAL_END], label="validation starts")
savefig(p5, joinpath(plotdir, "05_validation_residuals.png"))

p6 = bar(["R", "Tg"], [est.R_hat/R_TRUE, est.Tg_hat/TG_TRUE],
    ylabel="identified / true", title="Identified parameter ratios", legend=false)
hline!(p6, [1.0], label="truth")
savefig(p6, joinpath(plotdir, "06_parameter_identification.png"))

println("RESULT_DIR=", outdir)
println("PLOT_DIR=", plotdir)
println("FREKI_OBJECTIVE_01=PASS")
