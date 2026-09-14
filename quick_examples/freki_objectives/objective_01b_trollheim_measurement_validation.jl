# FREKI Objective 1B — Trollheim HPP measurement validation
#
# Uses the existing processed Trollheim measurement record to calibrate a
# nonlinear HydroPowerDynamics.jl Francis-turbine relation on one time window
# and validate it on held-out plant data. This is a measurement-based evidence
# layer for FREKI Objective 1; it is not an official FCR prequalification test.

using CSV
using DataFrames
using Statistics
using LinearAlgebra
using Printf
using Plots

const RHO = 1000.0
const G = 9.81
const ETA_GEN = 0.985
const D_RUNNER = 2.5

const CAL_START = 1200.0
const CAL_END = 1800.0
const VAL_START = 1801.0
const VAL_END = 2400.0
const STARTUP_START = 850.0
const STARTUP_END = 1200.0

base = @__DIR__
data_path = joinpath(base, "..", "trollheim_experimental", "01_preprocessing", "trollheim_processed.csv")
plot_dir = joinpath(base, "plots")
result_dir = joinpath(base, "results")
mkpath(plot_dir)
mkpath(result_dir)

df = CSV.read(data_path, DataFrame)

println("=== FREKI Objective 1B: Trollheim measurement validation ===")
println("Dataset: $(nrow(df)) samples, 1 s sampling")
@printf("Calibration window: %.0f-%.0f s\n", CAL_START, CAL_END)
@printf("Held-out validation window: %.0f-%.0f s\n", VAL_START, VAL_END)
@printf("Independent startup check: %.0f-%.0f s\n\n", STARTUP_START, STARTUP_END)

validrow(r) = isfinite(r.H) && isfinite(r.Q) && isfinite(r.tau_o) &&
              isfinite(r.P_elec_MW) && r.H > 1.0 && r.Q > 1.0 && r.tau_o > 0.02

cal = filter(r -> CAL_START <= r.t <= CAL_END && validrow(r), df)
val = filter(r -> VAL_START <= r.t <= VAL_END && validrow(r), df)
startup = filter(r -> STARTUP_START <= r.t <= STARTUP_END && validrow(r) && r.P_elec_MW > 5.0, df)

isempty(cal) && error("No usable calibration samples")
isempty(val) && error("No usable validation samples")
isempty(startup) && error("No usable startup samples")

# -----------------------------------------------------------------------------
# Nonlinear Francis model used in HydroPowerDynamics.jl experimental workflow
#
#   Q = tau_o * Kq * D^2 * sqrt(H)
#   eta = eta_max * [1 - c_eta (Q/Qr - 1)^2]
#   P_elec = eta_gen * rho*g*Q*H*eta
#
# Kq, eta_max and c_eta are identified only from the calibration window.
# -----------------------------------------------------------------------------

fcal = Float64.(cal.tau_o) .* D_RUNNER^2 .* sqrt.(Float64.(cal.H))
qcal = Float64.(cal.Q)
Kq = sum(fcal .* qcal) / sum(fcal .^ 2)
Qr = mean(qcal)

eta_cal = Float64.(cal.P_mech) ./ (RHO * G .* qcal .* Float64.(cal.H))
xcal = qcal ./ Qr .- 1.0
A = hcat(ones(length(xcal)), -(xcal .^ 2))
coef = A \ eta_cal
eta_max = coef[1]
c_eta = coef[2] / eta_max

function predict(seg)
    H = Float64.(seg.H)
    gate = Float64.(seg.tau_o)
    qhat = Kq .* gate .* D_RUNNER^2 .* sqrt.(H)
    etahat = eta_max .* (1.0 .- c_eta .* (qhat ./ Qr .- 1.0).^2)
    etahat = clamp.(etahat, 0.0, 1.0)
    pmech_hat = RHO * G .* qhat .* H .* etahat
    pelec_hat_MW = ETA_GEN .* pmech_hat ./ 1e6
    return qhat, etahat, pelec_hat_MW
end

qhat_cal, etahat_cal, phat_cal = predict(cal)
qhat_val, etahat_val, phat_val = predict(val)
qhat_start, etahat_start, phat_start = predict(startup)

rmse(a,b) = sqrt(mean((a .- b).^2))
fitpct(y, yhat) = 100 * (1 - norm(yhat .- y) / max(norm(y .- mean(y)), eps()))

q_val = Float64.(val.Q)
p_val = Float64.(val.P_elec_MW)
q_start = Float64.(startup.Q)
p_start = Float64.(startup.P_elec_MW)

rmse_q_val = rmse(q_val, qhat_val)
rmse_p_val = rmse(p_val, phat_val)
rel_q_val = 100 * rmse_q_val / mean(q_val)
rel_p_val = 100 * rmse_p_val / mean(p_val)
fit_q_val = fitpct(q_val, qhat_val)
fit_p_val = fitpct(p_val, phat_val)
bias_q_val = mean(qhat_val .- q_val)
bias_p_val = mean(phat_val .- p_val)

rmse_q_start = rmse(q_start, qhat_start)
rmse_p_start = rmse(p_start, phat_start)
rel_q_start = 100 * rmse_q_start / mean(q_start)
rel_p_start = 100 * rmse_p_start / mean(p_start)

function traffic(relq, relp)
    if relq < 1.0 && relp < 1.0
        return "GREEN"
    elseif relq < 3.0 && relp < 3.0
        return "AMBER"
    else
        return "RED"
    end
end

turbine_health = traffic(rel_q_val, rel_p_val)
hydraulic_health = "AMBER"
governor_health = "AMBER"
fcr_readiness = "MORE DATA REQUIRED"

println("Identified from calibration data:")
@printf("  Kq       = %.6f\n", Kq)
@printf("  eta_max  = %.6f\n", eta_max)
@printf("  c_eta    = %.6f\n", c_eta)
@printf("  Qr       = %.4f m3/s\n\n", Qr)

println("Held-out validation:")
@printf("  Q RMSE   = %.4f m3/s  (%.3f %%)\n", rmse_q_val, rel_q_val)
@printf("  P RMSE   = %.4f MW    (%.3f %%)\n", rmse_p_val, rel_p_val)
@printf("  Q FIT    = %.2f %%\n", fit_q_val)
@printf("  P FIT    = %.2f %%\n", fit_p_val)
@printf("  Q bias   = %.4f m3/s\n", bias_q_val)
@printf("  P bias   = %.4f MW\n\n", bias_p_val)

println("Independent startup check:")
@printf("  Q RMSE   = %.4f m3/s  (%.2f %%)\n", rmse_q_start, rel_q_start)
@printf("  P RMSE   = %.4f MW    (%.2f %%)\n\n", rmse_p_start, rel_p_start)

println("Model-health decision:")
println("  TURBINE MODEL   = $turbine_health")
println("  HYDRAULIC MODEL = $hydraulic_health")
println("  GOVERNOR MODEL  = $governor_health")
println("  FCR READINESS   = $fcr_readiness")

summary = DataFrame(
    quantity = [
        "K_q", "eta_max", "c_eta", "Q_rated",
        "RMSE_Q_validation", "RMSE_P_validation",
        "REL_RMSE_Q_validation", "REL_RMSE_P_validation",
        "FIT_Q_validation", "FIT_P_validation",
        "RMSE_Q_startup", "RMSE_P_startup",
        "turbine_model_health", "hydraulic_model_health",
        "governor_model_health", "FCR_readiness"
    ],
    value = Any[
        Kq, eta_max, c_eta, Qr,
        rmse_q_val, rmse_p_val,
        rel_q_val, rel_p_val,
        fit_q_val, fit_p_val,
        rmse_q_start, rmse_p_start,
        turbine_health, hydraulic_health,
        governor_health, fcr_readiness
    ],
    unit = [
        "-", "-", "-", "m3/s",
        "m3/s", "MW", "%", "%", "%", "%",
        "m3/s", "MW", "state", "state", "state", "state"
    ]
)
CSV.write(joinpath(result_dir, "objective_01b_summary.csv"), summary)

validation_ts = DataFrame(
    t = Float64.(val.t),
    H_measured_m = Float64.(val.H),
    gate_measured_pu = Float64.(val.tau_o),
    Q_measured_m3s = q_val,
    Q_model_m3s = qhat_val,
    P_measured_MW = p_val,
    P_model_MW = phat_val,
    Q_residual_m3s = qhat_val .- q_val,
    P_residual_MW = phat_val .- p_val,
    speed_measured_rpm = Float64.(val.n_rpm)
)
CSV.write(joinpath(result_dir, "objective_01b_validation_timeseries.csv"), validation_ts)

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
default(size=(900, 500), linewidth=1.6, legend=:best)

p1 = plot(Float64.(df.t), Float64.(df.P_elec_MW), label="Generator power", ylabel="MW")
plot!(p1, Float64.(df.t), 100 .* Float64.(df.tau_o), label="Gate x100", ylabel="MW / %")
vspan!(p1, [CAL_START, CAL_END], alpha=0.08, label="calibration")
vspan!(p1, [VAL_START, VAL_END], alpha=0.08, label="held-out")
xlabel!(p1, "Time [s]")
title!(p1, "Trollheim measurement record and validation windows")
savefig(p1, joinpath(plot_dir, "01b_01_measurement_windows.png"))

p2 = plot(Float64.(cal.t), qcal, label="Measured Q", linestyle=:dash, ylabel="Q [m3/s]")
plot!(p2, Float64.(cal.t), qhat_cal, label="Calibrated nonlinear model")
xlabel!(p2, "Time [s]")
title!(p2, "Calibration window: measured and modelled flow")
savefig(p2, joinpath(plot_dir, "01b_02_calibration_flow.png"))

p3 = plot(Float64.(val.t), q_val, label="Measured Q", linestyle=:dash, ylabel="Q [m3/s]")
plot!(p3, Float64.(val.t), qhat_val, label="HPD nonlinear turbine replay")
xlabel!(p3, "Time [s]")
title!(p3, @sprintf("Held-out flow validation — RMSE %.3f m3/s", rmse_q_val))
savefig(p3, joinpath(plot_dir, "01b_03_heldout_flow.png"))

p4 = plot(Float64.(val.t), p_val, label="Measured P", linestyle=:dash, ylabel="P_elec [MW]")
plot!(p4, Float64.(val.t), phat_val, label="HPD nonlinear turbine replay")
xlabel!(p4, "Time [s]")
title!(p4, @sprintf("Held-out power validation — RMSE %.3f MW", rmse_p_val))
savefig(p4, joinpath(plot_dir, "01b_04_heldout_power.png"))

p5 = plot(Float64.(val.t), qhat_val .- q_val, label="Q residual [m3/s]", ylabel="Residual")
plot!(p5, Float64.(val.t), phat_val .- p_val, label="P residual [MW]")
hline!(p5, [0.0], linestyle=:dash, label="zero")
xlabel!(p5, "Time [s]")
title!(p5, "Held-out validation residuals")
savefig(p5, joinpath(plot_dir, "01b_05_residuals.png"))

p6 = plot(Float64.(startup.t), p_start, label="Measured P", linestyle=:dash, ylabel="P_elec [MW]")
plot!(p6, Float64.(startup.t), phat_start, label="Model replay")
xlabel!(p6, "Time [s]")
title!(p6, @sprintf("Independent startup check — power RMSE %.2f MW", rmse_p_start))
savefig(p6, joinpath(plot_dir, "01b_06_startup_power.png"))

health_labels = ["Turbine", "Hydraulic", "Governor", "FCR readiness"]
health_score = [
    turbine_health == "GREEN" ? 3.0 : turbine_health == "AMBER" ? 2.0 : 1.0,
    2.0,
    2.0,
    1.0
]
p7 = bar(health_labels, health_score, ylim=(0,3.5), legend=false, ylabel="Evidence level")
hline!(p7, [1.0,2.0,3.0], linestyle=:dot)
title!(p7, "Objective 1B model-health evidence")
savefig(p7, joinpath(plot_dir, "01b_07_model_health.png"))

println("\nGenerated:")
println("  results/objective_01b_summary.csv")
println("  results/objective_01b_validation_timeseries.csv")
println("  plots/01b_01_measurement_windows.png ... 01b_07_model_health.png")
