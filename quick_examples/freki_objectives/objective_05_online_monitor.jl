# FREKI Objective 5 — online-style model health and FCR capability monitor
#
# Replays the existing one-hour Trollheim measurement record sample-by-sample.
# This is NOT a live SCADA connection. It is an offline streaming replay used to
# test the FREKI logic:
#
#   measurements -> model prediction -> rolling residuals -> model health
#                -> model-based FCR screening capacity
#
# A numerical FCR screening value is reported continuously, but it is marked
# trusted only when the local measurement/model residual test is GREEN.
# Formal Statnett/Nordic qualification is not claimed.

using CSV
using DataFrames
using Statistics
using Printf
using Plots

const RHO = 1000.0
const G = 9.81
const ETA_GEN = 0.985
const D_RUNNER = 2.5
const GATE_MAX = 1.0
const WINDOW_S = 60

base = @__DIR__
data_path = joinpath(base, "..", "trollheim_experimental", "01_preprocessing", "trollheim_processed.csv")
result_dir = joinpath(base, "results")
plot_dir = joinpath(base, "plots")
mkpath(result_dir)
mkpath(plot_dir)

obj1b = CSV.read(joinpath(result_dir, "objective_01b_summary.csv"), DataFrame)
obj2 = CSV.read(joinpath(result_dir, "objective_02_summary.csv"), DataFrame)

df = CSV.read(data_path, DataFrame)

function lookup_num(tbl, key)
    i = findfirst(==(key), tbl.quantity)
    i === nothing && error("Missing $key")
    parse(Float64, string(tbl.value[i]))
end

KQ = lookup_num(obj1b, "K_q")
ETA_MAX = lookup_num(obj1b, "eta_max")
ROBUST_FCR_MW = lookup_num(obj2, "robust_screening_capacity_MW")

# The online monitor intentionally uses the locally supported turbine map only.
# c_eta from Objective 1B was weakly identified, so a constant local efficiency
# is used rather than extrapolating the unphysical fitted curvature.
function local_prediction(H, gate)
    q = gate * KQ * D_RUNNER^2 * sqrt(max(H, 0.0))
    p = ETA_GEN * RHO * G * q * max(H, 0.0) * ETA_MAX / 1e6
    return q, p
end

function static_headroom(H, gate)
    _, pnow = local_prediction(H, gate)
    _, pmax = local_prediction(H, GATE_MAX)
    max(0.0, pmax - pnow)
end

n = nrow(df)
qhat = fill(NaN, n)
phat = fill(NaN, n)
rmse_q = fill(NaN, n)
rmse_p = fill(NaN, n)
rel_q = fill(NaN, n)
rel_p = fill(NaN, n)
health = fill("OFFLINE", n)
static_headroom_MW = fill(NaN, n)
screening_fcr_MW = fill(NaN, n)
trusted_fcr_MW = fill(NaN, n)

for i in 1:n
    H = Float64(df.H[i])
    gate = Float64(df.tau_o[i])
    q = Float64(df.Q[i])
    p = Float64(df.P_elec_MW[i])

    active = isfinite(H) && isfinite(gate) && isfinite(q) && isfinite(p) &&
             H > 1.0 && gate > 0.05 && q > 5.0 && p > 5.0

    if active
        qhat[i], phat[i] = local_prediction(H, gate)
        static_headroom_MW[i] = static_headroom(H, gate)
        screening_fcr_MW[i] = min(ROBUST_FCR_MW, static_headroom_MW[i])
    end

    # Rolling residual health only after enough active samples are available.
    i0 = max(1, i - WINDOW_S + 1)
    idx = [j for j in i0:i if isfinite(qhat[j]) && isfinite(phat[j])]

    if length(idx) >= 30
        qmeas = Float64.(df.Q[idx])
        pmeas = Float64.(df.P_elec_MW[idx])
        qpred = qhat[idx]
        ppred = phat[idx]

        rmse_q[i] = sqrt(mean((qpred .- qmeas).^2))
        rmse_p[i] = sqrt(mean((ppred .- pmeas).^2))
        rel_q[i] = 100 * rmse_q[i] / max(mean(abs.(qmeas)), eps())
        rel_p[i] = 100 * rmse_p[i] / max(mean(abs.(pmeas)), eps())

        if rel_q[i] < 1.0 && rel_p[i] < 1.0
            health[i] = "GREEN"
        elseif rel_q[i] < 3.0 && rel_p[i] < 3.0
            health[i] = "AMBER"
        else
            health[i] = "RED"
        end

        # Trusted here means trusted as a local model-based screening signal only.
        # It does NOT mean formally qualified FCR.
        if health[i] == "GREEN" && isfinite(screening_fcr_MW[i])
            trusted_fcr_MW[i] = screening_fcr_MW[i]
        end
    end
end

online = DataFrame(
    t = Float64.(df.t),
    H_m = Float64.(df.H),
    gate_pu = Float64.(df.tau_o),
    Q_measured_m3s = Float64.(df.Q),
    Q_model_m3s = qhat,
    P_measured_MW = Float64.(df.P_elec_MW),
    P_model_MW = phat,
    rolling_RMSE_Q_m3s = rmse_q,
    rolling_RMSE_P_MW = rmse_p,
    rolling_REL_RMSE_Q_pct = rel_q,
    rolling_REL_RMSE_P_pct = rel_p,
    model_health = health,
    static_upward_headroom_MW = static_headroom_MW,
    screening_fcr_MW = screening_fcr_MW,
    trusted_screening_fcr_MW = trusted_fcr_MW,
)
CSV.write(joinpath(result_dir, "objective_05_online_timeseries.csv"), online)

active_health = filter(!=("OFFLINE"), health)
green_n = count(==("GREEN"), active_health)
amber_n = count(==("AMBER"), active_health)
red_n = count(==("RED"), active_health)
health_n = length(active_health)

green_pct = health_n > 0 ? 100green_n / health_n : 0.0
trusted_vals = filter(isfinite, trusted_fcr_MW)
median_trusted = isempty(trusted_vals) ? NaN : median(trusted_vals)
max_trusted = isempty(trusted_vals) ? NaN : maximum(trusted_vals)

summary = DataFrame(
    quantity = [
        "stream_samples",
        "rolling_window_s",
        "green_samples",
        "amber_samples",
        "red_samples",
        "green_fraction_pct",
        "objective_02_robust_screening_MW",
        "median_trusted_screening_MW",
        "max_trusted_screening_MW",
        "formal_qualification_state",
        "streaming_mode"
    ],
    value = Any[
        n,
        WINDOW_S,
        green_n,
        amber_n,
        red_n,
        green_pct,
        ROBUST_FCR_MW,
        median_trusted,
        max_trusted,
        "NOT ESTABLISHED",
        "OFFLINE REPLAY OF REAL TROLLHEIM MEASUREMENTS"
    ]
)
CSV.write(joinpath(result_dir, "objective_05_summary.csv"), summary)

println("=== FREKI Objective 5: online-style monitor ===")
println("Replay samples: $n")
@printf("Rolling window: %d s\n", WINDOW_S)
@printf("GREEN fraction of evaluated samples: %.1f %%\n", green_pct)
@printf("Objective 2 robust screening reference: %.2f MW\n", ROBUST_FCR_MW)
if !isempty(trusted_vals)
    @printf("Median trusted screening signal: %.2f MW\n", median_trusted)
    @printf("Maximum trusted screening signal: %.2f MW\n", max_trusted)
else
    println("No GREEN interval produced a trusted screening signal")
end
println("Formal qualified capacity: NOT ESTABLISHED")
println("Mode: offline replay, not live SCADA")

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
default(size=(950, 500), linewidth=1.5, legend=:best)
t = Float64.(df.t)

p1 = plot(t, Float64.(df.P_elec_MW), label="Measured power", ylabel="MW", xlabel="Time [s]")
plot!(p1, t, phat, label="Online local-model prediction")
title!(p1, "Objective 5 — streaming replay: measured vs predicted power")
savefig(p1, joinpath(plot_dir, "05_01_online_power_prediction.png"))

p2 = plot(t, rel_q, label="Flow relative RMSE", ylabel="Rolling RMSE [%]", xlabel="Time [s]")
plot!(p2, t, rel_p, label="Power relative RMSE")
hline!(p2, [1.0], linestyle=:dash, label="GREEN threshold")
hline!(p2, [3.0], linestyle=:dot, label="RED threshold")
title!(p2, "Rolling 60 s model residuals")
savefig(p2, joinpath(plot_dir, "05_02_rolling_residuals.png"))

health_score = [h == "GREEN" ? 3.0 : h == "AMBER" ? 2.0 : h == "RED" ? 1.0 : 0.0 for h in health]
p3 = plot(t, health_score, label="Model health", yticks=([0,1,2,3],["OFF","RED","AMBER","GREEN"]),
          ylim=(-0.2,3.2), xlabel="Time [s]", ylabel="Health")
title!(p3, "Online model-health state")
savefig(p3, joinpath(plot_dir, "05_03_model_health.png"))

p4 = plot(t, screening_fcr_MW, label="Available screening estimate", ylabel="FCR [MW]", xlabel="Time [s]")
plot!(p4, t, trusted_fcr_MW, label="Trusted when model health is GREEN")
hline!(p4, [ROBUST_FCR_MW], linestyle=:dash, label="Objective 2 robust reference")
title!(p4, "Online model-based FCR screening signal")
savefig(p4, joinpath(plot_dir, "05_04_online_fcr_capacity.png"))

p5 = bar(["GREEN","AMBER","RED"], [green_n, amber_n, red_n], legend=false, ylabel="Evaluated samples")
title!(p5, "Objective 5 model-health occupancy")
savefig(p5, joinpath(plot_dir, "05_05_health_occupancy.png"))

println("Generated Objective 5 CSVs and plots.")
