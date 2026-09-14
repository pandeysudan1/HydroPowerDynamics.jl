# FREKI Objective 3 — identify practical measures that increase FCR screening capacity
#
# Starts from the executed Objective 2 robust screening result and asks:
#   What engineering change increases credible reserve capability most?
#
# This is a decision-support study, not a claim of formal Statnett qualification.

using CSV
using DataFrames
using Printf
using Plots

const RHO = 1000.0
const G = 9.81
const ETA_GEN = 0.985
const D_RUNNER = 2.5
const H0 = 384.9
const GATE0 = 0.928
const GATE_MIN = 0.05
const GATE_MAX = 1.00
const DT = 0.02
const T_END = 30.0
const T_CHECK = 5.0
const DELIVERY_MIN = 0.95
const DELIVERY_MAX = 1.20

base = @__DIR__
plot_dir = joinpath(base, "plots")
result_dir = joinpath(base, "results")
mkpath(plot_dir); mkpath(result_dir)

obj1b = CSV.read(joinpath(result_dir, "objective_01b_summary.csv"), DataFrame)
obj2  = CSV.read(joinpath(result_dir, "objective_02_summary.csv"), DataFrame)

function getnum(df, name)
    i = findfirst(==(name), df.quantity)
    i === nothing && error("Missing $name")
    parse(Float64, string(df.value[i]))
end

KQ = getnum(obj1b, "K_q")
ETA_MAX = getnum(obj1b, "eta_max")
QR = getnum(obj1b, "Q_rated")
BASELINE = getnum(obj2, "robust_screening_capacity_MW")

struct Case
    name::String
    description::String
    gate0::Float64
    gate_max::Float64
    dymax::Float64
    Tg::Float64
    H::Float64
    Kq::Float64
    eta_max::Float64
    c_eta::Float64
end

# The baseline mirrors the conservative Objective 2 uncertainty case.
baseline = Case(
    "baseline",
    "Current conservative screening case",
    GATE0, GATE_MAX, 0.12, 0.50, 0.98H0, 0.985KQ, 0.99ETA_MAX, 0.50)

# One-at-a-time engineering interventions. Values are study assumptions, not plant promises.
cases = Case[
    baseline,
    Case("faster_servo", "Reduce effective governor/servo time constant 0.50 -> 0.30 s",
         GATE0, GATE_MAX, 0.12, 0.30, 0.98H0, 0.985KQ, 0.99ETA_MAX, 0.50),
    Case("higher_gate_rate", "Increase permitted guide-vane rate 0.12 -> 0.18 pu/s",
         GATE0, GATE_MAX, 0.18, 0.50, 0.98H0, 0.985KQ, 0.99ETA_MAX, 0.50),
    Case("lower_operating_gate", "Pre-position unit at gate 0.88 pu before upward reserve service",
         0.88, GATE_MAX, 0.12, 0.50, 0.98H0, 0.985KQ, 0.99ETA_MAX, 0.50),
    Case("gate_travel_extension", "Increase usable upper gate from 1.00 -> 1.03 pu equivalent",
         GATE0, 1.03, 0.12, 0.50, 0.98H0, 0.985KQ, 0.99ETA_MAX, 0.50),
    Case("nominal_head", "Operate at nominal measured head instead of 2% low-head case",
         GATE0, GATE_MAX, 0.12, 0.50, H0, 0.985KQ, 0.99ETA_MAX, 0.50),
    Case("better_turbine_gain", "Remove conservative -1.5% turbine-gain derating",
         GATE0, GATE_MAX, 0.12, 0.50, 0.98H0, KQ, 0.99ETA_MAX, 0.50),
]

function turbine_power_MW(y, c::Case)
    q = y * c.Kq * D_RUNNER^2 * sqrt(c.H)
    eta = c.eta_max * (1 - c.c_eta * (q / QR - 1)^2)
    eta = clamp(eta, 0.0, 1.0)
    ETA_GEN * RHO * G * q * c.H * eta / 1e6
end

function deltaP_MW(y, c::Case)
    turbine_power_MW(y, c) - turbine_power_MW(c.gate0, c)
end

function gate_for_target(target, c::Case)
    target <= 0 && return c.gate0
    if target >= deltaP_MW(c.gate_max, c)
        return c.gate_max
    end
    lo, hi = c.gate0, c.gate_max
    for _ in 1:60
        mid = (lo + hi)/2
        if deltaP_MW(mid, c) < target
            lo = mid
        else
            hi = mid
        end
    end
    (lo + hi)/2
end

function simulate(cap, c::Case)
    t = collect(0.0:DT:T_END)
    y = similar(t); dp = similar(t)
    y[1] = c.gate0; dp[1] = 0.0
    ycmd = gate_for_target(cap, c)
    hit_rate = false
    hit_gate = false
    for k in 1:length(t)-1
        raw = (ycmd - y[k]) / c.Tg
        if abs(raw) > c.dymax + 1e-10
            hit_rate = true
        end
        rate = clamp(raw, -c.dymax, c.dymax)
        y[k+1] = clamp(y[k] + DT*rate, GATE_MIN, c.gate_max)
        if y[k+1] >= c.gate_max - 1e-6
            hit_gate = true
        end
        dp[k+1] = deltaP_MW(y[k+1], c)
    end
    i5 = argmin(abs.(t .- T_CHECK))
    d5 = dp[i5]; dend = dp[end]
    pass = d5 >= DELIVERY_MIN*cap && dend >= DELIVERY_MIN*cap && dend <= DELIVERY_MAX*cap
    (;t,y,dp,d5,dend,ycmd,pass,hit_rate,hit_gate)
end

function max_capacity(c::Case)
    best = 0.0
    besttr = simulate(0.0, c)
    for cap in 0.05:0.05:25.0
        tr = simulate(cap, c)
        if tr.pass
            best = cap
            besttr = tr
        else
            break
        end
    end
    best, besttr
end

rows = DataFrame(case=String[], description=String[], capacity_MW=Float64[], gain_MW=Float64[], gain_pct=Float64[], static_headroom_MW=Float64[], active_constraint=String[])
traces = Dict{String,Any}()

for c in cases
    cap, tr = max_capacity(c)
    gain = cap - BASELINE
    pct = BASELINE > 0 ? 100gain/BASELINE : 0.0
    static = deltaP_MW(c.gate_max, c)
    constraint = tr.hit_gate ? "gate headroom" : tr.hit_rate ? "gate-rate / servo dynamics" : "delivery envelope / static gain"
    push!(rows, (c.name, c.description, cap, gain, pct, static, constraint))
    traces[c.name] = tr
end

# Rank only intervention cases, not the baseline itself.
interventions = rows[rows.case .!= "baseline", :]
sort!(interventions, :gain_MW, rev=true)
best_case = interventions[1, :]

println("=== FREKI Objective 3: increase FCR ===")
@printf("Objective 2 robust baseline = %.2f MW\n", BASELINE)
@printf("Best single intervention = %s\n", best_case.case)
@printf("Screening capacity after intervention = %.2f MW\n", best_case.capacity_MW)
@printf("Increase = %.2f MW (%.1f %%)\n", best_case.gain_MW, best_case.gain_pct)
println("Limiting mechanism after intervention = $(best_case.active_constraint)")

CSV.write(joinpath(result_dir, "objective_03_interventions.csv"), rows)
CSV.write(joinpath(result_dir, "objective_03_ranking.csv"), interventions)

summary = DataFrame(
    quantity=["baseline_robust_MW","best_intervention","best_screening_capacity_MW","increase_MW","increase_pct","post_intervention_constraint","formal_qualification_state"],
    value=Any[BASELINE,best_case.case,best_case.capacity_MW,best_case.gain_MW,best_case.gain_pct,best_case.active_constraint,"NOT ESTABLISHED"])
CSV.write(joinpath(result_dir, "objective_03_summary.csv"), summary)

# Plot 1 — intervention ranking
p1 = bar(interventions.case, interventions.capacity_MW, legend=false,
         ylabel="Screening capacity [MW]", xlabel="Engineering intervention",
         xrotation=25, title="FREKI Objective 3 — which change increases FCR most?")
hline!(p1, [BASELINE], linestyle=:dash, label="Objective 2 robust baseline")
savefig(p1, joinpath(plot_dir, "03_01_intervention_ranking.png"))

# Plot 2 — incremental gain over baseline
p2 = bar(interventions.case, interventions.gain_MW, legend=false,
         ylabel="Capacity increase [MW]", xlabel="Engineering intervention",
         xrotation=25, title="Incremental FCR screening gain over baseline")
savefig(p2, joinpath(plot_dir, "03_02_capacity_gain.png"))

# Plot 3 — static headroom versus dynamic screening capacity
p3 = scatter(rows.static_headroom_MW, rows.capacity_MW,
             xlabel="Static upward headroom [MW]", ylabel="Dynamic screening capacity [MW]",
             title="Static headroom versus dynamically deliverable reserve", legend=false)
for r in eachrow(rows)
    annotate!(p3, r.static_headroom_MW, r.capacity_MW, text(r.case, 7))
end
savefig(p3, joinpath(plot_dir, "03_03_headroom_vs_dynamic.png"))

# Plot 4 — baseline versus best intervention response
bestname = String(best_case.case)
tr0 = traces["baseline"]
trb = traces[bestname]
p4 = plot(tr0.t, tr0.dp, label=@sprintf("baseline %.2f MW", BASELINE), xlabel="Time [s]", ylabel="ΔP [MW]")
plot!(p4, trb.t, trb.dp, label=@sprintf("%s %.2f MW", bestname, best_case.capacity_MW))
title!(p4, "Baseline and best-intervention reserve delivery")
savefig(p4, joinpath(plot_dir, "03_04_baseline_vs_best_response.png"))

# Plot 5 — decision table encoded as a bar plot of gain percentages
p5 = bar(interventions.case, interventions.gain_pct, legend=false,
         ylabel="Capacity increase [%]", xlabel="Intervention", xrotation=25,
         title="Engineering leverage for FCR screening capacity")
savefig(p5, joinpath(plot_dir, "03_05_engineering_leverage.png"))

println("Generated Objective 3 CSVs and plots.")
