# FREKI Objective 2 — robust FCR-capacity screening from a measurement-supported Trollheim model
#
# This script deliberately separates
#   (1) model-based engineering screening capacity, and
#   (2) formally qualified FCR capacity.
#
# Objective 1B supports the local turbine map near full load, but hydraulic and
# governor dynamics are still AMBER. Therefore this script does NOT claim formal
# prequalification. It computes a conservative screening capacity across a small
# uncertainty set and reports what limits the candidate.

using CSV
using DataFrames
using Statistics
using Printf
using Plots

const RHO = 1000.0
const G = 9.81
const ETA_GEN = 0.985
const D_RUNNER = 2.5
const F0 = 50.0
const DF_FULL = 0.10

# Measured operating point from Trollheim Objective 1B / experimental analysis
const H0 = 384.9
const Q0 = 36.36
const GATE0 = 0.928
const P0 = 128.26

# Plant-side screening limits. These are engineering assumptions for this study,
# not Statnett requirements.
const GATE_MIN = 0.05
const GATE_MAX = 1.00
const DY_MAX = 0.12           # pu/s
const T_END = 30.0
const DT = 0.02
const RESPONSE_CHECK_T = 5.0
const DELIVERY_MIN = 0.95     # screening acceptance band
const DELIVERY_MAX = 1.20

base = @__DIR__
plot_dir = joinpath(base, "plots")
result_dir = joinpath(base, "results")
mkpath(plot_dir)
mkpath(result_dir)

summary_path = joinpath(result_dir, "objective_01b_summary.csv")
obj1b = CSV.read(summary_path, DataFrame)

function getnum(name)
    row = findfirst(==(name), obj1b.quantity)
    row === nothing && error("Missing $name in Objective 1B summary")
    parse(Float64, string(obj1b.value[row]))
end

function getstr(name)
    row = findfirst(==(name), obj1b.quantity)
    row === nothing && error("Missing $name in Objective 1B summary")
    string(obj1b.value[row])
end

KQ_ID = getnum("K_q")
ETA_ID = getnum("eta_max")
QR_ID = getnum("Q_rated")
CETA_ID = getnum("c_eta")
TURBINE_HEALTH = getstr("turbine_model_health")
HYD_HEALTH = getstr("hydraulic_model_health")
GOV_HEALTH = getstr("governor_model_health")

# Negative c_eta from Objective 1B is not physically credible. Objective 2 therefore
# does not use that value outside the narrow calibration region. Instead, it brackets
# local efficiency curvature with physically non-negative values.

struct Scenario
    name::String
    Kq::Float64
    eta_max::Float64
    c_eta::Float64
    H::Float64
    Tg::Float64
end

scenarios = Scenario[
    Scenario("nominal",      KQ_ID,        ETA_ID,        0.25, H0,       0.30),
    Scenario("low_gain",     0.99KQ_ID,    0.995ETA_ID,   0.40, 0.99H0,   0.40),
    Scenario("slow_servo",   KQ_ID,        ETA_ID,        0.25, H0,       0.50),
    Scenario("low_head",     KQ_ID,        ETA_ID,        0.25, 0.98H0,   0.35),
    Scenario("conservative", 0.985KQ_ID,   0.99ETA_ID,    0.50, 0.98H0,   0.50),
]

function turbine_power_MW(y, s::Scenario)
    q = y * s.Kq * D_RUNNER^2 * sqrt(s.H)
    eta = s.eta_max * (1 - s.c_eta * (q / QR_ID - 1)^2)
    eta = clamp(eta, 0.0, 1.0)
    return ETA_GEN * RHO * G * q * s.H * eta / 1e6
end

# Re-anchor each scenario to the measured operating point. Capacity is based on
# incremental power, so this avoids treating small static parameter mismatch as
# an artificial reserve contribution.
function deltaP_MW(y, s::Scenario)
    turbine_power_MW(y, s) - turbine_power_MW(GATE0, s)
end

function gate_for_target(target_MW, s::Scenario)
    target_MW <= 0 && return GATE0
    maxinc = deltaP_MW(GATE_MAX, s)
    target_MW >= maxinc && return GATE_MAX
    lo, hi = GATE0, GATE_MAX
    for _ in 1:60
        mid = 0.5(lo + hi)
        if deltaP_MW(mid, s) < target_MW
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5(lo + hi)
end

function simulate_capacity(cap_MW, s::Scenario)
    t = collect(0.0:DT:T_END)
    y = similar(t)
    dp = similar(t)
    y[1] = GATE0
    dp[1] = 0.0
    ycmd = gate_for_target(cap_MW, s)

    for k in 1:length(t)-1
        # first-order servo plus explicit guide-vane rate limit
        raw_rate = (ycmd - y[k]) / s.Tg
        rate = clamp(raw_rate, -DY_MAX, DY_MAX)
        y[k+1] = clamp(y[k] + DT * rate, GATE_MIN, GATE_MAX)
        dp[k+1] = deltaP_MW(y[k+1], s)
    end

    i5 = argmin(abs.(t .- RESPONSE_CHECK_T))
    delivered5 = dp[i5]
    delivered_end = dp[end]
    lower = DELIVERY_MIN * cap_MW
    upper = DELIVERY_MAX * cap_MW

    pass_delivery = delivered5 >= lower && delivered_end >= lower && delivered_end <= upper
    pass_gate = maximum(y) <= GATE_MAX + 1e-9
    pass = pass_delivery && pass_gate

    return (; t, y, dp, delivered5, delivered_end, ycmd, pass)
end

# Search in 0.05 MW increments. This is intentionally transparent and reproducible.
capacity_grid = collect(0.05:0.05:20.0)
scenario_rows = DataFrame(
    scenario=String[], max_screen_MW=Float64[], static_headroom_MW=Float64[],
    delivered_5s_MW=Float64[], delivered_30s_MW=Float64[], gate_cmd_pu=Float64[])

best_by_scenario = Dict{String,Float64}()
trace_by_scenario = Dict{String,Any}()

for s in scenarios
    best = 0.0
    besttrace = simulate_capacity(0.0, s)
    for c in capacity_grid
        tr = simulate_capacity(c, s)
        if tr.pass
            best = c
            besttrace = tr
        else
            break
        end
    end
    best_by_scenario[s.name] = best
    trace_by_scenario[s.name] = besttrace
    push!(scenario_rows, (
        s.name,
        best,
        deltaP_MW(GATE_MAX, s),
        besttrace.delivered5,
        besttrace.delivered_end,
        besttrace.ycmd,
    ))
end

robust_capacity = minimum(values(best_by_scenario))
limiting_scenario = first([k for (k,v) in best_by_scenario if v == robust_capacity])
nominal_capacity = best_by_scenario["nominal"]

# Model-health gate: a numerical screening result can be produced, but formal
# qualified capacity is intentionally withheld until hydraulic/governor evidence is GREEN.
all_required_green = (HYD_HEALTH == "GREEN") && (GOV_HEALTH == "GREEN")
qualification_state = all_required_green ? @sprintf("%.2f MW candidate for formal test", robust_capacity) : "NOT ESTABLISHED"

println("=== FREKI Objective 2: FCR capacity screening ===")
@printf("Measured operating point: P0 = %.2f MW, gate = %.3f pu, H = %.1f m\n", P0, GATE0, H0)
@printf("Nominal screening capacity = %.2f MW\n", nominal_capacity)
@printf("Robust screening capacity  = %.2f MW\n", robust_capacity)
println("Limiting uncertainty scenario = $limiting_scenario")
println("Formal qualified FCR capacity = $qualification_state")
println("Hydraulic model health = $HYD_HEALTH")
println("Governor model health  = $GOV_HEALTH")

CSV.write(joinpath(result_dir, "objective_02_scenarios.csv"), scenario_rows)

summary = DataFrame(
    quantity=[
        "nominal_screening_capacity_MW",
        "robust_screening_capacity_MW",
        "limiting_scenario",
        "formal_qualified_capacity",
        "hydraulic_model_health",
        "governor_model_health"
    ],
    value=Any[
        nominal_capacity,
        robust_capacity,
        limiting_scenario,
        qualification_state,
        HYD_HEALTH,
        GOV_HEALTH
    ]
)
CSV.write(joinpath(result_dir, "objective_02_summary.csv"), summary)

# Plot 1: scenario capacity comparison
p1 = bar(scenario_rows.scenario, scenario_rows.max_screen_MW,
         ylabel="Screening capacity [MW]", xlabel="Uncertainty scenario",
         legend=false, title="FREKI Objective 2 — model-based FCR screening capacity",
         xrotation=20)
hline!(p1, [robust_capacity], linestyle=:dash, label="robust minimum")
savefig(p1, joinpath(plot_dir, "02_01_capacity_by_scenario.png"))

# Plot 2: static headroom by scenario
p2 = bar(scenario_rows.scenario, scenario_rows.static_headroom_MW,
         ylabel="Static upward headroom [MW]", xlabel="Uncertainty scenario",
         legend=false, title="Upward power headroom from measured gate position to full gate",
         xrotation=20)
savefig(p2, joinpath(plot_dir, "02_02_static_headroom.png"))

# Plot 3: robust-scenario response
srob = scenarios[findfirst(s -> s.name == limiting_scenario, scenarios)]
trrob = simulate_capacity(robust_capacity, srob)
p3 = plot(trrob.t, trrob.dp, label="Delivered incremental power", ylabel="ΔP [MW]", xlabel="Time [s]")
hline!(p3, [robust_capacity], linestyle=:dash, label="screening offer")
hline!(p3, [DELIVERY_MIN*robust_capacity], linestyle=:dot, label="95% lower screening bound")
hline!(p3, [DELIVERY_MAX*robust_capacity], linestyle=:dot, label="120% upper screening bound")
title!(p3, "Robust-scenario reserve delivery")
savefig(p3, joinpath(plot_dir, "02_03_robust_response.png"))

# Plot 4: gate trajectory
p4 = plot(trrob.t, trrob.y, label="Guide vane", ylabel="Gate [pu]", xlabel="Time [s]")
hline!(p4, [GATE_MAX], linestyle=:dash, label="full gate")
hline!(p4, [GATE0], linestyle=:dot, label="measured operating gate")
title!(p4, "Guide-vane use at robust screening capacity")
savefig(p4, joinpath(plot_dir, "02_04_gate_response.png"))

# Plot 5: decision bridge from model health to capacity claim
labels = ["Local turbine map", "Hydraulic dynamics", "Governor dynamics", "Formal FCR claim"]
score = [3.0, HYD_HEALTH == "GREEN" ? 3.0 : 2.0, GOV_HEALTH == "GREEN" ? 3.0 : 2.0, all_required_green ? 3.0 : 1.0]
p5 = bar(labels, score, ylim=(0,3.5), legend=false, ylabel="Evidence level", xrotation=20,
         title="Model-health gate before claiming qualified FCR")
savefig(p5, joinpath(plot_dir, "02_05_model_health_gate.png"))

println("Generated Objective 2 CSVs and plots.")
