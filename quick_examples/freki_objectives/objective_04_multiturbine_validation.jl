# FREKI Objective 4 — common validation workflow across Francis, Pelton and Kaplan
#
# Trollheim Francis case is anchored to measured-data work from Objective 1B.
# Pelton and Kaplan cases in this first implementation are canonical benchmark
# models only; they are not presented as plant-measurement validation.

using DataFrames
using CSV
using Statistics
using Printf
using Plots

base = @__DIR__
plot_dir = joinpath(base, "plots")
result_dir = joinpath(base, "results")
mkpath(plot_dir)
mkpath(result_dir)

# Common test horizon: a small commanded reserve step is applied at t = 2 s.
const DT = 0.01
const TEND = 20.0
const TSTEP = 2.0
const OFFER = 0.08   # pu reserve request relative to unit rating
const F0 = 50.0

abstract type TurbineBenchmark end

struct FrancisBenchmark <: TurbineBenchmark
    name::String
    Pbase_MW::Float64
    Tg::Float64
    Tw::Float64
    rate::Float64
    gain::Float64
end

struct PeltonBenchmark <: TurbineBenchmark
    name::String
    Pbase_MW::Float64
    Tneedle::Float64
    Tjet::Float64
    rate::Float64
    gain::Float64
end

struct KaplanBenchmark <: TurbineBenchmark
    name::String
    Pbase_MW::Float64
    Tgate::Float64
    Tblade::Float64
    rate_gate::Float64
    rate_blade::Float64
    gain::Float64
end

# Francis: Trollheim-inspired benchmark, linked conceptually to Objective 1B.
francis = FrancisBenchmark("Francis / Trollheim-inspired", 150.0, 0.35, 0.55, 0.12, 1.0)
# Pelton/Kaplan: canonical benchmark parameters for methodology comparison only.
pelton  = PeltonBenchmark("Pelton benchmark", 120.0, 0.20, 0.10, 0.20, 1.0)
kaplan  = KaplanBenchmark("Kaplan benchmark", 100.0, 0.40, 0.80, 0.10, 0.08, 1.0)

function satrate(dx, lim)
    clamp(dx, -lim, lim)
end

function simulate(m::FrancisBenchmark)
    t = collect(0.0:DT:TEND)
    gate = zeros(length(t)); flow = zeros(length(t)); power = zeros(length(t))
    for k in 1:length(t)-1
        cmd = t[k] >= TSTEP ? OFFER : 0.0
        dgate = satrate((cmd - gate[k])/m.Tg, m.rate)
        gate[k+1] = gate[k] + DT*dgate
        dflow = (gate[k] - flow[k])/m.Tw
        flow[k+1] = flow[k] + DT*dflow
        power[k+1] = m.gain * flow[k+1]
    end
    return (; t, actuator=gate, aux=flow, power)
end

function simulate(m::PeltonBenchmark)
    t = collect(0.0:DT:TEND)
    needle = zeros(length(t)); jet = zeros(length(t)); power = zeros(length(t))
    for k in 1:length(t)-1
        cmd = t[k] >= TSTEP ? OFFER : 0.0
        dneedle = satrate((cmd - needle[k])/m.Tneedle, m.rate)
        needle[k+1] = needle[k] + DT*dneedle
        djet = (needle[k] - jet[k])/m.Tjet
        jet[k+1] = jet[k] + DT*djet
        power[k+1] = m.gain * jet[k+1]
    end
    return (; t, actuator=needle, aux=jet, power)
end

function simulate(m::KaplanBenchmark)
    t = collect(0.0:DT:TEND)
    gate = zeros(length(t)); blade = zeros(length(t)); power = zeros(length(t))
    for k in 1:length(t)-1
        cmd = t[k] >= TSTEP ? OFFER : 0.0
        dgate = satrate((cmd - gate[k])/m.Tgate, m.rate_gate)
        gate[k+1] = gate[k] + DT*dgate
        # Coordinated blade command lags gate demand.
        dblade = satrate((cmd - blade[k])/m.Tblade, m.rate_blade)
        blade[k+1] = blade[k] + DT*dblade
        power[k+1] = m.gain * (0.6*gate[k+1] + 0.4*blade[k+1])
    end
    return (; t, actuator=gate, aux=blade, power)
end

function metrics(m, tr)
    target = OFFER
    idx5 = argmin(abs.(tr.t .- (TSTEP + 5.0)))
    idx10 = argmin(abs.(tr.t .- (TSTEP + 10.0)))
    p5 = tr.power[idx5]
    p10 = tr.power[idx10]
    pend = tr.power[end]
    rise_idx = findfirst(x -> x >= 0.9target, tr.power)
    t90 = isnothing(rise_idx) ? NaN : tr.t[rise_idx] - TSTEP
    overshoot = 100 * max(maximum(tr.power) - target, 0.0) / target
    delivery5 = 100*p5/target
    delivery10 = 100*p10/target
    return (; p5, p10, pend, t90, overshoot, delivery5, delivery10)
end

models = [francis, pelton, kaplan]
traces = Dict{String,Any}()
rows = DataFrame(
    turbine=String[], evidence=String[], Pbase_MW=Float64[],
    t90_s=Float64[], delivery_5s_pct=Float64[], delivery_10s_pct=Float64[],
    overshoot_pct=Float64[], steady_delivery_pct=Float64[])

for m in models
    tr = simulate(m)
    traces[m.name] = tr
    met = metrics(m, tr)
    evidence = m isa FrancisBenchmark ? "measurement-anchored methodology" : "canonical benchmark only"
    push!(rows, (m.name, evidence, m.Pbase_MW, met.t90, met.delivery5,
                 met.delivery10, met.overshoot, 100*met.pend/OFFER))
end

# A common pass screen used only to compare methodology consistency, not to claim qualification.
rows[!, :method_screen] = [
    (r.delivery_5s_pct >= 95 && r.delivery_10s_pct >= 95 && r.overshoot_pct <= 20) ? "PASS" : "REVIEW"
    for r in eachrow(rows)
]

CSV.write(joinpath(result_dir, "objective_04_multiturbine_summary.csv"), rows)

# Long-form trajectories for reproducibility.
ts = DataFrame(t=traces[francis.name].t)
for m in models
    tr = traces[m.name]
    key = replace(lowercase(first(split(m.name, " "))), "/"=>"")
    ts[!, Symbol(key*"_power_pu")] = tr.power
    ts[!, Symbol(key*"_actuator_pu")] = tr.actuator
    ts[!, Symbol(key*"_aux_pu")] = tr.aux
end
CSV.write(joinpath(result_dir, "objective_04_multiturbine_timeseries.csv"), ts)

# Plot 1: power response comparison
p1 = plot(title="FREKI Objective 4 — common reserve-step response", xlabel="Time [s]", ylabel="Incremental power [pu]")
for m in models
    tr = traces[m.name]
    plot!(p1, tr.t, tr.power, label=m.name)
end
hline!(p1, [OFFER], linestyle=:dash, label="requested reserve")
savefig(p1, joinpath(plot_dir, "04_01_power_response_comparison.png"))

# Plot 2: actuator dynamics
p2 = plot(title="Primary actuator response by turbine type", xlabel="Time [s]", ylabel="Actuator [pu]")
for m in models
    tr = traces[m.name]
    plot!(p2, tr.t, tr.actuator, label=m.name)
end
savefig(p2, joinpath(plot_dir, "04_02_actuator_comparison.png"))

# Plot 3: 90% rise time
p3 = bar(rows.turbine, rows.t90_s, legend=false, ylabel="t90 [s]", xrotation=15,
         title="Dynamic response speed across turbine benchmarks")
savefig(p3, joinpath(plot_dir, "04_03_t90_comparison.png"))

# Plot 4: delivery at 5 and 10 seconds
p4 = bar(rows.turbine, hcat(rows.delivery_5s_pct, rows.delivery_10s_pct),
         label=["5 s" "10 s"], ylabel="Delivered reserve [%]", xrotation=15,
         title="Common reserve-delivery metric across turbine types")
hline!(p4, [95.0], linestyle=:dash, label="95% screen")
savefig(p4, joinpath(plot_dir, "04_04_delivery_comparison.png"))

# Plot 5: evidence status — Francis measured anchor vs canonical benchmark cases.
score = [3.0, 1.0, 1.0]
p5 = bar(rows.turbine, score, ylim=(0,3.5), legend=false, ylabel="Evidence level", xrotation=15,
         title="Evidence status: measured anchor vs benchmark models")
savefig(p5, joinpath(plot_dir, "04_05_evidence_status.png"))

println("=== FREKI Objective 4 ===")
for r in eachrow(rows)
    @printf("%-28s t90=%5.2f s, delivery5=%6.2f%%, screen=%s, evidence=%s\n",
            r.turbine, r.t90_s, r.delivery_5s_pct, r.method_screen, r.evidence)
end
println("Generated Objective 4 CSVs and plots.")
