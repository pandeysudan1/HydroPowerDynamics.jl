# FREKI Objective 1: validate an HPD plant model from normal-operation data.
#
# Workflow:
# small normal load variations -> nonlinear HPD truth -> noisy measurements
# -> identify R, Tg and effective penstock f_D -> infer roughness
# -> second nonlinear HPD replay -> held-out validation.
#
# The identification is equation-based. It avoids a large black-box optimizer.

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
const ROUGHNESS_TRUE = 1.5e-5
const ETA0 = ETA_MAX * (1 - C_ETA * (Q0 / Q_RATED - 1)^2)
const T_END = 65.0
const T_CAL_END = 35.0
const DT = 0.05

# Small ambient multi-sine load motion. The two small faster terms add enough
# bandwidth to see the servo without creating a dedicated qualification event.
function normal_load(tt)
    tt < 0 && return P0
    P0 +
    0.55e6 * sin(2π * tt / 17.0) +
    0.30e6 * sin(2π * tt / 7.5) +
    0.15e6 * sin(2π * tt / 31.0) +
    0.10e6 * sin(2π * tt / 2.4) +
    0.06e6 * sin(2π * tt / 1.2)
end

function pipe_guess(L, Dpipe, roughness)
    A = π * Dpipe^2 / 4
    μ = 1e-3
    Re = DM0 * Dpipe / (μ * A)
    fD = 0.25 / log10(roughness/(3.7Dpipe) + 5.74/Re^0.9)^2
    dpf = fD * (L/Dpipe) * (DM0/(RHO*A))^2 * RHO/2
    (; Re, fD, dpf)
end

function build_model(; R=R_TRUE, Tg=TG_TRUE, roughness=ROUGHNESS_TRUE)
    @named upper = Reservoir(H=H_GROSS, rho=RHO, g=G, p_atm=P_ATM)
    @named tail = Reservoir(H=0.0, rho=RHO, g=G, p_atm=P_ATM)
    @named headrace = Penstock(L=L_HEADRACE, D_pipe=D_HEADRACE, rho=RHO, B=2e9, e=roughness)
    @named penstock = Penstock(L=L_PENSTOCK, D_pipe=D_PENSTOCK, rho=RHO, B=2e9, e=roughness)
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
            P0 + 0.55e6*sin(2π*t/17.0) + 0.30e6*sin(2π*t/7.5) +
            0.15e6*sin(2π*t/31.0) + 0.10e6*sin(2π*t/2.4) + 0.06e6*sin(2π*t/1.2)),
    ]
    raw = ODESystem(eqs, t; name=:FREKIObjective1,
        systems=[upper, tail, headrace, penstock, turbine, rotor, generator, speed, governor])
    sys = structural_simplify(raw)
    (; sys, headrace, penstock, turbine, rotor, generator, governor)
end

function simulate(; R=R_TRUE, Tg=TG_TRUE, roughness=ROUGHNESS_TRUE)
    m = build_model(R=R, Tg=Tg, roughness=roughness)
    hrg = pipe_guess(L_HEADRACE, D_HEADRACE, roughness)
    png = pipe_guess(L_PENSTOCK, D_PENSTOCK, roughness)
    u0 = Dict(
        m.headrace.dm => DM0, m.headrace.p_avg => P_UP,
        m.penstock.dm => DM0, m.penstock.p_avg => P_UP,
        m.rotor.omega => W0, m.rotor.theta => 0.0,
        m.governor.xi => 0.0, m.governor.gate => GATE0,
    )
    guesses = Dict(
        m.headrace.Re => hrg.Re, m.headrace.f_D => hrg.fD, m.headrace.dp_f => hrg.dpf,
        m.penstock.Re => png.Re, m.penstock.f_D => png.fD, m.penstock.dp_f => png.dpf,
        m.turbine.H => H_GROSS, m.turbine.Q => Q0, m.turbine.eta => ETA0,
        m.turbine.P_mech => P0, m.turbine.dm => DM0, m.turbine.tau_shaft => P0/W0,
    )
    prob = ODEProblem(m.sys, u0, (-300.0, T_END); guesses)
    sol = solve(prob, Rodas5P(); abstol=1e-8, reltol=1e-8, saveat=DT)
    sol.t[end] >= T_END - 1e-6 || error("Simulation stopped early: $(sol.retcode)")
    keep = findall(>=(0.0), sol.t)
    (; t=sol.t[keep],
       f=sol[m.rotor.omega][keep] ./ W0 .* F0,
       pm=sol[m.turbine.P_mech][keep] ./ 1e6,
       pe=sol[m.generator.P_elec][keep] ./ 1e6,
       gate=sol[m.governor.tau_o][keep],
       q=sol[m.turbine.Q][keep],
       p_pen_in=sol[m.penstock.port_a.p][keep],
       p_pen_out=sol[m.penstock.port_b.p][keep],
       fD=sol[m.penstock.f_D][keep],
       Re=sol[m.penstock.Re][keep])
end

function moving_average(x, halfwindow)
    y = similar(x)
    for i in eachindex(x)
        lo=max(firstindex(x),i-halfwindow); hi=min(lastindex(x),i+halfwindow)
        y[i]=mean(@view x[lo:hi])
    end
    y
end

function cumulative_trapezoid(tvec, x)
    z=zeros(length(x))
    for i in 2:length(x)
        z[i]=z[i-1]+0.5*(x[i]+x[i-1])*(tvec[i]-tvec[i-1])
    end
    z
end

# Step 1: estimate R. The derivative regression is retained only for the ratio
# a/b = R because the previous experiment showed that this ratio is robust even
# when the individual fast time constant estimate is not.
function estimate_R(tvec, fmeas, gatemeas)
    gs=moving_average(gatemeas,6)
    dg=zeros(length(gs))
    for i in 2:length(gs)-1
        dg[i]=(gs[i+1]-gs[i-1])/(tvec[i+1]-tvec[i-1])
    end
    dg[1]=dg[2]; dg[end]=dg[end-1]
    ef=(F0 .- fmeas)./F0
    xi=cumulative_trapezoid(tvec,ef)
    x1=GATE0 .+ KI.*xi .- gs
    idx=findall(i->2.0<=tvec[i]<=T_CAL_END,eachindex(tvec))
    a,b=hcat(x1[idx],ef[idx]) \ dg[idx]
    Rhat=a/b
    (;Rhat,a,b,ef,xi,gs)
end

# Step 2: with R fixed, estimate Tg from the integrated servo equation over
# many short windows. This avoids differentiating the gate to estimate Tg.
function estimate_Tg(tvec, fmeas, gatemeas, Rhat)
    gs=moving_average(gatemeas,4)
    fs=moving_average(fmeas,2)
    ef=(F0 .- fs)./F0
    xi=cumulative_trapezoid(tvec,ef)
    ucmd=GATE0 .+ ef./Rhat .+ KI.*xi
    drive=ucmd .- gs
    Idrive=cumulative_trapezoid(tvec,drive)
    w=max(2,round(Int,0.8/DT))
    xv=Float64[]; yv=Float64[]
    for i in eachindex(tvec)
        j=i+w
        j>length(tvec) && break
        (2.0<=tvec[i] && tvec[j]<=T_CAL_END) || continue
        push!(xv,Idrive[j]-Idrive[i])
        push!(yv,gs[j]-gs[i])
    end
    c=dot(xv,yv)/dot(xv,xv)
    (;Tghat=1/c, c, ucmd, gs)
end

# Step 3: identify effective Darcy friction directly from the HPD momentum law
#
# (L/A) d(dm)/dt = dp - fD*k*dm^2,
# k=L/(2*D*rho*A^2).
#
# Integrating over short windows gives a scalar linear regression for fD.
function estimate_friction(tvec, qmeas, pin, pout)
    A=π*D_PENSTOCK^2/4
    k=L_PENSTOCK/(2D_PENSTOCK*RHO*A^2)
    dm=RHO .* moving_average(qmeas,4)
    dp=moving_average(pin .- pout,3)
    Idp=cumulative_trapezoid(tvec,dp)
    Idm2=cumulative_trapezoid(tvec,dm.^2)
    w=max(2,round(Int,1.0/DT))
    xv=Float64[]; yv=Float64[]
    for i in eachindex(tvec)
        j=i+w
        j>length(tvec) && break
        (2.0<=tvec[i] && tvec[j]<=T_CAL_END) || continue
        push!(xv,k*(Idm2[j]-Idm2[i]))
        push!(yv,(Idp[j]-Idp[i])-(L_PENSTOCK/A)*(dm[j]-dm[i]))
    end
    fDhat=dot(xv,yv)/dot(xv,xv)
    Rehat=mean(abs.(dm).*D_PENSTOCK./(1e-3*A))
    term=10.0^(-0.5/sqrt(fDhat))-5.74/Rehat^0.9
    ehat=3.7D_PENSTOCK*term
    (;fDhat,Rehat,ehat)
end

rmse(a,b)=sqrt(mean((a.-b).^2))
function fit_percent(y,yhat)
    den=norm(y.-mean(y)); den<eps() && return NaN
    100*(1-norm(y.-yhat)/den)
end

println("=== FREKI Objective 1: controller + hydraulic model validation ===")
println("Solving nonlinear HPD truth model...")
truth=simulate()

Random.seed!(240914)
σf=0.0020; σp=0.020; σg=0.00020; σq=0.0030; σpressure=20.0
f_meas=truth.f .+ σf.*randn(length(truth.f))
pe_meas=truth.pe .+ σp.*randn(length(truth.pe))
gate_meas=truth.gate .+ σg.*randn(length(truth.gate))
q_meas=truth.q .+ σq.*randn(length(truth.q))
pin_meas=truth.p_pen_in .+ σpressure.*randn(length(truth.p_pen_in))
pout_meas=truth.p_pen_out .+ σpressure.*randn(length(truth.p_pen_out))

r_est=estimate_R(truth.t,f_meas,gate_meas)
tg_est=estimate_Tg(truth.t,f_meas,gate_meas,r_est.Rhat)
hyd_est=estimate_friction(truth.t,q_meas,pin_meas,pout_meas)
fD_ref=mean(truth.fD[findall(<=(T_CAL_END),truth.t)])

@printf("R: truth %.6f, identified %.6f, error %.2f%%\n",R_TRUE,r_est.Rhat,100*(r_est.Rhat/R_TRUE-1))
@printf("Tg: truth %.6f s, identified %.6f s, error %.2f%%\n",TG_TRUE,tg_est.Tghat,100*(tg_est.Tghat/TG_TRUE-1))
@printf("penstock fD: reference %.8f, identified %.8f, error %.2f%%\n",fD_ref,hyd_est.fDhat,100*(hyd_est.fDhat/fD_ref-1))
@printf("roughness: truth %.3e m, inferred %.3e m\n",ROUGHNESS_TRUE,hyd_est.ehat)

(0.2<r_est.Rhat<1.0) || error("R estimate outside physical range")
(0.05<tg_est.Tghat<1.5) || error("Tg estimate outside physical range")
(0.002<hyd_est.fDhat<0.05) || error("fD estimate outside physical range")
(1e-7<hyd_est.ehat<5e-3) || error("roughness estimate outside physical range")

println("Replaying nonlinear HPD with identified R, Tg and roughness...")
identified=simulate(R=r_est.Rhat,Tg=tg_est.Tghat,roughness=hyd_est.ehat)
ival=findall(>(T_CAL_END),truth.t)
metrics=Dict(
    "RMSE_frequency_Hz"=>rmse(f_meas[ival],identified.f[ival]),
    "RMSE_electrical_power_MW"=>rmse(pe_meas[ival],identified.pe[ival]),
    "RMSE_gate_pu"=>rmse(gate_meas[ival],identified.gate[ival]),
    "RMSE_flow_m3s"=>rmse(q_meas[ival],identified.q[ival]),
    "FIT_frequency_pct"=>fit_percent(f_meas[ival],identified.f[ival]),
    "FIT_electrical_power_pct"=>fit_percent(pe_meas[ival],identified.pe[ival]),
    "FIT_flow_pct"=>fit_percent(q_meas[ival],identified.q[ival]),
)
@printf("Held-out RMSE: f %.6f Hz, Pe %.6f MW, gate %.8f pu, Q %.6f m3/s\n",
    metrics["RMSE_frequency_Hz"],metrics["RMSE_electrical_power_MW"],metrics["RMSE_gate_pu"],metrics["RMSE_flow_m3s"])
@printf("Held-out FIT: f %.2f%%, Pe %.2f%%, Q %.2f%%\n",
    metrics["FIT_frequency_pct"],metrics["FIT_electrical_power_pct"],metrics["FIT_flow_pct"])

outdir=joinpath(@__DIR__,"results"); plotdir=joinpath(@__DIR__,"plots")
mkpath(outdir); mkpath(plotdir)
load_MW=normal_load.(truth.t)./1e6
res_f=f_meas.-identified.f; res_pe=pe_meas.-identified.pe

CSV.write(joinpath(outdir,"objective_01_timeseries.csv"),DataFrame(
    time_s=truth.t,load_MW=load_MW,
    measured_frequency_Hz=f_meas,model_frequency_Hz=identified.f,
    measured_electrical_power_MW=pe_meas,model_electrical_power_MW=identified.pe,
    measured_gate_pu=gate_meas,model_gate_pu=identified.gate,
    measured_flow_m3s=q_meas,model_flow_m3s=identified.q,
    measured_penstock_inlet_Pa=pin_meas,measured_penstock_outlet_Pa=pout_meas,
    frequency_residual_Hz=res_f,power_residual_MW=res_pe))

summary=DataFrame(
    quantity=["R","Tg","penstock_fD","roughness","RMSE_frequency","RMSE_electrical_power","RMSE_gate","RMSE_flow","FIT_frequency","FIT_electrical_power","FIT_flow"],
    reference=[R_TRUE,TG_TRUE,fD_ref,ROUGHNESS_TRUE,NaN,NaN,NaN,NaN,NaN,NaN,NaN],
    identified_or_metric=[r_est.Rhat,tg_est.Tghat,hyd_est.fDhat,hyd_est.ehat,
        metrics["RMSE_frequency_Hz"],metrics["RMSE_electrical_power_MW"],metrics["RMSE_gate_pu"],metrics["RMSE_flow_m3s"],
        metrics["FIT_frequency_pct"],metrics["FIT_electrical_power_pct"],metrics["FIT_flow_pct"]],
    unit=["pu/pu","s","-","m","Hz","MW","pu","m3/s","%","%","%"])
CSV.write(joinpath(outdir,"objective_01_summary.csv"),summary)

p1=plot(truth.t,load_MW,xlabel="Time [s]",ylabel="Load [MW]",title="Normal-operation excitation",legend=false)
savefig(p1,joinpath(plotdir,"01_normal_operation_input.png"))
p2=plot(truth.t,f_meas,label="measurement",xlabel="Time [s]",ylabel="Frequency [Hz]",title="Held-out frequency validation")
plot!(p2,truth.t,identified.f,label="identified HPD"); vline!(p2,[T_CAL_END],label="validation starts")
savefig(p2,joinpath(plotdir,"02_frequency_validation.png"))
p3=plot(truth.t,pe_meas,label="measurement",xlabel="Time [s]",ylabel="Electrical power [MW]",title="Held-out power validation")
plot!(p3,truth.t,identified.pe,label="identified HPD"); vline!(p3,[T_CAL_END],label="validation starts")
savefig(p3,joinpath(plotdir,"03_power_validation.png"))
p4=plot(truth.t,gate_meas,label="measured gate",xlabel="Time [s]",ylabel="Gate [pu]",title="Guide-vane validation")
plot!(p4,truth.t,identified.gate,label="identified HPD"); vline!(p4,[T_CAL_END],label="validation starts")
savefig(p4,joinpath(plotdir,"04_gate_validation.png"))
p5=plot(truth.t,res_f,label="frequency residual [Hz]",xlabel="Time [s]",ylabel="Residual",title="Validation residuals")
plot!(p5,truth.t,res_pe./10,label="power residual / 10 [MW]"); vline!(p5,[T_CAL_END],label="validation starts")
savefig(p5,joinpath(plotdir,"05_validation_residuals.png"))
p6=bar(["R","Tg","fD","roughness"],[r_est.Rhat/R_TRUE,tg_est.Tghat/TG_TRUE,hyd_est.fDhat/fD_ref,hyd_est.ehat/ROUGHNESS_TRUE],
    ylabel="identified / reference",title="Identified controller and hydraulic parameters",legend=false)
hline!(p6,[1.0]); savefig(p6,joinpath(plotdir,"06_parameter_identification.png"))
p7=plot(truth.t,q_meas,label="measured flow",xlabel="Time [s]",ylabel="Flow [m3/s]",title="Hydraulic held-out validation")
plot!(p7,truth.t,identified.q,label="identified HPD"); vline!(p7,[T_CAL_END],label="validation starts")
savefig(p7,joinpath(plotdir,"07_hydraulic_identification.png"))

println("RESULT_DIR=",outdir)
println("PLOT_DIR=",plotdir)
println("FREKI_OBJECTIVE_01=PASS")
