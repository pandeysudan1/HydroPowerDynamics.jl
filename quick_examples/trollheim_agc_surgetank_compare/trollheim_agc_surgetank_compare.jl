# Trollheim AGC comparison using HydroPowerDynamics.jl
# Same 75 -> 90 MW disturbance; only hydraulic topology differs.

using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq
using Printf

const PBASE=150e6, P0=75e6, P1=90e6, TSTEP=5.0, F0=50.0
const POLES=16, RPM0=120F0/POLES, W0=2π*RPM0/60
const H_GROSS=371.0, RHO=1000.0, G=9.81
const Q_RATED=37.0, D_RUNNER=2.5, ETA_MAX=0.97, C_ETA=0.25, TAU_RATED=0.90
const KQ=Q_RATED/(TAU_RATED*D_RUNNER^2*sqrt(H_GROSS))
const TAU0=0.5384797177203611, Q0=22.13749950628151, DM0=RHO*Q0
const H_INERTIA=1.83, JTOTAL=2H_INERTIA*PBASE/W0^2
const R_PU=0.50, R_DROOP=R_PU*W0/PBASE, TP=0.05, TGOV=0.30, KI=0.10/PBASE
const L_HEADRACE=500.0, D_HEADRACE=6.0, L_PENSTOCK=500.0, D_PENSTOCK=4.0
const ROUGHNESS=1.5e-5
const ATANK=π/4*3.4^2, ZSURGE0=69.97, D_RISER=3.4, L_RISER=87.0, E_RISER=1e-2

function build_trollheim(; with_surge::Bool)
    @named upper=Reservoir(H=H_GROSS,rho=RHO,g=G,p_atm=101325.0)
    @named tail=Reservoir(H=0.0,rho=RHO,g=G,p_atm=101325.0)
    @named headrace=Penstock(L=L_HEADRACE,D_pipe=D_HEADRACE,rho=RHO,B=2e9,e=ROUGHNESS)
    @named penstock=Penstock(L=L_PENSTOCK,D_pipe=D_PENSTOCK,rho=RHO,B=2e9,e=ROUGHNESS)
    @named turbine=FrancisTurbineAffinity(D=D_RUNNER,K_q=KQ,eta_max=ETA_MAX,c_eta=C_ETA,
                                         Q_rated=Q_RATED,rho=RHO,g=G)
    @named rotor=RotorInertia(J=JTOTAL,b_v=0.0,omega_0=W0)
    @named generator=SimpleGenerator(P_rated=P0,omega_rated=W0,omega_s=W0,D_d=0.0,eta_gen=1.0)
    @named speed=RotationalSpeedSensor()
    @named pshaft=MechanicalPowerSensor()
    @named governor=GGOV1Governor(R_droop=R_DROOP,T_p=TP,T_gov=TGOV,omega_ref=W0,P_set=P0,
                                  tau_min=0.05,tau_max=1.0,K_i=KI,tau_0=TAU0)
    systems=Any[upper,tail,headrace,penstock,turbine,rotor,generator,speed,pshaft,governor]
    eqs=Equation[
        connect(upper.port,headrace.port_a),
        connect(headrace.port_b,penstock.port_a),
        connect(penstock.port_b,turbine.port_in),
        connect(turbine.port_out,tail.port),
        connect(governor.gate_out,turbine.opening),
        connect(turbine.shaft,pshaft.flange_a),
        connect(pshaft.flange_b,rotor.flange_turbine),
        connect(rotor.flange_generator,generator.flange,speed.flange),
        connect(speed.w,governor.speed_in),
        connect(pshaft.P,governor.power_in),
    ]
    surge=nothing
    if with_surge
        @named surge_comp=SurgeTank(A_t=ATANK,Z_0=ZSURGE0,D_riser=D_RISER,L_riser=L_RISER,
                                    e_riser=E_RISER,rho=RHO,g=G,p_atm=101325.0)
        surge=surge_comp
        push!(systems,surge)
        deleteat!(eqs,2)
        insert!(eqs,2,connect(headrace.port_b,penstock.port_a,surge.port))
    end
    raw=ODESystem(eqs,t;name=with_surge ? :TrollheimWithSurge : :TrollheimNoSurge,systems=systems)
    sys=structural_simplify(raw)
    (;sys,raw,headrace,penstock,turbine,rotor,generator,governor,surge)
end

function simulate_case(;with_surge::Bool)
    m=build_trollheim(with_surge=with_surge)
    u0=Dict(m.headrace.dm=>DM0,m.penstock.dm=>DM0,m.rotor.omega=>W0,
            m.governor.P_meas=>P0,m.governor.x_gov=>0.0,m.governor.x_int=>TAU0)
    with_surge && (u0[m.surge.Z]=ZSURGE0)
    prob=ODEProblem(m.sys,merge(u0,Dict(m.generator.P_rated=>P0)),(-300.0,65.0))
    load_step=DiscreteCallback(
        (u,tt,integrator)->tt==TSTEP,
        integrator->(integrator.ps[m.generator.P_rated]=P1; nothing),
        save_positions=(true,true))
    sol=solve(prob,Rodas5P();callback=load_step,tstops=[TSTEP],abstol=1e-8,reltol=1e-8,saveat=0.01)
    (;m,sol)
end

function extract(run)
    m,sol=run.m,run.sol
    tv=sol.t
    f=sol[m.rotor.omega]./W0.*F0
    pm=sol[m.turbine.P_mech]./1e6
    gate=sol[m.governor.tau_o]
    q=sol[m.turbine.Q]
    pre=findall(tt->0<=tt<TSTEP,tv); post=findall(>=(TSTEP),tv)
    i=post[argmin(f[post])]
    iae=0.0
    for k in 2:length(tv)
        tv[k]>=TSTEP || continue
        iae += 0.5*(abs(f[k]-F0)+abs(f[k-1]-F0))*(tv[k]-tv[k-1])
    end
    settle=NaN
    for j in post
        if all(abs.(f[j:end].-F0).<=0.05); settle=tv[j]; break; end
    end
    base=(t=tv,f=f,pm=pm,gate=gate,q=q,nadir=f[i],t_nadir=tv[i],final_f=f[end],final_pm=pm[end],
          iae=iae,settling=settle,pre_f_range=maximum(f[pre])-minimum(f[pre]),
          pre_pm_range=maximum(pm[pre])-minimum(pm[pre]))
    m.surge===nothing ? base : merge(base,(Z=sol[m.surge.Z],Qs=sol[m.surge.port.dm]./RHO))
end

function write_csv(path,r;surge=false)
    open(path,"w") do io
        println(io,surge ? "time_s,frequency_Hz,mechanical_power_MW,gate_pu,turbine_flow_m3s,surge_level_m,surge_flow_m3s" :
                          "time_s,frequency_Hz,mechanical_power_MW,gate_pu,turbine_flow_m3s")
        for i in eachindex(r.t)
            r.t[i] >= 0 || continue
            if surge
                @printf(io,"%.6f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",r.t[i],r.f[i],r.pm[i],r.gate[i],r.q[i],r.Z[i],r.Qs[i])
            else
                @printf(io,"%.6f,%.10f,%.10f,%.10f,%.10f\n",r.t[i],r.f[i],r.pm[i],r.gate[i],r.q[i])
            end
        end
    end
end

println("=== HydroPowerDynamics.jl Trollheim AGC surge comparison ===")
@printf("Load %.1f -> %.1f MW at t=%.1f s; H=%.2f s\n",P0/1e6,P1/1e6,TSTEP,H_INERTIA)
println("Solving WITHOUT surge tank..."); r0=extract(simulate_case(with_surge=false))
println("Solving WITH surge tank...");    r1=extract(simulate_case(with_surge=true))

@printf("%-25s %14s %14s\n","Metric","No surge","With surge")
@printf("%-25s %14.6f %14.6f\n","Pre f range [Hz]",r0.pre_f_range,r1.pre_f_range)
@printf("%-25s %14.6f %14.6f\n","Pre Pm range [MW]",r0.pre_pm_range,r1.pre_pm_range)
@printf("%-25s %14.6f %14.6f\n","Nadir [Hz]",r0.nadir,r1.nadir)
@printf("%-25s %14.3f %14.3f\n","Nadir time [s]",r0.t_nadir,r1.t_nadir)
@printf("%-25s %14.6f %14.6f\n","f(65) [Hz]",r0.final_f,r1.final_f)
@printf("%-25s %14.6f %14.6f\n","Pm(65) [MW]",r0.final_pm,r1.final_pm)
@printf("%-25s %14.6f %14.6f\n","IAE [Hz s]",r0.iae,r1.iae)
@printf("%-25s %14.3f %14.3f\n","Settling [s]",r0.settling,r1.settling)
@printf("Nadir delta (with - without) = %.6f Hz\n",r1.nadir-r0.nadir)

outdir=joinpath(@__DIR__,"results"); mkpath(outdir)
write_csv(joinpath(outdir,"without_surgetank.csv"),r0)
write_csv(joinpath(outdir,"with_surgetank.csv"),r1;surge=true)
open(joinpath(outdir,"metrics.csv"),"w") do io
    println(io,"metric,without_surgetank,with_surgetank")
    @printf(io,"nadir_Hz,%.10f,%.10f\n",r0.nadir,r1.nadir)
    @printf(io,"nadir_time_s,%.10f,%.10f\n",r0.t_nadir,r1.t_nadir)
    @printf(io,"final_frequency_Hz,%.10f,%.10f\n",r0.final_f,r1.final_f)
    @printf(io,"final_mechanical_power_MW,%.10f,%.10f\n",r0.final_pm,r1.final_pm)
    @printf(io,"IAE_Hz_s,%.10f,%.10f\n",r0.iae,r1.iae)
    @printf(io,"settling_s,%.10f,%.10f\n",r0.settling,r1.settling)
end
println("RESULT_DIR=",outdir)
