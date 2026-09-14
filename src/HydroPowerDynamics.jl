"""
    HydroPowerDynamics.jl

Acausal, equation-based hydropower component library built on ModelingToolkit.jl.

Physical domains covered
  • Hydraulic  (pressure / mass-flow)
  • Mechanical rotational (angular velocity / torque)
  • Signal  (causal scalar)

Component families
  • Hydraulic  : Reservoir, Penstock, SurgeTank, DraftTube, GuideVane
  • Turbine    : FrancisTurbineAffinity, PeltonTurbine
  • Mechanical : RotorInertia, SimpleGenerator, VariableLoadGenerator
  • Governor   : PIDGovernor, GGOV1Governor, SimpleGovernorAGC
"""
module HydroPowerDynamics

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
import OrdinaryDiffEq
using DataInterpolations
using ModelingToolkitStandardLibrary.Blocks

export t, D
export Blocks

include("utils.jl")
include("connectors.jl")
include("hydraulic.jl")
include("turbine.jl")
include("mechanical.jl")
include("governor.jl")
include("agc.jl")

export HydraulicPort, RotationalPort, SignalInPort, SignalOutPort
export Reservoir, Penstock, SurgeTank, DraftTube, GuideVane
export FrancisTurbineAffinity, PeltonTurbine
export RotorInertia, SimpleGenerator, VariableLoadGenerator, RotationalSpeedSensor, MechanicalPowerSensor
export PIDGovernor, GGOV1Governor, SimpleGovernorAGC
export gross_head, net_head, darcy_head_loss, power_cascade, joukowsky_pressure, wave_speed
export critical_closure_time, unit_speed, unit_discharge, plant_efficiency, hydraulic_efficiency
export darcy_factor, swamee_jain

end # module HydroPowerDynamics
