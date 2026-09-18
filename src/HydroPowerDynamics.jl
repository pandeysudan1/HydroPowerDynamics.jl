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

include("Core/Utils.jl")
include("Interfaces/Connectors.jl")

include("Waterways/LegacyHydraulic.jl")
include("Waterways/Reservoirs/Reservoirs.jl")
include("Waterways/RigidPipe.jl")
include("Waterways/SurgeTank.jl")

include("Turbines/LegacyTurbines.jl")
include("Turbines/TurbineLookup.jl")

include("Electromechanical/Mechanical.jl")
include("Controls/Governors.jl")
include("Controls/AGC.jl")

export HydraulicPort, HydraulicPortZ
export RotationalPort, MechanicalRotationalPort
export TranslationalPort, MechanicalTranslationalPort
export ElectricalACPort, ElectricalPort, ACBusPort, ElectricalDQPort, ElectricalDCPort
export SignalInPort, SignalOutPort
export Reservoir, ReservoirBoundaryZ, DynamicReservoir
export RigidPipe, Penstock, SurgeTank, OpenHPLSurgeTank, DraftTube, GuideVane
export FrancisTurbineAffinity, PeltonTurbine, TurbineLookup, bilinear3x3_clamped
export RotorInertia, SimpleGenerator, VariableLoadGenerator, RotationalSpeedSensor, MechanicalPowerSensor
export PIDGovernor, GGOV1Governor, SimpleGovernorAGC
export gross_head, net_head, darcy_head_loss, power_cascade, joukowsky_pressure, wave_speed
export critical_closure_time, unit_speed, unit_discharge, plant_efficiency, hydraulic_efficiency
export darcy_factor, swamee_jain, model_structure_report

end # module HydroPowerDynamics
