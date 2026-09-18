# ============================================================================ #
# connectors.jl – physical-domain interfaces for hydropower + power networks
#
# Design rule
# -----------
# Every acausal connector contains:
#   • potential/across variables: equal at a connection set
#   • flow/through variables:     sum to zero at a connection set
#
# Positive flow variables are defined as flowing INTO a component.
#
# The connector layer intentionally contains no component equations. Energy,
# constitutive laws and kinematic relations belong to components.
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Hydraulic domain – compact backward-compatible port
#
# across : p   [Pa]
# through: dm  [kg/s]
#
# Use this connector for existing HydroPowerDynamics components and for systems
# where elevation is represented by component parameters.
# ---------------------------------------------------------------------------- #
@connector HydraulicPort begin
    p(t),  [description = "Pressure [Pa]"]
    dm(t), [description = "Mass flow rate into component [kg/s]", connect = Flow]
end

# ---------------------------------------------------------------------------- #
# Hydraulic domain – elevation-aware OpenHPL-style port
#
# across : p [Pa], z [m]
# through: dm [kg/s]
#
# z is an across/topology variable. Components that change elevation (pipe,
# penstock, reservoir boundary, etc.) must explicitly relate inlet and outlet z.
# This port is introduced alongside HydraulicPort so existing models are not
# broken during the OpenHPL refactor.
# ---------------------------------------------------------------------------- #
@connector HydraulicPortZ begin
    p(t),  [description = "Pressure [Pa]"]
    z(t),  [description = "Absolute elevation [m]"]
    dm(t), [description = "Mass flow rate into component [kg/s]", connect = Flow]
end

# ---------------------------------------------------------------------------- #
# Rotational mechanical domain
#
# across : phi [rad], omega [rad/s]
# through: tau [N*m]
#
# Components carrying inertia should impose D(phi) ~ omega. Keeping both phi
# and omega on the connector is useful for generators, shafts and sensors.
# ---------------------------------------------------------------------------- #
@connector RotationalPort begin
    phi(t),   [description = "Angular position [rad]"]
    omega(t), [description = "Angular velocity [rad/s]"]
    tau(t),   [description = "Torque into component [N*m]", connect = Flow]
end

# Explicit alias used in network-oriented code.
const MechanicalRotationalPort = RotationalPort

# ---------------------------------------------------------------------------- #
# Translational mechanical domain
#
# across : x [m], v [m/s]
# through: f [N]
#
# Not central to the water-to-grid chain, but useful for servos, actuators and
# electro-mechanical auxiliary models. Components carrying mass should impose
# D(x) ~ v.
# ---------------------------------------------------------------------------- #
@connector TranslationalPort begin
    x(t), [description = "Position [m]"]
    v(t), [description = "Velocity [m/s]"]
    f(t), [description = "Force into component [N]", connect = Flow]
end

const MechanicalTranslationalPort = TranslationalPort

# ---------------------------------------------------------------------------- #
# Electrical AC network domain – rectangular phasor connector
#
# across : vr, vi [V]
# through: ir, ii [A]
#
# Rectangular voltage/current variables preserve Kirchhoff voltage/current laws
# directly and avoid making active/reactive power connector flow variables.
#
# Component terminal complex power (positive INTO component):
#
#   S = V * conj(I)
#   P = vr*ir + vi*ii
#   Q = vi*ir - vr*ii
#
# This port is suitable for balanced positive-sequence network models, SMIB,
# multi-machine networks and reduced Nordic-grid studies.
# ---------------------------------------------------------------------------- #
@connector ElectricalACPort begin
    vr(t), [description = "Real voltage component [V]"]
    vi(t), [description = "Imaginary voltage component [V]"]
    ir(t), [description = "Real current into component [A]", connect = Flow]
    ii(t), [description = "Imaginary current into component [A]", connect = Flow]
end

# Short network-oriented aliases.
const ElectricalPort = ElectricalACPort
const ACBusPort = ElectricalACPort

# ---------------------------------------------------------------------------- #
# Electrical dq domain
#
# across : vd, vq [V]
# through: id, iq [A]
#
# Useful inside synchronous-machine, converter and excitation models. Network
# interfaces should normally use ElectricalACPort unless all connected models
# share the same rotating reference frame.
# ---------------------------------------------------------------------------- #
@connector ElectricalDQPort begin
    vd(t), [description = "d-axis voltage [V]"]
    vq(t), [description = "q-axis voltage [V]"]
    id(t), [description = "d-axis current into component [A]", connect = Flow]
    iq(t), [description = "q-axis current into component [A]", connect = Flow]
end

# ---------------------------------------------------------------------------- #
# Electrical DC domain
#
# across : v [V]
# through: i [A]
#
# Useful for excitation systems, DC links, batteries and UPS/BESS extensions.
# ---------------------------------------------------------------------------- #
@connector ElectricalDCPort begin
    v(t), [description = "Voltage [V]"]
    i(t), [description = "Current into component [A]", connect = Flow]
end

# ---------------------------------------------------------------------------- #
# Causal signal domain
#
# Delegated to ModelingToolkitStandardLibrary.Blocks. These are deliberately
# separate from physical acausal connectors.
# ---------------------------------------------------------------------------- #
const SignalInPort  = Blocks.RealInput
const SignalOutPort = Blocks.RealOutput
