# ============================================================================
# openhpl_surge_tank.jl – OpenHPL-inspired two-state surge tank
#
# This component is intentionally added alongside the reduced SurgeTank model.
# It follows the simple OpenHPL surge-shaft momentum/mass structure while using
# HydroPowerDynamics.jl's existing HydraulicPort convention.
#
# States:
#   h(t)     vertical water height above the branch [m]
#   Vdot(t)  volumetric flow into the surge shaft [m^3/s]
#
# Geometry:
#   A = pi*D^2/4
#   cos_theta = H/L
#   l = h/cos_theta
#
# Balances:
#   dm/dt = mdot
#   dM/dt = mdot*v + F_p - F_f - F_g
#
# with
#   m = rho*A*l
#   v = Vdot/A
#   M = m*v
#   mdot = rho*Vdot
#   F_p = (p_bottom-p_atm)*A
#   F_g = m*g*cos_theta
#
# Darcy wall friction is represented as
#   F_f = f*(l/D)*(rho*A/2)*v*abs(v)
#
# The connector mass flow is positive INTO the surge tank.
# ============================================================================

"""
    OpenHPLSurgeTank(; name, H, L, D_tank, h_0, Vdot_0,
                       rho, g, p_atm, f)

OpenHPL-inspired simple surge shaft with two dynamic states: water height
`h` and shaft volume flow `Vdot`.

Unlike the reduced `SurgeTank`, which imposes the port pressure
algebraically from water level, this model retains water-column momentum.
The connected hydraulic network therefore determines the bottom pressure,
which accelerates/decelerates the surge-shaft water column.

This is the first high-fidelity component in the OpenHPL refactor path.

# Parameters
- `H`: vertical component of full shaft length [m]
- `L`: full shaft length [m]
- `D_tank`: shaft diameter [m]
- `h_0`: initial vertical water height [m]
- `Vdot_0`: initial volume flow into shaft [m^3/s]
- `rho`: water density [kg/m^3]
- `g`: gravitational acceleration [m/s^2]
- `p_atm`: pressure at free surface [Pa]
- `f`: Darcy friction factor [-]

# Port
- `port`: hydraulic branch connection at the bottom of the shaft

# Notes
The model corresponds to the **simple/open surge tank** physics in OpenHPL.
Air-cushion, sharp-orifice and throttle variants are intentionally deferred.
"""
@mtkmodel OpenHPLSurgeTank begin
    @parameters begin
        H       = 87.0,      [description = "Vertical shaft component [m]"]
        L       = 87.0,      [description = "Shaft length [m]"]
        D_tank  = 3.4,       [description = "Surge shaft diameter [m]"]
        h_0     = 50.0,      [description = "Initial vertical water height [m]"]
        Vdot_0  = 0.0,       [description = "Initial shaft volume flow [m^3/s]"]
        rho     = 1000.0,    [description = "Water density [kg/m^3]"]
        g       = 9.81,      [description = "Gravity [m/s^2]"]
        p_atm   = 101_325.0, [description = "Free-surface pressure [Pa]"]
        f       = 0.02,      [description = "Darcy friction factor [-]"]
    end

    @variables begin
        h(t)    = h_0,    [description = "Vertical water height [m]"]
        Vdot(t) = Vdot_0, [description = "Volume flow into surge shaft [m^3/s]"]

        A(t),           [description = "Shaft cross-sectional area [m^2]"]
        cos_theta(t),   [description = "Shaft slope ratio H/L [-]"]
        l(t),           [description = "Water-column length along shaft [m]"]
        v(t),           [description = "Water velocity along shaft [m/s]"]
        m(t),           [description = "Water mass in shaft [kg]"]
        M(t),           [description = "Water-column momentum [kg*m/s]"]
        mdot(t),        [description = "Mass flow into shaft [kg/s]"]
        F_p(t),         [description = "Bottom pressure force [N]"]
        F_f(t),         [description = "Darcy friction force [N]"]
        F_g(t),         [description = "Gravity force along shaft [N]"]
    end

    @components begin
        port = HydraulicPort()
    end

    @equations begin
        # Geometry
        A ~ pi * D_tank^2 / 4
        cos_theta ~ H / L
        l ~ h / cos_theta

        # Kinematics / constitutive relations
        v ~ Vdot / A
        mdot ~ rho * Vdot
        m ~ rho * A * l
        M ~ m * v

        # Connector convention: positive flow enters this component.
        port.dm ~ mdot

        # Forces in the OpenHPL simple surge-shaft momentum balance.
        F_p ~ (port.p - p_atm) * A
        F_f ~ f * (l / D_tank) * (rho * A / 2) * v * abs(v)
        F_g ~ m * g * cos_theta

        # Mass balance. Since m = rho*A*h/cos(theta), this also gives dh/dt.
        D(m) ~ port.dm

        # Momentum balance.
        D(M) ~ port.dm * v + F_p - F_f - F_g
    end
end
