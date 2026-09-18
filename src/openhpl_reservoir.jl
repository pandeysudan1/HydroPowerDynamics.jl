# ============================================================================
# openhpl_reservoir.jl – elevation-aware reservoir models inspired by OpenHPL
#
# Two models are provided deliberately:
#
#   ReservoirBoundaryZ
#       Constant-level hydrostatic boundary for short-timescale dynamics.
#
#   DynamicReservoir
#       Variable-level storage reservoir with OpenHPL-style geometry and
#       mass balance. This represents the default OpenHPL reservoir mode.
#
# Both use HydraulicPortZ so downstream waterway components can propagate
# absolute elevation explicitly.
# ============================================================================

"""
    ReservoirBoundaryZ(; name, h, z_out, rho, g, p_atm)

Elevation-aware constant-level upstream reservoir.

Hydrostatic outlet pressure:

    p_out = p_atm + rho*g*h

Absolute outlet elevation is imposed as:

    z = z_out

This is the preferred reservoir for FCR and short transient studies where
reservoir level is effectively constant.
"""
@mtkmodel ReservoirBoundaryZ begin
    @parameters begin
        h      = 100.0,     [description = "Water level above outlet [m]"]
        z_out  = 0.0,       [description = "Absolute outlet elevation [m]"]
        rho    = 1000.0,    [description = "Water density [kg/m^3]"]
        g      = 9.81,      [description = "Gravity [m/s^2]"]
        p_atm  = 101_325.0, [description = "Atmospheric pressure [Pa]"]
    end

    @variables begin
        h_abs(t), [description = "Absolute free-surface elevation [m]"]
    end

    @components begin
        port = HydraulicPortZ()
    end

    @equations begin
        port.z ~ z_out
        h_abs ~ z_out + h
        port.p ~ p_atm + rho * g * h
    end
end


"""
    DynamicReservoir(; name, h_0, z_out, L, W, alpha, rho, g, p_atm)

Variable-level storage reservoir based on the default OpenHPL reservoir model.

Geometry follows OpenHPL:

    A_vertical = h * (W + h*tan(alpha))
    m          = rho * A_vertical * L

The single hydraulic outlet uses the HydroPowerDynamics sign convention:
positive connector mass flow is INTO the reservoir. Therefore reservoir mass
balance is simply:

    dm/dt = port.dm

For normal generation flow, water leaves the reservoir and `port.dm < 0`, so
the stored mass and water level decrease.

Outlet pressure remains hydrostatic:

    p_out = p_atm + rho*g*h

This component intentionally captures the robust default OpenHPL mode first.
Reservoir inflow/momentum/friction mode will be added as a higher-fidelity
variant after the basic waterway chain is validated.
"""
@mtkmodel DynamicReservoir begin
    @parameters begin
        h_0    = 50.0,      [description = "Initial level above outlet [m]"]
        z_out  = 0.0,       [description = "Absolute outlet elevation [m]"]
        L      = 500.0,     [description = "Reservoir length [m]"]
        W      = 100.0,     [description = "Reservoir bed width [m]"]
        alpha  = 0.0,       [description = "Side-wall angle [rad]"]
        rho    = 1000.0,    [description = "Water density [kg/m^3]"]
        g      = 9.81,      [description = "Gravity [m/s^2]"]
        p_atm  = 101_325.0, [description = "Atmospheric pressure [Pa]"]
    end

    @variables begin
        h(t) = h_0, [description = "Water level above outlet [m]"]
        A(t),       [description = "Vertical reservoir cross-sectional area [m^2]"]
        V(t),       [description = "Stored water volume [m^3]"]
        m(t),       [description = "Stored water mass [kg]"]
        h_abs(t),   [description = "Absolute free-surface elevation [m]"]
        Q_out(t),   [description = "Positive volumetric outflow [m^3/s]"]
    end

    @components begin
        port = HydraulicPortZ()
    end

    @equations begin
        # Geometry from OpenHPL/Waterway/Reservoir.mo
        A ~ h * (W + h * tan(alpha))
        V ~ A * L
        m ~ rho * V

        # Elevation and hydrostatic outlet pressure
        port.z ~ z_out
        h_abs ~ z_out + h
        port.p ~ p_atm + rho * g * h

        # Positive Q_out means water leaving the reservoir.
        Q_out ~ -port.dm / rho

        # Storage balance. Connector flow is positive INTO component.
        D(m) ~ port.dm
    end
end
