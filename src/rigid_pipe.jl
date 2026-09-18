# ============================================================================ #
# rigid_pipe.jl – incompressible elevation-aware hydraulic pipe
#
# Physics:
#   port_a.dm + port_b.dm = 0
#   port_a.dm = dm
#
# Momentum:
#   L/A * d(dm)/dt = p_a - p_b + rho*g*(z_a-z_b) - dp_f
#
# Darcy-Weisbach:
#   v    = dm/(rho*A)
#   Re   = |dm|*D/(mu*A)
#   dp_f = f*(L/D)*(rho/2)*v*|v|
#
# This is an incompressible rigid-column model. Elastic water hammer is handled
# by a later ElasticPenstock model.
# ============================================================================ #

"""
    RigidPipe(; name, L, D_pipe, rho, g, mu, roughness)

Elevation-aware incompressible hydraulic pipe with one water-column inertia
state. The component uses two HydraulicPortZ connectors and is intended for
headrace tunnels, short rigid conduits, and control-oriented waterway models.

Positive dm is defined from port_a to port_b.
"""
@mtkmodel RigidPipe begin
    @parameters begin
        L         = 500.0,   [description = "Pipe length [m]"]
        D_pipe    = 2.0,     [description = "Internal diameter [m]"]
        rho       = 1000.0,  [description = "Water density [kg/m^3]"]
        g         = 9.81,    [description = "Gravity [m/s^2]"]
        mu        = 1.0e-3,  [description = "Dynamic viscosity [Pa*s]"]
        roughness = 1.5e-5,  [description = "Absolute roughness [m]"]
    end

    @variables begin
        dm(t) = 0.0, [description = "Mass flow from port_a to port_b [kg/s]"]
        A(t),        [description = "Cross-sectional area [m^2]"]
        Q(t),        [description = "Volume flow rate [m^3/s]"]
        v(t),        [description = "Mean velocity [m/s]"]
        Re(t),       [description = "Reynolds number [-]"]
        f_D(t),      [description = "Darcy friction factor [-]"]
        dp_f(t),     [description = "Signed Darcy pressure drop [Pa]"]
        h_f(t),      [description = "Signed Darcy head loss [m]"]
        dp_drive(t), [description = "Pressure/elevation driving term [Pa]"]
    end

    @components begin
        port_a = HydraulicPortZ()
        port_b = HydraulicPortZ()
    end

    @equations begin
        A ~ pi * D_pipe^2 / 4
        Q ~ dm / rho
        v ~ Q / A

        port_a.dm + port_b.dm ~ 0
        port_a.dm ~ dm

        Re ~ abs(dm) * D_pipe / (mu * A)
        f_D ~ darcy_factor(Re, D_pipe, roughness)
        dp_f ~ f_D * (L / D_pipe) * (rho / 2) * v * abs(v)
        h_f ~ dp_f / (rho * g)

        dp_drive ~ port_a.p - port_b.p + rho * g * (port_a.z - port_b.z)
        D(dm) ~ (A / L) * (dp_drive - dp_f)
    end
end
