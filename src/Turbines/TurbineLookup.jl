# ============================================================================
# TurbineLookup.jl
#
# Data-driven Francis-style turbine model using a compact 3x3 hill chart.
#
# Independent lookup coordinates:
#   gate  : guide-vane opening [-]
#   n11   : unit speed = n_rpm * D / sqrt(H)
#
# Lookup outputs:
#   Q11   : unit discharge
#   eta   : hydraulic efficiency
#
# Scaling:
#   Q      = Q11 * D^2 * sqrt(H)
#   P_hyd  = rho*g*Q*H
#   P_mech = eta*P_hyd
#
# The default table is a smooth Trollheim-inspired demonstration table.
# It is NOT claimed to be measured Trollheim hill-chart data. Replace the
# scalar table parameters with identified/field values when available.
# ============================================================================

"""
    bilinear3x3_clamped(x, y,
                        x1,x2,x3, y1,y2,y3,
                        z11,z12,z13,z21,z22,z23,z31,z32,z33)

Clamped bilinear interpolation on a 3x3 rectangular grid.

Rows correspond to x-grid values and columns correspond to y-grid values.
The scalar-argument form is intentional: it is straightforward for
ModelingToolkit symbolic registration and parameterization.
"""
function bilinear3x3_clamped(x, y,
    x1,x2,x3, y1,y2,y3,
    z11,z12,z13,z21,z22,z23,z31,z32,z33)

    xc = min(max(x, x1), x3)
    yc = min(max(y, y1), y3)

    if xc <= x2
        xa, xb = x1, x2
        if yc <= y2
            ya, yb = y1, y2
            za, zb, zc, zd = z11, z12, z21, z22
        else
            ya, yb = y2, y3
            za, zb, zc, zd = z12, z13, z22, z23
        end
    else
        xa, xb = x2, x3
        if yc <= y2
            ya, yb = y1, y2
            za, zb, zc, zd = z21, z22, z31, z32
        else
            ya, yb = y2, y3
            za, zb, zc, zd = z22, z23, z32, z33
        end
    end

    tx = (xc - xa) / (xb - xa)
    ty = (yc - ya) / (yb - ya)

    zlow  = za + ty * (zb - za)
    zhigh = zc + ty * (zd - zc)
    return zlow + tx * (zhigh - zlow)
end

@register_symbolic bilinear3x3_clamped(
    x, y,
    x1,x2,x3, y1,y2,y3,
    z11,z12,z13,z21,z22,z23,z31,z32,z33
)

"""
    TurbineLookup(; name, ...)

Acausal hydraulic-to-rotational turbine based on a compact hill-chart lookup.

The model starts from conservation/energy conversion:

    port_a.dm + port_b.dm = 0
    P_hyd = rho*g*Q*H
    P_mech = eta*P_hyd

The hydraulic characteristic comes from lookup tables in the standard turbine
unit quantities n11 and Q11.

Default lookup axes:
- gate = [0.2, 0.6, 1.0]
- n11  = [40, 50, 60]

The default table is only a demonstration benchmark. For plant validation,
replace Q11 and eta entries with measured or manufacturer hill-chart data.
"""
@mtkmodel TurbineLookup begin
    @parameters begin
        rho      = 1000.0, [description = "Water density [kg/m^3]"]
        g        = 9.81,   [description = "Gravity [m/s^2]"]
        D_runner = 2.5,    [description = "Runner diameter [m]"]
        H_eps    = 1.0e-4, [description = "Head regularization [m]"]
        omega_eps = 1.0e-3,[description = "Speed regularization [rad/s]"]

        # Lookup axes
        gate_1 = 0.2
        gate_2 = 0.6
        gate_3 = 1.0
        n11_1  = 40.0
        n11_2  = 50.0
        n11_3  = 60.0

        # Q11 table: rows = gate, columns = n11
        q11_11 = 0.060
        q11_12 = 0.062
        q11_13 = 0.060
        q11_21 = 0.185
        q11_22 = 0.190
        q11_23 = 0.184
        q11_31 = 0.300
        q11_32 = 0.30735108593749
        q11_33 = 0.298

        # Efficiency table: rows = gate, columns = n11
        eta_11 = 0.72
        eta_12 = 0.76
        eta_13 = 0.70
        eta_21 = 0.89
        eta_22 = 0.93
        eta_23 = 0.88
        eta_31 = 0.94
        eta_32 = 0.97
        eta_33 = 0.92
    end

    @variables begin
        H(t),         [description = "Net hydraulic head [m]"]
        H_eff(t),     [description = "Regularized positive head [m]"]
        n_rpm(t),     [description = "Runner speed [rpm]"]
        n11(t),       [description = "Unit speed [rpm*m/sqrt(m)]"]
        Q11(t),       [description = "Unit discharge"]
        eta(t),       [description = "Hydraulic efficiency [-]"]
        Q(t),         [description = "Volume flow [m^3/s]"]
        mdot(t),      [description = "Mass flow [kg/s]"]
        P_hyd(t),     [description = "Hydraulic power [W]"]
        P_mech(t),    [description = "Mechanical shaft power [W]"]
        tau_mech(t),  [description = "Mechanical torque magnitude [N*m]"]
    end

    @components begin
        port_a = HydraulicPortZ()
        port_b = HydraulicPortZ()
        shaft  = RotationalPort()
        gate   = SignalInPort()
    end

    @equations begin
        # 1. MASS / ENERGY BALANCE
        port_a.dm + port_b.dm ~ 0
        P_hyd  ~ rho * g * Q * H_eff
        P_mech ~ eta * P_hyd

        # 2. ALGEBRAIC / LOOKUP RELATIONS
        H ~ (port_a.p - port_b.p) / (rho * g) + (port_a.z - port_b.z)
        H_eff ~ max(H, H_eps)

        n_rpm ~ abs(shaft.omega) * 60 / (2*pi)
        n11 ~ n_rpm * D_runner / sqrt(H_eff)

        Q11 ~ bilinear3x3_clamped(
            gate.u, n11,
            gate_1, gate_2, gate_3,
            n11_1, n11_2, n11_3,
            q11_11, q11_12, q11_13,
            q11_21, q11_22, q11_23,
            q11_31, q11_32, q11_33
        )

        eta ~ bilinear3x3_clamped(
            gate.u, n11,
            gate_1, gate_2, gate_3,
            n11_1, n11_2, n11_3,
            eta_11, eta_12, eta_13,
            eta_21, eta_22, eta_23,
            eta_31, eta_32, eta_33
        )

        Q ~ Q11 * D_runner^2 * sqrt(H_eff)
        mdot ~ rho * Q
        tau_mech ~ P_mech * shaft.omega / (shaft.omega^2 + omega_eps^2)

        # 3. CONNECTION EQUATIONS
        port_a.dm ~ mdot
        port_b.z ~ port_a.z

        # Turbine exports mechanical torque, hence negative torque INTO turbine.
        shaft.tau ~ -tau_mech
    end
end
