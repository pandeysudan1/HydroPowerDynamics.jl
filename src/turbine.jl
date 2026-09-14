# ============================================================================ #
# turbine.jl  –  Turbine models
#
# Section 3 of the HydroMTK Mathematical Reference.
#
# FrancisTurbineAffinity  – quadratic affinity-law model     (Section 3.3)
# PeltonTurbine           – jet-impact model                 (Section 3.4)
# ============================================================================ #

# ============================================================================ #
# FrancisTurbineAffinity  (Section 3.3)
# ============================================================================ #
"""
    FrancisTurbineAffinity(; name, D, K_q, eta_max, c_eta, Q_rated, rho, g)

Francis turbine – affinity-law performance model.

`RotationalPort.tau` is positive **into** a component. A turbine delivers
mechanical torque out of the turbine, so its shaft flow variable is negative.
"""
@mtkmodel FrancisTurbineAffinity begin
    @parameters begin
        D       = 2.5,    [description = "Runner diameter [m]"]
        K_q     = 0.35,   [description = "Flow coefficient at rated"]
        eta_max = 0.92,   [description = "Peak efficiency"]
        c_eta   = 0.25,   [description = "Parabolic loss coefficient"]
        Q_rated = 100.0,  [description = "Rated volumetric flow [m³/s]"]
        rho     = 1000.0, [description = "Water density [kg/m³]"]
        g       = 9.81,   [description = "Gravity [m/s²]"]
    end
    @variables begin
        H(t),         [description = "Net head across turbine [m]"]
        Q(t),         [description = "Volumetric flow rate [m³/s]"]
        eta(t),       [description = "Hydraulic efficiency [-]"]
        P_mech(t),    [description = "Mechanical power output [W]"]
        tau_shaft(t), [description = "Shaft torque magnitude [N·m]"]
        dm(t),        [description = "Mass flow rate [kg/s]"]
    end
    @components begin
        port_in  = HydraulicPort()
        port_out = HydraulicPort()
        opening  = SignalInPort()
        shaft    = RotationalPort()
    end
    @equations begin
        port_in.dm + port_out.dm ~ 0.0
        port_in.dm ~ dm
        H ~ (port_in.p - port_out.p) / (rho * g)
        Q ~ opening.u * K_q * D^2 * sqrt(abs(H) + 1e-6)
        dm ~ rho * Q
        eta ~ eta_max * (1.0 - c_eta * (Q / Q_rated - 1.0)^2)
        P_mech ~ rho * g * Q * abs(H) * eta
        tau_shaft ~ P_mech / (abs(shaft.omega) + 1e-6)

        # Positive connector torque means torque entering the turbine. Since the
        # turbine is a driver, mechanical power leaves this component.
        shaft.tau ~ -tau_shaft
    end
end

# ============================================================================ #
# PeltonTurbine  (Section 3.4)
# ============================================================================ #
"""
    PeltonTurbine(; name, R_mean, D_jet_max, C_v, phi, beta2_deg, rho, g)

Pelton turbine – jet-impact performance model. Mechanical shaft torque follows
the same connector convention as `FrancisTurbineAffinity`: output torque is
negative at the turbine port because positive flow is defined into a component.
"""
@mtkmodel PeltonTurbine begin
    @parameters begin
        R_mean    = 1.2,    [description = "Pitch circle radius [m]"]
        D_jet_max = 0.15,   [description = "Max jet diameter [m]"]
        C_v       = 0.98,   [description = "Velocity coefficient"]
        phi       = 0.85,   [description = "Relative velocity ratio"]
        beta2_deg = 170.0,  [description = "Bucket exit angle [°]"]
        rho       = 1000.0, [description = "Water density [kg/m³]"]
        g         = 9.81,   [description = "Gravity [m/s²]"]
    end
    @variables begin
        H_net(t),     [description = "Net head at nozzle [m]"]
        v_jet(t),     [description = "Jet velocity [m/s]"]
        A_jet(t),     [description = "Jet area [m²]"]
        Q(t),         [description = "Volumetric flow [m³/s]"]
        dm(t),        [description = "Mass flow rate [kg/s]"]
        u1(t),        [description = "Bucket peripheral velocity [m/s]"]
        tau_shaft(t), [description = "Shaft torque magnitude [N·m]"]
    end
    @components begin
        port_in  = HydraulicPort()
        port_out = HydraulicPort()
        opening  = SignalInPort()
        shaft    = RotationalPort()
    end
    @equations begin
        port_in.dm + port_out.dm ~ 0.0
        port_in.dm ~ dm
        H_net ~ (port_in.p - port_out.p) / (rho * g)
        v_jet ~ C_v * sqrt(2.0 * g * abs(H_net) + 1e-6)
        A_jet ~ opening.u^2 * (π / 4 * D_jet_max^2)
        Q  ~ A_jet * v_jet
        dm ~ rho * Q
        u1 ~ abs(shaft.omega) * R_mean
        tau_shaft ~ rho * Q * v_jet * u1 *
                    (1.0 + phi * cosd(beta2_deg)) / (abs(shaft.omega) + 1e-6)
        shaft.tau ~ -tau_shaft
    end
end
