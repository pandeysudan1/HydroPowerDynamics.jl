# ============================================================================ #
# mechanical.jl  –  Rotating mass and generator models
# ============================================================================ #

# ============================================================================ #
# RotorInertia
#
# RotationalPort convention: tau > 0 means torque entering a component.
# Therefore, for the rotor storage element,
#
#   J*dω/dt = τ_turbine_node + τ_generator_node - b_v*ω
#
# With a turbine driver, its own port has tau < 0 and the rotor-side port gets
# tau > 0 through the connection equation. With a generator load, its port has
# tau > 0 and the rotor-side port gets tau < 0.
# ============================================================================ #
"""
    RotorInertia(; name, J, b_v, omega_0)

Lumped rotating mass (turbine runner + shaft + generator rotor).
"""
@mtkmodel RotorInertia begin
    @parameters begin
        J       = 5000.0,  [description = "Polar inertia [kg·m²]"]
        b_v     = 5.0,     [description = "Viscous damping [N·m·s/rad]"]
        omega_0 = 0.0,     [description = "Initial angular velocity [rad/s]"]
    end
    @variables begin
        omega(t) = 0.0,  [description = "Angular velocity [rad/s]"]
        theta(t) = 0.0,  [description = "Angular position [rad]"]
    end
    @components begin
        flange_turbine   = RotationalPort()
        flange_generator = RotationalPort()
    end
    @equations begin
        flange_turbine.omega   ~ omega
        flange_generator.omega ~ omega
        flange_turbine.phi     ~ theta
        flange_generator.phi   ~ theta
        D(theta) ~ omega
        J * D(omega) ~ flange_turbine.tau + flange_generator.tau - b_v * omega
    end
end

# ============================================================================ #
# SimpleGenerator
# ============================================================================ #
"""
    SimpleGenerator(; name, P_rated, omega_rated, omega_s, D_d, eta_gen)

Simplified algebraic generator. `flange.tau` is positive because mechanical
shaft torque enters the generator/load component.
"""
@mtkmodel SimpleGenerator begin
    @parameters begin
        P_rated     = 50.0e6,  [description = "Rated power [W]"]
        omega_rated = 157.08,  [description = "Rated speed [rad/s]"]
        omega_s     = 157.08,  [description = "Synchronous speed [rad/s]"]
        D_d         = 1.5,     [description = "Damping coefficient"]
        eta_gen     = 0.97,    [description = "Generator efficiency"]
    end
    @variables begin
        tau_gen(t), [description = "Generator torque magnitude [N·m]"]
        P_elec(t),  [description = "Electrical output power [W]"]
    end
    @components begin
        flange = RotationalPort()
    end
    @equations begin
        tau_gen ~ (P_rated / omega_rated) *
                  (1.0 + D_d * (flange.omega - omega_s) / omega_s)
        P_elec ~ tau_gen * abs(flange.omega) * eta_gen
        flange.tau ~ tau_gen
    end
end

# ============================================================================ #
# RotationalSpeedSensor
# ============================================================================ #
"""
    RotationalSpeedSensor(; name)

Ideal non-invasive speed sensor.
"""
@mtkmodel RotationalSpeedSensor begin
    @components begin
        flange = RotationalPort()
        w      = SignalOutPort()
    end
    @equations begin
        flange.tau ~ 0
        w.u ~ flange.omega
    end
end

# ============================================================================ #
# MechanicalPowerSensor
# ============================================================================ #
"""
    MechanicalPowerSensor(; name)

Ideal in-line rotational power sensor. Positive `P.u` means power flows from
`flange_a` toward `flange_b`.
"""
@mtkmodel MechanicalPowerSensor begin
    @components begin
        flange_a = RotationalPort()
        flange_b = RotationalPort()
        P        = SignalOutPort()
    end
    @equations begin
        flange_b.tau ~ -flange_a.tau
        flange_b.omega ~ flange_a.omega
        flange_b.phi ~ flange_a.phi
        # Power entering through side a. For a turbine upstream, the sensor-side
        # torque is positive after the connection equation.
        P.u ~ flange_a.tau * flange_a.omega
    end
end
