# ============================================================================ #
# agc.jl — simple frequency droop + secondary integral AGC governor
# ============================================================================ #

"""
    SimpleGovernorAGC(; name, R, T_g, K_i, omega_ref, gate_bias, gate_min, gate_max)

Minimal hydro governor for load-frequency-control studies.

The component intentionally leaves dynamic-state initial conditions to the
assembled plant model. In particular, `tau_o` is an algebraic alias of `gate`
and therefore carries no independent initialization binding.
"""
@mtkmodel SimpleGovernorAGC begin
    @parameters begin
        R         = 0.50,   [description = "Permanent speed droop [pu/pu]"]
        T_g       = 0.30,   [description = "Guide-vane servo time constant [s]"]
        K_i       = 0.10,   [description = "Secondary integral gain [1/s]"]
        omega_ref = 157.08, [description = "Reference shaft speed [rad/s]"]
        gate_bias = 0.50,   [description = "Nominal gate bias [pu]"]
        gate_min  = 0.05,   [description = "Minimum gate opening [pu]"]
        gate_max  = 1.00,   [description = "Maximum gate opening [pu]"]
    end
    @variables begin
        e_f(t)
        xi(t)
        u_cmd(t)
        u_sat(t)
        gate(t)
        tau_o(t)
    end
    @components begin
        speed_in = Blocks.RealInput()
        gate_out = Blocks.RealOutput()
    end
    @equations begin
        e_f ~ (omega_ref - speed_in.u) / omega_ref
        D(xi) ~ e_f
        u_cmd ~ gate_bias + e_f / R + K_i * xi
        u_sat ~ min(gate_max, max(gate_min, u_cmd))
        T_g * D(gate) ~ u_sat - gate
        gate_out.u ~ gate
        tau_o ~ gate
    end
end

# ============================================================================ #
# VariableLoadGenerator — generator reaction torque driven by an input load
# ============================================================================ #
"""
    VariableLoadGenerator(; name, omega_rated, omega_s, D_d, eta_gen)

Simple algebraic generator/load model with an external active-power demand input.
`RotationalPort.tau` is positive into the component, so a generator/load absorbs
positive shaft torque.
"""
@mtkmodel VariableLoadGenerator begin
    @parameters begin
        omega_rated = 157.08, [description = "Rated shaft speed [rad/s]"]
        omega_s     = 157.08, [description = "Synchronous shaft speed [rad/s]"]
        D_d         = 0.0,    [description = "Speed-dependent load damping [-]"]
        eta_gen     = 1.0,    [description = "Generator efficiency [-]"]
    end
    @variables begin
        tau_gen(t)
        P_elec(t)
    end
    @components begin
        flange = RotationalPort()
        P_load = Blocks.RealInput()
    end
    @equations begin
        tau_gen ~ (P_load.u / omega_rated) *
                  (1.0 + D_d * (flange.omega - omega_s) / omega_s)
        P_elec ~ tau_gen * abs(flange.omega) * eta_gen
        flange.tau ~ tau_gen
    end
end
