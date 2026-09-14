# ============================================================================ #
# agc.jl — simple frequency droop + secondary integral AGC governor
#
# e_f      = (omega_ref - omega)/omega_ref
# d(xi)/dt = e_f
# u_cmd    = gate_bias + e_f/R + K_i*xi
# u_sat    = clamp(u_cmd, gate_min, gate_max)
# T_g*d(gate)/dt = u_sat - gate
# ============================================================================ #

"""
    SimpleGovernorAGC(; name, R, T_g, K_i, omega_ref, gate_bias, gate_min, gate_max)

Minimal hydro governor for load-frequency-control studies.

The primary response is permanent speed droop and the secondary response is an
integral of per-unit frequency error. `gate` is the servomotor state and is sent
directly to the turbine opening input.

Equations:

```math
 e_f = (\omega_{ref}-\omega)/\omega_{ref}
```

```math
 \dot\xi=e_f
```

```math
 u_{cmd}=u_0+e_f/R+K_i\xi
```

```math
 T_g\dot u = \operatorname{sat}(u_{cmd})-u
```
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
        e_f(t),            [description = "Per-unit frequency error"]
        xi(t) = 0.0,       [description = "Integral of frequency error [s]"]
        u_cmd(t),          [description = "Unsaturated gate command [pu]"]
        u_sat(t),          [description = "Saturated gate command [pu]"]
        gate(t) = 0.50,    [description = "Guide-vane servo state [pu]"]
        tau_o(t) = 0.50,   [description = "Backward-compatible gate alias [pu]"]
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
