# ============================================================================ #
# utils.jl  –  System-level equations and helper functions
#
# Equations from Section 6 & 7 of the HydroMTK Mathematical Reference.
# ============================================================================ #

# ── Friction factor models ───────────────────────────────────────────────────

"""
    darcy_factor(Re, D, e) -> f_D

Multi-regime Darcy-Weisbach friction factor (dimensionless).

Three regimes, C¹-continuous:
- Laminar     Re ≤ 2100 : f = 64/Re                    (Hagen-Poiseuille)
- Transition  2100<Re<2300 : cubic Hermite spline       (matched value + slope)
- Turbulent   Re ≥ 2300 : Swamee-Jain variant           (Colebrook approximation)

Registered with ModelingToolkit for use inside `@mtkmodel` equation blocks.
"""
function darcy_factor(Re, D, e)
    Re_l = 2100.0
    Re_t = 2300.0

    f_l(r) = 64.0 / r
    f_t(r) = 1.0 / (2.0 * log10(e / (3.7 * D) + 5.7 / (r^0.9)))^2

    Re <= 0    && return 0.0
    Re <= Re_l && return f_l(Re)
    Re >= Re_t && return f_t(Re)

    # Cubic Hermite spline — matches value and first derivative at both ends
    y1  = f_l(Re_l)
    y2  = f_t(Re_t)
    dy1 = -64.0 / Re_l^2
    δ   = 1e-3 * Re_t
    dy2 = (f_t(Re_t + δ) - f_t(Re_t - δ)) / (2δ)
    X   = [
        Re_l^3   Re_l^2  Re_l  1.0
        Re_t^3   Re_t^2  Re_t  1.0
        3Re_l^2  2Re_l   1.0   0.0
        3Re_t^2  2Re_t   1.0   0.0
    ]
    k = X \ [y1, y2, dy1, dy2]
    return k[1]*Re^3 + k[2]*Re^2 + k[3]*Re + k[4]
end
@register_symbolic darcy_factor(Re, D, e)

"""
    swamee_jain(Re, D, e) -> f_D

Swamee-Jain explicit approximation to the Colebrook-White equation.
Valid for: 5×10³ ≤ Re ≤ 10⁸, 10⁻⁶ ≤ e/D ≤ 10⁻².
Error < 3% vs Colebrook-White.
"""
swamee_jain(Re, D, e) =
    0.25 / (log10(e / (3.7 * D) + 5.74 / (Re^0.9 + 1e-10)))^2
@register_symbolic swamee_jain(Re, D, e)

# ── Fluid constants (defaults) ───────────────────────────────────────────────
const RHO_0  = 1000.0    # kg/m³  reference density
const B_BULK = 2.0e9     # Pa     bulk modulus of water
const MU_W   = 1.0e-3    # Pa·s   dynamic viscosity
const P_ATM  = 101_325.0 # Pa     atmospheric pressure
const G_GRAV = 9.81      # m/s²   gravitational acceleration

# ── 6.1  Gross and Net Head  [m] ─────────────────────────────────────────────
"""
    gross_head(Z_reservoir, Z_tailwater)

Return H_gross = Z_reservoir - Z_tailwater  [m].  (Eq. SL.1)
"""
gross_head(Z_res, Z_tail) = Z_res - Z_tail

"""
    net_head(H_gross, h_f_penstock, h_f_draft = 0.0)

Return H_net = H_gross - h_f_penstock - h_f_draft  [m].  (Eq. SL.2)
"""
net_head(H_gross, h_f_p, h_f_d = 0.0) = H_gross - h_f_p - h_f_d

"""
    darcy_head_loss(f_D, L, D, v, g = G_GRAV)

Darcy-Weisbach friction head loss:  h_f = f_D*(L/D)*v²/(2g)  [m].  (Eq. SL.3)
"""
darcy_head_loss(f_D, L, D, v, g = G_GRAV) = f_D * (L / D) * v^2 / (2g)

# ── 6.2  Power Cascade  [W] ──────────────────────────────────────────────────
"""
    power_cascade(rho, g, Q, H_gross, H_net, eta_t, eta_gen, eta_xfmr)

Return a `NamedTuple` with the full power cascade:
  `P_hydraulic`, `P_available`, `P_turbine`, `P_electrical`, `P_grid`
and overall `eta_plant`.  (Eqs. PC.1–PC.7)
"""
function power_cascade(rho, g, Q, H_gross, H_net, eta_t, eta_gen, eta_xfmr = 1.0)
    P_hydraulic  = rho * g * Q * H_gross
    P_available  = rho * g * Q * H_net
    P_turbine    = P_available * eta_t
    P_electrical = P_turbine  * eta_gen
    P_grid       = P_electrical * eta_xfmr
    eta_plant    = P_grid / P_hydraulic
    return (; P_hydraulic, P_available, P_turbine, P_electrical, P_grid, eta_plant)
end

"""
    plant_efficiency(eta_hyd, eta_mech, eta_gen, eta_xfmr = 1.0)

Overall plant efficiency as product of component efficiencies.  (Eq. PC.7)
"""
plant_efficiency(eta_hyd, eta_mech, eta_gen, eta_xfmr = 1.0) =
    eta_hyd * eta_mech * eta_gen * eta_xfmr

# ── 6.3  Water Hammer  ───────────────────────────────────────────────────────
"""
    wave_speed(B, rho, D, E_pipe, e_wall)

Pressure-wave celerity in a buried pipeline:
    a = sqrt(B / (rho * (1 + B*D/(E_pipe*e_wall))))  [m/s].  (Eq. WH.2)
"""
wave_speed(B, rho, D, E_pipe, e_wall) =
    sqrt(B / (rho * (1 + B * D / (E_pipe * e_wall))))

"""
    joukowsky_pressure(rho, a_wave, Delta_v)

Maximum instantaneous pressure rise from gate closure:
    Δp = ρ · a · Δv  [Pa].  (Eq. WH.1)
"""
joukowsky_pressure(rho, a_wave, Delta_v) = rho * a_wave * abs(Delta_v)

"""
    critical_closure_time(L, a_wave)

Critical gate-closure time  T_c = 2L/a  [s].  (Eq. WH.3)
Closures faster than T_c experience full Joukowsky pressure.
"""
critical_closure_time(L, a_wave) = 2L / a_wave

# ── 3.1  Dimensionless turbine parameters ────────────────────────────────────
"""
    unit_speed(n_rpm, D, H)

Unit speed  n₁₁ = n·D/√H  [rpm·m^(1/2)].  (Eq. T.1)
"""
unit_speed(n, D, H) = n * D / sqrt(H)

"""
    unit_discharge(Q, D, H)

Unit discharge  Q₁₁ = Q/(D²·√H)  [(m^(1/2))/s].  (Eq. T.2)
"""
unit_discharge(Q, D, H) = Q / (D^2 * sqrt(H))

"""
    hydraulic_efficiency(P_mech, rho, g, Q, H)

η = P_mech / (ρ·g·Q·H).  (Eq. T.4)
"""
hydraulic_efficiency(P_mech, rho, g, Q, H) = P_mech / (rho * g * Q * H)


# ── Model structure / well-posedness diagnostics ─────────────────────────────

"""
    model_structure_report(sys)

Return a compact structural count for a ModelingToolkit system.

Fields:
- equations   : number of equations
- unknowns    : number of unknown variables
- parameters  : number of parameters
- balanced    : true when equations == unknowns

This is a bookkeeping diagnostic, not a substitute for successful
structural compilation and initialization.
"""
function model_structure_report(sys)
    neq = length(equations(sys))
    nunk = length(unknowns(sys))
    npar = length(parameters(sys))
    return (
        equations = neq,
        unknowns = nunk,
        parameters = npar,
        balanced = neq == nunk,
    )
end
