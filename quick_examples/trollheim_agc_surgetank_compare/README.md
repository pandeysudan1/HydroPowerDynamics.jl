# Trollheim AGC — with vs without surge tank

This example reproduces the Trollheim active-power/frequency-control study using **HydroPowerDynamics.jl** and ModelingToolkit.jl.

## Study question

How much does a surge tank change the AGC/frequency response of the same hydropower plant?

The two cases are intentionally identical except for one hydraulic branch:

```text
Case A — no surge tank
Reservoir -> headrace -> penstock -> turbine -> rotor/generator -> load

Case B — with surge tank
Reservoir -> headrace -> junction -> penstock -> turbine -> rotor/generator -> load
                              |
                              +-> surge tank
```

Both cases use the same:

- Trollheim-inspired turbine parameters,
- headrace and penstock geometry,
- rotor inertia,
- GGOV1-style governor/AGC,
- 75 MW initial load,
- 90 MW final load,
- step time `t = 5 s`,
- simulation horizon through `t = 65 s`.

A hidden preconditioning interval `-300...0 s` is used so the visible interval begins from the nonlinear operating point rather than arbitrary initial guesses.

## HydroPowerDynamics.jl components

The example reuses package components directly:

- `Reservoir`
- `Penstock`
- `SurgeTank`
- `FrancisTurbineAffinity`
- `RotorInertia`
- `SimpleGenerator`
- `RotationalSpeedSensor`
- `MechanicalPowerSensor`
- `GGOV1Governor`

The surge-tank case uses a three-way acausal hydraulic connection:

```julia
connect(headrace.port_b, penstock.port_a, surge.port)
```

so mass continuity is enforced by the connector equations.

## Disturbance

The electrical load is represented by the `SimpleGenerator.P_rated` torque demand. At `t = 5 s` it is changed from

$$
P_L = 75\ \mathrm{MW}
$$

to

$$
P_L = 90\ \mathrm{MW}.
$$

The same callback is used in both cases.

## Metrics

The script evaluates:

- pre-step frequency drift,
- pre-step mechanical-power drift,
- frequency nadir,
- nadir time,
- frequency at 65 s,
- mechanical power at 65 s,
- integral absolute frequency error,
- ±0.05 Hz settling time,
- guide-vane trajectory,
- turbine-flow trajectory,
- surge level and surge branch flow for the surge-tank case.

## Plot convention

All comparison plots use:

- **dashed line** — without surge tank,
- **solid line** — with surge tank.

Generated plots:

```text
plots/frequency_compare.png
plots/power_compare.png
plots/gate_compare.png
plots/flow_compare.png
plots/surge_level.png
plots/surge_flow.png
```

## Run

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. quick_examples/trollheim_agc_surgetank_compare/trollheim_agc_surgetank_compare.jl
```

## Model equations

### Penstock momentum

HydroPowerDynamics.jl uses

$$
\frac{d\dot m}{dt}
=
\frac{A}{L}
\left(p_{in}-p_{out}-\Delta p_f\right),
$$

with Darcy-Weisbach friction.

### Surge tank

$$
A_t\frac{dZ}{dt}=\frac{\dot m_s}{\rho},
$$

and

$$
p_j=p_{atm}+\rho g Z-\Delta p_{riser}.
$$

### Francis turbine

$$
Q=u_g K_q D^2\sqrt{|H|},
$$

$$
P_m=\rho gQ|H|\eta(Q).
$$

### Rotor

$$
J\dot\omega=\tau_t-\tau_g.
$$

### GGOV1-style governor

$$
P_{ref}=P_{set}+\frac{\omega_{ref}-\omega}{R},
$$

$$
T_p\dot P_m^{filt}=P_m-P_m^{filt},
$$

$$
T_g\dot x_g=e_g-x_g,
$$

$$
\dot x_i=K_i e_g,
$$

$$
u_g=\operatorname{sat}(x_g+x_i).
$$

## Interpretation

The purpose is not to assume that the surge tank always improves the frequency nadir. Instead, the simulation isolates how surge storage changes the hydraulic trajectory seen by the governor. Final steady-state power should be similar; the differences should appear mainly in the transient flow, gate demand, turbine-power trajectory, nadir, and settling behavior.
