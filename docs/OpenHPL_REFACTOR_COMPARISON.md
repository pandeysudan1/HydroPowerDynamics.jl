# HydroPowerDynamics.jl → OpenHPL refactor comparison

Branch: `Refractor_OpenHPL`

Reference repositories:
- MTK implementation: `pandeysudan1/HydroPowerDynamics.jl`
- Modelica reference: `pandeysudan1/OpenHPL`

## Objective

Keep the composability and symbolic-analysis advantages of ModelingToolkit while progressively replacing simplified hydraulic/electromechanical models with equations that reproduce the physical structure of OpenHPL.

The first target is not a line-by-line Modelica translation. It is a component-by-component equivalence layer with matched ports, states, parameters, initialisation assumptions, and validation tests.

## 1. Component map

| HydroPowerDynamics.jl | OpenHPL reference | Current MTK fidelity | Main gap |
|---|---|---|---|
| `Reservoir` | `Waterway/Reservoir.mo` | Low/medium | MTK is fixed-head only; OpenHPL supports variable level, inflow, momentum, friction and topology/elevation |
| `Penstock` | `Waterway/Penstock.mo`, `PenstockKP.mo`, `Pipe.mo` | Medium for lumped studies | OpenHPL uses distributed compressible/elastic dynamics and gravity/elevation; MTK is essentially one lumped momentum state |
| `SurgeTank` | `Waterway/SurgeTank.mo` | Medium | MTK has level + empirical riser loss; OpenHPL includes water momentum, moving mass, gravity and several surge-tank types |
| `DraftTube` | `Waterway/DraftTube.mo` | Medium | MTK is algebraic Bernoulli/loss model; OpenHPL is more detailed and geometry-aware |
| `GuideVane` | `Waterway/Gate.mo` / turbine guide-vane equations | Medium | MTK uses a simple orifice; OpenHPL can represent gate geometry/nonlinear opening relationships |
| `FrancisTurbineAffinity` | `ElectroMech/Turbines/Francis.mo` | Low for nonlinear transient physics | OpenHPL uses Euler turbine equations, runner geometry, guide-vane kinematics and loss terms |
| `PeltonTurbine` | `ElectroMech/Turbines/Pelton.mo` | Medium | Needs parameter/equation alignment and validation |
| `RotorInertia` + `SimpleGenerator` | `Generators/SimpleGen.mo`, base classes | Medium | MTK splits inertia/generator; OpenHPL couples energy balance, friction and electrical load |
| `PIDGovernor` / `GGOV1Governor` | `Controllers/Governor.mo` | Different model family | Current models do not reproduce OpenHPL transient-droop governor topology |
| `HydraulicPort` | `Interfaces.Contact_i/o` | Good basic analogy | OpenHPL also propagates elevation/topology information |

## 2. Connector philosophy

Both libraries use an acausal hydraulic connector concept:

- potential variable: pressure `p`
- flow variable: mass flow `mdot`

This is already a strong basis for translation.

OpenHPL additionally propagates elevation through its connectors. For a closer translation, HydroPowerDynamics should add elevation explicitly, either:

1. directly to `HydraulicPort`, or
2. through a separate elevation connector/parameter layer.

For the first refactor pass, use a component parameter `Δz` / `H` to avoid making every existing model incompatible.

## 3. Reservoir

### Current MTK

The current model is a pressure boundary:

```text
p_out = p_atm + ρ g H
```

This is appropriate for short FCR simulations where reservoir level is effectively constant.

### OpenHPL

OpenHPL can represent:

- variable reservoir water level,
- inflow/outflow mass balance,
- changing cross-sectional area,
- reservoir momentum,
- wall/bed friction,
- hydrostatic pressure,
- fixed-level mode.

### Refactor direction

Do not remove the simple reservoir. Introduce two levels:

- `ReservoirBoundary`: current constant-head model
- `DynamicReservoir`: OpenHPL-equivalent mass/momentum model

This preserves fast FCR studies while allowing nonlinear waterway studies.

## 4. Penstock

### Current MTK

The present `Penstock` uses a single flow state and Darcy loss:

```text
d(mdot)/dt = A/L · (p_in - p_out - Δp_f)
```

It also defines `p_avg`, but because the two-port flow balance enforces
`mdot_in + mdot_out = 0`, the present compressibility equation does not create a true distributed water-hammer state.

### OpenHPL

OpenHPL's penstock formulation contains:

- multiple spatial cells,
- pressure states,
- flow states,
- water compressibility,
- pipe-wall elasticity through effective compressibility,
- gravity due to elevation drop,
- Darcy friction,
- pressure-dependent density and area.

Its older staggered-grid `Penstock.mo` is explicitly marked as less robust; `PenstockKP.mo` is the more important high-fidelity reference.

### Refactor direction

Create three models rather than forcing one model to do everything:

- `RigidPipe`: low-order incompressible inertance + friction
- `ElasticPenstockLumped`: 2–4 state water-hammer approximation
- `ElasticPenstockFV`: N-cell finite-volume / semi-discrete OpenHPL-equivalent model

FCR model identification can use the first two; pressure-wave/cavitation work can use the third.

## 5. Surge tank

### Current MTK

Current state:

```text
dZ/dt = mdot / (ρ A_t)
p_port = p_atm + ρ g Z - Δp_riser
```

This is a useful reduced model but has no explicit water-column momentum state.

### OpenHPL

OpenHPL models:

```text
M = m v
dM/dt = momentum_flux + pressure_force - friction - gravity
dm/dt = net_mass_flow
```

and supports:

- simple tank,
- air-cushion tank,
- sharp-orifice tank,
- throttle-valve tank,
- creek inflow,
- steady-state initialisation.

### Refactor direction

The first major physics upgrade should be an `OpenHPLSurgeTank` with two differential states:

- water level / mass,
- water-column momentum or riser flow.

This is especially important for transient droop and water-mass oscillation studies.

## 6. Francis turbine

### Current MTK

`FrancisTurbineAffinity` uses:

```text
Q = opening · Kq · D² · sqrt(|H|)
P = ρ g Q H η
τ = P / ω
```

This is compact and excellent for control development, but turbine dynamics are dominated by an empirical affinity relation.

### OpenHPL

The detailed Francis model contains:

- Euler turbine power,
- inlet/outlet runner geometry,
- guide-vane linkage and angle,
- blade velocity triangles,
- runner pressure drop,
- hydraulic loss terms,
- optional water compressibility,
- low-load behavior.

Representative structure:

```text
W_s =
  mdot·ω·R1·(Vdot/A1)·cot(α1)
  - mdot·ω·R2·(ωR2 + (Vdot/A2)·cot(β2))
```

### Refactor direction

Keep both abstraction levels:

- `FrancisTurbineAffinity`
- `FrancisTurbineOpenHPL`

The second should initially implement Euler power + guide-vane geometry + dominant loss terms. Add the remaining empirical corrections only after the basic model validates.

## 7. Generator / rotating mass

Current HydroPowerDynamics separates:

```text
Turbine → RotorInertia → SimpleGenerator
```

This is actually a useful architecture for grid studies.

OpenHPL's simple generator is more tightly formulated around mechanical/electrical energy balance and load power.

Recommended approach:

- retain separate MTK components,
- move bearing/friction terms into a clearly defined shaft/generator model,
- create an `OpenHPLSimpleGenerator` validation model,
- later add a synchronous-machine model rather than overloading `SimpleGenerator`.

## 8. Governor

The current MTK governors are PID/GGOV1-style.

OpenHPL's governor is structurally different and includes:

- power-to-guide-vane lookup,
- permanent droop,
- transient droop,
- pilot servo lag `T_p`,
- main servo integration `T_g`,
- transient droop block `T_r`,
- guide-vane opening/closing rate limits,
- guide-vane position limits.

For FREKI/FCR work this is a high-priority model because transient droop interacts directly with waterway oscillations.

Add a separate:

`OpenHPLGovernor`

rather than altering `PIDGovernor` or `GGOV1Governor`.

## 9. Proposed new source structure

```text
src/
  interfaces/
    hydraulic_port.jl
    rotational_port.jl

  hydraulic/
    reservoir_boundary.jl
    dynamic_reservoir.jl
    rigid_pipe.jl
    elastic_penstock_lumped.jl
    elastic_penstock_fv.jl
    surge_tank.jl
    draft_tube.jl
    gate.jl

  turbines/
    francis_affinity.jl
    francis_openhpl.jl
    pelton.jl

  mechanical/
    shaft.jl
    simple_generator.jl

  control/
    pid_governor.jl
    ggov1.jl
    openhpl_governor.jl
```

Do not perform this file split until equivalence tests exist; first implement the new models alongside the current source layout.

## 10. Refactor order

### Stage 1 — interfaces and baseline
- Freeze current examples/tests.
- Define common sign conventions.
- Establish units and nominal operating-point parameters.
- Add OpenHPL-equivalence tests.

### Stage 2 — waterway core
- Dynamic reservoir.
- Rigid pipe with elevation/gravity.
- Two-state surge tank.
- Lumped elastic penstock.

### Stage 3 — turbine
- Francis Euler-power model.
- Guide-vane geometry.
- Runner losses.

### Stage 4 — governor
- Reproduce OpenHPL permanent + transient droop.
- Servo lags.
- Opening/closing rate limits.
- Gate saturation.

### Stage 5 — distributed water hammer
- N-cell elastic penstock.
- Compare eigenmodes and pressure-wave propagation with OpenHPL PenstockKP.

### Stage 6 — integrated validation
Use the same operating point and disturbance in both implementations:

```text
reservoir
  → headrace/pipe
  → surge tank
  → penstock
  → Francis turbine
  → shaft/generator
  → grid/load
        ↑
     governor
```

Compare:

- flow,
- turbine inlet pressure/head,
- surge level,
- mechanical power,
- shaft speed/frequency,
- guide-vane position,
- dominant oscillation frequency,
- damping,
- steady-state error.

## 11. First implementation target

Start with the surge tank, not the detailed Francis turbine.

Reason: it gives a contained, physically meaningful two-state DAE upgrade and directly attacks the initialization/water-mass-oscillation problem relevant to the FREKI/FCR work.

Target states:

```text
x = [h, Vdot]
```

Core equations for the simple OpenHPL case:

```text
m = ρ A h
v = Vdot/A
M = m v

dm/dt = mdot_branch

dM/dt =
    mdot_branch·v
  + (p_bottom - p_atm)A
  - F_f(v,h)
  - m g

mdot_branch = ρ Vdot
```

After this model works, compare it directly against the current one-state `SurgeTank` under the same pressure/flow disturbance.

## 12. Validation principle

Every OpenHPL-derived MTK model should have:

1. equation-reference comments,
2. a minimal standalone test,
3. steady-state initialization test,
4. disturbance simulation,
5. comparison against the reduced MTK model,
6. comparison against OpenHPL output where available.

The goal of this branch is therefore:

> **OpenHPL physics + ModelingToolkit composability + control/FCR analysis tools.**
