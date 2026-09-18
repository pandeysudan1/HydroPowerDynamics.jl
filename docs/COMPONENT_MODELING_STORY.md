# Building Equation-Based Hydropower Components — A Modeling Story

This note describes the general method used in the `Refractor_OpenHPL` branch
to build reusable hydropower and power-system components with
ModelingToolkit.jl.

The central idea is simple:

> We do not start from code. We start from conservation laws, then add
> constitutive relations, then connection equations, then ask whether the
> assembled problem is structurally well posed.

The current examples are the reservoir, rigid pipe, surge tank, and a new
lookup-table turbine.

---

## 1. Source tree organization

The source tree is now organized by physical domain:

```text
src/
├── HydroPowerDynamics.jl
├── Core/
│   └── Utils.jl
├── Interfaces/
│   └── Connectors.jl
├── Waterways/
│   ├── LegacyHydraulic.jl
│   ├── RigidPipe.jl
│   ├── SurgeTank.jl
│   └── Reservoirs/
│       └── Reservoirs.jl
├── Turbines/
│   ├── LegacyTurbines.jl
│   └── TurbineLookup.jl
├── Electromechanical/
│   └── Mechanical.jl
└── Controls/
    ├── Governors.jl
    └── AGC.jl
```

This makes the package read like the plant:

```text
Interfaces
   ↓
Waterways
   ↓
Turbines
   ↓
Electromechanical
   ↓
Controls
```

Legacy components are retained while new OpenHPL-style models are developed
beside them.

---

# 2. The general problem-solving method

When we build a physical component, we follow the same sequence.

## Step 1 — choose the control volume

Ask:

- What enters the component?
- What leaves?
- What is stored?
- What physical domain is exchanged through each port?

For a rigid pipe, mass passes through but momentum is stored.

For a reservoir, mass is stored.

For a surge tank, both mass and momentum are stored.

For a turbine, hydraulic energy enters and mechanical energy leaves.

---

## Step 2 — write conservation equations first

Examples:

### Mass

```math
\frac{dm}{dt}
=
\sum \dot m.
```

### Linear momentum

```math
\frac{dM}{dt}
=
\sum F
+
\text{momentum flux}.
```

### Rotational momentum

```math
J\frac{d\omega}{dt}
=
\sum \tau.
```

### Energy conversion

For a turbine,

```math
P_h
=
\rho gQH
```

and

```math
P_m
=
\eta P_h.
```

These equations define the physics before numerical details are introduced.

---

## Step 3 — add constitutive and algebraic relations

Examples include:

```math
A=\frac{\pi D^2}{4},
```

```math
Q=\frac{\dot m}{\rho},
```

```math
v=\frac{Q}{A},
```

```math
\Delta p_f
=
f_D\frac{L}{D}\frac{\rho v|v|}{2}.
```

In a lookup turbine, the constitutive law is not an analytical formula for
efficiency. Instead, the constitutive relation comes from a hill chart.

---

## Step 4 — write connector equations

The connector equations tell us how the component participates in a network.

For a two-port hydraulic component:

```math
\dot m_a+\dot m_b=0.
```

For the refactored hydraulic interfaces:

```text
HydraulicPortZ:
    p   pressure
    z   elevation
    dm  mass flow
```

For a turbine shaft:

```text
RotationalPort:
    phi
    omega
    tau
```

The package convention is always:

> flow variables are positive into the component.

That means a turbine exporting mechanical power has negative shaft torque into
the turbine component.

---

# 3. Reservoir story

The reservoir starts from mass storage.

For a rectangular geometry,

```math
m=\rho LWh.
```

The mass balance is

```math
\frac{dm}{dt}=\dot m_{port}.
```

Then the hydrostatic relation is added:

```math
p=p_{atm}+\rho gh.
```

The final component therefore has:

```text
balance:
    dm/dt

algebra:
    geometry
    volume
    mass
    absolute level

connection:
    pressure
    elevation
    connector mass flow
```

Its step and ramp tests have exact analytical solutions, making the reservoir
an ideal first validation component.

---

# 4. Rigid-pipe story

The rigid pipe is a momentum-storage element.

The governing equation is

```math
\frac{L}{A}\frac{d\dot m}{dt}
=
p_a-p_b
+
\rho g(z_a-z_b)
-
\Delta p_f.
```

Then we define:

```math
Q=\frac{\dot m}{\rho},
```

```math
v=\frac{Q}{A},
```

```math
Re=\frac{|\dot m|D}{\mu A},
```

and the Darcy loss.

The Trollheim headrace benchmark provides:

```text
L = 4496.5 m
D = 6.3 m
Q = 37 m³/s
```

so the model can be checked against direct Darcy-Weisbach calculations.

---

# 5. Surge-tank story

The surge tank contains both mass and momentum storage.

Its states are

```text
h(t)
Vdot(t)
```

with

```math
m=\rho A\ell,
```

```math
M=mv.
```

Mass balance:

```math
\frac{dm}{dt}=\dot m_{port}.
```

Momentum balance:

```math
\frac{dM}{dt}
=
\dot m_{port}v+F_p-F_f-F_g.
```

This creates a physical restoring mechanism:

```text
junction pressure changes
        ↓
shaft water accelerates
        ↓
tank level changes
        ↓
gravity force changes
        ↓
flow reverses / oscillates
```

That is why a surge tank cannot be represented faithfully by only a static
pressure relation when studying hydraulic oscillations.

---

# 6. Why a lookup-table turbine?

A mechanistic turbine is valuable because it exposes the internal physics.

OpenHPL's main Francis model follows this philosophy and uses Euler turbine
relations, runner geometry, velocity triangles, guide-vane geometry, and loss
models.

However, real plant work often provides:

- measured hill charts,
- manufacturer efficiency maps,
- commissioning data,
- identified operating-point tables.

In those cases a lookup turbine can represent the known plant behavior more
directly.

Therefore the package should support two complementary turbine families:

```text
FrancisMechanistic
    OpenHPL / Euler physics
    useful for physical interpretation and extrapolation

TurbineLookup
    measured / manufacturer tables
    useful for plant-specific validation and FREKI workflows
```

---

# 7. TurbineLookup formulation

The new model uses the standard dimensionless turbine coordinates.

Unit speed:

```math
n_{11}
=
\frac{nD}{\sqrt{H}}.
```

Unit discharge:

```math
Q_{11}
=
\frac{Q}{D^2\sqrt{H}}.
```

The table coordinates are:

```text
x = gate opening
y = n11
```

and two surfaces are stored:

```text
Q11(gate,n11)
eta(gate,n11)
```

The component then reconstructs actual flow:

```math
Q
=
Q_{11}D^2\sqrt{H}.
```

Hydraulic power:

```math
P_h
=
\rho gQH.
```

Mechanical power:

```math
P_m
=
\eta P_h.
```

Mechanical torque:

```math
\tau_m
\approx
\frac{P_m}{\omega}
```

with a small regularization around zero speed.

---

# 8. Lookup table structure

The first implementation uses a compact 3x3 hill chart.

Gate axis:

```text
0.2   0.6   1.0
```

Unit-speed axis:

```text
40   50   60
```

The implementation performs clamped bilinear interpolation.

The default table is deliberately labeled as a demonstration table.

It is tuned so that the design point

```text
H      = 371 m
D      = 2.5 m
n11    = 50
gate   = 1.0
Q      = 37 m³/s
eta    = 0.97
```

matches the Trollheim benchmark already used in the repository.

This does **not** mean the complete table is measured Trollheim data.

When real hill-chart data becomes available, the table entries should be
replaced without changing the surrounding component equations.

---

# 9. Why use unit quantities?

Using (n_{11}) and (Q_{11}) separates the turbine map from one exact head.

A dimensional table such as

```text
gate -> Q
```

is valid only near one head.

A unit-quantity map allows the component to scale with changing head:

```text
changing H
   ↓
changing n11
   ↓
lookup Q11 and eta
   ↓
recover actual Q
```

This is much better suited to a hydropower transient model.

---

# 10. Equation and unknown consistency

A component connected to the outside world is intentionally an **open**
mathematical object.

For example, a turbine port exposes pressure, elevation, mass flow, speed,
torque, and gate input. Those external quantities are not all determined by
the turbine alone.

Therefore:

> Equation count should be checked on an assembled test system, not by requiring
> an isolated physical component with free ports to be square.

The package now provides:

```julia
model_structure_report(sys)
```

which returns:

```text
equations
unknowns
parameters
balanced
```

with

```text
balanced = equations == unknowns
```

The turbine test assembles:

```text
upstream hydraulic boundary
        ↓
TurbineLookup
        ↓
downstream hydraulic boundary

shaft speed boundary
        ↔
TurbineLookup shaft

gate source
        →
TurbineLookup gate
```

After `@mtkcompile`, the test requires

```julia
report.equations == report.unknowns
```

and prints the actual counts in the CI log.

The lookup turbine itself currently contains:

```text
29 scalar parameters
11 named internal algebraic variables
4 physical/signal connectors
```

The exact equation/unknown count of the assembled compiled model is obtained
from ModelingToolkit because connection expansion and structural
transformations change the raw count.

Compilation plus structural balance plus successful initialization is stronger
evidence than manually counting source lines.

---

# 11. Turbine validation

The design-point unit test chooses

```text
H = 371 m
gate = 1.0
n11 = 50
D = 2.5 m
```

The required rotational speed is

```math
n
=
\frac{n_{11}\sqrt{H}}{D}
\approx
385.23\;\mathrm{rpm}.
```

At this point the lookup table gives

```math
Q_{11}
=
\frac{37}{2.5^2\sqrt{371}},
```

so the reconstructed flow must be

```math
Q=37\;\mathrm{m^3/s}.
```

The efficiency table gives

```math
\eta=0.97.
```

Therefore

```math
P_m
=
1000\times9.81\times37\times371\times0.97
\approx
130.622\;\mathrm{MW}.
```

The unit test checks all of these quantities.

It also tests a point halfway between table rows to verify the bilinear
interpolation itself.

---

# 12. Validation hierarchy

For every component we should use the same sequence.

## Layer 1 — equation-level physics

Derive the result analytically.

## Layer 2 — component unit test

Use the smallest possible boundary system.

## Layer 3 — structural test

Check

```text
equations == unknowns
```

after assembling boundaries.

## Layer 4 — standard signal test

Use:

```text
step
ramp
sine
operating-point sweep
```

as appropriate.

## Layer 5 — Trollheim benchmark

Connect progressively:

```text
Reservoir
   ↓
RigidPipe
   ↓
SurgeTank
   ↓
Penstock
   ↓
TurbineLookup
   ↓
Shaft
   ↓
Generator
```

## Layer 6 — plant-data identification

For FREKI-style work:

```text
measurements
    ↓
identify / tune model
    ↓
validate model
    ↓
use model for FCR qualification studies
```

---

# 13. Current files

New organization:

```text
src/Interfaces/Connectors.jl
src/Waterways/Reservoirs/Reservoirs.jl
src/Waterways/RigidPipe.jl
src/Waterways/SurgeTank.jl
src/Turbines/TurbineLookup.jl
```

Turbine test:

```text
test/test_turbine_lookup.jl
```

Existing hydraulic validation:

```text
test/test_reservoir.jl
test/test_trollheim_reservoir_validation.jl
test/test_rigid_pipe.jl
test/test_openhpl_surge_tank.jl
test/test_standard_signals.jl
```

---

# 14. What comes next

The next useful integrated problem is no longer another isolated component.

It is the first complete nonlinear hydraulic chain:

```text
ReservoirBoundaryZ
        ↓
Trollheim RigidPipe
        ↓
surge junction
        ├── OpenHPLSurgeTank
        ↓
penstock
        ↓
TurbineLookup
```

For that assembled system we should:

1. derive the steady operating point,
2. count equations and unknowns,
3. compile with MTK,
4. initialize from physical states only,
5. apply gate step and ramp tests,
6. compare flow, pressure, surge level and mechanical power with analytical
   expectations,
7. then attach the shaft and generator.

That is the natural bridge from isolated component tests to a real
hydropower-system model.
