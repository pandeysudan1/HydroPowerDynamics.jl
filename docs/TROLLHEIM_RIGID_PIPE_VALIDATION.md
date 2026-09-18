# Trollheim Rigid Pipe Validation

This document validates the refactored `RigidPipe` component against
closed-form hydraulic relations and the Trollheim headrace geometry already
used in the repository examples.

## 1. Trollheim headrace benchmark

The detailed Trollheim example uses:

| Quantity | Value |
|---|---:|
| Headrace length | 4496.5 m |
| Headrace diameter | 6.3 m |
| Rated discharge | 37 m³/s |
| Water density | 1000 kg/m³ |
| Dynamic viscosity | 1.0e-3 Pa·s |
| Roughness | 1.5e-5 m |
| Gravity | 9.81 m/s² |

The refactored component is intended to reproduce the corresponding analytical
velocity, Reynolds number, Darcy friction factor, pressure loss, and transient
water-column acceleration.

## 2. Model

`RigidPipe` uses two elevation-aware hydraulic connectors:

```julia
port_a = HydraulicPortZ()
port_b = HydraulicPortZ()
```

The package convention is positive mass flow from `port_a` to `port_b`.

### Continuity

```math
\dot m_a + \dot m_b = 0
```

and

```math
\dot m_a = \dot m.
```

### Kinematics

```math
A = \frac{\pi D^2}{4}
```

```math
Q = \frac{\dot m}{\rho}
```

```math
v = \frac{Q}{A}.
```

### Reynolds number

```math
Re = \frac{|\dot m|D}{\mu A}.
```

### Darcy-Weisbach loss

```math
\Delta p_f
=
f_D\frac{L}{D}\frac{\rho}{2}v|v|.
```

The corresponding head loss is

```math
h_f = \frac{\Delta p_f}{\rho g}.
```

### Momentum balance

```math
\frac{L}{A}\frac{d\dot m}{dt}
=
p_a-p_b
+
\rho g(z_a-z_b)
-
\Delta p_f.
```

This makes the pressure and elevation contributions explicit and avoids hiding
gravity inside an arbitrary head parameter.

## 3. Analytical Trollheim values

For

```text
L = 4496.5 m
D = 6.3 m
Q = 37 m³/s
```

the area is

```math
A = \frac{\pi 6.3^2}{4}
  = 31.1724531\;\mathrm{m^2}.
```

The mean velocity is

```math
v = \frac{37}{31.1724531}
  = 1.18694541\;\mathrm{m/s}.
```

The Reynolds number is

```math
Re \approx 7.47776\times10^6.
```

Using the package Darcy-factor relation gives

```math
f_D \approx 0.00869835.
```

The corresponding pressure loss is

```math
\Delta p_f
\approx
4373.23\;\mathrm{Pa}.
```

and the head loss is

```math
h_f
=
\frac{4373.23}{1000\times9.81}
\approx
0.445793\;\mathrm{m}.
```

So at rated Trollheim flow, the long 6.3 m headrace still has a relatively
small steady head loss because of its large cross-sectional area.

## 4. Steady-state MTK unit test

The unit test applies exactly the analytical friction pressure drop across the
pipe:

```text
upstream pressure - downstream pressure = 4373.23 Pa
```

with equal elevations.

The initial mass flow is

```math
\dot m
=
\rho Q
=
37000\;\mathrm{kg/s}.
```

At this point,

```math
\frac{d\dot m}{dt}=0
```

because

```math
p_a-p_b=\Delta p_f.
```

The MTK solution is checked against:

| Variable | Analytical target |
|---|---:|
| Q | 37 m³/s |
| Re | 7.477756e6 |
| f_D | 0.00869835 |
| pressure loss | 4373.232 Pa |
| head loss | 0.445793 m |

This validates the steady constitutive physics of the component.

## 5. Transient acceleration test

A second unit test checks the water-column inertia directly.

At zero initial flow,

```math
v=0
```

so initially

```math
\Delta p_f=0.
```

Applying a pressure step

```math
\Delta p = 10000\;\mathrm{Pa}
```

gives the analytical initial mass-flow acceleration

```math
\frac{d\dot m}{dt}
=
\frac{A}{L}\Delta p.
```

For the Trollheim headrace,

```math
\frac{d\dot m}{dt}
=
\frac{31.1724531}{4496.5}\times10000
\approx
69.32\;\mathrm{kg/s^2}.
```

The MTK simulation is integrated over a very short interval and the numerical
slope is compared against this analytical value.

This checks the dynamic part of the model independently from the friction test.

## 6. Role in the refactored plant

The component now gives the first physically meaningful connection after the
reservoir:

```text
ReservoirBoundaryZ
       │
       ▼
    RigidPipe
       │
       ▼
   Surge Tank
```

For Trollheim, the intended use is:

```text
Reservoir
   ↓
4496.5 m / 6.3 m headrace
   ↓
surge tank
```

This model is deliberately incompressible. It captures rigid-column inertia
and Darcy friction but not elastic water hammer.

## 7. Model hierarchy

The hydraulic transport hierarchy should therefore become:

```text
RigidPipe
    ↓
incompressible / control-oriented

ElasticPenstockLumped
    ↓
single-mode compressibility

ElasticPenstockFV
    ↓
distributed water-hammer / OpenHPL PenstockKP equivalent
```

This lets FCR studies use the simplest model that retains the relevant
frequency range.

## 8. Files

Implementation:

```text
src/rigid_pipe.jl
```

Analytical unit test:

```text
test/test_rigid_pipe.jl
```

Related Trollheim references:

```text
quick_examples/sim_test/iso_surgetank_comparison.jl
quick_examples/trollheim_agc_surgetank_compare/
```

## 9. Next component

The next refactor target should be the surge-tank junction using
`HydraulicPortZ`.

The key analytical relations will be:

```math
A_s\frac{dZ}{dt}
=
Q_{headrace}-Q_{penstock}
```

and, for the momentum-retaining OpenHPL-style shaft,

```math
\frac{dM_s}{dt}
=
F_p-F_g-F_f+\text{momentum flux}.
```

The existing `OpenHPLSurgeTank` should therefore be upgraded from its current
single branch `HydraulicPort` interface to an elevation-aware junction
compatible with the new reservoir and `RigidPipe`.
