# Trollheim Reservoir Unit Validation

This document validates the refactored reservoir components against the
Trollheim benchmark data already used in `quick_examples/` and against
closed-form hydraulic relations.

The purpose is not merely to confirm that the MTK model integrates. The unit
tests compare simulated variables against values that can be calculated
analytically from conservation laws and hydrostatics.

> The Trollheim values in this repository are benchmark/reference values used
> by the existing examples. They should not be interpreted here as newly
> field-verified measurements.

## 1. Trollheim reference values

The detailed Trollheim example in

```text
quick_examples/sim_test/iso_surgetank_comparison.jl
```

uses:

| Quantity | Value |
|---|---:|
| Nominal head | 371 m |
| Rated discharge | 37 m³/s |
| Water density | 1000 kg/m³ |
| Gravity | 9.81 m/s² |
| Atmospheric pressure | 101325 Pa |
| Peak turbine efficiency | 0.97 |
| Headrace total length | 4496.5 m |
| Headrace diameter | 6.3 m |
| Surge-tank diameter | 3.4 m |
| Penstock sections | 363 m × 4.7 m; 145 m × 3.3 m |

These values form the hydraulic reference for the reservoir unit test.

## 2. Component under test

Two reservoir models are checked:

- `ReservoirBoundaryZ`
- `DynamicReservoir`

Both use the elevation-aware connector

```julia
HydraulicPortZ(p, z, dm)
```

with positive `dm` defined as flow into a component.

---

## 3. Test A — hydrostatic reservoir pressure

For a fixed reservoir with water head

```math
H = 371\;\mathrm{m},
```

the absolute outlet pressure follows directly from hydrostatics:

```math
p = p_{atm} + \rho g H.
```

Substituting the Trollheim benchmark values:

```math
p
=
101325
+
1000\times 9.81\times371
=
3\,740\,835\;\mathrm{Pa}.
```

Therefore:

```text
Gauge pressure     = 3.639510 MPa
Absolute pressure  = 3.740835 MPa
```

The unit test connects `ReservoirBoundaryZ` to a zero-flow probe and checks

```math
p_{MTK} \approx 3.740835\;\mathrm{MPa}.
```

It also verifies the elevation relation

```math
h_{abs}=z_{out}+H.
```

For `z_out=0`:

```math
h_{abs}=371\;\mathrm{m}.
```

This validates the two most important reservoir boundary quantities passed to
the next hydraulic component:

```text
pressure + absolute elevation
```

---

## 4. Test B — dynamic reservoir against an exact analytical solution

The dynamic reservoir geometry is

```math
A(h)=h(W+h\tan\alpha).
```

For vertical walls,

```math
\alpha=0,
```

so

```math
A(h)=hW.
```

With reservoir length `L`, stored volume and mass become

```math
V=LWh,
```

```math
m=\rho LWh.
```

The component mass balance is

```math
\frac{dm}{dt}=\dot m_{port}.
```

Under constant water withdrawal `\dot m_{out}>0`,

```math
\dot m_{port}=-\dot m_{out}.
```

Therefore

```math
\rho LW\frac{dh}{dt}=-\dot m_{out},
```

giving the exact solution

```math
h(t)
=
h_0
-
\frac{\dot m_{out}}{\rho LW}t.
```

### Numerical benchmark used in the unit test

The test chooses

```text
h0      = 371 m
L       = 500 m
W       = 100 m
rho     = 1000 kg/m³
dm_out  = 1000 kg/s
t       = 10 s
```

Hence

```math
\Delta h
=
\frac{1000}{1000\times500\times100}\times10
=
0.0002\;\mathrm{m}.
```

The exact final level is therefore

```math
h(10)
=
370.9998\;\mathrm{m}.
```

and the exact absolute pressure at that level is

```math
p(10)
=
p_{atm}+\rho g h(10)
=
3\,740\,833.038\;\mathrm{Pa}.
```

The MTK simulation is tested directly against these values.

This is stronger than a generic solver-success test because it verifies that
the component obeys conservation of mass quantitatively.

---

## 5. Test C — hydraulic power consistency

The Trollheim example also provides a convenient energy check.

The ideal hydraulic power at

```math
Q=37\;\mathrm{m^3/s}
```

and

```math
H=371\;\mathrm{m}
```

is

```math
P_h=\rho gQH.
```

Using the benchmark values:

```math
P_h
=
1000\times9.81\times37\times371
=
134.661\;\mathrm{MW}.
```

With peak turbine efficiency

```math
\eta=0.97,
```

the corresponding mechanical power is

```math
P_m
=
\eta P_h
=
130.622\;\mathrm{MW}.
```

This does not test turbine dynamics yet, but it creates an energy reference
that the later Francis-turbine refactor should reproduce at the same operating
point.

---

## 6. Relation to the existing Trollheim simulations

The compact AGC example uses:

```text
Gross head = 371 m
Rated flow = 37 m³/s
```

and an initial operating flow near

```text
Q0 = 22.1375 m³/s.
```

The more detailed surge-tank example additionally uses the Trollheim waterway
geometry:

```text
Reservoir
   ↓
4496.5 m headrace tunnel
   ↓
3.4 m surge tank
   ↓
two-section penstock
   ↓
Francis turbine
```

The present reservoir validation establishes the first boundary condition for
that full chain.

---

## 7. Validation hierarchy

The refactor should validate every component at three levels.

### Level 1 — analytical unit physics

Examples:

```text
reservoir hydrostatics
mass conservation
Darcy pressure loss
water-column inertia
shaft energy balance
AC power relation
```

### Level 2 — isolated MTK component simulation

Each component is connected to the smallest possible set of boundaries and
tested against the Level-1 result.

### Level 3 — Trollheim integrated benchmark

Components are connected into

```text
Reservoir
   ↓
Headrace
   ↓
Surge Tank
   ↓
Penstock
   ↓
Francis
   ↓
Shaft
   ↓
Generator
   ↓
Grid / Load
```

and compared with the existing Trollheim examples.

This prevents a numerically successful integrated simulation from hiding a
wrong sign convention or incorrect physical balance in an individual
component.

---

## 8. Unit-test implementation

The tests are located in

```text
test/test_trollheim_reservoir_validation.jl
```

and check:

1. fixed reservoir hydrostatic pressure,
2. fixed reservoir absolute elevation,
3. dynamic reservoir level under constant withdrawal,
4. dynamic pressure against the analytical solution,
5. dynamic volumetric outflow,
6. Trollheim hydraulic and mechanical power reference values.

They run as part of

```text
test/runtests.jl
```

and therefore form part of the MTK-11 CI gate.

---

## 9. Analytical benchmark summary

| Quantity | Analytical value |
|---|---:|
| Trollheim nominal head | 371 m |
| Reservoir gauge pressure | 3.639510 MPa |
| Reservoir absolute pressure | 3.740835 MPa |
| Dynamic level after 10 s test | 370.9998 m |
| Dynamic absolute pressure after 10 s | 3.740833038 MPa |
| Rated ideal hydraulic power | 134.661 MW |
| Rated mechanical power at η = 0.97 | 130.622014 MW |

The MTK unit tests should reproduce these values within their specified
numerical tolerances.

---

## 10. Next validation target

The next component should be the rigid headrace pipe.

Its primary analytical unit relation will be

```math
\frac{L}{A}\frac{d\dot m}{dt}
=
p_a-p_b
+
\rho g(z_a-z_b)
-
\Delta p_f,
```

with

```math
\Delta p_f
=
f\frac{L}{D}\frac{\rho v|v|}{2}.
```

The Trollheim headrace geometry already provides a direct benchmark:

```text
L = 4496.5 m
D = 6.3 m
```

so the next unit test can compare MTK pressure loss and water-column
acceleration against closed-form Darcy-Weisbach and momentum calculations.
