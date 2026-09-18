# Modeling Notes — OpenHPL-Style Notation, Units, Balances, and Validation

This note defines the modeling discipline for the `Refractor_OpenHPL` branch.

The aim is not only to translate OpenHPL equations into ModelingToolkit syntax.
The aim is to make every hydropower component readable as a physical story:

```text
physical conservation law
        ↓
constitutive / algebraic relations
        ↓
connector equations
        ↓
initialization
        ↓
standard input signal
        ↓
analytical expectation
        ↓
MTK simulation
        ↓
unit-test comparison
        ↓
integrated Trollheim validation
```

The first three hydraulic components following this pattern are:

- `ReservoirBoundaryZ` / `DynamicReservoir`
- `RigidPipe`
- `OpenHPLSurgeTank`

The same pattern should later be used for the penstock, Francis turbine,
shaft, synchronous generator, governor, transformer, line, and grid.

---

## 1. Notation philosophy

The notation should stay close to OpenHPL and standard fluid-mechanics practice.

### Primary hydraulic symbols

| Symbol | Code name | Meaning | SI unit |
|---|---|---|---|
| (p) | `p` | pressure | Pa |
| (z) | `z` | absolute elevation | m |
| (h) | `h` | local water height / head | m |
| (Q) | `Q`, `Vdot` | volume flow rate | m³/s |
| (dot m) | `dm`, `mdot` | mass flow rate | kg/s |
| (ho) | `rho` | water density | kg/m³ |
| (v) | `v` | mean fluid velocity | m/s |
| (A) | `A` | cross-sectional area | m² |
| (L) | `L` | length | m |
| (D) | `D_pipe`, `D_tank` | diameter | m |
| (g) | `g` | gravity | m/s² |
| (f_D) | `f_D` | Darcy friction factor | – |
| (Re) | `Re` | Reynolds number | – |
| (Delta p_f) | `dp_f` | friction pressure loss | Pa |
| (h_f) | `h_f` | friction head loss | m |

### Mechanical symbols

| Symbol | Code name | Meaning | SI unit |
|---|---|---|---|
| (phi) | `phi` | angular position | rad |
| (omega) | `omega` | angular speed | rad/s |
| (	au) | `tau` | torque | N·m |
| (J) | `J` | rotational inertia | kg·m² |

### Electrical symbols

| Symbol | Code name | Meaning | SI unit |
|---|---|---|---|
| (v_r,v_i) | `vr, vi` | rectangular voltage components | V |
| (i_r,i_i) | `ir, ii` | rectangular current components | A |
| (P) | `P` | active power | W |
| (Q_e) | `Q` or `Q_e` | reactive power | var |

The hydraulic notation deliberately distinguishes:

```text
Q      volumetric flow [m³/s]
dm     mass flow       [kg/s]
```

with

```math
\dot m = \rho Q.
```

This distinction is essential because the physical connectors use mass flow,
while many hydropower reference equations are written in volumetric flow.

---

## 2. Connector convention

### HydraulicPortZ

The preferred OpenHPL-style hydraulic connector is

```julia
@connector HydraulicPortZ begin
    p(t)
    z(t)
    dm(t), [connect = Flow]
end
```

The variables have the following roles:

```text
across variables : p, z
through variable : dm
```

At a hydraulic connection set, ModelingToolkit imposes:

```math
p_1=p_2=\ldots
```

```math
z_1=z_2=\ldots
```

and

```math
\sum_k \dot m_k = 0.
```

The package sign convention is:

> positive connector mass flow is positive into a component.

This one rule must be used in every hydraulic component.

---

## 3. Required equation order for every component

Each new component should be readable in the same physical order.

### 3.1 Balance laws first

The first equations should be conservation equations.

Typical hydraulic mass balance:

```math
\frac{dm}{dt}=\sum \dot m.
```

Typical hydraulic momentum balance:

```math
\frac{dM}{dt}=\sum F + \text{momentum flux}.
```

Typical rotational balance:

```math
J\frac{d\omega}{dt}=\sum\tau.
```

Typical electrical dynamic balance would similarly follow machine flux or
energy equations before algebraic network equations.

### 3.2 Algebraic and constitutive relations second

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
Re=\frac{|\dot m|D}{\mu A},
```

and

```math
\Delta p_f
=
f_D\frac{L}{D}\frac{\rho}{2}v|v|.
```

### 3.3 Connection equations third

Examples:

```julia
port_a.dm + port_b.dm ~ 0
port_a.dm ~ dm
```

or for a boundary:

```julia
port.p ~ p_atm + rho*g*h
port.z ~ z_out
```

This keeps topology separate from internal physics.

### 3.4 Initialization last

Only independent physical states should be initialized directly.

Avoid giving guesses or fixed initial values to every algebraic variable.

For example, initialize

```text
reservoir : h
pipe      : dm
surge tank: h, Vdot
```

and derive area, mass, Reynolds number, pressure, forces, etc. algebraically.

This reduces the risk of overdetermined MTK initialization systems.

---

# Part I — Reservoir

## 4. Reservoir as a conservation problem

A reservoir is first a storage element.

For the dynamic model,

```math
m=\rho V.
```

The governing balance is

```math
\frac{dm}{dt}=\dot m_{port}.
```

Because positive flow is defined into the component, normal turbine discharge
appears as

```math
\dot m_{port}<0.
```

The water level therefore decreases naturally.

For vertical walls,

```math
V=LWh
```

and therefore

```math
\rho LW\frac{dh}{dt}=\dot m_{port}.
```

The hydrostatic algebraic relation is

```math
p=p_{atm}+\rho gh.
```

The elevation relation is

```math
h_{abs}=z_{out}+h.
```

---

## 5. Reservoir step input

Consider a step in outflow at (t=t_0):

```math
Q_{out}(t)
=
\begin{cases}
0,&t<t_0\\
Q_s,&t\ge t_0.
\end{cases}
```

For the rectangular reservoir,

```math
\frac{dh}{dt}
=
-\frac{Q_s}{LW}.
```

Hence

```math
h(t)
=
h_0
-
\frac{Q_s}{LW}(t-t_0)
```

for (t\ge t_0).

This is an exact analytical solution.

The unit test applies this signal to `DynamicReservoir` and verifies the MTK
result against the closed-form value.

---

## 6. Reservoir ramp input

Now let the withdrawal ramp after (t_0):

```math
Q_{out}(t)=k_r(t-t_0).
```

Then

```math
\frac{dh}{dt}
=
-\frac{k_r(t-t_0)}{LW}.
```

Integrating gives

```math
h(t)
=
h_0
-
\frac{k_r}{2LW}(t-t_0)^2.
```

This is also exact.

So the reservoir gives us an ideal first example of the workflow:

```text
step input → linear level response
ramp input → quadratic level response
```

and both can be tested without relying on a numerical reference simulation.

---

# Part II — Rigid Pipe

## 7. Pipe as a momentum problem

The rigid pipe does not store mass in the incompressible approximation.

Its continuity equation is

```math
\dot m_a+\dot m_b=0.
```

Its dynamic state comes from water-column momentum.

The model uses

```math
\frac{L}{A}\frac{d\dot m}{dt}
=
p_a-p_b
+
\rho g(z_a-z_b)
-
\Delta p_f.
```

This equation tells a direct physical story:

```text
pressure force
+ gravity
- friction
= water-column acceleration
```

The constitutive relations are

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
Re=\frac{|\dot m|D}{\mu A},
```

and Darcy-Weisbach friction

```math
\Delta p_f
=
f_D\frac{L}{D}\frac{\rho}{2}v|v|.
```

---

## 8. Trollheim steady validation

The detailed Trollheim headrace example uses approximately

```text
L = 4496.5 m
D = 6.3 m
Q = 37 m³/s
```

giving approximately

```text
A      = 31.17245 m²
v      = 1.18695 m/s
Re     = 7.48e6
f_D    = 0.00870
dp_f   = 4.37 kPa
h_f    = 0.446 m
```

The steady unit test applies

```math
p_a-p_b=\Delta p_f
```

at equal elevation.

The momentum equation then predicts

```math
\frac{d\dot m}{dt}=0.
```

The simulated flow should remain at the analytical operating point.

---

## 9. Pipe pressure-step input

At zero initial flow,

```math
\Delta p_f=0.
```

Applying an instantaneous pressure difference (Delta p) gives

```math
\left.\frac{d\dot m}{dt}\right|_{0^+}
=
\frac{A}{L}\Delta p.
```

For the Trollheim headrace and a 10 kPa pressure step,

```math
\left.\frac{d\dot m}{dt}\right|_{0^+}
\approx 69.3\;\mathrm{kg/s^2}.
```

The unit test estimates the initial numerical slope over a very short time and
compares it with this analytical result.

---

## 10. Pipe pressure-ramp input

Suppose

```math
\Delta p(t)=rt
```

and the pipe starts from rest.

At (t=0), friction is initially zero and

```math
\frac{d^2\dot m}{dt^2}
=
\frac{A}{L}r.
```

Therefore the short-time analytical response is

```math
\dot m(t)
\approx
\frac12\frac{A}{L}rt^2.
```

The ramp test verifies this early-time quadratic response.

It is intentionally a local analytical validation because Darcy friction makes
the full finite-time ramp response nonlinear.

---

# Part III — Surge Tank

## 11. Surge tank as mass plus momentum storage

The OpenHPL-style surge tank stores both water mass and water-column momentum.

The states are

```text
h(t)      water height
Vdot(t)   shaft volume flow
```

with

```math
A=\frac{\pi D^2}{4},
```

```math
v=\frac{\dot V}{A},
```

```math
m=\rho A\ell,
```

and

```math
M=mv.
```

The shaft geometry is

```math
\cos\theta=\frac{H}{L},
```

```math
\ell=\frac{h}{\cos\theta}.
```

### Mass balance

```math
\frac{dm}{dt}=\dot m_{port}.
```

### Momentum balance

```math
\frac{dM}{dt}
=
\dot m_{port}v
+
F_p
-
F_f
-
F_g.
```

where

```math
F_p=(p_{bottom}-p_{atm})A,
```

```math
F_g=mg\cos\theta,
```

and

```math
F_f
=
f\frac{\ell}{D}
\frac{\rho A}{2}v|v|.
```

This is the main physical difference between the reduced one-state surge tank
and the OpenHPL-style model: the shaft water does not respond infinitely fast.

---

## 12. Surge-tank hydrostatic equilibrium

For a vertical shaft,

```math
H=L,
```

so

```math
\cos\theta=1.
```

At

```math
\dot V=0
```

and steady height (h_0), the force balance becomes

```math
(p_b-p_{atm})A
=
\rho Ah_0g.
```

Hence

```math
p_b
=
p_{atm}+\rho gh_0.
```

This is the first surge-tank unit test.

---

## 13. Surge-tank pressure-step input

Now apply a small pressure step

```math
\Delta p
```

above the hydrostatic equilibrium pressure.

At the instant after the step,

```math
v=0,
```

so friction is zero.

For a vertical shaft the initial volume-flow acceleration reduces to

```math
\left.\frac{d\dot V}{dt}\right|_{0^+}
=
\frac{\Delta p}{\rho\ell_0}.
```

For example, with

```text
h0 = 50 m
Delta p = 1000 Pa
rho = 1000 kg/m³
```

the initial acceleration is

```math
\frac{d\dot V}{dt}
=
\frac{1000}{1000\times50}
=
0.02\;\mathrm{m^3/s^2}.
```

The MTK test estimates the initial slope and compares it with this value.

A later test can add a pressure ramp and compare the short-time second
derivative in the same way.

---

# Part IV — Standard validation signals

## 14. Why use step and ramp signals?

A step is useful because it exposes:

- sign conventions,
- static gain,
- inertia,
- discontinuity handling,
- initial acceleration.

A ramp is useful because it exposes:

- integration behavior,
- accumulation,
- second-order local response,
- whether the state equation is physically consistent.

For the first hydraulic components:

| Component | Step input | Expected response | Ramp input | Expected response |
|---|---|---|---|---|
| DynamicReservoir | outflow | linear (h(t)) | outflow | quadratic (h(t)) |
| RigidPipe | pressure | initial linear (dot m) | pressure | initial quadratic (dot m) |
| OpenHPLSurgeTank | pressure | initial linear (dot V) | pressure | later extension |

The important point is that the expected result is derived before running the
simulation.

---

# Part V — Unit-test structure

## 15. Required test hierarchy

Every new physical component should have four layers of tests.

### A. Dimensional / algebraic test

Check geometry and static equations.

Examples:

```text
A = pi*D²/4
Q = dm/rho
P = rho*g*Q*H
```

### B. Conservation test

Check mass, momentum, torque, or current balance.

Examples:

```text
dm_a + dm_b = 0
tau_a + tau_b = J*domega/dt
ir_1 + ir_2 + ... = 0
```

### C. Standard-signal test

Apply step/ramp input and compare with analytical behavior.

### D. Integrated plant test

Connect the component into the Trollheim chain and compare operating-point and
transient quantities.

---

## 16. Current test files

The relevant tests are:

```text
test/test_reservoir.jl
test/test_trollheim_reservoir_validation.jl
test/test_rigid_pipe.jl
test/test_openhpl_surge_tank.jl
test/test_standard_signals.jl
```

The new standard-signal test contains:

```text
DynamicReservoir:
  step outflow
  ramp outflow

RigidPipe:
  pressure step
  pressure ramp

OpenHPLSurgeTank:
  hydrostatic equilibrium
  pressure-step initial acceleration
```

All tests are included from

```text
test/runtests.jl
```

and therefore belong to the MTK-11 CI gate.

---

# Part VI — A component template for future refactors

## 17. Recommended source pattern

Every new component should approximately follow:

```julia
@mtkmodel ComponentName begin

    @parameters begin
        # geometry
        # material/fluid constants
        # rated/reference values
    end

    @variables begin
        # dynamic states first
        # algebraic variables second
    end

    @components begin
        # physical ports
        # causal signal ports
    end

    @equations begin

        # ------------------------------------------------------------
        # 1. MASS / MOMENTUM / ENERGY BALANCE
        # ------------------------------------------------------------
        # D(state) ~ ...

        # ------------------------------------------------------------
        # 2. ALGEBRAIC / CONSTITUTIVE RELATIONS
        # ------------------------------------------------------------
        # geometry
        # friction
        # efficiency
        # pressure / torque / power relations

        # ------------------------------------------------------------
        # 3. CONNECTION EQUATIONS
        # ------------------------------------------------------------
        # port flow assignments
        # topology relations
    end
end
```

The exact equation ordering does not change the mathematics, but it makes code
review much easier: a reader first sees what is conserved, then how the
component behaves internally, then how it connects to the network.

---

# Part VII — Trollheim story

## 18. Building the plant one conservation law at a time

The refactor is now developing naturally from upstream to downstream.

### Chapter 1 — reservoir

The reservoir stores mass.

```text
water withdrawal
      ↓
mass decreases
      ↓
water level decreases
      ↓
hydrostatic pressure changes
```

### Chapter 2 — rigid headrace

The headrace stores momentum.

```text
pressure + elevation difference
            ↓
     accelerates water
            ↓
      Darcy friction
            ↓
     opposes the motion
```

### Chapter 3 — surge tank

The surge tank stores both mass and momentum.

```text
junction pressure
      ↓
shaft water accelerates
      ↓
water enters/leaves tank
      ↓
tank level changes
      ↓
hydrostatic restoring force changes
```

These three components already produce the first nonlinear hydraulic subsystem:

```text
ReservoirBoundaryZ
        │
        ▼
    RigidPipe
        │
        ├──── OpenHPLSurgeTank
        │
        ▼
 next Penstock
```

The next step is to make this exact three-component subsystem executable with
the same Trollheim operating point and validate its steady state before
attaching the turbine.

---

## 19. Definition of done for each component

A refactored component is not considered validated merely because it compiles.

It should satisfy:

```text
[ ] notation and SI units documented
[ ] balance equations visible in source
[ ] algebraic relations documented
[ ] connector equations documented
[ ] standalone compile test
[ ] equilibrium test
[ ] step or ramp response test
[ ] analytical comparison
[ ] Trollheim parameter comparison where applicable
[ ] CI pass on current MTK branch
```

This checklist should become the common standard for the remainder of the
OpenHPL refactor.
