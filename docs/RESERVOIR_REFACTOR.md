# Reservoir Refactor — OpenHPL to ModelingToolkit 11

This note documents the first hydraulic component refactor in the
`Refractor_OpenHPL` branch.

The goal is to preserve the useful physics from OpenHPL while rebuilding the
components around the new HydroPowerDynamics.jl physical connector layer and
ModelingToolkit 11.

## 1. Role in the hydropower system

The reservoir is the upstream hydraulic boundary of the plant.

```text
Reservoir
    │
    │ HydraulicPortZ
    │ p, z, mdot
    ▼
Headrace / Pipe
    ▼
Surge Tank
    ▼
Penstock
    ▼
Turbine
```

The reservoir provides the initial hydraulic energy through water level and
elevation.

## 2. Connector

The refactored reservoir models use `HydraulicPortZ`:

```julia
@connector HydraulicPortZ begin
    p(t)
    z(t)
    dm(t), [connect = Flow]
end
```

with:

- `p`: pressure [Pa]
- `z`: absolute elevation [m]
- `dm`: mass flow rate [kg/s]

The package convention is:

> Positive `dm` means flow into the component.

Therefore, during normal generation, water leaves the reservoir and the
reservoir connector sees

```text
dm < 0
```

This convention is important because all later mass and momentum balances are
built around it.

## 3. Refactor strategy

Instead of one reservoir model with many Boolean mode switches, the refactor
starts with two explicit components:

1. `ReservoirBoundaryZ`
2. `DynamicReservoir`

This keeps each model structurally simple and makes MTK initialization easier
to reason about.

The legacy `Reservoir` component remains available during the refactor for
backward compatibility.

---

## 4. ReservoirBoundaryZ

`ReservoirBoundaryZ` represents a large upstream reservoir whose level is
effectively constant over the timescale of interest.

This is the preferred model for FCR and other short transient studies.

### Parameters

- `h`: water level above the outlet [m]
- `z_out`: absolute outlet elevation [m]
- `rho`: water density [kg/m³]
- `g`: gravitational acceleration [m/s²]
- `p_atm`: atmospheric pressure [Pa]

### Equations

The outlet pressure is hydrostatic:

```math
p_{out} = p_{atm} + \rho g h
```

The connector elevation is

```math
z = z_{out}
```

The absolute free-surface elevation is

```math
h_{abs} = z_{out} + h
```

### MTK implementation

```julia
@mtkmodel ReservoirBoundaryZ begin
    @parameters begin
        h = 100.0
        z_out = 0.0
        rho = 1000.0
        g = 9.81
        p_atm = 101_325.0
    end

    @variables begin
        h_abs(t)
    end

    @components begin
        port = HydraulicPortZ()
    end

    @equations begin
        port.z ~ z_out
        h_abs ~ z_out + h
        port.p ~ p_atm + rho * g * h
    end
end
```

---

## 5. DynamicReservoir

`DynamicReservoir` represents a reservoir whose water level changes as water
is withdrawn or added.

It follows the basic storage geometry used in
`OpenHPL/Waterway/Reservoir.mo`.

### Geometry

OpenHPL defines the vertical cross-sectional area as

```math
A(h) = h\left(W + h\tan\alpha\right)
```

where:

- `W` is the bed width,
- `alpha` is the side-wall angle.

The stored volume is

```math
V = A(h)L
```

and the stored mass is

```math
m = \rho V.
```

### Mass balance

Using the HydroPowerDynamics connector sign convention,

```math
\frac{dm}{dt} = \dot m_{port}.
```

For normal generation discharge,

```math
\dot m_{port} < 0,
```

so reservoir mass and level decrease.

### Hydrostatic outlet pressure

```math
p_{out} = p_{atm} + \rho g h.
```

### Absolute elevation

```math
h_{abs} = z_{out} + h.
```

### Volumetric outflow

A positive generation outflow is defined as

```math
Q_{out} = -\frac{\dot m_{port}}{\rho}.
```

### MTK implementation

```julia
@mtkmodel DynamicReservoir begin
    @parameters begin
        h_0 = 50.0
        z_out = 0.0
        L = 500.0
        W = 100.0
        alpha = 0.0
        rho = 1000.0
        g = 9.81
        p_atm = 101_325.0
    end

    @variables begin
        h(t) = h_0
        A(t)
        V(t)
        m(t)
        h_abs(t)
        Q_out(t)
    end

    @components begin
        port = HydraulicPortZ()
    end

    @equations begin
        A ~ h * (W + h * tan(alpha))
        V ~ A * L
        m ~ rho * V

        port.z ~ z_out
        h_abs ~ z_out + h
        port.p ~ p_atm + rho * g * h

        Q_out ~ -port.dm / rho

        D(m) ~ port.dm
    end
end
```

---

## 6. Mapping to OpenHPL

The main OpenHPL reference is:

```text
OpenHPL/Waterway/Reservoir.mo
```

The current refactor maps the following OpenHPL ideas:

| OpenHPL concept | HydroPowerDynamics refactor |
|---|---|
| Constant reservoir level | `ReservoirBoundaryZ` |
| Variable level reservoir | `DynamicReservoir` |
| Outlet pressure | Hydrostatic equation |
| Outlet elevation | `HydraulicPortZ.z` |
| Stored mass | `m = rho*V` |
| Storage dynamics | `D(m) ~ port.dm` |
| OpenHPL default no-inflow mode | `DynamicReservoir` |

OpenHPL also contains a more detailed inflow/momentum mode with reservoir
friction. That is intentionally deferred until the basic hydraulic chain is
validated.

---

## 7. Why keep two models?

For frequency-control studies, reservoir level changes are usually negligible
over the short simulation horizon. In that case, the fixed-level model avoids
adding an unnecessary slow state.

For longer or strongly nonlinear simulations, the storage state becomes
important.

This gives two useful model fidelities:

```text
ReservoirBoundaryZ
    ↓
fast FCR / governor / small-signal studies

DynamicReservoir
    ↓
longer nonlinear / storage / operating-point studies
```

This follows the broader package philosophy:

> Keep reduced control-oriented models and higher-fidelity OpenHPL-derived
> models side by side.

---

## 8. Initialization

For the dynamic reservoir, the main physical state is

```math
x = h.
```

A typical initial condition is

```julia
u0 = [
    sys.res.h => 50.0,
]
```

The remaining quantities are generated algebraically from the physical
relations.

This is preferable to independently guessing `A`, `V`, `m`,
`p_out`, and `h_abs`, which can make the MTK initialization system
overdetermined.

---

## 9. Validation

The initial tests check:

1. `ReservoirBoundaryZ` compiles with the new `HydraulicPortZ` connector.
2. A fixed reservoir imposes hydrostatic pressure and elevation consistently.
3. `DynamicReservoir` compiles under MTK 11.
4. Constant water withdrawal causes `h(t)` to decrease.
5. Absolute free-surface elevation is consistent with `z_out + h`.

These tests are intended to become the standard pattern for every refactored
physical component:

```text
equations
   ↓
standalone compile test
   ↓
steady-state / initialization test
   ↓
disturbance simulation
   ↓
comparison with reduced model
   ↓
comparison with OpenHPL where practical
```

---

## 10. Next component

The next hydraulic component should be the headrace/rigid pipe.

It should use two `HydraulicPortZ` connectors:

```text
Reservoir
  │
  ▼
port_a
┌──────────────┐
│  RigidPipe   │
└──────────────┘
port_b
  │
  ▼
Surge Tank
```

The first pipe model should include:

- flow inertia,
- gravity/elevation difference,
- Darcy-Weisbach friction,
- mass conservation,
- explicit sign conventions.

A suitable first momentum equation is

```math
\frac{L}{A}\frac{d\dot m}{dt}
=
p_a - p_b
+
\rho g(z_a-z_b)
-
\Delta p_f.
```

That component will create the first physically meaningful connection between
the refactored reservoir and the OpenHPL-style surge tank.
