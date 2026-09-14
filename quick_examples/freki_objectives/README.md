# FREKI objectives with HydroPowerDynamics.jl

This branch extends the FCR prequalification study into a **FREKI-oriented research workflow** using `HydroPowerDynamics.jl` (HPD) as the main physics-based plant model.

The aim is not to reproduce SINTEF FREKI or claim official Statnett prequalification. The aim is to turn the five FREKI objectives into small, executable hydropower studies that can later be connected to real plant measurements.

```text
real plant data
     |
     v
parameter/state estimation
     |
     v
validated HPD model
     |
     +----> FCR capacity search
     |
     +----> limiting-constraint diagnosis
     |
     +----> Pelton / Francis / Kaplan comparison
     |
     +----> online model-vs-measurement monitoring
```

---

## 1. Research question

> Can a nonlinear hydropower model, validated from normal operating data, be used to estimate FCR capability, explain the limiting physics, and support online prequalification-oriented monitoring?

We treat the five FREKI objectives as five engineering experiments.

| FREKI objective | HPD question | Main output |
|---|---|---|
| 1. Validate models from real-time measurements | Does the HPD model reproduce measured power, gate, flow and frequency trajectories? | residual plots + fit metrics |
| 2. Determine qualified FCR | What reserve can the validated model deliver within dynamic and hydraulic limits? | MW capability + limiting constraint |
| 3. Increase qualified FCR | Which tuning or plant change increases feasible reserve? | sensitivity / Pareto plot |
| 4. Pelton, Francis, Kaplan | Does the same workflow work across turbine types? | cross-technology comparison |
| 5. Online demonstration | Can model-vs-measurement checks run continuously? | rolling residual / health indicators |

---

# Objective 1 — validate the power-plant simulation model

## Idea

Use measured operation rather than a dedicated disturbance test whenever the data contain enough excitation.

For measured signals

$$
y_m(t)=\{f,\;P,\;y_g,\;Q,\;H\},
$$

and HPD predictions

$$
\hat y(t,\theta),
$$

estimate model parameters by minimizing

$$
\boxed{
\theta^*=\arg\min_{\theta}
\sum_k
\left(y_m(k)-\hat y(k,\theta)\right)^T
W
\left(y_m(k)-\hat y(k,\theta)\right)
}
$$

where parameters may include

$$
\theta=\{T_g,\;R,\;K_f,\;f_{pipe},\;\eta_t,\;T_w,\ldots\}.
$$

### First implementation

Start with the existing Trollheim-inspired Francis model and create synthetic measurements from the nonlinear HPD model. Add realistic noise and small operating disturbances, then estimate a reduced parameter set.

Use three layers:

```text
truth HPD model
    -> synthetic measurements + noise
    -> candidate HPD model
    -> residuals / fitted parameters
```

### Metrics

For signal $i$,

$$
RMSE_i=
\sqrt{\frac{1}{N}\sum_{k=1}^{N}
(y_{i,k}-\hat y_{i,k})^2}
$$

and normalized fit

$$
FIT_i=
100\left(1-
\frac{\|y_i-\hat y_i\|_2}
{\|y_i-\bar y_i\|_2}
\right).
$$

The first plots should be:

1. measured vs simulated active power;
2. measured vs simulated guide-vane position;
3. measured vs simulated flow;
4. residuals versus time;
5. parameter estimate convergence.

### Decision

```text
small unbiased residuals -> model accepted for the tested operating region
structured residuals     -> model structure or parameter set must be improved
```

---

# Objective 2 — determine how much FCR can be qualified

This objective directly continues the existing `FCR_study` branch.

The current workflow is

```text
nonlinear HPD plant
      -> local linear model
      -> frequency-domain screening
      -> JuMP reserve search
      -> nonlinear HPD replay
      -> PASS / REDUCE / RETUNE
```

The reserve problem is

$$
\boxed{\max P_{FCR}}
$$

subject to governor, waterway, turbine and qualification-oriented response constraints.

The first executed study found

$$
P_{FCR}^*=2.70\ \mathrm{MW},
\qquad
R^*=5.56\%,
$$

with the **guide-vane opening rate** as the active proxy constraint.

The important next step is to make the reserve estimate depend on the **validated parameter uncertainty** obtained in Objective 1.

Instead of one deterministic model,

$$
P_{FCR}^*(\theta),
$$

use an uncertainty set

$$
\theta\in\Theta_{validated}
$$

and calculate a robust capability

$$
\boxed{
P_{FCR}^{robust}
=
\min_{\theta\in\Theta_{validated}}
P_{FCR}^*(\theta)
}
$$

This provides a clean link between **model validation** and **reserve declaration**.

### Existing evidence from the FCR study

![JuMP optimized response](../fcr_prequalification_study/plots/04_jump_linear_response.png)

![Guide-vane trajectory](../fcr_prequalification_study/plots/05_jump_gate.png)

![Nonlinear HPD verification](../fcr_prequalification_study/plots/06_nonlinear_power_verification.png)

![Linear versus nonlinear model](../fcr_prequalification_study/plots/09_linear_vs_nonlinear.png)

The existing result already shows why model validation matters: the reduced and nonlinear models have similar actuator motion but different delivered mechanical power.

---

# Objective 3 — identify measures that increase qualified FCR

Once a validated model exists, do not only ask *how much FCR is possible?*

Ask

> **What must change to increase it?**

Define

$$
P_{FCR}^*=F(R,T_g,\dot y_{max},H,P_0,T_w,\theta_h,\ldots).
$$

Then calculate sensitivities such as

$$
S_x=
\frac{\partial P_{FCR}^*}{\partial x}.
$$

The first study should sweep:

| Variable | Engineering meaning |
|---|---|
| governor droop $R$ | strength of primary response |
| actuator constant $T_g$ | governor/servo speed |
| gate-rate limit $\dot y_{max}$ | actuator capability |
| operating power $P_0$ | available headroom / local operating point |
| hydraulic starting time $T_w$ | water-column dynamics |
| gross/net head | turbine operating condition |

For each case record

$$
\{P_{FCR}^*,\;\Delta P_{5s},\;\max|\dot y|,\;\max Q/Q_0,\;\min H\}.
$$

The target figure is a capability map such as

```text
                 feasible FCR
                      ^
                      |
          hydraulic  |       actuator
           limited   |       limited
                      |
                      +------------------> operating point
```

This turns the model from a pass/fail tool into a **plant-improvement tool**.

The current FCR case already gives the first diagnosis:

```text
guide-vane rate  -> active / near-active
flow limit       -> inactive
gate travel      -> inactive
minimum delivery -> satisfied
```

So the immediate engineering hypothesis is:

$$
\boxed{
\text{larger allowable gate rate or improved governor tuning}
\Rightarrow
\text{larger feasible FCR}
}
$$

subject to hydraulic and stability limits.

---

# Objective 4 — validate on Pelton, Francis and Kaplan plants

The framework should stay the same while the turbine physics change.

```text
frequency trajectory
        |
        v
    governor
        |
        v
+--------------------+
| turbine technology |
+--------------------+
   |       |       |
 Pelton  Francis  Kaplan
   |       |       |
   +-------+-------+
           |
           v
 power / flow / pressure / actuator limits
```

## Common interface

Each model should expose approximately

$$
u=\{f_{grid},P_{set},H_{up},H_{down}\}
$$

and outputs

$$
y=\{P_m,Q,H_t,y_g,\omega\}.
$$

### Francis

Use the existing Trollheim-inspired HPD model first. This is the present reference implementation.

### Pelton

The important additional dynamics are nozzle/needle actuation, jet flow and possible deflector logic.

### Kaplan

The model should include coupled guide-vane and runner-blade control, making the FCR capability operating-point dependent in two actuator dimensions.

### Comparison table to generate

| Turbine | $P_{FCR}^*$ | active limit | 5 s response | dominant hydraulic/actuator issue |
|---|---:|---|---:|---|
| Francis | TBD from validated study | TBD | TBD | TBD |
| Pelton | TBD | TBD | TBD | TBD |
| Kaplan | TBD | TBD | TBD | TBD |

The goal is not to force identical models. The goal is to use one **qualification workflow** around technology-specific physics.

---

# Objective 5 — demonstrate the method online with real-time data

The online version is a repeated model-validation problem.

At time $k$ receive

$$
\{f_k,P_k,y_{g,k},Q_k,H_k\}
$$

and propagate the HPD model

$$
\hat x_{k+1}=F(\hat x_k,u_k,\hat\theta_k).
$$

Calculate residuals

$$
r_k=y_k-\hat y_k.
$$

Then maintain rolling metrics such as

$$
RMSE_{30min},\qquad
bias_{30min},\qquad
\max|r|,
$$

and a simple model-health state

```text
GREEN  -> model remains inside validation envelope
AMBER  -> model mismatch increasing
RED    -> do not rely on current model for FCR capability declaration
```

A later Railway service can ingest a CSV/API stream, run the model and expose a compact status payload:

```json
{
  "model_status": "GREEN",
  "fcr_capability_mw": 2.7,
  "active_constraint": "gate_rate",
  "power_rmse_mw": 0.0,
  "last_update": "..."
}
```

The first online demonstration should use **simulated streaming data**. Real plant data can be connected later without changing the model-validation architecture.

---

# One integrated FREKI-oriented workflow

The five objectives are not separate notebooks. They form one chain:

```text
[1] normal-operation measurements
              |
              v
     validate HPD parameters
              |
              v
[2] calculate feasible FCR capability
              |
              v
[3] identify the active limitation
              |
              v
     tune / redesign / change operating point
              |
              v
[4] repeat for Francis, Pelton and Kaplan
              |
              v
[5] run validation and capability checks online
```

Mathematically, the full research problem can be written compactly as

$$
\theta^*=\arg\min_{\theta}J_{validation}(\theta)
$$

followed by

$$
P_{FCR}^*
=
\max_{P_{FCR}}
\left\{P_{FCR}:g(x,u,\theta^*)\le0\right\}.
$$

The research value is therefore the link

$$
\boxed{
\text{measurement}
\rightarrow
\text{validated nonlinear model}
\rightarrow
\text{FCR capability}
\rightarrow
\text{physical explanation}
}
$$

rather than simulation alone.

---

# Implementation sequence on this branch

The branch will be developed in this order:

1. **Objective 1:** generate synthetic normal-operation data from HPD and recover selected parameters;
2. **Objective 2:** reuse the existing FCR optimizer with the validated model;
3. **Objective 3:** run parameter/operating-point sweeps and produce an FCR capability map;
4. **Objective 4:** introduce Pelton and Kaplan model interfaces after the Francis workflow is stable;
5. **Objective 5:** expose the validation + capability calculation as a lightweight online Railway workflow.

The immediate next executable target is therefore **Objective 1**, because FREKI's central idea depends on trusting the model before using it to declare FCR capability.

---

## Scope statement

This repository is a research implementation inspired by the published FREKI objectives and Statnett-oriented FCR prequalification questions. It is **not** an official SINTEF FREKI implementation and does not itself constitute Statnett prequalification.
