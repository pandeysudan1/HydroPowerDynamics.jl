# FREKI objectives with HydroPowerDynamics.jl

This study uses the existing Trollheim-inspired AGC model in `HydroPowerDynamics.jl` as a practical test bed for the five FREKI objectives defined by SINTEF.

FREKI aims to simplify FCR prequalification by using verified simulation models and normal operating data instead of relying only on dedicated field tests. The purpose here is not to reproduce the full FREKI project or claim official Statnett qualification. The goal is to build a transparent, executable research workflow around the same engineering questions.

Reference project: SINTEF FREKI, 2024-2026.

---

## 1. Existing plant model

The present branch already contains a nonlinear Trollheim-inspired hydropower model with AGC:

```text
Reservoir -> headrace -> [surge tank] -> penstock -> Francis turbine
                                                   |
                                                   v
                                             rotating mass
                                                   |
                                                   v
                                            generator / load
                                                   ^
                                                   |
                                             governor / AGC
                                                   ^
                                                   |
                                                speed
```

The controller implemented in `src/agc.jl` is

$$
e_f = \frac{\omega_{ref}-\omega}{\omega_{ref}},
$$

$$
\dot{\xi}=e_f,
$$

$$
u_{cmd}=u_0+\frac{e_f}{R}+K_i\xi,
$$

$$
T_g\dot u_g=\operatorname{sat}(u_{cmd})-u_g.
$$

The nonlinear hydraulic model retains waterway dynamics, friction, turbine flow, mechanical power and the optional surge tank. This is useful because FREKI is fundamentally a model-validation problem, not only a governor-tuning problem.

---

# Objective 1 — validate the plant model from normal operation

## FREKI question

Can a simulation model be validated from real-time measurements without disturbing normal plant operation?

## Study idea

Use ordinary load and frequency variations as excitation. For the first implementation, the nonlinear AGC model itself generates a synthetic `measured` data set. Noise and small parameter offsets are then added so that a second HPD model must recover the plant behaviour.

The measured vector is

$$
y_m(t)=
\begin{bmatrix}
f & P_e & P_m & u_g & Q & H\end{bmatrix}^T.
$$

The model prediction is

$$
\hat y(t,\theta)=\mathcal{M}_{HPD}(u(t),\theta),
$$

where a first parameter vector can be

$$
\theta=
\begin{bmatrix}
R & T_g & K_i & f_{pipe} & \eta_{max} & H_g
\end{bmatrix}^T.
$$

Estimate the parameters from

$$
\theta^*=\arg\min_{\theta}
\sum_k
\left(y_m(k)-\hat y(k,\theta)\right)^T
W
\left(y_m(k)-\hat y(k,\theta)\right).
$$

### First outputs

```text
plots/01_measured_vs_model_frequency.png
plots/02_measured_vs_model_power.png
plots/03_measured_vs_model_gate_flow.png
plots/04_validation_residuals.png
plots/05_parameter_convergence.png
```

Report at least

| Metric | Meaning |
|---|---|
| RMSE frequency | dynamic frequency mismatch |
| RMSE power | active-power mismatch |
| RMSE gate | actuator mismatch |
| RMSE flow | hydraulic mismatch |
| FIT % | overall model fit |
| parameter error | identified vs reference value |

The important point is that the validation data should come from small naturally occurring changes, not from a dedicated large disturbance.

---

# Objective 2 — determine qualified FCR from the verified model

## FREKI question

Once the plant model is trusted, how much FCR can the unit safely offer?

The existing FCR study already follows the useful engineering chain

```text
validated nonlinear plant
        -> reduced model
        -> candidate FCR
        -> governor droop
        -> nonlinear replay
        -> PASS / REDUCE / RETUNE
```

The decision variable is

$$
P_{FCR}\ge 0.
$$

A simple capacity formulation is

$$
\max P_{FCR}
$$

subject to

$$
y_{min}\le y(t)\le y_{max},
$$

$$
|\dot y(t)|\le \dot y_{max},
$$

$$
Q_{min}\le Q(t)\le Q_{max},
$$

plus the relevant Statnett/Nordic response envelope.

The important extension for FREKI is uncertainty. A validated model is never exact. Therefore define a validated parameter region

$$
\Theta_{val}=\{\theta: J(\theta)\le J_{max}\}.
$$

Then estimate a robust FCR capability

$$
\boxed{
P_{FCR}^{robust}
=
\min_{\theta\in\Theta_{val}}
P_{FCR}^{max}(\theta)
}
$$

rather than relying only on one nominal parameter set.

### Outputs

```text
plots/06_fcr_response_envelope.png
plots/07_fcr_capacity_vs_model_uncertainty.png
results/fcr_capacity_summary.csv
```

The final result should answer three questions:

- how much FCR can be offered;
- which physical/control constraint limits it;
- how sensitive that number is to model uncertainty.

---

# Objective 3 — find measures that increase qualified FCR

## FREKI question

If the plant is limited to a certain FCR capacity, what should be changed to obtain more?

Treat the FCR capability as

$$
P_{FCR}^{max}
=F(R,T_g,K_i,\dot y_{max},P_0,H,T_w,\theta_h,\ldots).
$$

The first sensitivity study should vary

- governor droop `R`,
- servo time constant `T_g`,
- guide-vane rate limit,
- operating power `P0`,
- gross head,
- waterway dynamics,
- surge-tank representation.

For each parameter $p_i$ calculate a local sensitivity

$$
S_i=
\frac{\partial P_{FCR}^{max}}{\partial p_i}.
$$

A practical decision table can then be generated:

| Change | FCR effect | Frequency effect | Hydraulic effect | Action |
|---|---:|---:|---:|---|
| lower `T_g` | TBD | faster response | higher gate activity | evaluate |
| change `R` | TBD | changes primary gain | moderate | tune |
| higher gate-rate limit | TBD | faster response | more actuator duty | check hardware |
| different operating point | TBD | changes headroom | changes flow/head | schedule |
| surge-tank dynamics | TBD | changes early transient | changes pressure/flow | model explicitly |

### Outputs

```text
plots/08_fcr_sensitivity_tornado.png
plots/09_fcr_capability_map.png
plots/10_operating_point_map.png
```

The purpose is to move from

```text
FAIL
```

to

```text
FAIL -> limiting mechanism -> engineering action -> new FCR capability
```

---

# Objective 4 — Pelton, Francis and Kaplan

## FREKI question

Can the same validation and FCR-capability method work across different turbine technologies?

Keep a common external interface:

```text
frequency / speed
      |
      v
governor + actuator
      |
      v
hydraulic turbine
      |
      v
mechanical power / torque
```

but change the internal turbine/actuator physics.

### Francis

Use the current Trollheim-inspired model first.

Relevant states/signals:

- guide-vane opening,
- penstock flow,
- head,
- turbine power,
- surge-tank level when present.

### Pelton

Future model extension:

```text
frequency -> governor -> needle / deflector -> jet -> runner torque
```

Important additional dynamics include needle travel, jet flow and deflector action.

### Kaplan

Future model extension:

```text
frequency -> governor -> guide vane
                    -> blade pitch
                    -> turbine power
```

Kaplan introduces coordinated guide-vane and runner-blade control.

The validation API should remain turbine-independent:

```julia
simulate(model, inputs, parameters)
validate(measurements, prediction)
estimate_fcr(validated_model, requirements)
```

### Output

| Turbine | Main actuator | Important hydraulic state | Expected FCR limitation |
|---|---|---|---|
| Francis | guide vane | penstock/surge dynamics | gate rate + water inertia |
| Pelton | needle/deflector | jet dynamics | actuator/jet response |
| Kaplan | guide vane + blade pitch | flow/head + coordinated control | multi-actuator dynamics |

---

# Objective 5 — online demonstration with real-time data

## FREKI question

Can model validity and FCR capability be monitored continuously?

Use the validated HPD model as an online observer/predictor:

$$
\hat x_{k+1}=F(\hat x_k,u_k,\hat\theta_k),
$$

$$
\hat y_k=G(\hat x_k,u_k,\hat\theta_k),
$$

and calculate residuals

$$
r_k=y_k-\hat y_k.
$$

Over a moving window calculate

$$
RMSE_f,
\quad
RMSE_P,
\quad
RMSE_Q,
\quad
RMSE_y.
$$

A simple model-health state can be

```text
GREEN  : model residuals inside validated band
AMBER  : model mismatch increasing
RED    : model no longer trusted for FCR estimation
```

The online output should show

```text
live measurements
      |
      +------> HPD model prediction
      |               |
      |               v
      +----------> residuals
                      |
                      v
                 model health
                      |
                      v
              estimated FCR capacity
```

### Outputs

```text
plots/11_online_model_residual.png
plots/12_online_fcr_capacity.png
results/model_health.csv
```

---

# Relationship to Statnett prequalification

Statnett and the Nordic TSOs require FCR-providing entities to be prequalified before market delivery. The current Nordic process includes verification of entity properties, prequalification tests, measurement/data requirements and technical documentation.

This repository does **not** replace that formal process.

The proposed role of HydroPowerDynamics.jl is instead

```text
Statnett / Nordic requirements
             |
             v
normal plant measurements
             |
             v
validate nonlinear HPD model
             |
             v
simulate required FCR behaviour
             |
             v
estimate capability + uncertainty
             |
             v
identify limiting mechanism
             |
             v
support prequalification engineering
```

---

# First implementation

The first executable study should address **Objective 1 only**.

Use the existing AGC Trollheim model to generate normal-operation data with small load variations rather than one large step. Treat those signals as measurements and identify a deliberately perturbed model.

Suggested disturbance:

$$
P_L(t)=P_0+\Delta P_{slow}(t)+\Delta P_{noise}(t),
$$

with changes small enough to represent normal operation.

Start with only three unknown parameters:

$$
\boxed{\theta=[R,T_g,f_{pipe}]}
$$

and measured signals

$$
\boxed{y=[f,P_e,u_g,Q]}.
$$

A good first success criterion is

```text
1. simulation runs from steady state;
2. noisy synthetic measurements are generated;
3. parameter estimation converges;
4. validation data are not reused for fitting;
5. measured and predicted traces are plotted;
6. RMSE/FIT metrics are reported;
7. the validated parameter uncertainty is passed to Objective 2.
```

---

# Proposed folder

```text
quick_examples/freki_objectives/
├── README.md
├── objective_01_model_validation.jl
├── objective_02_fcr_capacity.jl
├── objective_03_fcr_improvement.jl
├── objective_04_turbine_types.jl
├── objective_05_online_validation.jl
├── plots/
└── results/
```

For now, this README is the research specification. Numerical values should only be added after the corresponding scripts have been executed successfully.

---

# References

1. SINTEF, **FREKI**, Innovation Project for the Industrial Sector, 2024-2026.
2. Statnett / Nordic TSOs, **Technical Requirements for Frequency Containment Reserve Provision in the Nordic Synchronous Area**.
3. Statnett / Nordic TSOs, **Test Program for Prequalification of FCR in the Nordic Synchronous Area**.
4. `HydroPowerDynamics.jl`, branch `AGC_Trollheim`, existing Trollheim AGC and surge-tank studies.
