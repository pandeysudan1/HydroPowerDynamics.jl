# FREKI Objective 1B — Validation of HydroPowerDynamics.jl with Trollheim HPP Measurement Data

This note continues Objective 1 from the synthetic normal-operation experiment and moves the same philosophy to measured Trollheim hydropower-plant data already stored in this repository.

The aim is not to claim that the full FREKI validation problem is solved. The aim is to show how existing plant measurements can be used to check and identify the physical parts of `HydroPowerDynamics.jl` before the model is used for FCR-capability studies.

```text
synthetic Objective 1
        |
        v
check identification logic
        |
        v
Trollheim plant measurements
        |
        v
fit physical turbine parameters
        |
        v
validate measured startup / steady operation
        |
        v
identify controller + hydraulic parameters
        |
        v
nonlinear HPD replay on held-out data
        |
        v
model health for FCR studies
```

---

## 1. Measurement record available in the repository

The existing experimental-data workflow uses one hour of measured Trollheim HPP data sampled at 1 s.

Source and processed data:

```text
quick_examples/data/Trollheim_data.xlsx
quick_examples/trollheim_experimental/01_preprocessing/trollheim_processed.csv
```

The record contains the following measured channels.

| Plant measurement | Raw unit | HPD quantity |
|---|---:|---|
| turbine inlet pressure | kPa | `p_in` [Pa] |
| turbine outlet pressure | kPa | `p_out` [Pa] |
| tailwater level | m | `h_tailwater_m` |
| shaft speed | RPM | `omega` [rad/s] |
| guide-vane servo position | mm | `tau_o` [pu] |
| generator active power | MW | `P_elec` [W] |
| tunnel discharge | m3/s | `Q` [m3/s] |

The preprocessing code converts the measurements into the variables used by `HydroPowerDynamics.jl`.

$$
H(t)=\frac{p_{in}(t)-p_{out}(t)}{\rho g},
$$

$$
\dot m(t)=\rho Q(t),
$$

$$
P_{mech}(t)=\frac{P_{elec}(t)}{\eta_{gen}},
$$

$$
\eta(t)=\frac{P_{mech}(t)}{\rho g Q(t)H(t)}.
$$

The measured guide-vane signal is normalized using 99.134 mm as the full-open servo position in the existing preprocessing workflow.

![Measured Trollheim operating record](../trollheim_experimental/images/trollheim_overview.png)

---

## 2. What is contained in the one-hour record

The existing analysis divides the record into four useful operating regions.

| Time [s] | Operating condition | Why it is useful for model validation |
|---:|---|---|
| 0–850 | plant off | sensor offsets and zero-flow checks |
| 850–1200 | startup ramp | dynamic flow, gate, power and speed response |
| 1200–2400 | near steady full-load operation | turbine-parameter identification |
| 2400–3600 | shutdown ramp | independent dynamic validation opportunity |

The reported mean operating point in the steady-state window is approximately

| Quantity | Measured value |
|---|---:|
| net head | 384.9 m |
| discharge | 36.36 m3/s |
| guide-vane opening | 0.928 pu |
| generator power | 128.26 MW |
| hydraulic efficiency | 0.9484 |
| shaft speed | about 375 RPM |

This is already much stronger evidence than the synthetic experiment because the model is now being compared with a real hydraulic machine and real measurement channels.

---

## 3. Step A — validate the turbine physics first

Before estimating governor dynamics, isolate the turbine from the controller.

Feed the measured guide-vane opening and measured net head directly into the `FrancisTurbineAffinity` equations:

$$
Q=\tau_o K_q D^2\sqrt{H},
$$

$$
\eta=\eta_{max}\left[1-c_\eta\left(\frac{Q}{Q_r}-1\right)^2\right],
$$

$$
P_{mech}=\rho gQH\eta.
$$

Using the measured steady-state interval from 1200 to 2400 s, the existing parameter-fitting study gives

| Parameter | Earlier design value | Fitted from Trollheim measurements |
|---|---:|---:|
| $K_q$ | 0.3415 | **0.319521** |
| $\eta_{max}$ | 0.97 | **0.948541** |
| $c_\eta$ | 0.25 | **2.98336** |
| $Q_r$ | 37.0 m3/s | **36.36 m3/s** |

The reported fit errors are

| Output | RMSE |
|---|---:|
| discharge $Q$ | **0.07142 m3/s** |
| hydraulic efficiency $\eta$ | **0.001614** |
| mechanical power | **0.1421 MW** |

![Trollheim turbine parameter fit](../trollheim_experimental/03_parameter_fitting/images/parameter_fit.png)

### Interpretation

This result says that the simple HPD Francis turbine equations can reproduce the measured full-load operating region closely after the turbine parameters are identified from plant data.

It does **not** yet validate the governor, waterway inertia, surge-tank dynamics or FCR response. It validates the turbine conversion from measured head and gate position to flow and power around the measured operating point.

---

## 4. Step B — use the startup ramp as a dynamic validation case

The startup interval contains natural excitation without imposing an artificial FCR test.

The existing open-loop study uses approximately

```text
t = 800–1200 s
```

and feeds measured guide-vane motion into the identified turbine equations. This checks whether the fitted turbine parameters remain useful during a large transient rather than only at steady state.

![Measured versus modelled startup response](../trollheim_experimental/02_openloop_validation/images/openloop_comparison.png)

For FREKI Objective 1, this is important because the preferred sequence is

```text
identify from one operating region
        ->
replay another operating region
        ->
check residuals on data not used for fitting
```

A model should not be called validated only because it fits the same samples used for parameter estimation.

---

## 5. Step C — extend the measured-data study to controller and waterway parameters

The synthetic Objective 1 experiment currently works with controller and hydraulic parameters such as

$$
\theta=[R,T_g,f_D].
$$

The measured Trollheim record provides enough channels to begin the same exercise, but the identifiability of each parameter must be checked before trusting the estimate.

### Controller block

For the governor servo,

$$
T_g\dot u_g=u_0+K_i\xi-u_g+\frac{e_f}{R}.
$$

The measured record contains shaft speed and guide-vane position, so a controller-identification test is possible if the actual governor operating mode and controller reference signal during the record are known.

This point is important: shaft speed and guide-vane motion alone are not sufficient to claim a unique physical $R$ or $T_g$ unless the controller structure and relevant reference/set-point signals are consistent with the model.

### Hydraulic block

For a penstock segment in HPD,

$$
\frac{L}{A}\frac{d\dot m}{dt}
=
\Delta p
-
f_D\frac{L}{2D\rho A^2}\dot m|\dot m|.
$$

The measured data already provide pressure difference and discharge. Therefore an effective friction parameter can be estimated over sufficiently dynamic windows and then replayed in the nonlinear HPD model.

The useful physical parameter set becomes

$$
\boxed{
\theta_{plant}
=
[K_q,\eta_{max},c_\eta,Q_r,f_D,\ldots]
}
$$

and, when the governor signals are sufficiently known,

$$
\boxed{
\theta_{control}
=
[R,T_g,K_i,\ldots].
}
$$

---

## 6. Proposed identification / validation split for the real data

Do not use the full one-hour record for both fitting and evaluation.

A practical first split is

```text
1200–1800 s  parameter fitting / calibration
1800–2400 s  held-out steady-state validation
850–1050 s   dynamic identification / startup calibration
1050–1200 s  held-out startup validation
2400–3600 s  independent shutdown validation
```

The exact windows should be adjusted only after checking data quality and whether the controller mode changes inside a window.

For each held-out window calculate at least

$$
RMSE_Q,
\quad RMSE_P,
\quad RMSE_\omega,
\quad RMSE_{\tau_o},
$$

plus normalized FIT values and residual plots.

The final decision should not be a single number. It should report which part of the model is supported by the measurements.

| Model block | Evidence required | Example status |
|---|---|---|
| turbine conversion | measured $H,Q,\tau_o,P$ | GREEN if held-out residuals are small |
| waterway | measured pressure + flow transients | GREEN / AMBER / RED |
| governor | speed + gate + known controller mode/reference | GREEN / AMBER / RED |
| complete FCR model | all above + frequency-response test | not yet established by this one-hour record |

---

## 7. Why this matters for FREKI and FCR prequalification

The practical question is not whether a simulation can be made to look similar to one measurement trace. The question is whether the model is trustworthy enough to make an engineering decision before a formal plant test.

The intended chain is

```text
normal Trollheim measurements
        |
        v
identify physical parameters
        |
        v
held-out nonlinear HPD validation
        |
        v
model-health / uncertainty statement
        |
        v
model-based FCR capability screening
        |
        v
select plant-test operating points
        |
        v
formal Statnett / Nordic prequalification tests
```

This is where academia can reduce practical effort: the model is used to narrow the search space, expose weakly identified physics and estimate which test cases are worth taking to the plant. It does not replace the formal TSO prequalification test.

---

## 8. Current evidence and next implementation

### Already available

- one hour of Trollheim plant measurements at 1 s resolution;
- processed HPD-aligned CSV data;
- measured startup, steady-state and shutdown regions;
- fitted Francis turbine parameters;
- steady-state fit errors for $Q$, $\eta$ and $P_{mech}$;
- open-loop startup comparison;
- synthetic FREKI Objective 1 identification workflow.

### Next executable study

Create

```text
objective_01b_trollheim_measurement_validation.jl
```

with the workflow

```text
load trollheim_processed.csv
        ->
select calibration and held-out windows
        ->
fit turbine / hydraulic parameters
        ->
run nonlinear HPD replay
        ->
compute residuals and FIT metrics
        ->
write model-health table
        ->
generate measurement-vs-model plots
```

Recommended outputs:

```text
plots/
├── 08_trollheim_measurement_overview.png
├── 09_trollheim_flow_validation.png
├── 10_trollheim_power_validation.png
├── 11_trollheim_speed_validation.png
├── 12_trollheim_gate_validation.png
├── 13_trollheim_hydraulic_residuals.png
└── 14_trollheim_model_health.png

results/
├── objective_01b_trollheim_summary.csv
└── objective_01b_trollheim_timeseries.csv
```

The result should end with a statement such as

```text
TURBINE MODEL      GREEN
HYDRAULIC MODEL    AMBER
GOVERNOR MODEL     NOT YET IDENTIFIABLE FROM AVAILABLE SIGNALS
FCR CAPABILITY     DO NOT USE UNTIL REQUIRED BLOCKS ARE GREEN
```

The actual labels must come from the executed residual tests; they should not be assigned in advance.

---

## Existing repository material used by this note

- `quick_examples/trollheim_experimental/ANALYSIS.md`
- `quick_examples/trollheim_experimental/01_preprocessing/preprocess_trollheim.jl`
- `quick_examples/trollheim_experimental/01_preprocessing/trollheim_processed.csv`
- `quick_examples/trollheim_experimental/02_openloop_validation/openloop_validation.jl`
- `quick_examples/trollheim_experimental/03_parameter_fitting/fit_turbine_params.jl`
- `quick_examples/trollheim_experimental/03_parameter_fitting/fitted_params.toml`
- `quick_examples/freki_objectives/objective_01_model_validation.jl`
