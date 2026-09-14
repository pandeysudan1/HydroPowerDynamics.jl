# FREKI Objective 1B — Trollheim measurement validation results

This note records the first executed measurement-based validation result for the Trollheim HPP data already available in this repository. It continues `TROLLHEIM_MEASUREMENT_VALIDATION.md` and should be read together with `objective_01b_trollheim_measurement_validation.jl`.

The purpose is not to claim FCR prequalification. The purpose is to decide which parts of the HydroPowerDynamics.jl model are already supported by measured plant data, which parts remain weakly identified, and what additional measurements are needed before model-based FCR capability is trusted.

## 1. Executed study

The GitHub Actions workflow completed successfully on the processed one-hour Trollheim record.

```text
1200–1800 s   calibration
1801–2400 s   held-out steady-state validation
850–1200 s    independent startup check
```

The local Francis relation was

$$
Q=\tau_oK_qD^2\sqrt{H},
$$

$$
\eta=\eta_{max}\left[1-c_\eta\left(\frac{Q}{Q_r}-1\right)^2\right],
$$

$$
P_e=\eta_g\rho gQH\eta.
$$

No held-out samples were used for fitting.

## 2. Identified parameters

| Parameter | Identified value | Interpretation |
|---|---:|---|
| $K_q$ | 0.319751 | physically plausible local flow coefficient |
| $\eta_{max}$ | 0.948071 | consistent with the measured high-load efficiency level |
| $Q_r$ | 36.4445 m3/s | close to observed full-load discharge |
| $c_\eta$ | -1.84392 | **not physically credible as an efficiency-hill curvature** |

The negative $c_\eta$ is an identifiability warning. The high-load calibration window contains too little off-design excitation to determine the curvature of the efficiency hill reliably.

```text
good output fit != all physical parameters are identified
```

## 3. Held-out validation

| Metric | Result |
|---|---:|
| Flow RMSE | **0.07945 m3/s** |
| Relative flow RMSE | **0.2190 %** |
| Electrical-power RMSE | **0.15704 MW** |
| Relative power RMSE | **0.1227 %** |
| Flow FIT | 33.91 % |
| Power FIT | 54.61 % |

The relative RMSE values are very small. Around the measured full-load point, the local turbine map reproduces flow and power closely on data not used for fitting.

The normalized FIT values are lower because the held-out signals have little variation around their mean. In such a narrow operating window, absolute and relative RMSE are more informative than FIT alone.

![Held-out flow validation](plots/01b_03_heldout_flow.png)

![Held-out power validation](plots/01b_04_heldout_power.png)

![Held-out residuals](plots/01b_05_residuals.png)

## 4. Independent startup check

| Metric | Result |
|---|---:|
| Startup flow RMSE | **0.14864 m3/s** |
| Startup power RMSE | **4.31884 MW** |

![Startup power check](plots/01b_06_startup_power.png)

Flow remains reasonably reproduced, while power error increases markedly away from the calibration operating point. This shows where the simple local map stops being sufficient.

Possible causes include efficiency variation away from full load, transient hydraulic effects, generator-efficiency assumptions, dynamic actuator/turbine effects and signal timing. The current test does not yet separate these mechanisms.

## 5. Model-health interpretation

The first script labelled the turbine block `GREEN` because held-out relative flow and power errors were below 1 %. That is too optimistic if physical parameter plausibility is also considered.

| Model block | Evidence | Decision |
|---|---|---|
| Local turbine flow/power map near full load | excellent held-out RMSE | **GREEN locally** |
| Full turbine parameter set | $c_\eta<0$ from weak excitation | **AMBER** |
| Hydraulic dynamics | waterway states not yet independently validated | **AMBER** |
| Governor dynamics | controller mode/reference signals are missing | **AMBER / not uniquely identifiable** |
| FCR capability | complete dynamic chain not validated | **MORE DATA REQUIRED** |

```text
LOCAL TURBINE MAP      GREEN
TURBINE PARAMETER SET  AMBER
HYDRAULIC DYNAMICS     AMBER
GOVERNOR DYNAMICS      AMBER
FCR READINESS           MORE DATA REQUIRED
```

The important FREKI lesson is that residual quality and parameter identifiability must be reported separately.

## 6. Measurement sufficiency: an important hydraulic correction

The available pressure channels are

- `p_penstock_a_`: turbine inlet / spiral-casing pressure;
- `p_dt_b_`: turbine outlet / draft-tube pressure.

Their difference is useful for turbine net head,

$$
H_t=\frac{p_{in}-p_{out}}{\rho g},
$$

but these are **not the two end pressures of a penstock segment**. Therefore the current record must not be used to identify a penstock Darcy factor directly from

$$
\frac{L}{A}\frac{d\dot m}{dt}
=
\Delta p-f_D\frac{L}{2D\rho A^2}\dot m|\dot m|.
$$

Doing so would assign the turbine pressure drop to the penstock model and would be physically incorrect.

This changes the next FREKI task from "fit $f_D$ now" to **measurement sufficiency and hydraulic identifiability**.

## 7. What can be identified with the current Trollheim record?

| Quantity / model block | Current channels sufficient? | Comment |
|---|---|---|
| local $K_q$ | **Yes** | $H$, gate and $Q$ are measured |
| local efficiency level | **Yes** | $P$, $Q$ and $H$ are measured |
| efficiency curvature $c_\eta$ | **Weak** | full-load window gives too little off-design excitation |
| turbine net-head response | **Yes** | turbine inlet and outlet pressures are measured |
| penstock friction $f_D$ | **No** | penstock upstream pressure/head is missing |
| surge-tank dynamics | **No / incomplete** | surge level/branch-flow measurements are not in this record |
| governor $R,T_g$ | **Not uniquely** | controller reference/mode signals are missing |
| complete FCR dynamic model | **No** | required blocks are not all validated |

This table is a useful FREKI output by itself: it tells the engineer which extra sensors or logged control signals provide the highest value before another plant test is planned.

## 8. Next implementation — Objective 1C: measurement sufficiency and dynamic residuals

The next executable study should not invent a penstock parameter. It should use the existing data to quantify what the current measurements can support and what remains structurally unidentifiable.

The proposed next script is

```text
objective_01c_measurement_sufficiency.jl
```

with the workflow

```text
Trollheim measured channels
        |
        v
classify available physical equations
        |
        v
check excitation of each parameter
        |
        v
compute local sensitivities / condition numbers
        |
        v
rank identifiable vs weak vs unavailable parameters
        |
        v
recommend additional plant signals
```

A practical output table should look like

| Parameter | Sensitivity | Data support | Action |
|---|---:|---|---|
| $K_q$ | high | available | retain |
| $\eta_{max}$ | high near full load | available | retain |
| $c_\eta$ | low | weak excitation | collect part-load data |
| $f_D$ | not observable | missing penstock boundary pressure | add/log upstream pressure/head |
| surge parameters | not observable | missing surge states | log surge level/branch flow if available |
| $R,T_g$ | ambiguous | missing governor references/mode | export controller signals |

The recommended additional control channels are governor frequency reference, power/set-point reference, controller mode, droop/transient-droop settings and internal servo command. For waterway identification, at least one appropriate upstream hydraulic boundary measurement is required in addition to turbine-inlet pressure and flow.

## 9. FCR decision logic

```text
measured plant data
      -> local turbine validation
      -> measurement sufficiency / identifiability
      -> collect missing high-value signals
      -> hydraulic dynamic validation
      -> governor validation
      -> uncertainty region
      -> model-based FCR capacity screening
      -> selected physical plant tests
      -> Statnett / Nordic prequalification process
```

The role of HydroPowerDynamics.jl is to reduce unnecessary plant testing and expose missing physics before the formal test. It should not turn a locally fitted model into an automatic claim of qualification.

## Generated evidence

- `results/objective_01b_summary.csv`
- `results/objective_01b_validation_timeseries.csv`
- `plots/01b_01_measurement_windows.png`
- `plots/01b_02_calibration_flow.png`
- `plots/01b_03_heldout_flow.png`
- `plots/01b_04_heldout_power.png`
- `plots/01b_05_residuals.png`
- `plots/01b_06_startup_power.png`
- `plots/01b_07_model_health.png`
