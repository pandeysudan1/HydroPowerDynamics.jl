# FCR prequalification study

This folder turns the report workflow into one executable engineering study using `HydroPowerDynamics.jl`, `ControlSystemsBase.jl`, `JuMP.jl` and `Ipopt.jl`.

The engineering chain is

```text
nonlinear HPD plant
      -> local linear screening model
      -> frequency-domain analysis
      -> JuMP reserve-capacity search
      -> nonlinear HPD replay
      -> PASS / REDUCE / RETUNE
```

This is a **prequalification aid**, not an official Statnett qualification decision. The present numerical limits are deliberately transparent proxy limits. Product-specific Nordic/Statnett trajectories, operating points, acceptance envelopes, measurement-path treatment and reporting requirements should later be supplied as an external qualification-data layer.

---

## 1. Research question

The practical question is:

> **How much FCR can a hydropower unit offer at a given operating point while respecting actuator, hydraulic and dynamic-response limits, and can a cheap linear/JuMP screening result be trusted when it is replayed on the nonlinear HydroPowerDynamics.jl model?**

The motivation is to reduce repeated manual prequalification work. Instead of manually trying one droop value after another, we first formulate the reserve capability as an optimization problem. The optimizer returns a candidate reserve and the candidate is then checked on the nonlinear plant.

The method therefore separates two questions:

1. **Screening:** what is the largest reserve that appears feasible in a reduced local model?
2. **Verification:** does that candidate remain feasible when nonlinear hydraulic physics are restored?

A disagreement between the two models is itself a useful engineering result because it identifies the need for derating, controller retuning, or a better screening model.

---

## 2. Operating point and test input

The first case is a Trollheim-inspired 150 MW hydro unit operated around

- base power: `P_base = 150 MW`,
- initial mechanical power: `P0 = 75 MW`,
- nominal frequency: `f0 = 50 Hz`,
- initial guide vane: `y0 = 0.53848 pu`,
- initial flow: `Q0 = 22.1375 m3/s`,
- gross head: `H = 371 m`,
- governor actuator time constant: `Tg = 0.30 s`.

The proxy prequalification disturbance is a frequency step

$$
\Delta f(t)=
\begin{cases}
0, & t<0,\\
-0.10\ \mathrm{Hz}, & t\ge 0.
\end{cases}
$$

The current proxy requirement asks for at least 90% of the offered reserve after 5 s and approximately the full reserve by the end of the simulation. These values are used to exercise the workflow and are **not claimed as current Statnett limits**.

---

## 3. Reduced linear model

The local hydraulic screening model uses normalized flow deviation `q` and guide-vane position `y`:

$$
\dot{\Delta q}
=
\frac{1}{T_w}
\left(
\frac{\Delta y}{y_0}-\Delta q
\right),
$$

$$
\dot{\Delta y}
=
\frac{1}{T_g}
\left(-K_f\Delta f-\Delta y\right),
$$

with incremental power

$$
\Delta P=P_0\Delta q.
$$

The water starting time is estimated from the headrace and penstock geometry,

$$
T_w =
\frac{L_h Q_0}{g A_h H}
+
\frac{L_p Q_0}{g A_p H},
$$

which gives approximately

$$
T_w \approx 0.350\ \mathrm{s}
$$

for the present case.

In state-space form,

$$
\dot x = Ax+B\Delta f,
\qquad
\Delta P=Cx,
$$

where

$$
x=
\begin{bmatrix}
\Delta q\\
\Delta y
\end{bmatrix},
$$

$$
A=
\begin{bmatrix}
-1/T_w & 1/(T_w y_0)\\
0 & -1/T_g
\end{bmatrix},
\qquad
B=
\begin{bmatrix}
0\\
-K_f/T_g
\end{bmatrix},
\qquad
C=
\begin{bmatrix}
P_0 & 0
\end{bmatrix}.
$$

This model is intentionally cheap. It is used for poles, Bode/Nyquist inspection and the JuMP reserve search. It is not treated as the final physical truth.

---

## 4. FCR reserve optimization in JuMP

The decision variable is the offered reserve

$$
P_{\mathrm{FCR}} \ge 0.
$$

The optimization problem is

$$
\boxed{
\max P_{\mathrm{FCR}}
}
$$

subject to the discretized governor and hydraulic dynamics

$$
y_{k+1}
=
y_k+
\frac{\Delta t}{T_g}
\left(y_k^{cmd}-y_k\right),
$$

$$
q_{k+1}
=
q_k+
\frac{\Delta t}{T_w}
\left(\frac{y_k}{y_0}-q_k\right),
$$

where the FCR command is proportional to the test frequency deviation,

$$
y_k^{cmd}
=
y_0-
\frac{P_{\mathrm{FCR}}}{P_0}
\frac{\Delta f_k}{\Delta f_{full}}.
$$

The candidate must respect guide-vane bounds

$$
y_{min}\le y_k\le y_{max},
$$

flow bounds

$$
q_{min}\le q_k\le q_{max},
$$

guide-vane opening and closing rates

$$
-\dot y_{close}
\le
\frac{y_{k+1}-y_k}{\Delta t}
\le
\dot y_{open},
$$

and dynamic delivery requirements

$$
P_0(q_{5s}-1)
\ge
\alpha_{req}P_{\mathrm{FCR}},
$$

$$
P_0(q_{end}-1)
\ge
0.98P_{\mathrm{FCR}}.
$$

The current proxy bounds are

$$
y\in[0.05,1.0],
\qquad
q\in[0.70,1.45],
$$

and

$$
|\dot y|\le0.12\ \mathrm{pu/s}.
$$

### Linear optimization result

For the present equations the reserve search is dominated by the **guide-vane opening-rate constraint**. Immediately after the frequency step,

$$
\dot y(0^+)\approx
\frac{P_{\mathrm{FCR}}/P_0}{T_g}.
$$

Therefore

$$
\frac{P_{\mathrm{FCR}}/75}{0.30}\le0.12,
$$

which gives

$$
\boxed{P_{\mathrm{FCR}}^*\approx2.7\ \mathrm{MW}}.
$$

This is a useful sanity check for the JuMP solution: if the optimizer returns approximately 2.7 MW, the active constraint has a clear physical meaning rather than being a black-box numerical result.

For this candidate the reduced model predicts approximately

- maximum guide-vane opening: `0.5745 pu`,
- maximum normalized flow: `1.0669 pu`,
- 5 s incremental power: `5.01 MW`,
- maximum opening rate: `0.12 pu/s`.

The large difference between the offered reserve (2.7 MW) and the reduced-model incremental power is a warning that the present first-order screening normalization is deliberately simplified. This is exactly why the nonlinear HPD replay is required before interpreting the reserve as plant capability.

---

## 5. Mapping the optimized reserve to governor droop

The candidate reserve is converted into the primary governor droop used by the nonlinear HPD model:

$$
R^*
=
\frac{\Delta f_{full}/f_0}
{P_{\mathrm{FCR}}^*/P_0}.
$$

For the linear candidate above,

$$
R^*
\approx
\frac{0.10/50}{2.7/75}
\approx0.0556,
$$

or about

$$
\boxed{R^*\approx5.56\%}.
$$

Integral action is disabled in this test (`Ki = 0`) so that the response represents primary frequency action rather than slower secondary restoration.

---

## 6. Nonlinear HydroPowerDynamics.jl verification

The nonlinear verification restores the physical chain

```text
Reservoir -> headrace -> penstock -> Francis turbine -> prescribed shaft-frequency boundary
                                      ^
                                      |
                                  FCR governor
```

The frequency trajectory used in the screening test is imposed at the shaft boundary. The nonlinear model then computes the hydraulic flow, friction, turbine power and guide-vane response rather than assuming the linearized relationship remains exact.

The candidate is accepted only if the nonlinear trajectory satisfies all of the following checks:

$$
\Delta P(5s)\ge0.90P_{\mathrm{FCR}}^*,
$$

$$
y_{min}\le y(t)\le y_{max},
$$

$$
q_{min}\le Q(t)/Q_0\le q_{max},
$$

and

$$
-\dot y_{close}\le\dot y(t)\le\dot y_{open}.
$$

The notebook writes the final engineering decision as either

```text
PASS
```

or

```text
REDUCE / RETUNE
```

The nonlinear result, not the linear optimization alone, is the important prequalification result.

---

## 7. Plot results and interpretation

### 7.1 Linear Bode magnitude

![Linear Bode magnitude](plots/01_linear_bode_magnitude.png)

**What to read:** this plot shows how strongly a frequency disturbance is converted into incremental mechanical power as disturbance frequency changes.

**Interpretation:** the low-frequency gain represents the quasi-steady FCR sensitivity. The magnitude rolls off when the disturbance becomes faster than the governor/hydraulic dynamics. A plant can therefore have adequate steady-state reserve but still fail a fast dynamic requirement.

### 7.2 Linear Bode phase

![Linear Bode phase](plots/02_linear_bode_phase.png)

**What to read:** phase lag quantifies delay between frequency deviation and delivered hydropower response.

**Interpretation:** increasing phase lag means the turbine response arrives later relative to the grid-frequency disturbance. This is important because prequalification is not only an MW-capacity problem; timing matters.

### 7.3 Nyquist trajectory

![Linear Nyquist](plots/03_linear_nyquist.png)

**What to read:** the Nyquist trajectory is a compact frequency-domain representation of gain and phase together.

**Interpretation:** in this notebook it is primarily a screening diagnostic. When a more complete closed-loop grid/plant model is added, the trajectory can be used for explicit robustness and stability-margin assessment.

### 7.4 JuMP optimized linear response

![JuMP linear response](plots/04_jump_linear_response.png)

**What to read:** the solid curve is the predicted incremental power and the dashed level is the reserve offered by JuMP.

**Interpretation:** the optimized candidate should meet the dynamic-delivery constraint without violating plant constraints. For the present proxy case, the reserve is expected to be limited mainly by guide-vane opening rate rather than available gate travel or flow capacity.

### 7.5 Optimized guide-vane trajectory

![Optimized gate](plots/05_jump_gate.png)

**What to read:** this is the actuator motion required to provide the optimized reserve.

**Interpretation:** a trajectory touching the rate limit means that the actuator, not the turbine MW rating, is the binding capability. This result is operationally useful because it tells the engineer what should be improved or retuned if more FCR is desired.

### 7.6 Nonlinear HPD power verification

![Nonlinear power](plots/06_nonlinear_power_verification.png)

**What to read:** nonlinear turbine power is compared with the JuMP reserve offer.

**Interpretation:** if the nonlinear power response reaches the required fraction within the required time and remains physically well behaved, the linear candidate survives the first physics-based validation. If not, the offer should be reduced or the controller/screening model revised.

### 7.7 Nonlinear guide-vane response

![Nonlinear gate](plots/07_nonlinear_gate.png)

**What to read:** this plot verifies that the governor trajectory produced by the actual nonlinear simulation remains inside gate and gate-rate limits.

**Interpretation:** disagreement with the optimized linear gate trajectory indicates that the reduced model is missing important nonlinear or hydraulic interactions.

### 7.8 Nonlinear turbine flow

![Nonlinear flow](plots/08_nonlinear_flow.png)

**What to read:** this is the hydraulic cost of the frequency response in terms of turbine discharge.

**Interpretation:** FCR feasibility is not only a generator-power question. Excessive flow excursion can become the real constraint even if the generator has sufficient MW headroom.

### 7.9 Linear screening versus nonlinear physics

![Linear versus nonlinear](plots/09_linear_vs_nonlinear.png)

This is the central result of the study.

**Interpretation:**

- close agreement means the reduced model is suitable for cheap reserve screening near this operating point;
- moderate mismatch means a safety factor or model correction may be enough;
- large mismatch means the linear candidate should not be used directly and the workflow should return `REDUCE / RETUNE`;
- repeated mismatch across operating points is evidence that the screening model itself should be upgraded.

The purpose is therefore not to force the nonlinear plant to agree with the optimizer. The purpose is to use the nonlinear plant to determine when the cheap optimizer can be trusted.

---

## 8. Why JuMP helps the prequalification workflow

A conventional workflow can involve repeated manual tests:

```text
choose droop -> simulate -> inspect -> change droop -> simulate again -> inspect again
```

The proposed workflow changes this to

```text
plant limits + test definition
          |
          v
JuMP searches maximum candidate reserve
          |
          v
one nonlinear verification layer
          |
          v
PASS / REDUCE / RETUNE + limiting constraint
```

The benefit is not that JuMP replaces the physical qualification test. The benefit is that it reduces the search space before expensive/high-fidelity testing and explains **why** the capability is limited.

This can make prequalification preparation less repetitive and potentially cheaper, especially when many operating points, hydro units or controller settings must be screened.

---

## 9. Generated artifacts

After a successful CI run the study contains

```text
quick_examples/fcr_prequalification_study/
├── FCR_Prequalification_Study.jl
├── FCR_Prequalification_Study.ipynb
├── README.md
├── plots/
│   ├── 01_linear_bode_magnitude.png
│   ├── 02_linear_bode_phase.png
│   ├── 03_linear_nyquist.png
│   ├── 04_jump_linear_response.png
│   ├── 05_jump_gate.png
│   ├── 06_nonlinear_power_verification.png
│   ├── 07_nonlinear_gate.png
│   ├── 08_nonlinear_flow.png
│   └── 09_linear_vs_nonlinear.png
└── results/
    ├── fcr_candidate_summary.csv
    └── fcr_validation_trajectory.csv
```

The summary CSV records

- optimized reserve,
- optimized droop,
- linear 5 s delivery,
- nonlinear 5 s delivery,
- nonlinear final delivery,
- maximum gate,
- maximum normalized flow,
- maximum guide-vane rate,
- final decision.

---

## 10. Reproduce locally

```bash
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
pip install jupytext nbconvert jupyter
jupytext --to ipynb quick_examples/fcr_prequalification_study/FCR_Prequalification_Study.jl
jupyter nbconvert --to notebook --execute --inplace \
  quick_examples/fcr_prequalification_study/FCR_Prequalification_Study.ipynb \
  --ExecutePreprocessor.timeout=1800
```

---

## 11. Next step toward actual FCR prequalification

Keep the physics/optimization architecture fixed and replace the proxy test configuration with the applicable FCR-N/FCR-D requirements. Then repeat the workflow across several operating points:

$$
(P_0,H,Q_0)^{(1)},
(P_0,H,Q_0)^{(2)},\ldots,(P_0,H,Q_0)^{(n)}.
$$

This produces a capability surface rather than one number,

$$
P_{\mathrm{FCR,max}} =
\Phi(P_0,H,Q_0,\theta_{gov},\theta_{hyd}).
$$

That capability map is the more useful long-term result: it can support automatic prequalification preparation, controller tuning and a Freki-type decision-support layer around HydroPowerDynamics.jl.