# Prompt for Sonnet — Step 9: regress on the raw time-series measurements (NaCl)

You are implementing **Step 9**: redo the NaCl regression posed directly on the
**raw measured time series** — permeate conductivity $\kappa_p(t)$ and permeate
mass $m(t)$ — with **no transformation of the response**, and compare it to the
transformed-$\Delta c_{\mathrm m}$ fit of §2.1–2.2. Work in
`/Users/adowling/DowlingLab/Membranes/data3_regression`, env `data3-regression`.
This **populates the existing placeholder subsection** `sec:altform`
("Alternative formulation: regress on the raw time-series measurements") in the
NaCl section of `docs/reports/main.tex`, **explaining the mathematics inline**
(exploratory results section — math here is fine per the PI).

> **Run after** Step 8 (it reuses the AR(1) machinery) and the figure-compliance
> audit (finalized `plot_style`).

## Why (the errors-in-variables problem)
The current fit forms $\Delta c_{\mathrm m}=c_p J_w\,\ell/\Dm$ and regresses it —
so the noisy processed quantities $c_p$ and $J_w$ appear \emph{inside the
response}, and $c_p$ appears on both sides of the estimation. This mis-states the
uncertainty. The rigorous alternative places the error on the \emph{raw measured}
signal and predicts it from the model.

## Working forward model (Opus — state as working assumptions; exploratory)
Primary dataset: raw sheet `NF270_MC5.xlsx / 05.27.26_NaCl` (the raw twin of the
`(2)` sheet used in §2.1–2.2). Reuse `analysis/preprocessing/preprocess_nacl.py`
to load the raw time series and the embedded calibration.
- **Known inputs per time $t$:** retentate conductivity $\kappa_r(t)\to$ feed-bulk
  conc $c_f(t)$ via the linear calibration $c=(\kappa-a)/s$ (embedded MC5 values
  $a=-63.706,\ s=76.685$); permeate mass $m(t)\to$ water flux
  $J_w(t)=\frac{1}{\rho A_m}\frac{\mathrm dm}{\mathrm dt}$ (smoothed slope, as in
  `sec:preproc`). Feed interface $\cint(t)$ from $c_f$ (non-CP, matching the `(2)`
  convention — or note if you apply CP).
- **Model prediction of the measured $\kappa_p$:** the co-ion (Cl$^-$) flux balance
  ties $c_p$ to the feed side. With the corrected Donnan co-ion root
  $c_{\co}^{\mathrm m}(c^{\mathrm s};\theta)=\tfrac12(-|\chi|+\sqrt{|\chi|^2+4\dstar (c^{\mathrm s})^2})$,
  solute conservation gives the implicit equation for $c_p$:
  $$ J_w\,c_p \;=\; \frac{\Dm}{\ell}\Big[c_{\co}^{\mathrm m}\big(\cint;\theta\big)-c_{\co}^{\mathrm m}\big(c_p;\theta\big)\Big]. $$
  Solve this scalar equation for $c_p(\theta,t)$ at each $t$ (it has a unique
  physical root; use Brent/`scipy.optimize.brentq` on $[0,\cint]$). Then predict
  $\hat\kappa_p(t)=a+s\,c_p(\theta,t)$ (invert the calibration).
- **Response = the raw measured $\kappa_p(t)$**; residual $=\kappa_p^{\text{meas}}(t)-\hat\kappa_p(t)$.
  (Optionally also include a mass/flux residual channel; keep $\kappa_p$ as the
  primary. State clearly what you fit.)

## Statistical treatment
- Fit $\theta=(|\chi|,\dstar)$ by weighted least squares on the raw $\kappa_p$
  residuals with a **measurement-error model** for conductivity (heteroscedastic,
  e.g.\ $\sigma_{\kappa}\propto$ magnitude or a constant floor + proportional term,
  as in DATA~2.0). State the error model.
- Report estimates + covariance/Wald and profile CIs.
- **With and without autocorrelation:** repeat the CIs under i.i.d.\ and under the
  Step-8 AR(1) correction applied to the $\kappa_p$ residuals (reuse
  `analysis/step8_autocorrelation/`).

## Analysis / comparison (the point)
- Do the estimates $(|\chi|,\dstar)$ move relative to the transformed-$\Delta c_{\mathrm m}$
  fit (§2.1–2.2)? By how much?
- Does the uncertainty change — is the raw-measurement formulation more/less
  precise, and how does AR(1) affect it here vs.\ in Step 8?
- Discuss the impact of removing the feature/response transformation: is the more
  statistically rigorous formulation materially different, and which should the
  group prefer?

## Figures (finalized `plot_style`, PNG+PDF)
`nacl_raw_fit.{png,pdf}` — measured vs.\ predicted $\kappa_p(t)$ (and/or vs.\
$\cint$), plus residuals; `nacl_raw_compare.{png,pdf}` — $(|\chi|,\dstar)$ point
estimates + CIs for transformed vs.\ raw-measurement, i.i.d.\ vs.\ AR(1).

## Report write-up (populate `sec:altform`)
Replace the "Planned" placeholder content of `sec:altform` with the actual
analysis: (i) the errors-in-variables motivation, (ii) the forward-model
mathematics above (calibration, flux balance / implicit $c_p$, predicted
$\kappa_p$, the measurement-error model and likelihood) with precision, (iii) the
results and comparison table (transformed vs.\ raw, i.i.d.\ vs.\ AR(1)), the
figures, and (iv) a verdict on whether the rigorous formulation changes the
conclusions. Keep the working assumptions explicitly flagged (this is an
exploratory reduced forward model; the fully rigorous version is a DATA~2.0-style
dynamic DAE on all raw channels — note as future work).

## Guardrails
- Edit `main.tex` only within `sec:altform`; do not touch other prose/results/model
  math, `refs.bib`, `data/`, `prior_analysis/`, or the §2.1–2.2 transformed-fit
  numbers.
- New code under `analysis/step9_raw_measurement/` (reuse preprocessing + Step-2/3/8
  imports; do not modify those files' numerics).
- Compile clean (0 undefined, 0 overfull); **rasterize + eyeball** the figures;
  `latexmk -c` after.

## Report back
The raw-measurement $(|\chi|,\dstar)$ ± CI (i.i.d.\ and AR(1)) vs.\ the transformed
fit; how much the estimates/uncertainty move; the measurement-error model used;
the forward-model working assumptions; the `sec:altform` write-up added; and
confirmation of a clean compile with the §2.1–2.2 numbers unchanged.
