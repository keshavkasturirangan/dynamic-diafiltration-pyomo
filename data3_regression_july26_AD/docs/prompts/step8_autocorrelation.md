# Prompt for Sonnet — Step 8: autocorrelation-corrected confidence regions (AR(1)) for NaCl

You are implementing **Step 8**: correcting the NaCl confidence regions for the
autocorrelation that the i.i.d.-error assumption ignores. Work in
`/Users/adowling/DowlingLab/Membranes/data3_regression`, env `data3-regression`.
Add a **new subsection to the NaCl section** of `docs/reports/main.tex` that
**explains the mathematics inline** (the PI wants the math in this exploratory
results section; it can move to Methods later) and reports the results.

> **Run after** the figure-compliance audit (which finalizes `analysis/plot_style.py`
> — boxed legends, bold math). Use that finalized style for all figures.

## Motivation (already flagged in the report)
The regressions minimize $\SSE(\theta)=\sum_i r_i^2$, $r_i=y_i-f(\theta;x_i)$ with
$y_i=\Delta c_{\mathrm m}$, assuming the $r_i$ are i.i.d.\ $N(0,\sigma^2)$. But the
$\sim$692 points are a dense single-run **time series**, so the residuals are
strongly autocorrelated and every reported region (Wald, F-test, profile,
bootstrap) is **too narrow**. Quantify how much, and check whether the corrected
region is **more nonlinear**.

## Dataset & reuse
Primary NaCl dataset `NF270_MC5 05.27.26_NaCl (2)` (the one in §2.1–2.2),
corrected model $f_{\text{corr}}$. Reuse `analysis/step2_nacl/regress_nacl.py`
(data, model, fit) and `analysis/step3_nacl_uncertainty/scipy_uncertainty.py`
(F-test threshold, SSE surface, profile) — the AR(1) layer wraps these.

## Mathematics to implement AND write up (with precision)
1. **AR(1) residual model.** Order the OLS residuals by time and model
   $r_i=\phi\,r_{i-1}+a_i$, $a_i\sim\text{iid }N(0,\sigma_a^2)$, $|\phi|<1$.
   Estimate the lag-1 autocorrelation
   $\hat\phi=\dfrac{\sum_{i=2}^n r_i r_{i-1}}{\sum_{i=1}^n r_i^2}$.
2. **Effective sample size** $n_{\mathrm{eff}}=n\dfrac{1-\phi}{1+\phi}$, and the
   naive CI-inflation factor $\sqrt{n/n_{\mathrm{eff}}}=\sqrt{(1+\phi)/(1-\phi)}$.
3. **GLS via Prais–Winsten whitening.** The AR(1) covariance is
   $\Sigma_{ij}=\frac{\sigma_a^2}{1-\phi^2}\phi^{|i-j|}$. Whiten:
   $\tilde y_1=\sqrt{1-\phi^2}\,y_1$, $\tilde y_i=y_i-\phi y_{i-1}$ ($i\ge2$), and
   the model $\tilde f$ likewise; refit $\theta$ on the whitened residuals
   (Cochrane–Orcutt: iterate $\hat\theta\leftrightarrow\hat\phi$ to convergence, or
   joint MLE). Report whether the **point estimate moves** vs.\ OLS.
4. **Corrected uncertainty.** Wald covariance from the whitened Jacobian
   $\mathbf V=\hat\sigma_a^2(\tilde{\mathbf J}^\top\tilde{\mathbf J})^{-1}$; the
   **nonlinear F-test / profile region on the whitened (GLS) objective** using
   $n_{\mathrm{eff}}$ (or the AR(1) log-likelihood). These are the
   correctly-sized regions.
5. **Sufficiency of AR(1).** Compute the ACF of the OLS residuals and of the
   \emph{whitened} residuals; if whitening removes the correlation (whitened ACF
   $\approx0$ beyond lag 0) AR(1) suffices; if structure remains, note AR($p$)/ARMA
   as future work.

## Analysis
- Report $\hat\phi$, $n_{\mathrm{eff}}$ and $n_{\mathrm{eff}}/n$, and the inflation
  factor.
- Table: $(|\chi|,\dstar)$ and their 95% CIs **i.i.d.\ vs.\ AR(1)-corrected** (Wald
  and profile); state how many-fold wider the corrected CIs are.
- Compare region **shape** i.i.d.\ vs.\ corrected (overlay; note any added
  curvature/asymmetry).
- Whitened-residual ACF to judge AR(1) sufficiency.

## Figures (finalized `plot_style`)
Into `docs/reports/figures/` (PNG+PDF): `nacl_ar1_acf.{png,pdf}` (residual ACF,
raw vs.\ whitened, with the $\pm1.96/\sqrt n$ band) and
`nacl_ar1_regions.{png,pdf}` (i.i.d.\ vs.\ AR(1)-corrected confidence regions
overlaid on the $(|\chi|,\dstar)$ plane).

## Report write-up
Add a subsection **after** `sec:naclUQ` (e.g.\ `\subsection{Autocorrelation-corrected
uncertainty}\label{sec:naclAR1}`) that (i) states the math of items 1–5 with
precision, (ii) reports the results/table/figures, and (iii) gives a verdict:
how badly the i.i.d.\ regions understated uncertainty, whether the region became
more nonlinear, and whether AR(1) is adequate or richer models are future work.
Also update the existing "Planned correction" note in `sec:naclUQ` to point to this
new subsection (it is no longer just planned).

## Guardrails
- Edit `main.tex` only by adding `sec:naclAR1` and updating that one pointer note
  in `sec:naclUQ`; do not touch other prose/results/model math, `refs.bib`,
  `data/`, `prior_analysis/`, or the NaCl point estimates.
- New code under `analysis/step8_autocorrelation/` (reuse Step-2/3 imports).
- Compile clean (`latexmk -pdf`: 0 undefined, 0 overfull); **rasterize + eyeball**
  the two figures; `latexmk -c` after.

## Report back
$\hat\phi$, $n_{\mathrm{eff}}/n$, the i.i.d.-vs-corrected CI table and fold-widening,
whether the region grew more nonlinear, the AR(1)-sufficiency verdict, the
subsection added, and confirmation of a clean compile + unchanged point estimates.
