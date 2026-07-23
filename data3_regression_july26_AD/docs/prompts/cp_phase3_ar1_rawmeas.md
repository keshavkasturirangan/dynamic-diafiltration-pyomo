# Prompt for Sonnet — CP Phase 3: re-run Step 8 (AR(1)) and Step 9 (raw-measurement) concentration-polarization–corrected

You are finishing the CP redo. **Phase 2** made CP the standing convention for the
single-salt point-estimate/uncertainty regressions. **This is Phase 3**: re-run the
two NaCl uncertainty analyses that build on the primary fit —
**Step 8 (AR(1) autocorrelation, `sec:naclAR1`)** and
**Step 9 (raw-measurement reformulation, `sec:altform`)** — on the **CP-corrected**
NaCl fit, and update those two sections. Work in
`/Users/adowling/DowlingLab/Membranes/data3_regression`, env `data3-regression`.

## Gating
- **Run after the Phase 2 refinement (`cp_phase2_refine.md`) is committed.** (Both
  build on the step2 primary NaCl fit, which the refinement does not change — the
  refinement only touches the all-NaCl secondary table — but run after it to avoid a
  concurrent `main.tex` edit.)
- The CP-corrected primary NaCl fit is the (2)-sheet bulk + our CP: **$|\chi|=56.74$
  mM** (was 50.36 non-CP), $\dstar$ correspondingly shifted. Both Step 8 and Step 9
  must be recomputed against this fit.
- Keep both analyses' **non-CP results reproducible** for the CP-vs-non-CP comparison
  (the committed report currently carries the non-CP numbers).

## Part A — Step 8: AR(1) confidence regions, CP-corrected (`sec:naclAR1`)
`analysis/step8_autocorrelation/ar1_correction.py` reuses `regress_nacl as step2`
(`step2.f_corrected`, `step2.load_sheet`), so it inherits the `APPLY_CP` toggle.
- Re-run the full AR(1) machinery on the **CP-corrected** primary fit: $\hat\phi$,
  $n_{\mathrm{eff}}=n\frac{1-\phi}{1+\phi}$, **Method A** (Prais–Winsten/GLS whitening,
  Wald + $F$/profile region with the actual $n-p$ df) and **Method B** ($n_{\mathrm{eff}}$
  heuristic on the OLS surface). Report the widening factors **CP vs non-CP**.
- **Preserve the Step-8b framing and do NOT reintroduce the double-count:** whitening
  (Method A) and the $n_{\mathrm{eff}}$ deflation (Method B) are *alternative*
  corrections; the region uses $n-p$ df on the whitened SSE. Both should again give a
  similar $\sim$2.7–2.9$\times$ widening and profile$\approx$Wald. If CP materially
  changes $\hat\phi$ / $n_{\mathrm{eff}}$ / the widening, report it; otherwise state it
  is essentially unchanged (autocorrelation is a property of the dense time series, not
  of the CP shift).
- Update `sec:naclAR1`'s three-way CI table and numbers to CP-corrected (non-CP in a
  comparison column or a sentence), the region figure, and **remove the Phase-3-pending
  one-line note** Phase 2 added.

## Part B — Step 9: raw-measurement reformulation, CP-corrected (`sec:altform`)
`analysis/step9_raw_measurement/raw_regression.py` predicts the measured $\kappa_p(t)$
from the retentate conductivity $\kappa_r(t)$: bulk feed $c_f=(\kappa_r-a)/s$ (embedded
NaCl calibration $a=-63.706$, $s=76.685$), $\Jw$ from the smoothed mass slope, and
currently sets the **feed interface $c_{\mathrm{int}}(t)=c_f(t)$ directly (non-CP)**.
- **Insert the thin-film CP** into the $c_{\mathrm{int}}$ reconstruction, consistent
  with Phase 2: $c_{\mathrm{int}}=c_p+(c_f-c_p)\exp(\Jw/k_{\ce{NaCl}})$, permeate
  unpolarized, $k$ imported from `cp_authors_bulk`/Phase 1 (single source of truth).
- **Handle the coupling carefully (the real subtlety):** here $c_p$ is *model-predicted*
  (from the Donnan flux balance $\Jw c_p=(\Dm/\ell)[\cm{\co}(c_{\mathrm{int}})-\cm{\co}(c_p)]$),
  yet $c_p$ also appears in the CP formula for $c_{\mathrm{int}}$. So under CP,
  $(c_{\mathrm{int}}, c_p)$ must be solved **jointly** per time point (a coupled
  2-equation root-find / fixed-point), not the current single brentq for $c_p$ at fixed
  $c_{\mathrm{int}}=c_f$. Implement the coupled solve, or if you use an approximation
  (e.g. seeding $c_p$ from the measured $\kappa_p$), **document it explicitly** and note
  its effect. Then $\hat\kappa_p=a+s\,c_p$, residual $=\kappa_{p,\mathrm{meas}}-\hat\kappa_p$
  as before.
- Re-fit; report the CP-corrected $(|\chi|,\dstar)$ **vs** the non-CP raw fit
  ($|\chi|\approx115$, $\dstar\approx1.01$) **and vs** the transformed CP fit
  ($|\chi|=56.74$). The transformed-vs-raw gap under CP is the errors-in-variables
  sensitivity — quantify it.
- **Keep the honest, exploratory framing:** the near-unit-root AR breakdown
  ($\hat\phi\to0.999$, $n_{\mathrm{eff}}<1$) is a non-stationarity issue that CP does
  **not** fix — report whether $\hat\phi$ changed, and keep the "exploratory, not a
  replacement; needs a non-stationary/state-space noise model" verdict. Update
  `sec:altform`'s numbers/figure and **remove the Phase-3-pending note**.

## Report updates (`docs/reports/main.tex`)
- `sec:naclAR1` and `sec:altform`: CP-corrected primary, non-CP retained for comparison,
  pending-notes removed.
- **Discussion points 13 (AR(1) convention) and 14 (raw-measurement)** in
  `sec:questions`: update the quoted numbers (the $\sim$2.7–2.9$\times$ widening; the
  transformed$\to$raw shift, previously "$+128\%$" from $50.36\to115$ — recompute under
  CP from $56.74\to$ the new raw value).
- **Abstract:** update the "i.i.d.\ CIs $\sim$2.7–2.9$\times$ too narrow" and the
  "$+128\%$" raw-measurement figures if CP moves them.

## Guardrails
- Scope: `step8_autocorrelation/`, `step9_raw_measurement/`, and the `sec:naclAR1` /
  `sec:altform` sections (+ the abstract / discussion-points 13–14 numbers they feed).
  **Do not** change the CP mechanism, other salts, the model math, `refs.bib`, `data/`,
  or `prior_analysis/`.
- Keep non-CP reproducible for both (regression test); import $k$ from Phase 1.
- Compile clean (`latexmk -pdf`: 0 undefined, 0 overfull; `latexmk -c`); report page
  count; rasterize + eyeball the two regenerated figures.

## Report back
- Step 8: $\hat\phi$, $n_{\mathrm{eff}}$, and Method A / Method B widening — CP vs
  non-CP; confirmation the double-count is not reintroduced and profile$\approx$Wald.
- Step 9: CP-corrected $(|\chi|,\dstar)$ vs non-CP raw vs transformed-CP; how the
  $(c_{\mathrm{int}},c_p)$ coupling was solved; whether the near-unit-root behavior
  persists.
- Sections/figures/abstract/discussion-points updated; clean-compile confirmation
  (0 undefined, 0 overfull, page count).
