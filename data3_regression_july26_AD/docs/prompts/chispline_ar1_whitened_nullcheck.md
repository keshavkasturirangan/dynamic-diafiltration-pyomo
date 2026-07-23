# Prompt for Sonnet — AR(1)-whitened χ(c) null-check (decide if the CaCl2 drift is real or absorbed autocorrelation)

You are running the **cheap, decisive follow-up** to the χ(c)-spline first cut
(§\ref{sec:chispline}, `analysis/chi_spline/chispline_cacl2.py`). Work in
`/Users/adowling/DowlingLab/Membranes/data3_regression`, env `data3-regression`.

## The question
The first cut found that a $\chi(c)$ spline dramatically beats constant-$\chi$ on
\ce{CaCl2} (SSE $418\to52.6$, $F(4,1639)=2851$) --- **but the \ce{NaCl} null check
failed** (the same machinery finds significant, non-flat $\chi(c)$ on the
well-identified $1{:}1$ case that should have no adsorption drift). The suspected
cause is a **time--concentration confound within a single run**: under
concentrating-lag diafiltration $c$ rises monotonically with time, so
$\chi(c)\approx\chi(t)$ and an unpenalized spline can absorb \ce{NaCl}'s AR(1)
temporal drift ($\hat\phi\approx0.74$, \S\ref{sec:naclAR1}) rather than chemistry.
The i.i.d.\ $F$-test/AIC used in the first cut do not know about that
autocorrelation and are anti-conservative.

**The test:** redo the model comparison on **AR(1)-whitened residuals** (which
remove the lag-1 temporal correlation), for **both** datasets. Then:
- If the spurious \ce{NaCl} "drift" **evaporates** once whitened (constant-vs-spline
  $F$ no longer significant, $\hat\chi(c)$ flattens) **and \ce{CaCl2} survives**
  (still strongly favors the spline) $\Rightarrow$ **green light**: the \ce{CaCl2}
  drift is real concentration structure, autocorrelation was the confound.
- If the \ce{NaCl} improvement **persists** even whitened $\Rightarrow$ the
  time--concentration confound is not just AR(1); it must be broken by **pooling
  across differently-paced runs** (time and $c$ decorrelate *across* runs).

## Method (reuse existing machinery --- this is meant to be small)
Combine Step 8's whitening (`analysis/step8_autocorrelation/ar1_correction.py`:
Prais--Winsten/Cochrane--Orcutt $\hat\phi$ estimation and the
`whitened_resid_fn`/whitening operator) with Step 12's spline
(`analysis/chi_spline/chispline_cacl2.py`: `f_chispline`, `f_nacl_chispline`,
`fit_chispline`, the nested $F$-test/AIC, the constant baselines). Put the new code
in `analysis/chi_spline/` (a new script, e.g.\ `chispline_ar1_whitened.py`, importing
both --- do not modify the Step 8 or first-cut scripts).

For each dataset (\ce{NaCl} MC5 05.27.26 null check; \ce{CaCl2} MC5 05.27.26),
CP-corrected, matching the first cut's data exactly:
1. Fit **constant-$\chi$** (i.i.d.) as the baseline; from its time-ordered residuals
   estimate $\hat\phi$ (report it --- for \ce{CaCl2} this is a **new** number; Step 8
   only characterized \ce{NaCl}).
2. **Whiten** the residual sequence with the whitening operator $\mathbf W(\hat\phi)$
   (Prais--Winsten: $\tilde r_1=\sqrt{1-\hat\phi^2}\,r_1$,
   $\tilde r_i=r_i-\hat\phi r_{i-1}$). Refit **both** the constant model and the
   **1-knot** $\chi(c)$ spline by minimizing the **whitened** SSE
   $\lVert\mathbf W(\mathbf y-\hat{\mathbf y}(\boldsymbol\theta))\rVert^2$ (reuse Step 8's
   whitened-residual function, generalized to take the spline model's prediction).
3. **Nested $F$-test on the whitened SSE**, constant vs.\ 1-knot spline, with the
   actual $n-p$ df (Step 8's **Method A**, the GLS-consistent choice). Report the
   whitened SSE, $F$, $p$, $\Delta$AIC, and the fitted $\hat\chi(c)$ control-point
   spread (flat or not?).
4. **Fair-comparison note:** whiten both models with a **common $\hat\phi$** (from the
   constant baseline) so the $F$-test isolates the mean-model change under one noise
   model; also report the sensitivity of re-estimating $\hat\phi$ per model
   (Cochrane--Orcutt), as Step 8b did with its Method A/B. Keep the framing honest and
   do **not** reintroduce Step 8's double-count (whitening OR $n_{\mathrm{eff}}$, not
   both).

## Deliverables
- A small table: for \ce{NaCl} and \ce{CaCl2}, i.i.d.\ vs.\ whitened --- constant-vs-1-knot
  $F$, $p$, $\Delta$AIC, and $\hat\chi(c)$ spread; plus each dataset's $\hat\phi$.
- The **verdict** per the two branches above, stated plainly.
- A compact figure (`docs/reports/figures/cacl2_chispline_whitened.{png,pdf}`,
  plot_style-compliant): the whitened $\hat\chi(c)$ for \ce{NaCl} (does it flatten?)
  and \ce{CaCl2} (does it survive?), ideally beside the i.i.d.\ curves for contrast.
- **Report update (scoped):** update §\ref{sec:chispline}'s "Disentangle before
  interpreting" future-work bullet and the verdict to report the *outcome* of this
  test (it currently *recommends* the test as the next step). Keep the section's
  honest, preliminary framing; if the test is decisive in either direction, say so
  plainly. Flag the change for coordinator (Opus) review rather than rewriting the
  section's tone unilaterally.

## Guardrails
- Reuse Step 8 + Step 12 code; do not modify those scripts or any other analysis,
  `data/`, `refs.bib`, or other report sections. New code confined to
  `analysis/chi_spline/`.
- Do not change the constant-$\chi$, CP, or first-cut numbers already in the report ---
  this adds a whitened comparison alongside them.
- Compile clean if you edit the report (`latexmk -pdf`: 0 undefined, 0 overfull;
  `latexmk -c`); rasterize + eyeball the new figure.

## Report back
- The i.i.d.-vs-whitened comparison table (both datasets), each $\hat\phi$, and the
  common-$\phi$ vs.\ per-model-$\phi$ sensitivity.
- Whether the \ce{NaCl} spurious drift evaporates under whitening and whether \ce{CaCl2}
  survives --- and therefore which branch (green light / confound-persists) holds.
- The §\ref{sec:chispline} update made; clean-compile confirmation.
