# Prompt for Sonnet — Step 8 follow-up: AR(1) confidence regions computed BOTH ways, with the double-counting subtlety documented

You are revising the **`sec:naclAR1`** subsection (and its script/figures) in
`docs/reports/main.tex` (`/Users/adowling/DowlingLab/Membranes/data3_regression`,
env `data3-regression`). Documentation + calculation task — do not change the
model or the point estimates.

## The issue to resolve
The current `sec:naclAR1` computes the AR(1)-corrected confidence region on the
**whitened** SSE **and** sets the $F$-test denominator df to $n_{\mathrm{eff}}-p$.
Those are **two alternative ways to account for autocorrelation, not
complementary** — combining them double-counts the correlation:
- **GLS / whitening** already encodes the information loss in the whitened
  Jacobian $\tilde{\mathbf J}$; the correct region uses the \emph{actual} $n-p$ df
  (this is why the Wald result is $\sim\!2.7\times$).
- The **effective-sample-size ($n_{\mathrm{eff}}$) deflation** is the heuristic you
  use \emph{instead}, on the un-whitened OLS surface.

Combining them inflates the profile region by an extra $\sqrt{n/n_{\mathrm{eff}}}\approx2.8\times$,
which is almost certainly why the profile looks "$\sim\!8\times$ / much more
nonlinear." The PI wants **both consistent calculations shown on paper, the math
for each written out, and the discrepancy flagged for follow-up** — do not silently
pick one.

## Compute the region BOTH consistent ways
Report, for the primary NaCl dataset (corrected model), the 95% CIs three ways:
1. **i.i.d.** (baseline, already in §2.2).
2. **Method A — GLS / whitening (rigorous).** Prais–Winsten whitening +
   Cochrane–Orcutt (already implemented). Wald
   $\mathbf V=\hat\sigma_a^2(\tilde{\mathbf J}^\top\tilde{\mathbf J})^{-1}$,
   $\hat\sigma_a^2=\SSE_{\mathrm w}/(n-p)$; $F$-test/profile region on
   $\SSE_{\mathrm w}(\boldsymbol\theta)$ with **$n-p$** df. (This is the current
   code with the df changed from $n_{\mathrm{eff}}-p$ back to $n-p$.)
3. **Method B — effective-sample-size heuristic.** No whitening: keep the OLS fit
   and OLS SSE surface, but replace $n$ by $n_{\mathrm{eff}}=n\frac{1-\phi}{1+\phi}$
   — Wald SE $\times\sqrt{n/n_{\mathrm{eff}}}$; $F$-test/profile region on the
   \emph{OLS} $\SSE(\boldsymbol\theta)$ with **$n_{\mathrm{eff}}-p$** df.

Expectation to verify and state: Methods A and B should give **similar** widening
(both $\sim\!2.7$–$2.8\times$), and under each the profile should track its Wald
(mild nonlinearity, as in the i.i.d.\ case) — i.e.\ the region does **not** become
dramatically more nonlinear once the correction is applied \emph{consistently}.

## Report write-up (`sec:naclAR1`)
- Keep/So state the math already there (AR(1) model, $\hat\phi$, $n_{\mathrm{eff}}$,
  whitening, Cochrane–Orcutt). **Add the explicit math for Method B** (the
  $n_{\mathrm{eff}}$ deflation of the i.i.d.\ region) alongside Method A's GLS math,
  so both routes are written out.
- Replace the single CI table with a **three-column comparison** (i.i.d.\ | Method A
  GLS | Method B $n_{\mathrm{eff}}$), Wald and profile for $|\chi|$ and $\dstar$.
- **Correct the interpretation:** remove "the region grew substantially more
  nonlinear" / "$\sim\!8\times$ should be quoted going forward." State that both
  consistent methods give $\sim\!2.7\times$ and the profile tracks the Wald; the
  earlier $\sim\!8\times$ came from combining whitening with $n_{\mathrm{eff}}$ df.
- **Flag for the PI (a `\todo` / `\note`):** whitening (Method A) and
  $n_{\mathrm{eff}}$ deflation (Method B) are alternative corrections; combining
  them double-counts (the earlier draft did this). Open questions for the PI:
  which to adopt as the headline, and whether the $n_{\mathrm{eff}}$ heuristic is
  even appropriate given AR(1) is not fully sufficient (the whitened ACF still has
  structure, incl.\ a negative lag-1). Keep the AR(1)-insufficiency discussion and
  the "extend to CaCl2 / raw reformulation" future-work note.
- Update the region figure so it shows the i.i.d.\ vs.\ Method-A (GLS) regions
  (and, if clean, Method B) — with the corrected framing (no "8×"); keep the
  zoom+full two-panel style and the compliant plot style (boxed legends, bold
  math).

## Guardrails
- Edit only `sec:naclAR1` (+ its figures) in `main.tex`; touch nothing else
  (§2.2 pointer note may stay). Do not change point estimates, other sections,
  `refs.bib`, `data/`, `prior_analysis/`.
- `analysis/step8_autocorrelation/ar1_correction.py` gets Method B added and the
  Method A df fixed to $n-p$; regenerate figures; **rasterize + eyeball**.
- Compile clean (0 undefined, 0 overfull); `latexmk -c` after.

## Report back
The three-way CI table (i.i.d.\ / A / B) for $|\chi|$ and $\dstar$ (Wald + profile),
confirmation that A and B agree and that profile$\approx$Wald under each, the
corrected interpretation, the PI-follow-up flag added, and a clean-compile
confirmation with point estimates unchanged.
