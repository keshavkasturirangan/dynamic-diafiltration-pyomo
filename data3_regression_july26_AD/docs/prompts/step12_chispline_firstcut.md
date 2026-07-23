# Prompt for Sonnet — Step 12 (first cut): a SIMPLE concentration-dependent χ(c) spline regression on CaCl2

You are attempting the **first, deliberately simple** version of a
concentration-dependent membrane-charge model: replace the scalar $|\chi|$ with a
smooth $\chi(c)$ represented by a B-spline, regressed **jointly** with the lumped
constant $K$, on the \ce{CaCl2} data. The goal of this first cut is narrow: **find
out whether the apparent adsorption drift is a real, identifiable $\chi(c)$ or an
artifact of the constant-charge closure** — not to build the full penalized/UQ
machinery yet. Work in `/Users/adowling/DowlingLab/Membranes/data3_regression`, env
`data3-regression` (`/opt/anaconda3/envs/data3-regression/bin/python`).

## Gating / dependencies (do not start until these hold)
- **Run only after CP Phase 2 (`cp_phase2_regressions.md`) is committed and audited**
  (and Phase 3, if we do one). This analysis **must use the CP-corrected
  \ce{CaCl2} $\cint$** — a flexible spline will otherwise absorb CP/calibration
  misspecification and return a meaningless $\chi(c)$ (see the reasoning in the
  drafted methodology, next bullet). Use whatever CP-corrected data path Phase 2
  established in `analysis/cacl2/regress_cacl2.py` (the authors'-bulk+CP `c_int`).
- **The mathematics is already written up** in
  `docs/reports/draft_chispline_futurework.tex` (a standalone future-work section
  drafted by the coordinator). **Read it first** — this prompt implements a
  simplified subset of §"Estimation" / §"Proposed first study" there. Do not
  re-derive; follow that formulation.

## Scope of the first cut (keep it simple)
- **One salt, one wide-range dataset:** \ce{CaCl2} **MC5 05.27.26** (the widest range,
  $\cint$ to $\sim\SI{108}{mM}$, $n\approx1645$) — the range that makes $\chi(c)$
  estimable. (Pooling all three CaCl2 datasets is a later refinement; not now.)
- **scipy only** (`scipy.optimize.least_squares`), reusing the existing CaCl2 co-ion
  **cubic root solver** in `analysis/cacl2/regress_cacl2.py`. (The Pyomo
  nonlinear-constraint version is the eventual generalizable path — note it as future
  work, don't build it here.)
- **No roughness penalty yet** — an unpenalized few-knot regression spline is the
  simplest thing that can work. (Mention the P-spline penalty as the obvious
  robustness add-on.)

## Model
Replace the scalar $|\chi|$ in the CaCl2 co-ion cubic with
$\chi(c)=\sum_{k}\beta_k B_k(c)$, a **cubic B-spline** in $\cint$, evaluated **at each
face** (feed and permeate): the prediction is
$\cm{\co}(\cint;\chi(\cint),K)-\cm{\co}(\cp;\chi(\cp),K)$, exactly as in
Eqs. (2)–(3) of the drafted section. The response $\Delta c_m=J_{s}\ell/D_m$ is
**unchanged** (uses $\cp$, $\Jw$), so only the predictor shifts.
- **Basis:** build the B-spline basis with `scipy.interpolate.BSpline` /
  `splinelab`-style knots, or a small explicit design matrix. Interior knots at
  **quantiles of the observed $\cint$**; clamp/boundary knots at the data range
  endpoints (no extrapolation). Try **1 and 2 interior knots** (i.e. two model
  complexities), plus the constant-$\chi$ case (0 interior knots).
- **Positivity, simply:** cubic B-spline basis functions are non-negative, so
  **constrain all control points $\beta_k\ge0$** (`least_squares` bounds) — that alone
  guarantees $\chi(c)\ge0$ without an $\exp$ link or extra constraints.
- **Parameters:** $(\boldsymbol\beta, K)$, jointly. **Initialize** $\beta_k$ all equal
  to the current constant-$\chi$ estimate (flat spline) and $K$ at its constant-fit
  value — a warm, physically sensible start.

## Analysis / deliverables
1. **Fit** the 0- (constant), 1-, and 2-interior-knot models; report for each:
   $\hat K$, the $\beta_k$, SSE, a scale-free fit metric (pseudo-$R^2$), and the number
   of free parameters. Confirm multi-start/convergence robustness on the flexible fits.
2. **Nested model comparison:** constant-$\chi$ vs 1-knot vs 2-knot by an **$F$-test**
   (extra-sum-of-squares) **and AIC/BIC**. State plainly whether added flexibility is
   justified by the data or is fitting noise.
3. **The key cross-check (cheap and decisive):** overlay the regressed $\chi(\cint)$
   curve on the **model-free sliding-window "effective $|\chi|$" diagnostic** already
   computed for this dataset (`sec:cacl2diag` / the CaCl2 code). If the regressed
   spline tracks the independent windowed estimates, that is strong corroboration that
   $\chi(c)$ is real; if not, be honest about the discrepancy.
4. **Pointwise band (simple version):** from the Gauss–Newton covariance
   $\hat\sigma^2(\mathbf J^\top\mathbf J)^{-1}$, plot $\hat\chi(\cint)$ with a pointwise
   $\pm z\,\mathrm{s.e.}$ band, $\mathrm{Var}(\hat\chi(c))=\mathbf B(c)^\top
   \mathrm{Cov}(\hat{\boldsymbol\beta})\mathbf B(c)$ (Eq. (6) of the drafted section).
5. **Optional null-check (recommended, cheap):** run the same 1-knot spline on a
   well-identified **NaCl** dataset (MC5 05.27.26). Expect the spline to collapse to
   ≈constant ($\beta_k$ ≈ equal, no AIC improvement) — a control that the method does
   **not** manufacture spurious drift. Report the outcome briefly.

## Figure (plot_style compliant; PNG+PDF; rasterize-verify)
`docs/reports/figures/cacl2_chispline_firstcut.{png,pdf}`: (a) fit vs data (constant
vs spline); (b) $\hat\chi(\cint)$ with pointwise band, overlaid on the constant-$\chi$
line and the sliding-window effective-$\chi$ points.

## Report / integration
- **Do NOT edit `main.tex` in this pass.** Deliver the analysis + figure + report-back
  numbers; the coordinator will integrate the drafted future-work section *and* this
  first result together after reviewing (and will decide whether it stays "future
  work" or becomes "preliminary results").
- New code under **`analysis/chi_spline/`** (reuse the CaCl2 cubic root solver +
  `plot_style`; import the CP path from the Phase-2 CaCl2 driver). Do not change the
  existing CaCl2 constant-χ code or its committed numbers.

## Guardrails
- Use CP-corrected $\cint$ (Phase 2). Do not change the model math (the cubic, the
  transformed-response definition) — only $|\chi|\to\chi(\cint)$.
- Don't touch `data/`, `refs.bib`, `prior_analysis/`, or other salts' code.
- Keep it simple: scipy, unpenalized, ≤2 interior knots, $\beta_k\ge0$. Flag anything
  that looks non-identifiable (e.g. wild $\beta_k$, exploding CI band) rather than
  forcing a fit.

## Report back
- Table: constant vs 1-knot vs 2-knot — $\hat K$, SSE, pseudo-$R^2$, #params, and the
  $F$-test / AIC / BIC comparison; the verdict on whether $\chi(c)$ drift is real and
  identifiable.
- Whether the regressed $\chi(\cint)$ tracks the independent sliding-window
  effective-$\chi$ (the decisive cross-check).
- The NaCl null-check outcome (did the spline stay flat?).
- The figure; any identifiability red flags; and a recommendation on whether the full
  treatment (penalty, pooling, Pyomo/ParmEst, UQ, multi-salt) is worth pursuing.
