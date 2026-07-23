# LaCl3 (3:1) Donnan regression + adsorption-hypothesis diagnostics (report §4)

Extends the CaCl2 (§3) adsorption-hypothesis analysis to \ce{LaCl3} (3:1),
the last of the single-salt trilogy. The 3:1 co-ion Donnan equation is a
**quartic with no compact closed form**, so this task's methodological point
is to fit it **two independent ways** and cross-check, mirroring Step 3's
scipy-vs-ParmEst validation:

- **Method A** (`regress_lacl3.py`): solve the quartic numerically per point
  via a bracketed root-find (Brent) and fit with
  `scipy.optimize.least_squares` — fast, drives the diagnostics.
- **Method B** (`parmest_lacl3.py`): pose the quartic as a **Pyomo nonlinear
  constraint** (the membrane co-ion concentration as a bounded decision
  variable) and estimate via ParmEst/Ipopt — the formulation that
  **generalizes to the multicomponent case** (no closed form exists there
  either for any salt combination).

## How to run

```bash
conda activate data3-regression
python analysis/lacl3/regress_lacl3.py      # Method A + all figures
python analysis/lacl3/parmest_lacl3.py       # Method B cross-check
```

Both take well under a minute. Reuses `analysis/step2_nacl/regress_nacl.py`
(`FitResult`) and `analysis/plot_style.py` unmodified; `parmest_lacl3.py`
imports `regress_lacl3.py` for the data and Method A comparison.

## Dataset

One processed \ce{LaCl3} run: `NF270_MC2 05.21.24_LaCl3` in
`data/BoE Analysis.xlsx` (also in `Rejection_Analysis.xlsx`), 623 valid rows,
standard 15-column schema (B=`c_int`, H=`c_p`, J=`J_w`). **Narrow
concentration range**, `c_int` ≈ 2–28 mM. A raw-only **MC4
`07.11.24_SLaCl3`** run exists in the raw campaign workbook but is **out of
scope** here (would need the preprocessing pipeline; noted as future work,
like the raw-only NaCl runs before they were processed).

**Anticipated headwinds** (both realized, see Results): (i) the narrow
2–28 mM range gives weak `|χ|` identifiability, the same range effect
flagged for NaCl's MC2 and CaCl2's MC3 (07.11.24, S); (ii) La³⁺ rejection is
low, so the fit may be ill-conditioned in the `|χ|` direction specifically.

## Model

**Co-ion = Cl⁻**, stoichiometry `c_co_s = 3*c_LaCl3` (3 Cl⁻ per LaCl3), for
both `c_int` and `c_p`.

**Donnan co-ion quartic** (report Eq. 31co), lumped constant `K` (≥0)
absorbing `δ*(δ°_co)²` and stoichiometric factors:
```
c_co_m^4 + |chi|*c_co_m^3 - K*c_co_s^4 = 0
```
Descartes' rule of signs (coefficient pattern `+,+,0,0,-` for `t>0`)
guarantees exactly one positive real root.

**Transformed response**: `deltaC_m = J_s_Cl*l/D_m`, `J_s_Cl = 3*c_p*J_w`
(Cl⁻ flux = 3× the LaCl3 molar flux), `l = 80e-9` m, 3:1 ambipolar membrane
diffusivity `D_m = 4*D_La_m*D_Cl_m/(3*D_La_m + D_Cl_m)`,
`D_La_m = 0.626e-12`, `D_Cl_m = 2.03e-12` m²/s.

**Prediction** = `c_co_m(c_int) - c_co_m(c_p)`. Fit `(|χ|, K)` **bounded ≥0**
(as for CaCl2 — `χ` enters the quartic unsquared/asymmetrically, so an
unconstrained fit can wander negative with no physical meaning).

### Root-solve implementation

A per-point bracketed Brent search (`scipy.optimize.brentq` with an
expanding-bracket search for the upper endpoint), verified against
`numpy.roots` on a 3000-point random parameter scan
(`chi∈[1e-3,1e3]`, `K∈[1e-6,1e2]`, `c_s∈[1e-2,1e3]`) to **<1e-6 relative
error**, and **~2.5× faster** than a per-point `numpy.roots` call at this
dataset's scale (`n=623` — small enough that either approach is fast; Brent
was chosen for consistency with Step 9's implicit-equation approach and for
headroom if this script is later extended to a larger LaCl3 dataset).

### Bug caught during development

The shared `profile_likelihood()` grid-bound floor (`max(..., 1e-6)`,
inherited from the CaCl2/NaCl scripts where `chi`/`K` are always `≳1e-4`)
**broke** for LaCl3's `K≈8.8e-8` — for a parameter this small,
`beta_hat + span_se*se` (the natural upper end of the profile grid) can
itself be *below* the 1e-6 floor, producing a grid with `lo > hi` (reversed
direction) and a nonsensical inverted CI (initially printed as
`(0.000001, 0.000000)`, i.e. lower bound above the "upper" bound — a red
flag that shouldn't have been ignored). Fixed by capping the floor at
`min(natural_hi * 0.5, 1e-10)` instead of a fixed `1e-6`, so it can never
exceed the grid's natural span. This is local to `analysis/lacl3/
regress_lacl3.py` and does not affect the NaCl/CaCl2 scripts (their `K`/`δ*`
values are always well above `1e-6`, so the bug never triggers there — but
worth being aware of if either is ever extended to a comparably deep,
weakly-identified regime).

## Results

```
n = 623
|chi| = 0.0000 mM  (95% CI 0.0000, 94.3717)   -- pinned at the lower bound
K     = 8.758e-08  (95% CI 8.298e-08, 9.236e-08)
SSE = 9.6770   RMSE/range = 0.1246   pseudo-R2 = 0.8028
Multi-start: 100.0% global (100 starts, method trf)
```

**Method A vs. Method B cross-check** (ParmEst warm-started at Method A's
converged root, for a fair/strong comparison):
```
Method A (scipy):   |chi| = 8.313e-14 mM, K = 8.75775e-08, SSE = 9.6770
Method B (ParmEst): |chi| = 1.130e-03  mM, K = 8.75802e-08, SSE = 9.6775
  Relative diff (K):   3.1e-05
  Relative diff (SSE): 5.1e-05
```
**Confirmed agreement.** `K` matches to 5 significant figures and SSE to
<0.01% between the two independent numerical routes (bracketed-Brent scipy
fit vs. Pyomo/Ipopt nonlinear-constraint NLP). `|χ|`'s raw relative
difference looks large (`~1e9`) only because BOTH methods independently
converge to a value indistinguishable from zero (`8×10⁻¹⁴` vs. `1×10⁻³` mM,
both far inside numerical noise for an unidentified boundary-pinned
parameter) — dividing two near-zero numbers is not a meaningful relative
comparison. The substantive result — **`|χ|` is not identifiable from this
dataset and both methods agree it is effectively zero** — is exactly the
cross-validation this task set out to obtain.

## Adsorption-hypothesis diagnostics

**Residuals vs. concentration**: a strong, smooth, systematic wave — dips to
about −0.19 mM around `c_int`≈10 mM, rises steadily back through zero by
`c_int`≈20 mM, and up to +0.15 mM at the top of the range. Sign-runs test:
**2 runs observed vs. ~293 expected under randomness, z=−24.9** — about as
extreme a non-random signature as could be observed on `n=623` residuals
(essentially two long monotonic-sign stretches, not the alternating pattern
white noise would produce).

**Effective `|χ|` vs. concentration**: of 19 sliding windows, only 7 (all at
the high-concentration end, `c_int`≳21 mM) are reliably identified; the
low/mid-range windows are all unidentified (pinned at/near the search-box
bound) and excluded from the trend. Among the reliable windows, local `|χ|`
rises monotonically from ~0.06 to ~0.85 mM as `c_int` increases — i.e. even
in the one part of the range where a local estimate is possible at all,
`|χ|` stays two orders of magnitude below NaCl's ~18–50 mM and CaCl2's
already-small values, and the trend is monotonically *increasing* with
concentration (unlike CaCl2's non-monotonic inverted-U).

**Comparison across salts**: `|χ|→0` was already CaCl2's headline finding;
LaCl3 pushes this further — the global point estimate is pinned exactly at
the `|χ|≥0` bound (not just small), and only a fraction (0.06–0.85 mM) of
even the local windowed estimates ever exceeds numerical noise. The
`pseudo-R²=0.80` is noticeably worse than any CaCl2 dataset (0.84–0.95),
consistent with a poorer-conditioned, more strongly adsorption-affected 3:1
system, as anticipated. **Verdict: the adsorption hypothesis extends to
La³⁺, and more severely than for Ca²⁺** — consistent with the trend NaCl
(well-identified, no adsorption signature) → CaCl2 (`|χ|→0`, moderate
residual structure) → LaCl3 (`|χ|` pinned at exactly 0, strongest residual
structure of the three, worst fit quality) predicted by increasing cation
charge/adsorption affinity. This conclusion is honestly caveated by the
narrow concentration range's weak identifiability — a wider-range LaCl3
dataset (or the raw-only MC4 run, once preprocessed) would strengthen it.

## Figures (in `docs/reports/figures/`, `plot_style.py` conventions)

- `lacl3_fit.png` — half-width; fit vs. data.
- `lacl3_residuals_vs_conc.png` — half-width; the key diagnostic.
- `lacl3_effective_chi_vs_conc.png` — half-width; local `|χ|` vs. `c_int`
  (unreliable windows shown as faint `×`).
- `lacl3_appendix.png` — full-width, compact 4-panel (fit, residual,
  effective-χ, parameter summary) — the single-dataset analog of the CaCl2
  per-dataset appendix.

## Guardrails honored

- `data/`, `prior_analysis/`, `refs.bib` untouched.
- `analysis/step2_nacl/regress_nacl.py` and `analysis/plot_style.py` imported,
  not modified. NaCl/CaCl2 numbers unaffected (verified: no shared code
  changed; `analysis/cacl2/` and `analysis/nacl_all_datasets/` untouched).
- `main.tex` edited only within `sec:lacl3` and the new LaCl3 appendix
  subsection.
