# NaCl regression across all available datasets (report §2.4)

Fits the NaCl (1:1) Donnan model to all three processed NaCl datasets and
assesses reproducibility of `|χ|` (fixed charge, mM) and `δ*` (partition
factor) across membrane cuts (MC2 vs MC5) and a repeat run on the same cut
(MC5 `NaCl` vs `NaCl (2)`).

## How to run

```bash
conda activate data3-regression
python analysis/nacl_all_datasets/regress_all_nacl.py
```

Takes a few seconds. Reuses, unmodified: `analysis/step2_nacl/regress_nacl.py`
(`build_deltaC_m`, `f_matlab`, `f_corrected`, `fit_model`, `load_sheet`) and
`analysis/step3_nacl_uncertainty/scipy_uncertainty.py` (`f_test_threshold`,
`sse_surface`, `profile_likelihood`, `profile_ci_from_curve`). No
model/fitting logic is duplicated in this script.

## Datasets

All three sheets live in `data/Rejection_Analysis.xlsx` (identically in
`data/BoE Analysis.xlsx`), same 15-column schema (B=`c_int`, H=`c_p`,
J=`J_w`):

| Membrane | Sheet | Valid rows |
|---|---|---|
| MC2 | `NF270_MC2 05.07.24_NaCl` | 482 |
| MC5 | `NF270_MC5 05.27.26_NaCl` | 747 |
| MC5 | `NF270_MC5 05.27.26_NaCl (2)` | 747 (Step 2/3 primary dataset) |

The 6 raw-only NaCl runs in `data/NF270_MC{2,3,4,5}.xlsx` are out of scope
here (need interfacial-concentration/`J_w` preprocessing first) — noted as
future work.

## Windowing rule

Unlike Bill's ad hoc row window (57–748) for the MC5(2) sheet, MC2 and plain
MC5 have no equivalent hand-picked window, so for comparability **the full
valid data range is the primary fit** for all three datasets. A
**startup-trimmed** fit (drop the first `STARTUP_TRIM_FRAC = 7%` of rows —
the low-concentration startup transient) is reported as a sensitivity check.

## Results (primary = `f_corrected`, full range, 95% CIs = profile likelihood)

| Dataset | n | `\|χ\|` [mM] | 95% CI | `δ*` [–] | 95% CI | SSE | Multi-start % global |
|---|---|---|---|---|---|---|---|
| MC2 (05.07.24) | 482 | 18.46 | (18.06, 18.87) | 0.2888 | (0.2857, 0.2920) | 4.47 | 100% |
| MC5 (05.27.26) | 747 | 49.56 | (48.58, 50.55) | 0.3251 | (0.3232, 0.3270) | 130.07 | 100% |
| MC5 (05.27.26) (2) | 747 | 49.55 | (48.58, 50.54) | 0.3251 | (0.3232, 0.3270) | 129.62 | 100% |

Multi-start: 100 random starts per dataset, `|χ|₀ ∈ [1,200]` mM, `δ*₀`
log-uniform on `[10⁻³,10]`, `scipy.optimize.least_squares` method `lm`
(Step 3 found `lm`/`trf` reliably reach the global optimum here, unlike
`dogbox`, which can get trapped at `|χ|→0`). All three datasets converged to
the global optimum from **100% of 100 starts** — the same reliable behavior
found for the Step-3 primary dataset extends to MC2 and plain MC5.

Secondary (uncorrected `f_matlab`, for continuity with Bill's original
values, full range):

| Dataset | `\|χ\|` [mM] | `δ*` [–] | SSE |
|---|---|---|---|
| MC2 (05.07.24) | 9.228 ± 0.080 | 0.07220 ± 0.00031 | 4.47 |
| MC5 (05.27.26) | 24.780 ± 0.194 | 0.08127 ± 0.00019 | 130.07 |
| MC5 (05.27.26) (2) | 24.777 ± 0.194 | 0.08127 ± 0.00019 | 129.62 |

(As in Step 2, corrected = 2×|χ|, 4×δ* of uncorrected, exactly, for every
dataset — the factor-of-2 identity is model-structural, not
dataset-specific.)

### Startup-trimmed sensitivity (drop first 7% of rows, corrected model)

| Dataset | n dropped | `\|χ\|` full → trimmed | `δ*` full → trimmed |
|---|---|---|---|
| MC2 (05.07.24) | 34 | 18.46 → 18.50 (+0.25%) | 0.2888 → 0.2891 (+0.11%) |
| MC5 (05.27.26) | 52 | 49.56 → 50.21 (+1.32%) | 0.3251 → 0.3261 (+0.32%) |
| MC5 (05.27.26) (2) | 52 | 49.55 → 50.21 (+1.33%) | 0.3251 → 0.3261 (+0.32%) |

Estimates move by ≤1.3% under trimming for all three datasets — small
relative to the cross-dataset differences below, so the startup transient is
not the cause of the MC2-vs-MC5 disagreement.

## Reproducibility verdict

**The MC5 repeat is highly reproducible**: `|χ|` (49.56 vs 49.55 mM) and `δ*`
(0.3251 vs 0.3251) agree to within a few hundredths of a percent, with fully
overlapping 95% CIs — re-running the same membrane cut gives essentially
identical parameters. **MC2 does *not* agree with MC5**: `|χ|` is roughly
2.7× smaller (18.46 vs 49.56 mM) and `δ*` is measurably lower (0.289 vs
0.325), with non-overlapping 95% CIs — a real, not noise-level, discrepancy
across membrane cuts. Neither cut's 95% CI contains the independently-used
BoE reference of `|χ| ≈ 44 mM`; MC5's point estimate (≈49.5 mM) is
reasonably close (~13% high), while MC2's (≈18.5 mM) is off by more than
2×. Because MC2's data span only `c_int ∈ [1.6, 44.4]` mM versus MC5's
`[0.8, 182]` mM (a 4× narrower range with no data past the point where MC5's
curve is most informative about `|χ|` and `δ*` jointly), we cannot yet tell
whether this reflects genuine cut-to-cut membrane heterogeneity or a
weaker-identifiability artifact of MC2's narrower concentration range;
resolving that is future work (fit MC2 and MC5 jointly with a shared vs.
per-cut `|χ|,δ*` and compare via a likelihood-ratio test, once the raw-only
runs are preprocessed to add more narrow-range and wide-range replicates).

## Figures (`docs/reports/figures/`, PNG + PDF, shared `plot_style`)

- `nacl_all_estimates.{png,pdf}` — two stacked panels, `|χ|` and `δ*` point
  estimates ± 95% profile CI across the 3 datasets, one color per dataset,
  BoE reference line on the `|χ|` panel.
- `nacl_all_confidence_regions.{png,pdf}` — the three 95% nonlinear
  (F-test) confidence regions overlaid on the `(|χ|, δ*)` plane; the two MC5
  regions are drawn with different linestyles (solid/dashed) since they
  nearly perfectly overlap (visual confirmation of the repeat-run
  reproducibility above).
- `nacl_all_fits.{png,pdf}` — 3-panel `f_corrected` fit vs. `Δc_m` data, one
  panel per dataset, shared y-axis with `label_outer()`.

## Verification

- Ran in the `data3-regression` conda env; all outputs above are from that run.
- This script only *reads* `data/Rejection_Analysis.xlsx` and *writes* new
  figure files under `docs/reports/figures/` — it does not modify
  `prior_analysis/`, `data/`, `docs/reports/main.tex`, or
  `docs/reports/refs.bib`, and `latexmk` was not invoked.
- Restyled `analysis/step2_nacl/regress_nacl.py` and both
  `analysis/step3_nacl_uncertainty/*.py` (added the shared `plot_style` +
  PDF output, no logic changes) were re-run and reproduce the previously
  committed numeric estimates unchanged (verified by diffing stdout against
  the prior run): Step 2 `|χ|=25.1777/50.3555 mM`, Step 3 scipy multi-start
  100%/72.7%/100% global rates and bootstrap CIs, Step 3 ParmEst
  `θ̂=(50.355517, 0.326342)` and 100% multi-start — all identical to the
  committed values.
