# Salt-general raw-to-processed preprocessing pipeline (NaCl, CaCl2, LaCl3)

Converts RAW diafiltration campaign sheets (`data/NF270_MC{2,3,4,5}.xlsx`)
into the PROCESSED analysis columns the regressions consume (interfacial
concentration, ionic strength, rejection, permeate interfacial
concentration, water flux, solute permeability). Originally built for
single-salt NaCl (unlocking 6 previously raw-only NaCl runs for the
multi-dataset regression); **Phase 1 of the CP redo**
(`docs/prompts/cp_preprocessing_extend.md`) generalized it to CaCl2 (2:1)
and LaCl3 (3:1) and moved **concentration polarization (CP)** into this
layer, computed once from raw signals for every single-salt dataset (see
"Phase 1: salt generalization + CP" below). **Phase 2** (re-pointing the
regressions at the CP-corrected CSVs, the CP-vs-non-CP comparison, and the
results-section rewrites) is a later, separate task — not done here.

## Files

- `preprocess_nacl.py` — the salt-general pipeline (module name kept for
  backward compatibility with existing NaCl importers): raw-sheet loading,
  water flux (Jw) from a smoothed mass-vs-time slope, per-experiment
  conductivity calibration, the `Salt` dataclass (per-salt ICP molar mass,
  valence-aware ionic-strength multiplier, CP solution diffusivity),
  concentration polarization (CP, optional), derived columns.
- `validate.py` — validates the original NaCl pipeline against the two runs
  where we have both raw and already-processed (Excel) data.
- `process_raw_runs.py` — applies the pipeline to the 6 raw-only NaCl runs;
  writes both the original non-CP `processed/*.csv` (byte-for-byte
  unchanged by the Phase 1 refactor — regression-tested, see below) and new
  CP-corrected `processed/*_CP.csv` twins.
- `process_raw_multisalt.py` — **Phase 1**: applies the salt-general
  pipeline (both `apply_cp=True` and `apply_cp=False`) to the 3 CaCl2 runs,
  the 1 LaCl3 run, and the 2 NaCl "validation twin" runs, and validates each
  against the authors' processed sheets. See "Phase 1" below.
- `processed/*.csv` — output datasets, one row per raw-sheet time point,
  columns `c_int_mM, ionic_strength_mM, rejection, c_p_mM, Jw_m3_m2_s,
  B_um_s, time_s, ret_cond_uS_cm, perm_cond_uS_cm`, with the fitted
  calibration and (as of Phase 1) salt/CP provenance recorded in `#`-comment
  header lines.

## How to run

```bash
conda activate data3-regression
python analysis/preprocessing/validate.py            # ~4 min (5 large-workbook loads)
python analysis/preprocessing/process_raw_runs.py     # ~3 min (6 raw-workbook loads x2 CP states)
python analysis/preprocessing/process_raw_multisalt.py  # ~10 min (CaCl2/LaCl3/NaCl-twin x2 CP states + validation)
python analysis/nacl_all_datasets/regress_all_nacl.py  # picks up the (unchanged) non-CP NaCl CSVs automatically
```

## Recipe (equations)

Same equations as `docs/reports/main.tex` §1 "Data preprocessing" (label
`sec:preproc`):

1. **Water flux**: `Jw = (1/(rho*A_m)) * dm/dt`, `dm/dt` from a centered
   rolling-window linear-regression slope of mass(t) (window tuned to 51
   rows against the one reference `J` column we have — see Validation 1).
2. **Bulk concentration**: `c_bulk = (cond - intercept) / slope`, calibration
   fit per-experiment from that run's own vial (conductivity, ICP) pairs.
3. **Interfacial concentration** (thin-film CP, optional via `apply_cp`):
   `c_int = c_p + (c_bulk_ret - c_p) * exp(Jw/k)`, permeate assumed
   unpolarized (`c_p` = bulk permeate concentration directly). If
   `apply_cp=False`, bulk == interfacial on both sides (matches the "(2)"
   sheet convention Step 2/Step 3/§2.4 already use).
4. **Derived**: `I = c_int` (1:1 salt), `R = 1 - c_p/c_int`,
   `B = Jw*c_p/(c_int-c_p) * 1e6` [µm/s].

For CaCl2/LaCl3, step 4's `I` and step 3's CP diffusivity `D` generalize
per-salt — see "Phase 1" below.

## Phase 1: salt generalization + CP (`docs/prompts/cp_preprocessing_extend.md`)

**Refactor approach.** `preprocess_nacl.py` keeps its module name and
public API (`load_raw_sheet`, `fit_calibration`, `process_experiment`, ...)
unchanged for every existing NaCl caller — all of them call it with no
`salt=` argument. A new `Salt` frozen dataclass (`z_cation`,
`icp_molar_mass`, `d_cation_m2_s`) derives the valence-aware ionic-strength
multiplier and the CP-relevant ambipolar solution diffusivity from a single
number, `z_cation`:
```
I = 0.5 * z*(z+1) * c_salt            # 1 (NaCl), 3 (CaCl2), 6 (LaCl3)
D_s = (1+z) * D_cation * D_Cl / (z*D_cation + D_Cl)
```
`process_experiment(..., salt=SALT_NACL)` defaults to NaCl, so every
existing call site (unmodified) reproduces its previous behavior exactly.
**Regression test**: reprocessing all 6 raw-only NaCl runs with
`apply_cp=False` after the refactor is **byte-for-byte identical** to the
pre-refactor CSVs (`diff` clean on all 6 files) — confirmed before writing
any CaCl2/LaCl3 code.

**Per-salt constants:**

| Salt | z | ICP molar mass (elemental) | `D_cation` [m²/s] | `D_s` (ambipolar, solution) [m²/s] | `k` [m/s] | `I` rule |
|---|---|---|---|---|---|---|
| NaCl | 1 | 22.99 (Na⁺) | 1.33e-9 | 1.6071e-9 | 2.5473e-5 | `I = c` |
| CaCl2 | 2 | 40.08 (Ca²⁺) | 0.79e-9 | 1.3327e-9 | 2.2471e-5 | `I = 3c` |
| LaCl3 | 3 | 138.91 (La³⁺) | 0.626e-9 | 1.3007e-9 | 2.2107e-5 | `I = 6c` |

(Shared anion `D_Cl = 2.03e-9` m²/s. `D_cation`/`D_Cl` are the
**solution-phase** values used for the CP boundary layer — ×1000 smaller
than the **membrane-phase** `D`'s the regressions use for `B`,
`analysis/{step2_nacl,cacl2,lacl3}/regress_*.py`'s `D_NA_M`/`D_CA_M`/
`D_LA_M`/`D_CL_M` — one source of truth, see `preprocess_nacl.py`'s
`Salt.solution_D` docstring.)

**Datasets processed (all with CP on AND off; `process_raw_multisalt.py`,
plus the 6 raw-only NaCl runs' new `_CP.csv` twins from `process_raw_runs.py`):**

| Dataset | Calibration source | `(a, s, n_vial, R²)` |
|---|---|---|
| CaCl2 MC3 07.11.24_SCaCl2 | fit, own vials | `(96.406, 187.577, 14, 0.9756)` |
| CaCl2 MC5 05.27.26_CaCl2 | **given** (embedded formula, see below) | `(0, 149.73, —, —)` |
| CaCl2 MC5 06.27.26_CaCl2 | **approximated** (05.27.26's formula; flagged PROVISIONAL) | `(0, 149.73, —, —)` |
| LaCl3 MC2 05.21.24_LaCl3 | fit, own vials | `(93.837, 260.084, 14, 0.9991)` |
| NaCl MC5 05.27.26_NaCl (twin) | given (pre-existing embedded formula) | `(-63.706, 76.685, —, —)` |
| NaCl MC2 05.07.24_NaCl | approximated (MC5 formula; flagged PROVISIONAL, pre-existing) | `(-63.706, 76.685, —, —)` |

**Deviation from the task prompt, flagged.** The prompt states the CaCl2/
LaCl3 raw sheets "DO record" per-vial conductivity. True for
`07.11.24_SCaCl2` and `05.21.24_LaCl3` (used to fit their own calibrations
above). **False** for both MC5 CaCl2 sheets (`05.27.26_CaCl2`,
`06.27.26_CaCl2` — vial "Conductivity @ Temp" is entirely blank), the same
situation the two NaCl validation-twin sheets were already in. For
`05.27.26_CaCl2` we found a **live-formula fallback** by inspecting cell
formulas (`data_only=False`), exactly as was previously done for the NaCl
"(2)" sheet: `=N{row}/149.73` (intercept 0, slope 149.73) for the sheet's
first several rows. For `06.27.26_CaCl2` there is **no live formula
anywhere** in its column either — its calibration is genuinely
unrecoverable from this workbook, so it is approximated with
`05.27.26_CaCl2`'s formula (same membrane cut/campaign, closest date) and
flagged PROVISIONAL, mirroring the pre-existing MC2 05.07.24_NaCl
precedent.

**CP shift per dataset** (median / max % change in `c_int`, non-CP → CP):

| Dataset | median | max |
|---|---|---|
| CaCl2 MC3 07.11.24_SCaCl2 | +14.9% | +30.7% |
| CaCl2 MC5 05.27.26_CaCl2 | +7.8% | +37.8% |
| CaCl2 MC5 06.27.26_CaCl2 | +6.6% | +41.1% |
| LaCl3 MC2 05.21.24_LaCl3 | +28.7% | +57.8% |
| NaCl MC5 05.27.26_NaCl (twin) | +10.2% | +28.3% |
| NaCl MC2 05.07.24_NaCl | +14.3% | +29.3% |

The NaCl twin's +10.2% median is in the right ballpark of, if a bit below,
the task prompt's expected "~12–22%" — consistent with genuine CP physics,
not a stray implementation difference (see the CP-on cross-check below,
which reproduces the pre-existing validate.py numbers to 2 decimal places).

**Validation — CP-off vs. the authors' non-CP processed sheets** (median /
max relative error; `B`=c_int, `H`=c_p, `D`=ionic strength, `J`=Jw,
`F`=rejection, `K`=B coefficient):

| Dataset | B | H | J | F | K |
|---|---|---|---|---|---|
| CaCl2 MC3 07.11.24_SCaCl2 | 23.4% / 68.0% | 11.2% / 96.8% | 4.1% / 99.5% | 13.3% / 5573% | 36.8% / 69.0% |
| CaCl2 MC5 05.27.26_CaCl2 | 1.5% / 9.9% | 15.0% / 596% | 6.3% / 99.5% | 8.0% / 5863% | 30.5% / 100% |
| CaCl2 MC5 06.27.26_CaCl2 | 3.4% / 29.0% | 15.0% / 93.1% | 7.9% / 101% | 6.1% / 186% | 25.8% / 100% |
| LaCl3 MC2 05.21.24_LaCl3 | 26.4% / 55.6% | 33.5% / 297% | 4.7% / 99.5% | 0.4% / 31.7% | 7.0% / 248% |
| NaCl MC5 05.27.26_NaCl (2), twin (pre-existing, for comparison) | 1.0% / 38.4% | 19.9% / 96.7% | 3.8% / 99.7% | n/a (blank) | 71.8% / 100% |

**Validation — CP-on vs. the one authors' already-CP-corrected sheet**
(NaCl only; `B`/`H`, median / max relative error):

| Dataset | B (c_int, CP) | H (c_p) |
|---|---|---|
| NaCl MC5 05.27.26_NaCl (main) | 11.9% / 208% | 19.9% / 96.7% |
| NaCl MC2 05.07.24_NaCl (approx. calib.) | 21.9% / 34.1% | 28.7% / 122% |

These two rows reproduce the pre-existing `validate.py` numbers (documented
above in "Validation 2") to within rounding — the Phase 1 refactor changed
nothing about the NaCl CP path itself; this is a consistency check, not new
information.

**Do not hide failures — CaCl2/LaCl3 CP-off reconstruction does NOT reach
NaCl's validation quality. Root cause (revised after Opus audit): a
calibration-SOURCE mismatch, NOT nonlinearity.**

NaCl's CP-off `B` (c_int) reproduces the authors' sheet to **~1% median**
(both the twin and the "(2)" sheet). CaCl2's three datasets range from
**1.5% to 23.4% median**, and LaCl3 is **26.4% median**.

The original hypothesis here was a "nonlinear conductivity→concentration
relationship." **That is wrong.** Recovering the authors' *own* implied
calibration from their processed sheets (their c_int vs. their inline
retentate conductivity) shows it is **linear**: `kappa = -132.6 + 152.1*c`
for CaCl2 MC3 and `kappa = -302.6 + 216.3*c` for LaCl3, both R²≈0.999 (a
quadratic term barely improves it). The real cause is that we fit the
calibration from the **vial (conductivity, ICP) pairs** (`kappa = +96.4 +
187.6*c` for CaCl2 MC3), whereas the authors used a different **inline**
conductivity calibration — different slope (~20%) and opposite-sign
intercept. The shrinking reference/ours ratio (~3.1× at low c → ~1.3× at
high c) is exactly the signature of two *linear* calibrations with
different intercepts (intercept-dominated at low c; approaching the slope
ratio 187.6/152.1 ≈ 1.23 at high c), which is what first (mis)suggested
nonlinearity. Rejection `F` agrees far better than `B`/`H` for LaCl3
(0.4% vs. 26%) because `R=1-c_p/c_int` cancels the shared multiplicative
calibration mismatch — consistent with a systematic calibration
difference, not noise.

For NaCl the raw workbook carries the authors' embedded inline calibration
(so we reproduce their sheet); the raw **CaCl2/LaCl3 workbooks contain no
embedded conductivity calibration**, and no vial-fit convention
(elemental/molecular molar mass, with/without dilution) reproduces the
authors' coefficients — their calibration is external to the data we have.
This is exactly the **vial-vs-inline instrument question (report discussion
point 9)**, now quantified at 20–26% for multivalent salts. **Decision (PI,
after audit):** the CaCl2/LaCl3 CP analysis (Phase 2) applies CP *on top of
the authors' bulk c_int* (their linear calibration) to isolate the CP
effect and reproduce prior work; the raw-derived vial-ICP CSVs here are
retained as a labeled calibration-sensitivity artifact.

## Working choices / flags

1. **`A_m` (≈3.95 cm²) and `rho`** are not documented upstream; unchanged
   from the value already in `docs/reports/main.tex` (inferred by an earlier
   session by back-solving the Jw equation). Still a flag for the team to
   confirm against the actual stirred-cell spec sheet.
2. **Jw smoothing**: centered rolling-window linear-regression slope,
   **window = 51 rows**, tuned by minimizing median relative error against
   the one available reference `J` column (05.27.26_NaCl "(2)" sheet):
   3.8% median error. This is noticeably worse than the ~1e-6 we can get for
   calibration-formula columns — Jw is fundamentally a derivative of a noisy
   signal and 51-row smoothing is a genuine approximation, not a bug.
3. **Conductivity calibration — fit from vial data, with a real bug found
   and fixed along the way:**
   - **The "ICP Salt 1 (mg/L)" column is the concentration of the DILUTED
     ICP sample, not the original vial.** It must be scaled up by the
     dilution factor `(sample_vol_mL + acid_vol_mL) / sample_vol_mL` (from
     the same vial-data block) before converting mg/L → mM. Missing this
     (an earlier version of this pipeline did) silently understates
     concentration by 40–200x and produces a physically impossible
     calibration slope (~12,000 µS/cm per mM, vs. real NaCl solution
     conductivity of ~O(100) µS/cm per mM) with poor R² (0.32–0.76). After
     the fix, the Feed/Diafiltrate ICP-derived concentrations land within
     ~30% of the values stated in each run's own notes (e.g. "1 mM NaCl
     Feed, 150 mM NaCl Diafiltrate"), and calibration R² jumps to
     0.995–0.9995 — this cross-check is what caught the bug.
   - **ICP molar mass: elemental Na⁺ (22.99 g/mol), not molecular NaCl
     (58.44 g/mol)**, since ICP-OES is an elemental technique. Confirmed
     by the same Feed/Diafiltrate sanity check above (using 58.44 would
     put concentrations off by 2.54x in the wrong direction).
   - **The two validation-target sheets ("05.27.26_NaCl" and
     "05.07.24_NaCl") have NO per-vial conductivity recorded at all**
     (column P entirely blank for every vial) — unlike the 6 raw-only runs,
     which all have full 14-point (13 for one MC4 run) vial conductivity+ICP
     data. This was not anticipated in the original task plan (which
     expected to fit-and-compare a per-experiment calibration for
     "05.27.26" specifically). We could not test "does 05.27.26's own
     vial-fit calibration match the embedded formula" directly for that
     reason; instead we used the embedded formula
     (`cond = -63.706 + 76.685*c`) as given/ground-truth for that sheet, and
     validated the *general* fitting approach (dilution correction + Na
     molar mass) on the 6 runs that do have vial conductivity — see the
     calibration table below, where the fitted slopes (93–98) land close to
     (about 25% above) the embedded formula's 76.685, a reasonable
     cross-check given they're different membrane cuts/sessions/probes.
4. **CP vs non-CP**: implemented as the `apply_cp` flag. The 6 new CSVs use
   `apply_cp=False` for direct comparability with the existing Step
   2/Step 3/§2.4 results (which used the non-CP "(2)" sheet convention).
5. **H/permeate-side smoothing**: approximated with a short floor (first 20
   rows at 0.04 mM) + rolling-median smooth. Investigating the reference
   sheet directly suggests a physical explanation for the real flattening
   artifact: columns O (permeate conductivity) are themselves blank for the
   first ~14–25 rows of the reference sheets, consistent with sensor
   dead-volume before permeate first reaches the conductivity cell — not
   arbitrary smoothing. Not reproduced exactly; flagged as approximate.

## Validation 1 (exact target): NF270_MC5 "05.27.26_NaCl (2)"

**Important finding, not anticipated going in**: inspecting the sheet's
cell *formulas* (not just values) shows column B (`c_int`) is a **live
formula for only the first 22 of 747 rows** (`=(N{row}+63.706)/76.685`); the
remaining 725 rows (97% of the sheet) are **hardcoded pasted numbers** that
do *not* satisfy that formula exactly (they regress against the verbatim
retentate conductivity with R²=0.999 but a visibly different
effective slope/intercept, and residuals up to ~4 mM) — i.e. most of this
"exact target" sheet was actually produced by some undocumented downstream
step (a refined/different calibration, or one computed from smoothed rather
than instantaneous conductivity), not by the formula visible in the cell.
We validate the two regimes separately:

| Region | n | median rel. err | max rel. err |
|---|---|---|---|
| B, rows 2–23 (live formula) | 22 | 0.0 | 1.8e-16 (machine precision) |
| B, rows 24–748 (hardcoded/pasted) | 725 | 1.0e-2 | 3.8e-1 |

Full-sheet column match (using our own smoothed Jw, window=51):

| col | meaning | median rel err | max rel err |
|---|---|---|---|
| B | c_int (calibration) | 9.8e-3 | 3.8e-1 |
| D | ionic strength | 9.8e-3 | 3.8e-1 |
| M | time (verbatim) | 0.0 | 0.0 |
| N | ret. cond. (verbatim) | 0.0 | 2.2e-16 |
| J | Jw (smoothed) | 3.8e-2 | 1.0 |
| K | B coefficient (our Jw, ref H) | 7.2e-1 | 1.0 |

K formula check (our B, but the *reference* J and H — isolates the K
formula from Jw/H-smoothing uncertainty): median 2.3e-2, max 5.5e-1. F
(rejection) is blank throughout this sheet — no reference to check.

## Validation 2 (approximate, CP-corrected targets)

| Run | col | median rel err | max rel err |
|---|---|---|---|
| MC5 05.27.26_NaCl (main, CP-corrected; same raw source+calibration) | B | 11.9% | 208% |
| " | H | 19.9% | 96.7% |
| MC2 05.07.24_NaCl (different cut/session; calibration UNKNOWN, approximated with MC5 constants) | B | 21.9% | 34.1% |
| " | F | 6.0% | 29.1% |
| " | H | 28.7% | 122% |

As expected, MC2's mismatch is dominated by the unknown/approximated
calibration (different membrane cut and session, so almost certainly a
different probe/cell constant), not by the CP model itself — we have no way
to independently verify MC2's true calibration since it has neither vial
conductivity nor a live formula.

## Per-experiment calibrations, 6 raw-only runs

Non-CP (`apply_cp=False`), fit from each run's own 13–14 vial
(conductivity, dilution-corrected ICP) pairs:

| Run | intercept | slope | n | R² |
|---|---|---|---|---|
| MC3 07.09.24_NaCl | 145.02 | 98.02 | 14 | 0.9951 |
| MC3 07.22.24_SNaCl | 88.79 | 96.60 | 14 | 0.9995 |
| MC4 07.11.24_SNaCl | 17.35 | 95.36 | 13 | 0.9971 |
| MC5 07.23.24_NaCl | 111.33 | 93.00 | 14 | 0.9962 |
| MC5 07.23.24_SNaCl | 177.02 | 93.21 | 14 | 0.9965 |
| MC5 07.23.24_S2NaCl | 136.69 | 98.42 | 14 | 0.9991 |
| *(reference)* embedded MC5 05.27.26 formula | -63.71 | 76.69 | — | — |

Slopes are remarkably consistent across all 6 runs (93–98) despite spanning
3 membrane cuts and 2 months — a good internal-consistency check on the
dilution-correction fix — and land within ~25% of the embedded formula's
slope, a reasonable cross-probe/session agreement. **MC4's missing "Final
Tube" vial conductivity (13/14 points) did not degrade its fit** (R²=0.997,
in line with the other 5).

## Expanded multi-dataset regression (9 NaCl runs total)

Re-ran `analysis/nacl_all_datasets/regress_all_nacl.py` across the original
3 Excel-sourced datasets plus these 6 new CSVs (non-CP throughout, for
comparability). **Bottom line: adding the wide-range raw-only runs does
*not* sharpen the MC2-vs-MC5 identifiability picture — it adds a third,
noisier tier of scatter on top of it:**

| Dataset | n | `\|χ\|` [mM] | `δ*` | Note |
|---|---|---|---|---|
| MC2 (05.07.24) | 482 | 18.46 | 0.2888 | Excel, clean |
| MC5 (05.27.26) | 747 | 49.56 | 0.3251 | Excel, clean |
| MC5 (05.27.26) (2) | 747 | 49.55 | 0.3251 | Excel, clean (repeat of above) |
| MC3 (07.09.24) | 391 | — | — | **Ill-conditioned** (diverges to nonphysical `\|χ\|`~2.1e6 mM) |
| MC3 (07.22.24, S) | 1312 | 20.33 | 0.3000 | 95% CI (12.4, 31.9) — wide |
| MC4 (07.11.24, S) | 317 | 106.63 | 1.931 | 95% CI (94.5, 120.1) — far from all others |
| MC5 (07.23.24) | 496 | 61.99 | 0.7035 | 95% CI (60.5, 63.5) — tighter but still 25% off MC5 (05.27.26) |
| MC5 (07.23.24, S) | 500 | 169.15 | 0.7147 | 95% CI (88, 535) — barely identified |
| MC5 (07.23.24, S2) | 1324 | 35.63 | 0.3434 | 95% CI (24.4, 52.1) — overlaps MC5 (05.27.26) |

Three findings:
1. **The two clean Excel-sourced MC5 runs still agree beautifully with each
   other** (as before); the new preprocessed runs do not reproduce that
   level of agreement even among themselves (three different MC5 07.23.24
   runs give 62, 169, and 36 mM) — this is much more scatter than the
   Excel-sourced repeat showed, and is best attributed to the preprocessing
   pipeline's own approximations (calibration R²~0.995–0.9995 rather than
   near-machine precision, 51-row Jw smoothing, heuristic CP-off/H
   handling) compounding into the fit, not to genuine membrane-to-membrane
   variability.
2. **MC3 07.09.24_NaCl is ill-conditioned**: the reconstructed data shows a
   near-zero apparent rejection (`c_int` and `c_p` converge to within a few
   mM of each other over most of the run), which starves the 2-parameter
   Donnan model of curvature to fit and lets `least_squares` wander to a
   nonphysical optimum. `regress_all_nacl.py` now detects this
   (`\|χ\|`/`δ*` > 1e4) and skips multi-start/profile-CI for that dataset
   rather than crashing or reporting garbage confidence intervals.
3. **The `S`/`S2` (sequential/reverse-diluting) runs show an anomalous
   "hook"** at the high-concentration end of the corrected-model fit (see
   `nacl_all_fits.png`) and are extremely sensitive to the startup-trim
   sensitivity check (e.g. MC5 07.23.24 S2's `\|χ\|` moves +17,000% under a
   7% row trim) — these dilution-reversal runs likely need a different
   preprocessing/trimming treatment than the standard NaCl runs, which is
   future work, not attempted here.

**Practical recommendation**: treat the 3 original Excel-sourced datasets as
the reliable basis for the reproducibility conclusion (MC5 repeat
excellent, MC2 genuinely different); treat the 6 newly-unlocked raw-only
datasets as provisional pending either (a) a more precise Jw-smoothing/CP
treatment, or (b) obtaining/restoring the actual live calibration formulas
these older runs were originally processed with (analogous to the "(2)"
sheet's embedded formula), if such a record exists.
