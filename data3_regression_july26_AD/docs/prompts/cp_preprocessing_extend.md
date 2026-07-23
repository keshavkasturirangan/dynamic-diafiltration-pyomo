# Prompt for Sonnet — Phase 1 of the CP redo: extend the preprocessing pipeline to CaCl2/LaCl3 and recompute concentration polarization FROM RAW for every single-salt dataset

You are moving **concentration polarization (CP)** into the **data-preprocessing
layer** so it is computed once, from the raw signals, in a single validated place —
the most robust architecture (PI's decision). This is **Phase 1**: extend the
already-validated NaCl preprocessing pipeline to **CaCl2 (2:1)** and **LaCl3 (3:1)**,
**reprocess every analyzed single-salt dataset from its raw campaign workbook** with
CP applied, and **validate the reconstruction against the authors' processed
sheets**. The regression re-runs and report results-section rewrites are **Phase 2**
(a later prompt) — do **not** touch the regression drivers or the salt results
sections here.

> **Phase 2 scope, for context (a later prompt — do NOT do it here):** re-point all
> the single-salt regressions at the CP-corrected CSVs and re-run them, **including
> the Step 8 AR(1) autocorrelation correction (`sec:naclAR1`) and the Step 9
> raw-measurement reformulation (`sec:altform`)** — both build on the NaCl fit, so
> CP shifts their inputs and they must be recomputed CP-corrected (Step 9's
> `c_int(κ_r)` reconstruction should use the same CP mass-transfer coefficient this
> phase defines). Phase 2 also does the CP-vs-non-CP comparison, the results
> sections (§2–§4), the abstract, and the discussion-point resolution. Everything in
> Phase 1 must therefore keep the non-CP path reproducible so Phase 2 can compare.

Work in `/Users/adowling/DowlingLab/Membranes/data3_regression`, env
`data3-regression` (`/opt/anaconda3/envs/data3-regression/bin/python`).

**Why (the discussion point being resolved).** The preliminary Excel analysis
skipped CP for speed; the report fits **non-CP** interfacial concentrations
throughout and flags "CP vs. non-CP" as an open team decision (§1 discussion list;
the `\todo` at `sec:preproc` follow-up item 5). Recomputing CP inside preprocessing
resolves it robustly: every dataset flows through one pipeline, CP is documented in
one place, and the regressions (Phase 2) just consume CP-corrected processed CSVs.

## What already exists (build on it — do not rewrite from scratch)

`analysis/preprocessing/preprocess_nacl.py` is the **validated NaCl pipeline**:
- `load_raw_sheet(path, sheet)` parses the raw campaign schema (time, mass,
  pressure, retentate/permeate temp+conductivity, vial-swap; a Vial Data block with
  per-vial conductivity, sample/acid volumes, ICP mg/L; and ICP calibration blocks).
- `fit_calibration(vials, icp_molar_mass=...)` fits `cond = a + s*c_mM` from a run's
  own vial (conductivity, ICP) pairs — with the **ICP dilution fix** already in place
  (`c_mM = icp_mgL * (sample_vol+acid_vol)/sample_vol / molar_mass`; skipping the
  dilution silently understates concentration ~40–200×).
- `mass_transfer_coefficient(D=...)` = `0.23 v0^0.57 D^0.67 nu^-0.24 b^-0.43`,
  `v0=omega*b/2`, `rpm=350`, `b=0.0254`, `nu=1.003e-6`.
- `compute_interfacial(c_bulk_ret, c_bulk_perm, Jw, apply_cp=, k=)` applies the
  thin-film CP (`c_int = c_p + (c_bulk_ret - c_p) exp(Jw/k)`; permeate unpolarized).
- `process_experiment(raw, apply_cp=, ...)` ties it together → `c_int_mM`,
  `c_p_mM`, `ionic_strength_mM`, `rejection`, `B_um_s`, `Jw_m3_m2_s`, ...
- `process_raw_runs.py` drives the 6 raw-only NaCl runs (currently `apply_cp=False`).

**Confirmed:** the raw sheet schema is **identical across NaCl/CaCl2/LaCl3** (same
27-column layout), so `load_raw_sheet` already parses all three. The generalization
is per-salt **constants and derived-quantity formulas**, not parsing.

## The generalization to 2:1 and 3:1 salts

Refactor `preprocess_nacl.py` into a **salt-general** module (e.g. add a `Salt`
config/dataclass, or a `preprocess.py` that parametrizes the NaCl one) carrying, per
salt:

1. **ICP metal molar mass** (ICP-OES measures the elemental metal; moles metal =
   moles salt, so metal molarity = salt molarity):
   - Na⁺ = 22.99 (unchanged), **Ca²⁺ = 40.08**, **La³⁺ = 138.91** g/mol.
   The dilution correction is identical; only the molar mass changes in
   `fit_calibration`.

2. **Valence-aware ionic strength** `I = ½ Σ z_i² c_i` for a fully-dissociated salt
   at salt concentration `c`:
   - NaCl (1:1): `I = c` (unchanged).  CaCl2 (2:1): `I = 3c`.  LaCl3 (3:1): `I = 6c`.

3. **Ambipolar (Nernst–Hartley) SOLUTION diffusivity** for the CP mass-transfer
   coefficient `k` (same functional form the regressions use for the *membrane* `D_m`,
   but with **solution** ion diffusivities — the codebase's membrane constants × 1000):
   - `D_Na=1.33e-9`, `D_Ca=0.79e-9`, `D_La=0.626e-9`, `D_Cl=2.03e-9` m²/s.
   - NaCl: `D_s = 2 D_Na D_Cl/(D_Na+D_Cl)` (= existing `D_NACL_M2_S`).
   - CaCl2: `D_s = 3 D_Ca D_Cl/(2 D_Ca + D_Cl)`.
   - LaCl3: `D_s = 4 D_La D_Cl/(3 D_La + D_Cl)`.
   Compute and **report the resulting `k` for each salt**. (Cross-check: the cacl2/
   lacl3 regression scripts define `D_CA_M=0.79e-12`, `D_LA_M=0.626e-12`,
   `D_CL_M=2.03e-12` — these are the ×1000-smaller membrane values; reuse them ×1000
   for the solution `D_s` so there is one source of truth, with a clear comment.)

4. **Rejection / B**: unchanged formulas on salt concentrations
   (`R = 1 - c_p/c_int`, `B = Jw c_p/(c_int-c_p)`).

The **co-ion stoichiometry** (Cl⁻ = 2c or 3c) belongs to the *regression* model
(Phase 2), NOT preprocessing — preprocessing emits **salt** concentrations
`c_int_mM`, `c_p_mM`.

**Regression-test NaCl:** the refactor must leave NaCl outputs **byte-for-byte
unchanged** (`apply_cp=False` reproduces the current processed CSVs exactly). Verify
before proceeding.

## Datasets to (re)process from raw (all with CP on; also emit non-CP for validation/Phase-2 comparison)

Raw campaign workbooks in `data/` (all present, verified):
- **NaCl:** `NF270_MC5.xlsx::05.27.26_NaCl` (the validation twin — see below),
  `NF270_MC2.xlsx::05.07.24_NaCl`, and the 6 raw-only runs already in
  `process_raw_runs.py`.
- **CaCl2:** `NF270_MC3.xlsx::07.11.24_SCaCl2`, `NF270_MC5.xlsx::05.27.26_CaCl2`,
  `NF270_MC5.xlsx::06.27.26_CaCl2`.
- **LaCl3:** `NF270_MC2.xlsx::05.21.24_LaCl3`.

**Calibration source per sheet:** fit from the run's own vial (conductivity, ICP)
pairs where per-vial conductivity is recorded — **the CaCl2 and LaCl3 raw sheets DO
record it** (Vial Data "Conductivity @ Temp" is populated). Where it is absent (the
NaCl `05.27.26`/`05.07.24` sheets record no per-vial conductivity), keep the existing
fallback (embedded/stated calibration), exactly as the current pipeline does. Report
each dataset's calibration `(intercept, slope, n_vial_points, R²)`.

Emit **CP-corrected** processed CSVs to `analysis/preprocessing/processed/` (the new
standing data), and also write the **non-CP** version (suffix or a flag) so Phase 2
can build the CP-vs-non-CP comparison. Each CSV header must record provenance:
source sheet, calibration, ICP molar mass, salt, valence-aware `I` rule, CP `D_s`
and `k`, and `apply_cp`.

## Validation — this is the gate; do it before trusting anything downstream

1. **CP-off reproduces the authors' non-CP processed sheets.** Reprocess with
   `apply_cp=False` and compare column-by-column against the authors' processed
   sheets (which we determined are **non-CP**): LaCl3 `NF270_MC2 05.21.24_LaCl3`, the
   3 CaCl2 sheets, and NaCl `NF270_MC5 05.27.26_NaCl (2)` (all in `BoE Analysis.xlsx`/
   `Rejection_Analysis.xlsx`). Report **median and max relative error** on `c_int`,
   `c_p`, `Jw`, `I`, `R`, `B`, plus parity plots. (Diagnostic already run: for a
   non-CP sheet `c_int/κ_retentate` has CV ≈ 0.05–0.09; for the CP sheet ≈ 0.25.)
2. **CP-on reproduces the authors' already-CP NaCl sheet.** The non-"(2)" sheet
   `NF270_MC5 05.27.26_NaCl` is **already CP-corrected** by the authors. Reprocess
   `05.27.26_NaCl` with `apply_cp=True` and confirm it reproduces that sheet's
   `c_int` — this validates the CP implementation end-to-end. Report the agreement
   and the **CP shift magnitude** (report expects a median ~12–22% rise in `c_int`).
3. **Report the CP shift per salt** (median/max % change in `c_int` from CP-off to
   CP-on) for CaCl2 and LaCl3 too.
4. **Do not hide failures.** If a dataset's CP-off reconstruction does not reproduce
   its authors' sheet (reconstruction noise — as happened for the 6 raw-only NaCl
   runs), say so explicitly with the error magnitude and a hypothesis; flag it for
   the PI rather than silently accepting it.

## Report write-up — Phase 1 scope only (methods + validation, NOT the salt results)

Edit `docs/reports/main.tex` **only** in `sec:preproc`/`sec:cp` (and the
`sec:preproc` follow-up items). Document the salt-general pipeline: per-salt ICP
molar mass, valence-aware ionic strength, and the per-salt CP `k`/`D_s`; state that
CP is now recomputed from raw in preprocessing for all single salts; add the
validation results (the CP-off parity against non-CP sheets, and the CP-on
reproduction of the authors' CP NaCl sheet, with the shift magnitude). **Leave the
salt results sections (§2 `sec:naclAll`, §3 `sec:cacl2`, §4 `sec:lacl3`), the
regression drivers, the abstract headline numbers, and the discussion-point
resolution for Phase 2** — but you may downgrade the `sec:preproc` follow-up item 5
`\todo` to a note that CP is now implemented and the results-level resolution is
pending Phase 2.

## Guardrails
- **Do not touch** the regression drivers (`analysis/step2_nacl/`,
  `step3_nacl_uncertainty/`, `cacl2/`, `lacl3/`, `nacl_all_datasets/`), the salt
  results sections, the abstract, `refs.bib`, `data/`, or `prior_analysis/`.
- NaCl `apply_cp=False` outputs must stay **byte-for-byte identical** (regression
  test) — the NaCl multi-dataset analysis currently reads these CSVs.
- Extend `analysis/preprocessing/README.md` with the salt generalization, the
  per-salt constants, and the validation results.
- If you edit the report, compile clean (`latexmk -pdf`: 0 undefined, 0 overfull;
  `latexmk -c` after) and report the page count.

## Report back
- The refactor approach and confirmation NaCl `apply_cp=False` outputs are unchanged.
- Per-salt: ICP molar mass, valence-aware `I` rule, CP `D_s` and `k`, and each
  dataset's calibration `(a, s, n_vials, R²)`.
- **Validation:** CP-off median/max relative error vs each authors' non-CP sheet;
  the CP-on reproduction of the authors' CP NaCl sheet; the CP shift magnitude per
  salt; and any dataset that fails to reproduce (with magnitude + hypothesis).
- The list of CP-corrected (and non-CP) CSVs written, with their provenance headers.
- `sec:preproc`/`sec:cp` changes made; clean-compile confirmation if the report was
  edited. (Note that Phase 2 — regression re-runs, results sections, discussion-point
  resolution, abstract — is deliberately not done here.)
