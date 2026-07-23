# Prompt for Sonnet — CP Phase 2: re-run all single-salt regressions concentration-polarization–corrected (on the authors' bulk), with a CP-vs-non-CP comparison

You are completing the concentration-polarization (CP) redo. **Phase 1** put a
validated, salt-general CP correction in the preprocessing layer and computed the
per-salt mass-transfer coefficient `k`. **This is Phase 2**: re-point the single-salt
**point-estimate and uncertainty** regressions at CP-corrected interfacial
concentrations, make CP the standing convention, and present a transparent
CP-vs-non-CP comparison. Work in
`/Users/adowling/DowlingLab/Membranes/data3_regression`, env `data3-regression`
(`/opt/anaconda3/envs/data3-regression/bin/python`).

**Step 8 (AR(1), `sec:naclAR1`) and Step 9 (raw-measurement, `sec:altform`) are
Phase 3 — do NOT re-run or edit them here** (they build on the *finalized* CP NaCl
fit this phase produces). You may add a one-line note in each that a CP-corrected
re-run is Phase 3.

## The CP mechanism for Phase 2 — apply CP on the authors' bulk (not the raw reconstruction)

The Phase-1 audit established (see `sec:preproc` validation paragraph + `sec:condmodels`):
the raw→processed reconstruction reproduces the authors' NaCl sheets but is 20–26%
off for \ce{CaCl2}/\ce{LaCl3} because the authors' inline conductivity calibration is
external to the raw data (discussion point 9). **PI decision:** to isolate the CP
effect cleanly and stay comparable with all prior DATA3 results, Phase 2 applies CP
**on top of the authors' (non-CP) bulk** $\cint$, i.e. the col-B value the regressions
already read — **not** the Phase-1 raw/vial-ICP reconstruction.

**Mechanism (per dataset with an authors' sheet):**
```
c_int_CP = c_p + (c_int_authors - c_p) * exp(Jw / k_salt)      # permeate unpolarized: c_p unchanged
```
- `c_int_authors`, `c_p`, `Jw` are the columns the drivers already load (B, H, J).
- `k_salt` = the Phase-1 salt-specific mass-transfer coefficient — **import it from the
  Phase-1 code** (`analysis/preprocessing/preprocess_nacl.py`: the `Salt` config +
  `mass_transfer_coefficient(D=salt.D_s)`), the single source of truth. Do not
  recompute `k` independently. (Values for reference: NaCl 2.55e-5, CaCl2 2.25e-5,
  LaCl3 2.21e-5 m/s.)
- **CP shifts only the predictor.** The response $\Delta c_{\mathrm m}=J_{s}\ell/D_m$
  uses $\cp$ and $\Jw$ (both unchanged), so CP only raises the feed-side $\cint$
  argument of the model prediction $c_m(\cint)-c_m(\cp)$ — a clean re-fit, no change to
  the model math.

**Two data-path details:**
- **NaCl "(2)" sheet is non-CP** (its col B is bulk) → apply the mechanism. Phase 1
  validated that this reproduces the authors' *already-CP* main MC5 sheet — so for the
  all-NaCl set, the MC5 and MC5(2) rows **coincide once CP is applied**; keep a single
  MC5 (CP) entry and note the duplicate now agrees (don't plot the same run twice).
- **The 6 raw-only NaCl runs** have no authors' sheet, so "CP on authors' bulk" does
  not apply; use the Phase-1 raw-reprocessed **`*_CP.csv`** twins for them (vial-ICP
  calibration + CP). Flag this asymmetry explicitly in the write-up: the 3 processed
  NaCl + all CaCl2 + LaCl3 use authors'-bulk+CP, while the 6 raw-only NaCl use
  vial-ICP+CP (already reconstruction-noisy, as documented).

## Implementation (add a CP toggle; keep non-CP reproducible)

- Add a shared helper (e.g. `analysis/cp_authors_bulk.py`) `apply_cp(c_int_bulk, c_p, Jw, salt)`
  importing `k_salt` from Phase-1. Use it from every driver.
- Give each driver an `APPLY_CP` flag. **Default the primary run to CP=on**, but keep
  `APPLY_CP=False` runnable and **byte-for-byte reproducing the current committed
  numbers** (regression test — verify before committing).
- Produce **both** CP and non-CP results for the comparison tables.

## Regressions to re-run (Phase 2 = point estimates + uncertainty only)
1. **NaCl single** — `analysis/step2_nacl/regress_nacl.py`.
2. **NaCl uncertainty (scipy + ParmEst)** — `analysis/step3_nacl_uncertainty/*.py`.
3. **NaCl all datasets** — `analysis/nacl_all_datasets/regress_all_nacl.py` (handle the
   MC5/MC5(2) coincidence; raw-only runs via `*_CP.csv`).
4. **CaCl2** — `analysis/cacl2/regress_cacl2.py` (3 sheets).
5. **LaCl3** — `analysis/lacl3/regress_lacl3.py` (1 sheet).
Regenerate all their figures (finalized `analysis/plot_style.py`: boxed legends, bold
math; PNG+PDF; rasterize-and-eyeball each).

## Deliverable: CP-vs-non-CP comparison + does CP change the conclusion?
For each salt, a compact comparison (non-CP vs CP): $(|\chi|, K/\dstar)$ ± CI, SSE, the
scale-free fit metric, multi-start %. Report the **CP shift in $\cint$** used (Phase-1
found median +10–14% NaCl, +6.6–14.9% CaCl2, +28.7% LaCl3) and the resulting shift in
the parameters. **Answer per salt:** does CP change the *scientific conclusion*?
- NaCl: does $|\chi|\approx\SI{50}{mM}$ move materially?
- CaCl2 / LaCl3: do they still drive $|\chi|\to0$? Does the NaCl→CaCl2→LaCl3 adsorption
  trend survive? (Expect qualitative conclusions to hold, since they are
  residual-shape / effective-$\chi$-drift diagnostics, not absolute-scale — tie this to
  the "working assumptions bias absolute values, not diagnostic shapes" cross-cut,
  `sec:crosscuts`.)

## Report changes (`docs/reports/main.tex`)
- **Results sections** (§2 `sec:naclAll` + the NaCl single/UQ subsections, §3
  `sec:cacl2`, §4 `sec:lacl3`): make the **CP-corrected estimates primary**, add the
  CP-vs-non-CP comparison per salt, and state whether the conclusion changes. Keep the
  existing honest caveats (narrow-range identifiability, boundary solutions, etc.).
- **Abstract / headline numbers:** update any quoted $|\chi|$/$\dstar$ that move under
  CP, and the surrounding text.
- **Discussion point 11 ("CP vs non-CP", `sec:questions`): mark RESOLVED** — CP is now
  the standing convention throughout, computed in preprocessing (Phase 1) and applied
  in the regressions (Phase 2). Keep one sentence on the shift magnitude and that
  non-CP is retained for comparison.
- **Discussion point 9 (calibration source):** do **not** relitigate — it is already
  adjudicated in `sec:condmodels` (physics favors the vial-ICP calibration over the
  authors' inline). State plainly that Phase 2's **primary analysis uses the authors'
  bulk $\cint$ for comparability with prior DATA3 results and clean CP isolation**,
  that the vial-ICP calibration is physically preferable (a ~13–24% absolute-slope
  question, `sec:condmodels`), and that the **qualitative conclusions are robust to
  this because they are shape-based** — so the calibration-source question affects
  absolute $|\chi|$/$K$ magnitudes, not the adsorption story. (This is the third
  sensitivity axis alongside CP-vs-non-CP and transformed-vs-raw.)
- Keep the "response unchanged, only predictor shifts" note where CP is introduced.

## Guardrails
- **Do not change the model math** (Donnan polynomials, the $\tfrac12$ $f_{\text{corr}}$
  reparameterization, the transformed-response definition) — only the $\cint$ fed in.
- **Do not touch** Step 8/9 code or their sections (beyond the one-line Phase-3 note),
  `refs.bib`, `data/`, `prior_analysis/`, or `analysis/conductivity_models/`.
- Non-CP outputs must stay byte-for-byte reproducible (`APPLY_CP=False`) — verify
  before committing.
- Import `k_salt` from Phase-1 (single source of truth); do not hardcode/recompute.
- Compile clean (`latexmk -pdf`: 0 undefined, 0 overfull — watch the new comparison
  tables for overfull; `latexmk -c` after); report the page count. Rasterize + eyeball
  every regenerated figure.

## Report back
- Per salt: non-CP vs CP $(|\chi|, K/\dstar)$ ± CI, SSE, fit metric, multi-start %; the
  $\cint$ CP shift; and a one-line verdict on whether CP changes the conclusion.
- Confirmation the MC5/MC5(2) coincidence was handled and the raw-only asymmetry noted.
- Confirmation `APPLY_CP=False` reproduces the committed non-CP numbers.
- Sections/figures updated; discussion point 11 marked resolved; point 9 referenced
  (not relitigated); abstract updated.
- Clean-compile confirmation (0 undefined, 0 overfull, page count).
