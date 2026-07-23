# Prompt for Sonnet — NaCl raw→processed preprocessing pipeline (validated), + document the methods in report Section 1

You are building a **validated preprocessing pipeline** that converts RAW
diafiltration campaign data into the PROCESSED analysis columns our regression
consumes, for single-salt NaCl. Work in
`/Users/adowling/DowlingLab/Membranes/data3_regression`, in the `data3-regression`
conda env. Two goals of equal importance: **(A)** produce & validate the pipeline
and process the raw-only runs, and **(B)** document the preprocessing methods, with
mathematical precision, in **Section 1 of the report** (`docs/reports/main.tex`).

Why this matters: it unlocks 6 more NaCl experiments (helping the MC2 narrow-range
identifiability question), and it is the forward model for the later raw
time-series reformulation (report §2.3 `sec:altform`).

## Raw sheet schema (per-experiment tab in `data/NF270_MC{2,3,4,5}.xlsx`)
Row 1: `Datapoints: N | Notes: <text>`. Row 3: time-series headers, cols **A–H** =
`Time (s)`, `Mass (g)` (cumulative permeate/balance mass), `Pressure (psi)`,
`Retentate Temp`, **`Retentate Cond @ Temp (uS/cm)`**, `Permeate Temp`,
**`Permeate Cond @ Temp (uS/cm)`**, `Vial Swap` (0/1). A side metadata block has
initial/final stirred-cell weights. A **Vial Data block (cols ~O–V)** has per-vial
`Conductivity @ Temp`, `ICP Salt 1 (mg/L)`, and ICP calibration blocks — i.e. each
experiment carries its own conductivity↔concentration pairs (Feed / Diafiltrate /
Retentate / Final Tube / Vial 1…10).

## Processed target columns (`data/BoE Analysis.xlsx`, e.g. sheet `NF270_MC5 05.27.26_NaCl`)
15 cols; for single salt the populated ones are **B** interfacial conc (mM),
**D** ionic strength (mM), **F** rejection, **H** permeate interfacial conc (mM,
smoothed), **J** Jw (m³/m²·s), **K** B coefficient (µm/s), **M** time, **N**
retentate cond, **O** permeate cond (smoothed). Direct copies: M=raw Time, N=raw
Retentate Cond (verbatim). O = smoothed raw Permeate Cond. B, D, H, J, K computed.

## Recipe (equations + constants) — implement these
- **(a) Water flux** `Jw = (dMass/dt)/(ρ·A_m)` [m³/m²·s], from a **smoothed** slope
  of `Mass(g)` vs `Time(s)` (the reference Jw is smooth/monotone, not raw
  finite-difference). Constants: `ρ ≈ 1000 g/L`, `A_m ≈ 3.95e-4 m²`.
- **(b) Bulk concentration** from raw conductivity via a **linear** calibration
  `c_bulk = (cond − intercept)/slope`.
- **(c) Interfacial concentration** via the thin-film concentration-polarization
  model: `c_int = c_p + (c_bulk − c_p)·exp(Jw·δ/D)`, `δ = D/k`, with
  `k = 0.23·v0^0.57·D^0.67/(ν^0.24·b^0.43)`, `b = 0.0254 m`, `ν = 1.003e-6 m²/s`,
  350 RPM (`ω = 2π·350/60`, `v0 = ω·b/2`), `D_Na = 1.33e-9`, `D_Cl = 2.03e-9`
  (use the appropriate NaCl-salt effective D — document your choice).
- **(d)** `I = ½ Σ zᵢ² cᵢ` (=c for 1:1 NaCl); `R = 1 − H/B = 1 − c_{p,int}/c_int`;
  `B_coeff = Jw·c_{p,int}/(c_int − c_{p,int})` ×1e6 [µm/s].

## WORKING CHOICES + FLAGS (state these in code, README, and the §1 write-up)
These are unresolved in the source material — make a sensible choice, **flag it**,
and let validation adjudicate where possible:
1. **A_m (≈3.95 cm²) and ρ are inferred/undocumented** — parameterize them; note
   the k-correlation's `b=2.54 cm` cell bore is a *different* area than A_m.
2. **Jw smoothing** is undocumented — pick a method (e.g. rolling-window linear
   slope, or fit a smooth monotone mass(t)); tune to match the reference `J` column.
3. **Calibration source.** The embedded MC5 curve is `c=(cond+63.706)/76.685`, but
   the standalone `Multicomponent Calibration Curves.xlsx` (Probe #1) differs
   (slope ~92.7–104). **Preferred approach:** fit a **per-experiment** linear
   calibration from that sheet's own vial (conductivity, ICP) pairs — self-contained
   and probe/session-correct. ICP is `mg/L`; convert to mM with the correct molar
   mass (Na⁺ 22.99 vs NaCl 58.44 — **determine and flag which the ICP column is**).
   Validate: does MC5 05.27.26's vial-fit calibration ≈ `(cond+63.706)/76.685`?
4. **CP vs no-CP.** The main BoE sheet is CP-corrected; the `(2)` sheet skips CP
   (bulk = interfacial). **Note: our Step-2 regression read the `(2)` (non-CP)
   sheet.** Implement **both** paths (a flag `apply_cp: True/False`).
5. H/O "smoothing" is undocumented (H is flattened to ~0.04 early in 05.27.26) —
   reproduce approximately and flag; don't over-engineer.

## (A) Validation-first plan
Validate against the runs where we have BOTH raw and processed:
1. **Exact target — the `(2)` sheet** `NF270_MC5 05.27.26_NaCl (2)` (live formulas,
   **non-CP**, calibration `/76.685`): with `apply_cp=False` and that calibration,
   reproduce **B, D, F, K to ~1e-6 relative**, and Jw as closely as the smoothing
   allows. This is the primary numeric validation of calibration + I/R/B + Jw.
2. **CP targets (approximate)** — main `NF270_MC5 05.27.26_NaCl` and
   `NF270_MC2 05.07.24_NaCl` (CP-corrected, pasted): with `apply_cp=True`, reproduce
   B/H and report the % deviation per column; document residual mismatch (expected,
   from undocumented smoothing/D).
Report a per-column match table for each validation run.

## (A) Process the 6 raw-only runs
Apply the validated pipeline to: MC3 `07.09.24_NaCl`, MC3 `07.22.24_SNaCl`, MC4
`07.11.24_SNaCl`, MC5 `07.23.24_NaCl`, MC5 `07.23.24_SNaCl`, MC5 `07.23.24_S2NaCl`
(all raw sheets confirmed present). Use each run's **own vial calibration**; flag
any run whose vials are too sparse to calibrate (fall back + mark provisional).
Write processed outputs as **CSV** to `analysis/preprocessing/processed/` (one per
run; do **not** modify `data/` or the Box workbooks). Note the `S`/`S2` prefixes are
sequential/reverse (diluting) runs.

## (A) Re-run the multi-dataset regression on the expanded set
Re-run `analysis/nacl_all_datasets/regress_all_nacl.py` (extend it to also read the
new CSVs) across all available NaCl runs. **For comparability with the prior
results, use the non-CP interfacial concentrations** (matching the `(2)` sheet our
Step-2/§2.4 used); note this choice. Report the expanded `(|χ|, δ*)` table — does
adding wide-range runs change the MC2-vs-MC5 picture / sharpen identifiability?

## (B) Refine the report Section 1 preprocessing subsection (equally important)
A subsection **`\subsection{Data preprocessing: from raw signals to the analysis
columns}`** (label `sec:preproc`) **already exists** in Section 1 of
`docs/reports/main.tex` — Opus drafted it with equations (a)–(d), the constants,
and a flagged "Follow-up items (data reduction)" list (A_m/ρ inferred; Jw smoothing;
calibration source + ICP molar mass; CP vs non-CP incl. that §2 used the non-CP
values; H/O smoothing). **Refine it in place — do not rewrite from scratch, and do
not change its label or placement:**
- Confirm the equations match your implementation; correct anything that differs.
- **Update the flagged follow-up items with your validated values**: the fitted
  per-experiment calibration coefficients (and whether MC5's matches `/76.685`),
  the Jw smoothing method you chose, the confirmed-or-still-inferred `A_m`, the ICP
  molar mass used, and whether the CP decision is resolved or still open.
- Keep notation consistent (§1.3 conventions; `\si{}` units, the `\Jw`/`\cint`/`\cp`
  macros). Keep it a methods description (no results tables).
- **Guardrail:** edit only this subsection; do **not** alter other prose, equations,
  tables, numbers, or `\todo`/`\note` flags elsewhere. Compile with `latexmk -pdf`,
  confirm **0 undefined citations, 0 overfull hboxes**, run `latexmk -c` after.
  I (Opus) will review it for notation consistency.

## Deliverables & guardrails
- `analysis/preprocessing/preprocess_nacl.py` (load raw sheet → Jw, calibration,
  optional CP, I/R/B → processed columns), a validation script, and
  `analysis/preprocessing/README.md` (recipe + validation results + flags).
- `analysis/preprocessing/processed/*.csv` for the 6 raw-only runs.
- Updated `analysis/nacl_all_datasets/regress_all_nacl.py` (expanded set).
- The new §1 subsection in `main.tex`.
- Use `analysis/plot_style.py` for any figures (fig_width, below-axes/shared legends).
- Do **not** modify `data/`, `prior_analysis/`, `refs.bib`, or the Box folder.

## Report back
1. Validation: per-column match table for the `(2)` sheet (should be ~exact), and %
   deviations for the CP-corrected main + MC2 runs; what the residual mismatch is.
2. The per-experiment calibrations fit from vial data (and whether MC5's matches the
   embedded `/76.685`); the ICP molar-mass choice.
3. Expanded multi-dataset `(|χ|, δ*)` table across all NaCl runs and what it says
   about reproducibility/identifiability.
4. Confirmation: scripts run in `data3-regression`; the §1 subsection compiles clean
   (0 undefined, 0 overfull); `data/` and `refs.bib` untouched. List every
   working-choice flag you made.
