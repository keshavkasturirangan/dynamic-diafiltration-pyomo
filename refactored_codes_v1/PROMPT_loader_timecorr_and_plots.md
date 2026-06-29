# PROMPT — Generalize the DATA2 time correction + DATA2-styled prediction plots

You are working in the `dynamic-diafiltration-pyomo` repo. Three deliverables, all in
`refactored_codes_v1/`. **Read before you write**: skim the files referenced in each task
and confirm the actual function names / line numbers (they drift) before editing.

## Repository orientation
- Core library: `refactored_codes_v1/refactored_ucb_library.py`
  - Excel/legacy loader: `_load_legacy_data_stru_from_excel(...)` (~line 4780) — builds `data_stru`.
  - Model solve: `solve_model(...)` → returns `(fit, sim, ...)`.
  - Existing comparison plots: `plot_nf270_data_comparison(...)` (~13486) and
    `plot_data_comparison(...)` (~5631).
  - Permeate time-correction metadata already lives on `data_stru["permeate_time_correction"]`.
- Existing prediction-plot generator: `refactored_codes_v1/_make_predictions.py`.
- Conductivity model: `refactored_codes_v1/conductivity_paper.py`
  (`variant_shedlovsky(...)` → single-salt specific conductivity in mS/cm; `msa_transport(...)`
  → multi-salt) and `refactored_codes_v1/CONDUCTIVITY_PAPER_EXPLAINER.md`.
- DATA2 `.mat` time-correction reference: `UnifiedFramework/Trials/Data_loader/load_data.m`
  (MATLAB — see eqns around lines 215–233 for the vial-close / shared-axis convention).
- Python envs: `~/miniforge3/envs/my-idaes-env/bin/python` (pyomo + IPOPT, for solving) and
  `~/miniforge3/bin/python` (plotting). Plotting tasks do not need IPOPT.

## Hard constraints (do not violate)
1. **Do NOT modify any MATLAB `.m` file.** `load_data.m` even has unresolved merge-conflict
   markers — read it for the convention only; never edit it.
2. **Do NOT change DATA1 or DATA2 behavior.** All new behavior must be opt-in / config-driven
   so existing campaigns are byte-for-byte unaffected.
3. Match the **DATA2 plotting convention** (colors + markers) — do not invent a palette.
4. Verify, don't assume: keep the DATA1/DATA2 pytest suite green and visually inspect every
   new/changed figure.

---

## Task 1 — Make the DATA2 time correction a reusable, campaign-level feature

**Problem.** The DATA2 `.mat` workflow applies a time adjustment to the mass, retentate, and
permeate streams (shared time axis; the single ICP permeate value compared at vial close,
`m.cV[vial, tau.last()]`). DATA3/NF270 needs the same correction, but it currently has to be
re-applied/patched by hand each time.

**Do this.** Implement the DATA2 time-correction logic as a **single reusable routine in
`refactored_ucb_library.py`** (e.g. `apply_campaign_time_correction(data_stru, convention=...)`,
or a config flag honored inside `_load_legacy_data_stru_from_excel`), driven by **per-campaign
configuration** so any campaign can enable it declaratively — DATA3 turned on now, others later
— with no manual editing per run.
- Source the convention from `load_data.m` (read-only) and from the existing DATA3 vial-close
  work already in the loader; consolidate them into one code path.
- Record what was applied on `data_stru` metadata (method, anchor, which streams), as is done
  today for `permeate_time_correction`.

**Acceptance.**
- A DATA3 sheet loads with the correction applied automatically via config (no hand-patching).
- DATA1 and DATA2 loads are unchanged (the correction is off / identical for them).
- `python -m pytest pytest_refactored_codes_v1/tests/test_data1_smoke.py
  pytest_refactored_codes_v1/tests/test_data2_paper_comparison.py -q` stays green.
- A short note added to `Architecture.md` describing the new config knob.

---

## Task 2 — Plot-cleanup utility (legend never covers the data; DATA2 colors/symbols)

**Problem.** In the mass-vs-time and concentration-vs-time comparison plots the legend overlaps
the trends being analyzed.

**Do this.** Add a small, reusable styling/cleanup helper (new module
`refactored_codes_v1/_plot_style.py`, or a helper in the library) and apply it to the existing
mass-vs-time and concentration-vs-time plots so:
- The **legend is moved off the data** — place it outside the axes (e.g. `bbox_to_anchor`) or
  auto-select the emptiest quadrant; never occlude the curves.
- A single shared **`DATA2_STYLE`** dict defines colors + markers for: retentate (measurements),
  permeate (measurements), vial (measurements), retentate ICP, permeate/retentate predictions
  (model lines). **Extract these exact colors/markers from the DATA2 plotting code**
  (`plot_data_comparison` and/or the DATA2 reproduction notebook/figures) — do not guess.
- Reuse the same `DATA2_STYLE` everywhere (Task 2 and Task 3) for consistency.

**Acceptance.** Regenerate a representative mass and concentration plot; confirm by rendering to
PNG that the legend sits outside/clear of the data and colors/markers match DATA2.

---

## Task 3 — New per-experiment prediction-fit figure set (5 quantities)

For each experiment (per sheet/run), produce a set of time-series figures comparing
model/prediction to data, reusing `DATA2_STYLE` and the Task-2 legend handling:

1. **Mass per vial vs time** — keep the existing plot as-is (it's good).
2. **ICP concentrations vs time** — retentate and permeate ICP points (separate figure).
3. **Conductivity-probe signal vs time** — raw probe trace for retentate and permeate
   (separate figure). Find the raw conductivity time series in the loader/`data_stru`.
4. **Conductivity-derived concentration vs time** — convert the probe conductivity to
   concentration using `conductivity_paper.py` (**invert** `variant_shedlovsky` per point with a
   root-find such as `scipy.optimize.bisect`, using the salt-specific parameters), plotted for
   retentate and permeate (separate figure). Note in code how this compares to the loader's
   existing conductivity→concentration conversion.
5. **Applied pressure (ΔP) vs time** — separate figure.

Implementation notes:
- Build this as a new generator (e.g. `_make_prediction_panels.py`) or extend
  `_make_predictions.py`; keep it driven by the same run registry used elsewhere.
- One figure per quantity per experiment (separate files), saved under the predictions output
  dir; **also** emit a combined multi-panel overview per experiment if cheap.
- Single-salt inversion is sufficient for the current DATA3 campaign; leave a clear hook for
  multi-salt via `msa_transport` if needed later.

**Acceptance.** Render the 5 figures for one representative DATA3 sheet (e.g. a concentrating
CaCl₂ run); verify each is readable, legend-clear, DATA2-styled, and that the
conductivity-derived concentration tracks the ICP points sensibly.

---

## Task 4 — Per-sheet concentration-range table

For **each experimental data sheet**, generate a table summarizing the **range (min–max)** of
the key concentrations [mM] over the run:
- **Diafiltrate** `c_D`
- **Feed-side retentate** `c_F`
- **Interfacial** `c_in` (membrane surface)
- **Permeate** `c_H`

Notes:
- `c_D`, `c_F`, `c_H` come from the loaded data (`data_stru` / ICP / conductivity-derived
  series); `c_in` is a **model quantity** — take it from the fitted forward-sim
  (`solve_model(...)` → `sim[vial]["cIn"]`). Label each column with its source (measured vs model).
- Report **min and max** per stream (mean optional); keep units consistent (mM).
- Emit **one table per sheet** AND a **combined campaign summary** (one row per sheet) as CSV
  *and* a human-readable markdown table.

**Acceptance.** A range table is produced for every DATA3 sheet; spot-check that `c_F` spans the
dilution/concentration trajectory and that `c_in` sits on the expected side of `c_F` given the
polarization direction (concentrating: `c_in > c_F`; diluting: `c_in < c_F`).

---

## Task 5 — WSSE identifiability contours for the new response channels

Extend the existing per-response **WSSE contour** analysis (σ×B log-WSSE surfaces at fixed Lₚ —
the same method and visual style as the current mass / retentate / permeate contours) to the
**new response channels** from Task 3. Produce one contour per channel:
- **Mass per vial** (as today).
- **ICP retentate + permeate concentration.**
- **Conductivity-probe signal** (retentate + permeate) — WSSE computed **in conductivity space**:
  map the model concentration → predicted conductivity via
  `conductivity_paper.variant_shedlovsky(...)` and compare to the measured probe trace.
- **Conductivity-derived concentration** (retentate + permeate) — WSSE computed **in
  concentration space**: invert the probe conductivity → concentration (Task 3 inversion) and
  compare to the model concentration.

Implementation:
- **Reuse the existing contour machinery** — `_run_Bsigma_fixedLp.py`, `calc_contour_2d_py` /
  `calc_contour_3d_py`, `calc_FIM`, and the per-response WSSE objective branches in
  `refactored_ucb_library.py` — plus the established contour styling (the "v7" look with the
  optimum-marker box, cf. `_replot_Bsigma_titled.py` / `_rerender_contours_with_optimum_box.py`).
- The conductivity channels add **two new residual definitions** (conductivity-space and
  concentration-space); wire them in as additional per-response objective channels without
  altering the existing mass/retentate/permeate channels.
- **ΔP (applied pressure) is a control input, not a fitted response — no WSSE contour for it.**

**Acceptance.** A log-WSSE σ×B contour rendered per response channel for a representative DATA3
sheet, in the same style as the existing contours, with the optimum marker shown; confirm the
per-channel identifiability pattern (narrow valley in B, ~flat in σ) is consistent with the
established single-experiment finding.

---

## Task 6 — Rejection vs log(Pe) diagnostic (σ-family) for choosing σ per experiment (analysis phase)

Build an analysis-phase plot of **rejection R (y-axis) vs log(Pe) (x-axis)**, overlaying **one
trajectory per σ value** so the best σ to use for a given experiment can be read off the figure.
- Each trajectory is the Spiegler–Kedem rejection curve at a fixed σ:
  `R(Pe; σ) = σ(1 − e^−Pe) / (1 − σ·e^−Pe)`, swept over log(Pe) for a grid of σ (e.g. 0.1 … 0.99).
  Each curve rises from `R → 0` at `Pe → 0` to a plateau `R → σ` at `Pe ≫ 1`.
- **Overlay the experiment's observed/fitted rejection at its operating Pe** (use Pe / Pe* from
  the Péclet analysis); the σ-curve passing through that point is the best-supported σ — mark it.
- Purpose: a visual, analysis-phase way to select σ per experiment and to see how tightly σ is
  constrained (ties into the established single-experiment σ–B identifiability finding — when the
  operating Pe is large, many σ-curves bunch near their plateau, so σ is poorly distinguished).
- Reuse the Pe definition/sweep from Tasks 1/5; use the established plot styling + the Task-2
  legend handling (legend off the data; a labeled σ per curve).

**Acceptance.** An R-vs-log(Pe) figure with the σ-family of curves overlaid for a representative
DATA3 sheet, the experiment's operating point marked, and the best-fit σ identifiable.

---

## Task 7 — Fix the diafiltrate feed rates (S0, S) from recorded values for DATA3 (lag mode)

**Context.** All DATA3/NF270 runs are **lag mode**. In lag mode the cell mass balance
(`ode_mF_rule`, ~line 2029 of `refactored_ucb_library.py`) is
`dmF/dt = −S0 − area·ρ·Jw` on the active concentration legs and `dmF/dt = S/3600` on the final
closed leg. Today `S0` is a fixed `Param` but **seeded from an `M_O`-based guess**, and `S` is a
**free `Var`** seeded from `M_O/Σ(tf−ti)·3600` (~lines 1810, 1956–1958). So the apparatus feed
terms are guessed/estimated rather than taken from the experiment — even though these flow rates
are **recorded per run**. They should be fixed inputs, not fitted.

**Solution (implement this):**
1. **Ingest the recorded feed rates** S0 (g/s) and S (g/hr) per experiment into `data_config`
   (from the run logs / Excel) under explicit keys, e.g. `feed_rate_S0_g_per_s`,
   `feed_rate_S_g_per_hr`.
2. In `model_construct_inter`, **for the DATA3 lag path only**, set BOTH `m.S0` and `m.S` as
   fixed `Param`s from those `data_config` values — replace the `M_O`-based seed for `S0` and the
   `m.S = Var(...)` for `S` (~lines 1950–1958). Guard with `workflow_family == "DATA3"`; leave the
   DATA1/DATA2 and Overflow branches byte-identical.
3. **Fall back** to the current behavior if a recorded value is missing (don't break loading),
   and record on `data_stru` metadata which source was used (recorded vs fallback).

**Why.** Removes the only free apparatus degree of freedom (`S`), so the fit estimates just the
transport parameters of interest (Lₚ, the B(c) coefficients, σ). It also eliminates a phantom
free variable on any run with no final closed leg (where `S` never enters the balance and would
otherwise drift unconstrained).

**Caveat.** This is a correctness / parameter-count cleanup, **not** an σ-identifiability fix — σ
stays flat for the transport-regime reason (high-Pe plateau). Expect fitted Lₚ/B/σ to be ~unchanged
(S was already pinned by the mass channel), just with one fewer free parameter.

**Acceptance.** A DATA3 lag sheet fits with S0/S fixed from recorded values (metadata shows the
source); the free-parameter count drops by one; DATA1/DATA2 results are byte-identical; pytest green.

---

## Task 8 — Square (forward-solve) WSSE contour surfaces for every B(c) form × every DATA3 single-salt experiment

**Objective.** For **each DATA3 single-salt experiment** (all sheets) and **each B(c) correlation**
— constant / linear / quadratic / cubic / saturating-exp / donnan — generate the WSSE contour
surface for the **response channels of Task 5**, computed as a **square forward solve** at every
grid node (no per-node nonlinear optimization).

**How every node stays square (the decision rule):**
1. Fit each B-form **once** per experiment (or read its `result_<form>.json`) → the best-fit
   (maximum-likelihood) parameter set θ\*.
2. Hold **fixed across the whole grid**: `Lp` (its optimum), `S0` and `S` (recorded feed rates,
   per Task 7), and **every B(c) coefficient that is NOT on a contour axis**, pinned at θ\*.
3. Axes = `σ` × `[one B parameter]` (e.g. the leading coefficient / in-window B level). With two
   axes setting two unknowns and everything else pinned, the model has **zero free parameters →
   forward-solve the DAE (square) and *evaluate* the WSSE** at each node.
4. This is the **slice** convention — a conditional landscape *through the optimum* (its minimum
   coincides with the fitted (σ, B)). Label it as conditional, not a profile.

**Response channels (per Task 5), one contour panel each:** mass; ICP retentate conc; ICP permeate
conc; conductivity-probe (retentate + permeate, conductivity space); conductivity-derived
concentration (retentate + permeate, concentration space). ΔP is excluded (control input).

**Implementation.**
- Reuse `calc_contour_2d_py` / the `method="slice"` path of `_run_profile_contours.py` plus the
  established contour styling (optimum marker); extend the channel set to the conductivity ones.
- Every node is a CPU-capped forward solve (ideally with a *bounded* CasADi/IDAS init —
  `integrator_options={'max_num_steps': …}` — or `skip_sim_init`), so it runs fast and avoids the
  per-node optimization holes the profile campaign hit.
- Optional robustness check: regenerate a slice at the MLE ± standard error of the pinned
  coefficients to confirm the conclusion is insensitive to the pin.

**Deliverable.** One multi-panel contour figure per **(experiment × B-form)** — a panel per
response channel — covering all DATA3 single-salt sheets × all six forms, saved under
`bform_study/` with a small index.

**Acceptance.** Square forward-solve contours produced for every (experiment, B-form, channel);
spot-check that each surface's minimum sits at the fitted (σ, B) and that σ is ~flat (the expected
non-identifiability), consistent across forms and salts.

---

## Task 9 — Profile-likelihood WSSE contour surfaces (rigorous counterpart of Task 8)

**Objective.** Same scope as Task 8 — every DATA3 single-salt experiment × every B(c) form × the
Task-5 response channels — but compute the **profile-likelihood** surface: at each grid node, **do
NOT pre-fix the non-axis coefficients — re-optimize them.**

**Method (per node).** Clamp the two axes (σ and the chosen B parameter), hold `Lp`, `S0`, `S`
fixed (`Lp` optimum; `S0`/`S` recorded, per Task 7), then solve a sub-optimization over the
remaining free parameters λ (the non-axis B coefficients):

```
Phi_profile(psi)  =  min over lambda  of  WSSE(psi, lambda | Lp, S0, S fixed)
```

Record `Phi_profile` as the pixel value, AND store `lambda*(psi)` (the profiled nuisance values).

**Why (vs Task 8).** Task 8 pins λ at the MLE → a *conditional* (slice) surface, which is **not** a
valid confidence region. The profile re-optimizes λ → the surface supports **likelihood-ratio
confidence regions**:

```
CR  =  { psi :  Phi_profile(psi) - Phi_min  <=  chi-square threshold }
```

Use it when the goal is honest identifiability / confidence intervals on σ and the B parameters.

**By-product to plot.** The profiled trajectories `lambda*(psi)` (e.g. β₀ as σ varies) directly
visualize the parameter correlations (the β₀–β₁ trade-off; whether β₁ tracks σ).

**Cost & machinery (NOT square).** This is a nonlinear optimization at every node — you cannot keep
it a forward solve. Reuse the existing `_run_profile_contours.py` `method="profile"` campaign:
- warm-start continuation (seed each node from a converged neighbour — a *guess*, not a fix),
- per-node IPOPT CPU cap + a bounded / `skip_sim_init` CasADi init,
- watchdog kill-restart + JSONL resume for stiff nodes (expect occasional holes).
Task 7 (fixed S0/S) shrinks λ, so each per-node solve is faster and more robust.

**Deliverable.** One multi-panel profile contour per (experiment × B-form) — a panel per response
channel — with the LR confidence regions overlaid (optionally the profiled-λ trajectories); all
DATA3 single-salt sheets × all six forms.

**Acceptance.** Profile surfaces produced for every (experiment, B-form, channel); each minimum
coincides with the fit; LR confidence regions drawn; the flat-σ non-identifiability is confirmed
**rigorously** (not just conditionally), and any difference from the Task-8 slice is noted.

---

## Task 10 — Dimensional-consistency (units) audit + unit-change sensitivity test

**Objective.** Verify the whole DATA3 model + analysis pipeline is **dimensionally consistent**
(no silent unit mismatches), document the unit of every variable/parameter/constant, and add a
small **sensitivity test** that re-expresses inputs in different units (with the matching
conversion) to confirm the physics is invariant.

**Part A — units audit.**
1. Build a **units table** for every model quantity: mass (g), `S0` (g/s), `S` (g/hr),
   concentrations `cF/cV/cD/cIn/cH` (mM), `Jw` (state which — LMH from `Lp·ΔP` vs velocity cm/s),
   `Lp` (L·m⁻²·h⁻¹·bar⁻¹), `ΔP` (bar), `Δπ` (bar, van 't Hoff `ni·R·T·Δc`), `B`, `σ` (–),
   `D` (cm²/s), the mass-transfer coeff `k` (cm/s), `Pe` (–), area `Am`, density `rho`, time
   (min vs s), and the conductivity domain (mS/cm, µS/cm; M vs mM in Shedlovsky/MSA).
2. Check each equation **term-by-term**, paying special attention to the cross-unit interfaces
   that are easy to get wrong:
   - **`Jw` consistency:** it appears as a velocity in `exp(Jw/k)`, in `Am·rho·Jw` (a mass rate),
     and in `Pe = Jw·L/Dm`. Confirm `Jw` is converted to the *same* velocity units everywhere
     (`Lp·ΔP` yields LMH — flag any spot that uses LMH where cm/s is required).
   - **`B` has two conventions:** single uses `Js = B·(cIn−cH)`; the B(c) forms use
     `Js = Jw·B·(cIn−cH)` → **different units for `B`**. Confirm the apparent-µm/s conversion
     (`×Jw·1e4`) is applied consistently and the contour/overlay axes use one convention.
   - **`ΔP` vs `Δπ`** pressure units in `Jw = Lp(ΔP − σ·Δπ)` (van 't Hoff mM→pressure; bar vs Pa).
   - **`S0` (g/s) vs `S` (g/hr)** — the `/3600`; confirm the time-scaling `(TF−TI)/Tauf` is
     dimensionless and the time axis (min vs s) is consistent.
   - **Conductivity:** Shedlovsky/MSA input concentration units (M vs mM), outputs (mS/cm vs
     µS/cm), and the EC25 temperature compensation.
3. Document findings; fix any genuine mismatch (guard DATA3-only if it touches the shared model;
   never silently change DATA1/DATA2 numerics — flag those).

**Part B — unit-change sensitivity test.**
- For a representative DATA3 sheet: run the fit, then **re-run with one input in different units +
  the correct conversion** (e.g. `ΔP` in Pa with `Lp` converted; or concentrations in M not mM;
  or `Jw`/area in SI). The fitted `Lp/B/σ` (back in consistent units) and the WSSE **must be
  invariant.** Any change exposes a hidden hard-coded unit assumption (a bug).
- Add a **dimensional scaling check**: multiply an input by a known factor and confirm the outputs
  scale exactly as dimensional analysis predicts.
- Optionally wire 1–2 of these as **pytest invariants** (DATA3-guarded) so unit regressions are
  caught automatically.

**Acceptance.** A units table + term-by-term consistency report committed (e.g. a `UNITS.md` or a
section in `Architecture.md`); the unit-change test passes (results invariant under correct
conversion) for ≥1 DATA3 sheet; any mismatch found is documented and (if DATA3-only) fixed.

---

## Verification checklist (run at the end)
- `pytest` DATA1 + DATA2 smoke/paper-comparison green (DATA1/DATA2 untouched).
- Per-sheet concentration-range table (Task 4) emitted as CSV + markdown for each DATA3 sheet.
- WSSE σ×B contours (Task 5) rendered per response channel for a representative sheet, matching
  the established contour style.
- Rejection-vs-log(Pe) σ-family diagnostic (Task 6) rendered for a representative sheet, with the
  experiment's operating point marked and best-fit σ identifiable.
- Diafiltrate feed rates S0/S fixed from recorded values for DATA3 lag fits (Task 7); free-param
  count reduced by one; DATA1/DATA2 byte-identical.
- Square forward-solve WSSE contour surfaces (Task 8) generated for every DATA3 single-salt
  experiment × every B(c) form × every Task-5 response channel, with non-axis coefficients pinned
  at the per-form MLE; each surface's minimum coincides with the fit.
- Profile-likelihood WSSE contour surfaces (Task 9) generated per experiment × B-form × channel
  with non-axis coefficients re-optimized per node; LR confidence regions drawn; flat-σ confirmed
  rigorously vs the Task-8 slice.
- Units audit (Task 10): units table + term-by-term consistency report committed; unit-change
  sensitivity test passes (fit + WSSE invariant under correct conversion) for ≥1 DATA3 sheet.
- A DATA3 sheet loads with the time correction on via config; metadata recorded.
- All regenerated/new figures rendered to PNG and visually inspected (legend clear of data,
  DATA2 colors/markers).
- `Architecture.md` updated with the new time-correction config + plotting helpers.

## Assumptions baked in (change if wrong)
- Task 3 figures are **separate files per quantity per experiment** (+ optional combined panel).
- Conductivity→concentration uses **single-salt** `variant_shedlovsky` inversion for now.
- "Concentration vs time" plots cover both retentate and permeate on shared axes where it aids
  comparison; split if it gets crowded.
