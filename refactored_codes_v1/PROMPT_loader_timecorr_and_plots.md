# DATA3 PROMPT DOC — loader, plots, contours, and the physics push

This is the working task list for the DATA3/NF270 diafiltration analysis in the
`dynamic-diafiltration-pyomo` repo (branch `refactored_unified_codebase`).
It is written to be read **one card at a time** — every task is self-contained.

---

## How to read this file (30 seconds)

- **Every task is one card** with the same five slots:
  **Goal** (one line) → **Why** → **Do** (checklist) → **Done when** → **Files**.
- **Status chips:** ✅ = done · 🔶 = partly done · ⬜ = not started.
  The dashboard below tells you where everything stands at a glance.
- **Part A** (Tasks 1–11) = the original loader/plots/contours work.
  **Part B** (Tasks 12–23) = added 2026-07-22. Task 11+k = point k of the
  2026-07-22 request, so the numbering maps 1:1 to that message.
- Line numbers drift. **Skim the referenced file and confirm names/lines before
  editing.**
- **This `.md` is the ONLY copy** — the `.docx`/`.pdf` exports were retired
  2026-07-22. Edit here (LaTeX `$…$` math and `![](path.png)` images welcome);
  need a Word/PDF to share? See **Sharing** at the very bottom.

---

## Status dashboard

| # | Task (short name) | Status | Output / where it lives |
|--|--------------------|----|----------------------------------|
| 1 | Config-driven DATA2 time correction | ✅ | `apply_campaign_time_correction` in the library; `Architecture.md` note |
| 2 | Legend-off-data + shared DATA2 style | 🔶 | applied in new figures; no single reusable `DATA2_STYLE` helper yet |
| 3 | Per-experiment prediction figures (8 plots) | 🔶 | 5-panel combined version in library (`plot_meas_vs_pred_conductivity`); 8 separate files + ΔP plot pending |
| 4 | Per-sheet concentration-range table | ⬜ | (`concentration_range.png` exists as a one-off) |
| 5 | 7-channel σ×B WSSE slice contours | ✅ | `_run_Bsigma_channels.py` → `bform_study/Bsigma_channels/` |
| 6 | Rejection R vs log Pe (σ-family) | ⬜ | superseded/extended by Task 14 |
| 7 | Fix S0/S from recorded feed rates | ⬜ | — |
| 8 | Square slice contours: all B-forms × all sheets | 🔶 | `_run_Bsigma_coeffs.py` → `bform_study/Bsigma_coeffs/` + A4 PDFs (NaCl set + representatives done) |
| 9 | Profile-likelihood contours (rigorous) | 🔶 | `_run_profile_contours.py` machinery; full campaign + LR regions pending |
| 10 | Units audit + unit-change sensitivity | ⬜ | — |
| 11 | Selectable cond↔conc converter (3 backends) | ✅ | `concentration_to_conductivity` / `conductivity_to_concentration`, runfile branch `v` |
| 12 | Pull Alex's regression repo + notes PDF | 🔶 | imported → `data3_regression_july26_AD/`; the *reading/reconciling* is pending |
| 13 | Consolidate ALL plotting into the library | 🔶 | meas-vs-pred ported (the worked example); contour drivers et al. pending |
| 14 | B(c_int) vs log Pe rejection plots (σ-trajectories) | ⬜ | — |
| 15 | Derive B–c_int correlation per basis function | ⬜ | — |
| 16 | The "beyond-Pe" evidence plot (6 representative runs) | ⬜ | — |
| 17 | Derive the replacement dimensionless number | ⬜ | — |
| 18 | Analysis with the new number × B(c_int) | ⬜ | — |
| 19 | Js/Jw vs Pe and vs the new number (2×2) | ⬜ | — |
| 20 | Electromigration: hammer vs scalpel | ⬜ | — |
| 21 | Diffusivity-coefficient audit (cations/anions) | ⬜ | — |
| 22 | Integrate the BoE coupled-diffusivity equation | 🔶 | Box archive imported → `box_multicomponent_diafiltration/`; the port is pending |
| 23 | Sensitivity-analysis plots as library functions | ⬜ | (`sigma_sensitivity-*.png` exist as one-offs) |

**Suggested order for Part B** (cheap → heavy):
1. **Task 12** (read Alex's notes — it reframes everything else),
2. **Tasks 20 + 21** (paper-and-pencil physics + a value audit; no big compute),
3. **Tasks 16 → 17 → 18 → 19** (the Pe story, in that order — each feeds the next),
4. **Tasks 14 + 15** (rejection plots + derivations),
5. **Tasks 13 + 23** (library consolidation + sensitivity functions),
6. then the heavy campaigns: Tasks 7, 4, 10, 8, 9.

---

## Ground rules (never violate)

1. **Do NOT modify any MATLAB `.m` file.** (`load_data.m` even has unresolved
   merge-conflict markers — read for the convention only.)
2. **Do NOT change DATA1 or DATA2 behavior.** Everything new is opt-in /
   config-driven; existing campaigns stay byte-for-byte identical.
3. Match the **DATA2 plotting convention** (colors + markers) — never invent a palette.
4. **Verify, don't assume:** keep the DATA1/DATA2 pytest suite green and visually
   inspect every new/changed figure.
5. Do not edit `conductivity_paper.py` — wrap it from the library instead.
6. `box_multicomponent_diafiltration/` and `data3_regression_july26_AD/` are
   **read-only archives** — copy out what you need.

## Key files & environments

- **Core library (put reusable code HERE):** `refactored_codes_v1/refactored_ucb_library.py`
  — loader `_load_legacy_data_stru_from_excel`, solver `solve_model(...)`
  (`sim_opt=True` = square forward solve, `False` = fit), time correction
  `apply_campaign_time_correction`, converters
  `concentration_to_conductivity` / `conductivity_to_concentration`,
  meas-vs-pred `plot_meas_vs_pred_conductivity` + `run_nf270_meas_vs_pred`.
- **Interactive driver:** `refactored_codes_v1/refactored_ucb_runfile.py`
  (branch letters; `v` = meas-vs-pred is the wiring pattern to copy).
- **Conductivity physics (do not edit):** `refactored_codes_v1/conductivity_paper.py`
  (`variant_shedlovsky`, `msa_transport`) + `CONDUCTIVITY_PAPER_EXPLAINER.md`.
- **Contour drivers (to be absorbed by Task 13):** `_run_Bsigma_channels.py`,
  `_run_Bsigma_coeffs.py`, `_run_Bsigma_fixedLp.py`, `_run_profile_contours.py`,
  `_run_all_pairs_5channel.py`, `_run_all_11_sheets_5channel.py`.
- **Figure output root:** `UnifiedFramework/DATA3/results/paper_artifacts/nf270/`.
- **Alex's regression snapshot:** `data3_regression_july26_AD/`
  (notes = `DATA3_regression_notes_AD_2026-07-07.pdf`; see its `PROVENANCE.md`).
- **Box archive:** `box_multicomponent_diafiltration/` (BoE derivation + MATLAB;
  see its `PROVENANCE.md`).
- **Envs:** `~/miniforge3/envs/my-idaes-env/bin/python` (pyomo + IPOPT — anything
  that solves) · `~/miniforge3/bin/python` (plot-only).
- **DATA2 time-correction reference (read-only):**
  `legacy/data1_matlab/functions/load_data.m` + DATA2 notebooks.

---

# PART A — original tasks (1–11)

---

## Task 1 — Config-driven DATA2 time correction  ✅ DONE

**Goal.** One reusable, campaign-level routine for the DATA2 permeate time
correction — no more hand-patching per run.

**What was built.**
- `apply_campaign_time_correction(data_stru, convention=None, *, workflow_family="DATA3")`
  in the library (+ helper `_permeate_anchor_index`).
- Conventions: `"vial_close"` (default; byte-identical to the old inline logic),
  `"tube_transit"` (retired V_tube model), `"none"`.
- Config: global `CAMPAIGN_TIME_CORRECTION = {"DATA3": "vial_close"}` +
  per-run override `data_config['permeate_time_correction_convention']`.
- The Excel loader calls it automatically and stamps
  `data_stru["permeate_time_correction"]` (method, anchor, streams).
- DATA1/DATA2 (`.mat`) never reach it → untouched. `Architecture.md` documents
  the knob ("Permeate time-correction convention").

**Verified by.** DATA1 smoke + DATA2 paper-comparison pytest green; DATA3 sheets
load with the correction on via config.

---

## Task 2 — Legend never covers data + one shared DATA2 style  🔶 PARTLY DONE

**Goal.** A single reusable styling helper so every figure uses the DATA2
palette and the legend never occludes the trends.

**Why.** Legends were sitting on top of the curves; colors were re-derived
ad hoc per script.

**Done so far.** Legend-off-data placement applied to the coefficient-contour
and meas-vs-pred figures. The DATA2 palette is encoded inside
`plot_meas_vs_pred_conductivity`: mass `r.` (meas) / `b-` (model); retentate
`ms` (meas) / `g-` (model); permeate vial `cs` + `r^` (cV) / `r-` (cH);
backends by line style (— Shedlovsky, ‑‑ MSA, ·· DATA2).

**Do (remaining).**
- [ ] Extract ONE shared `DATA2_STYLE` dict (colors + markers per stream) into
  the library and route the existing plots through it.
- [ ] A small `place_legend_clear(ax)` helper (outside axes or emptiest corner).

**Done when.** A representative mass + concentration figure regenerates with the
shared dict; nothing overlaps the data; colors match DATA2 exactly.

---

## Task 3 — Per-experiment prediction-figure set (8 separate plots)  🔶 PARTLY DONE

**Goal.** For each DATA3 sheet: **8 separate** measurement-vs-prediction
time-series figures (retentate and permeate always split), DATA2-styled, on the
DATA2-corrected time axis.

1. Mass per vial (keep as-is — per-vial reset)
2. ICP retentate concentration
3. ICP permeate concentration
4. Conductivity probe — retentate (raw trace)
5. Conductivity probe — permeate
6. Conductivity-derived concentration — retentate (via Task-11 converter, per point)
7. Conductivity-derived concentration — permeate
8. Applied pressure ΔP (control input)

**Done so far.** The library's `plot_meas_vs_pred_conductivity` +
`run_nf270_meas_vs_pred` (runfile branch `v`) produce a combined **5-panel**
version of plots 1/4/5/6/7 with all three converter backends overlaid, with the
full DATA2 time convention (shared `t_delay` origin; model drawn only for
collection vials `i ≥ n_v0−1` — the startup-holdup vial is *not* predicted;
per-vial mass reset; permeate ICP anchored at vial close).

**Do (remaining).**
- [ ] Emit the **8 individual figure files** per sheet (add the ICP-only panels
  and the ΔP plot), reusing the same library internals — do it inside the
  library, not a new script (Task 13 rule).
- [ ] Keep the combined panel as the cheap overview.

**Done when.** All 8 files render for one representative sheet (e.g. a
concentrating CaCl₂ run), legend-clear, DATA2-styled, holdup handled.

---

## Task 4 — Per-sheet concentration-range table  ⬜ NOT STARTED

**Goal.** One table per sheet + one campaign summary: **min–max [mM]** of
diafiltrate `c_D`, feed/retentate `c_F`, interfacial `c_in`, permeate `c_H`.

**Why.** The B(c) story needs to know what concentration window each experiment
actually explores.

**Do.**
- [ ] `c_D`, `c_F`, `c_H` from the loaded data; `c_in` from the fitted forward
  sim (`solve_model(...)` → `sim[vial]["cIn"]`) — label measured vs model.
- [ ] Emit per-sheet CSV **and** markdown, plus one combined campaign table
  (a row per sheet).

**Done when.** Every DATA3 sheet has a table; spot-check `c_in > c_F` while
concentrating and `c_in < c_F` while diluting.

---

## Task 5 — 7-channel σ×B WSSE slice contours  ✅ DONE

**Goal.** Extend the σ×B log-WSSE contour (fixed Lₚ, square forward solve at
every node — the *slice*, NOT profile likelihood) to all seven response
channels, one contour each.

**What was built.** `_run_Bsigma_channels.py` (`exp <rid>` / `all` / `plot <rid>`):
mass · ICP retentate · ICP permeate · conductivity-probe ret/perm (σ-space) ·
cond-derived concentration ret/perm (c-space). ΔP excluded (control input).
Reuses the `Bsigma_fixedLp_nofloor` grid; channels go through the Task-11
converter. **Gotcha handled:** the permeate-ICP channel needs absolute
measurement-error floors (`FLOOR_MM=0.5`, `FLOOR_US=50`) or the near-zero
startup-vial ICP (~0.08 mM) dominates the WSSE.

**Output.** `bform_study/Bsigma_channels/<rid>/single/{contour7.png,
contourdata-7ch.csv, nodes.jsonl, grid.json}` for the three representative
sheets (CaCl₂ MC2.05.07.24, LaCl₃ MC2.05.21.24, NaCl MC2.05.07.24).

---

## Task 6 — Rejection R vs log Pe with the σ-family  ⬜ (extended by Task 14)

**Goal.** Analysis-phase figure: Spiegler–Kedem rejection curves
`R(Pe; σ) = σ(1 − e^−Pe) / (1 − σ·e^−Pe)` for a grid of σ (0.1…0.99) vs
log Pe, with each experiment's operating point overlaid — read the
best-supported σ off the figure.

**Why.** Visualizes *why* σ is weakly identified: at large operating Pe all
σ-curves sit at their plateau `R → σ`, so tiny R noise moves σ a lot.

**Do.** See **Task 14**, which wraps this and the B(c_int)-vs-log-Pe variant
into one deliverable (same machinery, one card).

---

## Task 7 — Fix S0/S from recorded feed rates (DATA3 lag mode)  ⬜ NOT STARTED

**Goal.** Stop guessing the apparatus feed terms: set `S0` (g/s) and `S` (g/hr)
from the **recorded per-run values**, as fixed Params — not seeded guesses, not
a free Var.

**Why.** Today `S0` is seeded from an `M_O`-based guess and `S` is a free `Var`
(~lines 1810, 1956–58; lag balance `ode_mF_rule` ~line 2029:
`dmF/dt = −S0 − area·ρ·Jw` on active legs, `dmF/dt = S/3600` on the closed
leg). That's one phantom degree of freedom — and unconstrained on runs with no
closed leg. This is a correctness cleanup, **not** a σ-identifiability fix
(σ stays flat for the high-Pe-plateau reason).

**Do.**
- [ ] Ingest recorded rates into `data_config` under explicit keys
  (`feed_rate_S0_g_per_s`, `feed_rate_S_g_per_hr`).
- [ ] In `model_construct_inter`, DATA3-lag path only
  (`workflow_family == "DATA3"`): make BOTH `m.S0` and `m.S` fixed Params from
  those values. DATA1/DATA2/Overflow branches byte-identical.
- [ ] Fall back to current behavior when a recorded value is missing; stamp the
  source (recorded vs fallback) on `data_stru` metadata.

**Done when.** A DATA3 lag sheet fits with S0/S fixed (metadata shows source);
free-parameter count drops by one; DATA1/DATA2 byte-identical; pytest green.
Expect fitted Lₚ/B/σ ≈ unchanged (mass channel already pinned S).

---

## Task 8 — Square slice contours: every B(c) form × every sheet  🔶 PARTLY DONE

**Goal.** For each DATA3 single-salt sheet × each B(c) form (constant / linear /
quadratic / cubic / saturating-exp / donnan): WSSE contours for the Task-5
channels, every node a **square forward solve** (no per-node optimization).

**The decision rule that keeps every node square.**
1. Fit each form once (or read `result_<form>.json`) → θ\*.
2. Pin across the grid: `Lp` (optimum), `S0`/`S` (recorded — Task 7), and every
   B-coefficient NOT on an axis, at θ\*.
3. Axes = σ × one B-parameter → zero free parameters → forward-solve the DAE and
   *evaluate* WSSE.
4. Label it a **conditional slice through the optimum** — its minimum coincides
   with the fit; it is *not* a profile.

**Done so far.** `_run_Bsigma_coeffs.py`: per-coefficient (β_k × σ) slices for
linear/quad/cubic/sat (B built on interfacial conc `cIn`, `B_form_rule3`),
subprocess worker + watchdog + JSONL resume + a cheap B(cIn)-physicality
pre-filter; β₀×β₁ correlation slice for NaCl; A4 multi-page PDFs
(`<rid>_coeff_contours_A4.pdf`). Library toggle **`NF270_BFORM_SLICE_UNBOUNDED_B`**
(default False = byte-identical) because the poly `sim_opt` `m.B` is a bounded
Var that's infeasible off-optimum. **Honest findings:** slices are feasible only
near each optimum (off-optimum B(c) leaves the physical region → infeasible
holes, themselves informative); **cubic is essentially un-sliceable for NaCl**
(fitted B(cIn) goes negative over the run) but fine for CaCl₂/LaCl₃.

**Do (remaining).**
- [ ] Extend to **all** DATA3 single-salt sheets × **all six** forms × all 7
  channels; save under `bform_study/` with a small index.
- [ ] Optional robustness: re-slice at MLE ± SE of the pinned coefficients.

**Done when.** Every (sheet, form, channel) has a slice; each minimum sits at
the fitted (σ, B); the flat-σ pattern is consistent across forms and salts.

---

## Task 9 — Profile-likelihood contours (rigorous counterpart)  🔶 PARTLY DONE

**Goal.** Same scope as Task 8, but at each node **re-optimize the non-axis
coefficients** (nuisances λ):
`Φ_profile(ψ) = min_λ WSSE(ψ, λ | Lp, S0, S fixed)` — the surface that supports
likelihood-ratio confidence regions
`CR = {ψ : Φ_profile(ψ) − Φ_min ≤ χ² threshold}`.

**Why.** The Task-8 slice is conditional — NOT a valid confidence region. The
profile is the honest identifiability statement. By-product: the profiled
trajectories `λ*(ψ)` directly visualize parameter correlations (β₀–β₁
trade-off; whether β₁ tracks σ).

**Done so far.** `_run_profile_contours.py` has the `method="profile"` campaign
machinery (warm-start continuation from converged neighbours, per-node IPOPT CPU
cap, `skip_sim_init` / bounded CasADi init, watchdog kill-restart, JSONL
resume); earlier per-form profile contours exist for representative sheets.

**Do (remaining).**
- [ ] Full campaign: every sheet × form × channel, with LR confidence regions
  overlaid (+ optionally the λ*(ψ) trajectories).
- [ ] Note differences vs the Task-8 slice. (Task 7 shrinks λ → faster nodes.)

**Done when.** Profile surfaces exist per (sheet, form, channel); minima
coincide with fits; LR regions drawn; flat-σ confirmed **rigorously**.

---

## Task 10 — Units audit + unit-change sensitivity test  ⬜ NOT STARTED

**Goal.** Prove the pipeline is dimensionally consistent; document every unit;
add a test that re-expresses inputs in different units and demands identical
physics.

**Part A — audit.**
- [ ] Units table for every quantity: mass g · `S0` g/s · `S` g/hr ·
  `cF/cV/cD/cIn/cH` mM · `Jw` (LMH vs cm/s — state which where) ·
  `Lp` L·m⁻²·h⁻¹·bar⁻¹ · ΔP bar · Δπ bar (van 't Hoff `n_i·R·T·Δc`) · `B` ·
  σ (–) · `D` cm²/s · `k` cm/s · Pe (–) · `Am`, ρ, time (min vs s) ·
  conductivity mS/cm vs µS/cm, M vs mM.
- [ ] Term-by-term checks at the risky interfaces:
  `Jw` as velocity in `exp(Jw/k)` vs mass rate in `Am·ρ·Jw` vs `Pe = Jw·L/D`;
  **B's two conventions** (`Js = B·Δc` for single vs `Js = Jw·B·Δc` for the
  B(c) forms — different units; the `×Jw·1e4` apparent-µm/s conversion);
  ΔP vs Δπ pressure units; the `S/3600`; `(TF−TI)/Tauf` dimensionless;
  Shedlovsky/MSA input/output units + EC25 compensation.
- [ ] Fix genuine mismatches DATA3-guarded; never silently change DATA1/DATA2.

**Part B — sensitivity test.**
- [ ] Re-run one fit with an input in different units + correct conversion
  (ΔP in Pa with Lp converted; conc in M; SI area). Fitted Lp/B/σ and WSSE
  **must be invariant** — any drift is a hidden hard-coded unit (a bug).
- [ ] A dimensional-scaling check (multiply input by known factor → outputs
  scale exactly as predicted). Wire 1–2 as DATA3-guarded pytest invariants.

**Done when.** `UNITS.md` (or an `Architecture.md` section) with the table +
report; invariance test passes on ≥1 sheet.

---

## Task 11 — Selectable conductivity ↔ concentration converter  ✅ DONE

**Goal.** One converter pair with a selectable backend so every
conductivity-derived quantity can be computed under — and compared across —
three models.

**What was built.** In the library (wrapping `conductivity_paper.py`, unedited):
- `concentration_to_conductivity(c_mM, salt, T=25, *, method=...)` and
  `conductivity_to_concentration(cond_uS, salt, T=25, *, method=...)`,
  `method ∈ {"shedlovsky", "msa", "data2"}`.
- `"shedlovsky"` = `variant_shedlovsky` (inverse via per-point bisection;
  NaN-safe where the map folds at high conc);
  `"msa"` = `msa_transport` (inverse NaN at very low conc — bisection
  limitation → permeate-side gaps);
  `"data2"` = the DATA2 linear calibration
  `conc[mM] = 0.008813·cond[µS/cm] − 0.6949`
  (`DATA2_CONDUCTIVITY_CALIBRATION`).
- Round-trip exact for shedlovsky/data2; all conductivity channels (Tasks 3/5/8/9)
  route through these.

**Key finding.** The three backends **agree closely** on this campaign → the
probe-vs-ICP concentration gap is calibration/model, *not* the conversion
equation.

---

# PART B — new tasks (12–23), added 2026-07-22

*(Task 11+k = point k of the 2026-07-22 request.)*

---

## Task 12 — Pull in Alex's regression work + his notes PDF  🔶 IMPORT DONE, READING PENDING

**Goal.** Alex (Dowling, PI) built a regression analysis on a **different
process model** — an extended **Nernst–Planck / Donnan-equilibrium** membrane
description regressing the membrane **fixed charge χ** and a thermodynamic
**partition factor δ\*** — physics our DATA2-style Spiegler–Kedem (Lp, B, σ)
workflow does not have. Bring it into this branch and reconcile it with ours.

**Already done (2026-07-22).**
- Repo snapshot imported: `data3_regression_july26_AD/`
  (from https://github.com/dowlinglab/data3_regression, branch `main`,
  commit `5b8b1b7` of 2026-07-17, `.git` stripped; see `PROVENANCE.md` there).
- His write-up imported and renamed:
  `data3_regression_july26_AD/DATA3_regression_notes_AD_2026-07-07.pdf`
  ("Towards Multicomponent Diafiltration: Mathematical Modeling and Regression
  Notes", 62 pp). ⚠ The repo HEAD is 10 days **newer** than the PDF
  ("Added feedback from meeting with collaborators") — trust the repo where
  they disagree.

**Do (the real work).**
- [ ] Read the notes; map his model onto ours: χ, δ\* (Donnan) ↔ our B(c_int),
  σ; his CP-corrected `c_int` ↔ our film-model `cIn`.
- [ ] His statistics are directly reusable on OUR fits: autocorrelation-corrected
  confidence regions (i.i.d. CIs are ~2.5–2.7× too narrow on dense time series
  — this hits our WSSE contours too), effective sample size n_eff, profile
  likelihood, bootstrap, raw-measurement regression (§2.3–2.4 of the notes).
- [ ] Walk his §1.2 "open decisions" list; flag which ones touch our workflow.
- [ ] Short reconciliation note (md) in `data3_regression_july26_AD/` or
  `Architecture.md`: what we adopt, what stays separate.

**Done when.** The reconciliation note exists and at least the
autocorrelation/n_eff correction has a concrete plan (or implementation) for our
contour/CI story.

---

## Task 13 — Consolidate ALL plotting into `refactored_ucb_library.py`  🔶 PATTERN ESTABLISHED

**Goal.** Every plotting capability developed in standalone
`refactored_codes_v1/_*.py` scripts becomes a **library function** (+ a runfile
branch where it's a per-sheet campaign). **Especially the 5-panel contour plots
that are our norm for DATA3 single-salt sheets.** No orphan scripts.
*(This absorbs the scattered "new module or library helper" phrasing that was
in Tasks 2/3 — this card is now the single consolidation task.)*

**The worked example to copy (already done).** meas-vs-pred:
`plot_meas_vs_pred_conductivity` (figure fn) + `run_nf270_meas_vs_pred`
(per-sheet batch: load → fit → forward-sim → plot) + runfile branch letter `v`
+ the standalone `_make_meas_pred_nacl.py` **deleted** after the port.

**Do — port these, one at a time, same pattern (figure fn + `run_nf270_*`
dispatcher + runfile letter + delete/thin the script):**
- [ ] The **5-channel σ×B contour norm** (`_run_all_pairs_5channel.py`,
  `_run_all_11_sheets_5channel.py`, `_run_5channel_demo.py`,
  `_contour_batch_sigma_B.py`).
- [ ] The **7-channel** slice contours (`_run_Bsigma_channels.py`, Task 5).
- [ ] The **coefficient contours + A4 renderer** (`_run_Bsigma_coeffs.py`, Task 8)
  — keep the watchdog/resume machinery; it can live in the library too.
- [ ] Fixed-Lp contours (`_run_Bsigma_fixedLp.py`) and profile contours
  (`_run_profile_contours.py`, Task 9).
- [ ] Prediction generator (`_make_predictions.py`) → becomes the Task-3
  8-figure set.
- [ ] The stirred-cell mass-balance / interfacial-conc time plot
  (outputs in `bform_study/retentate_massbalance/`).
- [ ] Anything else `_make_*.py` / `_run_*.py` that renders a figure the
  campaign should be able to regenerate.

**Rules.** Solver-touching code guarded `workflow_family="DATA3"`; heavy
campaigns get a `run_nf270_*` batch entry point; scripts that remain (if any)
are thin CLI shims that import the library. `Architecture.md` gets one line per
new branch letter.

**Done when.** The runfile can regenerate every figure class end-to-end with no
standalone script logic; pytest green; deleted scripts leave no importers.

---

## Task 14 — Rejection plots: B(c_int) vs log Pe with σ-trajectories  ⬜ NOT STARTED
*(extends Task 6; do both variants with the same machinery)*

**Goal.** Per experiment: plot **B(c_interfacial) on Y vs log Pe on X**, with
**one trajectory per σ value**, so the σ-consistency of each run is visible.

**Careful — there are TWO Peclet numbers in play; label which one:**
- film/CP Peclet `Pe_film = Jw/k` (drives `c_in`), and
- membrane Peclet `Pe_m = Jw(1−σ)/B` (drives rejection).
For the σ-trajectories on a B-vs-Pe_m plot, note `B = Jw(1−σ)/Pe_m` — at a
given instantaneous `Jw`, each σ is a straight line of slope −1 in
log B–log Pe space; the run's per-vial (Pe, B(c_int)) points trace a trajectory
across them as `c_int` (and `Jw`) evolve.

**Do.**
- [ ] Per sheet: compute per-vial `Jw`, `c_int` (from the fitted forward sim),
  `B(c_int)` for the fitted form, both Peclets; plot B vs log Pe with the
  σ-family overlaid (σ grid 0.1…0.99); mark time/vial order (arrow or
  colorbar) so the trajectory reads as a path.
- [ ] Also render the **Task-6 classic**: `R(Pe; σ) = σ(1−e^−Pe)/(1−σe^−Pe)`
  σ-family vs log Pe with the observed per-vial rejection
  `R_obs = 1 − c_p/c_F` overlaid — best-supported σ read off the figure.
- [ ] **Suggested extra (pick up if useful):** the most physically direct
  variant is **B(c_int) vs c_int** (it IS the correlation we fit) with vial
  order marked, and `R_obs` vs `Jw` as the operating-curve companion. Offer
  all; recommend the (a) R-vs-logPe + (b) B-vs-logPe pair as the deliverable.

**Done when.** Both figures render per single-salt sheet (library functions —
Task 13 pattern), σ-family legend clear of data, Peclet definition stated on
the axis label.

---

## Task 15 — Derive the B–c_int correlation for each basis function  ⬜ NOT STARTED

**Goal.** A short derivation (as in the DATA2 paper,
`UnifiedFramework/DATA3/published_works/DATA2_main.pdf`) of what each basis we
fit — constant, linear, quadratic, cubic, saturating-exponential, Donnan —
**implies physically** for B(c_interfacial), and what each coefficient means.

**Why.** Right now the forms are curve-fitting choices. The DATA2 paper's route
(solution-diffusion permeability + concentration-dependent partitioning, Taylor
expansion about the operating window) turns each βₖ into a statement about
partition-coefficient derivatives — that's what makes the AIC/order-selection
result mean something.

**Do.**
- [ ] Follow the DATA2 derivation: start from `Js = B(c)·(c_in − c_H)` with
  `B = D_m·K(c)/L` (solution–diffusion), expand `K(c)` (Taylor for the
  polynomial family; Langmuir-type saturation for sat-exp; Donnan exclusion
  `K ∝ (c/χ)^{|z_co|/…}`-style for the charged form).
- [ ] Map: β₀ = dilute-limit permeability; β₁ ∝ ∂K/∂c; β₂ ∝ ∂²K/∂c²; sat-exp
  amplitude/scale = site saturation; Donnan exponent = valence/charge screening.
- [ ] Write it as a 2–3 page note (md or LaTeX; native equations), with the
  DATA2 equations cited by number; cross-link to the AIC order-selection
  result and to Alex's Donnan mechanistic model (Task 12) — his χ, δ\* should
  appear as the parameters of the Donnan-form K(c).

**Done when.** The note exists, each fitted coefficient has a physical reading,
and the Donnan-form derivation connects our B(c) to Alex's (χ, δ\*).

---

## Task 16 — The "beyond-Pe" evidence plot  ⬜ NOT STARTED

**Goal.** Pe-based (diffusion + convection only) reasoning assumes those are
the only dominant effects. Our campaign suggests otherwise. Produce **the one
plot that shows this with no room for doubt**, using one representative
experiment per salt (NaCl, CaCl₂, LaCl₃) × per regime (diluting,
concentrating) — 6 runs.

**Candidate figures (build (a); (b) is the kill shot; (c) optional).**
- [ ] **(a) Constant-coefficient SK failure, per run:** fit each run with
  constant (Lp, B, σ); plot observed per-vial rejection vs log Pe on top of the
  fitted `R(Pe; σ)` curve, and (inset or second row) the residuals vs `c_int`.
  If diffusion+convection sufficed, residuals are noise; instead they trend
  systematically with `c_int` — same membrane, same salt, coefficient must
  drift with concentration.
- [ ] **(b) Cross-salt collapse failure:** the 6 runs' `R_obs` vs log Pe on ONE
  axis. Pure diffusion–convection predicts each salt's points collapse onto a
  single σ-plateau curve (differences only via D in Pe). Show the valence
  ordering (NaCl > CaCl₂ > LaCl₃ rejection at comparable Pe) is far larger
  than the D-spread can produce — quantify: recompute Pe with each salt's D
  (Task 21 values) and show the curves STILL don't collapse.
- [ ] (c) `B_apparent` per vial (from `Js/(c_in − c_H)`) vs `c_int` across the
  6 runs — orders-of-magnitude drift at similar Pe = a "constant" that isn't.

**Done when.** The chosen figure(s) exist for the 6 representative runs with a
caption that states the argument in two sentences; the D-spread quantification
is in the caption or an SI table.

---

## Task 17 — Derive the replacement dimensionless number  ⬜ NOT STARTED

**Goal.** If Pe (convection/diffusion) is not the whole story, derive — don't
guess — the dimensionless group(s) that should join or replace it for our
charged-membrane system.

**Route.** Nondimensionalize the **extended Nernst–Planck + Donnan** problem
(the model Alex uses, Task 12; also
`box_multicomponent_diafiltration/Manuscript/JAO/Transport, MembranePhase -
Diffusion, Convection, Electromagration.docx`). The natural groups that fall
out:
- `Pe = Jw·δ/D` — convection vs diffusion (what we already use),
- **`ξ = χ / (z·c_int)`** — membrane fixed charge vs external ionic strength:
  the Donnan/Teorell–Meyer–Sievers exclusion parameter. Controls co-ion
  exclusion, hence rejection, hence the *apparent* B(c). Large ξ = strong
  exclusion; ξ → 0 = neutral membrane, Pe story recovered.
- `φ_D = z·F·Δψ_Donnan/(R·T)` — Donnan potential vs thermal voltage
  (equivalent information to ξ through the Donnan isotherm).

**Recommendation to verify in the derivation:** **ξ** is the right second axis —
it is built from Alex's regressed χ and our model's `c_int`, it is
per-vial-computable, and its c_int-dependence is exactly why B looks
concentration-dependent.

**Do.**
- [ ] 1–2 page derivation note: nondimensionalize, exhibit Pe and ξ (and φ_D),
  state the limits (ξ→0 neutral; Pe→∞ convective).
- [ ] Library helper computing per-vial ξ (χ from Alex's per-salt fits;
  `c_int` from our forward sim), DATA3-guarded.

**Done when.** The note + helper exist; ξ values per vial are available for
Tasks 18–19.

---

## Task 18 — What does (new number) × B(c_int) tell us?  ⬜ NOT STARTED

**Goal.** Re-do the Task-14 analysis with ξ in place of (or alongside) Pe and
answer: **is it useful?**

**The test is collapse.** Plot `B(c_int)` (and `R_obs`) vs ξ for the 6
representative runs (and then all single-salt sheets):
- If the salt-to-salt and regime-to-regime spread that Pe could **not** explain
  collapses onto one master curve in ξ — the charge physics is the missing
  driver, and ξ is the right correlating variable (answer: useful).
- If not, quantify what remains (partitioning non-ideality δ\*, dielectric
  exclusion) and say so honestly.

**Do.**
- [ ] B(c_int) vs ξ and R_obs vs ξ, per salt and pooled; same 6 representative
  runs as Task 16, then the full campaign.
- [ ] Quantify collapse: compare a pooled fit's WSSE/R² in ξ vs in Pe vs in
  c_int; report the ranking.

**Done when.** The collapse comparison (ξ vs Pe vs c_int) is plotted +
quantified, with a one-paragraph verdict on usefulness.

---

## Task 19 — Js and Jw vs Pe and vs the new number (2×2)  ⬜ NOT STARTED

**Goal.** Four panels per (representative) experiment set: `Js` vs Pe, `Jw` vs
Pe, `Js` vs ξ, `Jw` vs ξ. State up front what each can and cannot show.

**Expected information content (write this into the caption).**
- `Js` vs Pe — the diffusive→convective transition; curvature/slope reads the
  effective B; the most diagnostic of the four.
- `Jw` vs Pe — **near-tautological** (Pe ∝ Jw when δ/D is ~constant): expect a
  straight line; only *deviations* from linearity are informative (they mean k
  or D changed — CP regime shift). Plot it as the control, say so.
- `Js` vs ξ — charge regulation of solute flux: does Js drop as exclusion
  strengthens (ξ ↑ as c_int ↓)? This is the Donnan signature.
- `Jw` vs ξ — should be ~flat unless osmotic/charge effects feed back on water
  flux (σ·Δπ coupling); a trend here implicates the reflection term.

**Do.**
- [ ] Per-vial Js (= d(m_V·c_V)/dt / area, or from the model), Jw, Pe, ξ; 2×2
  figure per run + one pooled version colored by salt; library function
  (Task 13 pattern).

**Done when.** The 2×2 renders for the 6 representative runs + pooled; captions
state the expected/observed information content; a short verdict on which
panels carry signal.

---

## Task 20 — Electromigration: hammer or scalpel?  ⬜ NOT STARTED

**Goal.** You've been told explicit electromigration modeling is "a hammer for
a job that needs a scalpel." **Is that true?** Give the simple physics-based
answer with equations.

**The answer to write up (verify each step).** It is TRUE for single-salt
fitting, FALSE-but-tractable for multicomponent:
1. Extended Nernst–Planck flux:
   `J_i = −D_i∇c_i − z_i c_i (D_i F/RT) ∇φ + c_i v`.
2. Zero net current: `Σ z_i J_i = 0` → solve for `∇φ` and substitute back.
3. **Single binary electrolyte:** the migration term is eliminated *exactly*,
   leaving Fickian diffusion with the ambipolar (Nernst–Hartley) coefficient
   `D_± = (z_+ − z_−) D_+ D_− / (z_+ D_+ − z_− D_−)` — i.e. one scalar D per
   salt. Simulating the electric field explicitly adds **zero** physics for a
   single salt (the hammer); the scalpel is (a) using the correct ambipolar
   `D_±` (Task 21) and (b) Donnan partitioning at the interfaces (Alex's χ, δ\*).
4. **Multicomponent (2 cations + common anion):** migration coupling does NOT
   vanish — but it reduces to a concentration-dependent **coupled diffusion
   matrix** `D_ij(c₁,c₂)` (off-diagonals = one ion dragging another via the
   shared field). That is *exactly* what the BoE MATLAB does
   (`calculate_D` in `box_multicomponent_diafiltration/Analyses/MATLAB codes/
   BackoftheEnvelope_Calculations.m`, lines 98–122) — still no explicit field
   solve (still a scalpel, Task 22).

**Do.**
- [ ] Write the 1-page note with the three equations above + the D_± values for
  NaCl/CaCl₂/LaCl₃ (from Task 21); cross-ref the JAO transport docx variants
  (diffusion / +convection / +electromigration) in the Box archive and Alex's
  notes §1.5/§1.8.

**Done when.** The note exists and gives the advisor-ready two-sentence answer:
"for single salts, zero-current NP reduces exactly to Fick with D_± — explicit
migration adds nothing; for mixtures it becomes a D_ij matrix, which is the BoE
route, still no field solve."

---

## Task 21 — Are our diffusivity coefficients right?  ⬜ NOT STARTED

**Goal.** Audit every diffusion coefficient in our pipeline (film model `k`,
Pe, any D in the loader/analysis): which value, which convention (single-ion vs
salt/ambipolar), which source — and are they right for each cation/anion pair?

**Reference values to audit against.**
- Per-ion (BoE, `BackoftheEnvelope_Calculations.m` lines 21–24, m²/s):
  `D_Na = 1.33e−9`, `D_Ca = 0.79e−9`, `D_La = 6.26e−10`, `D_Cl = 2.03e−9`.
- Ambipolar salt values (Nernst–Hartley from those ions, 25 °C, infinite
  dilution): NaCl ≈ 1.61e−9, CaCl₂ ≈ 1.33e−9, LaCl₃ ≈ 1.29e−9 m²/s.
  **Watch for the classic error:** using the *cation* D for the *salt*
  (underestimates NaCl by ~20%, LaCl₃ by ~2×).
- Concentration dependence: Robinson & Stokes chapters are in the Box archive
  (`Conductivity Analyses/`) if we need D(c) beyond infinite dilution.

**Do.**
- [ ] Find every place a D (or k built from D) enters the library/loader; table
  of used-vs-literature values per salt, with convention labeled.
- [ ] Check the mass-transfer correlation too — the BoE stirred-cell form is
  `k = 0.23·v₀^0.57·D^0.67 / (ν^0.24·b^0.43)` (`.m` lines 124–130); compare
  against what our film model uses.
- [ ] Re-fit one sheet per salt with corrected D: report the effect on
  (Lp, B, σ) and on Pe (this feeds Tasks 14/16).

**Done when.** The audit table exists; any wrong value is fixed (DATA3-guarded)
or justified; the sensitivity of the fits to the D choice is quantified.

---

## Task 22 — Integrate the BoE coupled-diffusivity equation  🔶 ARCHIVE IMPORTED, PORT PENDING

**Goal.** Port the diffusivity-coefficient machinery derived in the Box
"Back of the Envelope" documents into `refactored_ucb_library.py` as opt-in
functions, so the multicomponent film model can use it.

**Already done (2026-07-22).** The Box folder is in the branch:
`box_multicomponent_diafiltration/` (see its `PROVENANCE.md`). The key sources:
- `Manuscript/Back of the envelope calculations for SI.docx` — the derivation
  (the boxed diffusivity-coefficient equation lives here).
- `Analyses/MATLAB codes/BackoftheEnvelope_Calculations.m` — reference
  implementation: `calculate_D(c1,c2)` builds the 2×2 coupled NP diffusion
  matrix `D_ij` (two cations + common anion, zero current), then the film
  problem is solved by eigen-decomposition
  (`Δc_a = T·diag(exp(Jw·δ/σ_eig))·T⁻¹·Δc_0`) → interfacial concentrations.
- `Analyses/Back of the Envelope Calculations_Coefficients Calculations.pptx`
  + `BoE Analysis.xlsx` — the worked numbers.

**Port carefully — known quirks in the `.m` (do NOT copy blindly, and do NOT
edit the `.m` itself):**
- Line 48: `Delta_La = params.DCa_s/k(params.DLa_s)` — numerator uses **DCa**
  (looks like a typo for DLa).
- `Delta_avg` averages Na/La/Cl only (sheet-specific: NaCl+LaCl₃ mixture);
  `calculate_D` hard-codes D1=Na, D2=La and `chi = 0` (placeholder).

**Do.**
- [ ] Library functions (DATA3-guarded, opt-in): `boe_coupled_D_matrix(c1, c2,
  ions=...)` (general z_i, D_i — not hard-coded), `boe_stirred_cell_k(D)`
  (the 0.23·v₀^0.57 correlation with rig constants b=0.0254 m, ν=1.003e−6 m²/s,
  350 RPM as defaults), and the eigen-decomposed film solve →
  `c_int` per ion.
- [ ] **Single-salt limit test:** with one salt, `D_ij` must collapse to the
  ambipolar D_± and the film solve to our current `exp(Jw/k)` CP model —
  pytest that.
- [ ] Cross-validate against `BoE Analysis.xlsx` numbers for one mixture vial.

**Done when.** The functions exist + tests pass (single-salt reduction, xlsx
cross-check); a mixture c_int can be computed end-to-end from the library.

---

## Task 23 — Sensitivity-analysis plots as library functions  ⬜ NOT STARTED

**Goal.** Standard sensitivity figures answering "are model + data sensitive to
our parameters?" — as **functions in `refactored_ucb_library.py`** (Task 13
pattern), not scripts.

**Do.**
- [ ] **Local, forward:** perturb each parameter (Lp, σ, each βₖ, S0/S, D, k)
  ±X% about the fit; forward-solve (square — no optimization); plot
  (a) trajectory ribbons (mass, cF, cV bands under the perturbation) and
  (b) a tornado chart of ΔWSSE per parameter (total + per response channel).
- [ ] **Local, derivative-based:** reuse `calc_FIM` — report per-parameter
  sensitivities `∂y/∂θ` scaled, eigen-spectrum of the FIM (ties to the
  contour/identifiability story; cite the flat-σ finding).
- [ ] Optional if cheap: a coarse Morris screening across all sheets.
- [ ] Entry points: `plot_parameter_sensitivity(...)` +
  `run_nf270_sensitivity(...)` + a runfile branch letter; precedent one-offs
  (`sigma_sensitivity-*.png`) get regenerated through the new path.

**Done when.** One command produces the ribbon + tornado + FIM-spectrum set for
a representative sheet (and batch for all sheets); figures land under
`paper_artifacts/nf270/sensitivity/`; pytest green.

---

## Verification checklist (run at the end of any working session)

- [ ] `pytest` DATA1 + DATA2 smoke/paper-comparison green (they stay untouched).
- [ ] Every new/changed figure rendered to PNG and **looked at** (legend clear,
  DATA2 colors, time correction visible in the mass panel).
- [ ] New library entry points reachable from the runfile; `Architecture.md`
  updated (branch letters, config knobs).
- [ ] Anything ported out of a standalone script: script deleted or thinned to
  a shim; no orphan importers.
- [ ] Heavy campaign outputs land under
  `UnifiedFramework/DATA3/results/paper_artifacts/nf270/` with a small index
  (JSON/md) per campaign.
- [ ] Derivation notes (Tasks 15/17/20) committed alongside the figures they
  explain.

## Assumptions baked in (change if wrong)

- Task 3 = **8 separate figure files** per experiment (+ the combined panel as
  an overview).
- Single-salt conversions via `variant_shedlovsky` inversion remain the default;
  MSA is the multi-salt hook.
- ξ (dimensionless fixed charge, Task 17) is the working candidate for "the new
  dimensionless number" until the derivation says otherwise.
- The BoE port (Task 22) generalizes the hard-coded Na/La pair to arbitrary
  (z_i, D_i) — the `.m` stays untouched.

---

## Sharing (export on demand)

This file is the single source of truth — no `.docx`/`.pdf` copies are kept in
the repo. To produce a shareable Word + PDF version when needed:

```bash
cd refactored_codes_v1 \
  && pandoc PROMPT_loader_timecorr_and_plots.md -f markdown+task_lists+pipe_tables+emoji \
       -t docx -o /tmp/PROMPT_loader_timecorr_and_plots.docx \
  && /opt/homebrew/bin/soffice --headless --convert-to pdf --outdir /tmp \
       /tmp/PROMPT_loader_timecorr_and_plots.docx
```

(`-f markdown`, not `gfm` — gfm ignores the dashboard table's column-width
hints. `$…$` math becomes native Word equations; images come along if their
relative paths resolve.) Don't commit the exports.
