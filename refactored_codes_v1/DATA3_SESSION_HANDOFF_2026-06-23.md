# DATA3 Session Handoff — 2026-06-23 (what's left to do)

Companion to `Architecture.md §20` (which records **what we did**). Everything here is
**DATA3 / NF270 only**. DATA1, DATA2, and the legacy MATLAB are untouched.

Environment for every command below:
```bash
cd refactored_codes_v1
export PATH="$HOME/miniforge3/envs/my-idaes-env/bin:$PATH"   # so pyomo finds ipopt
```

---

## 1. Regenerate DATA3 fits → figures → decks after the vial-close change  **(HEAVY — run on your machine)**

The vial-close permeate change (§20.1) moves the permeate residual from an interior
timestamp to vial close, so fitted **B/σ shift** and every DATA3 fit figure + deck must be
regenerated. Nothing has been regenerated yet.

- **Stage 1 — refit (hours, CPU-heavy):**
  ```bash
  python3 -u _run_bform_campaign.py all all 6
  python3 -u _run_Bsigma_fixedLp.py campaign
  ```
- **Stage 2 — re-render figures:**
  ```bash
  python3 _make_predictions.py all single
  python3 _replot_Bsigma_titled.py
  python3 _make_peclet.py
  python3 _make_donnan_reconciliation.py
  python3 _make_taylor_vs_donnan_lacl3.py
  python3 _make_taylor_vs_donnan.py MC3.07.11.24_SCaCl2 MC5.07.23.24_NaCl
  python3 _make_conc_flux_timeseries.py all
  ```
- **Stage 3 — rebuild deck:** `python3 _build_figure_pitch_deck.py`
  - Caption Lp/B/σ **auto-fill** from `result_single.json` (already wired in
    `_build_figure_pitch_deck.py` via `_fit_params`), so the slide numbers track the new fits.
- **Note:** Figs 9–11 aggregate across salts (Péclet regimes, Donnan master curve) → run the
  **full 11-sheet** campaign in Stage 1, not a subset.

## 2. Curated holdup boundaries — review / finalize

- The curated values in `nf270_holdup_boundaries.csv` come from **collection-line extrapolation
  to mass = 0** (`_curate_holdup_boundaries.py`), not experimenter annotation.
- Review `holdup_diagnostics/curated/<run_id>.png`; if you disagree, edit `holdup_end_s` (seconds,
  sheet time axis) — blank keeps the mass-rise auto-detect. The loader picks up edits on the next run.
- **Decision needed:** is the extrapolation curation final, or should lab-notebook boundaries replace it?

## 3. MATLAB-port analysis functions — run + open items

Run from the runfile: `python3 refactored_ucb_runfile.py` → ROOT=DATA3 → BRANCHES include **`x`**
(MATLAB ports; 2D+3D contours by default, DoE/σ-sensitivity heatmaps opt-in). Or call the library
directly: `lib.run_nf270_matlab_ports(...)`, or the individual `*_py` functions.

Open items:
- **Per-cell solver CPU cap (robustness).** Contour cells at pathological (B, Lp, σ) corners can
  grind thousands of ipopt restoration iterations before `calc_ind_objectives_py` returns NaN —
  no `solver_max_cpu_time` is threaded through `_nf270_contour_objectives_at_theta` → `solve_model`.
  This is inherited from the existing `calc_contour_2d_py`. Consider adding an optional CPU cap so
  3D sweeps don't stall on bad corners.
- **DoE / σ heatmaps are condition-sweep MBDoE tools** (hypothetical C_F0/C_D × ΔP), heavy. Validate
  outputs against expectations before trusting. The `C_D` Jw feasibility filter uses the data's
  initial retentate cF as the osmotic proxy (MATLAB used `max(cF)` from a forward sim) — refine if
  needed for diafiltration sheets.
- **Non-DATA3 B-bound divergence (documented, no DATA3 impact):** `calc_contour_2d_py` non-DATA3
  default `(1e-6, 50)` vs `calc_contour_3d_py` `(1e-6, 2)` — both mirror their MATLAB sources; DATA3
  uses `NF270_B_BOUNDS_PER_SALT` in both. Unify if desired.

## 4. DATA3-spec toggle — "we will get back to this"

Implemented (§20.5): `NF270_FORCE_DATA3_SPEC` (default `False`) + `run_with_data3_spec()`, honored by
`_nf270_conc_scales`, wired into the runfile NF270 dispatch. DATA1/DATA2 untouched. Possible follow-ups:
- Decide whether the toggle should also govern **B-bounds** and **σ interior bounds** (today it only
  forces the **weight scales**; B-bounds/σ-bounds come from `workflow_family` via `solve_model`).
- Add a **pytest invariant** asserting the DATA3 weights (cF 2% / cV 3%) are active during the DATA3
  workflow, so a future edit can't silently revert them.
- The shared `calc_FIM` default `workflow_family` was left as `'DATA1'` (changing it risks DATA1/DATA2
  regressions). If a cleaner contract is wanted, consider an explicit assert rather than a default flip.
- Decide whether any **other DATA3 entry points** (outside the runfile dispatch) should also wrap with
  `run_with_data3_spec`.

## 5. Other / pre-existing

- The legacy MATLAB loader `UnifiedFramework/Trials/Data_loader/load_data.m` still has **unresolved git
  merge-conflict markers** (`<<<<<<<` / `=======` / `>>>>>>>`). Not touched this session (per "don't
  change MATLAB"); flag for whoever owns that file.

---

## Files touched this session (for review/commit)

**Modified**
- `refactored_ucb_library.py` — vial-close permeate loader; curated-holdup override + helpers;
  `calc_contour_3d_py` / `doe_heatmap_py` / `heatmap_sigma_sensitivity_py`; `run_nf270_matlab_ports`;
  `_apply_sweep_condition` / `_initial_retentate_cF`; per-sheet mode propagation; `NF270_FORCE_DATA3_SPEC`
  + `run_with_data3_spec` + `_nf270_conc_scales` toggle.
- `refactored_ucb_runfile.py` — branch `x`; `run_with_data3_spec` wiring on the NF270 dispatch.
- `_build_figure_pitch_deck.py` — caption auto-fill from `result_single.json`.
- `_preflight_audit.py`, `_audit_loader_fidelity.py`, `Architecture.md` — doc/diagnostic updates.

**New**
- `nf270_holdup_boundaries.csv`, `_make_holdup_boundary_table.py`, `_curate_holdup_boundaries.py`,
  `holdup_diagnostics/` (+ `curated/`).
- `DATA3_SESSION_HANDOFF_2026-06-23.md` (this file).

---

## Update — 2026-06-24: collaborator-update deck + resume script

- **One-stop resume:** `RESUME_SESSION_2026-06-24.sh` — `bash RESUME_SESSION_2026-06-24.sh
  summary | env | verify | deck | figures | holdup`. Captures the whole session and reproduces it.
- **Deck built:** `DATA3_collaborator_update_2026-06-24.pptx` (8 slides) +
  `DATA3_collaborator_update_2026-06-24_backup.pptx` (13). Builder
  `_build_collab_update_deck.py`; assets in `_collab_update_build/`. White bg / dark text;
  **native editable OMML equations** (pandoc via `_eq2omml.py`); figures kept as crisp images
  (edit at source). Featured = figure_pitch slides 13–16 (Péclet regimes, Donnan master curve,
  Taylor-vs-Donnan LaCl₃/CaCl₂); 5-core by dropping the 3 explanatory slides.
- **Two envs:** `~/miniforge3/envs/my-idaes-env/bin/python` (fits/ports/loader; IPOPT) and
  `~/miniforge3/bin/python` (deck build; python-pptx + pandoc).
- **Also edited:** `_make_taylor_vs_donnan.py` / `_make_taylor_vs_donnan_lacl3.py` — softened the
  "truth" wording in the figure legend/suptitle.
- **Open item:** confirm in PowerPoint that the deck's equations are editable (a malformed
  double-`<a:pPr>` had broken this; fixed, OOXML validator now passes — but unverifiable outside
  PowerPoint here).
