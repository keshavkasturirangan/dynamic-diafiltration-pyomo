# DATA3 NF270 follow-up — resume checklist (paused 2026-06-10)

State at pause. Tasks 1–3 of the §18 follow-up are **implemented, validated, and
committed** (commits `431faf3` code+docs, `0b29ad4` scripts). DATA1/DATA2 stay
byte-identical: `python3 refactored_codes_v1/tests/regression_data1_data2_guards.py`
(6/6 PASS) + `pytest refactored_codes_v1/tests/test_data1_regression.py
refactored_codes_v1/tests/test_data2_regression.py` (2 passed, 2 skipped).

## Done
- **Task 1 — β contours:** validator + `sweep_pairs` β-slice gating wired.
- **Task 2 — B=f(I):** `NF270_B_USE_IONIC_STRENGTH`, `k_I·cIn` in `B_form_rule2/3`.
  Rescale identity validated to machine precision; cross-salt unification by ionic
  strength **refuted** (honest negative — `b_ionic_strength/b_ionic_validation.json`).
- **Task 3 — per-band θ:** `band_vials` mask (solve_model/B_fix/calc_FIM/contour),
  `solve_model_per_concentration_band` + per-band FIM. Value shown on GOOD/OKAY/POOR:
  high-cF bands recover interior σ; B rises with cF in MC3 SNaCl (corr +0.93);
  POOR LaCl3 bands rail to opposite σ bounds (flat surface). FIMs non-singular.
- Architecture §18.3/§18.4/§18.6 updated.
- **Contour restyle:** band overlays + per-vial trends re-rendered in the canonical
  figure_s5 line style (`_rerender_band_style.py`).

## To run when you return (ordered)

1. **Finish the figure_s5 restyle** (band overlays + per-vial already done):
   - The band-MASKED contour (`_run_band_contour.py`) and β contour
     (`_run_beta_contour_fast.py`) renderers still use `contourf` fill. Route them
     through `lib._plot_heatmap_frame(df, x, y, z, ax=...)` like
     `_rerender_band_style.py` does, then re-run:
     ```
     python3 refactored_codes_v1/_run_band_contour.py MC3.07.22.24_SNaCl 16
     python3 refactored_codes_v1/_run_beta_contour_fast.py MC3.07.22.24_SNaCl 18
     ```
   - (The fast β-contour run was started then killed to free CPU — never finished;
     it needs the restyle + a clean run. ~5 min.)

2. **Deck update** — add slides to `_build_followup_deck.py` for: B(I) rescale
   identity + cross-salt negative; per-band σ-recovery table (GOOD/OKAY/POOR);
   band-masked contour (flat-σ vs curved-σ basins); per-vial B-vs-cF trend.
   Rebuild → `DATA3_followup_analysis_2026-06-11.pptx`.

3. **Large-binary tracking decision (BLOCKING for a clean repo).** NOT committed:
   - decks `refactored_codes_v1/*.pptx,*.pdf` (~204 MB)
   - artifacts `UnifiedFramework/DATA3/results/paper_artifacts/nf270/{animations3d 65M,
     contour3d 24M, meeting_prep_2026-06-10 5.9M, matlab_ports* , band_value_test 1.1M, …}`
   - `pytest_refactored_codes_v1/` (35 MB test cache), the 19 modified `data1/notebook_figures/*.png`
   Decide: git-lfs, external store, or `.gitignore` (then `git rm --cached`). The
   small deliverables (`band_value_test/`, `b_ionic_strength/`) you may want tracked.

4. **Task 2 deeper (optional):** cross-salt was negative — B is more salt-specific
   than k_I explains. Next: add a salt/valence exclusion term, or a JOINT multi-salt
   fit sharing β in I-space with per-salt Lp/σ, and re-test transferability.

5. **Per-vial joint model (optional escalation):** the per-vial trend used rolling
   windows. The heavier `m.Lp[n]/m.sigma[n]` joint model (mirroring `B_form='pervial'`)
   is the next step IF the trend warrants — gate carefully to preserve DATA1/DATA2.

6. **Broaden band analysis:** extend to all ≥6-vial sheets; try 3 bands on ≥9-vial
   sheets (`solve_model_per_concentration_band(n_bands=3, max_bands=3)`).

## Re-run cheatsheet
```
python3 refactored_codes_v1/_run_band_value_test.py     # per-band fits + overlays
python3 refactored_codes_v1/_run_band_wrapper.py        # formal wrapper + per-band FIM
python3 refactored_codes_v1/_run_per_vial_trend.py      # per-vial trend + vial↔contour corr
python3 refactored_codes_v1/_run_b_ionic_validation.py  # B(I) identity + cross-salt clustering
python3 refactored_codes_v1/_rerender_band_style.py     # restyle overlays (figure_s5)
```
Gold/value-test sheets: GOOD `MC3.07.22.24_SNaCl`, OKAY `MC2.05.07.24_NaCl`,
POOR `MC2.05.21.24_LaCl3`. Set `lib.NF270_CF_RESIDUAL_FLOOR_MM=1.0` (campaign cfg).
