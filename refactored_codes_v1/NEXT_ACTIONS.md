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
- **Contour restyle COMPLETE (2026-06-12):** band-masked (`_run_band_contour.py`) and
  β (`_run_beta_contour_fast.py`) renderers now route through `lib._plot_heatmap_frame`
  + save grid CSVs; both reran. β driver now derives the β₁ sweep range from
  feasibility (B = β₀+β₁·cIn within bounds over the cF span) so the σ×β₁ panel
  isn't blank.
- **Deck DONE (2026-06-12):** 6 result slides added to `_build_followup_deck.py`
  (per-band σ-recovery table + overlay, band-masked contours, per-vial B-vs-cF
  trend, campaign-wide σ-recovery, B(I) identity/cross-salt, β contours); rebuilt
  (42 slides) and PDF-verified via LibreOffice. `.pptx` NOT committed (item 1).
- **Campaign band run DONE (2026-06-12, `_run_band_campaign.py`):** all 11 single-salt
  sheets. HONEST finding — σ-recovery (high-cF band interior, full railed) is
  **sheet-specific: 3/11**; the robust generalizable signal is **B rising with cF
  on 9/11**. Tempers the clean 3-sheet value-test; deck slide 7 + Architecture §18.4
  + memory updated. (This was old to-run item 4.)

## To run when you return (ordered)

1. **Large-binary tracking decision (BLOCKING for a clean repo).** NOT committed:
   - decks `refactored_codes_v1/*.pptx,*.pdf` (~204 MB)
   - artifacts `UnifiedFramework/DATA3/results/paper_artifacts/nf270/{animations3d 65M,
     contour3d 24M, meeting_prep_2026-06-10 5.9M, matlab_ports* , band_value_test 1.1M, …}`
   - `pytest_refactored_codes_v1/` (35 MB test cache), the 19 modified `data1/notebook_figures/*.png`
   Decide: git-lfs, external store, or `.gitignore` (then `git rm --cached`). The
   small deliverables (`band_value_test/`, `b_ionic_strength/`) you may want tracked.

2. **Task 2 deeper (optional):** cross-salt was negative — B is more salt-specific
   than k_I explains. Next: add a salt/valence exclusion term, or a JOINT multi-salt
   fit sharing β in I-space with per-salt Lp/σ, and re-test transferability.

3. **Per-vial joint model (optional escalation):** the per-vial trend used rolling
   windows. The heavier `m.Lp[n]/m.sigma[n]` joint model (mirroring `B_form='pervial'`)
   is the next step IF the trend warrants — gate carefully to preserve DATA1/DATA2.

4. **3-band test (optional):** the campaign used 2 bands. Try 3 bands on ≥9-vial
   sheets (`solve_model_per_concentration_band(n_bands=3, max_bands=3)`) to resolve
   finer B–cF structure.

5. **Deck polish (optional):** the 6 new result-slides (deck slides 4–9) are
   PDF-verified for layout; review wording/emphasis with the collaborators and
   re-order vs. the existing narrative if desired.

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
