# Action items — tomorrow evening (2026-06-19)

## 1. Re-fit the single-salt constant-B fits at pure-2% weighting (full end-to-end consistency)

**Why.** The contour / 3D-grid / animation pipeline centers on
`bform_study/<run_id>/result_single.json`, but those fits were computed at an EARLIER cF
weighting (legacy 0.3% / mixed), NOT the canonical Lilonfe 2% no-floor. Lp* is mass-driven
and robust, so the Lp-identifiability story holds; but B* and σ* (and the whole fit) should
be redone at 2% so the ANCHOR and the swept contour use identical weights. (Recall the
σ-flip finding: σ's minimizer moves with the cF weight.)

**How.**
- Pin weighting: `NF270_CF_RESIDUAL_SCALE_FRACTION = 0.02`, `NF270_CF_RESIDUAL_FLOOR_MM = None`,
  `workflow_family="DATA3"`.
- Back up the existing `result_single.json` (→ `result_single.0p3pct.json`) for all 11 sheets.
- Re-fit constant-B for all 11 single-salt sheets (`_run_bform_campaign.py` single, or
  `solve_model(sim_opt=False, B_form='single', …)` with multistart). Write new `result_single.json`.
- Regenerate downstream on the 2%-anchored fits: 3D grids (`_run_3d_contour_concentrating.py`),
  B×σ@Lp* (`_run_Bsigma_fixedLp.py single`), animations (`_make_3d_contour_animations.py`).
- Compare new vs old Lp*/B*/σ* (expect Lp* stable, σ* may shift). Note any sheet whose σ
  flips interior↔rail under the honest weight.

**Scope.** DATA3 single-salt only. DATA1/DATA2 must stay byte-identical (regression guards).

## 2. (after re-fit) Refresh the figure-pitch deck figures with the 2%-anchored versions.
