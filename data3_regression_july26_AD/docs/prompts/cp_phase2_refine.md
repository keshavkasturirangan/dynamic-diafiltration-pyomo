# Prompt for Sonnet — CP Phase 2 refinement: fix the MC2 CP mis-classification and the mixed-convention all-NaCl table

A coordinator audit of CP Phase 2 (commit `8936050`) found **one correctness bug and
one convention cleanup**, both confined to the **all-NaCl multi-dataset analysis**
(`analysis/nacl_all_datasets/regress_all_nacl.py` + its report table/figures/synthesis
text). The primary NaCl fit (step2), \ce{CaCl2}, and \ce{LaCl3} are correct and must
**not** change. Work in `/Users/adowling/DowlingLab/Membranes/data3_regression`, env
`data3-regression`.

## Bug 1 (correctness): MC2 05.07.24 is non-CP but was marked `already_cp: True`
In `regress_all_nacl.py::DATASETS`, `NF270_MC2 05.07.24_NaCl` has
`"already_cp": True`, so the Phase-2 CP correction was **skipped** for it and its
"CP-corrected" row is actually the non-CP value. **It is in fact non-CP.**
Verification (coordinator; you should re-confirm): the $\cint/\kappa_{\mathrm{ret}}$
coefficient of variation is **0.096** for MC2 — in the non-CP cluster (MC5(2) = 0.088)
and far from the genuinely-CP MC5-main (0.255). Applying the Phase-2
authors'-bulk CP mechanism to MC2 shifts its $\cint$ by a **median 18.8%**
(14–49%).
- **Fix:** set MC2's `already_cp: False` so `cp_authors_bulk.apply_cp` is applied to
  its bulk $\cint$ (exactly as for MC5(2)). Re-run; MC2's CP $|\chi|$ will rise from
  18.46\,mM (as NaCl's did, $\sim$+13\% or more given the larger shift).
- Re-confirm the CP-status of every Excel NaCl dataset with the CV diagnostic and
  make `already_cp` match reality: MC2 → non-CP (False); MC5-main → CP (True, correct);
  MC5(2) → non-CP (False, correct). (The 6 raw-only CSV runs correctly use their
  Phase-1 `*_CP.csv` twins.)

## Bug 2 (convention): the all-NaCl CP table mixes our-CP and authors'-CP
After fixing MC2, the datasets carrying a CP correction are MC2, MC5(2), and the 6
raw-only runs --- all using **our** thin-film CP. The lone exception is **MC5-main**,
whose column is the **authors' own** CP. Presenting MC5-main (authors'-CP,
$|\chi|\approx49.6$) beside MC5(2) (our-CP, $56.07$) currently reads as an
unexplained "MC5 $\ne$ MC5(2)".
- **Fix / reframe:** make the primary CP set **uniformly our-CP** --- use
  **MC5(2)+our-CP as the MC5 05.27.26 entry** in the primary CP table/figures, and
  present **MC5-main (authors'-CP) as an explicitly-labeled comparison** (not a second
  primary row). Do **not** apply our CP to MC5-main (keep `already_cp: True` there ---
  double-applying would be wrong); just relabel its role.
- Frame the $\sim$13\% difference honestly: it is the **our-CP vs.\ authors'-CP
  implementation discrepancy** --- a *third* CP-related sensitivity axis alongside
  CP-vs-non-CP (point 11) and calibration-source (point 9, \S\ref{sec:condmodels}).
  It arises because our Leveque-type thin-film $k$ differs from whatever CP the
  authors applied (Phase 1 already found our-CP reproduces the authors' CP MC5 sheet
  only to $\sim$12\% median). Add one or two sentences noting this as a minor open
  item; do not over-dramatize it.

## Report fixes (`docs/reports/main.tex`, all-NaCl area only)
- The `tab:naclall` MC2 row: remove the "(already CP)" annotation; MC2 is CP-corrected
  by our mechanism like the others. Update its $|\chi|$/$\dstar$/SSE to the new
  CP values.
- The synthesis line stating "MC2 $|\chi|=18.46$\,mM, unaffected by CP" is now wrong ---
  correct it to MC2's new CP value and note it *was* CP-corrected.
- Update any other place that cites the MC2 CP number or the MC5/MC5(2) coincidence
  framing (the earlier text said applying CP to MC5(2) does not reproduce MC5's
  already-CP value --- reframe per Bug 2 as the our-CP-vs-authors'-CP discrepancy).
- Refresh the all-NaCl table/figures (`nacl_all_*`) and the MC5-repeat/MC2-vs-MC5
  discussion accordingly.

## Guardrails
- **Do not change** the step2 primary NaCl fit, \ce{CaCl2}, \ce{LaCl3}, the
  conductivity-models section, `refs.bib`, `data/`, or `prior_analysis/`. Scope is
  `regress_all_nacl.py` + the all-NaCl table/figures/synthesis text in `main.tex`.
- Keep `APPLY_CP=False` byte-for-byte reproducing the committed non-CP numbers
  (regression test) --- verify after the change.
- Recompile clean (`latexmk -pdf`: 0 undefined, 0 overfull; `latexmk -c`); report page
  count. Rasterize + eyeball the regenerated all-NaCl figures.

## Report back
- The corrected `already_cp` classification (with the CV diagnostic per NaCl dataset),
  and MC2's non-CP vs CP $|\chi|$ (old 18.46 vs new).
- The reframed MC5 primary (MC5(2)+our-CP) vs authors'-CP comparison, and the stated
  our-CP-vs-authors'-CP $\sim$13\% discrepancy.
- Confirmation the step2/CaCl2/LaCl3 numbers are untouched and `APPLY_CP=False` still
  reproduces committed values.
- Clean-compile confirmation (0 undefined, 0 overfull, page count).
