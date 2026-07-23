# Prompt for Sonnet — NaCl report refinement, round 2: show ALL the data + fix specific issues

You are doing a second refinement pass on `docs/reports/main.tex` in
`/Users/adowling/DowlingLab/Membranes/data3_regression` (11pt article; build with
`latexmk -pdf` in `docs/reports/`; env `data3-regression`). **Documentation/figure
task — do not change analysis results or the settled model mathematics;** verify
every number against the code / `analysis/preprocessing/README.md` / the processed
CSVs and report faithfully.

## Guiding principle (read first)
This is a **preliminary internal report**, not a paper. The PI and collaborators
want to **understand what is happening in every experiment**, so favor
**completeness and transparency over brevity** — it is fine, and desired, for the
report to be long and somewhat of a "data dump." The current report reads too much
like a summarized paper results section. **Show every NaCl dataset individually:
its preprocessing outputs and its fit/residuals — not just aggregate summaries.**

The tasks below come from specific PI comments; #5 is the big one.

## 1. Fix the dataset inventory table (`tab:nacldata`, currently "Table 6") — it is stale
Its status column still says several runs are "raw only --- needs preprocessing,"
but **all six have since been preprocessed and regressed.** Also the status
values are confusing (why is only one MC5 "analyzed"?). Fix by making the table's
purpose explicit and current. Recommended: reframe it as the **data provenance /
initial-state inventory** (which runs arrived already processed vs. raw), and add
a clear **current-status** indication that all nine are now processed and
regressed (cross-referencing where each is analyzed). A reader must be able to
tell, per run, what data we started with and what we have done with it.

## 2. Make it unmistakable that ALL nine datasets were regressed (`tab:naclall` + `tab:naclexp`, "Table 7"+)
Right now the results are split (3 originally-processed in one table, 6
reconstructed in another), so it reads as if only 3 were regressed. **Present all
nine in one place** — either a single comprehensive results table (with a column
marking "originally processed" vs. "reconstructed from raw"), or clearly linked
tables with an explicit statement of the count ("all nine NaCl runs were
regressed; three ... and six ..."). Remove the ambiguity.

## 3. Deep-dive the confidence-regions figure (`nacl_all_confidence_regions`, "Fig 4b")
The PI flagged an apparent issue with the MC5 confidence region. Diagnose and fix:
- **Confirm there is no numerical problem** with the well-identified runs' regions
  (MC2, MC5, MC5-repeat) — re-check the F-test/profile region computation for
  these; the Step-3 result had MC5 well-identified (|χ| CI ≈ 48.6–50.5), so the
  regions should be small and tight.
- The real problem is **visualization**: the axis is dominated by the poorly
  identified outlier runs (MC5-S at |χ|≈169, MC4 at δ*≈2), so the three
  well-identified runs collapse into an unreadable overlapping blob in the corner
  and their regions render as dots. **Add a zoomed panel/inset** for the
  well-identified runs showing the actual (small) F-test confidence regions, and
  **plot regions as regions, not just markers**. A two-panel figure (well- vs.
  poorly-identified, or full + zoom) is a good solution.
- **Fix the color collision**: MC2 and MC5 (07.23.24, S2) are both dark blue.
  Give every dataset a distinct, consistent color (reuse across all figures).

## 4. Fix Figure 5 (`nacl_all_fits`) — vertical compression
The 3×3 fit grid is too short: each panel is wide and vertically squished. Increase
the per-panel height (taller overall figure / larger aspect ratio) so the fits are
legible. Keep the shared "data / fit" legend below.

## 5. Add comprehensive per-dataset detail for ALL NaCl datasets (the main request)
The PI wants to see, for **every** NaCl run, both the preprocessing results and the
fit — even though this makes the report long. Add a **per-dataset appendix** (e.g.
"Appendix: per-dataset NaCl detail," one subsection or block per run) containing,
for each of the nine runs:
- **Preprocessing outputs** (for the 6 reconstructed runs especially): the
  per-experiment **conductivity calibration fit** (conductivity vs.
  dilution-corrected ICP, with the fitted line, slope/R²), the **Jw(t)** used, and
  the reconstructed interfacial concentrations $\cint(t),\cp(t)$. For the 3
  originally-processed runs, note they came pre-processed (no calibration to show).
- **Fit + residuals**: the corrected-model fit vs. $\Delta c_{\mathrm m}$ and a
  **residual-vs-$\cint$** panel, plus the run's $(|\chi|,\dstar)$ with CI and
  fit quality — so the reader sees which runs fit well and where the reconstructed
  ones break down (e.g. the reverse/diluting S/S2 runs diverge at high $\cint$;
  MC3 07.09.24 does not converge).
Generate these with new figures (use `analysis/plot_style.py`: `fig_width`,
below-axes/shared legends, PNG+PDF; add generating code under
`analysis/`). Keep the main §2.4 as the summary and reference the appendix for the
full detail. It is fine for this appendix to add many pages.

## Guardrails
- Do **not** alter analysis results/numbers, the settled model math
  (`sec:model`), notation, `\author`/preamble, or the CaCl2/LaCl3 stubs.
- Verify all numbers against code / `analysis/preprocessing/README.md` / CSVs.
- Use the shared `analysis/plot_style.py` for every figure; keep one consistent
  color per dataset across all figures.
- Compile clean: `latexmk -pdf` → **0 undefined citations, 0 overfull hboxes**;
  **rasterize (`pdftoppm`) and visually check** every new/changed figure and table
  for legibility and layout (this is why the PI flagged Figs 4b/5 — verify yours
  actually render well); `latexmk -c` after; do not commit `.bbl/.aux`; `main.pdf`
  is gitignored.
- May edit `main.tex` and add figures under `docs/reports/figures/` (+ code under
  `analysis/`). Do **not** modify `data/`, `prior_analysis/`, or `refs.bib`.

## Report back
- A change log (tables reframed, figure fixes, the new per-dataset appendix).
- What the Fig-4b deep dive found (numerical check result + the visualization fix).
- Any dataset whose per-dataset detail revealed something notable (e.g. the S/S2
  divergence, the non-converging run) — worth a sentence in the report.
- Confirmation: clean compile (0 undefined, 0 overfull), final page count, results/
  numbers unchanged, `data/`/`refs.bib` untouched.
