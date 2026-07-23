# Prompt for Sonnet — audit ALL report figures for publication-guideline compliance and fix them

You are auditing **every figure in the report** `docs/reports/main.tex`
(`/Users/adowling/DowlingLab/Membranes/data3_regression`, env `data3-regression`,
build with `latexmk -pdf` in `docs/reports/`) against the group's publication
guidelines, and **fixing any figure that is out of compliance.** This is a
figure/style task — do **not** change analysis results, numbers, or the model math.

Guidelines: <https://ndcbe.github.io/data-and-computing/notebooks/01/Publication-Quality-Figures.html>
plus the project conventions already in `analysis/plot_style.py`. Most fixes should
be made **centrally in `plot_style.py`** (so all figures and future CaCl2/LaCl3 work
inherit them), then figures regenerated; a few will need per-figure tweaks.

## Compliance checklist (apply to every figure)
1. **Font size matches the report body (~11pt).** This is a *combination* of (a) the
   figure's saved `figsize` and (b) the `\includegraphics` width in `main.tex`: the
   two must agree so the on-page scale ≈ 1 and the ~11pt rcParams fonts render at
   ~11pt. For **each** embedded figure, check its `\includegraphics[width=...]`
   against the `figsize` it is saved at (`plot_style.fig_width("full")`≈6.5in /
   `"half"`≈3.19in). If a figure is saved at, say, 6.5in but included at
   `0.6\linewidth`, its text is shrunk ~40% — fix by matching the figsize to the
   display width (or the width to the figsize). List any mismatches you find/fix.
2. **Bold labels/titles ⇒ bold math.** Axis labels and titles are bold, but
   matplotlib does **not** bold mathtext just because the text is bold — so a bold
   `[mM]` sits next to a non-bold `$|\chi|$`. Make the math bold too. Note
   `\mathbf{}` does **not** bold lowercase Greek (χ, δ); use **`\boldsymbol{}`**
   (e.g. `r"$\boldsymbol{|\chi|}$ [mM]"`, `r"$\boldsymbol{\delta^*}$ [--]"`,
   `r"$\Delta \boldsymbol{c_m}$ [mM]"`). Sweep every axis label, title, and any bold
   annotation across all figures. Verify by eye that symbols render bold.
3. **Box the legend.** Use `frameon=True` with a light, visible frame
   (`edgecolor` a light gray, `framealpha=1`) — a box around every legend. This
   **updates** the previous `frameon=False` convention; set it as the default in
   `plot_style.apply_style()` (e.g. `rcParams["legend.frameon"]=True`,
   `legend.edgecolor`, `legend.framealpha`) so it applies everywhere, and update the
   `plot_style` docstring accordingly.
4. **No text overlap.** Check title vs. y-axis label, tick labels vs. neighbors,
   legend vs. axes/other text, and (multi-panel) inter-panel labels. Fix with
   `constrained_layout`/`tight_layout`, `supylabel`/`label_outer`, adjusted
   `bbox_to_anchor`, or a taller/wider figure. (Recall Fig. `nacl_all_fits` was
   fixed for vertical compression and `tab:sources`-adjacent figures for
   title/ylabel collisions — hold that standard everywhere.)
5. **Legend vs. data / placement.** If a legend overlaps the data, move it: prefer
   **below the axes** (boxed) or **outside** the plot; and where 2+ subplots share
   the same legend entries, use **one shared figure-level legend** instead of
   repeating it per panel (as `nacl_multistart` / `nacl_all_fits` already do). Ensure
   no legend hides data points.
6. (Keep the existing standards: colorblind-distinct color per dataset reused across
   figures, inward ticks top+right, units in bracketed labels, save **PNG + PDF**
   with `bbox_inches='tight'`.)

## Scope — every figure embedded in the report
Audit all figures `\includegraphics`'d in `main.tex`. They are produced by these
scripts (regenerate via them after editing `plot_style.py`):
- `analysis/step2_nacl/regress_nacl.py` → `nacl_fit`, `nacl_residual_contour`
- `analysis/step3_nacl_uncertainty/scipy_uncertainty.py` → `nacl_multistart`,
  `nacl_confidence_scipy`, `nacl_profile_vs_wald`, `nacl_bootstrap_scipy`
- `analysis/step3_nacl_uncertainty/parmest_uncertainty.py` (~6 min) →
  `nacl_multistart_parmest`, `nacl_confidence_parmest`, `nacl_bootstrap_parmest`
- `analysis/nacl_all_datasets/regress_all_nacl.py` → `nacl_all_estimates`,
  `nacl_all_confidence_regions`, `nacl_all_fits`
- `analysis/nacl_all_datasets/per_dataset_appendix.py` → the 9 `nacl_appendix_*`
- `analysis/preprocessing/make_figures.py` → `nacl_preproc_calibration`,
  `nacl_preproc_reconstruction`
- `analysis/cacl2/regress_cacl2.py` → `cacl2_fits`, `cacl2_estimates`,
  `cacl2_residuals_vs_conc`, `cacl2_effective_chi_vs_conc`
- `analysis/cacl2/per_dataset_appendix.py` → the 3 `cacl2_appendix_*`
Confirm numeric outputs are unchanged after regenerating (this is a styling pass).
**Note:** the CaCl2 figures were generated *before* the boxed-legend / bold-math
conventions were adopted, so they especially need those fixes; treat them on equal
footing with the NaCl figures.

## Method (do this, don't skip the visual check)
1. Update `plot_style.py` centrally (bold-math approach, boxed-legend rcParams,
   any font/placement defaults) + its docstring.
2. Regenerate all figures; fix per-figure issues (labels needing `\boldsymbol`,
   overlaps, legend placement, figsize/`\includegraphics` width mismatches).
3. `latexmk -pdf`; then **rasterize every figure page** (`pdftoppm -png -r 150
   main.pdf ...`) and **visually inspect each one** against the checklist — the PI
   has repeatedly caught rendering issues that looked fine in code. Iterate until
   all pass. `latexmk -c` after.

## Guardrails
- Edit only `plot_style.py`, the figure-generating scripts (styling/labels/legend/
  figsize only — not the numerics), and `\includegraphics` widths in `main.tex`.
  Do **not** change analysis results, model math, prose, tables, `\author`/preamble,
  or `refs.bib`; do not modify `data/` or `prior_analysis/`.
- Compile must end clean: **0 undefined citations, 0 overfull hboxes**. `main.pdf`
  is gitignored; don't commit `.bbl/.aux`.

## Report back
- A **per-figure compliance table**: figure name → issues found (font-size/scale,
  bold-math, legend box, overlap, legend-vs-data) → fixed? A reader should see at a
  glance that every figure now complies.
- The `plot_style.py` changes made (boxed legend, bold-math handling, etc.).
- Confirmation: all figures regenerated (PNG+PDF) and **rasterize-verified**, report
  compiles clean (0 undefined, 0 overfull), page count, and analysis numbers
  unchanged.
