# Prompt for Sonnet — collaborator-readability pass on the report (clarity, consistency, attribution)

You are polishing `docs/reports/main.tex` for a specific collaborator audience. This
is a **content-preserving** pass: **do not change any numbers, results, model
equations, or scientific conclusions** — the goal is that five specific readers can
follow it. Work in `/Users/adowling/DowlingLab/Membranes/data3_regression`, env
`data3-regression`.

**Run now** (the report is quiet — the Step 12 χ(c)-spline results have been
integrated). The coordinator (Opus) audited the report through these readers' eyes;
the findings below are line-anchored to save you rediscovery.

**Report state (updated):** the report is now **59 pages** and includes a new
**§6 "Toward a concentration-dependent membrane charge $\chi(c)$ (preliminary)"**
(`sec:chispline`, inserted after the cross-salt synthesis) --- so the Claude-Code
section is now **§7** and the appendices **§8--10**. The line anchors below at or
before ~line 2349 (Tasks B and F items in `sec:questions`, `sec:preproc`,
`sec:nacl`, and the earlier `17×` occurrences) are **unchanged**; anything after
~line 2349 shifted by roughly **+90 lines** (re-locate by section, not raw line
number, in the later material). The new `sec:chispline` introduced **no** new
"Phase/Step N" jargon and is already attributed (names Estrada's dataset) --- it
mainly needs the Task-E term-glossing (see below).

## The audience (write for all five)
1. **Full professor, membrane expert, modest applied-stats familiarity** (Bill Phillip).
2. **Two grad students** learning membrane science, no recent stats.
3. **A grad student** in detailed membrane modeling + numerical models in Pyomo.
4. **A grad student** in parameter estimation / optimal experiment design, moderate
   membrane familiarity.
5. **A postdoc** with a strong controls / dynamic-modeling background, still learning
   membrane science and applied statistics.
So: membrane physics can be assumed for (1); statistics can be assumed only for (4);
**(2),(3),(5) need plain-language on-ramps for whichever half is outside their field.**

## Task A — De-jargon the internal project scaffolding (highest priority)
`Phase 1/2/3` appears **~25 times** and `Step 8 / Step 9 / Step 14` **~6 times** in the
prose. **Collaborators do not know this internal work-plan**, so these read as
undefined references. Replace them with **descriptive** language:
- "Phase 1" → "the raw-signal reprocessing (moving concentration polarization
  upstream)"; "Phase 2" → "re-running the regressions CP-corrected"; "Phase 3" → "the
  autocorrelation and raw-measurement re-runs". Or define the three once, compactly,
  and then minimize repetition.
- "Step 8/9" → "the autocorrelation (AR(1)) analysis" / "the raw-measurement
  reformulation"; "Step 14" → "the multicomponent extension". Keep section
  cross-references (`\S\ref{...}`) — it's the bare "Phase/Step N" tokens that must go.

## Task B — Stale-note sweep (real inconsistencies: the report contradicts itself)
The concentration-polarization work is **complete** and its results are now primary,
but several notes still say it is pending / not done. Fix each to reflect completion:
- **Line ~272** (discussion point 11): "Step 8 (AR(1)) and Step 9 (raw-measurement)
  still use the pre-CP fit; re-running them CP-corrected is Phase 3" — **now done**;
  §\ref{sec:naclAR1}/§\ref{sec:altform} are CP-corrected.
- **Lines ~718–730** (`sec:preproc` follow-up item 5): "results-level resolution
  pending … re-running … is Phase 2, deliberately not done here … remains discussion
  point 11" — **now done**.
- **Line ~951** (`sec:nacl` opening `\note`): "Re-fitting to the CP-corrected values
  is an open decision" — **now done / resolved**; the section is CP-corrected primary.
- Sweep for any other "deliberately not done / pending / open decision" that CP
  completion has overtaken (check ~1664, ~1730, ~1760 — some are legitimately open
  future work; keep those, fix only the CP-superseded ones).

## Task C — Name the people (per PI request: name students, not just Dowling/Phillip)
`the authors'` appears **~35 times** (e.g. "the authors' CaCl2/LaCl3 processed
sheets", "authors' bulk `c_int`", "authors' inline calibration") and is ambiguous
(could mean the report's authors). These refer to the collaborators who **ran the
experiments and produced the processed sheets**. Attribution (from the raw workbook
experiment names + source-file authorship):
- **Laurianne Estrada** ran NaCl MC5 05.27.26 and CaCl2 MC5 05.27.26 (raw files
  `LE_...`); wrote the multi-salt MATLAB (`Donnan_Multiple_Salts.m`).
- **Jonathan Ouimet** ran NaCl MC2 05.07.24, CaCl2 MC3 07.11.24, and LaCl3 MC2
  05.21.24 (raw files `JAO_...`).
- The soft-sensor conductivity models (§\ref{sec:condmodels}) are **Lilonfe, Estrada,
  Singh, Ouimet**, Phillip, Dowling \cite{lilonfe2025softsensors}.
- "Bill's MATLAB" (already named ~11×) is Phillip's NaCl regression.
**Add a brief "data provenance / who ran what" note early** (e.g. end of
§\ref{sec:scope} or a footnote on the dataset table) naming Estrada and Ouimet per
dataset, then **replace the vague "the authors'" with the specific name(s)** where it
means their processed data/calibration (e.g. "Estrada's and Ouimet's processed
sheets", "the processed-sheet inline calibration"). Keep it collegial and factual.
Do **not** rename the report author line unless the PI directs.

## Task D — Cut redundancy in the calibration-source story
The vial-vs-inline / calibration-source thread is told at near-full length in **four
places**: discussion point 9 (~50 lines, §\ref{sec:questions}), `sec:preproc`
follow-up item 3 (~35-line `\note`), the Phase-1 validation paragraph, and
§\ref{sec:condmodels}. Consolidate: keep **one** full treatment (the `sec:preproc`
validation paragraph + §\ref{sec:condmodels} is the natural home), and **shorten
discussion point 9 and follow-up item 3 to a few sentences + a cross-reference**.
Same for the "our-CP vs authors'-CP ~13% / third sensitivity axis" framing — state it
once well, cross-reference elsewhere. Preserve every number; just stop repeating the
full argument.

## Task E — Reader on-ramps (add brief, plain-language scaffolding)
1. **Statistical-methods glosses (`sec:stats`)** for readers (1),(2),(5): add **one
   plain sentence per method** saying what it is *for* and how to read it — Wald
   ("simplest error bars; assume the model is locally straight; can be optimistic if
   curved"), nonlinear/$F$-test region, profile likelihood, bootstrap ("model-free
   resampling check"), multi-start. Keep the equations; add the intuition.
2. **A one-paragraph membrane/diafiltration physical primer** for readers (4),(5) and
   the students — early (start of §\ref{sec:model} or §\ref{sec:scope}): what
   diafiltration is, what is physically measured (flux, conductivity, ICP → ion
   concentrations feed vs permeate), and what $\chi$ (membrane fixed charge) and
   $\dstar$ (partition factor) *mean physically* and why they matter.
3. **Define key terms at first use, briefly**: "identifiability", "confounding",
   "autocorrelation", "effective sample size $n_{\mathrm{eff}}$", "concentration
   polarization" — one clause each. The controls postdoc will recognize the
   near-unit-root / non-stationarity discussion (§\ref{sec:altform}); a half-sentence
   connecting it to state-space/time-series framing they know is welcome.
4. **New §\ref{sec:chispline} ($\chi(c)$ spline) terms** for readers (1),(2),(5): the
   modeling/param-est students (3),(4) are fine, but gloss for the others —
   "varying-coefficient model" (a parameter becomes a smooth function of
   concentration), "B-spline / control points" (a flexible curve built from a few
   local pieces), "nested submodel" (constant-$\chi$ is the spline with no wiggle),
   "roughness penalty / P-spline", "pointwise band", and "adsorption isotherm" — one
   clause each at first use. The section's honest verdict (constant-$\chi$ is
   inadequate, but the $\chi(c)$ shape is not yet interpretable — the failed \ce{NaCl}
   null check) is important; keep that framing intact, just make the vocabulary
   accessible. It is a *preliminary/foreshadowing* section — ensure it still reads as
   a possible next step, not a settled result.

## Task F — Specific fixes
- **"$\gg17\times$ below NaCl"** appears ~4× (lines ~264, ~1896, ~1991, ~2227). The
  actual ratio is $56.74/2.95\approx19\times$. Either make it "$\approx19\times$" or,
  if the $\gg17$ lower-bound phrasing is intentional, ensure it's consistent and
  briefly justified in all four spots.
- Consider whether the abstract (very dense) can open with one plainer sentence
  framing the study for a mixed audience before the jargon; keep it short.

## Task G — Figure audit (math rendering + legend placement)
Rasterize and eyeball **every figure the report actually includes** (there are ~25;
pay special attention to the ones generated or regenerated since the last
figure-compliance pass: `cacl2_chispline_firstcut`, `cacl2_chispline_whitened`,
`conductivity_model_calibration`, the CP-regenerated `nacl_all_*`, and the AR(1)/raw
figures). Audit each for three things, and **fix the offenders**:

1. **Correct, consistent rendering of math symbols.** Symbols must render as real
   math, not spelled-out or garbled: Greek $\chi$, $\delta^{*}$; subscripted
   $c_{\mathrm{int}}$, $c_{\mathrm p}$; $\Delta c_{\mathrm m}$; super/subscripts; and
   the salts as \ce{CaCl2}/\ce{LaCl3} (or mathtext `CaCl$_2$`), **not** plain "chi",
   "delta", "CaCl2", "c_int" in titles/labels. Known offenders to start from: several
   figure \emph{titles} spell things out (e.g.\ `cacl2_chispline_firstcut`'s title
   "CaCl2 chi(c) spline…" and its axis uses proper $\chi$/$c_{int}$ — make the title
   match). Follow the established `analysis/plot_style.py` convention: bold axis
   titles with bold math via `\boldsymbol{}`. No raw LaTeX showing literally, no
   mojibake.
2. **Reasonable legend placement.** No legend should **overlap the plotted data** or
   **clip an axis label or title**. Known offenders to check: `cacl2_chispline_firstcut`
   and `cacl2_chispline_whitened` have a below-axis legend that appears to clip the
   "$c$ [mM]" x-label; the first-cut panel (a) legend sits "upper left" and may sit on
   the data. Follow the plotting convention (boxed legend, `frameon=True`; placed
   below the axes or outside; a single shared legend for multi-panel figures; move
   inside-the-axes legends that cover data). Use `tight_layout`/`bbox_inches` /
   `constrained_layout` so nothing is clipped.
3. **Font size $\ge$ the body text ($11$\,pt).** Every text element in a figure ---
   axis labels, tick labels, legend, title, annotations --- must render on the page at
   **at least the 11\,pt main-text size**, so a reader never squints at figure text
   smaller than the prose around it. The on-page size is the figure's native font size
   scaled by (display width $/$ figure's saved width):
   $\text{pt}_{\text{page}}=\text{pt}_{\text{native}}\times W_{\text{includegraphics}}/W_{\text{figsize}}$.
   Two failure modes to hunt for and fix:
   - **Scaled-down embedding.** A figure saved at `fig_width('full')` but
     `\includegraphics`-ed at, e.g., `0.6\linewidth` (or two separate
     `minipage{0.49\textwidth}` figures each at `width=\linewidth`) shrinks the font
     to $\sim$60\%/$\sim$50\% of native --- below body size. **Fix by matching the
     saved `figsize` width to the on-page display width** (`analysis/plot_style.py`'s
     `fig_width('full')`$=6.5$in for a `\linewidth` figure, `fig_width('half')`$\approx3.19$in
     for a genuinely half-width one) so the scale factor is $\approx1$.
   - **Explicit small-font overrides in the plotting code.** Some scripts hard-code
     `fontsize=7`/`fontsize=10` (e.g.\ the `chispline` figures use a 7-pt legend and
     10-pt titles) --- below body size even at scale 1. Remove these and inherit
     `plot_style.py`'s $\sim$11-pt defaults.
   When **side-by-side horizontal panels** force each panel (and its fonts) too narrow
   to reach 11\,pt at the page width, **stack the panels vertically instead** --- switch
   the script's `plt.subplots(1,2)` to `subplots(2,1)` and re-flow, or replace two
   side-by-side `minipage`-embedded figures with two full-width `\includegraphics`
   stacked vertically. A taller single-column figure that is legible beats a wide one
   that is not.

**Fixing requires editing the figure-generating scripts** (under `analysis/…`) and,
for font-size/stacking fixes, the figure's `\includegraphics`/`minipage` embedding in
`main.tex` — but all edits are **presentation-only**: titles, axis/tick/legend text,
legend location, on-page font size, `figsize`\,↔\,embed-width matching, panel
orientation (side-by-side\,$\to$\,stacked), and the `\includegraphics` width. Do
**not** change any plotted data, model curve, fitted value, axis \emph{range} that
would hide/rescale data, or numeric annotation. After a fix, regenerate both PNG and
PDF and **rasterize-verify** the figure renders correctly \emph{at its on-page size}
(check the font is $\ge$ body text). If a figure is already compliant, leave it
untouched.

## Guardrails
- **No changes to numbers, results tables, model equations, or conclusions.** This is
  clarity/consistency/attribution only. If you find a genuine numerical error, flag it
  in the report-back — do not silently "fix" a result.
- Edit `main.tex` for Tasks A–F. For **Task G only**, you may make
  **presentation-only** edits to the figure-generating scripts under `analysis/` and
  regenerate the affected figures (PNG+PDF) — no change to plotted data, curves,
  fitted values, or numbers. Do not touch `data/` or `refs.bib` (attribution uses
  names in prose, not new bib entries).
- Keep the `\todo`/`\note` flag mechanism; you're editing their content, not removing
  the flagging convention.
- Compile clean (`latexmk -pdf`: 0 undefined, 0 overfull; `latexmk -c`); report the
  page count.

## Report back
- Task A: how "Phase/Step N" was handled (defined-once vs replaced), with counts.
- Task B: each stale note fixed (line + old→new intent).
- Task C: the provenance note added, and how many "the authors'" were renamed.
- Task D: what was consolidated and where the single home now is.
- Task E: the primer + glosses added.
- Task F: the 17× resolution.
- Task G: which figures were audited, which had math-rendering, legend-placement, or
  **font-size (< 11 pt on page)** issues, and what presentation-only fix each got
  (including any side-by-side\,$\to$\,vertically-stacked panel changes and
  `figsize`/embed-width adjustments), with a rasterize-verify note; confirm no plotted
  data/curves/numbers changed.
- Clean-compile confirmation (0 undefined, 0 overfull, page count) and confirmation no
  numbers/results/equations changed.
