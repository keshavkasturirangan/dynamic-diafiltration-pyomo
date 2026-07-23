# Prompt for Sonnet — NaCl report refinement: make the report a self-contained, precise handoff document

You are refining the report `docs/reports/main.tex` in
`/Users/adowling/DowlingLab/Membranes/data3_regression` (11pt article; build with
`latexmk -pdf` in `docs/reports/`; env `data3-regression`). This is a
**documentation/audit task, not new analysis.**

## Goal
The report must become the **standalone handoff and discussion document** for the
collaborators (Laurianne, Bill, Molly): a reader should understand **all of the
NaCl analysis we have done — especially the preprocessing automation — without
opening any code, CSV, or the `analysis/preprocessing/README.md`.** Audit the
report extensively for content, fill gaps, and spell out the preprocessing
pipeline precisely. Do **not** change any numerical results or the settled model
mathematics — verify everything against the existing code/README/outputs and
report faithfully.

## Sources of truth (read these; do not modify them)
- `analysis/preprocessing/README.md` (validation tables, per-experiment
  calibrations, the (2)-sheet finding, CP deviations, all flags),
- `analysis/preprocessing/{preprocess_nacl,validate,process_raw_runs}.py`,
- `analysis/preprocessing/processed/*.csv`,
- `analysis/step2_nacl/regress_nacl.py`, `analysis/step3_nacl_uncertainty/*.py`,
  `analysis/nacl_all_datasets/regress_all_nacl.py`.
Cross-check every number you add against these; flag (do not invent) anything you
cannot source.

## Part 1 — spell out the preprocessing automation precisely (highest priority)
The methods + flags are in `\S`~`sec:preproc`, but the **validation results and the
step-by-step automation** currently live only in the code/README. Bring them into
the report so it stands alone:
1. **Raw→processed data-flow.** Add a short table (or itemized map) listing each
   processed column (B `c_int`, D `I`, F `R`, H `c_p`, J `Jw`, K `B`, N/O
   conductivities) → how it is produced (direct copy / calibration / CP / derived
   formula) → which raw input(s) it uses. A collaborator should see the whole
   pipeline at a glance.
2. **Validation results as report tables/prose.** Summarize, in the report:
   - the `(2)`-sheet check — exact (machine precision) on the 22 live-formula rows,
     and the finding that the other ~97% are pasted values deviating up to ~38%;
   - the CP-corrected-sheet deviations (main MC5, MC2) and what drives them;
   - the per-experiment calibration table (6 runs: fitted slope, R²) with the
     93–98 range and agreement with the embedded 76.685;
   - the Jw smoothing accuracy (51-row rolling slope, ~3.8% median error).
3. **Illustrative figure(s)** (use `analysis/plot_style.py`: `fig_width`,
   below-axes/shared legends, PNG+PDF): e.g. a per-experiment **calibration fit**
   (conductivity vs concentration with the fitted line) and/or a
   **reconstructed-vs-reference** panel (pipeline output vs the `(2)`-sheet on its
   formula rows). Keep it honest — show the scatter where it exists. Add the
   generating code to `analysis/preprocessing/` if you make new figures.
4. Ensure every **flag** in `sec:preproc` is stated as a crisp, self-contained
   open question (a reader shouldn't need context we haven't written down).

## Part 2 — audit the rest of the NaCl analysis for completeness
Section by section (`\S`~`sec:reproMATLAB`, `sec:naclUQ`, `sec:naclAll`), confirm a
reader can follow **what we did, why, and what we found** without the code:
- Step 2 (MATLAB reproduction + factor-of-2): is the setup, the exact-reparameterization
  finding, and the recommendation fully stated? 
- Step 3 (uncertainty): are the multi-start convergence diagnosis, nonlinear vs.
  Wald regions, profile vs. Wald CIs, bootstrap, and the ParmEst cross-validation
  each explained in prose (not only in figure captions)? Several Step-3 companion
  figures are referenced but not embedded — either embed the important ones or make
  sure the text conveys their content.
- Multi-dataset + expanded set (§2.4): is the MC2-vs-MC5 story and the
  reconstruction-dominated expanded-set caveat clear?
Fill genuine gaps with concise, accurate prose. Do not pad.

## Part 3 — facilitate the conversation
Make the open decisions easy to find and discuss. Ensure the modeling questions
(`sec:questions`), the `sec:preproc` follow-up items, and the recommendations
(`sec:naclDisc`) are consistent and cross-referenced. If it helps, add a short
consolidated **"Open decisions for the team"** list (or ensure `sec:questions`
covers the preprocessing items too) — the report should read as an agenda for the
collaborative meeting. Keep the existing `\todo`/`\note` color flags.

## Consistency pass
Verify numbers in prose match the tables, figures, and code outputs (e.g.
$|\chi|\approx50$, the expanded-run values, calibration slopes, deviations). Fix
mismatches toward the code/README values.

## Guardrails
- **Do not alter** analysis results/numbers, the settled model math (`sec:model`
  transport/Donnan/statistics), notation/conventions, or the `\author`/preamble.
- Leave the `\ce{CaCl2}`/`\ce{LaCl3}` sections as the clearly-marked planned stubs.
- Compile clean: `latexmk -pdf` → **0 undefined citations, 0 overfull hboxes**;
  rasterize (`pdftoppm`) and eyeball any new/changed figures and tables for
  layout; run `latexmk -c` after (do not commit `.bbl/.aux`; `main.pdf` is
  gitignored).
- May edit `main.tex` and add figures under `docs/reports/figures/` (+ their
  generating code under `analysis/preprocessing/`). Do **not** modify `data/`,
  `prior_analysis/`, or `refs.bib`.

## Report back
- A **section-by-section change log** (what you added/clarified and why), so Opus
  can review efficiently.
- The new tables/figures added and where.
- Any content you could not source or that needs a human decision (flagged).
- Confirmation of a clean compile (0 undefined, 0 overfull), final page count, and
  that results/numbers and `data/`/`refs.bib` are unchanged.
