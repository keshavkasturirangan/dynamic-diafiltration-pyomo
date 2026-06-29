# Prompt: Reproducibility Audit — DATA1 + DATA2 (Refactored Codebase)

> Authored 2026-06-09 from the operator's explicit framing. Locked here so
> the scope can be audited later. Companion artifact:
> [`REPRODUCIBILITY_AUDIT.md`](./REPRODUCIBILITY_AUDIT.md) (the actual report).

---

## 1. Why this audit exists

The operator's stated goal:

> "Use only one code to produce results for new experiments or any old
> experiment, published or unpublished. The point of the pytest is to
> ensure that irrespective of the changes we make to the codebase moving
> forward, nothing changes our capability to reproduce our prior results.
> This fear comes from wanting to build a code-based system that is as
> light as possible, where we activate different features of the process
> model or analysis to generate different results."

In one line: **prove that `refactored_codes_v1/` alone can reproduce every
published DATA1 + DATA2 figure and table**, so that future feature flags +
process-model toggles never silently break that capability.

## 2. Authoritative source-of-truth = the four published PDFs

The audit's "targets" are NOT the `target_id` strings in
`UnifiedFramework/DATA3/docs/validation/target_notebook_source_map.csv`.
They are every Figure (with all sub-panels) and every Table in:

- `UnifiedFramework/DATA3/published_works/DATA1_main.pdf`  (Ouimet et al. 2022, J. Membr. Sci.)
- `UnifiedFramework/DATA3/published_works/DATA1_SI.pdf`     (Supporting Information)
- `UnifiedFramework/DATA3/published_works/DATA2_main.pdf`   (Liu et al., NF270 / MBDoE)
- `UnifiedFramework/DATA3/published_works/DATA2_SI.pdf`     (Supporting Information)

The locked CSVs (`target_notebook_source_map.csv`,
`simulation_validation_manifest.csv`) are USEFUL CONTEXT for the audit but
are not authoritative — if they disagree with the PDFs, the PDFs win.

## 3. Authoritative target codebase

The audit checks reproducibility **only against**:

- `refactored_codes_v1/refactored_ucb_runfile.py`  (orchestrator — runfile, tree-shaped dispatcher)
- `refactored_codes_v1/refactored_ucb_library.py`  (all science + plotting lives here)
- `refactored_codes_v1/conductivity_paper.py`      (read-only — DO NOT MODIFY)

No other Python source counts toward "reproducible from current code." The
goal is a single self-contained codebase.

## 4. Legacy material — context only, read-only

These predate the refactor and are referenced ONLY as context for what the
target should produce. The audit does not consume them at test time, does
not modify them, and any port to Python must end up in
`refactored_ucb_library.py`, not in these files.

### MATLAB scripts in `legacy/data1_matlab/functions/`

Used to produce DATA1 (and some DATA2) figures and `.mat` artifacts:

- `run_data_analysis.m` — fits Lp, B, σ; simulates; optional contour exploration
- `doe_heatmap_diafiltration.m` — A/D/E/modified-E heatmaps for diafiltration design
- `doe_heatmap_filtration.m` — same metrics for filtration design
- `run_sigma_sensitivity.m` — σ-sensitivity for filtration + diafiltration
- `heatmap_sigma_sensitivity_diafiltration.m` — σ heatmaps (mass / retentate / permeate)
- `heatmap_sigma_sensitivity_filtration.m` — σ heatmaps for filtration

### Python references

- `legacy/data1_matlab/functions/diafiltration_plots.py` — pre-refactor plotting functions
- `legacy/data1_matlab/DiafiltrationPaperPlots.ipynb` — DATA1 paper-figure reproductions
- `legacy/data1_matlab/PSE2021Plots.ipynb` — DATA-MBDoE figure reproductions
- `utility.py` (repo root) — pre-refactor Python utilities (parameter estimation, FIM, visualization)
- `DATA1_model_demo.ipynb`, `DATA2_model_demo.ipynb`, `DATA2_visualization.ipynb` (repo root) — pre-refactor notebooks that produced the published numbers

**Rule: catalog these for context; do not modify.**

## 5. Per-target audit method

For every Figure + Table in the four PDFs:

1. **Trace** which `refactored_ucb_library.py` function would produce it,
   via `materialize_all` (or `run_pipeline`, `run_data1_direct_contour_branch`,
   etc.).
2. **Check** whether that function exists in the current library.
3. **Resolve** to one of four verdicts:

   | Verdict | Meaning |
   |---|---|
   | `PRESENT` | `refactored_codes_v1` can produce this from scratch with NO `.mat` dependency. |
   | `PRESENT_WITH_MAT_DEPENDENCY` | A Python function exists but it consumes a pre-computed `.mat` from MATLAB (e.g. `legacy/data1_matlab/data/*sigma_sensitivity*.mat`). |
   | `MISSING_NEEDS_MATLAB_PORT` | No Python equivalent; a named MATLAB script must be ported. |
   | `MISSING_OUT_OF_SCOPE` | The figure isn't model-derived (e.g. apparatus schematic in Fig. 1) and shouldn't be auto-regenerated. |

4. If `MISSING_NEEDS_MATLAB_PORT` or `PRESENT_WITH_MAT_DEPENDENCY`:
   name the MATLAB source (from §4) AND the closer Python reference (notebook
   cell, `utility.py` function) so the port has a concrete model to follow.

## 6. Output of the audit

A single Markdown table per published paper at
`pytest_refactored_codes_v1/docs/REPRODUCIBILITY_AUDIT.md`:

```
| paper | figure / table | description | verdict | python function (existing) | mat dependency | matlab source | python reference (legacy) | notes |
|---|---|---|---|---|---|---|---|---|
| DATA1_main | Fig. 1 | apparatus schematic | MISSING_OUT_OF_SCOPE | — | — | — | — | static SVG; not model-derived |
| DATA1_main | Fig. 2A | mass vs time, dataset 501.1 | PRESENT | plot_sim_comparison | — | — | DiafiltrationPaperPlots.ipynb cell 4 | |
| DATA1_main | Fig. 4 (3 panels) | σ-sensitivity heatmaps | PRESENT_WITH_MAT_DEPENDENCY | render_sigma_sensitivity (reads .mat) | sigma_sensitivity_*.mat | heatmap_sigma_sensitivity_filtration.m + heatmap_sigma_sensitivity_diafiltration.m | diafiltration_plots.py:render_sigma_sensitivity | |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
```

Plus a summary block:

```
Verdict counts per paper:
  DATA1_main:   PRESENT=X  PRESENT_WITH_MAT=Y  NEEDS_PORT=Z  OUT_OF_SCOPE=W
  DATA1_SI:     ...
  DATA2_main:   ...
  DATA2_SI:     ...
  TOTAL:        ...

Distinct MATLAB scripts requiring port (sorted by figure-count impact):
  1. run_sigma_sensitivity.m + heatmap_sigma_sensitivity_*.m  → blocks N figures
  2. doe_heatmap_diafiltration.m + doe_heatmap_filtration.m   → blocks M figures
  3. run_data_analysis.m contour branch                       → blocks K figures
```

## 7. What this audit does NOT do

- Does NOT modify any code in `refactored_codes_v1/`, `legacy/`, or anywhere else.
- Does NOT run `materialize_all` (solver-heavy; 3D fill currently using all 10 cores).
- Does NOT propose ports inline — each `MISSING_NEEDS_MATLAB_PORT` gets a
  separate discussion before any change. (Operator's explicit instruction:
  "Do not change the DATA1 and DATA2 workflow without talking to me clearly
  about it and what you want to change.")
- Does NOT update `target_notebook_source_map.csv` or any other locked
  mapping. The PDFs are the new source of truth.
- Does NOT touch `conductivity_paper.py` under any circumstance.

## 8. Acceptance criteria

The audit is "done" when `REPRODUCIBILITY_AUDIT.md` exists and:

1. Lists every Figure and Table from each of the four PDFs.
2. Assigns each a verdict from the four-element vocabulary.
3. Names a Python function for every `PRESENT` row.
4. Names a MATLAB source AND a Python reference for every
   `MISSING_NEEDS_MATLAB_PORT` row.
5. Names the `.mat` file consumed for every `PRESENT_WITH_MAT_DEPENDENCY` row.
6. Includes the summary block from §6.
7. Has been reviewed by the operator before any port discussion begins.

## 9. Where future work lives

After this audit lands and the operator reviews it, each
`MISSING_NEEDS_MATLAB_PORT` and `PRESENT_WITH_MAT_DEPENDENCY` row spawns a
separate proposal under `pytest_refactored_codes_v1/docs/ports/`:

- `ports/PORT_<MATLAB_script>.md` — what it does, proposed Python signature,
  acceptance test against the original `.mat` output.

No port proposal triggers a code change until the operator explicitly
approves it.
