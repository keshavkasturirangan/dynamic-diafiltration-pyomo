# DATA1 Panel-by-Panel Checklist (Notebook Regeneration)

> **Canonical code path (source of truth):** `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/`  
> Runner entrypoint: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`  
> Core pipeline: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`  
> Legacy root-level scripts (for example `utility.py`, `run_*.py`) are compatibility/reference paths, not primary development paths.


Date: 2026-02-21  
Run ID: `20260221-data1-notebook`  
Artifacts: `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260221-data1-notebook/figures/data1_notebook_regen`

## DATA1 Main

| Target | Expected panels | Regenerated artifact(s) | Status | Notes |
|---|---|---|---|---|
| Fig. 2 | A-D | `mass-dat501.1.png`, `concentration-dat501.1.png`, `mass-dat511.12.png`, `concentration-dat511.12.png` | Generated | Panel mapping consistent with notebook section headers. |
| Table 1 | N/A | N/A (table values from fit outputs; not a PNG artifact) | Separate numeric workflow | Covered via reproduction tables, not notebook plotting. |
| Fig. 3 | Single | `concentration_range.png` | Generated | Matches notebook “Experiment space Figure 3”. |
| Table 2 | N/A | N/A (table values from fit outputs; not a PNG artifact) | Separate numeric workflow | Covered via reproduction tables, not notebook plotting. |
| Fig. 4 | A-C | `sigma_sensitivity-mass.png`, `sigma_sensitivity-reten_conc.png`, `sigma_sensitivity-perme_conc.png` | Generated | Sigma sensitivity set from notebook logic. |
| Fig. 5 | A-C | `contour_fixsig-mass.png`, `contour_fixsig-retentate_conc.png`, `contour_fixsig-permeate_conc.png` | Generated | Contour set from `x_sigma-y_Lp` path. |
| Fig. 6 | A-C | `contour_fixB-mass.png`, `contour_fixB-retentate_conc.png`, `contour_fixB-permeate_conc.png` | Generated | Contour set from `x_B-y_Lp` path. |

## DATA1 SI

| Target | Expected panels | Regenerated artifact(s) | Status | Notes |
|---|---|---|---|---|
| Fig. S2 | all | `cr_measure-dat511.12.png`, `cr_measure-dat501.1.png` | Locked | Locked by visual comparison against extracted SI panel images + notebook section intent. |
| Fig. S3 | all | `mass-dat511.11.png`, `concentration-dat511.11.png` | Locked | Locked by visual comparison against extracted SI panel images + notebook section intent. |
| Fig. S4 | all | `contour_preface_fixsig-mass.png`, `contour_preface_fixsig-retentate_conc.png`, `contour_preface_fixsig-permeate_conc.png` | Locked | Locked by visual comparison against extracted SI panel images + notebook section intent. |
| Fig. S5 | all | `mass-dat501.11.png`, `concentration-dat501.11.png` | Locked | Locked by visual comparison against extracted SI panel images + notebook section intent. |
| Fig. S6 | all | `contour_preface_fixB-mass.png`, `contour_preface_fixB-retentate_conc.png`, `contour_preface_fixB-permeate_conc.png` | Locked | Locked by visual comparison against extracted SI panel images + notebook section intent. |
| Fig. S7 | all | `contour_sensitivity-mass.png`, `contour_sensitivity-retentate_conc.png`, `contour_sensitivity-permeate_conc.png`, `contour_sensitivity-endvial1-*.png`, `contour_sensitivity-endvial5-*.png`, `contour_sensitivity-endvial10-*.png` | Locked | Includes filtration + diafiltration end-vial sensitivity maps; locked by SI page-image comparison. |

## Additional generated notebook artifacts
- `calib_curve.png`
- `colorbar-horizontal.png`
- `colorbar-vertical.png`
- `mass_preface-dat511.12.png`
- `concentration_preface-dat511.12.png`
- `cr_measure-dat301.1.png`

## Lock evidence
- SI figure image extracts used for direct visual lock:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260221-data1-notebook/pdf_extract/data1_si/`
- Regeneration manifest:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260221-data1-notebook/figures/data1_notebook_regen/manifest_data1_notebook_regen.json`

## Remaining check before tolerance sign-off
- Numeric tolerance pass/fail is still separate from figure lock and tracked via side-by-side table outputs.
