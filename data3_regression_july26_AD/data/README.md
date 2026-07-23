# Data provenance

Experimental datasets for the DATA3 regression, copied **2026-07-05** from the
shared Box folder unless noted. Box paths below are relative to the shared-folder
root `Diafiltration Modeling and Experiments (Phillip and Dowling groups)/`.

These are working copies; the Box folder remains the system of record. If the Box
files are updated, re-copy and note the date here.

| File | Source | Notes |
|------|--------|-------|
| `Rejection_Analysis.xlsx` | From **Bill Phillip via email** (attached alongside his MATLAB script `B_regression_concentration_NAClanal.m`). | The file the analysis code currently reads (`analysis/step2_nacl/regress_nacl.py`). Content twin of `BoE Analysis.xlsx`; analysis uses sheet `NF270_MC5 05.27.26_NaCl (2)`. Previously stored under `prior_analysis/`; consolidated here on 2026-07-05. |
| `BoE Analysis.xlsx` | Box `Multicomponent Diafiltration/Analyses/BoE Analysis.xlsx` | Authoritative working analysis workbook (27 sheets: single-salt NaCl/CaCl₂/LaCl₃ + multicomponent `E#` experiments). Processed columns B/H/J/K/M/N/O per experiment. |
| `NF270_MC2.xlsx` | Box `Multicomponent Diafiltration/Diafiltration Campaign and Raw Data/NF270_MC2.xlsx` | Raw campaign data, membrane cut MC2 (28-col raw schema; `Jw` not pre-computed). |
| `NF270_MC3.xlsx` | Box `.../Diafiltration Campaign and Raw Data/NF270_MC3.xlsx` | Raw campaign data, membrane cut MC3. |
| `NF270_MC4.xlsx` | Box `.../Diafiltration Campaign and Raw Data/NF270_MC4.xlsx` | Raw campaign data, membrane cut MC4. |
| `NF270_MC5.xlsx` | Box `.../Diafiltration Campaign and Raw Data/NF270_MC5.xlsx` | Raw campaign data, membrane cut MC5 (contains the analyzed `05.27.26_NaCl` runs). |
| `NF270 - Experimental Campaign.xlsx` | Box `.../Diafiltration Campaign and Raw Data/NF270 - Experimental Campaign.xlsx` | Campaign index / experiment log across membrane cuts. |

## Notes
- Membrane cuts `MC#` are different coupons cut from the same NF270 flat sheet.
- Available single-salt NaCl experiments and their status (processed vs. raw-only)
  are tabulated in `docs/reports/main.tex` (Table `tab:nacldata`). One raw tab,
  `NF270_MC5 05.26.26_NaCl`, is **mislabeled** — its own note identifies it as a
  CaCl₂ run.
- The raw `NF270_MC#.xlsx` workbooks need the interfacial-concentration / water-flux
  processing step before the current pipeline (which expects processed columns
  B/H/J) can read them.
