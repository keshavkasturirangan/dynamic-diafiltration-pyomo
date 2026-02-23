# DATA1 Notebook vs Published Figures (Initial Comparison)

Date: 2026-02-21
Notebook: `DATA1_matlab/DiafiltrationPaperPlots.ipynb`
Published refs: `/Users/kkasturi/Downloads/DATA1_main.pdf`, `/Users/kkasturi/Downloads/DATA1_SI.pdf`

## Scope
This is an initial artifact-level comparison between:
1. Figures the notebook is intended to generate (from notebook section headers and plotting calls).
2. PNG files currently present in `DATA1_matlab/`.
3. Required DATA1 targets from the validation plan.

This pass is coverage/status only (not yet pixel-perfect visual matching against PDF panels).

## Notebook-declared figure intent
From markdown sections in the notebook:
- Filtration M1: `Figure 2, S5, S6`
- Filtration M2/M3/M4: `Figure S5, S6`
- Sigma sensitivity (filtration M1): `Figure 4`
- Sigma sensitivity across conditions: `Figure S7`
- Diafiltration M1: `Figure 2, 5, 6, S3, S4`
- Diafiltration M2/M3/M4: `Figure 5, 6, S3, S4`
- Sigma sensitivity (diafiltration M1): `Figure 4`
- Experiment space: `Figure 3`
- Concentration-ratio measurement plot: `Figure S1, S2`

## PNG artifacts currently present
Found 16 PNGs in `DATA1_matlab/`:
- `mass-dat501.1.png`
- `concentration-dat501.1.png`
- `mass-dat501.11.png`
- `concentration-dat501.11.png`
- `mass-dat511.12.png`
- `concentration-dat511.12.png`
- `contour_fixB-mass.png`
- `contour_fixB-permeate_conc.png`
- `contour_fixB-retentate_conc.png`
- `contour_fixsig-mass.png`
- `contour_fixsig-permeate_conc.png`
- `contour_fixsig-retentate_conc.png`
- `sigma_sensitivity-mass.png`
- `sigma_sensitivity-perme_conc.png`
- `sigma_sensitivity-reten_conc.png`
- `concentration_range.png`

## Initial target mapping status
- `Figure 2 (main)`: likely covered by `mass-dat501.1.png`, `concentration-dat501.1.png`, `mass-dat511.12.png`, `concentration-dat511.12.png`.
- `Figure 3 (main)`: likely covered by `concentration_range.png`.
- `Figure 4 (main)`: likely covered by 3 sigma-sensitivity plots.
- `Figure 5 (main)`: likely covered by `contour_fixsig-*` set (3 panels).
- `Figure 6 (main)`: likely covered by `contour_fixB-*` set (3 panels).
- `Figure S1/S2/S3/S4/S5/S6/S7 (SI)`: partial at best in current folder; several expected SI-specific outputs are missing from the checked-in PNG set.

## Notable missing outputs (expected from plotting code)
Plotting functions define additional output filenames that are not present now, including examples:
- `cr_measure-dat*.png`
- `calib_curve.png`
- `contour_sensitivity*-mass.png`
- `contour_sensitivity*-retentate_conc.png`
- `contour_sensitivity*-permeate_conc.png`
- `colorbar*-horizontal.png`
- `colorbar*-vertical.png`

## Main conclusion
- The current notebook artifacts appear to cover DATA1 main Figures 2–6 at a first-pass level.
- SI reproduction is incomplete from the currently available PNG outputs.
- Next step should be a deterministic rerun of notebook logic (or script-equivalent) with normalized paths, then a panel-by-panel comparison against DATA1 main/SI PDFs.
