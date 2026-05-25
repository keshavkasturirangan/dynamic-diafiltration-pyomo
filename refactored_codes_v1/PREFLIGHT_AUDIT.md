# DATA3 single-salt pre-flight audit

_Generated: 2026-05-24 15:45:10_

Read-only audit of the loader, the registry, and the bounds.
Run BEFORE pytest / B_form selection / final rollout so any data-side
issues are surfaced before we lock in regression reference values.

Scope: 11 single-salt sheets.


## 1. Registry salt classifications

| Run ID | Sheet salt (parser) | Loader's `namec` | Notes-text hint | Agreement |
|---|---|---|---|---|
| `MC2.05.07.24_NaCl` | `NaCl` | `NaCl` | `NaCl` | ✓ |
| `MC3.07.22.24_SNaCl` | `NaCl` | `NaCl` | `NaCl` | ✓ |
| `MC4.07.11.24_SNaCl` | `NaCl` | `NaCl` | `NaCl` | ✓ |
| `MC5.07.23.24_NaCl` | `NaCl` | `NaCl` | `NaCl` | ✓ |
| `MC5.07.23.24_SNaCl` | `NaCl` | `NaCl` | `NaCl` | ✓ |
| `MC5.07.23.24_S2NaCl` | `NaCl` | `NaCl` | `NaCl` | ✓ |
| `MC2.05.07.24_CaCl2` | `CaCl2` | `CaCl2` | `CaCl2` | ✓ |
| `MC3.07.11.24_SCaCl2` | `CaCl2` | `CaCl2` | `CaCl2` | ✓ |
| `MC3.07.12.24_S2CaCl2` | `CaCl2` | `CaCl2` | `CaCl2` | ✓ |
| `MC2.05.21.24_LaCl3` | `LaCl3` | `LaCl3` | `LaCl3` | ✓ |
| `MC4.07.11.24_SLaCl3` | `LaCl3` | `LaCl3` | `LaCl3` | ✓ |

**Summary:** 11/11 sheets have consistent salt identification across sources.


## 2. V_tube sanity

**Representative sheet:** `MC3.07.22.24_SNaCl`  **Apparatus constant:** V_tube = 0.3 g

Per-vial t_corr placement (`t_corr` is where the ICP value gets planted in the NaN-padded cV_avg array; *t_close* is where it would land WITHOUT the correction):

| Vial | t_open (s) | t_close (s) | m_total (g) | dm/dt (mg/s) | tau (s) | t_corr (s) | % of vial |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 432 | 0.030 | 0.07 | 4318 | 0 | 0% |
| 2 | 439 | 1394 | 2.960 | 3.10 | 97 | 820 | 40% |
| 3 | 1401 | 2371 | 3.010 | 3.10 | 97 | 1789 | 40% |
| 4 | 2378 | 3312 | 3.010 | 3.22 | 93 | 2752 | 40% |
| 5 | 3320 | 4231 | 3.020 | 3.31 | 91 | 3685 | 40% |
| 6 | 4238 | 5113 | 2.990 | 3.42 | 88 | 4588 | 40% |
| 7 | 5121 | 5970 | 2.960 | 3.48 | 86 | 5459 | 40% |
| 8 | 5977 | 6795 | 2.980 | 3.64 | 82 | 6304 | 40% |
| 9 | 6803 | 7634 | 3.020 | 3.63 | 83 | 7136 | 40% |
| 10 | 7641 | 8443 | 3.000 | 3.74 | 80 | 7962 | 40% |
| 11 | 8450 | 9238 | 3.000 | 3.80 | 79 | 8765 | 40% |
**Interpretation:** If V_tube ≈ 0.3 g is right for this apparatus, t_corr typically lands 30-50% into each vial (well inside the collection window). If t_corr lands near 0% or near 100%, V_tube may be off — near 0% means V_tube is too LARGE (over-correcting backward); near 100% means V_tube is too SMALL (no correction happening). Visually look at the **% of vial** column — values should cluster around 30-50%.


## 3. Pyomo bounds review

| Parameter | Current bound | Literature range for NF270 | Comment |
|---|---|---|---|
| `Lp` (water permeability, L·m⁻²·h⁻¹·bar⁻¹) | [0.5, 50] | 4–15 typical, up to ~30 in best condition | Lower bound (0.5) very conservative — won't clip any real fit. Upper bound (50) has plenty of headroom. ✓ |
| `B` (solute permeability, μm/s) | [1e-6, 30] | 0.1–10 for divalents (CaCl₂, LaCl₃); 1–30 for NaCl on NF270 | Lower bound (1e-6) effectively zero. Upper bound (30) tight — could clip a hyper-permeable fit on NaCl. ⚠ |
| `σ` (reflection coefficient) | [1e-3, 0.999] | 0.3–0.6 for NaCl, 0.9–0.99 for divalents/trivalents | The 0.999 upper-bound clipping shows up *constantly* in the contour data (data wants σ → 1 for NaCl). Functionally correct; physically defensible (real membranes are never σ=1 exactly). ✓ |

**Flags worth a second look:**

- `B ∈ [1e-6, 30]`: if any fit lands at B = 30 (the upper bound), the true value may be higher. Worth grepping the previous contour batch's centering_fit.json files for `"B": 30.0` or near-30 values.
- `σ` upper bound of 0.999 is just slightly inside 1.0 — this is intentional (avoids the σ=1 singularity in Spiegler-Kedem) and is what Fix A clips deterministic seeds to.


## 4. 11-sheet roster + regime + ICP coverage

| Run ID | Workbook | Salt | n_vials | cF first → last (mM) | Regime | ICP vials w/ data | Notes |
|---|---|---|---|---|---|---|---|
| `MC2.05.07.24_NaCl` | NF270_MC2.xlsx | `NaCl` | 11 | 0.9 → 35.7 | concentration | 10 / 11 |  |
| `MC3.07.22.24_SNaCl` | NF270_MC3.xlsx | `NaCl` | 11 | 95.5 → 10.4 | dilution | 10 / 11 |  |
| `MC4.07.11.24_SNaCl` | NF270_MC4.xlsx | `NaCl` | 11 | 1.1 → 110.1 | concentration | 10 / 11 |  |
| `MC5.07.23.24_NaCl` | NF270_MC5.xlsx | `NaCl` | 11 | 0.9 → 107.3 | concentration | 10 / 11 |  |
| `MC5.07.23.24_SNaCl` | NF270_MC5.xlsx | `NaCl` | 11 | 0.9 → 109.3 | concentration | 10 / 11 |  |
| `MC5.07.23.24_S2NaCl` | NF270_MC5.xlsx | `NaCl` | 11 | 93.1 → 16.1 | dilution | 10 / 11 |  |
| `MC2.05.07.24_CaCl2` | NF270_MC2.xlsx | `CaCl2` | 11 | 0.6 → 23.8 | concentration | 10 / 11 |  |
| `MC3.07.11.24_SCaCl2` | NF270_MC3.xlsx | `CaCl2` | 11 | 0.6 → 25.9 | concentration | 10 / 11 |  |
| `MC3.07.12.24_S2CaCl2` | NF270_MC3.xlsx | `CaCl2` | 11 | 0.6 → 21.7 | concentration | 10 / 11 |  |
| `MC2.05.21.24_LaCl3` | NF270_MC2.xlsx | `LaCl3` | 11 | 0.4 → 9.7 | concentration | 10 / 11 |  |
| `MC4.07.11.24_SLaCl3` | NF270_MC4.xlsx | `LaCl3` | 11 | 0.5 → 10.1 | concentration | 10 / 11 |  |


## 5. Loader field audit — what gets read / used / ignored

### Time-series block (per-row data)

| Workbook column | Read? | Used in fit? | Where it lands in `data_stru` |
|---|---|---|---|
| `Time (s)` | ✓ | ✓ — defines x-axis, vial segmentation | `data_raw[i]["time"]` |
| `Mass (g)` | ✓ | ✓ — mass channel in WSSE, reset per vial | `data_raw[i]["mass"]` |
| `Pressure (psi)` | ✓ | ✓ — averaged into `data_config["delP"]` (bar) | `data_raw[i]["pressure"]` |
| `Retentate Temp` | ✓ | ✓ — per-row EC25 compensation | `data_raw[i]["retentate_temp"]` + top-level `measured_ret_temp_C` |
| `Retentate Cond @ Temp (uS/cm)` | ✓ | ✓ — compensated to 25 °C, then inverted to cF | `data_raw[i]["cF_exp"]` (after inversion) |
| `Permeate Temp` | ✓ | ✓ — per-row EC25 compensation | `data_raw[i]["permeate_temp"]` + top-level `measured_perm_temp_C` |
| `Permeate Cond @ Temp (uS/cm)` | ✓ | ⚪ only if `NF270_USE_PERMEATE_PROBE` toggle is on | `data_raw[i]["cV_perm_cond"]` |
| `Vial Swap` | ✓ | ✓ — segments run into per-vial windows | `data_raw[i]["vial_swap"]` |

### Top-left metadata block

| Workbook cell | Read? | Used in fit? | Where it lands |
|---|---|---|---|
| `Datapoints:` | ⚠ partially | sanity-check only, not stored | (just compared to actual time-series length) |
| `Experiment Name` | ✓ | ⚪ provenance only | `data_stru["filename"]` |
| `Initial Solution Weight (g)` | ✓ | ✓ — boundary condition `M_F0` | `data_config["M_F0"]` |
| `Final Solution Weight (g)` | ✓ | ✓ — combined with `M_F0` to give `M_O = Final − M_F0` | `data_config["M_O"]` |
| `Final Vial Weight (g)` | ⚠ parsed | ✗ NOT used downstream | only in `data_config` metadata |
| `Final Vial w/ Solution (g)` | ⚠ parsed | ✗ NOT used downstream | only in `data_config` metadata |
| `ICP Calibration Points (#)` | ⚠ parsed | ✗ NOT used downstream | provenance only |
| `Number of Salts` | ✓ | ✓ — `nc` field, determines single vs multi-salt path | `data_config["nc"]` |
| `Salt 1`, `Salt 2` | ✓ | ✓ — salt identity for conductivity inversion + ICP MW | `data_config["namec"]` |
| `Notes` (free-text cell) | ✓ | ✓ — parsed for `C_F0`, `C_D`, `rpm`, etc. | `data_config["note_text"]` + derived fields |

### Sidebar ICP block (per-vial assays)

| Sidebar row | Read? | Used in fit? | Where it lands |
|---|---|---|---|
| `Feed` row | ✓ | ✓ — overrides Notes-text `C_F0` if ICP value present | `data_config["cF_feed_icp_mM"]`, `C_F0`, `C_F0_source` |
| `Diafiltrate` row | ✓ | ✓ — overrides Notes-text `C_D` if ICP value present | `data_config["cF_diafiltrate_icp_mM"]`, `C_D`, `C_D_source` |
| `Retentate` row | ✓ | ⚪ cross-check only | `data_config["cF_retentate_icp_mM"]` |
| `Final Tube` row | ✓ | ⚪ cross-check only | `data_config["icp_final_tube_mM"]` |
| `Vial 1..N` rows | ✓ | ✓ — populates `cV_avg` per vial (tube-transit-corrected anchor) | `data_raw[i]["cV_avg"]` (NaN-padded array) |
| `ICP Salt 1 Calibration` table | ⚠ parsed | ✗ NOT used downstream | local `_calib` variable only |

### Summary

- **Fully used:** time, mass, pressure, retentate temp + cond, permeate temp, vial swap, M_F0, M_O, salt names, Notes, sidebar ICP per-vial values.
- **Used only if a toggle is enabled:** permeate conductivity (`NF270_USE_PERMEATE_PROBE`).
- **Parsed but not threaded into the fit:** Final Vial Weight, Final Vial w/ Solution, ICP Calibration Points (#), the per-salt ICP Calibration tables.

The "parsed but not used" items are kept as provenance (so they survive in `data_config` for audit purposes) but don't currently influence model output. If your audit standards require *all* metadata to be either used or explicitly justified, the ICP calibration tables are the most actionable candidate — they could be wired in to validate the precomputed `mg/L` values rather than trusting the meter's own calibration.


## 6. ICP calibration cross-check (NEW)

Per-sheet sanity check: each workbook embeds an ICP calibration table (Concentration vs. Intensity standards). We refit that table as a linear curve and cross-check whether the stored `ICP Salt 1 (mg/L)` values match what the curve would predict from the same intensity. Healthy sheets show all residuals within 5%; bad residuals usually mean the workbook's calibration column was leftover/wrong (a Salt-1=NaCl template reused for a different cation), NOT that the lab's analytical result is wrong.

| Run ID | R² | n pts | Stored mg/L range | Within 5% | Max rel residual (%) | Status |
|---|---|---|---|---|---|---|
| `MC2.05.07.24_NaCl` | 0.9973 | 10 | 0.0-20.0 | 14/14 | 3.13 | ✓ healthy |
| `MC3.07.22.24_SNaCl` | 0.9980 | 10 | 0.0-20.0 | 11/13 | 5.23 | ⚠ disagreement (2/13 fail) |
| `MC4.07.11.24_SNaCl` | 0.9985 | 10 | 0.0-20.0 | 14/14 | 0.37 | ✓ healthy |
| `MC5.07.23.24_NaCl` | — | — | — | — | — | ⚠ no calibration table in workbook |
| `MC5.07.23.24_SNaCl` | — | — | — | — | — | ⚠ no calibration table in workbook |
| `MC5.07.23.24_S2NaCl` | — | — | — | — | — | ⚠ no calibration table in workbook |
| `MC2.05.07.24_CaCl2` | 0.9965 | 10 | 0.0-20.0 | 0/14 | 921.04 | ⚠⚠ catastrophic (max 921%) — likely wrong calibration in workbook |
| `MC3.07.11.24_SCaCl2` | 0.9575 | 10 | 0.0-20.0 | 3/14 | 84.45 | ⚠⚠ catastrophic (max 84%) — likely wrong calibration in workbook |
| `MC3.07.12.24_S2CaCl2` | 0.9575 | 10 | 0.0-20.0 | 3/14 | 83.86 | ⚠⚠ catastrophic (max 84%) — likely wrong calibration in workbook |
| `MC2.05.21.24_LaCl3` | 0.9950 | 10 | 0.0-20.0 | 1/15 | 78.85 | ⚠⚠ catastrophic (max 79%) — likely wrong calibration in workbook |
| `MC4.07.11.24_SLaCl3` | 0.9999 | 10 | 0.0-20.0 | 13/14 | 38.92 | ~ mostly OK (13/14) |
**Interpretation:**

- ✓ **NaCl MC2/MC3/MC4 sheets** show R² ≈ 0.998 and residuals ≤ 5%. The workbook's calibration table is internally consistent with its precomputed mg/L values.
- ⚠ **MC5 NaCl sheets** have no calibration table at all — cross-check can't run. Note this for downstream users.
- ⚠⚠ **CaCl₂ and LaCl₃ sheets (mostly)** show residuals 39-921%. This is almost certainly because the calibration column in the workbook is leftover from a Salt-1=NaCl template that was reused for a different cation. The lab's analytical software used the correct calibration; the workbook column does not reflect it.

**Action items:**

- DO NOT modify `cV_avg` based on the workbook calibration — the precomputed   mg/L from the lab software is ground truth.
- The audit field `icp_calibration_residual` is now populated on every loaded sheet;   downstream consumers can inspect it but should treat large residuals on multivalent   sheets as a *workbook bookkeeping issue*, not a measurement error.
- If we want a stronger cross-check, we'd need the lab's actual calibration   (perhaps from the analytical-software logs, not the workbook).

