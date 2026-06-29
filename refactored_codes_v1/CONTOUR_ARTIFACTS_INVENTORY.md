# DATA3 contour artifacts & scripts — locations

All absolute paths share two roots:
- **`REPO_ROOT`** = `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo`
- **`SCRIPTS/`** = `REPO_ROOT/refactored_codes_v1/`
- **`nf270/`** = `REPO_ROOT/UnifiedFramework/DATA3/results/paper_artifacts/nf270/`

The 11 campaign run-ids: `MC2.05.07.24_CaCl2`, `MC2.05.07.24_NaCl`, `MC2.05.21.24_LaCl3`, `MC3.07.11.24_SCaCl2`, `MC3.07.12.24_S2CaCl2`, `MC3.07.22.24_SNaCl`, `MC4.07.11.24_SLaCl3`, `MC4.07.11.24_SNaCl`, `MC5.07.23.24_NaCl`, `MC5.07.23.24_S2NaCl`, `MC5.07.23.24_SNaCl`.

---

## 1. TL;DR map — pipeline stage → script → output

| Stage | Script (`SCRIPTS/`) | Output location (under `nf270/`) | Kind |
|---|---|---|---|
| 3D objective grid (Lp×B×σ) | `_run_3d_contour_concentrating.py` | `contour3d/<run_id>/grid3d.csv` | CSV (6-col) |
| σ×B at fixed Lp (canonical, no floor) | `_run_Bsigma_fixedLp.py` (`2pct_nofloor`) | `bform_study/Bsigma_fixedLp_nofloor/<run_id>/single/` | CSV + JSONL |
| σ×B at fixed Lp (floored variant) | `_run_Bsigma_fixedLp.py` (`2pct_floor1mm`) | `bform_study/Bsigma_fixedLp/<run_id>/<form>/` | JSONL only |
| Profile/slice contours (re-fit B(c)) | `_run_profile_contours.py` | `bform_study/profile_contours/…`, `bform_study/stack_contours/…` | JSONL only |
| 5-channel all-pairs slice | `_run_all_pairs_5channel.py` | `matlab_ports_5ch_allpairs/<run_id>/` | CSV (7-col) |
| β / ionic-strength contour (library) | `_run_beta_contour.py` | `beta_contours/<run_id>/` | CSV (5-col) |
| β / ionic-strength contour (fast MP) | `_run_beta_contour_fast.py` | `beta_contours/` (flat) | CSV (3-col) |
| Band-masked σ×Lp | `_run_band_contour.py` | `band_value_test/` | CSV (3-col) |
| Constant-B pairwise panels (library) | `_run_bform_contours.py` | `bform_study/contours/<run_id>/` | CSV (5-col) |
| 2-D coarse panels (legacy) | (library / earlier driver) | `contour_panels/<run_id>/` | CSV (5-col) |
| MATLAB-replica 5-ch full sweep | `_run_all_pairs_5channel.py` / drivers | `matlab_ports_5ch_full/<run_id>/` | CSV (7-col) |
| 3D contour GIFs | `_make_3d_contour_animations.py` | `animations3d/<run_id>/` | GIF |
| DATA1 identifiability panels | `_render_data1_panels.py` | `bform_study/Lp_identifiability_panels/<run_id>_identifiability.png` | PNG |
| Decks | `_build_*deck*.py` | `SCRIPTS/*.pptx` | PPTX |
| Fit cache (warm-start θ*) | — | `warm_start_fits/summary.json` | JSON |
| Fit cache (per-form θ*) | `_run_bform_campaign.py` | `bform_study/<run_id>/result_<form>.json` | JSON |

---

## 2. Generation scripts (`SCRIPTS/`)

| Script | Purpose | Command | Output path(s) | Plane / method |
|---|---|---|---|---|
| `_run_3d_contour_concentrating.py` | 30×30×10 Lp×B×σ 3-D sweep, 11 sheets; port of `calc_contour_3d.m` | `python3 SCRIPTS/_run_3d_contour_concentrating.py [--smoke-test] [--workers N] [--sheets …] [--grid-lp/b/sigma …] [--nfe 80] [--no-render]` | `nf270/contour3d/<run_id>/grid3d.csv` (+`_meta.json`, `_failures.log`, `<run_id>_grid3d_slices.png`) | Lp×B×σ **slice** (forward-sim per cell, no re-fit) |
| `_run_Bsigma_fixedLp.py` | σ×B identifiability with Lp pinned at Lp* | `python3 _run_Bsigma_fixedLp.py {exp\|stack\|campaign\|single} …`; `BSIGMA_VARIANT=2pct_floor1mm …` | `nf270/bform_study/Bsigma_fixedLp_nofloor/<run_id>/<form>/` (canonical, CSV) or `…/Bsigma_fixedLp/<run_id>/<form>/` (floored, JSONL); pooled `…/_stack/<group>/<form>/` | σ×B at **fixed Lp\*** **slice** (B swept by scaling form coeff `k`); `SLICE_FORMS={single,sat,donnan}` only |
| `_run_profile_contours.py` | PROFILE contours (re-solve B(c)+S at each node) + slice + stacks | `python3 _run_profile_contours.py {plan\|run\|plot\|campaign\|stack\|bformsweep\|stackplot} …` | `nf270/bform_study/profile_contours/<rid>/<form>/<plane>/`; stacks `…/stack_contours/<group>/<form>/<plane>/<method>/` | sigLp / BLp / sigB; **profile** (re-fit, `solve_model_B_fix`) **and** slice |
| `_run_all_pairs_5channel.py` | All 3 pairwise sweeps, 5-channel kernel | `python3 SCRIPTS/_run_all_pairs_5channel.py [--sheet …] [--grid N] [--workers N] [--nfe 80]` | `nf270/matlab_ports_5ch_allpairs/<run_id>/contourdata-x_*-y_*.csv` (+`objcontour-*.png`, `_meta.json`) | (σ,Lp)/(B,Lp)/(σ,B) **slice**; Y-priority Lp>B>σ |
| `_run_beta_contour.py` | β-form (B=β₀+β₁·cIn) contour via library | `python3 _run_beta_contour.py [RUN_ID] [GRID=18] [NFE=80]` | `nf270/beta_contours/<RUN_ID>/contourdata-x_{beta_0\|sigma}-y_{Lp\|beta_1}.csv` (+`objcontour-*.png`, `centering_fit.json`) | (β₀,Lp)/(σ,Lp)/(σ,β₁) **slice** |
| `_run_beta_contour_fast.py` | Capped-solver MP β contour | `python3 _run_beta_contour_fast.py [RUN_ID] [GRID=18]` | `nf270/beta_contours/contourdata-{beta_0-Lp\|sigma-beta_1}-<RUN_ID>.csv` (flat, +`.png`/`.json`) | (β₀,Lp) & (σ,β₁) **slice**; β₁ feasibility-bounded |
| `_run_band_contour.py` | Band-masked σ×Lp (per-band θ, band vials only) | `python3 _run_band_contour.py [RUN_ID] [GRID=16]` | `nf270/band_value_test/band_masked_contourdata-<RUN_ID>-band{a}_{b}.csv` (+`.png`/`.json`) | σ×Lp **slice**, per cF-band |
| `_run_bform_contours.py` | 3 constant-B pairwise WSSE maps (seeding workflow) | `python3 _run_bform_contours.py [grid=30] [workers=6] [sheets=all\|<rid,…>]` | `nf270/bform_study/contours/<run_id>/contourdata-x_*-y_*.csv` (+`objcontour-*.png`) | (B,Lp)/(σ,Lp)/(σ,B) **slice** |

**Not grid generators** (orchestrate fits/figures only): `_run_band_campaign.py`, `_run_band_wrapper.py`, `_run_band_value_test.py` (only *reads* an existing `matlab_ports_5ch_full/…/contourdata-x_sigma-y_Lp.csv`), `_run_bform_campaign.py` (fits the B-forms; exports `load(rid)` + `OUT`).

**Common weighting in drivers:** `_run_3d_contour_concentrating.py` sets `_lib.NF270_CF_RESIDUAL_SCALE_FRACTION=0.02` + `NF270_CF_RESIDUAL_FLOOR_MM=None` per worker (canonical 2%/no-floor). `_run_profile_contours.py`, `_run_beta_contour*.py`, `_run_band_contour.py`, `_run_bform_contours.py` set `NF270_CF_RESIDUAL_FLOOR_MM=1.0` (2%+floor). `_run_Bsigma_fixedLp.py` switches on `BSIGMA_VARIANT` (default `2pct_nofloor`).

---

## 3. Library entry points — `SCRIPTS/refactored_ucb_library.py`

Section header documenting the DATA1 `calc_contour_2d.m` / `plot_contour.m` workflow: **lines 8020–8043**.

| Function | Line | Role |
|---|---|---|
| `_nf270_contour_objectives_at_theta(data_stru, theta, *, mode="Lag", B_form="single", workflow_family="DATA3", nfe=80, band_vials=None)` | **8045** | DATA1 `calc_contour_2d` **cell kernel**: forward-sim one sheet at fixed θ via `solve_model(sim_opt=True)`; returns `(obj_m, obj_cv, obj_cr)` **raw** (caller logs), failed solve → `(nan,nan,nan)`. Called by `_run_3d_contour_concentrating.py._eval_cell`. |
| `_nf270_contour_grid_dataframe(... x_var, y_var, *, grid_density=20, …, band_vials=None)` | **8093** | Builds one 2-D **slice** grid (double loop, third axis held at `theta_fit`). Schema `[x_var, y_var, Obj_mass, Obj_concentration, Obj_retentate_concentration]`, objectives **log10** (≤0/NaN→NaN). Axes: `Lp/B/sigma/beta_0/beta_1/beta_2`. |
| `run_nf270_contour_for_sheet(run_id, *, save_dir, grid_density=20, theta_fit=None, …, B_form="single", band_vials=None)` | **8200** | Top-level: optional centering fit (`centering_fit.json`), 3 slices → `contourdata-x_{x}-y_{y}.csv` + `objcontour-*.png` + paper-style plots; per-CSV checkpoint/resume. |

**Axis ranges** (`_axis_values`, lines 8132–8150): `Lp→linspace(0.1·Lp_fit, 2.0·Lp_fit, n)`; `sigma→linspace(0,1,n)` (never log); `B→linspace(0.1·|B_fit|, 3.0·|B_fit|, n)` (`B_abs=max(|B_fit|,1e-4)`); `beta_0→linspace(0.1·β₀, 3.0·β₀, n)`; `beta_1/2→linspace(βk−2·span, βk+2·span, n)`, `span=max(|βk|,1)`.

**`sweep_pairs`** (lines 8264–8275), Y-priority **Lp > B(β) > σ**:
- lumped-B (`single`/string/fractional): `(B,Lp)`, `(sigma,Lp)`, `(sigma,B)`.
- numeric `B_form ≥ 1`: `(beta_0,Lp)`, `(sigma,Lp)`, `(sigma,beta_1)`.

**Objective / weighting constants:**

| Constant | Line | Current value | Legacy fallback |
|---|---|---|---|
| `NF270_CF_RESIDUAL_SCALE_FRACTION` | **299** | `0.02` (2%) | `0.003` (0.3%) |
| `NF270_CP_RESIDUAL_SCALE_FRACTION` | **300** | `None` | `0.03` (3%) |
| `NF270_CF_RESIDUAL_FLOOR_MM` | **271** | `None` (disabled) | n/a (DATA3 guard) |

- `_nf270_conc_scales(workflow_family)` — **lines 303–313**. DATA3 → cF=2%, cP=3% (cP falls back since its knob is `None`).
- `obj_rule(m)` — **line 2548** (in `solve_model`); duplicate at **line 3440** (in `solve_model_B_fix`). Mass `mass_scale=data_config["mass_scale_g"]` default **0.01 g** (lines 2554–2556, residual line 2625); permeate `(res_cp/(_cp_frac·cp_meas))²` (lines 2638, 2654); retentate `relative_scale=_cf_frac·cf_meas` (line 2685, residual 2691), with optional `max(_cf_frac·cf_meas, floor)` (lines 2674–2689) when `NF270_CF_RESIDUAL_FLOOR_MM` set and DATA3 (inert at `None`). Same fractions feed `calc_FIM` (lines 3113, 3130–3131).
- `solve_model(...)` — **line 2426**: `sim_opt=False` = fit (minimized `Objective`, line 2842); `sim_opt=True` = forward eval (dummy `Obj_1=Objective(expr=1)`, `obj_rule` as `Expression`, lines 2838–2841). Related: `solve_model_per_concentration_band` (**line 3238**), `solve_model_B_fix` (**line 3374**, B held fixed for slices).
- Original MATLAB port (not the NF270 port): `calc_contour_2d_py` (**line 14424**, CSV writer **14553**); renderers `plot_contour` (**6348**), `_plot_contour_paper_style` (**6475**), `_plot_heatmap_frame` (**6295**), `plot_contour_py` (**14558**).
- Shared v7 line-contour renderer `_contour_panel(...)` — **line 328 of `SCRIPTS/_run_profile_contours.py`** (called line 295); mirrored by `_render_data1_panels.py._panel` (line 45).

---

## 4. CSV data files

### 4.1 3-D objective grids — `contour3d/<run_id>/grid3d.csv`
Path glob: `nf270/contour3d/<run_id>/grid3d.csv` (+`_meta.json`, `<run_id>_grid3d_slices.png`)
**Schema (6 cols, objectives log10):** `Lp,B,sigma,Obj_mass,Obj_concentration,Obj_retentate_concentration`

| run_id | grid (Lp×B×σ) | rows | status |
|---|---|---|---|
| MC2.05.07.24_CaCl2 | 30×30×10 | 9000 | full |
| MC2.05.07.24_NaCl | 15×12×10 | 1800 | reduced (has cf fields) |
| MC2.05.21.24_LaCl3 | 30×30×10 | 9000 | full |
| MC3.07.11.24_SCaCl2 | 30×30×10 | 9000 | full |
| MC3.07.12.24_S2CaCl2 | 30×30×10 | 9000 | full |
| MC3.07.22.24_SNaCl | 3×3×3 | 27 | **smoke test** |
| MC4.07.11.24_SLaCl3 | 30×30×10 | 9000 | full |
| MC4.07.11.24_SNaCl | 30×30×10 | 9000 | full |
| MC5.07.23.24_NaCl | 30×30×10 | 9000 | full |
| MC5.07.23.24_SNaCl | 30×30×10 | 9000 | full |
| MC5.07.23.24_S2NaCl | — | — | **empty dir, no files** |

**Weighting:** Only `contour3d/MC2.05.07.24_NaCl/_meta.json` carries `cf_residual_scale_fraction=0.02`, `cf_residual_floor_mM=null`, `weighting="retentate cF = 2% relative, NO floor (Lilonfe et al.); permeate cV 3%; mass 0.01 g"`. All other `_meta.json` lack these keys (implicitly legacy 0.3%). `MC5.07.23.24_S2NaCl` has no `_meta.json`.

### 4.2 Fixed-Lp σ×B — no-floor (CSV) vs floored (JSONL only)
- **No-floor (canonical, CSV present):** `nf270/bform_study/Bsigma_fixedLp_nofloor/<run_id>/single/`
  - Only **2 run-ids**: `MC2.05.07.24_CaCl2` (full: `grid.json, nodes.jsonl, contour.png, WORKER_DONE,` **`contourdata-x_sigma-y_B.csv`**) and `MC3.07.22.24_SNaCl` (`grid.json, nodes.jsonl` only — CSV not yet emitted).
  - **`contourdata-x_sigma-y_B.csv`** schema (3-channel, log10): `sigma,B,Obj_mass,Obj_concentration,Obj_retentate_concentration`; **49 rows (7×7)**; `grid.json` has `axis_sigma[7], axis_B[7]`.
- **Floored (JSONL only, NO CSV):** `nf270/bform_study/Bsigma_fixedLp/<run_id>/<form>/{grid.json,nodes.jsonl,contour.png,WORKER_DONE}`
  - `<form>`∈`{single,sat,donnan}`; most runs `single`+`sat`; `MC5.07.23.24_S2NaCl` adds `donnan`; `MC3.07.22.24_SNaCl, MC4.07.11.24_SLaCl3, MC5.07.23.24_NaCl, MC5.07.23.24_SNaCl` are `single` only.
  - `grid.json`: `rid, form, Lp_fixed, Bfeed_fit, sigma_fit, axis_sigma[9], axis_B[9]` (9×9; donnan=69).
  - `nodes.jsonl` (81 nodes): `i,j,sigma,B,status,Bfeed,obj_m,obj_cv,obj_cr,WSSE3,dt` → **WSSE3 = pooled 3-channel objective**.
  - Pooled: `…/_stack/<group>/<form>/stacked.{json,png}`; `<group>`∈`{CaCl2,LaCl3,NaCl_concentrating,NaCl_diluting}`, `<form>`∈`{single,sat}`; `stacked.json`: `group, form, Lp, pooled_min{sigma,B,WSSE3}, per_exp_min[]`.

### 4.3 Profile / stack / per-experiment contours (JSONL only, NO CSV)
- **Stack contours:** `nf270/bform_study/stack_contours/<group>/<form>/<plane>/{profile|slice}/<run_id>/{grid.json,nodes.jsonl[,.skip],WORKER_DONE}`; pooled `…/<profile|slice>/stacked.{json,png}` one level up.
  - `<group>`∈`{CaCl2,LaCl3,NaCl_concentrating,NaCl_diluting}`; `<form>`∈`{single,sat,donnan,poly1,poly2,poly3}`; `<plane>`∈`{sigLp,BLp,sigB}`; method∈`{profile,slice}`.
  - `grid.json`: `rid,form,plane,xname,yname,axisA[9],axisB[9],anchor{Lp,beta_0,beta_1,sigma,S0,S},anchor_ij,out_dir,skip_sim_init,method`.
  - `nodes.jsonl`: `i,j,x,y,status,obj_m,obj_cv,obj_cr,WSSE3,params{…},dt`.
  - `stacked.json`: `group,form,plane,method,stacked_min{Lp,sigma,WSSE3},per_exp_min[]`.
- **Profile contours (per-experiment):** `nf270/bform_study/profile_contours/<run_id>/<form>/<plane>/{grid.json,nodes.jsonl[,.skip,.cur],contour.png,run.log,WORKER_DONE}`
  - **4 run-ids**: `MC2.05.07.24_CaCl2, MC2.05.07.24_NaCl, MC2.05.21.24_LaCl3, MC3.07.22.24_SNaCl`. `<form>`∈`{single,sat,donnan,poly1,poly2,poly3}`; `<plane>`∈`{sigLp,BLp,sigB}` (BLp/sigB only under `single`). Same schema as stack contours.

### 4.4 DATA1 per-salt panels
- **Coarse 2-D panels (10×10):** `nf270/contour_panels/<run_id>/contourdata-x_<X>-y_<Y>.csv` — **11 run-ids × 3 CSVs = 33**.

  | filename | schema (5-col, log10) | rows |
  |---|---|---|
  | `contourdata-x_B-y_Lp.csv` | `B,Lp,Obj_mass,Obj_concentration,Obj_retentate_concentration` | 100 (10×10) |
  | `contourdata-x_sigma-y_B.csv` | `sigma,B,Obj_mass,Obj_concentration,Obj_retentate_concentration` | 100 |
  | `contourdata-x_sigma-y_Lp.csv` | `sigma,Lp,Obj_mass,Obj_concentration,Obj_retentate_concentration` | 100 |

  Co-located per run: `centering_fit.json` (`theta_fit`: `Lp,B,sigma,S0,S` + `obj_m/obj_cv/obj_cr`), `objcontour-*.png`, `paper-*-{mass,permeate_conc,retentate_conc}.png`. No per-CSV weighting field. Archive siblings at `nf270/` root: `contour_panels.zip`, `contour_panels_cF_floor_1mM.zip` (a 1 mM-floor variant was archived).
- **Fine 2-D panels (30×30):** `nf270/bform_study/contours/<run_id>/contourdata-x_<X>-y_<Y>.csv` — same 3 filenames/schema as above but **900 rows (30×30)**; **11 run-ids × 3 = 33**. Co-located `objcontour-*.png`, `paper-*-*.png`. No per-run weighting JSON.

### 4.5 β / ionic-strength contours — `beta_contours/`
- **Library driver (`_run_beta_contour.py`), per-rid subdir:** `nf270/beta_contours/<RUN_ID>/contourdata-x_{beta_0|sigma}-y_{Lp|beta_1}.csv` — 3-channel 5-col schema; planes (β₀,Lp)/(σ,Lp)/(σ,β₁). Co-located `objcontour-*.png`, `centering_fit.json`.
- **Fast MP driver (`_run_beta_contour_fast.py`), flat:** distinct naming `contourdata-<X>-<Y>-<run_id>.csv` —

  | filename | schema (single-channel) | rows |
  |---|---|---|
  | `contourdata-beta_0-Lp-MC3.07.22.24_SNaCl.csv` | `beta_0,Lp,log10_obj_cr` | 324 |
  | `contourdata-sigma-beta_1-MC3.07.22.24_SNaCl.csv` | `sigma,beta_1,log10_obj_cr` | 324 |

  Single-channel `log10_obj_cr` (retentate-conductivity log objective). Co-located `beta_contours-MC3.07.22.24_SNaCl.{json,png}`. Only run-id present: `MC3.07.22.24_SNaCl` (its `MC3.07.22.24_SNaCl/` subdir is empty).

### 4.6 Band-masked σ×Lp — `band_value_test/`
Path glob: `nf270/band_value_test/band_masked_contourdata-<RUN_ID>-band{first}_{last}.csv` (per band) — schema **`sigma,Lp,Obj_retentate_concentration`** (single-channel, log10 obj_cr). Co-located `band_masked_contour-<RUN_ID>.png`, `band_masked_contour-<RUN_ID>.json`.

### 4.7 `matlab_ports*` — MATLAB-replica sweeps
3-channel CSVs = 5 cols (`…,Obj_mass,Obj_concentration,Obj_retentate_concentration`); 5-channel CSVs = 7 cols (adds `Obj_permeate_conductivity,Obj_retentate_conductivity`). Co-located `objcontour-*.png`.

| variant (path glob under `nf270/`) | run-ids | CSV filename(s) | channels | grid (rows) |
|---|---|---|---|---|
| `matlab_ports/<run_id>/` | 1: `MC2.05.07.24_NaCl` | `contourdata-x_B-y_Lp.csv`, `contourdata-x_sigma-y_Lp.csv` | 3-ch (5 col) | 25 (5×5); +`sigma_sensitivity.png` |
| `matlab_ports_5ch/<run_id>/` | 3: CaCl2, NaCl, LaCl3 | `contourdata-x_sigma-y_Lp.csv` | 5-ch (7 col) | 400 (20×20) |
| `matlab_ports_5ch_allpairs/<run_id>/` | 1: `MC2.05.07.24_NaCl` | `-x_B-y_Lp`, `-x_B-y_sigma`, `-x_sigma-y_B`, `-x_sigma-y_Lp` | 5-ch (7 col) | 400 (20×20) |
| `matlab_ports_5ch_full/<run_id>/` | **11 (all)** + `_index.json` | `contourdata-x_sigma-y_Lp.csv` (+`_meta.json`) | 5-ch (7 col) | 225 (15×15) |

**5-channel schema** (e.g. `matlab_ports_5ch_full`): `sigma,Lp,Obj_mass,Obj_concentration,Obj_retentate_concentration,Obj_permeate_conductivity,Obj_retentate_conductivity`.
**Sweep metadata (`matlab_ports_5ch_full`):** root `_index.json` (`sweep_type:"5channel_sigma_x_Lp"`, `grid:15`, `n_workers:6`, `nfe:80`, `sheets[]`×11 with `warm_start{Lp,B,sigma,S0,S}, Lp_range, sigma_range, B_fixed, argmin_per_channel{…}`); per-run `_meta.json` mirrors per-sheet fields. **B held fixed (`B_fixed`)** while sweeping σ×Lp. No explicit cF weighting field (conductivity-channel sweep; weighting set in driver).

---

## 5. Fit caches — where θ* lives

### 5.1 Warm-start converged fits — `warm_start_fits/summary.json`
`nf270/warm_start_fits/summary.json` — JSON array of 11 sheets. Each carries `warm_start_theta_init`, `baseline_theta`, and the converged `warm_start_theta_fit` (the θ*), decomposing `obj_total` into `warm_start_obj_m/obj_cv/obj_cr`. Read by `_run_3d_contour_concentrating.py` and `_run_all_pairs_5channel.py` (key `warm_start_theta_fit`, matched by `run_id`). **These values predate the current weighting; several σ pin to 0/1 bounds.**

| run_id | Lp | B | sigma | obj_total |
|---|---|---|---|---|
| MC2.05.07.24_NaCl | 9.95339 | 9.73225 | 1.0 | 4357.27 |
| MC3.07.22.24_SNaCl | 8.21586 | 14.67912 | 0.45058 | 555.75 |
| MC4.07.11.24_SNaCl | 7.78596 | 2.02359 | 0.0 | 9536.22 |
| MC5.07.23.24_NaCl | 10.02301 | 3.69828 | 1.0 | 5793.63 |
| MC5.07.23.24_SNaCl | 9.44521 | 3.06517 | 1.0 | 4763.59 |
| MC5.07.23.24_S2NaCl | 7.38062 | 10.14469 | 0.0 | 829.62 |
| MC2.05.07.24_CaCl2 | 4.00043 | 0.50937 | 1.0 | 6717.14 |
| MC3.07.11.24_SCaCl2 | 5.53322 | 9.10985 | 1.0 | 4079.67 |
| MC3.07.12.24_S2CaCl2 | 4.63008 | 30.0 (bound) | 1.0 | 9251.67 |
| MC2.05.21.24_LaCl3 | 2.49193 | 0.31335 | 1.0 | 17711.17 |
| MC4.07.11.24_SLaCl3 | 2.77582 | 0.17151 | 1.0 | 21345.14 |

### 5.2 Per-run B-form fits — `bform_study/<run_id>/result_<form>.json`
Glob: `nf270/bform_study/<run_id>/result_<form>.json` — **all 11 run-ids × all 6 forms** (`single,sat,donnan,poly1,poly2,poly3`) = **66 files**. Sibling dirs (`contours/, model_selection/, peclet/, pooling/, predictions/, profile_contours/, stack_contours/, Bsigma_fixedLp[_nofloor]/, Lp_identifiability_panels/`) hold **no** `result_*.json`.
**Schema:** `parameters{Lp,B,sigma,S0,S}, WSSE3, Obj, obj_m, obj_cv, obj_cr, n_Bparams, n_data, AIC, AICc, BIC, multistart, FIM{std,eig_min,eig_max,cond,det}, run_id, salt, stage, form`. These are the anchor fits read by `_run_Bsigma_fixedLp.py`, `_run_profile_contours.py`, and `_run_bform_contours.py` (`result_single.json` → `parameters`).

`single`-form θ* + per-run best-AIC form:

| run_id | salt | Lp | B | sigma | best-AIC form |
|---|---|---|---|---|---|
| MC2.05.07.24_CaCl2 | CaCl2 | 7.539 | 5.711 | 1.00 | sat |
| MC2.05.07.24_NaCl | NaCl | 9.913 | 9.997 | 1.00 | poly3 |
| MC2.05.21.24_LaCl3 | LaCl3 | 3.536 | 0.316 | 0.00 | poly3 |
| MC3.07.11.24_SCaCl2 | CaCl2 | 6.778 | 8.722 | 1.00 | poly3 |
| MC3.07.12.24_S2CaCl2 | CaCl2 | 5.525 | 14.034 | 1.00 | poly3 |
| MC3.07.22.24_SNaCl | NaCl | 8.984 | 13.876 | 1.00 | poly3 |
| MC4.07.11.24_SLaCl3 | LaCl3 | 0.500 | 0.133 | 1.00 | poly2 |
| MC4.07.11.24_SNaCl | NaCl | 11.003 | 5.675 | 0.72 | poly3 |
| MC5.07.23.24_NaCl | NaCl | 9.312 | 10.138 | 1.00 | poly3 |
| MC5.07.23.24_S2NaCl | NaCl | 7.881 | 11.722 | 0.38 | poly3 |
| MC5.07.23.24_SNaCl | NaCl | 8.863 | 8.195 | 1.00 | poly3 |

### 5.3 Other anchor/centering fits
- `nf270/contour_panels/<run_id>/centering_fit.json` and `nf270/bform_study/contours/<run_id>/centering_fit.json` (when cold-fit) — `theta_fit{Lp,B,sigma,S0,S}` + `obj_m/obj_cv/obj_cr`.
- `nf270/beta_contours/<RUN_ID>/centering_fit.json` (β-form centering: `single` fit mapped to `(β₀=B, β₁=0)`).

---

## 6. Weighting status callout

| Convention | What it means | Where it lives |
|---|---|---|
| **2% / no-floor (CANONICAL)** | `CF_FRAC=0.02`, `CF_FLOOR=None` (`cf_residual_scale_fraction=0.02`, `cf_residual_floor_mM=null`) — the Lilonfe et al. weighting; library defaults at lines 271/299/300 | `contour3d/<run_id>/grid3d.csv` (driver sets it per worker; only `MC2.05.07.24_NaCl/_meta.json` records it); `bform_study/Bsigma_fixedLp_nofloor/…` (dir name + `BSIGMA_VARIANT=2pct_nofloor` default + 2026-06-18 timestamps) |
| **2% + 1 mM floor** | `CF_FRAC=0.02` with `NF270_CF_RESIDUAL_FLOOR_MM=1.0` set in-driver | `bform_study/Bsigma_fixedLp/…` (`BSIGMA_VARIANT=2pct_floor1mm`); `bform_study/profile_contours/…` & `stack_contours/…`; `beta_contours/…` (both drivers); `band_value_test/…`; `bform_study/contours/…`. Archived 1 mM-floor panel set: `nf270/contour_panels_cF_floor_1mM.zip` |
| **Legacy 0.3% (STALE)** | `cf_residual_scale_fraction` absent from `_meta.json` → implicitly the older 0.3% convention (gave WSSE ~44× too tight on cF); these grids predate the reweighting | All `contour3d/<run_id>/_meta.json` **except** `MC2.05.07.24_NaCl`; `warm_start_fits/summary.json` θ* (predate current weighting); the non-`_nofloor` `contour_panels/` set (live CSVs; weighting only inferable from sibling `centering_fit.json`) |
| **Unrecorded / conductivity-weighted** | cF weighting not in JSON; conductivity-channel sweep, weights set in driver | all `matlab_ports*` (`_meta.json`/`_index.json` carry `B_fixed`+`argmin_per_channel`, no `cf_residual_scale_fraction`) |

**Single most authoritative weighting record in the whole tree:** `nf270/contour3d/MC2.05.07.24_NaCl/_meta.json` → `cf_residual_scale_fraction=0.02`, `cf_residual_floor_mM=null`, `weighting="retentate cF = 2% relative, NO floor (Lilonfe et al.); permeate cV 3%; mass 0.01 g"`.

**Quick rule of thumb:**
- Want canonical (pure 2%, no floor): use `contour3d/grid3d.csv` and `bform_study/Bsigma_fixedLp_nofloor/`.
- Floored (2%+1 mM): everything under `profile_contours/`, `stack_contours/`, `Bsigma_fixedLp/`, `beta_contours/`, `band_value_test/`, `bform_study/contours/`.
- Treat `warm_start_fits/summary.json` θ* and any `contour3d/_meta.json` without `cf_residual_scale_fraction` as **stale 0.3%-era** provenance.
