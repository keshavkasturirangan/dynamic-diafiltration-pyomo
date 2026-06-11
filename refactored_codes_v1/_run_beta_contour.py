#!/usr/bin/env python3
"""Drive a beta-vs-Lp/sigma objective-contour run for a DATA3 NF270 sheet.

Exercises the numeric-B_form (B_form>=1) branch of ``run_nf270_contour_for_sheet``:
the centering fit produces beta_0/beta_1 (B = beta_0 + beta_1*cIn) instead of a
lumped B, and ``sweep_pairs`` switches to the beta slices:

    (beta_0, Lp)      Lp on Y   (beta_0 plays the role of B)
    (sigma,  Lp)      Lp on Y
    (sigma,  beta_1)  beta_1 on Y

The single-shot B_form=1 centering fit from the default init is fragile on the
gold sheet, so we WARM-START it: first do a well-behaved B_form='single' fit,
map (Lp, B, sigma, S0, S) -> (Lp, beta_0=B, beta_1=0, sigma, S0, S), then refit
with B_form=1 from that seed.  That fitted theta centers the contour grid.

Usage:
    python3 _run_beta_contour.py [RUN_ID] [GRID] [NFE]
Defaults: RUN_ID=MC3.07.22.24_SNaCl, GRID=18, NFE=80.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib  # noqa: E402

RUN_ID = sys.argv[1] if len(sys.argv) > 1 else "MC3.07.22.24_SNaCl"
GRID = int(sys.argv[2]) if len(sys.argv) > 2 else 18
NFE = int(sys.argv[3]) if len(sys.argv) > 3 else 80

lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0   # match the warm-anchor campaign config

OUT = (HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts"
       / "nf270" / "beta_contours" / RUN_ID)
OUT.mkdir(parents=True, exist_ok=True)

fam = lib.NF270_RUN_REGISTRY[RUN_ID]
wb = lib.NF270_DEFAULT_ROOT / fam["workbook"]
ds = lib.loadxlsx(wb, sheet=fam["sheet"])["data_stru"]
mode = ds["mode"]

print(f"[beta-contour] run_id={RUN_ID} mode={mode} grid={GRID} nfe={NFE}")
print(f"[beta-contour] out={OUT}")

# --- Step 1: well-behaved lumped-B fit to seed the numeric fit ---------------
print("\n[beta-contour] step 1: B_form='single' seed fit ...")
fit_s, _, _ = lib.solve_model(ds, mode, sim_opt=False, B_form="single",
                              workflow_family="DATA3", nfe=NFE, LOUD=False,
                              solver_max_cpu_time=180)
ps = dict(fit_s["parameters"])
print(f"[beta-contour]   single fit: Lp={ps['Lp']:.4g} B={ps['B']:.4g} "
      f"sigma={ps['sigma']:.4g}  (WSSE_cr={fit_s.get('obj_cr')})")

seed = {"Lp": ps["Lp"], "sigma": ps["sigma"], "beta_0": ps["B"], "beta_1": 0.0}
for k in ("S0", "S"):
    if k in ps:
        seed[k] = ps[k]

# --- Step 2: numeric B_form=1 centering fit, warm-started from the seed -------
print("\n[beta-contour] step 2: B_form=1 centering fit (warm start) ...")
theta_fit = None
try:
    fit_1, _, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=False, B_form=1,
                                  workflow_family="DATA3", nfe=NFE, LOUD=False,
                                  solver_max_cpu_time=300)
    if isinstance(fit_1, dict) and "parameters" in fit_1:
        theta_fit = dict(fit_1["parameters"])
        print(f"[beta-contour]   B_form=1 fit: Lp={theta_fit.get('Lp'):.4g} "
              f"beta_0={theta_fit.get('beta_0'):.4g} beta_1={theta_fit.get('beta_1'):.4g} "
              f"sigma={theta_fit.get('sigma'):.4g}  (WSSE_cr={fit_1.get('obj_cr')})")
except Exception as exc:
    print(f"[beta-contour]   B_form=1 fit raised: {exc!r}")

if theta_fit is None:
    print("[beta-contour]   numeric fit failed; centering on the constructed seed.")
    theta_fit = dict(seed)

# --- Step 3: beta contour sweeps centered on theta_fit -----------------------
print("\n[beta-contour] step 3: beta contour sweeps ...")
artifact = lib.run_nf270_contour_for_sheet(
    RUN_ID, save_dir=OUT, grid_density=GRID, theta_fit=theta_fit,
    nfe=NFE, mode=mode, B_form=1,
)

print("\n[beta-contour] DONE")
print(f"[beta-contour] theta_fit center = {artifact.get('theta_fit')}")
for p in artifact.get("csv_paths", []):
    print(f"   CSV {p}")
for p in artifact.get("png_paths", []):
    print(f"   PNG {p}")
