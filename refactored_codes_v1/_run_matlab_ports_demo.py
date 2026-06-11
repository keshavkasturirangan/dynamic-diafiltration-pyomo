#!/usr/bin/env python3
"""Demo driver for the four MATLAB→Python ports against a DATA3 NF270 sheet.

For each sheet (default: the 3 concentrating-regime sheets used by
_run_3d_contour_concentrating.py), runs the full MATLAB-equivalent pipeline:

    1. Load data_stru + warm-start θ.
    2. calc_contour_2d_py(σ, Lp)   →  contourdata-x_sigma-y_Lp.csv
       plot_contour_py(...)        →  objcontour-x_sigma-y_Lp.png
    3. calc_contour_2d_py(B, Lp)   →  contourdata-x_B-y_Lp.csv
       plot_contour_py(...)        →  objcontour-x_B-y_Lp.png
    4. sigma_sensitivity_py(...)   →  sigma_sensitivity.png

Output lands in:
    paper_artifacts/nf270/matlab_ports/{run_id}/
        contourdata-x_sigma-y_Lp.csv      ← matches MATLAB filename
        contourdata-x_B-y_Lp.csv          ← matches MATLAB filename
        objcontour-x_sigma-y_Lp.png       ← matches MATLAB filename
        objcontour-x_B-y_Lp.png           ← matches MATLAB filename
        sigma_sensitivity.png             ← matches MATLAB filename

Usage:
    python3 refactored_codes_v1/_run_matlab_ports_demo.py                    # all 3 sheets
    python3 refactored_codes_v1/_run_matlab_ports_demo.py --grid 20          # faster
    python3 refactored_codes_v1/_run_matlab_ports_demo.py --sheets MC2.05.07.24_NaCl
    python3 refactored_codes_v1/_run_matlab_ports_demo.py --no-sigma-sens    # skip σ sweep
"""

import os
import sys
import json
import time
import argparse
from pathlib import Path

REPO_ROOT = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo")
sys.path.insert(0, str(REPO_ROOT / "refactored_codes_v1"))
os.chdir(REPO_ROOT)

import refactored_ucb_library as lib


SHEETS_DEFAULT = [
    "MC2.05.07.24_NaCl",
    "MC2.05.07.24_CaCl2",
    "MC2.05.21.24_LaCl3",
]


def _load_warm_start_theta(run_id):
    summary_path = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                    / "paper_artifacts" / "nf270" / "warm_start_fits"
                    / "summary.json")
    summ = json.load(open(summary_path))
    for entry in summ:
        if entry.get("run_id") == run_id:
            return dict(entry["warm_start_theta_fit"])
    raise KeyError(f"{run_id} not in {summary_path}")


def _load_sheet(run_id):
    family = lib.NF270_RUN_REGISTRY[run_id]
    workbook_path = lib.NF270_DEFAULT_ROOT / family["workbook"]
    sheet_name = family["sheet"]
    loaded = lib.loadxlsx(workbook_path, sheet=sheet_name)
    return loaded["data_stru"]


def run_one_sheet(run_id, *, grid_density, do_sigma_sens, n_sigma_points):
    """Run all four MATLAB-port ports against one sheet."""
    print(f"\n{'='*68}")
    print(f"=== {run_id} ===")

    base_theta = _load_warm_start_theta(run_id)
    data_stru  = _load_sheet(run_id)
    salt = str(data_stru.get("data_config", {}).get("namec", "?"))
    print(f"  salt: {salt}")
    print(f"  warm-start theta:  Lp={base_theta['Lp']:.3f}  "
          f"B={base_theta['B']:.3f}  sigma={base_theta['sigma']:.3f}  "
          f"S0={base_theta.get('S0', 'NA')}  S={base_theta.get('S', float('nan')):.4f}")

    save_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                / "paper_artifacts" / "nf270" / "matlab_ports" / run_id)
    save_dir.mkdir(parents=True, exist_ok=True)
    print(f"  save_dir: {save_dir.relative_to(REPO_ROOT)}")

    # ---- Sweep 1: (sigma, Lp) — MATLAB run_data_analysis.m line 75-86 ----
    print(f"\n  [1/4] calc_contour_2d_py(sigma, Lp)  grid={grid_density}x{grid_density}")
    t0 = time.time()
    df_sigLp = lib.calc_contour_2d_py(
        data_stru, base_theta,
        x_var="sigma", y_var="Lp",
        grid_density=grid_density, save_dir=save_dir,
    )
    print(f"        elapsed: {time.time()-t0:.1f} s   "
          f"n_finite/total: {df_sigLp['Obj_mass'].notna().sum()}/{len(df_sigLp)}")

    # ---- Sweep 2: (B, Lp) — MATLAB run_data_analysis.m line 88-98 ----
    print(f"\n  [2/4] calc_contour_2d_py(B, Lp)      grid={grid_density}x{grid_density}")
    t0 = time.time()
    df_BLp = lib.calc_contour_2d_py(
        data_stru, base_theta,
        x_var="B", y_var="Lp",
        grid_density=grid_density, save_dir=save_dir,
    )
    print(f"        elapsed: {time.time()-t0:.1f} s   "
          f"n_finite/total: {df_BLp['Obj_mass'].notna().sum()}/{len(df_BLp)}")

    # ---- Render: plot_contour_py for both sweeps ----
    print(f"\n  [3/4] plot_contour_py — sigma×Lp")
    p1 = lib.plot_contour_py(df_sigLp, "sigma", "Lp",
                              save_dir=save_dir, title=f"{run_id} — σ × Lp")
    print(f"        wrote: {p1.relative_to(REPO_ROOT)}")
    p2 = lib.plot_contour_py(df_BLp, "B", "Lp",
                              save_dir=save_dir, title=f"{run_id} — B × Lp")
    print(f"        wrote: {p2.relative_to(REPO_ROOT)}")

    # ---- Sigma sensitivity sweep ----
    sigsen_out = None
    if do_sigma_sens:
        import numpy as np
        print(f"\n  [4/4] sigma_sensitivity_py  n_sigma={n_sigma_points}")
        # MATLAB convention: 0, 0.5, 1.0 (3-point cycle), or finer
        sigma_values = np.linspace(0.0, 1.0, n_sigma_points)
        t0 = time.time()
        _, sigsen_out = lib.sigma_sensitivity_py(
            data_stru, sigma_values, base_theta,
            save_dir=save_dir, LOUD=True,
        )
        print(f"        elapsed: {time.time()-t0:.1f} s")
        if sigsen_out is not None:
            print(f"        wrote: {sigsen_out.relative_to(REPO_ROOT)}")
    else:
        print(f"\n  [4/4] sigma_sensitivity_py — SKIPPED (--no-sigma-sens)")

    return {
        "run_id":   run_id,
        "save_dir": save_dir,
        "outputs":  [save_dir / "contourdata-x_sigma-y_Lp.csv",
                     save_dir / "contourdata-x_B-y_Lp.csv",
                     p1, p2,
                     sigsen_out] if sigsen_out else [
                     save_dir / "contourdata-x_sigma-y_Lp.csv",
                     save_dir / "contourdata-x_B-y_Lp.csv",
                     p1, p2],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheets", nargs="*", default=None)
    ap.add_argument("--grid",   type=int, default=20,
                    help="Grid density for calc_contour_2d_py (50 = MATLAB default).")
    ap.add_argument("--no-sigma-sens", action="store_true",
                    help="Skip sigma_sensitivity_py (faster).")
    ap.add_argument("--n-sigma", type=int, default=3,
                    help="Number of sigma values for sigma_sensitivity_py.")
    args = ap.parse_args()

    sheets = args.sheets or SHEETS_DEFAULT
    print(f"\nSheets:  {sheets}")
    print(f"Grid:    {args.grid}×{args.grid} per sweep")
    print(f"σ-sens:  {'OFF' if args.no_sigma_sens else f'{args.n_sigma} points'}")

    overall_t0 = time.time()
    results = []
    for run_id in sheets:
        try:
            results.append(run_one_sheet(
                run_id,
                grid_density=args.grid,
                do_sigma_sens=not args.no_sigma_sens,
                n_sigma_points=args.n_sigma,
            ))
        except Exception as e:
            print(f"\n  *** {run_id} failed: {type(e).__name__}: {e}")

    overall_dt = time.time() - overall_t0
    print(f"\n{'='*68}")
    print(f"DONE.  Total: {overall_dt/60:.2f} min")
    for r in results:
        print(f"  {r['run_id']}: {len(r['outputs'])} artifacts → "
              f"{r['save_dir'].relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
