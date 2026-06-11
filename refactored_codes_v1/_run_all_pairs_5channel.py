#!/usr/bin/env python3
"""Run all three 2D parameter-pair sweeps (σ × Lp, B × Lp, B × σ) for a single
sheet, each with the 5-channel kernel (mass + 2 conc + 2 cond).

Output per sweep:
    paper_artifacts/nf270/matlab_ports_5ch_allpairs/{run_id}/
        ├── contourdata-x_sigma-y_Lp.csv         (σ × Lp at fixed B)
        ├── contourdata-x_B-y_Lp.csv             (B × Lp at fixed σ)
        ├── contourdata-x_B-y_sigma.csv          (B × σ at fixed Lp)
        ├── objcontour-x_sigma-y_Lp.png          (5-panel 1×5)
        ├── objcontour-x_B-y_Lp.png              (5-panel 1×5)
        ├── objcontour-x_B-y_sigma.png           (5-panel 1×5)
        └── _meta.json

Default: MC2.05.07.24_NaCl, 20×20 grid per sweep, 4 workers.

Usage:
  python3 refactored_codes_v1/_run_all_pairs_5channel.py                       # MC2 NaCl, 20×20
  python3 refactored_codes_v1/_run_all_pairs_5channel.py --sheet MC3.07.22.24_SNaCl
  python3 refactored_codes_v1/_run_all_pairs_5channel.py --grid 30 --workers 6
"""

import os
import sys
import json
import time
import argparse
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import pandas as pd

REPO_ROOT = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo")
sys.path.insert(0, str(REPO_ROOT / "refactored_codes_v1"))
os.chdir(REPO_ROOT)


SHEET_DEFAULT = "MC2.05.07.24_NaCl"
GRID = 20
NFE = 80
N_WORKERS = 4

# Order: (x_var, y_var) — y_var is the Y axis; third parameter held at warm-start.
# AXIS CONVENTION (rule-of-thumb, 2026-06-11): Y-axis priority L_p > B > sigma.
SWEEPS = [
    ("sigma", "Lp"),     # Lp on Y  (L_p > sigma)
    ("B",     "Lp"),     # Lp on Y  (L_p > B)
    ("sigma", "B"),      # B  on Y  (B > sigma)  — was ("B","sigma") (sigma on Y); flipped per rule
]


_DATA_STRU = None
_BASE_THETA = None


def _load_warm_start_theta(run_id):
    p = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
         / "paper_artifacts" / "nf270" / "warm_start_fits" / "summary.json")
    for e in json.load(open(p)):
        if e.get("run_id") == run_id:
            return dict(e["warm_start_theta_fit"])
    raise KeyError(f"{run_id} not in warm_start_fits/summary.json")


def _worker_init(run_id, base_theta):
    global _DATA_STRU, _BASE_THETA
    import refactored_ucb_library as _lib
    fam = _lib.NF270_RUN_REGISTRY[run_id]
    loaded = _lib.loadxlsx(_lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])
    _DATA_STRU = loaded["data_stru"]
    _BASE_THETA = dict(base_theta)


def _eval_cell(args):
    import refactored_ucb_library as _lib
    x_var, x_val, y_var, y_val, nfe = args
    theta = dict(_BASE_THETA)
    theta[x_var] = float(x_val)
    theta[y_var] = float(y_val)

    def _l10(v):
        if v is None or not np.isfinite(v) or v <= 0:
            return float("nan")
        return float(np.log10(v))

    try:
        ind = _lib.calc_ind_objectives_5channel_py(
            theta, _DATA_STRU,
            mode="Lag", B_form="single",
            workflow_family="DATA3", nfe=nfe,
        )
        return (x_val, y_val,
                _l10(ind.m), _l10(ind.cp), _l10(ind.cr),
                _l10(ind.cp_cond), _l10(ind.cr_cond),
                None)
    except Exception as e:
        return (x_val, y_val,
                float("nan"), float("nan"), float("nan"),
                float("nan"), float("nan"),
                f"{type(e).__name__}: {e}")


def _bounds_for(var, theta, salt):
    import refactored_ucb_library as _lib
    if var == "Lp":
        lp = float(theta["Lp"])
        return (0.1 * lp, 2.0 * lp)
    if var == "sigma":
        return (0.0, 1.0)
    if var == "B":
        return _lib.NF270_B_BOUNDS_PER_SALT.get(salt, _lib.NF270_B_BOUNDS_DEFAULT)
    raise ValueError(var)


def _salt_from_run_id(run_id):
    if "NaCl" in run_id:  return "NaCl"
    if "CaCl2" in run_id: return "CaCl2"
    if "LaCl3" in run_id: return "LaCl3"
    return None


def run_sweep(run_id, x_var, y_var, base_theta, save_dir, *,
              grid, nfe, n_workers, pool):
    """Run one 2D sweep, write CSV, render PNG. Reuses an existing Pool."""
    import refactored_ucb_library as _lib
    salt = _salt_from_run_id(run_id)
    x_lb, x_ub = _bounds_for(x_var, base_theta, salt)
    y_lb, y_ub = _bounds_for(y_var, base_theta, salt)

    fixed = {k: base_theta[k] for k in ("Lp", "B", "sigma") if k not in (x_var, y_var)}

    print(f"\n  ── sweep: {x_var} × {y_var} ──")
    print(f"     {x_var} range: [{x_lb:.3g}, {x_ub:.3g}]")
    print(f"     {y_var} range: [{y_lb:.3g}, {y_ub:.3g}]")
    print(f"     fixed: {fixed}")

    x_grid = np.linspace(x_lb, x_ub, grid)
    y_grid = np.linspace(y_lb, y_ub, grid)
    cells = [(x_var, x, y_var, y, nfe) for x in x_grid for y in y_grid]
    n_cells = len(cells)

    t0 = time.time()
    results = []
    chunksize = max(1, n_cells // (n_workers * 10))
    for i, res in enumerate(pool.imap_unordered(_eval_cell, cells, chunksize=chunksize), 1):
        results.append(res)
        if i % max(1, n_cells // 5) == 0 or i == n_cells:
            dt = time.time() - t0
            rate = i / max(dt, 1e-9)
            eta = (n_cells - i) / max(rate, 1e-9)
            print(f"     [{i:>4d}/{n_cells}]  ({100*i/n_cells:5.1f}%)  "
                  f"rate={rate:.2f}/s  eta={eta/60:.1f}min")
    dt = time.time() - t0
    print(f"     elapsed: {dt/60:.2f} min")

    df = pd.DataFrame(results, columns=[
        x_var, y_var,
        "Obj_mass", "Obj_concentration", "Obj_retentate_concentration",
        "Obj_permeate_conductivity", "Obj_retentate_conductivity",
        "_err",
    ])
    n_failed = df["_err"].notna().sum()
    print(f"     failures: {n_failed} / {n_cells}")

    csv_path = save_dir / f"contourdata-x_{x_var}-y_{y_var}.csv"
    df.drop(columns=["_err"]).to_csv(csv_path, index=False)

    try:
        png = _lib.plot_contour_py(
            df.drop(columns=["_err"]),
            x_var, y_var,
            save_dir=save_dir,
            title=f"{run_id}  ·  {x_var} × {y_var}  ·  5-channel  "
                  f"(fixed: {', '.join(f'{k}={v:.3g}' for k,v in fixed.items())})",
        )
        print(f"     rendered: {png.relative_to(REPO_ROOT)}")
    except Exception as e:
        print(f"     RENDER FAILED: {type(e).__name__}: {e}")

    return {
        "x_var": x_var, "y_var": y_var,
        "x_range": [float(x_lb), float(x_ub)],
        "y_range": [float(y_lb), float(y_ub)],
        "fixed": fixed,
        "elapsed_s": dt,
        "n_cells": n_cells,
        "n_failed": int(n_failed),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet",   default=SHEET_DEFAULT)
    ap.add_argument("--grid",    type=int, default=GRID)
    ap.add_argument("--workers", type=int, default=N_WORKERS)
    ap.add_argument("--nfe",     type=int, default=NFE)
    args = ap.parse_args()

    run_id = args.sheet
    base_theta = _load_warm_start_theta(run_id)
    save_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                / "paper_artifacts" / "nf270" / "matlab_ports_5ch_allpairs" / run_id)
    save_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== {run_id} — all 3 pairs × 5 channels ===")
    print(f"  warm-start θ: Lp={base_theta['Lp']:.3f}  B={base_theta['B']:.3f}  "
          f"sigma={base_theta['sigma']:.3f}  S0={base_theta.get('S0', 'NA')}  "
          f"S={base_theta.get('S', float('nan')):.4f}")
    print(f"  grid: {args.grid}×{args.grid} per sweep, 3 sweeps  →  "
          f"{3 * args.grid * args.grid} total cells")
    print(f"  workers: {args.workers}    nfe: {args.nfe}")
    print(f"  save_dir: {save_dir.relative_to(REPO_ROOT)}")

    overall_t0 = time.time()
    sweep_meta = []
    with Pool(args.workers, initializer=_worker_init,
              initargs=(run_id, base_theta)) as pool:
        for x_var, y_var in SWEEPS:
            try:
                meta = run_sweep(run_id, x_var, y_var, base_theta, save_dir,
                                 grid=args.grid, nfe=args.nfe,
                                 n_workers=args.workers, pool=pool)
                sweep_meta.append(meta)
            except Exception as e:
                print(f"  *** sweep ({x_var}, {y_var}) FAILED: "
                      f"{type(e).__name__}: {e}")

    overall_dt = time.time() - overall_t0
    print(f"\n{'='*60}")
    print(f"DONE.  Total: {overall_dt/60:.2f} min")

    meta = {
        "run_id":       run_id,
        "grid":         args.grid,
        "n_workers":    args.workers,
        "nfe":          args.nfe,
        "warm_start":   base_theta,
        "sweeps":       sweep_meta,
        "total_min":    overall_dt / 60.0,
    }
    with open(save_dir / "_meta.json", "w") as f:
        json.dump(meta, f, indent=2, default=float)
    print(f"meta:  {(save_dir / '_meta.json').relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
