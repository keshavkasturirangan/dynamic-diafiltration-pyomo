#!/usr/bin/env python3
"""Run the 5-channel MATLAB-port contour sweep (σ × Lp) on the 3 concentrating-regime sheets.

Grid:
  σ  ∈ linspace(0.0, 1.0, 20)
  Lp ∈ [0.1·Lp_warm, 2·Lp_warm]   linspace, 20 points
  Fixed: B = warm-start B, S0/S = warm-start values

Output per sheet:
  paper_artifacts/nf270/matlab_ports_5ch/{run_id}/
    ├── contourdata-x_sigma-y_Lp.csv          ← 7-col tidy long-form (5 channel SSRs)
    ├── objcontour-x_sigma-y_Lp.png            ← 1×5 panel composite (figure_s5 style)
    ├── _meta.json                             ← provenance
    └── _failures.log                          ← failed grid cells

Per-cell wall time ≈ 3 s on 8 workers; 400 cells/sheet × 3 sheets = ~7-10 min total.

Usage:
  python3 refactored_codes_v1/_run_5channel_demo.py                          # full run
  python3 refactored_codes_v1/_run_5channel_demo.py --grid 10                # half resolution
  python3 refactored_codes_v1/_run_5channel_demo.py --sheets MC2.05.07.24_NaCl
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


SHEETS_DEFAULT = [
    "MC2.05.07.24_NaCl",
    "MC2.05.07.24_CaCl2",
    "MC2.05.21.24_LaCl3",
]
GRID = 20
NFE = 80
N_WORKERS = 8


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
    sigma, Lp, nfe = args
    theta = dict(_BASE_THETA)
    theta["Lp"]    = float(Lp)
    theta["sigma"] = float(sigma)
    # B kept at warm-start value (so this is a σ × Lp sweep at fixed B)

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
        return (sigma, Lp,
                _l10(ind.m), _l10(ind.cp), _l10(ind.cr),
                _l10(ind.cp_cond), _l10(ind.cr_cond),
                None)
    except Exception as e:
        return (sigma, Lp,
                float("nan"), float("nan"), float("nan"),
                float("nan"), float("nan"),
                f"{type(e).__name__}: {e}")


def sweep_sheet(run_id, *, grid_density, nfe, n_workers):
    import refactored_ucb_library as _lib

    base_theta = _load_warm_start_theta(run_id)
    print(f"\n=== {run_id} ===")
    print(f"  warm-start theta: Lp={base_theta['Lp']:.3f}  "
          f"B={base_theta['B']:.3f}  sigma={base_theta['sigma']:.3f}")

    Lp_warm = float(base_theta["Lp"])
    Lp_grid = np.linspace(0.1 * Lp_warm, 2.0 * Lp_warm, grid_density)
    sigma_grid = np.linspace(0.0, 1.0, grid_density)

    cells = [(sg, lp, nfe) for sg in sigma_grid for lp in Lp_grid]
    n_cells = len(cells)

    save_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                / "paper_artifacts" / "nf270" / "matlab_ports_5ch" / run_id)
    save_dir.mkdir(parents=True, exist_ok=True)

    print(f"  grid:  sigma={grid_density}  Lp={grid_density}  →  {n_cells} cells")
    print(f"  Lp range: [{Lp_grid[0]:.2f}, {Lp_grid[-1]:.2f}]")
    print(f"  workers: {n_workers}  nfe: {nfe}")
    print(f"  save_dir: {save_dir.relative_to(REPO_ROOT)}")

    t0 = time.time()
    with Pool(n_workers, initializer=_worker_init,
              initargs=(run_id, base_theta)) as pool:
        results = []
        chunksize = max(1, n_cells // (n_workers * 10))
        for i, res in enumerate(pool.imap_unordered(_eval_cell, cells, chunksize=chunksize), 1):
            results.append(res)
            if i % max(1, n_cells // 10) == 0 or i == n_cells:
                dt = time.time() - t0
                rate = i / max(dt, 1e-9)
                eta = (n_cells - i) / max(rate, 1e-9)
                print(f"  [{i:>4d}/{n_cells}]  ({100*i/n_cells:5.1f}%)  "
                      f"rate={rate:.2f}/s  eta={eta/60:.1f}min")
    dt = time.time() - t0
    print(f"  elapsed: {dt/60:.2f} min  ({dt/n_cells*1000:.0f} ms/cell)")

    df = pd.DataFrame(results, columns=[
        "sigma", "Lp",
        "Obj_mass", "Obj_concentration", "Obj_retentate_concentration",
        "Obj_permeate_conductivity", "Obj_retentate_conductivity",
        "_err",
    ])
    failed = df[df["_err"].notna()]
    print(f"  failures: {len(failed)} / {n_cells}")

    csv_path = save_dir / "contourdata-x_sigma-y_Lp.csv"
    df.drop(columns=["_err"]).to_csv(csv_path, index=False)

    if len(failed) > 0:
        with open(save_dir / "_failures.log", "w") as f:
            for _, r in failed.iterrows():
                f.write(f"sigma={r.sigma:.4f} Lp={r.Lp:.4f} :: {r._err}\n")

    meta = {
        "run_id":        run_id,
        "grid_density":  grid_density,
        "n_cells":       n_cells,
        "n_failures":    int(len(failed)),
        "wall_time_min": dt / 60.0,
        "Lp_range":      [float(Lp_grid[0]), float(Lp_grid[-1])],
        "sigma_range":   [0.0, 1.0],
        "fixed":         {"B": base_theta["B"], "S0": base_theta.get("S0"),
                          "S": base_theta.get("S")},
        "warm_start":    base_theta,
    }
    with open(save_dir / "_meta.json", "w") as f:
        json.dump(meta, f, indent=2, default=float)

    # Render 1×5 panel
    try:
        png = _lib.plot_contour_py(
            df.drop(columns=["_err"]),
            "sigma", "Lp",
            save_dir=save_dir,
            title=f"{run_id}  ·  σ × Lp  ·  5-channel (mass + 2 conc + 2 cond)",
        )
        print(f"  rendered: {png.relative_to(REPO_ROOT)}")
    except Exception as e:
        print(f"  RENDER FAILED: {type(e).__name__}: {e}")

    return csv_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheets", nargs="*", default=None)
    ap.add_argument("--grid",   type=int, default=GRID)
    ap.add_argument("--workers", type=int, default=N_WORKERS)
    ap.add_argument("--nfe",     type=int, default=NFE)
    args = ap.parse_args()

    sheets = args.sheets or SHEETS_DEFAULT
    print(f"5-CHANNEL CONTOUR SWEEP")
    print(f"  sheets:  {sheets}")
    print(f"  grid:    {args.grid}×{args.grid} σ × Lp")
    print(f"  workers: {args.workers}")
    print(f"  nfe:     {args.nfe}")

    overall_t0 = time.time()
    for run_id in sheets:
        try:
            sweep_sheet(run_id,
                        grid_density=args.grid,
                        nfe=args.nfe,
                        n_workers=args.workers)
        except Exception as e:
            print(f"\n  *** {run_id} FAILED: {type(e).__name__}: {e}")

    dt = time.time() - overall_t0
    print(f"\n{'='*60}")
    print(f"DONE.  Total: {dt/60:.2f} min")


if __name__ == "__main__":
    main()
