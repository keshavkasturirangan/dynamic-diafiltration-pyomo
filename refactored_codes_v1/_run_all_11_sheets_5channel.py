#!/usr/bin/env python3
"""5-channel σ × Lp sweep for ALL 11 NF270 single-salt sheets in DATA3.

This is the comprehensive scan that lets the advisor deck show every sheet
classified into one of three stages of identifiability:

    Stage 1 — well-identified         (4-5 channels argmin near warm-start)
    Stage 2 — inconsistent / minimized (channels argmin at different points)
    Stage 3 — not estimable            (no consistent basin)

Output per sheet — paper_artifacts/nf270/matlab_ports_5ch_full/{run_id}/
    ├── contourdata-x_sigma-y_Lp.csv     (400 rows × 7 cols)
    ├── objcontour-x_sigma-y_Lp.png      (1×5 panels with optimum text boxes)
    ├── _meta.json                       (provenance + per-channel argmin)
    └── _failures.log                    (if any)

Cost (4 workers / sheet, sheets sequential):
    ~5-10 min per sheet × 11 = ~75-110 min total.

Existing sheets that already have data in matlab_ports_5ch/ (MC2.05.07.24_NaCl,
MC2.05.07.24_CaCl2, MC2.05.21.24_LaCl3) are RE-RUN here so all 11 results live
in one consistent directory tree with the new text-box-styled PNGs.
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


# All 11 NF270 single-salt sheets — grouped by salt for log readability
SHEETS_ALL = [
    # NaCl (6 sheets)
    "MC2.05.07.24_NaCl",
    "MC3.07.22.24_SNaCl",
    "MC4.07.11.24_SNaCl",
    "MC5.07.23.24_NaCl",
    "MC5.07.23.24_SNaCl",
    "MC5.07.23.24_S2NaCl",
    # CaCl2 (3 sheets)
    "MC2.05.07.24_CaCl2",
    "MC3.07.11.24_SCaCl2",
    "MC3.07.12.24_S2CaCl2",
    # LaCl3 (2 sheets)
    "MC2.05.21.24_LaCl3",
    "MC4.07.11.24_SLaCl3",
]

GRID    = 20
NFE     = 80
WORKERS = 4

_DATA_STRU  = None
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
    loaded = _lib.loadxlsx(_lib.NF270_DEFAULT_ROOT / fam["workbook"],
                            sheet=fam["sheet"])
    _DATA_STRU = loaded["data_stru"]
    _BASE_THETA = dict(base_theta)


def _eval_cell(args):
    import refactored_ucb_library as _lib
    sigma, Lp, nfe = args
    theta = dict(_BASE_THETA)
    theta["Lp"]    = float(Lp)
    theta["sigma"] = float(sigma)

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
        return (sigma, Lp, *([float("nan")] * 5),
                f"{type(e).__name__}: {e}")


def sweep_sheet(run_id, *, grid, nfe, n_workers):
    import refactored_ucb_library as _lib

    base = _load_warm_start_theta(run_id)
    Lp_warm = float(base["Lp"])
    Lp_grid = np.linspace(0.1 * Lp_warm, 2.0 * Lp_warm, grid)
    sigma_grid = np.linspace(0.0, 1.0, grid)
    cells = [(sg, lp, nfe) for sg in sigma_grid for lp in Lp_grid]
    n_cells = len(cells)

    save_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                / "paper_artifacts" / "nf270" / "matlab_ports_5ch_full" / run_id)
    save_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== {run_id} ===")
    print(f"  warm-start: Lp={base['Lp']:.3f}  B={base['B']:.3f}  "
          f"sigma={base['sigma']:.3f}    B fixed during sweep.")
    print(f"  grid: {grid}×{grid} = {n_cells} cells   workers: {n_workers}")

    t0 = time.time()
    with Pool(n_workers, initializer=_worker_init,
              initargs=(run_id, base)) as pool:
        results = []
        for i, res in enumerate(
                pool.imap_unordered(_eval_cell, cells,
                                    chunksize=max(1, n_cells // (n_workers * 8))), 1):
            results.append(res)
            if i % max(1, n_cells // 5) == 0 or i == n_cells:
                dt = time.time() - t0
                rate = i / max(dt, 1e-9)
                eta = (n_cells - i) / max(rate, 1e-9)
                print(f"  [{i:>4d}/{n_cells}]  ({100*i/n_cells:5.1f}%)  "
                      f"rate={rate:.2f}/s  eta={eta/60:.1f}min")
    dt = time.time() - t0
    print(f"  elapsed: {dt/60:.2f} min")

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

    # Per-channel argmin summary
    cols = ["Obj_mass", "Obj_concentration", "Obj_retentate_concentration",
            "Obj_permeate_conductivity", "Obj_retentate_conductivity"]
    short = {"Obj_mass":"mass","Obj_concentration":"perm_conc",
             "Obj_retentate_concentration":"ret_conc",
             "Obj_permeate_conductivity":"perm_cond",
             "Obj_retentate_conductivity":"ret_cond"}
    argmins = {}
    for c in cols:
        s = df[c]
        if s.notna().sum() == 0:
            argmins[short[c]] = None
            continue
        i = s.idxmin()
        argmins[short[c]] = {
            "sigma":  float(df.at[i, "sigma"]),
            "Lp":     float(df.at[i, "Lp"]),
            "log10":  float(s.iloc[i]),
        }

    meta = {
        "run_id":          run_id,
        "grid":            grid,
        "n_cells":         n_cells,
        "n_failures":      int(len(failed)),
        "wall_time_min":   dt / 60.0,
        "warm_start":      base,
        "Lp_range":        [float(Lp_grid[0]), float(Lp_grid[-1])],
        "sigma_range":     [0.0, 1.0],
        "B_fixed":         float(base["B"]),
        "argmin_per_channel": argmins,
    }
    with open(save_dir / "_meta.json", "w") as f:
        json.dump(meta, f, indent=2, default=float)

    # Render with text-box-style optima (uses the updated plot_contour_py)
    try:
        png = _lib.plot_contour_py(
            df.drop(columns=["_err"]),
            "sigma", "Lp",
            save_dir=save_dir,
            title=f"{run_id}  ·  σ × Lp  ·  5-channel  "
                  f"(B fixed at warm-start = {base['B']:.2f})",
        )
        print(f"  rendered: {png.relative_to(REPO_ROOT)}")
    except Exception as e:
        print(f"  RENDER FAILED: {type(e).__name__}: {e}")

    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheets",  nargs="*", default=None)
    ap.add_argument("--grid",    type=int, default=GRID)
    ap.add_argument("--workers", type=int, default=WORKERS)
    ap.add_argument("--nfe",     type=int, default=NFE)
    args = ap.parse_args()

    sheets = args.sheets or SHEETS_ALL
    print(f"5-CHANNEL σ × Lp SWEEP — ALL NF270 SINGLE-SALT SHEETS")
    print(f"  sheets:  {len(sheets)}  ({sheets})")
    print(f"  grid:    {args.grid}×{args.grid}")
    print(f"  workers: {args.workers}    nfe: {args.nfe}")

    overall_t0 = time.time()
    metas = []
    for run_id in sheets:
        try:
            meta = sweep_sheet(run_id, grid=args.grid, nfe=args.nfe,
                                n_workers=args.workers)
            metas.append(meta)
        except Exception as e:
            print(f"\n  *** {run_id} FAILED: {type(e).__name__}: {e}")
            metas.append({"run_id": run_id, "error": f"{type(e).__name__}: {e}"})

    overall_dt = time.time() - overall_t0

    # Master index file at the campaign root
    index_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                 / "paper_artifacts" / "nf270" / "matlab_ports_5ch_full")
    index_dir.mkdir(parents=True, exist_ok=True)
    with open(index_dir / "_index.json", "w") as f:
        json.dump({
            "sweep_type":       "5channel_sigma_x_Lp",
            "grid":             args.grid,
            "n_workers":        args.workers,
            "nfe":              args.nfe,
            "total_min":        overall_dt / 60.0,
            "sheets":           metas,
        }, f, indent=2, default=float)

    print(f"\n{'='*68}")
    print(f"DONE.  Total: {overall_dt/60:.2f} min  for {len(sheets)} sheets")
    print(f"Index: {(index_dir / '_index.json').relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
