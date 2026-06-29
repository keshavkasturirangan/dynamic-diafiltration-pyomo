#!/usr/bin/env python3
"""Run 30×30×10 3D objective-contour sweep for 3 concentrating-regime DATA3 NF270 sheets.

Mirrors legacy/data1_matlab/functions/calc_contour_3d.m + plot_grid_contour.m,
adapted for the DATA3 NF270 single-salt workflow.

Grid:
  L_p ∈ linspace(0.5, 50, 30)
  B   ∈ linspace(1e-6, B_max_per_salt, 30)
        where B_max = 30 (NaCl) / 15 (CaCl₂) / 10 (LaCl₃) — from NF270_B_BOUNDS_PER_SALT
  σ   ∈ linspace(0.0, 1.0, 10)

Sheets (all concentrating-regime, same MC2 workbook coupon):
  MC2.05.07.24_NaCl    — concentrating NaCl,  B upper = 30 µm/s
  MC2.05.07.24_CaCl2   — concentrating CaCl₂, B upper = 15 µm/s
  MC2.05.21.24_LaCl3   — concentrating LaCl₃, B upper = 10 µm/s

Output per sheet:
  paper_artifacts/nf270/contour3d/{run_id}/
    ├── grid3d.csv                  ← tidy long-form
    ├── _meta.json                  ← provenance
    ├── _failures.log               ← any failed grid cells
    └── {run_id}_grid3d_slices.png  ← 10-row × 3-col composite

Multiprocessing: configurable N_WORKERS per sheet (default 8), sheets processed
sequentially so each sheet gets the full pool.

Usage:
  python3 refactored_codes_v1/_run_3d_contour_concentrating.py                 # full run
  python3 refactored_codes_v1/_run_3d_contour_concentrating.py --smoke-test     # 3×3×3, 1 sheet
  python3 refactored_codes_v1/_run_3d_contour_concentrating.py --workers 4
  python3 refactored_codes_v1/_run_3d_contour_concentrating.py --sheets MC2.05.07.24_NaCl
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


# Defaults overridable on the command line.
# ALL 11 DATA3 NF270 single-salt sheets (the only scope we work in now).
SHEETS_DEFAULT = [
    "MC2.05.07.24_NaCl", "MC3.07.22.24_SNaCl", "MC4.07.11.24_SNaCl",
    "MC5.07.23.24_NaCl", "MC5.07.23.24_SNaCl", "MC5.07.23.24_S2NaCl",
    "MC2.05.07.24_CaCl2", "MC3.07.11.24_SCaCl2", "MC3.07.12.24_S2CaCl2",
    "MC2.05.21.24_LaCl3", "MC4.07.11.24_SLaCl3",
]
GRID_LP = 30
GRID_B  = 30
GRID_SIGMA = 10
NFE = 80
N_WORKERS = 8

# Measurement-error weighting (Lilonfe et al., ChemRxiv 2026): retentate cF = 2%
# relative, NO absolute floor.  Set EXPLICITLY in every worker so the 3D grids are
# byte-consistent with the canonical fits and the DATA1 σ×Lp/B×σ panels, regardless
# of any module-level drift.  (The legacy 0.3% weight gave WSSE ~44× too tight on cF.)
CF_RESIDUAL_SCALE_FRACTION = 0.02
CF_RESIDUAL_FLOOR_MM = None


# Per-worker globals — set by _worker_init in each subprocess
_DATA_STRU  = None
_RUN_ID     = None
_BASE_THETA = None   # warm-start theta (provides S0, S — fixed across the sweep)


def _load_warm_start_theta(run_id):
    """Centering θ* for the grid sweep.

    Prefer the CANONICAL per-form fit `bform_study/<run_id>/result_single.json` so the 3-D
    grid is centered on the SAME Lp*/B*/σ* as the σ×B-at-fixed-Lp slices, the AIC analysis,
    and the per-form campaign.  This matters: warm_start_fits/summary.json diverges badly
    on several sheets (e.g. CaCl₂ warm_start B=0.51 vs canonical 5.71), which would center
    the contour — and its fit-relative B range — on the wrong basin.  Fall back to
    warm_start_fits only if the per-form fit is missing.  S0/S (held fixed across the sweep)
    are present in both.
    """
    result = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts"
              / "nf270" / "bform_study" / run_id / "result_single.json")
    if result.exists():
        r = json.load(open(result))
        if isinstance(r, dict) and "parameters" in r:
            p = dict(r["parameters"])
            if all(k in p for k in ("Lp", "B", "sigma", "S0", "S")):
                return p
    summary_path = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                    / "paper_artifacts" / "nf270" / "warm_start_fits"
                    / "summary.json")
    summ = json.load(open(summary_path))
    for entry in summ:
        if entry.get("run_id") == run_id:
            return dict(entry["warm_start_theta_fit"])
    raise KeyError(f"{run_id}: no result_single.json and not in {summary_path}")


def _worker_init(run_id, base_theta):
    """Each worker imports the library and loads its sheet's data_stru ONCE.

    Subsequent _eval_cell calls reuse this in-process state, so we pay the
    ~5-10 s loadxlsx cost once per worker rather than once per cell.
    """
    global _DATA_STRU, _RUN_ID, _BASE_THETA
    import refactored_ucb_library as _lib
    # Pure-2% no-floor retentate weighting, set explicitly per worker (see header).
    _lib.NF270_CF_RESIDUAL_SCALE_FRACTION = CF_RESIDUAL_SCALE_FRACTION
    _lib.NF270_CF_RESIDUAL_FLOOR_MM = CF_RESIDUAL_FLOOR_MM
    family = _lib.NF270_RUN_REGISTRY[run_id]
    workbook_path = _lib.NF270_DEFAULT_ROOT / family["workbook"]
    sheet_name = family["sheet"]
    loaded = _lib.loadxlsx(workbook_path, sheet=sheet_name)
    _DATA_STRU = loaded["data_stru"]
    _RUN_ID = run_id
    _BASE_THETA = dict(base_theta)


def _eval_cell(args):
    """Forward-simulate one grid cell, return (Lp, B, sigma, log10 per-channel objs, err_or_None)."""
    import refactored_ucb_library as _lib

    Lp, B, sigma, nfe = args
    theta = dict(_BASE_THETA)
    theta["Lp"]    = float(Lp)
    theta["B"]     = float(B)
    theta["sigma"] = float(sigma)
    try:
        obj_m, obj_cv, obj_cr = _lib._nf270_contour_objectives_at_theta(
            _DATA_STRU, theta,
            mode="Lag", B_form="single",
            workflow_family="DATA3", nfe=nfe,
        )

        def _log10_or_nan(v):
            try:
                if v is None or not np.isfinite(v) or v <= 0:
                    return float("nan")
                return float(np.log10(v))
            except (TypeError, ValueError):
                return float("nan")

        return (Lp, B, sigma,
                _log10_or_nan(obj_m),
                _log10_or_nan(obj_cv),
                _log10_or_nan(obj_cr),
                None)
    except Exception as e:
        return (Lp, B, sigma,
                float("nan"), float("nan"), float("nan"),
                f"{type(e).__name__}: {e}")


def _salt_from_run_id(run_id):
    if "NaCl"  in run_id: return "NaCl"
    if "CaCl2" in run_id or "CaCl₂" in run_id: return "CaCl2"
    if "LaCl3" in run_id or "LaCl₃" in run_id: return "LaCl3"
    return None


def sweep_sheet(run_id, *, grid_lp, grid_b, grid_sigma, nfe, n_workers):
    import refactored_ucb_library as _lib

    salt = _salt_from_run_id(run_id)
    if salt not in _lib.NF270_B_BOUNDS_PER_SALT:
        raise KeyError(
            f"Salt {salt!r} (from run_id {run_id!r}) not in NF270_B_BOUNDS_PER_SALT. "
            f"Available: {list(_lib.NF270_B_BOUNDS_PER_SALT)}"
        )
    B_lower, B_upper = _lib.NF270_B_BOUNDS_PER_SALT[salt]

    # DATA1 calc_contour_2d convention: fit-RELATIVE axis ranges, so the well-identified
    # Lp* is resolved sharply.  The old fixed 0.5–50 / 0–B_upper grid snapped Lp* to a
    # coarse node (e.g. NaCl true 9.95 → grid 11.11) and spent most nodes on irrelevant
    # high Lp.  Lp ∈ [0.1·Lp*, 2·Lp*], B ∈ [0.1·|B*|, 3·|B*|], σ ∈ [0,1].
    base_theta = _load_warm_start_theta(run_id)
    Lp_fit = float(base_theta["Lp"]); B_abs = max(abs(float(base_theta["B"])), 1e-3)
    Lp_grid    = np.linspace(max(0.2, 0.1 * Lp_fit), 2.0 * Lp_fit, grid_lp)
    B_grid     = np.linspace(max(1e-6, 0.1 * B_abs), 3.0 * B_abs, grid_b)
    sigma_grid = np.linspace(0.0, 1.0, grid_sigma)

    cells = [(Lp, B, sg, nfe)
             for Lp in Lp_grid
             for B  in B_grid
             for sg in sigma_grid]
    n_cells = len(cells)

    save_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                / "paper_artifacts" / "nf270" / "contour3d" / run_id)
    save_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== {run_id} ===")
    print(f"  salt: {salt}    B upper: {B_upper}")
    print(f"  base theta from warm-start: Lp={base_theta['Lp']:.3f}  "
          f"B={base_theta['B']:.3f}  sigma={base_theta['sigma']:.3f}  "
          f"S0={base_theta.get('S0', 'NA')}  S={base_theta.get('S', float('nan')):.4f}")
    print(f"  grid: Lp={grid_lp}  B={grid_b}  sigma={grid_sigma}  →  {n_cells} cells")
    print(f"  workers: {n_workers}    nfe: {nfe}")
    print(f"  save_dir: {save_dir.relative_to(REPO_ROOT)}")

    t0 = time.time()
    with Pool(n_workers, initializer=_worker_init, initargs=(run_id, base_theta)) as pool:
        results = []
        chunksize = max(1, n_cells // (n_workers * 10))
        for i, res in enumerate(pool.imap_unordered(_eval_cell, cells, chunksize=chunksize), 1):
            results.append(res)
            if i % max(1, n_cells // 20) == 0 or i == n_cells:
                dt = time.time() - t0
                rate = i / max(dt, 1e-9)
                eta = (n_cells - i) / max(rate, 1e-9)
                print(f"  [{i:>5d} / {n_cells}]  ({100*i/n_cells:5.1f}%)  "
                      f"rate={rate:.2f} cells/s  eta={eta/60:.1f} min")
    dt = time.time() - t0
    print(f"  elapsed: {dt/60:.2f} min  ({dt/n_cells*1000:.0f} ms/cell)")

    df = pd.DataFrame(results, columns=[
        "Lp", "B", "sigma",
        "Obj_mass", "Obj_concentration", "Obj_retentate_concentration",
        "_err",
    ])
    failed = df[df["_err"].notna()]
    print(f"  failures: {len(failed)} / {n_cells}  ({100*len(failed)/n_cells:.2f}%)")

    csv_path = save_dir / "grid3d.csv"
    df.drop(columns=["_err"]).to_csv(csv_path, index=False)

    if len(failed) > 0:
        with open(save_dir / "_failures.log", "w") as f:
            for _, row in failed.iterrows():
                f.write(f"Lp={row.Lp:.4f} B={row.B:.4f} sigma={row.sigma:.4f} "
                        f":: {row._err}\n")

    meta = {
        "run_id":         run_id,
        "salt":           salt,
        "B_upper":        B_upper,
        "grid_Lp":        grid_lp,
        "grid_B":         grid_b,
        "grid_sigma":     grid_sigma,
        "total_cells":    n_cells,
        "failed_cells":   int(len(failed)),
        "wall_time_min":  dt / 60.0,
        "ms_per_cell":    dt / n_cells * 1000,
        "nfe":            nfe,
        "n_workers":      n_workers,
        "Lp_range":       [float(Lp_grid[0]), float(Lp_grid[-1])],
        "B_range":        [float(B_grid[0]),  float(B_grid[-1])],
        "sigma_range":    [float(sigma_grid[0]), float(sigma_grid[-1])],
        "cf_residual_scale_fraction": CF_RESIDUAL_SCALE_FRACTION,
        "cf_residual_floor_mM":       CF_RESIDUAL_FLOOR_MM,
        "weighting":      "retentate cF = 2% relative, NO floor (Lilonfe et al.); permeate cV 3%; mass 0.01 g",
    }
    with open(save_dir / "_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    return csv_path


def render_composite(run_id, csv_path):
    """Render 10-row × 3-col composite: one row per σ-slice, three channels."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(csv_path)
    save_dir = csv_path.parent

    channels = [
        ("Obj_mass",                    "Mass residual"),
        ("Obj_concentration",           "Permeate concentration residual"),
        ("Obj_retentate_concentration", "Retentate concentration residual"),
    ]

    sigma_slices = np.sort(df["sigma"].unique())
    nrows, ncols = len(sigma_slices), len(channels)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(4.5 * ncols, 3.0 * nrows),
        squeeze=False,
    )

    for i, sg in enumerate(sigma_slices):
        sub = df[np.isclose(df["sigma"], sg)]
        for j, (col, title) in enumerate(channels):
            ax = axes[i, j]
            piv = sub.pivot_table(index="B", columns="Lp", values=col, aggfunc="mean")
            X, Y = np.meshgrid(piv.columns.to_numpy(float), piv.index.to_numpy(float))
            Z = piv.to_numpy(float)
            if np.isfinite(Z).any():
                cp = ax.contour(X, Y, Z, 10, linewidths=1.5, cmap="viridis")
                ax.clabel(cp, cp.levels[::2], inline=True, fontsize=7,
                          colors="k", fmt="%1.1f")
                flat = np.nanargmin(Z)
                iy, ix = np.unravel_index(flat, Z.shape)
                x_min = float(X[iy, ix])
                y_min = float(Y[iy, ix])
                z_min = float(Z[iy, ix])
                ax.plot(x_min, y_min, "^", markersize=10,
                        markeredgecolor="red", markerfacecolor=[1, 0.6, 0.6],
                        clip_on=False, zorder=10)
                # Optimum text box (slice-level argmin: Lp, B at fixed σ).
                opt_str = (f"(L_p={x_min:.2g},  B={y_min:.2g},  "
                           f"log₁₀={z_min:.2f})")
                ax.text(0.5, -0.30, opt_str,
                        transform=ax.transAxes,
                        ha="center", va="top", fontsize=7,
                        bbox=dict(boxstyle="round,pad=0.25",
                                  facecolor="white",
                                  edgecolor=[0.6, 0.2, 0.2], linewidth=0.6))
            ax.set_xlabel("L_p [L/(m²·hr·bar)]", fontsize=9)
            ax.set_ylabel("B [µm/s]",             fontsize=9)
            ax.tick_params(direction="in", labelsize=8)
            if i == 0:
                ax.set_title(f"log₁₀  {title}", fontsize=10, fontweight="bold")
            if j == 0:
                ax.text(-0.30, 0.5, f"σ = {sg:.2f}",
                        transform=ax.transAxes, rotation=90,
                        va="center", fontsize=11, fontweight="bold")

    plt.suptitle(f"{run_id}   ·   3D objective contour   ·   "
                 f"L_p × B × σ slices",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0.02, 0, 1, 0.99])

    out_path = save_dir / f"{run_id}_grid3d_slices.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  rendered: {out_path.relative_to(REPO_ROOT)}")
    return out_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke-test", action="store_true",
                    help="Run 3×3×3 grid on first sheet only (~1 min). Calibrates timing.")
    ap.add_argument("--sheets", nargs="*", default=None,
                    help="Subset of sheets to run.")
    ap.add_argument("--workers", type=int, default=N_WORKERS)
    ap.add_argument("--grid-lp",    type=int, default=GRID_LP)
    ap.add_argument("--grid-b",     type=int, default=GRID_B)
    ap.add_argument("--grid-sigma", type=int, default=GRID_SIGMA)
    ap.add_argument("--nfe",        type=int, default=NFE)
    ap.add_argument("--no-render",  action="store_true",
                    help="Skip the per-sheet PNG render (CSV + meta only).")
    args = ap.parse_args()

    if args.smoke_test:
        grid_lp = grid_b = grid_sigma = 3
        sheets  = (args.sheets or SHEETS_DEFAULT)[:1]
        print(f"SMOKE TEST  ·  {grid_lp}×{grid_b}×{grid_sigma} = {grid_lp*grid_b*grid_sigma} cells "
              f"on {sheets[0]}")
    else:
        grid_lp    = args.grid_lp
        grid_b     = args.grid_b
        grid_sigma = args.grid_sigma
        sheets     = args.sheets or SHEETS_DEFAULT

    print(f"\nSheets:  {sheets}")
    print(f"Workers: {args.workers}")
    print(f"Grid:    Lp={grid_lp}  B={grid_b}  sigma={grid_sigma}  "
          f"= {grid_lp*grid_b*grid_sigma} cells per sheet")
    print(f"nfe:     {args.nfe}")

    overall_t0 = time.time()
    csvs = []
    for run_id in sheets:
        csv_path = sweep_sheet(
            run_id,
            grid_lp=grid_lp, grid_b=grid_b, grid_sigma=grid_sigma,
            nfe=args.nfe, n_workers=args.workers,
        )
        csvs.append(csv_path)
        if not args.no_render:
            try:
                render_composite(run_id, csv_path)
            except Exception as e:
                print(f"  RENDER FAILED: {type(e).__name__}: {e}")

    overall_dt = time.time() - overall_t0
    print(f"\n{'='*64}")
    print(f"DONE.  Total wall time: {overall_dt/60:.2f} min "
          f"({overall_dt:.0f} s)")
    print("Outputs:")
    for p in csvs:
        print(f"  {p.parent.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
