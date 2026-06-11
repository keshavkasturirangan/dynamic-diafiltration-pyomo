#!/usr/bin/env python3
"""Re-render every existing contour PNG with the new 'optimum text box' beneath each panel.

Walks the three contour-output directories and re-renders the corresponding
PNG from the existing CSV — no recompute, no forward sims.  Adds the text
box (x_var, y_var, log10_min) below each panel, matching the MATLAB
``plot_contour.m`` legend convention but as a bordered box for visibility.

Directories handled:
  - paper_artifacts/nf270/matlab_ports/       (3-channel σ × Lp + B × Lp PNGs)
  - paper_artifacts/nf270/matlab_ports_5ch/   (5-channel σ × Lp PNGs)
  - paper_artifacts/nf270/contour3d/          (10×3 slice composites)

Run with:
    python3 refactored_codes_v1/_rerender_contours_with_optimum_box.py
"""

import os
import sys
import json
from pathlib import Path

REPO_ROOT = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo")
sys.path.insert(0, str(REPO_ROOT / "refactored_codes_v1"))
os.chdir(REPO_ROOT)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import refactored_ucb_library as lib


def rerender_matlab_port(csv_path: Path):
    """Re-render a 2D contour CSV (3-channel or 5-channel auto-detected by columns)."""
    df = pd.read_csv(csv_path)
    # Parse x_var, y_var from filename: contourdata-x_{X}-y_{Y}.csv
    stem = csv_path.stem
    parts = stem.replace("contourdata-", "").split("-")
    x_var = parts[0].replace("x_", "")
    y_var = parts[1].replace("y_", "")
    save_dir = csv_path.parent
    run_id = save_dir.name
    n_ch = "5-channel" if "Obj_permeate_conductivity" in df.columns else "3-channel"
    title = f"{run_id}  ·  {x_var} × {y_var}  ·  {n_ch}"
    out = lib.plot_contour_py(df, x_var, y_var, save_dir=save_dir, title=title)
    return out


def rerender_3d_contour(csv_path: Path):
    """Re-render a 10×3 σ-slice composite from grid3d.csv with the new text box.

    Mirrors the renderer in _run_3d_contour_concentrating.py.render_composite
    but adds the per-panel optimum text box.
    """
    df = pd.read_csv(csv_path)
    save_dir = csv_path.parent
    run_id = save_dir.name

    channels = [
        ("Obj_mass",                    "Mass residual"),
        ("Obj_concentration",           "Permeate concentration residual"),
        ("Obj_retentate_concentration", "Retentate concentration residual"),
    ]

    sigma_slices = np.sort(df["sigma"].unique())
    nrows, ncols = len(sigma_slices), len(channels)

    # Meta for the slide title
    meta_path = save_dir / "_meta.json"
    meta = json.load(open(meta_path)) if meta_path.exists() else {}
    gL = meta.get("grid_Lp", "?")
    gB = meta.get("grid_B",  "?")
    gS = meta.get("grid_sigma", "?")

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(4.5 * ncols, 3.0 * nrows),
        squeeze=False,
    )

    for i, sg in enumerate(sigma_slices):
        sub = df[np.isclose(df["sigma"], sg)]
        for j, (col, t) in enumerate(channels):
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
                # Per-slice optimum text box
                opt = f"(L_p={x_min:.2g},  B={y_min:.2g},  log₁₀={z_min:.2f})"
                ax.text(0.5, -0.30, opt,
                        transform=ax.transAxes,
                        ha="center", va="top", fontsize=7,
                        bbox=dict(boxstyle="round,pad=0.25",
                                  facecolor="white",
                                  edgecolor=[0.6, 0.2, 0.2], linewidth=0.6))
            ax.set_xlabel("L_p [L/(m²·hr·bar)]", fontsize=9)
            ax.set_ylabel("B [µm/s]",             fontsize=9)
            ax.tick_params(direction="in", labelsize=8)
            if i == 0:
                ax.set_title(f"log₁₀  {t}", fontsize=10, fontweight="bold")
            if j == 0:
                ax.text(-0.30, 0.5, f"σ = {sg:.2f}",
                        transform=ax.transAxes, rotation=90,
                        va="center", fontsize=11, fontweight="bold")

    plt.suptitle(f"{run_id}   ·   3D objective contour   ·   "
                 f"L_p × B × σ slices  ({gL}×{gB}×{gS} grid)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0.02, 0, 1, 0.99])

    out_path = save_dir / f"{run_id}_grid3d_slices.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    base = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"

    # ---- 2D MATLAB-port + 5-channel re-renders ----
    for subdir in ("matlab_ports", "matlab_ports_5ch"):
        root = base / subdir
        if not root.exists():
            continue
        print(f"\n=== {subdir} ===")
        for csv in sorted(root.rglob("contourdata-*.csv")):
            try:
                out = rerender_matlab_port(csv)
                print(f"  re-rendered: {out.relative_to(REPO_ROOT)}")
            except Exception as e:
                print(f"  FAILED {csv.relative_to(REPO_ROOT)}: {type(e).__name__}: {e}")

    # ---- 3D contour slice composites ----
    print(f"\n=== contour3d ===")
    cd_root = base / "contour3d"
    if cd_root.exists():
        for csv in sorted(cd_root.glob("*/grid3d.csv")):
            # Skip the legacy ~485 KB pre-fix CSVs (they're stale: all-NaN data
            # from yesterday's failed runs that hadn't been overwritten yet).
            # The fresh ones from the current run are larger (~950 KB).
            if csv.stat().st_size < 700_000:
                print(f"  SKIP stale: {csv.relative_to(REPO_ROOT)}  "
                      f"({csv.stat().st_size:,} bytes — pre-fix junk)")
                continue
            try:
                out = rerender_3d_contour(csv)
                print(f"  re-rendered: {out.relative_to(REPO_ROOT)}")
            except Exception as e:
                print(f"  FAILED {csv.relative_to(REPO_ROOT)}: {type(e).__name__}: {e}")

    print("\nDONE.")


if __name__ == "__main__":
    main()
