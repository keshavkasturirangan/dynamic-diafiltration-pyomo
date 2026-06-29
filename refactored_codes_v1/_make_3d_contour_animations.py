#!/usr/bin/env python3
"""Generate animated GIFs from existing 3D L_p × B × σ contour data.

For each sheet's grid3d.csv, render three animations:

   σ-scrub:  vary σ, each frame shows B vs L_p contour (10 frames @ 30×30)  ← matches A3/A4
   Lp-scrub: vary Lp, each frame shows B vs σ contour (30 frames @ 30×10)
   B-scrub:  vary B, each frame shows L_p vs σ contour (30 frames @ 30×10)

Each GIF is a 1×3 composite (3 channels: mass, perm-conc, ret-conc).

Output:
   paper_artifacts/nf270/animations3d/{run_id}/
       ├── scrub-sigma_BvsLp.gif        ← σ-scrub
       ├── scrub-Lp_Bvssigma.gif        ← Lp-scrub
       └── scrub-B_Lpvssigma.gif        ← B-scrub

PowerPoint plays embedded GIFs natively, so these drop into appendix slides A3/A4
as replacements for the static slice composites.

Pure rendering — no compute. ~30-60 sec per sheet.
"""

import os
import sys
import io
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

REPO = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo")
BASE3D = REPO / "UnifiedFramework/DATA3/results/paper_artifacts/nf270/contour3d"
OUT_BASE = REPO / "UnifiedFramework/DATA3/results/paper_artifacts/nf270/animations3d"

CHANNELS = [
    ("Obj_mass",                    "Mass residual"),
    ("Obj_concentration",           "Permeate concentration residual"),
    ("Obj_retentate_concentration", "Retentate concentration residual"),
]


def render_frame(df_slice, *, x_var, y_var, scrub_var, scrub_val, run_id, n_levels=10,
                 figsize=(14, 4.5), dpi=110):
    """Render one frame: 1×3 panels (3 channels) at a fixed slice."""
    labels = {
        "Lp":    r"$L_p$  [L / m² / hr / bar]",
        "B":     r"$B$  [µm / s]",
        "sigma": r"$\sigma$  [dimensionless]",
    }
    fig, axes = plt.subplots(1, 3, figsize=figsize, dpi=dpi)
    for ax, (col, title) in zip(axes, CHANNELS):
        piv = df_slice.pivot_table(index=y_var, columns=x_var, values=col, aggfunc="mean")
        X = piv.columns.to_numpy(float)
        Y = piv.index.to_numpy(float)
        XX, YY = np.meshgrid(X, Y)
        Z = piv.to_numpy(float)
        if np.isfinite(Z).any():
            cp = ax.contour(XX, YY, Z, n_levels, linewidths=1.5, cmap="viridis")
            ax.clabel(cp, inline=True, fontsize=8, fmt="%.1f")
            flat = np.nanargmin(Z)
            iy, ix = np.unravel_index(flat, Z.shape)
            ax.plot(XX[iy, ix], YY[iy, ix], "^", markersize=10,
                    markeredgecolor="red", markerfacecolor=[1, 0.6, 0.6])
        ax.set_xlabel(labels[x_var], fontsize=10)
        ax.set_ylabel(labels[y_var], fontsize=10)
        ax.set_title(rf"$\log_{{10}}$  {title}", fontsize=11, fontweight="bold")
        ax.tick_params(labelsize=9)

    # Format scrub-variable label.  Only sigma has a LaTeX symbol; Lp and B are plain text.
    scrub_label = {
        "sigma": rf"$\sigma$ = {scrub_val:.3g}",
        "Lp":    f"L_p = {scrub_val:.2f}",
        "B":     f"B = {scrub_val:.3g} µm/s",
    }[scrub_var]
    fig.suptitle(f"{run_id}   ·   {x_var} × {y_var}  at  {scrub_label}",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0.02, 1, 0.93])

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def make_animation(df, *, scrub_var, x_var, y_var, run_id, out_path,
                   frame_duration_ms=400, fps_note=None):
    """Make one animation by scrubbing through unique scrub_var values."""
    scrub_vals = np.sort(df[scrub_var].unique())
    frames = []
    print(f"   ⏳ {scrub_var}-scrub  ({len(scrub_vals)} frames, {x_var} × {y_var})")
    for sv in scrub_vals:
        sub = df[np.isclose(df[scrub_var], sv)]
        if sub[CHANNELS[0][0]].notna().sum() == 0:
            continue
        frame = render_frame(sub, x_var=x_var, y_var=y_var,
                             scrub_var=scrub_var, scrub_val=sv, run_id=run_id)
        frames.append(frame)
    if not frames:
        print(f"   ✗ no valid frames")
        return None

    # Ping-pong: forward + reversed (excluding endpoints) for smoother loop
    pingpong = frames + frames[-2:0:-1]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        out_path,
        save_all=True,
        append_images=pingpong[1:],
        duration=frame_duration_ms,
        loop=0,
        optimize=True,
    )
    size_kb = out_path.stat().st_size / 1024
    print(f"   ✓ wrote: {out_path.relative_to(REPO)}  ({len(pingpong)} frames, {size_kb:.0f} KB)")
    return out_path


def animate_sheet(run_id, *, frame_ms=400):
    csv = BASE3D / run_id / "grid3d.csv"
    if not csv.exists():
        print(f"  SKIP {run_id} — no grid3d.csv")
        return []

    # Skip stale (tiny) CSVs from failed/smoke runs.  Threshold lowered to 100 KB so the
    # 2026-06-18 pure-2% regen grids (15×12×10 ≈ 180 KB) are NOT mistaken for junk; a
    # 3×3×3 smoke test is ~2 KB so it's still excluded.
    if csv.stat().st_size < 100_000:
        print(f"  SKIP {run_id} — stale CSV ({csv.stat().st_size:,} bytes, likely smoke/junk)")
        return []

    df = pd.read_csv(csv)
    out_dir = OUT_BASE / run_id
    print(f"\n=== {run_id} ===")
    print(f"  3D grid: Lp={df.Lp.nunique()}, B={df.B.nunique()}, σ={df.sigma.nunique()}, "
          f"{df[CHANNELS[0][0]].notna().sum()}/{len(df)} finite")

    outputs = []
    # AXIS CONVENTION (rule-of-thumb, 2026-06-11): Y-axis priority L_p > B > sigma.
    # The higher-priority variable goes on Y; the other on X.  All three scrubs
    # below comply: L_p×B (L_p on Y), B×sigma (B on Y), L_p×sigma (L_p on Y).
    #
    # sigma-scrub: vary sigma, contour L_p vs B (L_p on Y per the rule)
    p = make_animation(df, scrub_var="sigma", x_var="B", y_var="Lp",
                        run_id=run_id, out_path=out_dir / "scrub-sigma_LpvsB.gif",
                        frame_duration_ms=frame_ms)
    if p: outputs.append(p)

    # Lp-scrub: vary Lp, contour B vs sigma (B on Y per the rule)
    p = make_animation(df, scrub_var="Lp", x_var="sigma", y_var="B",
                        run_id=run_id, out_path=out_dir / "scrub-Lp_Bvssigma.gif",
                        frame_duration_ms=frame_ms // 2)  # more frames, snappier
    if p: outputs.append(p)

    # B-scrub: vary B, contour L_p vs sigma (L_p on Y per the rule)
    p = make_animation(df, scrub_var="B", x_var="sigma", y_var="Lp",
                        run_id=run_id, out_path=out_dir / "scrub-B_Lpvssigma.gif",
                        frame_duration_ms=frame_ms // 2)
    if p: outputs.append(p)

    return outputs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheets", nargs="*", default=None,
                    help="Subset of sheets. Default: every sheet with valid 3D data.")
    ap.add_argument("--frame-ms", type=int, default=400,
                    help="σ-scrub frame duration in ms (Lp/B scrubs use half).")
    args = ap.parse_args()

    # Find sheets with valid 3D data
    if args.sheets:
        sheets = args.sheets
    else:
        sheets = []
        for d in sorted(BASE3D.iterdir()):
            if not d.is_dir():
                continue
            csv = d / "grid3d.csv"
            if csv.exists() and csv.stat().st_size >= 100_000:
                sheets.append(d.name)

    print(f"Sheets with valid 3D data ({len(sheets)}):  {sheets}")
    total = 0
    for sid in sheets:
        outs = animate_sheet(sid, frame_ms=args.frame_ms)
        total += len(outs)
    print(f"\n{'='*60}")
    print(f"DONE.  Wrote {total} animations across {len(sheets)} sheets.")
    print(f"Index: {OUT_BASE}")


if __name__ == "__main__":
    main()
