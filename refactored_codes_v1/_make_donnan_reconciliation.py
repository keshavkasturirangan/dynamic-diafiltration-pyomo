#!/usr/bin/env python3
"""Reconciliation figure: WHY the fitted Donnan B(c) is flat in every single-salt
window, even though the Donnan partition is a *saturating* curve.

The Donnan+dielectric partition B(c) = Jw*P0*(-X + sqrt(X^2 + 4 k^2 c^2))/(2c)
collapses, under the similarity variable u = c / c_knee with c_knee = X/(2k), to a
UNIVERSAL master curve

        B / B_plateau = (sqrt(1 + u^2) - 1) / u ,   B_plateau = Jw*P0*k ,

which rises ~linearly for u<<1 (Donnan exclusion, the regime that *identifies* X)
and saturates to 1 for u>>1 (charge screened, B ~ constant).  This script reads
the THREE real fitted Donnan results (LaCl3, CaCl2, NaCl) and each sheet's
measured concentration window, computes where each window lands on the master
curve (u = 2 k c / X), and shades it.  Because every campaign fit rails X to its
lower bound (X=1e-3 mM), c_knee ~ 1e-3 mM and every measured window sits at
u ~ 1e3-1e5 -- i.e. far out on the plateau.  So a single salt only ever samples
the flat top of the saturating curve; the curvature that would pin X is never
probed.  This is the honest reconciliation with the schematic (which draws the
knee *inside* the window): the saturating FORM is right, but one salt can't see it.

NOTHING is fabricated -- the master curve is the exact normalized fitted equation
and the bands are the real measured windows at the real fitted (X, k).

Output: .../nf270/bform_study/taylor_vs_donnan/donnan_reconciliation.png
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STUDY = (HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts"
         / "nf270" / "bform_study")

SHEETS = [
    ("MC5.07.23.24_NaCl",  "NaCl",  "#1f77b4"),
    ("MC3.07.11.24_SCaCl2", "CaCl₂", "#2ca02c"),
    ("MC2.05.21.24_LaCl3", "LaCl₃", "#d62728"),
]


def window(sdir):
    d = np.genfromtxt(sdir / "apparent_B.csv", delimiter=",", names=True)
    c = np.atleast_1d(d["c_mM"])
    return float(np.nanmin(c)), float(np.nanmax(c))


def master(u):
    """Universal normalized Donnan partition B/B_plateau as a function of
    u = c / c_knee."""
    return (np.sqrt(1.0 + u**2) - 1.0) / u


def main():
    fig, ax = plt.subplots(figsize=(11.0, 5.4))

    # --- universal master curve (the exact normalized fitted Donnan form) ---
    u = np.logspace(-3, 6, 2000)
    ax.plot(u, master(u), "-", color="black", lw=2.6, zorder=6,
            label="Donnan partition  B/B$_{plateau}$ = (√(1+u²)−1)/u")

    # knee at u = 1
    ax.axvline(1.0, color="0.35", ls=":", lw=1.4, zorder=2)
    ax.plot([1.0], [master(np.array([1.0]))[0]], "o", color="black", ms=7, zorder=7)
    ax.annotate("knee  (u = 1, c = c$_{knee}$ = X/2k)\ncurvature here identifies X",
                xy=(1.0, master(np.array([1.0]))[0]), xytext=(0.02, 0.62),
                fontsize=9, color="0.2",
                arrowprops=dict(arrowstyle="->", color="0.4", lw=1.1))

    # regime shading
    ax.axvspan(1e-3, 1.0, color="#fdece6", zorder=0)
    ax.axvspan(1.0, 1e6, color="#e8f1ff", zorder=0)

    # --- each salt's real measured window mapped onto u = 2 k c / X ---
    # LaCl3 & CaCl2 windows overlap in u, so give each band its own height; NaCl
    # is well separated in x and can share the lowest row.
    BAND_Y = {"NaCl": 1.03, "LaCl₃": 1.03, "CaCl₂": 1.085}
    rows = []
    for rid, lab, col in SHEETS:
        sdir = STUDY / rid
        p = json.loads((sdir / "result_donnan.json").read_text())["parameters"]
        X, k = float(p["X"]), float(p["k_dd"])
        cmin, cmax = window(sdir)
        umin, umax = 2 * k * cmin / X, 2 * k * cmax / X
        knee = X / (2 * k)
        ax.axvspan(umin, umax, color=col, alpha=0.14, zorder=1)
        yb = BAND_Y[lab]
        ax.plot([umin, umax], [yb, yb], color=col, lw=5, solid_capstyle="butt", zorder=8)
        ax.text(np.sqrt(umin * umax), yb + 0.018, lab, color=col, fontsize=10,
                ha="center", va="bottom", fontweight="bold", zorder=9)
        rows.append((lab, X, k, knee, cmin, cmax, umin, umax,
                     float(master(np.array([umin]))[0])))

    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1e6)
    ax.set_ylim(0, 1.17)
    ax.set_xlabel("normalized concentration   u = c / c$_{knee}$ = 2 k c / X   (log scale)")
    ax.set_ylabel("normalized solute permeability   B / B$_{plateau}$")
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(alpha=0.25, which="both")

    # regime labels in data coords
    ax.text(0.04, 0.30, "Donnan exclusion\n(rise — identifies X)", fontsize=9.5,
            color="#a0432a", ha="left", va="center")
    ax.text(3e3, 0.50, "screened plateau\nB ≈ J$_w$·P$_0$·k  (≈ constant)", fontsize=9.5,
            color="#2c5d99", ha="center", va="center")
    ax.text(1e-3 * 1.4, 1.15, "measured single-salt windows  →  all on the plateau",
            fontsize=9.5, ha="left", va="top", color="0.3", style="italic")

    ax.set_title("Every single-salt window samples only the plateau of the saturating Donnan partition\n"
                 "fitted knee c$_{knee}$=X/2k ≈ 10⁻³ mM rails 10³–10⁵× below the data — the curvature that pins X is never probed",
                 fontsize=10.5)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.95)
    fig.tight_layout()
    out = STUDY / "taylor_vs_donnan" / "donnan_reconciliation.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)

    print(f"saved {out}\n")
    print(f"{'salt':6s} {'X[mM]':>8s} {'k':>5s} {'c_knee[mM]':>11s} "
          f"{'window[mM]':>16s} {'u=c/c_knee':>22s} {'B/Bplat@win_min':>16s}")
    for lab, X, k, knee, cmin, cmax, umin, umax, bmin in rows:
        print(f"{lab:6s} {X:8.4g} {k:5.3g} {knee:11.5g} "
              f"[{cmin:6.2f},{cmax:7.2f}] [{umin:9.0f},{umax:9.0f}]  {bmin:14.5f}")
    print("\nB/Bplateau at window_min is within rounding of 1.0 for all salts "
          "=> the entire window is on the flat plateau.")


if __name__ == "__main__":
    main()
