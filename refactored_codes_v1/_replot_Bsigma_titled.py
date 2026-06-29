#!/usr/bin/env python3
"""Re-plot the B-vs-sigma (fixed-Lp) identifiability contours from the saved
contourdata CSVs (NO re-solving), with the three panels titled simply
'mass', 'permeate', 'retentate' (the WSSE meaning is explained in the caption)
and the fixed Lp value in the figure title. v7 labelled-line-contour style.
Output: <Bsigma_fixedLp_nofloor>/<rid>/single/contour_titled.png
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
BS = (HERE.parent / "UnifiedFramework/DATA3/results/paper_artifacts/nf270"
      / "bform_study" / "Bsigma_fixedLp_nofloor")
REPS = {
    "diluting NaCl": "MC3.07.22.24_SNaCl",
    "concentrating NaCl": "MC2.05.07.24_NaCl",
    "concentrating CaCl₂": "MC2.05.07.24_CaCl2",
    "concentrating LaCl₃": "MC2.05.21.24_LaCl3",
}
PANELS = [("mass", "Obj_mass"), ("permeate", "Obj_concentration"),
          ("retentate", "Obj_retentate_concentration")]


def panel(ax, S, B, Z, title):
    levels = np.linspace(np.nanmin(Z), np.nanmax(Z), 14)
    cs = ax.contour(S, B, Z, levels=levels, cmap="turbo", linewidths=1.1)
    ax.clabel(cs, inline=True, fontsize=6, fmt="%.2f")
    j, i = np.unravel_index(np.nanargmin(Z), Z.shape)
    ax.plot(S[j, i], B[j, i], "^", color="red", ms=11, mec="k", mew=0.7, zorder=5)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("σ")
    ax.set_facecolor("white")


def main():
    for lab, rid in REPS.items():
        d = BS / rid / "single"
        csv = d / "contourdata-x_sigma-y_B.csv"
        if not csv.exists():
            print("MISSING", csv); continue
        Lp = json.load(open(d / "grid.json")).get("Lp_fixed", float("nan"))
        arr = np.genfromtxt(csv, delimiter=",", names=True)
        sig = np.unique(arr["sigma"]); Bv = np.unique(arr["B"])
        S, B = np.meshgrid(sig, Bv)
        fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3))
        for ax, (ttl, col) in zip(axes, PANELS):
            Z = arr[col].reshape(len(Bv), len(sig)) if arr[col].size == len(Bv) * len(sig) \
                else _pivot(arr, sig, Bv, col)
            panel(ax, S, B, Z, ttl)
        axes[0].set_ylabel("B  [µm s⁻¹]")
        fig.suptitle(f"{lab} ({rid.split('_')[0]}) — σ×B identifiability, "
                     f"hydraulic permeability fixed at its identified optimum Lₚ = {Lp:.2f} "
                     f"L m⁻² h⁻¹ bar⁻¹", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        out = d / "contour_titled.png"
        fig.savefig(out, dpi=160); plt.close(fig)
        print("WROTE", out)


def _pivot(arr, sig, Bv, col):
    Z = np.full((len(Bv), len(sig)), np.nan)
    si = {v: k for k, v in enumerate(sig)}; bi = {v: k for k, v in enumerate(Bv)}
    for row in arr:
        Z[bi[row["B"]], si[row["sigma"]]] = row[col]
    return Z


if __name__ == "__main__":
    main()
