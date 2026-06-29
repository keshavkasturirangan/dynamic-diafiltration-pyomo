#!/usr/bin/env python3
"""DATA3 §7a — stacking improves parameter estimability.

Per salt×regime group, compose the headline side-by-side:
   LEFT  = each individual experiment's σ×B contour (combined WSSE₃, on the COMMON grid)
   RIGHT = the single POOLED contour (summed = summed Fisher information) — a sharper basin.

Reads the canonical no-floor B×σ-at-fixed-Lp slices (Bsigma_fixedLp_nofloor/) already on disk.
Output: bform_study/Bsigma_fixedLp_nofloor/_stack/<group>/single/stacking_improvement.png
"""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _run_profile_contours as pc

VARIANT_DIR = "Bsigma_fixedLp_nofloor"


def _study():
    return pc._study()


def _wsse_grid(rid):
    od = _study() / VARIANT_DIR / rid / "single"
    g = json.loads((od / "grid.json").read_text())
    sig = np.array(g["axis_sigma"], float); B = np.array(g["axis_B"], float)
    Z = np.full((len(B), len(sig)), np.nan)
    for ln in (od / "nodes.jsonl").read_text().splitlines():
        try:
            r = json.loads(ln)
        except Exception:
            continue
        if r.get("status") == "ok":
            Z[r["j"], r["i"]] = r["WSSE3"]
    return sig, B, Z


def _basin_width(Z):
    """Count grid nodes within Δlog10 < 0.5 of the minimum — a coarse identifiability proxy
    (fewer = sharper basin = better-determined)."""
    if not np.any(~np.isnan(Z)):
        return None
    Zl = np.log10(np.clip(Z, 1e-6, None))
    return int(np.sum(Zl <= np.nanmin(Zl) + 0.5))


def make(group):
    sheets = [s for s in pc.POOL_GROUPS[group]
              if (_study() / VARIANT_DIR / s / "single" / "grid.json").exists()]
    if not sheets:
        print(f"[stack-improve] {group}: no slices"); return
    sig, B, _ = _wsse_grid(sheets[0]); X, Y = np.meshgrid(sig, B)
    Zs = {s: _wsse_grid(s)[2] for s in sheets}
    stack = np.zeros_like(X, float); cov = np.zeros_like(X, float)
    for s in sheets:
        m = ~np.isnan(Zs[s]); stack[m] += Zs[s][m]; cov += m
    stack[cov < len(sheets)] = np.nan

    npan = len(sheets) + 1
    fig, axes = plt.subplots(1, npan, figsize=(4.3 * npan, 4.5))
    if npan == 1:
        axes = [axes]

    def panel(ax, Z, title, pooled=False):
        if not np.any(~np.isnan(Z)):
            ax.text(0.5, 0.5, "(no overlap)", transform=ax.transAxes, ha="center"); return
        Zl = np.log10(np.clip(Z, 1e-6, None))
        cs = ax.contour(X, Y, Zl, levels=10, cmap="turbo", linewidths=1.3)
        ax.clabel(cs, inline=True, fontsize=5.5, fmt="%.2f")
        jm, im = np.unravel_index(np.nanargmin(Z), Z.shape)
        mn = (sig[im], B[jm])
        ax.plot(mn[0], mn[1], "^", ms=12, color="red", mec="white", mew=1.1, zorder=6)
        ax.text(0.5, -0.16, f"min σ={mn[0]:.2f}, B={mn[1]:.2f}",
                transform=ax.transAxes, ha="center", fontsize=7.5, color="red")
        ax.set_xlabel("$\\sigma$ [-]", fontsize=9); ax.set_ylabel("B@feed [$\\mu$m s$^{-1}$]", fontsize=9)
        ax.set_title(title, fontsize=9.5, fontweight="bold" if pooled else "normal")
        ax.grid(alpha=0.2, lw=0.4)
        if pooled:
            ax.set_facecolor((1.0, 0.96, 0.88))
        return mn

    mins = []
    for ax, s in zip(axes[:-1], sheets):
        mins.append(panel(ax, Zs[s], f"individual · {s.split('_')[0]}"))
    pmin = panel(axes[-1], stack, f"POOLED · {len(sheets)} stacked", pooled=True)
    # overlay each experiment's own minimum (○) on the pooled panel → scatter vs localization
    for (sm, bm) in mins:
        if sm is not None:
            axes[-1].plot(sm, bm, "o", ms=7, mfc="white", mec="0.15", mew=1.2, zorder=7)

    spread = "→".join(f"σ={m[0]:.2f}" for m in mins if m) if mins else ""
    fig.suptitle(f"{group} — pooling localizes the σ–B estimate "
                 f"(▲ each sheet's min: {spread};  ○ same mins on pooled panel → one consistent optimum)",
                 fontsize=8.5, y=0.98)
    fig.tight_layout(rect=[0.01, 0, 0.99, 0.90])
    od = _study() / VARIANT_DIR / "_stack" / group / "single"; od.mkdir(parents=True, exist_ok=True)
    out = od / "stacking_improvement.png"; fig.savefig(out, dpi=150); plt.close(fig)
    pooled_w = _basin_width(stack); indiv_w = [(_basin_width(Zs[s]), s.split('_')[0]) for s in sheets]
    print(f"[stack-improve] {group}: pooled basin={pooled_w}  vs individual={indiv_w}  -> {out}")


if __name__ == "__main__":
    groups = sys.argv[1:] or list(pc.POOL_GROUPS)
    for g in groups:
        make(g)
