#!/usr/bin/env python3
"""Curate NF270 holdup boundaries via DATA2-style collection-line extrapolation.

For each single-salt sheet, fit a line to the steady early collection region
(cumulative mass in [0.15, 0.75] g, after the auto boundary) and extrapolate it
back to the baseline (mass = 0) crossing. That crossing time is the true onset
of permeate collection — the same idea as DATA2's mass extrapolation — and is
written into nf270_holdup_boundaries.csv `holdup_end_s` as the curated boundary.

This handles both clean sheets (curated ~ auto) and gradual-onset sheets (the
linear extrapolation skips the sub-threshold creep). Comparison plots with the
auto (red) and curated (green) boundaries are saved to holdup_diagnostics/curated/.

Usage: python3 _curate_holdup_boundaries.py
"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib

CSV = HERE / "nf270_holdup_boundaries.csv"
DIAG = HERE / "holdup_diagnostics" / "curated"
DIAG.mkdir(parents=True, exist_ok=True)

FIT_LO, FIT_HI = 0.15, 0.75   # g — steady early-collection window for the line fit


def _full_vial1(ds):
    """Reconstruct the pre-split vial-1 mass/time and the auto boundary index."""
    cfg = ds["data_config"]; dr = ds["data_raw"]
    v1 = dr[0]
    t1 = np.asarray(v1["time"], float).reshape(-1)
    m1 = np.asarray(v1["mass"], float).reshape(-1)
    if cfg.get("vial1_split_status") == "split" and len(dr) > 1:
        v2 = dr[1]
        t2 = np.asarray(v2["time"], float).reshape(-1)
        m2 = np.asarray(v2["mass"], float).reshape(-1)
        t = np.concatenate([t1, t2])
        m = np.concatenate([m1, m2 + (m1[-1] if m1.size else 0.0)])
        i_b = len(t1)
    else:
        t, m, i_b = t1, m1, int(cfg.get("n_holdup_samples", 0) or 0)
    return t, m, i_b


def _extrapolated_onset(t, m, i_b):
    """Time where the steady collection line crosses mass = 0."""
    hi = max(FIT_HI, 0.0)
    sel = (m >= FIT_LO) & (m <= hi)
    sel[:i_b] = False                      # only the collection side
    if sel.sum() < 5:                      # widen if the window is thin
        sel = np.arange(len(m)) >= i_b
        sel &= np.isfinite(m)
    if sel.sum() < 3:
        return None
    a, b = np.polyfit(t[sel], m[sel], 1)   # mass = a*t + b
    if not np.isfinite(a) or a <= 0:
        return None
    return float(-b / a)                    # mass = 0 crossing


def main():
    df = pd.read_csv(CSV)
    rows = []
    for _, r in df.iterrows():
        rid = r["run_id"]
        if rid not in lib.NF270_RUN_REGISTRY:
            rows.append(r.to_dict()); continue
        fam = lib.NF270_RUN_REGISTRY[rid]
        ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
        t, m, i_b = _full_vial1(ds)
        auto_s = float(t[i_b]) if 0 <= i_b < len(t) else float("nan")
        cur = _extrapolated_onset(t, m, i_b)
        d = r.to_dict()
        if cur is not None:
            cur = max(float(t[0]), min(float(t[-1]), cur))   # clamp into the vial
            d["holdup_end_s"] = round(cur, 1)
            d["notes"] = "collection-line extrapolation to mass=0"
        rows.append(d)
        # comparison plot
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(t, m, ".", ms=4, color="tab:blue", label="vial-1 mass")
        ax.axhline(0.0, color="k", lw=0.8)
        ax.axvline(auto_s, color="red", ls="--", lw=1.5, label=f"auto = {auto_s:.0f} s")
        if cur is not None:
            ax.axvline(cur, color="green", ls="-", lw=1.5, label=f"curated = {cur:.0f} s")
            tl = np.linspace(cur, t[min(len(t) - 1, i_b + 40)], 30)
            a, b = np.polyfit(t[(m >= FIT_LO) & (m <= FIT_HI) & (np.arange(len(m)) >= i_b)],
                              m[(m >= FIT_LO) & (m <= FIT_HI) & (np.arange(len(m)) >= i_b)], 1)
            ax.plot(tl, a * tl + b, "g:", lw=1, label="collection line")
        ax.set_xlim(0, (auto_s if np.isfinite(auto_s) else t[-1]) * 2.2)
        ax.set_ylim(-0.1, 0.9)
        ax.set_xlabel("Time [s]"); ax.set_ylabel("Cumulative mass [g]")
        ax.set_title(f"{rid}  auto vs curated")
        ax.legend(fontsize=8, loc="best")
        fig.tight_layout(); fig.savefig(DIAG / f"{rid}.png", dpi=130); plt.close(fig)
        dshift = (cur - auto_s) if cur is not None else float("nan")
        print(f"  {rid:24s} auto={auto_s:7.1f}  curated={cur if cur is None else round(cur,1):>7}  Δ={dshift:+6.1f} s")

    out = pd.DataFrame(rows, columns=list(df.columns))
    out.to_csv(CSV, index=False)
    print(f"\nUPDATED {CSV}")
    print(f"Comparison plots in {DIAG}/")


if __name__ == "__main__":
    main()
