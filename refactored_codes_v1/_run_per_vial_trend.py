#!/usr/bin/env python3
"""Per-vial theta TREND + vial<->contour correlation (Architecture.md §18.4 steps 2-3).

Rather than a joint per-vial-parameter model (m.Lp[n]/m.sigma[n], which would be
invasive model surgery), we trace how theta moves with vial concentration by
fitting theta to a ROLLING WINDOW of vials (reusing the validated objective band
mask, regularized by warm-starting from the full-sheet fit).  A 3-vial window
keeps each fit determined while still resolving a per-vial-position trend.

Outputs, per sheet (>=6 fitted vials):
  - theta (sigma, Lp, B) vs window-mean terminal cF  (the per-vial trend)
  - Pearson correlation of window cF with sigma and with B
  - the window optima overlaid on the sheet's sigma x Lp contour, colored by cF
    (the vial<->contour-group correlation)

Run: python3 _run_per_vial_trend.py
"""
import sys, json
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib  # noqa: E402

ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
OUT = ART / "band_value_test"
OUT.mkdir(parents=True, exist_ok=True)

SHEETS = [("GOOD", "MC3.07.22.24_SNaCl"),
          ("OKAY", "MC2.05.07.24_NaCl"),
          ("POOR", "MC2.05.21.24_LaCl3")]
WINDOW = 3


def _contour_csv(rid):
    p = ART / "matlab_ports_5ch_full" / rid / "contourdata-x_sigma-y_Lp.csv"
    return p if p.exists() else None


def run(arg):
    role, rid = arg
    lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
    fam = lib.NF270_RUN_REGISTRY[rid]
    ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
    mode = ds["mode"]
    n = ds["data_config"]["n"]; n_v0 = ds["data_config"].get("n_v0", 1)
    fitted = [i for i in range(1, n + 1) if i >= n_v0]
    if len(fitted) < 6:
        return {"role": role, "rid": rid, "skip": f"only {len(fitted)} fitted vials"}
    # order fitted vials by terminal cF (low->high) so windows are concentration-contiguous
    order = sorted(fitted, key=lambda i: lib._vial_terminal_cf(ds, i))
    windows = [order[k:k + WINDOW] for k in range(0, len(order) - WINDOW + 1)]
    # explicit-bands wrapper call (regularized: each window seeded from full fit)
    res = lib.solve_model_per_concentration_band(
        ds, mode, bands=windows, B_form="single", workflow_family="DATA3",
        nfe=80, compute_fim=False, min_vials=0, max_bands=len(windows),
        solver_max_cpu_time=120)
    pts = []
    for e in res["per_band"]:
        cf = float(np.mean([lib._vial_terminal_cf(ds, i) for i in e["band_vials"]]))
        t = e["theta"]
        pts.append({"vials": e["band_vials"], "cf_mean": cf,
                    "sigma": t.get("sigma"), "Lp": t.get("Lp"), "B": t.get("B"),
                    "obj_cr": e.get("obj_cr")})
    cfs = np.array([p["cf_mean"] for p in pts])
    sig = np.array([p["sigma"] for p in pts], dtype=float)
    Bs = np.array([p["B"] for p in pts], dtype=float)

    def _pearson(a, b):
        if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
            return None
        return float(np.corrcoef(a, b)[0, 1])

    return {"role": role, "rid": rid, "full": res["full"]["theta"],
            "points": pts, "corr_cf_sigma": _pearson(cfs, sig),
            "corr_cf_B": _pearson(cfs, Bs)}


def plot_sheet(r):
    if "skip" in r:
        return
    rid = r["rid"]
    pts = r["points"]
    cfs = [p["cf_mean"] for p in pts]
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    # (a) trend: sigma, B vs window cF
    ax0 = ax[0]
    ax0.plot(cfs, [p["sigma"] for p in pts], "o-", color="#d62728", label="σ")
    ax0.set_xlabel("window mean terminal cF (mM)"); ax0.set_ylabel("σ", color="#d62728")
    ax0.set_ylim(-0.05, 1.05)
    ax0b = ax0.twinx()
    ax0b.plot(cfs, [p["B"] for p in pts], "s--", color="#1f77b4", label="B")
    ax0b.set_ylabel("B (μm/s)", color="#1f77b4")
    ax0.set_title(f"{r['role']} {rid}\nper-window θ trend  "
                  f"(corr cF·σ={r['corr_cf_sigma']}, cF·B={r['corr_cf_B']})", fontsize=9)
    # (b) Lp vs cF
    ax[1].plot(cfs, [p["Lp"] for p in pts], "^-", color="#2ca02c")
    ax[1].set_xlabel("window mean terminal cF (mM)"); ax[1].set_ylabel("Lp")
    ax[1].set_title("Lp vs concentration window")
    # (c) vials<->contour: overlay window optima on σ×Lp contour, colored by cF
    axc = ax[2]
    csv = _contour_csv(rid)
    if csv is not None:
        import pandas as pd
        df = pd.read_csv(csv)
        col = "Obj_retentate_concentration" if "Obj_retentate_concentration" in df.columns else df.columns[2]
        piv = df.pivot_table(index="Lp", columns="sigma", values=col)
        X, Y = np.meshgrid(piv.columns.values, piv.index.values)
        axc.contourf(X, Y, piv.values, levels=16, cmap="Greys", alpha=0.6)
    sc = axc.scatter([p["sigma"] for p in pts], [p["Lp"] for p in pts],
                     c=cfs, cmap="plasma", s=120, edgecolors="k", linewidths=1.0, zorder=5)
    fig.colorbar(sc, ax=axc, label="window cF (mM)")
    axc.set_xlabel("σ"); axc.set_ylabel("Lp")
    axc.set_title("window optima on σ×Lp contour\n(color = concentration)")
    fig.tight_layout()
    out = OUT / f"per_vial_trend-{rid}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[per-vial] wrote {out}")


if __name__ == "__main__":
    with Pool(3) as pool:
        results = pool.map(run, SHEETS)
    with open(OUT / "per_vial_trend_summary.json", "w") as fh:
        json.dump(results, fh, indent=2, default=float)
    print("\n=== PER-VIAL θ TREND (rolling 3-vial windows, low→high cF) ===")
    for r in results:
        if "skip" in r:
            print(f"  {r['role']:4s} {r['rid']}: SKIP ({r['skip']})")
            continue
        print(f"\n  {r['role']:4s} {r['rid']}  corr(cF,σ)={r['corr_cf_sigma']}  corr(cF,B)={r['corr_cf_B']}")
        for p in r["points"]:
            print(f"     vials {str(p['vials']):16s} cF̄={p['cf_mean']:5.1f}mM  σ={p['sigma']:.3f} Lp={p['Lp']:.2f} B={p['B']:.2f}")
        plot_sheet(r)
    print("\n[per-vial] DONE")
