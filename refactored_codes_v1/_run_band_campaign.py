#!/usr/bin/env python3
"""Campaign-wide per-band θ (NEXT_ACTIONS item 4): run the validated
solve_model_per_concentration_band on EVERY NF270 single-salt sheet and test
whether the GOOD/OKAY/POOR finding generalizes —
  (a) does the high-cF band recover an interior σ the whole-sheet fit rails away?
  (b) does B drift with concentration (low-band B vs high-band B) per sheet?
with per-band FIM non-singularity recorded.

Run: python3 _run_band_campaign.py
Outputs JSON + a campaign σ-recovery figure under .../nf270/band_value_test/.
"""
import sys, json, traceback
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib  # noqa: E402

OUT = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270" / "band_value_test"
OUT.mkdir(parents=True, exist_ok=True)


def _salt_of(rid):
    for s in ("S2NaCl", "SNaCl", "S2CaCl2", "SCaCl2", "NaCl", "CaCl2", "LaCl3", "SLaCl3"):
        if rid.endswith(s):
            return s.lstrip("S2").lstrip("S") or s
    return rid.split("_")[-1]


def _railed(x):
    return x is not None and (x < 0.02 or x > 0.98)


def _interior(x):
    return x is not None and 0.02 <= x <= 0.98


def run(rid):
    try:
        lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
        fam = lib.NF270_RUN_REGISTRY[rid]
        ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
        res = lib.solve_model_per_concentration_band(
            ds, ds["mode"], n_bands=2, B_form="single", workflow_family="DATA3",
            nfe=80, compute_fim=True, solver_max_cpu_time=120)
        pb = res["per_band"]
        if len(pb) < 2:
            return {"rid": rid, "salt": _salt_of(rid), "n_fitted": res["n_fitted"], "skip": "needs 2 bands"}
        lo, hi = pb[0], pb[1]   # low-cF, high-cF (partition is low->high)
        full_s = res["full"]["theta"].get("sigma")
        lo_s, hi_s = lo["theta"].get("sigma"), hi["theta"].get("sigma")
        lo_B, hi_B = lo["theta"].get("B"), hi["theta"].get("B")
        recovers = _railed(full_s) and (_interior(hi_s) or _interior(lo_s))
        return {
            "rid": rid, "salt": _salt_of(rid), "n_fitted": res["n_fitted"],
            "cf_lo": lo["cf_range"], "cf_hi": hi["cf_range"],
            "full_sigma": full_s, "lo_sigma": lo_s, "hi_sigma": hi_s,
            "lo_B": lo_B, "hi_B": hi_B,
            "B_drift": (hi_B - lo_B) if (lo_B is not None and hi_B is not None) else None,
            "recovers_interior_sigma": bool(recovers),
            "fim_singular": {"lo": lo.get("fim", {}).get("singular"), "hi": hi.get("fim", {}).get("singular")},
            "warnings": res["warnings"],
        }
    except Exception as exc:
        return {"rid": rid, "salt": _salt_of(rid), "error": repr(exc), "tb": traceback.format_exc()}


def figure(rows):
    ok = [r for r in rows if "full_sigma" in r and r.get("n_fitted", 0) >= 6]
    ok.sort(key=lambda r: (r.get("cf_hi") or [0, 0])[1])
    if not ok:
        return
    xs = np.arange(len(ok))
    fig, ax = plt.subplots(figsize=(max(8, 0.9 * len(ok)), 5.0))
    salt_col = {"NaCl": "#1f77b4", "CaCl2": "#2ca02c", "LaCl3": "#d62728"}
    for i, r in enumerate(ok):
        c = salt_col.get(r["salt"], "#777")
        ax.plot([i, i], [r["full_sigma"] or 0, r["hi_sigma"] or 0], color=c, lw=1.0, alpha=0.5, zorder=1)
        ax.scatter([i], [r["full_sigma"]], marker="o", s=80, facecolor="white", edgecolor=c, linewidths=1.6, zorder=3)
        ax.scatter([i], [r["hi_sigma"]], marker="*", s=190, facecolor=c, edgecolor="k", linewidths=1.0, zorder=4)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{r['rid'].split('_')[0]}\n{r['salt']}" for r in ok], fontsize=7, rotation=0)
    ax.set_ylabel("σ", fontweight="bold"); ax.set_ylim(-0.05, 1.05)
    ax.axhspan(0.98, 1.05, color="0.9", zorder=0); ax.axhspan(-0.05, 0.02, color="0.9", zorder=0)
    ax.set_title("Campaign per-band σ:  ○ whole-sheet fit   ★ high-cF band  (shaded = railed bounds)\n"
                 "star pulled off the wall = the high-cF band recovers an interior σ", fontsize=10)
    # legend
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", color="w", markerfacecolor="w", markeredgecolor="k", label="full-sheet σ"),
                       Line2D([], [], marker="*", color="w", markerfacecolor="0.5", markeredgecolor="k", markersize=12, label="high-cF band σ")],
              loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "band_campaign_sigma.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    sheets = [r for r in lib.NF270_SINGLE_SALT_RUNS if r in lib.NF270_RUN_REGISTRY]
    print(f"[campaign] {len(sheets)} single-salt sheets: {sheets}")
    with Pool(min(6, len(sheets))) as pool:
        rows = pool.map(run, sheets)
    with open(OUT / "band_campaign_summary.json", "w") as fh:
        json.dump(rows, fh, indent=2, default=float)
    figure(rows)
    print("\n=== CAMPAIGN PER-BAND σ-RECOVERY ===")
    n_rec = 0
    for r in sorted(rows, key=lambda r: r.get("rid", "")):
        if "error" in r:
            print(f"  {r['rid']:22s} ERROR {r['error'][:60]}"); continue
        if "skip" in r:
            print(f"  {r['rid']:22s} SKIP ({r['skip']}, n_fitted={r.get('n_fitted')})"); continue
        rec = "✓ RECOVERS interior σ" if r["recovers_interior_sigma"] else ""
        n_rec += int(r["recovers_interior_sigma"])
        bd = r["B_drift"]
        print(f"  {r['rid']:22s} {r['salt']:6s} nv={r['n_fitted']:2d}  "
              f"full σ={r['full_sigma']:.2f}  lo σ={r['lo_sigma']:.2f}  hi σ={r['hi_sigma']:.2f}  "
              f"B {r['lo_B']:.1f}→{r['hi_B']:.1f} (Δ{bd:+.1f})  {rec}")
    nq = sum(1 for r in rows if "full_sigma" in r and r.get("n_fitted", 0) >= 6)
    print(f"\n  {n_rec}/{nq} qualifying sheets: high-cF band recovers an interior σ the whole-sheet fit railed.")
    print("[campaign] DONE")
