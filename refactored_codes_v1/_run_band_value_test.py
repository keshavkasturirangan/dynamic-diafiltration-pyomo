#!/usr/bin/env python3
"""Per-concentration-band VALUE TEST (Architecture.md §18.4, value-first gate).

For one GOOD, one OKAY, one POOR sheet: split the fitted vials into two
concentration bands by median terminal cF, fit one theta per band by MASKING the
objective to that band's vials (solve_model(..., band_vials=...)), and overlay
each band's (sigma, Lp) on that sheet's existing sigma x Lp contour.

The question this answers: do the two bands' optima land in DISTINGUISHABLE
basins (systematic sigma/B drift with concentration)?  If yes, per-band/per-vial
seeding is worth building out; the POOR sheet should show the bands smearing
into the flat (unidentified-sigma) region.

Run:  python3 _run_band_value_test.py
Outputs JSON + one overlay PNG per sheet under .../nf270/band_value_test/.
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

REPO = HERE.parent
ART = REPO / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
OUT = ART / "band_value_test"
OUT.mkdir(parents=True, exist_ok=True)

SHEETS = [
    ("GOOD", "MC3.07.22.24_SNaCl"),
    ("OKAY", "MC2.05.07.24_NaCl"),
    ("POOR", "MC2.05.21.24_LaCl3"),
]
NFE = 80
CPU = 200  # per-fit CPU-time cap (s)


def _terminal_cf(ds, n):
    a = np.asarray(ds["data_raw"][n - 1]["cF_exp"], dtype=float)
    a = a[np.isfinite(a)]
    return float(a[-1]) if a.size else np.nan


def _theta_row(fit):
    p = dict(fit.get("parameters", {}))
    return {
        "Lp": p.get("Lp"), "B": p.get("B"), "sigma": p.get("sigma"),
        "obj_m": fit.get("obj_m"), "obj_cv": fit.get("obj_cv"), "obj_cr": fit.get("obj_cr"),
    }


def _contour_csv(rid):
    p = ART / "matlab_ports_5ch_full" / rid / "contourdata-x_sigma-y_Lp.csv"
    return p if p.exists() else None


def run_sheet(arg):
    role, rid = arg
    try:
        lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0  # campaign config (per worker process)
        fam = lib.NF270_RUN_REGISTRY[rid]
        ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
        mode = ds["mode"]
        n = ds["data_config"]["n"]
        n_v0 = ds["data_config"].get("n_v0", 1)
        fitted = [i for i in range(1, n + 1) if i >= n_v0]
        tcf = {i: _terminal_cf(ds, i) for i in fitted}
        vals = sorted((tcf[i], i) for i in fitted if np.isfinite(tcf[i]))
        med = float(np.median([v for v, _ in vals]))
        low = sorted(i for v, i in vals if v <= med)     # low-cF band (more dilute)
        high = sorted(i for v, i in vals if v > med)     # high-cF band (more concentrated)

        # full-sheet seed fit
        full, _, _ = lib.solve_model(ds, mode, sim_opt=False, B_form="single",
                                     workflow_family="DATA3", nfe=NFE, LOUD=False,
                                     solver_max_cpu_time=CPU)
        seed = dict(full["parameters"])

        res = {
            "role": role, "rid": rid, "mode": mode, "n_fitted": len(fitted),
            "median_cf_mM": med, "terminal_cf_mM": {str(i): tcf[i] for i in fitted},
            "bands": {"low": low, "high": high},
            "cf_range": {
                "low": [min(tcf[i] for i in low), max(tcf[i] for i in low)],
                "high": [min(tcf[i] for i in high), max(tcf[i] for i in high)],
            },
            "full": _theta_row(full),
        }
        for label, band in (("low", low), ("high", high)):
            fit, _, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=False, B_form="single",
                                        workflow_family="DATA3", nfe=NFE, LOUD=False,
                                        band_vials=band, solver_max_cpu_time=CPU)
            res[label] = _theta_row(fit)
            res[label]["band_vials"] = band
        # quick identifiability signal: did sigma move between bands? did it rail?
        sl, sh = res["low"]["sigma"], res["high"]["sigma"]
        res["sigma_low"], res["sigma_high"] = sl, sh
        res["sigma_split"] = (abs(sh - sl) if (sl is not None and sh is not None) else None)
        res["sigma_railed"] = {
            "low": (sl is not None and (sl < 1e-3 or sl > 1 - 1e-3)),
            "high": (sh is not None and (sh < 1e-3 or sh > 1 - 1e-3)),
        }
        _overlay_plot(res, rid)
        return res
    except Exception as exc:
        return {"role": role, "rid": rid, "error": repr(exc), "tb": traceback.format_exc()}


def _overlay_plot(res, rid):
    csv = _contour_csv(rid)
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    if csv is not None:
        import pandas as pd
        df = pd.read_csv(csv)
        col = "Obj_retentate_concentration"
        if col not in df.columns:
            col = df.columns[2]
        piv = df.pivot_table(index="Lp", columns="sigma", values=col)
        X, Y = np.meshgrid(piv.columns.values, piv.index.values)
        Z = piv.values
        cs = ax.contourf(X, Y, Z, levels=18, cmap="viridis", alpha=0.85)
        ax.contour(X, Y, Z, levels=10, colors="k", linewidths=0.3, alpha=0.4)
        fig.colorbar(cs, ax=ax, label=f"log10 {col}")
    pts = [
        ("full sheet", res["full"], "white", "o"),
        (f"low-cF band {res['cf_range']['low'][0]:.0f}-{res['cf_range']['low'][1]:.0f}mM",
         res["low"], "#ff5555", "^"),
        (f"high-cF band {res['cf_range']['high'][0]:.0f}-{res['cf_range']['high'][1]:.0f}mM",
         res["high"], "#55aaff", "s"),
    ]
    for label, row, color, mk in pts:
        if row.get("sigma") is not None and row.get("Lp") is not None:
            ax.scatter([row["sigma"]], [row["Lp"]], s=160, marker=mk, c=color,
                       edgecolors="k", linewidths=1.4, zorder=5,
                       label=f"{label}: σ={row['sigma']:.2f} Lp={row['Lp']:.2f} B={row['B']:.2f}")
    ax.set_xlabel("σ (reflection coefficient)")
    ax.set_ylabel("Lp (L m⁻² h⁻¹ bar⁻¹)")
    ax.set_title(f"{res['role']} · {rid}\nper-band θ on the σ×Lp contour (B,Lp pinned at fit)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=7, frameon=False)
    fig.tight_layout()
    out = OUT / f"band_overlay-{rid}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    print(f"[band-test] running {len(SHEETS)} sheets -> {OUT}")
    with Pool(min(3, len(SHEETS))) as pool:
        results = pool.map(run_sheet, SHEETS)
    summary = OUT / "band_value_test_summary.json"
    with open(summary, "w") as fh:
        json.dump(results, fh, indent=2, default=float)
    print(f"\n[band-test] wrote {summary}")
    print("\n=== PER-BAND theta SUMMARY ===")
    for r in results:
        if "error" in r:
            print(f"  {r['role']:4s} {r['rid']}: ERROR {r['error']}")
            continue
        f, lo, hi = r["full"], r["low"], r["high"]
        print(f"\n  {r['role']:4s} {r['rid']}  (median cF {r['median_cf_mM']:.1f} mM)")
        print(f"     full : σ={f['sigma']:.3f} Lp={f['Lp']:.2f} B={f['B']:.2f}  WSSE_cr={f['obj_cr']:.3g}")
        print(f"     low  : σ={lo['sigma']:.3f} Lp={lo['Lp']:.2f} B={lo['B']:.2f}  WSSE_cr={lo['obj_cr']:.3g}  cF {r['cf_range']['low'][0]:.0f}-{r['cf_range']['low'][1]:.0f}mM railed={r['sigma_railed']['low']}")
        print(f"     high : σ={hi['sigma']:.3f} Lp={hi['Lp']:.2f} B={hi['B']:.2f}  WSSE_cr={hi['obj_cr']:.3g}  cF {r['cf_range']['high'][0]:.0f}-{r['cf_range']['high'][1]:.0f}mM railed={r['sigma_railed']['high']}")
        print(f"     σ split |high-low| = {r['sigma_split']:.3f}")
    print("\n[band-test] DONE")
