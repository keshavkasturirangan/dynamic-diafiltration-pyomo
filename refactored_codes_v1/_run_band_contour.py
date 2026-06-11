#!/usr/bin/env python3
"""Band-MASKED sigma x Lp contours: re-run the contour sweep with the objective
restricted to each concentration band's vials, to confirm the two bands occupy
genuinely different basins (Architecture.md §18.4, rigor step 3).

For a sheet, partition the fitted vials into 2 cF bands, fit one theta per band,
then sweep sigma x Lp (B pinned at the band fit) using ONLY that band's vials in
the objective.  Render the two band contours side by side with each band's
optimum marked, plus the full-sheet optimum for reference.

Run: python3 _run_band_contour.py [RUN_ID] [GRID]
Default RUN_ID=MC3.07.22.24_SNaCl (GOOD), GRID=16.
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

RUN_ID = sys.argv[1] if len(sys.argv) > 1 else "MC3.07.22.24_SNaCl"
GRID = int(sys.argv[2]) if len(sys.argv) > 2 else 16
NFE = 80

_DS = None
_BAND = None
_BASE = None


def _winit(rid, band, base):
    global _DS, _BAND, _BASE
    lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
    fam = lib.NF270_RUN_REGISTRY[rid]
    _DS = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
    _BAND = band
    _BASE = base


def _cell(arg):
    sigma, Lp = arg
    theta = dict(_BASE)
    theta["Lp"] = float(Lp)
    theta["sigma"] = float(sigma)
    try:
        fit, _, _ = lib.solve_model(_DS, _DS["mode"], theta=theta, sim_opt=True,
                                    B_form="single", workflow_family="DATA3", nfe=NFE,
                                    band_vials=_BAND, solver_max_iter=300,
                                    solver_max_cpu_time=15, LOUD=False)
        v = fit.get("obj_cr") if isinstance(fit, dict) else None
        if v is None or not np.isfinite(v) or v <= 0:
            return (sigma, Lp, float("nan"))
        return (sigma, Lp, float(np.log10(v)))
    except Exception:
        return (sigma, Lp, float("nan"))


def sweep_band(rid, band, base, grid):
    Lp_fit = float(base["Lp"])
    sig = np.linspace(0.0, 1.0, grid)
    lp = np.linspace(0.1 * Lp_fit, 2.0 * Lp_fit, grid)
    cells = [(s, l) for l in lp for s in sig]
    with Pool(6, initializer=_winit, initargs=(rid, band, base)) as pool:
        rows = pool.map(_cell, cells)
    return rows


def main():
    fam = lib.NF270_RUN_REGISTRY[RUN_ID]
    lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
    ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
    mode = ds["mode"]
    bands = lib.partition_vials_by_terminal_cf(ds, n_bands=2)
    print(f"[band-contour] {RUN_ID} bands={bands} grid={GRID}")

    # full-sheet fit (reference marker)
    full, _, _ = lib.solve_model(ds, mode, sim_opt=False, B_form="single",
                                 workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=200)
    full_t = dict(full["parameters"])

    panels = []
    for bi, band in enumerate(bands):
        cfs = [lib._vial_terminal_cf(ds, i) for i in band]
        fit, _, _ = lib.solve_model(ds, mode, theta=full_t, sim_opt=False, B_form="single",
                                    workflow_family="DATA3", nfe=NFE, band_vials=band,
                                    solver_max_cpu_time=200)
        bt = dict(fit["parameters"])
        print(f"[band-contour] band{bi} cF[{min(cfs):.0f}-{max(cfs):.0f}] θ: "
              f"σ={bt['sigma']:.3f} Lp={bt['Lp']:.2f} B={bt['B']:.2f}; sweeping...")
        rows = sweep_band(RUN_ID, band, bt, GRID)
        panels.append({"band": band, "cf": [min(cfs), max(cfs)], "theta": bt, "rows": rows})

    # render
    fig, axes = plt.subplots(1, len(panels), figsize=(6.2 * len(panels), 5.0), squeeze=False)
    for ax, p in zip(axes[0], panels):
        arr = np.array(p["rows"])
        sig = np.unique(arr[:, 0]); lp = np.unique(arr[:, 1])
        Z = arr[:, 2].reshape(len(lp), len(sig))
        X, Y = np.meshgrid(sig, lp)
        cs = ax.contourf(X, Y, Z, levels=18, cmap="viridis")
        ax.contour(X, Y, Z, levels=10, colors="k", linewidths=0.3, alpha=0.4)
        fig.colorbar(cs, ax=ax, label="log10 band WSSE_cr")
        t = p["theta"]
        ax.scatter([t["sigma"]], [t["Lp"]], s=180, marker="*", c="#ffe14d",
                   edgecolors="k", linewidths=1.4, zorder=6,
                   label=f"band θ: σ={t['sigma']:.2f} Lp={t['Lp']:.2f}")
        ax.scatter([full_t["sigma"]], [full_t["Lp"]], s=120, marker="o", c="white",
                   edgecolors="k", linewidths=1.2, zorder=5,
                   label=f"full θ: σ={full_t['sigma']:.2f} Lp={full_t['Lp']:.2f}")
        ax.set_xlabel("σ"); ax.set_ylabel("Lp")
        ax.set_title(f"band cF {p['cf'][0]:.0f}-{p['cf'][1]:.0f} mM (vials {p['band']})\nobjective MASKED to this band")
        ax.legend(loc="lower left", fontsize=7, frameon=True)
    fig.suptitle(f"{RUN_ID} · band-masked σ×Lp contours (B pinned per band)", y=1.02)
    fig.tight_layout()
    out_png = OUT / f"band_masked_contour-{RUN_ID}.png"
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close(fig)
    with open(OUT / f"band_masked_contour-{RUN_ID}.json", "w") as fh:
        json.dump({"run_id": RUN_ID, "full_theta": full_t,
                   "panels": [{"band": p["band"], "cf": p["cf"], "theta": p["theta"]} for p in panels]},
                  fh, indent=2, default=float)
    print(f"[band-contour] wrote {out_png}")


if __name__ == "__main__":
    main()
