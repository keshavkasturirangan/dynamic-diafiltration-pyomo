#!/usr/bin/env python3
"""FAST beta-vs-Lp/sigma contours (Task 1) — multiprocessing + capped solver.

The library contour path (run_nf270_contour_for_sheet, B_form>=1) is correct but
each forward-sim cell can burn 3000 IPOPT iterations on extreme beta grid points.
This driver does the same two beta slices but with a 300-iter / 15 s cap per cell
and a worker pool, so a sheet finishes in minutes:

    (beta_0, Lp)      Lp on Y   (beta_0 plays B's role; L_p > B)
    (sigma,  beta_1)  beta_1 on Y  (beta > sigma)

Centering: warm-start B_form=1 from a 'single' fit (B -> beta_0, beta_1=0).
Run: python3 _run_beta_contour_fast.py [RUN_ID] [GRID]
Default RUN_ID=MC3.07.22.24_SNaCl, GRID=18.
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
OUT = ART / "beta_contours"
RUN_ID = sys.argv[1] if len(sys.argv) > 1 else "MC3.07.22.24_SNaCl"
GRID = int(sys.argv[2]) if len(sys.argv) > 2 else 18
NFE = 80

_DS = None
_BASE = None


def _winit(rid, base):
    global _DS, _BASE
    lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
    fam = lib.NF270_RUN_REGISTRY[rid]
    _DS = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
    _BASE = base


def _cell(arg):
    xname, xval, yname, yval = arg
    theta = dict(_BASE)
    theta[xname] = float(xval)
    theta[yname] = float(yval)
    try:
        fit, _, _ = lib.solve_model(_DS, _DS["mode"], theta=theta, sim_opt=True,
                                    B_form=1, workflow_family="DATA3", nfe=NFE,
                                    solver_max_iter=300, solver_max_cpu_time=15, LOUD=False)
        v = fit.get("obj_cr") if isinstance(fit, dict) else None
        if v is None or not np.isfinite(v) or v <= 0:
            return (xval, yval, float("nan"))
        return (xval, yval, float(np.log10(v)))
    except Exception:
        return (xval, yval, float("nan"))


def sweep(rid, base, xname, xvals, yname, yvals):
    cells = [(xname, xv, yname, yv) for yv in yvals for xv in xvals]
    with Pool(6, initializer=_winit, initargs=(rid, base)) as pool:
        return pool.map(_cell, cells)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fam = lib.NF270_RUN_REGISTRY[RUN_ID]
    lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
    ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
    mode = ds["mode"]

    # center: single -> B_form=1 warm fit
    fs, _, _ = lib.solve_model(ds, mode, sim_opt=False, B_form="single",
                               workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=180)
    ps = dict(fs["parameters"])
    seed = {"Lp": ps["Lp"], "sigma": ps["sigma"], "beta_0": ps["B"], "beta_1": 0.0}
    for k in ("S0", "S"):
        if k in ps:
            seed[k] = ps[k]
    f1, _, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=False, B_form=1,
                               workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=300)
    base = dict(f1["parameters"])
    print(f"[beta-fast] center θ: Lp={base['Lp']:.3f} beta_0={base['beta_0']:.4g} "
          f"beta_1={base['beta_1']:.4g} sigma={base['sigma']:.3f}")

    Lp, b0, b1 = float(base["Lp"]), float(base["beta_0"]), float(base["beta_1"])
    # beta_1 range derived from feasibility of B = beta_0 + beta_1*cIn over the
    # sheet's concentration span: B must stay in (0, B_max] up to cIn ~ cmax, so
    # beta_1 in (-beta_0/cmax, (B_max - beta_0)/cmax].  A naive ±2 span sweeps
    # mostly infeasible cells (B<0 or B over its Var bound) and renders blank.
    cmax = max(lib._vial_terminal_cf(ds, i)
               for i in range(1, ds["data_config"]["n"] + 1)
               if np.isfinite(lib._vial_terminal_cf(ds, i)))
    b_max = lib.NF270_B_BOUNDS_PER_SALT.get(
        str(ds["data_config"].get("namec") or "").strip(), lib.NF270_B_BOUNDS_DEFAULT)[1]
    b1_lo = -0.9 * b0 / cmax            # keep B > 0.1*beta_0 at cmax
    b1_hi = (b_max - b0) / cmax         # keep B <= bound at cmax
    sweeps = {
        "beta_0-Lp": dict(xname="beta_0", xvals=np.linspace(0.1 * abs(b0) + 1e-6, 3.0 * abs(b0) + 1e-3, GRID),
                          yname="Lp", yvals=np.linspace(0.1 * Lp, 2.0 * Lp, GRID)),
        "sigma-beta_1": dict(xname="sigma", xvals=np.linspace(0.0, 1.0, GRID),
                             yname="beta_1", yvals=np.linspace(b1_lo, b1_hi, GRID)),
    }
    panels = {}
    for key, s in sweeps.items():
        print(f"[beta-fast] sweeping {key} ...")
        rows = sweep(RUN_ID, base, s["xname"], s["xvals"], s["yname"], s["yvals"])
        panels[key] = {"rows": rows, **s}
        # persist CSV
        import csv
        with open(OUT / f"contourdata-{key}-{RUN_ID}.csv", "w", newline="") as fh:
            w = csv.writer(fh); w.writerow([s["xname"], s["yname"], "log10_obj_cr"])
            for r in rows:
                w.writerow(r)

    # render in the canonical figure_s5 contour-LINE style via lib._plot_heatmap_frame
    import pandas as pd
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    for ax, (key, p) in zip(axes, panels.items()):
        arr = np.array(p["rows"])
        df = pd.DataFrame({p["xname"]: arr[:, 0], p["yname"]: arr[:, 1],
                           "Obj_retentate_concentration": arr[:, 2]})
        lib._plot_heatmap_frame(df, p["xname"], p["yname"], "Obj_retentate_concentration",
                                ax=ax, show_title=False)
        ax.plot([base[p["xname"]]], [base[p["yname"]]], "*", markersize=18,
                markerfacecolor="#ffe14d", markeredgecolor="k", markeredgewidth=1.4,
                clip_on=False, zorder=7, label="fit θ")
        ax.set_title(f"{key}  ({p['yname']} on Y)", fontsize=10)
        ax.legend(loc="lower left", fontsize=8)
    fig.suptitle(f"{RUN_ID} · β contours (figure_s5 style; B_form=1: B=β₀+β₁·cIn)", y=1.02)
    fig.tight_layout()
    out_png = OUT / f"beta_contours-{RUN_ID}.png"
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.close(fig)
    with open(OUT / f"beta_contours-{RUN_ID}.json", "w") as fh:
        json.dump({"run_id": RUN_ID, "center_theta": base}, fh, indent=2, default=float)
    print(f"[beta-fast] wrote {out_png}")


if __name__ == "__main__":
    main()
