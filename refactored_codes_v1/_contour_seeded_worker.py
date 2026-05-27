"""Per-sheet contour-seeded fit worker.

Runs as a subprocess so the parent can hard-kill the entire process group
(worker + IPOPT child) via os.killpg on timeout. Reads its config from
argv, writes a small JSON summary to disk, and renders DATA3 paper-style
plots in-process (avoids the sim_stru JSON-serialization issue that the
warm-start worker hit).

v5 (slide-29) RECIPE — restored 2026-05-25 with corrections:
  B_form                              = 2          # power-law (β₀ + β₁·cF^p, p=2)
  NF270_CF_RESIDUAL_FLOOR_MM          = 3.0        # matches slide 26 noise floor
  NF270_USE_PERMEATE_PROBE            = None       # OFF (slide-29 recipe omits it)
  NF270_MULTISTART_USE_CONTOUR_SEEDS  = True       # per-sheet contour minima
  NF270_MULTISTART_INCLUDE_CROSS_SALT = True       # MC3.SNaCl gold-standard anchor

B_form CORRECTION (user directive 2026-05-25 morning):
  "B_form to be used is B_form=1 or 2 and not be held const."
  My earlier pivot to B_form='single' was chasing a v5-log artifact and
  was wrong — the convection-diffusion power-law (B = β₀ + β₁·cF^p) is
  the correct model.  B_form=2 (quadratic in cF) is the more flexible
  choice and matches the slide-29 documented setup.

CF UNCERTAINTY (slide 26):
  At cF = 100 mM, actual probe noise ≈ 1-2 mM, but 0.003·cF = 0.3 mM
  (3-6× over-weighted).  At cF = 200 mM, noise ≈ 2-5 mM, but
  0.003·cF = 0.6 mM (3-8× over-weighted).  Floor of 3 mM brings the
  residual scale into agreement with measurement reality everywhere
  below cF ≈ 1000 mM.
"""
import json
import os
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
os.chdir(str(ROOT.parent))

import refactored_ucb_library as lib

# Configure slide-29 recipe with corrections (B_form=2, cF floor=3 mM)
lib.NF270_CF_RESIDUAL_FLOOR_MM      = 3.0    # bumped 1 → 3 (slide 26 noise floor)
lib.NF270_USE_PERMEATE_PROBE        = None   # OFF
lib.NF270_MULTISTART_USE_CONTOUR_SEEDS = True
lib.NF270_MULTISTART_INCLUDE_CROSS_SALT = True  # MC3.SNaCl gold standard


def main():
    run_id   = sys.argv[1]
    wb_path  = sys.argv[2]
    sheet    = sys.argv[3]
    out_json = sys.argv[4]
    save_dir = sys.argv[5]
    multistart_iterations = int(sys.argv[6]) if len(sys.argv) > 6 else 2

    out = {"run_id": run_id, "multistart_iterations": multistart_iterations}
    t0 = time.time()

    # ─── load sheet ──────────────────────────────────────────────────────
    try:
        ds = lib.loadxlsx(wb_path, sheet=sheet)["data_stru"]
    except Exception as exc:
        out["error"] = f"load: {type(exc).__name__}: {exc}"
        Path(out_json).write_text(json.dumps(out))
        return

    out["mode"]    = ds.get("mode")
    out["n_vials"] = ds["data_config"]["n"]
    out["M_F0"]    = float(ds["data_config"]["M_F0"])
    out["M_O"]     = float(ds["data_config"]["M_O"])

    # ─── verify cF-floor contour panel exists for this sheet ─────────────
    panel = lib._resolve_nf270_panel_dir(ds)
    out["contour_panel"] = str(panel) if panel else None
    if panel is None:
        out["error"] = "no contour panel found (neither cF-floor nor wall-geometry)"
        Path(out_json).write_text(json.dumps(out))
        return

    # ─── fit ─────────────────────────────────────────────────────────────
    try:
        fit_stru, sim_stru, _ = lib.solve_model(
            ds, ds["mode"], None,        # theta=None → library multistart picks seeds
            sim_opt=False,
            B_form=2,                    # power-law (β₀ + β₁·cF^p, p=2) per user directive
            LOUD=False,
            workflow_family="DATA3",
            multistart=True,
            multistart_iterations=multistart_iterations,
            solver_max_cpu_time=90,      # per-trial IPOPT cap
        )
        out["fit_time_s"] = round(time.time() - t0, 1)
        params = fit_stru.get("parameters", {})
        out["parameters"] = {
            k: float(v) for k, v in params.items()
            if isinstance(v, (int, float))
        }
        out["obj_m"]  = float(fit_stru.get("obj_m",  0))
        out["obj_cv"] = float(fit_stru.get("obj_cv", 0))
        out["obj_cr"] = float(fit_stru.get("obj_cr", 0))
        out["status"] = "ok"

        # ─── render plots in-process (no sim_stru serialization needed) ─
        try:
            # Build fit_meta for the new init-guess + cF-uncertainty annotations
            # (Phase B, 2026-05-25)
            fit_meta = {
                "recipe_name": "slide-29 v5 (B_form=2, cF floor 3 mM)",
                "b_form": 2,
                "theta_init_summary": (
                    "contour seeds: σ-Lp top + B-Lp top (per-slice)\n"
                    "cross-salt anchor: MC3.SNaCl (Lp=8.22, B=14.7, σ=0.45)\n"
                    f"LHS: {multistart_iterations} samples (±20% box)\n"
                    f"total trials: {2 + 1 + multistart_iterations}"
                ),
                "cf_uncertainty_str": (
                    f"σ_cF = max(0.003·cF, {lib.NF270_CF_RESIDUAL_FLOOR_MM} mM)"
                ),
                "winning_theta": out.get("parameters", {}),
            }
            results_dict = {
                "data": ds,
                "sim_stru": sim_stru,
                "model_settings": {"workflow_family": "DATA3"},
                "fit_meta": fit_meta,
            }
            plot_paths = lib.run_data3_time_series_plots(
                results_dict, save_dir=save_dir, show=False
            )
            # 2026-05-25: also render pressure + osmotic-pressure plots
            try:
                pressure_paths = lib.run_data3_pressure_plots(
                    results_dict, save_dir=save_dir, show=False
                )
                plot_paths = list(plot_paths) + list(pressure_paths)
            except Exception as exc:
                out["pressure_plot_error"] = f"{type(exc).__name__}: {exc}"
            out["plot_paths"] = [str(p) for p in plot_paths]
            out["fit_meta"] = fit_meta
        except Exception as exc:
            out["plot_error"] = f"{type(exc).__name__}: {exc}"
            out["plot_traceback"] = traceback.format_exc()

    except Exception as exc:
        out["error"] = f"fit: {type(exc).__name__}: {exc}"
        out["fit_time_s"] = round(time.time() - t0, 1)
        out["traceback"] = traceback.format_exc()

    Path(out_json).write_text(json.dumps(out))


if __name__ == "__main__":
    main()
