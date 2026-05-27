"""Per-sheet warm-anchor fit worker for sheets WITHOUT cF-floor contour panels.

Runs K internal trials, each calling solve_model(multistart=False) with
a perturbed copy of the MC2.05.07.24_NaCl seed θ. Picks the trial with
the lowest combined objective and renders DATA3 paper-style plots for
the winner.

Why orchestrate trials in the worker rather than relying on the
library's multistart?
  The library's `_theta_variants_for_multistart` hard-overrides Lp to
  5.0 for DATA3 and uses fixed LHS ranges (Lp ∈ [4.8, 5.2], B log-scale
  [1e-6, 5.0]) that DO NOT include the MC2 NaCl seed (Lp=9.11, B=8.31).
  So we cannot ask the library to "LHS around MC2 seed" — that capability
  is only in `build_seeded_multistart_starts`, which requires a panel
  dir these sheets don't have.

Toggles set:
  NF270_CF_RESIDUAL_FLOOR_MM         = 1.0    # interior σ-basin (1/2 of recipe)
  NF270_USE_PERMEATE_PROBE           = 0.03   # 3 % perm-probe channel (2/2 of recipe)
  NF270_MULTISTART_USE_CONTOUR_SEEDS = False  # no contour seeds

WHY BOTH KNOBS: see _contour_seeded_worker.py header — cF-floor alone
is insufficient above cF ≈ 333 mM, perm-probe anchors σ interior.
"""
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
os.chdir(str(ROOT.parent))

import refactored_ucb_library as lib

# Configure warm-anchor mode for this worker
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
lib.NF270_USE_PERMEATE_PROBE = 0.03      # 3 % perm-probe; required for interior σ
lib.NF270_MULTISTART_USE_CONTOUR_SEEDS = False

# MC2.05.07.24_NaCl converged θ (v11 fit; interior σ) — the anchor for the
# concentration-regime NaCl sheets per user table 2026-05-24.
MC2_NACL_SEED = {
    "Lp": 9.11,
    "B": 8.31,
    "beta_0": 8.31,
    "beta_1": 1.0,
    "sigma": 0.66,
    "S0": 0,
}

# Per-parameter relative jitter (LHS ±frac around seed)
JITTER_FRAC = {
    "Lp":     0.15,   # ±15 %  → [7.74, 10.48]
    "B":      0.20,   # ±20 %  → [6.65, 9.97]
    "beta_0": 0.20,   # ±20 %  (mirrors B for B_form=1)
    "sigma":  0.15,   # ±15 %  → [0.561, 0.759]
}
# Hard bounds (must match model_construct_inter)
BOUNDS = {
    "Lp":     (0.5,  50.0),
    "B":      (1e-6, 50.0),
    "beta_0": (1e-6, 50.0),
    "sigma":  (1e-3, 0.999),
}


def build_perturbed_seeds(base_theta, k_trials, rng_seed=13):
    """Return list of K theta dicts: first is the exact seed, rest are LHS jitter."""
    seeds = [dict(base_theta)]
    if k_trials <= 1:
        return seeds
    rng = np.random.default_rng(rng_seed)
    keys = list(JITTER_FRAC.keys())
    # Latin-hypercube via permuted uniform
    n = k_trials - 1
    u = (rng.permutation(np.arange(n)).astype(float) + rng.uniform(size=(n, len(keys)))) / n
    for row in u:
        trial = dict(base_theta)
        for i, key in enumerate(keys):
            if key not in base_theta:
                continue
            v0 = float(base_theta[key])
            frac = JITTER_FRAC[key]
            low = max(BOUNDS[key][0], v0 * (1.0 - frac))
            high = min(BOUNDS[key][1], v0 * (1.0 + frac))
            trial[key] = float(low + row[i] * (high - low))
        seeds.append(trial)
    return seeds


def _combined_obj(out):
    """Score = obj_m + obj_cv + obj_cr (lower is better). NaN → +inf."""
    vals = [out.get("obj_m", float("nan")),
            out.get("obj_cv", float("nan")),
            out.get("obj_cr", float("nan"))]
    s = 0.0
    for v in vals:
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return float("inf")
        s += float(v)
    return s


def main():
    run_id   = sys.argv[1]
    wb_path  = sys.argv[2]
    sheet    = sys.argv[3]
    out_json = sys.argv[4]
    save_dir = sys.argv[5]
    k_trials = int(sys.argv[6]) if len(sys.argv) > 6 else 4

    out = {"run_id": run_id, "k_trials": k_trials, "seed_source": "MC2.05.07.24_NaCl"}
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

    # ─── build K perturbed seeds ────────────────────────────────────────
    seeds = build_perturbed_seeds(MC2_NACL_SEED, k_trials, rng_seed=13)
    out["seeds_used"] = seeds

    # ─── run K trials sequentially ──────────────────────────────────────
    trials = []
    best_idx = None
    best_score = float("inf")
    best_sim_stru = None
    for idx, theta_init in enumerate(seeds):
        ti = time.time()
        trial = {"idx": idx, "theta_init": theta_init}
        try:
            fit_stru, sim_stru, _ = lib.solve_model(
                ds, ds["mode"], theta_init,
                sim_opt=False,
                B_form=1,
                LOUD=False,
                workflow_family="DATA3",
                multistart=False,
                solver_max_cpu_time=90,
            )
            trial["fit_time_s"] = round(time.time() - ti, 1)
            params = fit_stru.get("parameters", {})
            trial["parameters"] = {
                k: float(v) for k, v in params.items()
                if isinstance(v, (int, float))
            }
            trial["obj_m"]  = float(fit_stru.get("obj_m",  0))
            trial["obj_cv"] = float(fit_stru.get("obj_cv", 0))
            trial["obj_cr"] = float(fit_stru.get("obj_cr", 0))
            score = _combined_obj(trial)
            trial["combined_obj"] = score
            if score < best_score:
                best_score = score
                best_idx = idx
                best_sim_stru = sim_stru
        except Exception as exc:
            trial["error"] = f"{type(exc).__name__}: {exc}"
            trial["fit_time_s"] = round(time.time() - ti, 1)
            trial["traceback"] = traceback.format_exc()
        trials.append(trial)

    out["trials"] = trials
    out["total_time_s"] = round(time.time() - t0, 1)
    out["best_idx"] = best_idx

    if best_idx is None:
        out["error"] = "all trials failed"
        out["status"] = "failed"
        Path(out_json).write_text(json.dumps(out))
        return

    best_trial = trials[best_idx]
    out["parameters"]  = best_trial["parameters"]
    out["obj_m"]       = best_trial["obj_m"]
    out["obj_cv"]      = best_trial["obj_cv"]
    out["obj_cr"]      = best_trial["obj_cr"]
    out["combined_obj"] = best_trial["combined_obj"]
    out["status"] = "ok"

    # ─── render plot for the winning trial ──────────────────────────────
    try:
        results_dict = {
            "data": ds,
            "sim_stru": best_sim_stru,
            "model_settings": {"workflow_family": "DATA3"},
        }
        plot_paths = lib.run_data3_time_series_plots(
            results_dict, save_dir=save_dir, show=False
        )
        out["plot_paths"] = [str(p) for p in plot_paths]
    except Exception as exc:
        out["plot_error"] = f"{type(exc).__name__}: {exc}"
        out["plot_traceback"] = traceback.format_exc()

    Path(out_json).write_text(json.dumps(out))


if __name__ == "__main__":
    main()
