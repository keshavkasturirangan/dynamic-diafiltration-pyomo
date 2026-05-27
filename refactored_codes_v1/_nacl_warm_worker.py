\
"""Per-sheet NaCl warm-start fit worker.

Runs as a subprocess so the parent can kill the whole process group
(worker + IPOPT child) via os.killpg on timeout. Reads its config from
argv, writes result to a JSON file.
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

# Enable cf-floor; disable contour seeds (May 22 wall-geometry incompatible
# with cf-floor's interior basin).
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
lib.NF270_MULTISTART_USE_CONTOUR_SEEDS = False

def main():
    run_id    = sys.argv[1]
    wb_path   = sys.argv[2]
    sheet     = sys.argv[3]
    theta_json = sys.argv[4]
    out_json  = sys.argv[5]

    theta = json.loads(theta_json)
    out = {"run_id": run_id, "seed_theta": theta}
    t0 = time.time()

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

    try:
        fit_stru, sim_stru, _ = lib.solve_model(
            ds, ds["mode"], theta,
            sim_opt=False,
            B_form=1,
            LOUD=False,
            workflow_family="DATA3",
            multistart=False,       # known-good seed, no LHS perturbation
            solver_max_cpu_time=240,
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
        out["sim_stru"] = [
            {
                "time": list(map(float, v.get("time", []))),
                "mF":   list(map(float, v.get("mF",   []))) if v.get("mF") is not None else None,
                "mV":   list(map(float, v.get("mV",   []))),
                "cF":   list(map(float, v.get("cF",   []))),
                "cH":   list(map(float, v.get("cH",   []))),
                "cV":   list(map(float, v.get("cV",   []))),
            }
            for v in sim_stru
        ]
    except Exception as exc:
        out["error"] = f"fit: {type(exc).__name__}: {exc}"
        out["fit_time_s"] = round(time.time() - t0, 1)
        out["traceback"] = traceback.format_exc()

    Path(out_json).write_text(json.dumps(out))

if __name__ == "__main__":
    main()
