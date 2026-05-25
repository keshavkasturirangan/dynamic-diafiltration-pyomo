
"""One-sheet fit worker. Invoked as a subprocess so the parent can
hard-kill it via OS timeout regardless of what IPOPT is doing."""
import json, os, sys, time, traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
os.chdir(str(ROOT.parent))

import refactored_ucb_library as lib

def main():
    run_id  = sys.argv[1]
    wb_path = sys.argv[2]
    sheet   = sys.argv[3]
    out_json = sys.argv[4]

    out = {"run_id": run_id}
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
            ds, ds["mode"], None,        # theta=None → per-sheet seeds
            sim_opt=False,
            B_form=1,                    # DATA2 recipe
            LOUD=False,
            workflow_family="DATA3",
            solver_max_cpu_time=90,
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
        # Also dump sim_stru so the parent can render plots without re-fitting
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
