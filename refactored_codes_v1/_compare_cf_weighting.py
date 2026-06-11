#!/usr/bin/env python3
"""DATA3 before/after: does relaxing the retentate cF weight 0.3% -> 2% relieve
the sigma-wall pinning and improve channel balance?

For each single-salt sheet, re-fit from an INTERIOR sigma=0.5 seed (warm-start
Lp/B/S0/S) under each cF weighting and report where sigma lands + per-channel
WSSE. If 0.3% pulls sigma to a bound while 2% lets it sit interior, the
physically-justified weight demonstrably helps.
"""
import csv, json, sys, traceback
from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import os; os.chdir(str(ROOT.parent))
import refactored_ucb_library as lib

REPO = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo")
NF = REPO / "UnifiedFramework/DATA3/results/paper_artifacts/nf270"
WARM = json.loads((NF / "warm_start_fits/summary.json").read_text())
warm = {e["run_id"]: dict(e["warm_start_theta_fit"]) for e in WARM}
reg = lib.NF270_RUN_REGISTRY

CONFIGS = [("baseline_cF0.3%", 0.003), ("new_cF2%", 0.02)]
sheets = list(warm.keys())
print(f"sheets={len(sheets)}  configs={[c[0] for c in CONFIGS]}\n")

rows = []
for sheet in sheets:
    fam = reg.get(sheet)
    if not fam:
        continue
    try:
        ds = lib.loadxlsx(str(lib.NF270_DEFAULT_ROOT / fam["workbook"]), sheet=fam["sheet"])["data_stru"]
    except Exception as exc:
        print(f"LOAD FAIL {sheet}: {exc}"); continue
    seed = dict(warm[sheet]); seed["sigma"] = 0.5      # interior start
    rec = {"sheet": sheet}
    for name, cf_frac in CONFIGS:
        lib.NF270_CF_RESIDUAL_SCALE_FRACTION = cf_frac   # cV stays 3%
        lib.NF270_CP_RESIDUAL_SCALE_FRACTION = None
        try:
            fit, _, _ = lib.solve_model(ds, ds["mode"], dict(seed), sim_opt=False,
                                        B_form="single", workflow_family="DATA3",
                                        multistart=False, solver_max_cpu_time=120, LOUD=False)
            p = fit["parameters"]
            rec[name] = {"Lp": round(float(p["Lp"]), 3), "B": round(float(p["B"]), 3),
                         "sigma": round(float(p["sigma"]), 4),
                         "obj_m": round(float(fit.get("obj_m", 0)), 1),
                         "obj_cv": round(float(fit.get("obj_cv", 0)), 1),
                         "obj_cr": round(float(fit.get("obj_cr", 0)), 1)}
        except Exception as exc:
            rec[name] = {"error": f"{type(exc).__name__}: {exc}"}
            traceback.print_exc()
    rows.append(rec)
    b, n = rec.get("baseline_cF0.3%", {}), rec.get("new_cF2%", {})
    def wall(s):
        try: return "WALL" if (s["sigma"] >= 0.98 or s["sigma"] <= 0.02) else "interior"
        except Exception: return "err"
    print(f"{sheet:24s}  baseline sigma={b.get('sigma','?'):<7} [{wall(b)}]   "
          f"new sigma={n.get('sigma','?'):<7} [{wall(n)}]")

lib.NF270_CF_RESIDUAL_SCALE_FRACTION = 0.02
(NF / "meeting_prep_2026-06-10/cf_weighting_compare.json").write_text(json.dumps(rows, indent=2, default=str))
print("\nsaved -> meeting_prep_2026-06-10/cf_weighting_compare.json")
