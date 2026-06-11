"""Rigorous check: under cF 0.3% vs 2%, is the gold-standard MC3 SNaCl sigma-min
interior or at the wall? Fit from sigma seeds {0.2,0.5,0.8} under each weighting;
report the best (lowest 3-channel WSSE) and its sigma. If interior wins under 2%,
the single-fit wall result was a local min."""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent; sys.path.insert(0, str(ROOT))
import os; os.chdir(str(ROOT.parent))
import refactored_ucb_library as lib
NF = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/paper_artifacts/nf270")
warm = {e["run_id"]: dict(e["warm_start_theta_fit"]) for e in json.loads((NF/"warm_start_fits/summary.json").read_text())}
for sheet in ["MC3.07.22.24_SNaCl", "MC5.07.23.24_S2NaCl"]:
    fam = lib.NF270_RUN_REGISTRY[sheet]
    ds = lib.loadxlsx(str(lib.NF270_DEFAULT_ROOT/fam["workbook"]), sheet=fam["sheet"])["data_stru"]
    print(f"\n===== {sheet} =====")
    for name, cf in [("0.3%", 0.003), ("2%", 0.02)]:
        lib.NF270_CF_RESIDUAL_SCALE_FRACTION = cf; lib.NF270_CP_RESIDUAL_SCALE_FRACTION = None
        best = None
        for s0 in (0.2, 0.5, 0.8):
            seed = dict(warm[sheet]); seed["sigma"] = s0
            try:
                fit,_,_ = lib.solve_model(ds, ds["mode"], seed, sim_opt=False, B_form="single",
                                          workflow_family="DATA3", multistart=False, solver_max_cpu_time=120, LOUD=False)
                w3 = float(fit.get("obj_m",0))+float(fit.get("obj_cv",0))+float(fit.get("obj_cr",0))
                sg = float(fit["parameters"]["sigma"])
                if best is None or w3 < best[0]: best = (w3, sg, s0)
            except Exception as e:
                print(f"   seed {s0}: FAIL {e}")
        if best:
            print(f"  cF={name:4s}  best WSSE3={best[0]:.1f}  sigma={best[1]:.3f}  (from seed {best[2]})  -> {'WALL' if best[1]>=0.98 or best[1]<=0.02 else 'INTERIOR'}")
