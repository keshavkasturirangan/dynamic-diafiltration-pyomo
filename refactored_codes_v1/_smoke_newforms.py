#!/usr/bin/env python3
"""Smoke test: confirm the new DATA3 B-forms ('sat','donnan') build, solve, and
return parameters on the GOOD NaCl sheet.  Also exercises the forward-sim path.
Run: python3 _smoke_newforms.py
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib
import numpy as np

lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
RID = "MC3.07.22.24_SNaCl"
fam = lib.NF270_RUN_REGISTRY[RID]
ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
mode = ds["mode"]
NFE = 80

# 1) constant-B seed fit
fs, _, _ = lib.solve_model(ds, mode, sim_opt=False, B_form="single",
                           workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=180)
ps = dict(fs["parameters"])
print(f"[single]  Lp={ps['Lp']:.3f}  B={ps['B']:.3f}  sigma={ps['sigma']:.3f}  WSSE={fs['Obj']:.1f}")

seed_common = {k: ps[k] for k in ("Lp", "sigma")}
for k in ("S0", "S"):
    if k in ps:
        seed_common[k] = ps[k]

# 2) saturating-exp fit
seed_sat = dict(seed_common); seed_sat.update({"B_inf": 1.0, "c_star": 20.0})
try:
    fsat, ssat, _ = lib.solve_model(ds, mode, theta=seed_sat, sim_opt=False, B_form="sat",
                                    workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=240)
    p = fsat["parameters"]
    print(f"[sat]     Lp={p['Lp']:.3f}  B_inf={p['B_inf']:.4g}  c_star={p['c_star']:.4g}  "
          f"sigma={p['sigma']:.3f}  WSSE={fsat['Obj']:.1f}")
    # forward-sim path at the fitted theta
    lib.solve_model(ds, mode, theta=p, sim_opt=True, B_form="sat",
                    workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=120)
    print("[sat]     forward-sim OK")
except Exception as e:
    import traceback; traceback.print_exc(); print("SAT FAILED:", e)

# 3) donnan-dielectric fit
seed_don = dict(seed_common); seed_don.update({"P0": 2.0, "X": 8.0, "k_dd": 0.5})
try:
    fdon, sdon, _ = lib.solve_model(ds, mode, theta=seed_don, sim_opt=False, B_form="donnan",
                                    workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=240)
    p = fdon["parameters"]
    print(f"[donnan]  Lp={p['Lp']:.3f}  P0={p['P0']:.4g}  X={p['X']:.4g}  k={p['k_dd']:.4g}  "
          f"sigma={p['sigma']:.3f}  WSSE={fdon['Obj']:.1f}")
    lib.solve_model(ds, mode, theta=p, sim_opt=True, B_form="donnan",
                    workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=120)
    print("[donnan]  forward-sim OK")
except Exception as e:
    import traceback; traceback.print_exc(); print("DONNAN FAILED:", e)

print("SMOKE DONE")
