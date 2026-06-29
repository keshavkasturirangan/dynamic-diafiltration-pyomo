#!/usr/bin/env python3
"""Single fixed-k Donnan DAE fit for one sheet, run as an isolated subprocess so
the parent (_recover_donnan.py) can hard-timeout it.  Writes result_donnan.json
ONLY on success (a timeout-kill therefore leaves any existing file intact).

Usage: python3 _donnan_fit_worker.py <run_id> <k_fixed> <salt>
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _run_bform_campaign as camp
import refactored_ucb_library as lib

camp.NFE = 80   # coarser mesh: the setup (single+sim+linear) must fit in the recovery timeout

rid = sys.argv[1]
kfix = float(sys.argv[2])
salt = sys.argv[3] if len(sys.argv) > 3 else ""
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
lib.NF270_DONNAN_FIX_K = kfix   # -> m.k_dd becomes a fixed Param (2-param P0,X fit)

ds = camp.load(rid)
mode = ds["mode"]
fs, _, _ = lib.solve_model(ds, mode, sim_opt=False, B_form="single",
                           workflow_family="DATA3", nfe=camp.NFE, solver_max_cpu_time=200, LOUD=False)
base = dict(fs["parameters"])
c, bp, bs = camp.apparent_partition(ds, base)
sl = {k: base[k] for k in ("Lp", "sigma") if k in base}
for k in ("S0", "S"):
    if k in base:
        sl[k] = base[k]
sl["beta_0"] = base.get("B", 1.0); sl["beta_1"] = 0.0
fl, _, _ = lib.solve_model(ds, mode, theta=sl, sim_opt=False, B_form=1,
                           workflow_family="DATA3", nfe=camp.NFE, solver_max_cpu_time=150, LOUD=False)
lin = (float(fl["parameters"]["beta_0"]), float(fl["parameters"]["beta_1"]))

# with the toggle set, fit_form's donnan solve is the 2-param (P0,X) fit at fixed k
res = camp.fit_form(ds, "donnan", base, c, lin, headline=False)
res["k_fixed"] = kfix
res["n_Bparams"] = 2
res.update({"run_id": rid, "salt": salt, "stage": None, "form": "donnan"})
out = camp.OUT / rid / "result_donnan.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(res, indent=2, default=float))
print(f"RECOVERED {rid} k={kfix} WSSE3={res.get('WSSE3')}", flush=True)
