import sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import _run_bform_campaign as camp
import refactored_ucb_library as lib
camp.NFE = 80
RID = sys.argv[1] if len(sys.argv) > 1 else "MC3.07.22.24_SNaCl"
ds = camp.load(RID); mode = ds["mode"]
t0 = time.time()
fs, _, _ = lib.solve_model(ds, mode, sim_opt=False, B_form="single", workflow_family="DATA3",
                           nfe=80, solver_max_cpu_time=120, LOUD=False)
base = dict(fs["parameters"])
print(f"[single] {time.time()-t0:.0f}s  B={base['B']:.3f} Lp={base['Lp']:.3f} sig={base['sigma']:.3f}", flush=True)
c, bp, bs = camp.apparent_partition(ds, base)
sl = {k: base[k] for k in ("Lp", "sigma") if k in base}
for k in ("S0", "S"):
    if k in base:
        sl[k] = base[k]
sl["beta_0"] = base.get("B", 1.0); sl["beta_1"] = 0.0
fl, _, _ = lib.solve_model(ds, mode, theta=sl, sim_opt=False, B_form=1, workflow_family="DATA3",
                           nfe=80, solver_max_cpu_time=120, LOUD=False)
lin = (float(fl["parameters"]["beta_0"]), float(fl["parameters"]["beta_1"]))
print(f"[linear] beta=({lin[0]:.4g},{lin[1]:.4g})", flush=True)
print("[seed sat]", camp.seed_sat(c, lin), flush=True)
print("[seed donnan]", camp.seed_donnan(c, lin), flush=True)
for form in ("2", "3", "sat", "donnan"):
    fv = int(form) if form.isdigit() else form
    t1 = time.time()
    try:
        res = camp.fit_form(ds, fv, base, c, lin, headline=False)
        pp = {k: round(v, 4) for k, v in res["parameters"].items() if k not in ("S0", "S")}
        print(f"[{form}] {time.time()-t1:.0f}s  WSSE3={res.get('WSSE3', float('nan')):.3f}  "
              f"AICc={res.get('AICc', float('nan')):.1f}  {pp}", flush=True)
    except Exception as e:
        print(f"[{form}] FAILED after {time.time()-t1:.0f}s: {e}", flush=True)
print("TEST DONE", flush=True)
