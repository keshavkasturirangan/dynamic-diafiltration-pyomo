#!/usr/bin/env python3
"""Exercise the formal solve_model_per_concentration_band wrapper + per-band FIM
on the GOOD/OKAY/POOR sheets, and report per-band theta + FIM non-singularity.

This is the formalized, FIM-guarded version of the value test.
Run: python3 _run_band_wrapper.py
"""
import sys, json, traceback
from pathlib import Path
from multiprocessing import Pool

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib  # noqa: E402

ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
OUT = ART / "band_value_test"
OUT.mkdir(parents=True, exist_ok=True)

SHEETS = [("GOOD", "MC3.07.22.24_SNaCl"),
          ("OKAY", "MC2.05.07.24_NaCl"),
          ("POOR", "MC2.05.21.24_LaCl3")]


def run(arg):
    role, rid = arg
    try:
        lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
        fam = lib.NF270_RUN_REGISTRY[rid]
        ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
        res = lib.solve_model_per_concentration_band(
            ds, ds["mode"], n_bands=2, B_form="single",
            workflow_family="DATA3", nfe=80, compute_fim=True, solver_max_cpu_time=200)
        res["role"] = role
        res["rid"] = rid
        return res
    except Exception as exc:
        return {"role": role, "rid": rid, "error": repr(exc), "tb": traceback.format_exc()}


if __name__ == "__main__":
    with Pool(3) as pool:
        results = pool.map(run, SHEETS)
    with open(OUT / "band_wrapper_summary.json", "w") as fh:
        json.dump(results, fh, indent=2, default=float)
    print("\n=== solve_model_per_concentration_band + per-band FIM ===")
    for r in results:
        if "error" in r:
            print(f"  {r['role']:4s} {r['rid']}: ERROR {r['error']}")
            continue
        f = r["full"]["theta"]
        print(f"\n  {r['role']:4s} {r['rid']}  n_fitted={r['n_fitted']}  bands={r['bands']}")
        if r["warnings"]:
            for w in r["warnings"]:
                print(f"     ! {w}")
        print(f"     full : σ={f.get('sigma'):.3f} Lp={f.get('Lp'):.2f} B={f.get('B'):.2f}")
        for e in r["per_band"]:
            t = e["theta"]; fim = e.get("fim", {})
            cf = e["cf_range"]
            fim_s = (f"min_eig={fim.get('min_eig'):.3g} cond={fim.get('cond'):.3g} singular={fim.get('singular')}"
                     if "min_eig" in fim else f"FIM:{fim}")
            print(f"     band{e['band_index']} cF[{cf[0]:.0f}-{cf[1]:.0f}]mM vials={e['band_vials']}")
            print(f"           σ={t.get('sigma'):.3f} Lp={t.get('Lp'):.2f} B={t.get('B'):.2f}  WSSE_cr={e.get('obj_cr'):.3g}")
            print(f"           {fim_s}")
    print("\n[band-wrapper] DONE")
