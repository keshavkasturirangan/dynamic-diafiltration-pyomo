#!/usr/bin/env python3
"""Profile-likelihood in the Donnan fixed-charge scale X for a single DATA3 sheet.

The campaign Donnan fits rail X to its lower bound (1e-3 mM); X->0 collapses the
Donnan B(c) to the constant Jw*P0*k (a FLAT curve).  This script pins X at a grid
of values (rail -> values large enough to bend the curve inside the measured
window) and re-optimizes the rest (P0, sigma, S, Lp; k free unless k_fixed given),
recording the objective cost of forcing the bend.  It NEVER overwrites the
canonical result_donnan.json — it writes only Xprofile.json.

Verdict logic (printed): if forcing X off the rail to a curve-bending value costs
only a small fraction of WSSE3, X is unidentified (degeneracy); if it costs a
large fraction, the data genuinely prefers X->0 (flat B).

Usage: python3 _profile_donnan_X.py <run_id> <salt> [k_fixed] [X1,X2,...]
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

camp.NFE = 80
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0

rid = sys.argv[1]
salt = sys.argv[2] if len(sys.argv) > 2 else ""
k_fixed = float(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] not in ("", "free", "None") else None
if len(sys.argv) > 4:
    XGRID = [float(x) for x in sys.argv[4].split(",")]
else:
    XGRID = [1e-3, 2.0, 5.0, 10.0, 20.0, 40.0]


def fit_at_X(ds, mode, base, c, lin, xval):
    """One Donnan fit with X pinned at xval (k fixed iff k_fixed given)."""
    lib.NF270_DONNAN_FIX_X = float(xval)
    lib.NF270_DONNAN_FIX_K = float(k_fixed) if k_fixed is not None else None
    try:
        res = camp.fit_form(ds, "donnan", base, c, lin, headline=False)
    except Exception as e:
        return {"X": xval, "error": repr(e)}
    finally:
        lib.NF270_DONNAN_FIX_X = None
        lib.NF270_DONNAN_FIX_K = None
    p = res.get("parameters") or {}
    return {
        "X": xval,
        "P0": p.get("P0"), "k_dd": p.get("k_dd"), "sigma": p.get("sigma"),
        "Lp": p.get("Lp"), "S": p.get("S"),
        "WSSE3": res.get("WSSE3"), "Obj": res.get("Obj"),
        "obj_m": res.get("obj_m"), "obj_cv": res.get("obj_cv"),
        "obj_cr": res.get("obj_cr"),
    }


def main():
    ds = camp.load(rid)
    mode = ds["mode"]
    # single-form warm start (same anchor as the donnan worker)
    fs, _, _ = lib.solve_model(ds, mode, sim_opt=False, B_form="single",
                               workflow_family="DATA3", nfe=camp.NFE,
                               solver_max_cpu_time=200, LOUD=False)
    base = dict(fs["parameters"])
    c, bp, bs = camp.apparent_partition(ds, base)
    sl = {k: base[k] for k in ("Lp", "sigma") if k in base}
    for k in ("S0", "S"):
        if k in base:
            sl[k] = base[k]
    sl["beta_0"] = base.get("B", 1.0); sl["beta_1"] = 0.0
    fl, _, _ = lib.solve_model(ds, mode, theta=sl, sim_opt=False, B_form=1,
                               workflow_family="DATA3", nfe=camp.NFE,
                               solver_max_cpu_time=150, LOUD=False)
    lin = (float(fl["parameters"]["beta_0"]), float(fl["parameters"]["beta_1"]))

    rows = []
    for xv in XGRID:
        r = fit_at_X(ds, mode, base, c, lin, xv)
        rows.append(r)
        print(f"  X={xv:<8.4g} -> WSSE3={r.get('WSSE3')}  P0={r.get('P0')}  "
              f"k={r.get('k_dd')}  obj_m={r.get('obj_m')} obj_cv={r.get('obj_cv')} "
              f"obj_cr={r.get('obj_cr')}", flush=True)

    ok = [r for r in rows if r.get("WSSE3") is not None]
    base_w = min((r["WSSE3"] for r in ok), default=None)
    out = {
        "run_id": rid, "salt": salt, "k_fixed": k_fixed,
        "c_window_mM": [float(c.min()), float(c.max())],
        "Xgrid": XGRID, "rows": rows,
        "WSSE3_min": base_w,
    }
    odir = camp.OUT / rid
    odir.mkdir(parents=True, exist_ok=True)
    (odir / "Xprofile.json").write_text(json.dumps(out, indent=2, default=float))
    print(f"\nsaved {odir/'Xprofile.json'}")
    if base_w:
        print(f"\n=== {rid} ({salt}) X-profile  (WSSE3_min={base_w:.4g}) ===")
        for r in ok:
            dw = 100.0 * (r["WSSE3"] - base_w) / base_w
            print(f"  X={r['X']:<8.4g}  WSSE3={r['WSSE3']:<12.5g}  "
                  f"(+{dw:6.2f}% vs min)   P0={r['P0']:.4g} k={r['k_dd']:.4g}")


if __name__ == "__main__":
    main()
