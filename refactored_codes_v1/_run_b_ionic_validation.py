#!/usr/bin/env python3
"""Validate the B = f(ionic strength) reparameterization (Architecture.md §18.3).

PART A — single-salt RESCALE IDENTITY (deterministic, no optimizer):
  For each salt, fit B(c) (B_form=1, toggle OFF) -> theta_c=(beta_0, beta_1).
  Then forward-sim B(c) at theta_c and B(I) at theta_I=(beta_0, beta_1/k_I) and
  confirm the per-channel objectives are IDENTICAL — proving the I-reparam is a
  pure rescaling beta_1 -> beta_1/k_I (k_I=cIn power).

PART B — cross-salt CLUSTERING (the point of the change):
  beta_0 is k_I-independent; the slope beta_1 differs by salt in c-space but, if B
  is governed by ionic strength, beta_1_I = beta_1_c / k_I should be MORE
  consistent across NaCl/CaCl2/LaCl3 than beta_1_c.  We report the coefficient of
  variation of beta_1_c vs beta_1_I across salts.

Run: python3 _run_b_ionic_validation.py
"""
import sys, json
from pathlib import Path
from multiprocessing import Pool

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib  # noqa: E402

ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
OUT = ART / "b_ionic_strength"
OUT.mkdir(parents=True, exist_ok=True)

SALTS = [("NaCl", "MC2.05.07.24_NaCl"),
         ("CaCl2", "MC2.05.07.24_CaCl2"),
         ("LaCl3", "MC2.05.21.24_LaCl3")]
NFE = 80


def _fit_bform1(ds, mode, seed):
    """Robust B_form=1 fit: seeded, then a multistart retry if the seeded fit
    returns no parameters.  Returns the parameters dict or None."""
    lib.NF270_B_USE_IONIC_STRENGTH = None
    fc, _, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=False, B_form=1,
                               workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=200)
    if isinstance(fc, dict) and isinstance(fc.get("parameters"), dict):
        return dict(fc["parameters"])
    fc, _, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=False, B_form=1,
                               workflow_family="DATA3", nfe=NFE, multistart=True,
                               multistart_iterations=6, solver_max_cpu_time=120)
    if isinstance(fc, dict) and isinstance(fc.get("parameters"), dict):
        return dict(fc["parameters"])
    return None


def run(arg):
    salt, rid = arg
    try:
        lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
        fam = lib.NF270_RUN_REGISTRY[rid]
        ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
        mode = ds["mode"]
        k_I = lib.nf270_ionic_strength_factor(salt)

        # well-behaved 'single' fit -> seed + the theta for the deterministic identity
        lib.NF270_B_USE_IONIC_STRENGTH = None
        fs, _, _ = lib.solve_model(ds, mode, sim_opt=False, B_form="single",
                                   workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=150)
        ps = dict(fs["parameters"])
        base = {"Lp": ps["Lp"], "sigma": ps["sigma"], "beta_0": ps["B"]}
        for k in ("S0", "S"):
            if k in ps:
                base[k] = ps[k]

        # ---- PART A: deterministic rescale identity (no optimizer) ----
        # Pick a fixed nonzero beta_1 so B = beta_0 + beta_1*cIn truly varies, then
        # forward-sim B(c) at (beta_0, beta_1) and B(I) at (beta_0, beta_1/k_I).
        # The two must match channel-for-channel (pure rescaling).
        b1_fixed = 0.01
        theta_c = dict(base); theta_c["beta_1"] = b1_fixed
        theta_I = dict(base); theta_I["beta_1"] = b1_fixed / k_I
        lib.NF270_B_USE_IONIC_STRENGTH = None
        oc = lib._nf270_contour_objectives_at_theta(ds, theta_c, mode=mode, B_form=1,
                                                    workflow_family="DATA3", nfe=NFE)
        lib.NF270_B_USE_IONIC_STRENGTH = True
        oI = lib._nf270_contour_objectives_at_theta(ds, theta_I, mode=mode, B_form=1,
                                                    workflow_family="DATA3", nfe=NFE)
        lib.NF270_B_USE_IONIC_STRENGTH = None

        def _rel(a, b):
            if a is None or b is None or not np.isfinite(a) or not np.isfinite(b):
                return None
            return abs(a - b) / max(abs(a), 1e-12)
        identity = {
            "obj_m": [oc[0], oI[0], _rel(oc[0], oI[0])],
            "obj_cv": [oc[1], oI[1], _rel(oc[1], oI[1])],
            "obj_cr": [oc[2], oI[2], _rel(oc[2], oI[2])],
        }
        max_rel = max([v[2] for v in identity.values() if v[2] is not None], default=None)

        # ---- PART B: fitted B-slope per salt (for cross-salt clustering) ----
        seed = dict(base); seed["beta_1"] = 0.0
        theta_fit = _fit_bform1(ds, mode, seed)
        if theta_fit is not None:
            b0f, b1c = theta_fit.get("beta_0"), theta_fit.get("beta_1")
            b1I = (b1c / k_I) if b1c is not None else None
        else:
            b0f = b1c = b1I = None

        return {"salt": salt, "rid": rid, "k_I": k_I,
                "beta_0": b0f, "beta_1_c": b1c, "beta_1_I": b1I,
                "fit_ok": theta_fit is not None,
                "identity": identity, "identity_max_rel": max_rel}
    except Exception as exc:
        import traceback
        return {"salt": salt, "rid": rid, "error": repr(exc), "tb": traceback.format_exc()}


def _cv(xs):
    xs = [x for x in xs if x is not None and np.isfinite(x)]
    if len(xs) < 2:
        return None
    m = np.mean(xs)
    return float(np.std(xs) / abs(m)) if m != 0 else None


if __name__ == "__main__":
    with Pool(3) as pool:
        results = pool.map(run, SALTS)
    ok = [r for r in results if "error" not in r and r.get("fit_ok")]
    b0 = [r["beta_0"] for r in ok]
    b1c = [r["beta_1_c"] for r in ok]
    b1I = [r["beta_1_I"] for r in ok]
    summary = {
        "per_salt": results,
        "cluster": {
            "n_salts": len(ok),
            "beta_0_CV": _cv(b0),
            "beta_1_c_CV": _cv(b1c),
            "beta_1_I_CV": _cv(b1I),
        },
    }
    with open(OUT / "b_ionic_validation.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)

    print("\n=== PART A: single-salt rescale identity (B(c)@(β0,β1)  vs  B(I)@(β0,β1/kI)) ===")
    for r in results:
        if "error" in r:
            print(f"  {r['salt']:6s} ERROR: {r['error']}")
            continue
        mr = r["identity_max_rel"]
        verdict = "IDENTICAL" if (mr is not None and mr < 1e-6) else ("CHECK" if mr is not None else "n/a")
        print(f"  {r['salt']:6s} k_I={r['k_I']:.0f}  max rel-diff across channels = "
              f"{mr if mr is None else f'{mr:.2e}'}  ({verdict})")
    print("\n=== PART B: cross-salt clustering of the fitted B-slope ===")
    print(f"  {'salt':6s} {'k_I':>4s} {'beta_0':>10s} {'beta_1_c':>12s} {'beta_1_I=β1c/kI':>16s}")
    for r in results:
        if "error" in r or not r.get("fit_ok"):
            print(f"  {r['salt']:6s}  (no B_form=1 fit)")
            continue
        print(f"  {r['salt']:6s} {r['k_I']:>4.0f} {r['beta_0']:>10.4g} {r['beta_1_c']:>12.4g} {r['beta_1_I']:>16.4g}")
    c = summary["cluster"]
    print(f"\n  ({c['n_salts']} salts fitted)")
    print(f"  CV(beta_0)   = {c['beta_0_CV']}")
    print(f"  CV(beta_1_c) = {c['beta_1_c_CV']}   (concentration-space slope)")
    print(f"  CV(beta_1_I) = {c['beta_1_I_CV']}   (ionic-strength-space slope)")
    if c["beta_1_c_CV"] and c["beta_1_I_CV"]:
        verdict = "TIGHTER in I-space (ionic strength unifies the slope)" if c["beta_1_I_CV"] < c["beta_1_c_CV"] else "NOT tighter in I-space"
        print(f"  -> beta_1 clustering: {verdict}")
    print("\n[b-ionic] DONE")
