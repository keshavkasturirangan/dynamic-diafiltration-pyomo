#!/usr/bin/env python3
"""DATA3 NF270 B-form campaign — fits every B(c) correlation on every single-salt
sheet and records theta, WSSE-by-channel, model-selection stats, FIM uncertainty,
and the apparent-B(c) curve.  Writes ONLY to paper_artifacts/nf270/bform_study/ —
the legacy constant / B_form=1 / B_form=2 artifacts are never touched.

Forms:  single (constant), 1 (linear), 2 (quadratic), 3 (cubic),
        sat (saturating exponential), donnan (Donnan-dielectric).

Rigor:  headline GOOD/OKAY/POOR sheets get LHS multistart (N=10) + FIM covariance;
        the rest get a single warm-started solve.  New forms (sat/donnan) are
        seeded by an offline scipy fit of the closed form to the single-fit
        apparent partition B_p(c) = Js/(Jw*(cIn-cH)), so IPOPT starts near the
        basin instead of cold.

Usage:  python3 _run_bform_campaign.py [sheets=all|<rid,rid>] [forms=all|<f,f>]
Resumable: a (sheet,form) whose result.json exists is skipped.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, math, traceback
from pathlib import Path
from multiprocessing import Pool
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib

try:
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
OUT = ART / "bform_study"
OUT.mkdir(parents=True, exist_ok=True)

NFE = 150
N_MULTISTART = 10
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0  # campaign cfg (NEXT_ACTIONS)

# 11 single-salt sheets: (run_id, salt, stage-label or None).  Headline = labelled.
SHEETS = [
    ("MC3.07.22.24_SNaCl",  "NaCl",  "GOOD"),
    ("MC5.07.23.24_S2NaCl", "NaCl",  None),
    ("MC4.07.11.24_SNaCl",  "NaCl",  None),
    ("MC2.05.07.24_NaCl",   "NaCl",  "OKAY"),
    ("MC5.07.23.24_NaCl",   "NaCl",  None),
    ("MC5.07.23.24_SNaCl",  "NaCl",  None),
    ("MC2.05.07.24_CaCl2",  "CaCl2", "OKAY"),
    ("MC3.07.11.24_SCaCl2", "CaCl2", None),
    ("MC3.07.12.24_S2CaCl2","CaCl2", None),
    ("MC2.05.21.24_LaCl3",  "LaCl3", "POOR"),
    ("MC4.07.11.24_SLaCl3", "LaCl3", None),
]
FORMS = ["single", 1, 2, 3, "sat", "donnan"]
NB = {"single": 1, 1: 2, 2: 3, 3: 4, "sat": 2, "donnan": 3}  # # B-shape params


def load(rid):
    fam = lib.NF270_RUN_REGISTRY[rid]
    return lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]


def apparent_partition(ds, theta_single):
    """Forward-sim at the single fit; return (c[mM], B_partition) over interior vials.
    B_partition = Js/(Jw*(cIn-cH)) is the Jw-scaled partition the polynomial / sat /
    donnan forms parameterize (Js = Jw*B*(cIn-cH))."""
    _, sim, _ = lib.solve_model(ds, ds["mode"], theta=theta_single, sim_opt=True,
                                B_form="single", workflow_family="DATA3", nfe=NFE,
                                solver_max_cpu_time=120, LOUD=False)
    nv0 = int(ds["data_config"].get("n_v0", 1))
    cs, bp, bs = [], [], []
    for k, v in sim.items():
        if (k + 1) < nv0:
            continue
        cIn = np.asarray(v["cIn"], float); cH = np.asarray(v["cH"], float)
        Jw = np.asarray(v["Jw"], float); Js = np.asarray(v["Js"], float)
        dc = cIn - cH
        m = (dc > 1e-6) & (Jw > 1e-9)
        if not m.any():
            continue
        Bp = Js[m] / (Jw[m] * dc[m])            # partition (Jw-scaled)
        Bmu = Js[m] * 10000.0 / dc[m]           # apparent µm/s
        cs.append(float(np.median(cIn[m]))); bp.append(float(np.median(Bp))); bs.append(float(np.median(Bmu)))
    o = np.argsort(cs)
    return np.array(cs)[o], np.array(bp)[o], np.array(bs)[o]


def _lin_curve(lin, c):
    """Linear partition B(c)=beta_0+beta_1*c over the data c-grid (partition units)."""
    b0, b1 = lin
    cg = np.asarray(c, float)
    if cg.size < 3:
        cg = np.linspace(1.0, 100.0, 12)
    y = np.maximum(b0 + b1 * cg, 1e-6)
    return cg, y


def seed_sat(c, lin):
    """Seed B_inf,c_star by fitting the saturating form to the converged linear
    partition curve (so the amplitude is in the model's native partition units)."""
    cg, y = _lin_curve(lin, c)
    Binf = max(1e-3, float(np.nanmax(y)) * 1.3)
    cstar = max(1.0, float(np.nanmedian(cg)))
    if HAVE_SCIPY:
        try:
            f = lambda x, Bi, cs: Bi * (1 - np.exp(-x / cs))
            p, _ = curve_fit(f, cg, y, p0=[Binf, cstar],
                             bounds=([1e-4, 1.0], [50.0, 500.0]), maxfev=20000)
            Binf, cstar = float(p[0]), float(p[1])
        except Exception:
            pass
    return {"B_inf": Binf, "c_star": cstar}


def seed_donnan(c, lin):
    """Seed P0,X,k by fitting the Donnan-dielectric form to the converged linear
    partition curve (units-consistent, near-basin start for IPOPT)."""
    cg, y = _lin_curve(lin, c)
    P0 = max(1e-3, float(np.nanmax(y)) * 1.6)
    X, k = 8.0, 0.5
    if HAVE_SCIPY:
        try:
            f = lambda x, P0, X, k: P0 * (-X + np.sqrt(X**2 + 4 * k**2 * x**2)) / (2 * x)
            p, _ = curve_fit(f, cg, y, p0=[P0, X, k],
                             bounds=([1e-4, 1e-2, 1e-2], [100.0, 200.0, 5.0]), maxfev=20000)
            P0, X, k = float(p[0]), float(p[1]), float(p[2])
        except Exception:
            pass
    return {"P0": P0, "X": X, "k_dd": k}


def fit_form(ds, form, base_single, c, lin, headline):
    """Return result dict for one (sheet, form).  `lin`=(beta_0,beta_1) from the
    converged linear fit, used to seed the mechanistic forms in partition units."""
    mode = ds["mode"]
    seed = {kk: base_single[kk] for kk in ("Lp", "sigma") if kk in base_single}
    for kk in ("S0", "S"):
        if kk in base_single:
            seed[kk] = base_single[kk]
    cpu = 200
    if form == "single":
        seed["B"] = base_single.get("B", 1.0)
    elif isinstance(form, int):
        # warm-start polynomials from the converged LINEAR fit (partition units),
        # not from the raw µm/s B — the latter is ~10x off and fails to converge.
        seed["beta_0"] = lin[0]
        seed["beta_1"] = lin[1]
        seed["beta_2"] = 0.0
        seed["beta_3"] = 0.0
    elif form == "sat":
        seed.update(seed_sat(c, lin)); cpu = 300
    elif form == "donnan":
        seed.update(seed_donnan(c, lin)); cpu = 500

    # multistart only for the well-conditioned forms; quad/cubic (ill-scaled c^k
    # term) and donnan (3 correlated params) are single-shot from their data seed.
    ms = bool(headline) and form in ("single", 1, "sat")
    fit, sim, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=False, B_form=form,
                                  workflow_family="DATA3", nfe=NFE,
                                  multistart=ms, multistart_iterations=N_MULTISTART,
                                  solver_max_cpu_time=cpu, LOUD=False)
    # NB: donnan fix-k recovery is done out-of-band by _recover_donnan.py with a
    # process-level timeout — an inline multi-k loop is unsafe because a failing
    # fix-k value can grind for 10-18 min in the uncapped casadi Simulator.
    k_fixed = None
    if not isinstance(fit, dict):
        raise RuntimeError("solver returned no feasible fit")
    wsse3 = sum(float(fit.get(kk, 0.0) or 0.0) for kk in ("obj_m", "obj_cv", "obj_cr"))
    res = {
        "parameters": fit["parameters"],
        "WSSE3": wsse3,                      # 3-channel weighted SSE (deck metric)
        "Obj": float(fit.get("Obj")),        # full objective (carries truncate term)
        "obj_m": float(fit.get("obj_m")), "obj_cv": float(fit.get("obj_cv")),
        "obj_cr": float(fit.get("obj_cr")),
        "n_Bparams": NB[form] if k_fixed is None else NB[form] - 1,
    }
    if k_fixed is not None:
        res["k_fixed"] = k_fixed
    rs = fit.get("res_std", {})
    n_data = sum(len(rs.get(kk, [])) for kk in ("res_m", "res_cp", "res_cf"))
    res["n_data"] = int(n_data)
    # AICc / BIC from the weighted 3-channel WSSE
    k = NB[form] + 2  # + Lp + sigma  (S is a nuisance state)
    sse = max(wsse3, 1e-12)
    if n_data > k + 1:
        res["AIC"] = n_data * math.log(sse / n_data) + 2 * k
        res["AICc"] = res["AIC"] + 2 * k * (k + 1) / (n_data - k - 1)
        res["BIC"] = n_data * math.log(sse / n_data) + k * math.log(n_data)
    if "multistart" in fit:
        res["multistart"] = {kk: fit["multistart"].get(kk) for kk in ("iterations", "n_trials", "best_score")}
    # FIM uncertainty (+ sensitivity conditioning) on headline sheets
    if headline:
        try:
            doe = lib.calc_FIM(ds, mode, theta=fit["parameters"], B_form=form,
                               workflow_family="DATA3", nfe=NFE)
            ev = [float(x) for x in (doe.get("eig_val") or [])]
            res["FIM"] = {
                "std": [float(x) for x in (doe.get("std") or [])],
                "eig_min": min(ev) if ev else None,
                "eig_max": max(ev) if ev else None,
                "cond": (max(ev) / min(ev) if ev and min(ev) > 0 else None),
                "det": float(doe.get("det")) if doe.get("det") is not None else None,
            }
        except Exception as e:
            res["FIM_error"] = str(e)
    return res


def process_sheet(task):
    """Fit every requested form on one sheet; write per-fit result files.
    Returns the list of summary rows for this sheet."""
    rid, salt, stage, forms = task
    rows = []
    headline = stage is not None
    try:
        ds = load(rid)
    except Exception as e:
        print(f"[LOAD FAIL] {rid}: {e}", flush=True)
        return rows
    sdir = OUT / rid
    sdir.mkdir(parents=True, exist_ok=True)
    # base single fit + apparent partition for seeding
    if True:
        try:
            # seed fit only (single-shot); the rigorous multistart 'single' result
            # is produced by the form loop below.
            fs, _, _ = lib.solve_model(ds, ds["mode"], sim_opt=False, B_form="single",
                                       workflow_family="DATA3", nfe=NFE,
                                       solver_max_cpu_time=200, LOUD=False)
            base_single = dict(fs["parameters"])
            c, bp, bs = apparent_partition(ds, base_single)
            np.savetxt(sdir / "apparent_B.csv",
                       np.column_stack([c, bp, bs]), delimiter=",",
                       header="c_mM,B_partition,B_apparent_ums", comments="")
            # linear fit (partition units) to seed the mechanistic forms consistently
            try:
                sl = {kk: base_single[kk] for kk in ("Lp", "sigma") if kk in base_single}
                for kk in ("S0", "S"):
                    if kk in base_single:
                        sl[kk] = base_single[kk]
                sl["beta_0"] = base_single.get("B", 1.0); sl["beta_1"] = 0.0
                fl, _, _ = lib.solve_model(ds, ds["mode"], theta=sl, sim_opt=False, B_form=1,
                                           workflow_family="DATA3", nfe=NFE,
                                           solver_max_cpu_time=150, LOUD=False)
                lin = (float(fl["parameters"]["beta_0"]), float(fl["parameters"]["beta_1"]))
            except Exception:
                lin = (float(base_single.get("B", 1.0)), 0.0)
        except Exception as e:
            print(f"[BASE FAIL] {rid}: {e}", flush=True)
            traceback.print_exc()
            return rows

        for form in forms:
            tag = form if isinstance(form, str) else f"poly{form}"
            rpath = sdir / f"result_{tag}.json"
            if rpath.exists():
                print(f"[skip] {rid} {tag} (exists)", flush=True)
                res = json.loads(rpath.read_text())
            else:
                try:
                    res = fit_form(ds, form, base_single, c, lin, headline)
                    res.update({"run_id": rid, "salt": salt, "stage": stage, "form": str(form)})
                    rpath.write_text(json.dumps(res, indent=2, default=float))
                    p = res["parameters"]
                    extra = {kk: round(p[kk], 4) for kk in p if kk not in ("S0", "S")}
                    print(f"[ok] {rid:22s} {tag:7s} WSSE3={res.get('WSSE3', float('nan')):.2f} "
                          f"AICc={res.get('AICc', float('nan')):.1f} {extra}", flush=True)
                except Exception as e:
                    print(f"[FIT FAIL] {rid} {tag}: {e}", flush=True)
                    traceback.print_exc()
                    res = {"run_id": rid, "salt": salt, "stage": stage, "form": str(form),
                           "error": str(e)}
                    rpath.write_text(json.dumps(res, indent=2, default=float))
            row = {"run_id": rid, "salt": salt, "stage": stage or "", "form": str(form),
                   "WSSE3": res.get("WSSE3"), "Obj": res.get("Obj"), "obj_m": res.get("obj_m"),
                   "obj_cv": res.get("obj_cv"), "obj_cr": res.get("obj_cr"),
                   "n_Bparams": res.get("n_Bparams"), "n_data": res.get("n_data"),
                   "AICc": res.get("AICc"), "BIC": res.get("BIC"),
                   "FIM_cond": (res.get("FIM") or {}).get("cond")}
            rows.append(row)
    return rows


def main():
    argsheets = sys.argv[1] if len(sys.argv) > 1 else "all"
    argforms = sys.argv[2] if len(sys.argv) > 2 else "all"
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    sheets = SHEETS if argsheets == "all" else [s for s in SHEETS if s[0] in argsheets.split(",")]
    forms = FORMS if argforms == "all" else [
        (int(f) if f.isdigit() else f) for f in argforms.split(",")]
    tasks = [(rid, salt, stage, forms) for rid, salt, stage in sheets]
    print(f"[campaign] {len(tasks)} sheets x {len(forms)} forms, workers={workers}", flush=True)

    all_rows = []
    if workers <= 1 or len(tasks) == 1:
        for t in tasks:
            all_rows += process_sheet(t)
            _write_summary(all_rows)
    else:
        with Pool(min(workers, len(tasks))) as pool:
            for rows in pool.imap_unordered(process_sheet, tasks):
                all_rows += rows
                _write_summary(all_rows)
    _write_summary(all_rows)
    print("CAMPAIGN DONE", flush=True)


def _write_summary(rows):
    import csv
    cols = ["run_id", "salt", "stage", "form", "WSSE3", "Obj", "obj_m", "obj_cv", "obj_cr",
            "n_Bparams", "n_data", "AICc", "BIC", "FIM_cond"]
    with open(OUT / "_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in cols})


if __name__ == "__main__":
    main()
