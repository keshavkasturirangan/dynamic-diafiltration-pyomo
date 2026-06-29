#!/usr/bin/env python3
"""DATA3 §8 — mass & concentration fits: individual fit vs pooled estimate.

For each poorly-identified rep, forward-sim the constant-B fit with (a) its OWN individual
theta and (b) the POOLED estimate (group-shared σ, B from the stacking; own Lp), and render
the OFFICIAL mass / permeate-conc / retentate-conc time series via the library's
run_data3_time_series_plots (NO custom plotting — same call as _make_predictions.py).

Honest read: the pooled (more identifiable) estimate reproduces the data essentially as well
as the individual best-fit → pooling improves identifiability at ~no fit cost.

Output: bform_study/fit_improvement/<rid>/{individual,pooled}/{mass,concentration}-*.png
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json
from pathlib import Path
import matplotlib; matplotlib.use("Agg")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib

ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
STUDY = ART / "bform_study"
OUT = STUDY / "fit_improvement"
lib.NF270_CF_RESIDUAL_SCALE_FRACTION = 0.02     # canonical 2% no-floor
lib.NF270_CF_RESIDUAL_FLOOR_MM = None

# poorly-identified reps → their salt×regime group (diluting NaCl excluded: well-identified)
REP_GROUP = {
    "MC2.05.07.24_NaCl": "NaCl_concentrating",
    "MC2.05.07.24_CaCl2": "CaCl2",
    "MC2.05.21.24_LaCl3": "LaCl3",
}


def _ds(rid):
    fam = lib.NF270_RUN_REGISTRY[rid]
    return lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]


def _render(rid, theta, label):
    ds = _ds(rid)
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        _, sim, _ = lib.solve_model(ds, ds["mode"], theta=dict(theta), sim_opt=True,
                                    B_form="single", workflow_family="DATA3", nfe=150,
                                    solver_max_cpu_time=150, LOUD=False)
    results = {
        "model_settings": {"workflow_family": "DATA3"},
        "data": ds, "sim_stru": sim,
        "fit_meta": {"recipe_name": label,
                     "winning_theta": {k: v for k, v in theta.items() if isinstance(v, (int, float))}},
    }
    saved = lib.run_data3_time_series_plots(results, save_dir=OUT / rid / label, show=False)
    return len(saved)


def main():
    for rid, group in REP_GROUP.items():
        indiv = json.loads((STUDY / rid / "result_single.json").read_text())["parameters"]
        pooled_min = json.loads((STUDY / "Bsigma_fixedLp_nofloor" / "_stack" / group / "single" / "stacked.json").read_text())["pooled_min"]
        pooled = dict(indiv)
        pooled["sigma"] = float(pooled_min["sigma"])
        pooled["B"] = float(pooled_min["B"])     # single-form: B@feed == constant B
        n1 = _render(rid, indiv, "individual")
        n2 = _render(rid, pooled, "pooled")
        print(f"[fit-improve] {rid} ({group}): individual σ={indiv['sigma']:.2f} B={indiv['B']:.2f} | "
              f"pooled σ={pooled['sigma']:.2f} B={pooled['B']:.2f}  -> {n1}+{n2} figs", flush=True)


if __name__ == "__main__":
    main()
