#!/usr/bin/env python3
"""DATA3 prediction-vs-measured plots per B-form.  For each fitted (sheet, form)
in bform_study/, forward-sim at the fitted theta and render the official mass /
permeate-conc / retentate-conc time series via run_data3_time_series_plots
(startup vial excluded by n_v0; built-in legend; [mM] axis).

Output: bform_study/predictions/<run_id>/<form>/{mass,concentration}-*.png
Usage:  python3 _make_predictions.py [sheets=all|<rid,..>] [forms=all|<f,..>]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib

ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
STUDY = ART / "bform_study"
PRED = STUDY / "predictions"
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0

FORMS = ["single", 1, 2, 3, "sat", "donnan"]
SHEET_IDS = [
    "MC3.07.22.24_SNaCl", "MC5.07.23.24_S2NaCl", "MC4.07.11.24_SNaCl",
    "MC2.05.07.24_NaCl", "MC5.07.23.24_NaCl", "MC5.07.23.24_SNaCl",
    "MC2.05.07.24_CaCl2", "MC3.07.11.24_SCaCl2", "MC3.07.12.24_S2CaCl2",
    "MC2.05.21.24_LaCl3", "MC4.07.11.24_SLaCl3",
]


def load_ds(rid):
    fam = lib.NF270_RUN_REGISTRY[rid]
    return lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]


def one(rid, form):
    tag = form if isinstance(form, str) else f"poly{form}"
    rp = STUDY / rid / f"result_{tag}.json"
    if not rp.exists():
        return f"{rid}/{tag}: no result"
    r = json.loads(rp.read_text())
    if "parameters" not in r or "error" in r:
        return f"{rid}/{tag}: failed fit"
    ds = load_ds(rid)
    try:
        _, sim, _ = lib.solve_model(ds, ds["mode"], theta=r["parameters"], sim_opt=True,
                                    B_form=form, workflow_family="DATA3", nfe=150,
                                    solver_max_cpu_time=150, LOUD=False)
    except Exception as e:
        return f"{rid}/{tag}: sim fail {e}"
    p = r["parameters"]
    results = {
        "model_settings": {"workflow_family": "DATA3"},
        "data": ds,
        "sim_stru": sim,
        "fit_meta": {"b_form": str(form),
                     "winning_theta": {k: v for k, v in p.items() if isinstance(v, (int, float))},
                     "recipe_name": f"B_form={form}"},
    }
    out = PRED / rid / tag
    saved = lib.run_data3_time_series_plots(results, save_dir=out, show=False)
    return f"{rid}/{tag}: {len(saved)} figs"


def main():
    sel = sys.argv[1] if len(sys.argv) > 1 else "all"
    fsel = sys.argv[2] if len(sys.argv) > 2 else "all"
    sheets = SHEET_IDS if sel == "all" else sel.split(",")
    forms = FORMS if fsel == "all" else [(int(f) if f.isdigit() else f) for f in fsel.split(",")]
    for rid in sheets:
        for form in forms:
            try:
                print("[pred]", one(rid, form), flush=True)
            except Exception as e:
                print(f"[pred FAIL] {rid} {form}: {e}", flush=True)
    print("PREDICTIONS DONE", flush=True)


if __name__ == "__main__":
    main()
