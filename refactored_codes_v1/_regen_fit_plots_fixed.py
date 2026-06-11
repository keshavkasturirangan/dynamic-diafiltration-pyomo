#!/usr/bin/env python3
"""Regenerate the per-sheet mass/concentration fit plots CORRECTLY.

Bug being fixed: the ad-hoc consistent_initial_guess plots drew the leading
startup/hold-up vial's predicted mV, so the mass trace showed an initial rise
with no corresponding collected-permeate measurement.  The sanctioned plotter
run_data3_time_series_plots() skips the startup vial's prediction
(`(i+1) >= n_v0`), matching the DATA1/DATA2 convention.

For each single-salt sheet we forward-simulate at the contour-consistent θ
(from summary_contour_consistent_theta.csv) and render via the official plotter.
Output: meeting_prep_2026-06-10/fits_fixed/<prefix>_{mass,concentration}*.png
"""
import csv
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import os
os.chdir(str(ROOT.parent))
import refactored_ucb_library as lib

REPO = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo")
SUMMARY = REPO / "UnifiedFramework/DATA3/results/paper_artifacts/nf270/meeting_prep_2026-06-10/consistent_initial_guess/summary_contour_consistent_theta.csv"
OUT = REPO / "UnifiedFramework/DATA3/results/paper_artifacts/nf270/meeting_prep_2026-06-10/fits_fixed"
OUT.mkdir(parents=True, exist_ok=True)

reg = getattr(lib, "NF270_RUN_REGISTRY", {})
print(f"registry keys: {len(reg)}")

WARM = REPO / "UnifiedFramework/DATA3/results/paper_artifacts/nf270/warm_start_fits/summary.json"
warm_by_sheet = {}
for e in json.loads(WARM.read_text()):
    warm_by_sheet[e["run_id"]] = dict(e["warm_start_theta_fit"])
print(f"warm-start sheets: {len(warm_by_sheet)}")

theta_by_sheet = {}
with open(SUMMARY) as f:
    for r in csv.DictReader(f):
        theta_by_sheet[r["sheet"]] = r

results = {}
for sheet, r in theta_by_sheet.items():
    try:
        fam = reg.get(sheet)
        if fam is None:
            results[sheet] = {"error": "not in NF270_RUN_REGISTRY"}
            print(f"SKIP {sheet}: not in registry")
            continue
        wb = lib.NF270_DEFAULT_ROOT / fam["workbook"]
        ds = lib.loadxlsx(str(wb), sheet=fam["sheet"])["data_stru"]
        # contour-consistent θ = warm-start θ with Lp, σ overridden from the contour argmin
        theta = dict(warm_by_sheet.get(sheet, {}))
        theta["Lp"] = float(r["Lp_contour"])
        theta["sigma"] = float(r["sigma_contour"])
        theta.setdefault("B", float(r["B_warm"]))
        theta.setdefault("S0", 0.0)
        fit_stru, sim_stru, _ = lib.solve_model(
            ds, ds["mode"], theta, sim_opt=True, B_form="single",
            LOUD=False, workflow_family="DATA3", multistart=False,
            solver_max_cpu_time=120)
        res = {"data": ds, "sim_stru": sim_stru,
               "model_settings": {"workflow_family": "DATA3"}}
        paths = lib.run_data3_time_series_plots(res, save_dir=str(OUT / sheet), show=False)
        results[sheet] = {"ok": True, "n_v0": ds["data_config"].get("n_v0"),
                          "plots": [str(p) for p in paths]}
        print(f"OK {sheet}: n_v0={ds['data_config'].get('n_v0')}  {len(paths)} plots")
    except Exception as exc:
        results[sheet] = {"error": f"{type(exc).__name__}: {exc}"}
        print(f"FAIL {sheet}: {exc}")
        traceback.print_exc()

(OUT / "_regen_index.json").write_text(json.dumps(results, indent=2, default=str))
print(f"\nDONE. index -> {OUT/'_regen_index.json'}")
