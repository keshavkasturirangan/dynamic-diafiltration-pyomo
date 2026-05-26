#!/usr/bin/env python3
"""
Legacy DATA2-style workflow, pointed at DATA3 / NF270 Excel sheets.

The goal is deliberately narrow:
    - load the Excel sheets through the refactored loader
    - convert conductivity to concentration before the legacy fit
    - fit with the old utility.py solver using B_form = 1
    - save mass-vs-time and concentration-vs-time plots

This is the compatibility bridge we can use to answer:
    "What changes if we force the DATA2 workflow onto the DATA3 campaign?"
"""
from __future__ import annotations

import json
import os
import sys
import time
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(HERE.parent) not in sys.path:
    sys.path.insert(0, str(HERE.parent))

import utility as legacy  # noqa: E402
import refactored_ucb_library as lib  # noqa: E402


NF270_ROOT = Path(
    os.environ.get(
        "DIAFILTRATION_NF270_ROOT",
        REPO_ROOT / "UnifiedFramework" / "ExperimentalDataFiles",
    )
).expanduser().resolve()

SINGLE_SALT_RUN_IDS = (
    "MC2.05.07.24_NaCl",
    "MC3.07.22.24_SNaCl",
    "MC4.07.11.24_SNaCl",
    "MC5.07.23.24_NaCl",
    "MC5.07.23.24_SNaCl",
    "MC5.07.23.24_S2NaCl",
    "MC2.05.07.24_CaCl2",
    "MC3.07.11.24_SCaCl2",
    "MC3.07.12.24_S2CaCl2",
    "MC2.05.21.24_LaCl3",
    "MC4.07.11.24_SLaCl3",
)

OUTPUT_ROOT = HERE / "outputs"
FIGURE_ROOT = OUTPUT_ROOT / "figures"
SUMMARY_PATH = OUTPUT_ROOT / "data2_workflow_for_data3_summary.json"
COMPARISON_CSV = OUTPUT_ROOT / "data2_workflow_for_data3_comparison.csv"
COMPARISON_JSON = OUTPUT_ROOT / "data2_workflow_for_data3_comparison.json"
B_FORM = 1
REFERENCE_SUMMARY_CANDIDATES = (
    REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270" / "regime_gated_fits" / "summary.json",
    REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270" / "regime_gated_fits_v2" / "summary.json",
    REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270" / "warm_start_fits" / "summary.json",
    REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270" / "warm_start_fits_with_perm" / "summary.json",
)


def _load_run(run_id: str):
    family = lib.NF270_RUN_REGISTRY[run_id]
    wb_path = NF270_ROOT / family["workbook"]
    if not wb_path.exists():
        raise FileNotFoundError(f"missing workbook: {wb_path}")
    return legacy.loadxlsx(wb_path, sheet=family["sheet"])["data_stru"]


def _param_snapshot(fit_stru):
    params = {}
    for key, value in (fit_stru.get("parameters", {}) or {}).items():
        if isinstance(value, (int, float)):
            params[key] = float(value)
        elif isinstance(value, dict):
            params[key] = {
                str(k): float(v)
                for k, v in value.items()
                if isinstance(v, (int, float))
            }
    return params


def _first_last_finite(value):
    arr = np.asarray(value, dtype=float).reshape(-1)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan"), float("nan")
    return float(arr[0]), float(arr[-1])


def _load_reference_summary():
    override = os.environ.get("DATA3_REFERENCE_SUMMARY")
    if override:
        path = Path(override).expanduser().resolve()
        if path.exists():
            return path, json.loads(path.read_text(encoding="utf-8"))
    for path in REFERENCE_SUMMARY_CANDIDATES:
        if path.exists():
            return path, json.loads(path.read_text(encoding="utf-8"))
    return None, None


def _compare_to_reference(legacy_summary: list[dict]) -> dict:
    ref_path, ref_summary = _load_reference_summary()
    report = {
        "reference_path": str(ref_path) if ref_path else None,
        "rows": [],
        "matched": 0,
        "unmatched": 0,
    }
    if not ref_summary:
        return report

    ref_index = {
        str(row.get("run_id")): row
        for row in ref_summary
        if isinstance(row, dict) and row.get("run_id")
    }
    for row in legacy_summary:
        if not isinstance(row, dict) or row.get("error"):
            continue
        run_id = str(row.get("run_id"))
        ref = ref_index.get(run_id)
        if ref is None:
            report["unmatched"] += 1
            continue

        legacy_params = row.get("parameters", {}) or {}
        ref_theta = ref.get("theta_fit", {}) or {}
        legacy_obj = row.get("obj", float("nan"))
        ref_obj = ref.get("total_obj", float("nan"))
        legacy_lp = legacy_params.get("Lp", float("nan"))
        ref_lp = ref_theta.get("Lp", float("nan"))
        legacy_sigma = legacy_params.get("sigma", float("nan"))
        ref_sigma = ref_theta.get("sigma", float("nan"))
        legacy_beta_0 = legacy_params.get("beta_0", float("nan"))
        legacy_beta_1 = legacy_params.get("beta_1", float("nan"))

        comparison_row = {
            "run_id": run_id,
            "reference_path": str(ref_path),
            "legacy_obj": float(legacy_obj),
            "reference_obj": float(ref_obj),
            "delta_obj": float(legacy_obj - ref_obj),
            "legacy_Lp": float(legacy_lp),
            "reference_Lp": float(ref_lp),
            "delta_Lp": float(legacy_lp - ref_lp),
            "legacy_sigma": float(legacy_sigma),
            "reference_sigma": float(ref_sigma),
            "delta_sigma": float(legacy_sigma - ref_sigma),
            "legacy_beta_0": float(legacy_beta_0),
            "legacy_beta_1": float(legacy_beta_1),
            "reference_B": float(ref_theta.get("B", float("nan"))),
            "legacy_cF_first_mM": float(row.get("cF_first_mM", float("nan"))),
            "reference_C_F0_mM": float((ref.get("mass_balance") or {}).get("C_F0_mM", float("nan"))),
            "legacy_cF_last_mM": float(row.get("cF_last_mM", float("nan"))),
            "reference_cF_final_icp_mM": float((ref.get("mass_balance") or {}).get("cF_final_icp_mM", float("nan"))),
            "legacy_fit_time_s": float(row.get("fit_time_s", float("nan"))),
            "reference_wall_time_s": float(ref.get("wall_time_s", float("nan"))),
        }
        report["rows"].append(comparison_row)
        report["matched"] += 1

    if report["rows"]:
        COMPARISON_CSV.parent.mkdir(parents=True, exist_ok=True)
        with COMPARISON_CSV.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(report["rows"][0].keys()))
            writer.writeheader()
            writer.writerows(report["rows"])
        COMPARISON_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def run_one(run_id: str) -> dict:
    data_stru = _load_run(run_id)
    t0 = time.time()
    fit_stru, sim_stru, sim_inter = legacy.solve_model(
        data_stru,
        data_stru["mode"],
        theta=None,
        sim_opt=False,
        B_form=B_FORM,
        LOUD=False,
    )
    fit_time_s = time.time() - t0
    legacy.plot_sim_comparison(
        data_stru,
        sim_stru,
        plot_pred=True,
        lg=False,
        LOUD=True,
    )
    return {
        "run_id": run_id,
        "workbook": lib.NF270_RUN_REGISTRY[run_id]["workbook"],
        "sheet": lib.NF270_RUN_REGISTRY[run_id]["sheet"],
        "mode": data_stru.get("mode"),
        "dataset": data_stru.get("dataset"),
        "sheet_name": data_stru.get("sheet_name"),
        "cF_first_mM": _first_last_finite(data_stru["data_raw"][0]["cF_exp"])[0],
        "cF_last_mM": _first_last_finite(data_stru["data_raw"][-1]["cF_exp"])[1],
        "fit_time_s": round(fit_time_s, 3),
        "obj": float(fit_stru.get("Obj", float("nan"))),
        "obj_m": float(fit_stru.get("obj_m", float("nan"))),
        "obj_cv": float(fit_stru.get("obj_cv", float("nan"))),
        "obj_cr": float(fit_stru.get("obj_cr", float("nan"))),
        "parameters": _param_snapshot(fit_stru),
    }


def run_all() -> list[dict]:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    os.chdir(OUTPUT_ROOT)

    print(f"[DATA2-for-DATA3] workbook root: {NF270_ROOT}")
    print(f"[DATA2-for-DATA3] B_form        : {B_FORM}")
    print(f"[DATA2-for-DATA3] output root   : {OUTPUT_ROOT}")
    print()

    summary = []
    for run_id in SINGLE_SALT_RUN_IDS:
        print(f"=== {run_id} ===")
        try:
            record = run_one(run_id)
        except Exception as exc:
            print(f"  ERROR: {type(exc).__name__}: {exc}")
            record = {
                "run_id": run_id,
                "error": f"{type(exc).__name__}: {exc}",
            }
        summary.append(record)
        if "error" not in record:
            params = record.get("parameters", {})
            print(
                f"  ok  fit={record['fit_time_s']:.1f}s  "
                f"obj={record['obj']:.3g}  "
                f"Lp={params.get('Lp', float('nan')):.3g}  "
                f"sigma={params.get('sigma', float('nan')):.3g}"
            )
        print()

    SUMMARY_PATH.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[DATA2-for-DATA3] wrote summary: {SUMMARY_PATH}")
    comparison = _compare_to_reference(summary)
    if comparison.get("reference_path"):
        print(f"[DATA2-for-DATA3] comparison reference: {comparison['reference_path']}")
        print(f"[DATA2-for-DATA3] comparison rows   : {comparison['matched']}")
        print(f"[DATA2-for-DATA3] comparison CSV    : {COMPARISON_CSV}")
    print(f"[DATA2-for-DATA3] figures live in: {FIGURE_ROOT}")
    return summary


if __name__ == "__main__":
    run_all()
