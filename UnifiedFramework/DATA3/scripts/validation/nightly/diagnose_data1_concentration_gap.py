#!/usr/bin/env python3
"""Diagnose the DATA1 concentration-gap between unified fits and legacy MATLAB baselines."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
from scipy.io import loadmat

REPO_ROOT = Path(__file__).resolve().parents[5]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from UnifiedFramework.DATA3.tests.regression.test_unified_codebase_pytest_validation import (
    _compare_series,
    _extract_unified_data1_predictions,
    _run_data1_pipeline,
    runfile,
)
SIM_ROOT = (
    REPO_ROOT
    / "UnifiedFramework"
    / "DATA3"
    / "docs"
    / "validation"
    / "simulation_validation_data_files"
    / "data1_main"
    / "fig2"
)
DATA1_ROOT = REPO_ROOT / "DATA1_matlab" / "data"
OUT_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "pytest_validation"


def _load_legacy_fit_params(dataset: str, fit_dir: str) -> dict[str, float]:
    fit_stru = loadmat(DATA1_ROOT / fit_dir / "fit_stru.mat", squeeze_me=True, struct_as_record=False)["fit_stru"]
    return {
        "Lp": float(fit_stru.Lp),
        "B": float(fit_stru.B),
        "sigma": float(fit_stru.sigma),
    }


def _endpoint_report(unified_df: pd.DataFrame, baseline_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for vial in sorted(unified_df["vial"].unique()):
        uni = unified_df[unified_df["vial"] == vial].sort_values("time_min")
        ret = baseline_df[baseline_df["series_id"] == f"retentate_prediction_vial_{vial}"].sort_values("time_min")
        perm = baseline_df[baseline_df["series_id"] == f"permeate_prediction_vial_{vial}"].sort_values("time_min")
        vialp = baseline_df[baseline_df["series_id"] == f"vial_prediction_vial_{vial}"].sort_values("time_min")
        rows.append(
            {
                "vial": int(vial),
                "unified_cF_end_mM": float(uni["retentate_concentration_mM"].iloc[-1]),
                "baseline_cF_end_mM": float(ret["concentration_mM"].iloc[-1]),
                "delta_cF_end_mM": float(uni["retentate_concentration_mM"].iloc[-1] - ret["concentration_mM"].iloc[-1]),
                "unified_cH_end_mM": float(uni["permeate_concentration_mM"].iloc[-1]),
                "baseline_cH_end_mM": float(perm["concentration_mM"].iloc[-1]),
                "delta_cH_end_mM": float(uni["permeate_concentration_mM"].iloc[-1] - perm["concentration_mM"].iloc[-1]),
                "unified_cV_end_mM": float(uni["vial_concentration_mM"].iloc[-1]),
                "baseline_cV_end_mM": float(vialp["concentration_mM"].iloc[-1]),
                "delta_cV_end_mM": float(uni["vial_concentration_mM"].iloc[-1] - vialp["concentration_mM"].iloc[-1]),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    dataset_specs = [
        ("501.1", runfile.PRESET_DATA1_PATHS[0], "data1_main_fig2_b.csv", "501.1 concpolar"),
        ("511.12", runfile.PRESET_DATA1_PATHS[3], "data1_main_fig2_d.csv", "511.12 concpolar"),
    ]

    summary_rows: list[dict[str, object]] = []
    endpoint_frames: list[pd.DataFrame] = []
    series_frames: list[pd.DataFrame] = []

    for dataset, dataset_path, baseline_name, fit_dir in dataset_specs:
        result = _run_data1_pipeline(dataset_path)
        unified = _extract_unified_data1_predictions(result)
        baseline = pd.read_csv(SIM_ROOT / baseline_name)

        unified_theta = {k: float(v) for k, v in result["parmest"]["theta"].to_dict().items()}
        legacy_theta = _load_legacy_fit_params(dataset, fit_dir)

        endpoint_df = _endpoint_report(unified, baseline)
        endpoint_df.insert(0, "dataset", dataset)
        endpoint_frames.append(endpoint_df)

        series_rows: list[dict[str, object]] = []
        for matcher, ucol in [
            ("retentate_prediction", "retentate_concentration_mM"),
            ("permeate_prediction", "permeate_concentration_mM"),
            ("vial_prediction", "vial_concentration_mM"),
        ]:
            part = baseline[baseline["series_id"].str.startswith(matcher)].copy()
            compared = pd.DataFrame(
                _compare_series(
                    part,
                    unified,
                    panel_id=f"{dataset}-{matcher}",
                    quantity_col="concentration_mM",
                    unified_col=ucol,
                )
            )
            compared.insert(0, "dataset", dataset)
            compared.insert(1, "matcher", matcher)
            series_rows.append(compared)
            summary_rows.append(
                {
                    "dataset": dataset,
                    "matcher": matcher,
                    "max_mae": float(compared["mae"].max()),
                    "max_abs_err": float(compared["max_abs_err"].max()),
                    "unified_Lp": unified_theta["Lp"],
                    "legacy_Lp": legacy_theta["Lp"],
                    "unified_B": unified_theta["B"],
                    "legacy_B": legacy_theta["B"],
                    "unified_sigma": unified_theta["sigma"],
                    "legacy_sigma": legacy_theta["sigma"],
                }
            )
        series_frames.append(pd.concat(series_rows, ignore_index=True))

    pd.DataFrame(summary_rows).to_csv(OUT_ROOT / "data1_concentration_gap_summary.csv", index=False)
    pd.concat(endpoint_frames, ignore_index=True).to_csv(OUT_ROOT / "data1_concentration_gap_endpoints.csv", index=False)
    pd.concat(series_frames, ignore_index=True).to_csv(OUT_ROOT / "data1_concentration_gap_series.csv", index=False)

    note = {
        "observation": (
            "The DATA1 concentration gap is driven by unified parameter estimates that differ materially "
            "from the legacy MATLAB fit, especially sigma for dataset 501.1."
        ),
        "artifacts": [
            "UnifiedFramework/DATA3/results/pytest_validation/data1_concentration_gap_summary.csv",
            "UnifiedFramework/DATA3/results/pytest_validation/data1_concentration_gap_endpoints.csv",
            "UnifiedFramework/DATA3/results/pytest_validation/data1_concentration_gap_series.csv",
        ],
    }
    (OUT_ROOT / "data1_concentration_gap_note.json").write_text(json.dumps(note, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
