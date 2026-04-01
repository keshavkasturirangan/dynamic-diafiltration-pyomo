#!/usr/bin/env python3
"""Compare unified DATA1 estimated outputs directly against legacy fit_stru.mat outputs."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat

REPO_ROOT = Path(__file__).resolve().parents[5]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from UnifiedFramework.DATA3.tests.regression.test_unified_codebase_pytest_validation import (  # noqa: E402
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


def _flatten_legacy_fit(dataset: str, fit_dir: str) -> pd.DataFrame:
    data_stru = loadmat(
        DATA1_ROOT / f"data_stru-dataset{dataset}.mat",
        squeeze_me=True,
        struct_as_record=False,
    )["data_stru"]
    fit_stru = loadmat(
        DATA1_ROOT / fit_dir / "fit_stru.mat",
        squeeze_me=True,
        struct_as_record=False,
    )["fit_stru"]
    t_delay = float(np.atleast_1d(data_stru.data_raw[0].time)[0])

    rows: list[dict[str, float | str]] = []
    for vial_idx, sim in enumerate(np.atleast_1d(fit_stru.sim_stru), start=1):
        time_min = (np.atleast_1d(sim.time).astype(float) - t_delay) / 60.0
        for t_val, y_val in zip(time_min, np.atleast_1d(sim.cF).astype(float)):
            rows.append(
                {
                    "series_id": f"retentate_prediction_vial_{vial_idx}",
                    "dataset": dataset,
                    "time_min": float(t_val),
                    "concentration_mM": float(y_val),
                }
            )
        for t_val, y_val in zip(time_min, np.atleast_1d(sim.cH).astype(float)):
            rows.append(
                {
                    "series_id": f"permeate_prediction_vial_{vial_idx}",
                    "dataset": dataset,
                    "time_min": float(t_val),
                    "concentration_mM": float(y_val),
                }
            )
        rows.append(
            {
                "series_id": f"vial_prediction_vial_{vial_idx}",
                "dataset": dataset,
                "time_min": float(time_min[-1]),
                "concentration_mM": float(np.atleast_1d(sim.cV).astype(float)[-1]),
            }
        )
    return pd.DataFrame(rows).sort_values(["series_id", "time_min"]).reset_index(drop=True)


def _compare_flat_exports(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    a = a.copy()
    b = b.copy()
    a["dataset"] = a["dataset"].astype(str)
    b["dataset"] = b["dataset"].astype(str)
    merged = a.merge(
        b,
        on=["series_id", "dataset", "time_min"],
        suffixes=("_direct", "_csv"),
        how="outer",
        indicator=True,
    )
    merged["delta_concentration_mM"] = merged["concentration_mM_direct"] - merged["concentration_mM_csv"]
    grouped = []
    for series_id, part in merged.groupby("series_id", dropna=False):
        deltas = part["delta_concentration_mM"].dropna().to_numpy(dtype=float)
        grouped.append(
            {
                "series_id": series_id,
                "rows_direct": int((part["_merge"] != "right_only").sum()),
                "rows_csv": int((part["_merge"] != "left_only").sum()),
                "all_rows_matched": bool((part["_merge"] == "both").all()),
                "max_abs_err": float(np.max(np.abs(deltas))) if len(deltas) else np.nan,
                "mae": float(np.mean(np.abs(deltas))) if len(deltas) else np.nan,
            }
        )
    return pd.DataFrame(grouped).sort_values("series_id").reset_index(drop=True)


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    dataset_specs = [
        ("501.1", runfile.PRESET_DATA1_PATHS[0], "data1_main_fig2_b.csv", "501.1 concpolar"),
        ("511.12", runfile.PRESET_DATA1_PATHS[3], "data1_main_fig2_d.csv", "511.12 concpolar"),
    ]

    summary_rows: list[dict[str, object]] = []
    series_frames: list[pd.DataFrame] = []
    export_equivalence_frames: list[pd.DataFrame] = []

    for dataset, dataset_path, baseline_name, fit_dir in dataset_specs:
        unified = _extract_unified_data1_predictions(_run_data1_pipeline(dataset_path))
        legacy_direct = _flatten_legacy_fit(dataset, fit_dir)
        legacy_csv = pd.read_csv(SIM_ROOT / baseline_name)

        direct_prediction = legacy_direct[legacy_direct["series_id"].str.contains("prediction")].copy()
        csv_prediction = legacy_csv[legacy_csv["series_id"].str.contains("prediction")].copy()
        equivalence = _compare_flat_exports(direct_prediction, csv_prediction)
        equivalence.insert(0, "dataset", dataset)
        export_equivalence_frames.append(equivalence)

        series_rows: list[dict[str, object]] = []
        for matcher, ucol in [
            ("retentate_prediction", "retentate_concentration_mM"),
            ("permeate_prediction", "permeate_concentration_mM"),
            ("vial_prediction", "vial_concentration_mM"),
        ]:
            part = legacy_direct[legacy_direct["series_id"].str.startswith(matcher)].copy()
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
                }
            )
        series_frames.append(pd.concat(series_rows, ignore_index=True))

    pd.DataFrame(summary_rows).to_csv(OUT_ROOT / "data1_legacy_fit_vs_unified_estimate_summary.csv", index=False)
    pd.concat(series_frames, ignore_index=True).to_csv(
        OUT_ROOT / "data1_legacy_fit_vs_unified_estimate_series.csv", index=False
    )
    pd.concat(export_equivalence_frames, ignore_index=True).to_csv(
        OUT_ROOT / "data1_legacy_fit_export_equivalence.csv", index=False
    )

    note = {
        "observation": (
            "Legacy fit_stru.mat predictions are exported faithfully into the simulation_validation_data_files "
            "DATA1 Figure 2 concentration CSVs; the remaining DATA1 concentration gap is between the unified "
            "estimated output and the legacy fit output itself."
        ),
        "artifacts": [
            "UnifiedFramework/DATA3/results/pytest_validation/data1_legacy_fit_vs_unified_estimate_summary.csv",
            "UnifiedFramework/DATA3/results/pytest_validation/data1_legacy_fit_vs_unified_estimate_series.csv",
            "UnifiedFramework/DATA3/results/pytest_validation/data1_legacy_fit_export_equivalence.csv",
        ],
    }
    (OUT_ROOT / "data1_legacy_fit_vs_unified_estimate_note.json").write_text(
        json.dumps(note, indent=2),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
