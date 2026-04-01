"""Pytest coverage for published-paper simulation validation figures."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "UnifiedFramework"
    / "DATA3"
    / "scripts"
    / "validation"
    / "nightly"
    / "validate_simulation_validation_data_files.py"
)
MANIFEST_CSV = (
    REPO_ROOT
    / "UnifiedFramework"
    / "DATA3"
    / "docs"
    / "validation"
    / "simulation_validation_data_files"
    / "simulation_validation_manifest.csv"
)
REPORT_CSV = (
    REPO_ROOT
    / "UnifiedFramework"
    / "DATA3"
    / "results"
    / "simulation_validation"
    / "simulation_validation_figure_report.csv"
)
SIM_VALIDATION_ROOT = MANIFEST_CSV.parent


MANIFEST_DF = pd.read_csv(MANIFEST_CSV) if MANIFEST_CSV.exists() else pd.DataFrame()
MANIFEST_CASES = [
    pytest.param(
        row,
        id=f"{row.paper}-{row.figure_id}-{row.panel_id}-{Path(row.relative_csv).stem}",
    )
    for row in MANIFEST_DF.itertuples(index=False)
]
REPORT_TARGET_IDS = sorted(pd.read_csv(REPORT_CSV)["target_id"].tolist()) if REPORT_CSV.exists() else []


def _load_validation_module():
    spec = importlib.util.spec_from_file_location("simulation_validation_figures", SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load validation script from {SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _assert_manifest_asset_shape(row, asset_path: Path) -> None:
    if row.artifact_kind == "error_note":
        payload = asset_path.read_text().strip()
        assert payload, f"Expected non-empty error note: {asset_path}"
        return

    df = pd.read_csv(asset_path)
    assert not df.empty, f"Expected non-empty artifact CSV: {asset_path}"

    if row.artifact_kind == "timeseries":
        assert len(df.columns) >= 3, f"Expected rich timeseries columns in {asset_path}"
        assert any(
            col in df.columns
            for col in ("series_id", "time_min", "vial_index", "cIn_mM", "interface_concentration_mM")
        )
        return

    if row.artifact_kind == "contour_grid":
        assert {"objective_value", "is_optimum"}.issubset(df.columns), asset_path
        assert bool(df["is_optimum"].astype(bool).any()), f"Expected flagged optimum in {asset_path}"
        return

    if row.artifact_kind == "scatter":
        assert len(df.columns) >= 3, f"Expected scatter coordinates in {asset_path}"
        return

    if row.artifact_kind in {"parameter_sweep", "bar_values", "residual_distribution"}:
        assert len(df.columns) >= 2, f"Expected structured data columns in {asset_path}"
        return

    raise AssertionError(f"Unhandled artifact_kind={row.artifact_kind} for {asset_path}")


@pytest.mark.regression
def test_simulation_validation_manifest_exists() -> None:
    """One-time generated simulation validation manifest should exist in the committed baseline folder."""
    assert MANIFEST_CSV.exists(), f"Missing simulation validation manifest: {MANIFEST_CSV}"
    df = pd.read_csv(MANIFEST_CSV)
    assert not df.empty
    assert {"paper", "figure_id", "relative_csv"}.issubset(df.columns)


@pytest.mark.regression
@pytest.mark.parametrize("row", MANIFEST_CASES)
def test_simulation_validation_manifest_row_has_valid_artifact(row) -> None:
    """Each committed panel/plot artifact should get its own pass/fail pytest outcome."""
    asset_path = SIM_VALIDATION_ROOT / row.relative_csv
    assert asset_path.exists(), f"Missing simulation validation artifact: {asset_path}"
    _assert_manifest_asset_shape(row, asset_path)


@pytest.fixture(scope="module")
def simulation_validation_report() -> pd.DataFrame:
    """Run the composite paper validation once and share the report across per-figure tests."""
    module = _load_validation_module()
    report = module.run_validation()
    assert REPORT_CSV.exists(), f"Expected comparison report at {REPORT_CSV}"
    assert not report.empty
    return report


@pytest.mark.regression
@pytest.mark.slow
@pytest.mark.nightly
@pytest.mark.parametrize("target_id", REPORT_TARGET_IDS)
def test_simulation_validation_figure_target_matches_paper_extract(
    simulation_validation_report: pd.DataFrame,
    target_id: str,
) -> None:
    """Each paper-figure composite gets its own pass/fail pytest outcome."""
    row = simulation_validation_report.loc[simulation_validation_report["target_id"] == target_id]
    assert not row.empty, f"Missing figure-validation row for {target_id}"
    status = row.iloc[0]["status"]
    assert status == "PASS", (
        f"Figure validation mismatch for {target_id}: "
        f"status={status}, best_score={row.iloc[0]['best_score']}, best_layout={row.iloc[0]['best_layout']}"
    )
