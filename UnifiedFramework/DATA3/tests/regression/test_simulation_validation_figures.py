"""Pytest coverage for published-paper simulation validation figures."""

from __future__ import annotations

import importlib.util
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


def _load_validation_module():
    spec = importlib.util.spec_from_file_location("simulation_validation_figures", SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load validation script from {SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.regression
def test_simulation_validation_manifest_exists() -> None:
    """One-time generated simulation validation manifest should exist in the committed baseline folder."""
    assert MANIFEST_CSV.exists(), f"Missing simulation validation manifest: {MANIFEST_CSV}"
    df = pd.read_csv(MANIFEST_CSV)
    assert not df.empty
    assert {"paper", "figure_id", "relative_csv"}.issubset(df.columns)


@pytest.mark.regression
@pytest.mark.slow
@pytest.mark.nightly
def test_simulation_validation_figures_match_paper_extracts() -> None:
    """Composite legacy figure artifacts should continue to match the extracted published-paper figures."""
    module = _load_validation_module()
    report = module.run_validation()
    assert REPORT_CSV.exists(), f"Expected comparison report at {REPORT_CSV}"
    assert not report.empty
    failing = report.loc[report["status"] != "PASS", ["target_id", "status", "best_score", "best_layout"]]
    assert failing.empty, "Figure validation mismatches:\n" + failing.to_string(index=False)
