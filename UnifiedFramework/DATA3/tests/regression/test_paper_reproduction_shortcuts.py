"""Regression tests for the local DATA1/DATA2 paper reproduction shortcuts."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest
from unittest.mock import patch


os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

REPO_ROOT = Path(__file__).resolve().parents[4]
REFAC_ROOT = REPO_ROOT / "refactored_codes_v1"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REFAC_ROOT) not in sys.path:
    sys.path.insert(0, str(REFAC_ROOT))

from refactored_ucb_library import run_paper_reproduction  # noqa: E402
from refactored_ucb_runfile import _paper_workflow_mode  # noqa: E402
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode import conductivity_paper as conductivity_paper_module  # noqa: E402


@pytest.mark.regression
def test_runfile_paper_mode_menu_is_simple() -> None:
    """The runfile should keep the paper shortcut menu small and obvious."""
    with patch("builtins.input", return_value="1"):
        assert _paper_workflow_mode() == "DATA1"
    with patch("builtins.input", return_value="2"):
        assert _paper_workflow_mode() == "DATA2"


@pytest.mark.regression
def test_conductivity_paper_module_is_available() -> None:
    """The conductivity-paper helper should stay importable and numerically stable."""
    cond = conductivity_paper_module.variant_shedlovsky(
        [0.01],
        temp=298.15,
        epsilon=78.5,
        eta=0.0089,
        lambda_0=126.5,
        a=1.0e-8,
        z_1=1,
        z_2=-1,
        lambda_0_cation=50.0,
        lambda_0_anion=76.5,
    )
    assert len(cond) == 1
    assert np.isfinite(cond[0])


@pytest.mark.slow
@pytest.mark.regression
def test_data1_paper_reproduction_shortcut_writes_artifacts(tmp_path: Path) -> None:
    """DATA1 paper mode should create the local reproduction outputs and metadata."""
    result = run_paper_reproduction("DATA1", repo_root=REPO_ROOT, output_root=tmp_path / "repro")
    out_dir = Path(result["output_dir"])

    assert out_dir.exists()
    meta_path = out_dir / "run_metadata.json"
    assert meta_path.exists()

    meta = json.loads(meta_path.read_text())
    assert meta["paper_family"] == "DATA1"
    assert "DATA1_511.12" in meta["datasets"]

    stage_a = meta.get("stage_a_artifacts", {})
    stage_b = meta.get("stage_b_artifacts", {})
    assert stage_a, "DATA1 stage A artifacts were not generated."
    assert stage_b, "DATA1 stage B artifacts were not generated."

    for artifact_path in stage_a.get("figures", []):
        assert Path(artifact_path).exists(), f"Missing DATA1 figure: {artifact_path}"
    for artifact_path in stage_a.get("tables", []):
        assert Path(artifact_path).exists(), f"Missing DATA1 table: {artifact_path}"

    if "manifest" in stage_b:
        assert Path(stage_b["manifest"]).exists(), f"Missing DATA1 manifest: {stage_b['manifest']}"


@pytest.mark.slow
@pytest.mark.regression
def test_data2_paper_reproduction_shortcut_writes_artifacts(tmp_path: Path) -> None:
    """DATA2 paper mode should create the local reproduction outputs and metadata."""
    result = run_paper_reproduction("DATA2", repo_root=REPO_ROOT, output_root=tmp_path / "repro")
    out_dir = Path(result["output_dir"])

    assert out_dir.exists()
    meta_path = out_dir / "run_metadata.json"
    assert meta_path.exists()

    meta = json.loads(meta_path.read_text())
    assert meta["paper_family"] == "DATA2"
    assert "DATA2_270511.123" in meta["datasets"]
    assert "DATA2_270611.123" in meta["datasets"]

    data2_artifacts = meta.get("data2_artifacts", {})
    assert data2_artifacts, "DATA2 artifacts were not generated."
    for dataset_id, artifact_group in data2_artifacts.items():
        for artifact_path in artifact_group.get("tables", []):
            assert Path(artifact_path).exists(), f"Missing DATA2 table for {dataset_id}: {artifact_path}"
        for artifact_path in artifact_group.get("figures", []):
            assert Path(artifact_path).exists(), f"Missing DATA2 figure for {dataset_id}: {artifact_path}"
