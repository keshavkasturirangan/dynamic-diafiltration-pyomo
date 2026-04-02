"""Regression coverage for the object-oriented unified workflow module."""

from __future__ import annotations

import os
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

REPO_ROOT = Path(__file__).resolve().parents[4]
UNICODE_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "ExperimentalDataAnalysis" / "UnifiedCode"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(UNICODE_ROOT) not in sys.path:
    sys.path.insert(0, str(UNICODE_ROOT))

from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode import unified_codebase_runfile as runfile  # noqa: E402
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library import (  # noqa: E402
    ExperimentMode,
    ModelOptions,
    ParameterTreatmentMode,
    ProcessModelProfile,
)
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.workflow.unified_workflow import (  # noqa: E402
    DataLoader,
    DatasetRequest,
    UQEngine,
)
from UnifiedFramework.DATA3.scripts.reproduce import reproduce_data1_data2 as reproduction_script  # noqa: E402


def _run_reproduction(tmp_path: Path, *datasets: str) -> Path:
    """Run the reproduction script for one or more datasets and return the latest output dir."""
    out_root = tmp_path / "repro"
    cmd = [
        sys.executable,
        str(REPO_ROOT / "UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py"),
        "--repo-root",
        str(REPO_ROOT),
        "--out-root",
        str(out_root),
        "--nfe",
        "10",
        "--skip-curve-metrics",
        "--datasets",
        *datasets,
    ]
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)
    run_dirs = sorted(out_root.glob("*"))
    assert run_dirs
    return run_dirs[-1]


@pytest.mark.regression
def test_object_oriented_loader_supports_mat_and_xlsx_requests() -> None:
    """The new DataLoader should treat MAT and XLSX as first-class requests."""
    loader = DataLoader()
    requests = [
        DatasetRequest(file_path=runfile.PRESET_DATA1_PATHS[0], dataset_id="DATA1_501_1"),
        DatasetRequest(
            file_path=runfile.PRESET_DATA3_PATH,
            selector=runfile.PRESET_DATA3_SELECTOR,
            dataset_id="DATA3_MC2",
        ),
    ]

    datasets, validations = loader.load_many(requests)

    assert len(datasets) == 2
    assert validations[0][0] is True
    assert validations[1][0] is True
    assert datasets[0].source.value == "mat"
    assert datasets[1].filename.endswith(".xlsx")
    assert datasets[1].sheet_name == runfile.PRESET_DATA3_SELECTOR
    assert len(datasets[1].vials) >= 1


@pytest.mark.regression
def test_uqengine_compare_models_returns_information_criteria() -> None:
    """The OOP UQEngine should expose AIC-style model comparison results."""
    loader = DataLoader()
    dataset, validation = loader.load(
        DatasetRequest(file_path=runfile.PRESET_DATA1_PATHS[0], dataset_id="DATA1_501_1")
    )
    assert validation[0] is True

    engine = UQEngine(data_loader=loader)
    comparison = engine.compare_models(
        [dataset],
        [
            ModelOptions(
                mode=ExperimentMode.DATA,
                parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
                b_form="single",
                nfe=10,
                use_sigma_logit_transform=False,
                process_model_profile=ProcessModelProfile.DATA1,
                paper_profile="DATA1_PAPER",
            )
        ],
    )

    assert not comparison.empty
    assert {"candidate_name", "objective", "aic", "aicc", "bic", "delta_aic"}.issubset(comparison.columns)
    assert np.isfinite(float(comparison.loc[0, "objective"]))


@pytest.mark.regression
def test_reproduction_run_dataset_routes_through_uqengine() -> None:
    """The DATA1/DATA2 reproduction path should use UQEngine as its main orchestrator."""
    engine = UQEngine(data_loader=DataLoader())
    spec = reproduction_script.CANONICAL_DATASETS[0]

    row, exp, options, estimation = reproduction_script.run_dataset(
        engine=engine,
        repo_root=REPO_ROOT,
        spec=spec,
        solver="ipopt",
        nfe=10,
        run_doe=False,
        calc_cov=False,
        skip_curve_metrics=True,
        ipopt_max_cpu_time=None,
    )

    assert row["workflow_class"] == "UQEngine"
    assert row["dataset_id"] == spec.dataset_id
    assert exp.dataset_id == spec.dataset_id
    assert options.process_model_profile == spec.process_model_profile
    assert isinstance(estimation, dict)


@pytest.mark.regression
def test_reproduction_script_writes_data2_artifact_metadata(tmp_path: Path) -> None:
    """A DATA2 reproduction smoke run should emit workflow-owned DATA2 artifact metadata."""
    latest = _run_reproduction(tmp_path, "DATA2_270611.123")
    meta = json.loads((latest / "run_metadata.json").read_text())
    assert "DATA2_270611.123" in meta.get("data2_artifacts", {})
    artifacts = meta["data2_artifacts"]["DATA2_270611.123"]
    assert artifacts["tables"]
    for artifact_path in artifacts["tables"]:
        assert Path(artifact_path).exists(), f"Missing DATA2 artifact: {artifact_path}"


@pytest.mark.regression
def test_reproduction_script_writes_data1_stage_artifact_files(tmp_path: Path) -> None:
    """A DATA1 reproduction smoke run should emit concrete Stage A/B artifact files."""
    latest = _run_reproduction(tmp_path, "DATA1_511.12")
    meta = json.loads((latest / "run_metadata.json").read_text())

    stage_a = meta.get("stage_a_artifacts", {})
    stage_b = meta.get("stage_b_artifacts", {})
    assert stage_a, "Missing DATA1 Stage A artifact metadata."
    assert stage_b, "Missing DATA1 Stage B artifact metadata."

    for artifact_path in stage_a.get("figures", []):
        assert Path(artifact_path).exists(), f"Missing DATA1 Stage A figure: {artifact_path}"
    for artifact_path in stage_a.get("tables", []):
        assert Path(artifact_path).exists(), f"Missing DATA1 Stage A table: {artifact_path}"

    manifest_path = Path(stage_b["manifest"])
    assert manifest_path.exists(), f"Missing DATA1 Stage B manifest: {manifest_path}"
    manifest = np.genfromtxt(manifest_path, delimiter=",", dtype=str, skip_header=1)
    assert manifest.size > 0


@pytest.mark.regression
def test_reproduction_script_writes_data2_artifact_files_with_expected_kinds(tmp_path: Path) -> None:
    """DATA2 workflow artifacts should include both measurement-comparison and curve-metric evidence."""
    latest = _run_reproduction(tmp_path, "DATA2_270611.123")
    meta = json.loads((latest / "run_metadata.json").read_text())
    artifacts = meta["data2_artifacts"]["DATA2_270611.123"]
    table_paths = [Path(path) for path in artifacts["tables"]]
    assert any("fig3_measurement" in path.name for path in table_paths)
    assert any("fig6_curve_metrics" in path.name for path in table_paths)
    for path in table_paths:
        assert path.exists(), f"Missing DATA2 artifact file: {path}"
