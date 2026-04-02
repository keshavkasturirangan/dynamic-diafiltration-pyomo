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
    ParameterScope,
    UQEngine,
)
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.workflow import unified_workflow as workflow_module  # noqa: E402
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
    """The OOP UQEngine should expose AIC-style comparison for one dataset + many model options."""
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
            ),
            ModelOptions(
                mode=ExperimentMode.DATA,
                parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
                b_form="pervial",
                nfe=10,
                use_sigma_logit_transform=False,
                process_model_profile=ProcessModelProfile.DATA1,
                paper_profile="DATA1_PAPER",
            ),
        ],
        candidate_names=["single_b", "per_vial_b"],
    )

    assert not comparison.empty
    assert {"candidate_name", "objective", "aic", "aicc", "bic", "delta_aic"}.issubset(comparison.columns)
    assert set(comparison["candidate_name"]) == {"single_b", "per_vial_b"}
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


@pytest.mark.regression
def test_uqengine_builds_explicit_multi_dataset_regression_plan() -> None:
    """The engine should expose shared vs dataset-specific regression intent explicitly."""
    loader = DataLoader()
    datasets, validations = loader.load_many(
        [
            DatasetRequest(file_path=runfile.PRESET_DATA2_PATHS[0], dataset_id="DATA2_A"),
            DatasetRequest(file_path=runfile.PRESET_DATA2_PATHS[1], dataset_id="DATA2_B"),
        ]
    )
    assert all(ok for ok, _ in validations)

    engine = UQEngine(data_loader=loader)
    plan = engine.build_regression_plan(
        datasets,
        ParameterScope(
            shared_parameters=("Lp", "sigma", "beta_0"),
            dataset_specific_parameters=("beta_1",),
        ),
    )

    assert plan.mode == "joint_shared_then_dataset_specific_refinement"
    assert plan.dataset_ids == ("DATA2_A", "DATA2_B")
    assert "Lp" in plan.shared_parameters
    assert "beta_1" in plan.dataset_specific_parameters


@pytest.mark.regression
def test_uqengine_multi_dataset_fit_returns_dataset_summaries(monkeypatch: pytest.MonkeyPatch) -> None:
    """A many-datasets + one-model workflow should report both global and per-dataset outputs."""
    loader = DataLoader()
    datasets, validations = loader.load_many(
        [
            DatasetRequest(file_path=runfile.PRESET_DATA2_PATHS[0], dataset_id="DATA2_A"),
            DatasetRequest(file_path=runfile.PRESET_DATA2_PATHS[1], dataset_id="DATA2_B"),
        ]
    )
    assert all(ok for ok, _ in validations)

    engine = UQEngine(data_loader=loader)
    options = ModelOptions(
        mode=ExperimentMode.LAG,
        parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
        b_form="convection",
        nfe=10,
        use_sigma_logit_transform=False,
        fix_sigma_in_estimation=True,
        beta_0_lb=0.5,
        beta_0_ub=5.0,
        beta_1_lb=0.0,
        beta_1_ub=0.1,
        process_model_profile=ProcessModelProfile.DATA2,
        paper_profile="DATA2_PAPER",
    )

    def _fake_estimate_parameters_with_parmest_v24(
        *,
        experiments,
        options,
        guess=None,
        calc_cov=True,
        cov_n=None,
        solver="ipopt",
        solver_options=None,
        tee=False,
    ):
        if len(experiments) > 1:
            return {
                "objective": 12.34,
                "theta": {"Lp": 1.2, "sigma": 0.9, "beta_0": 2.1, "beta_1": 0.03},
            }
        dataset_id = str(experiments[0].dataset_id)
        objective = 5.0 if dataset_id == "DATA2_A" else 7.34
        return {
            "objective": objective,
            "theta": {"Lp": 1.2, "sigma": 0.9, "beta_0": 2.1, "beta_1": 0.03},
        }

    monkeypatch.setattr(
        workflow_module,
        "estimate_parameters_with_parmest_v24",
        _fake_estimate_parameters_with_parmest_v24,
    )

    result = engine.fit_parameters(
        datasets,
        options,
        parameter_scope=ParameterScope(shared_parameters=("Lp", "sigma", "beta_0", "beta_1")),
        calc_cov=False,
        solver="ipopt",
        solver_options={"max_iter": 5000, "tol": 1e-6, "acceptable_tol": 1e-5},
        tee=False,
    )

    assert result["parameter_scope"]["mode"] == "joint_shared"
    assert len(result["dataset_summaries"]) == 2
    assert {row["dataset_id"] for row in result["dataset_summaries"]} == {"DATA2_A", "DATA2_B"}
    assert all("objective" in row for row in result["dataset_summaries"])
    assert "regression_plan" in result


@pytest.mark.regression
def test_uqengine_run_full_workflow_reports_many_dataset_many_model_metadata() -> None:
    """The full workflow should describe the many-datasets + many-models case explicitly."""
    engine = UQEngine(data_loader=DataLoader())
    requests = [
        DatasetRequest(file_path=runfile.PRESET_DATA2_PATHS[0], dataset_id="DATA2_A"),
        DatasetRequest(file_path=runfile.PRESET_DATA2_PATHS[1], dataset_id="DATA2_B"),
    ]
    candidate_options = [
        ModelOptions(
            mode=ExperimentMode.LAG,
            parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
            b_form="convection",
            nfe=5,
            use_sigma_logit_transform=False,
            fix_sigma_in_estimation=True,
            process_model_profile=ProcessModelProfile.DATA2,
            paper_profile="DATA2_PAPER",
        ),
        ModelOptions(
            mode=ExperimentMode.LAG,
            parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
            b_form="single",
            nfe=5,
            use_sigma_logit_transform=False,
            fix_sigma_in_estimation=True,
            process_model_profile=ProcessModelProfile.DATA2,
            paper_profile="DATA2_PAPER",
        ),
    ]

    workflow = engine.run_full_workflow(
        requests,
        candidate_options,
        parameter_scope=ParameterScope(shared_parameters=("Lp", "sigma"), dataset_specific_parameters=("beta_0",)),
        run_estimation=False,
        run_parameter_estimation_doe=False,
        run_model_discrimination_doe=False,
    )

    assert workflow.metadata["dataset_count"] == 2
    assert workflow.metadata["candidate_model_count"] == 2
    assert workflow.metadata["regression_plan"]["mode"] == "joint_shared_then_dataset_specific_refinement"


@pytest.mark.regression
def test_uqengine_parameter_estimation_doe_reports_explicit_workflow(monkeypatch: pytest.MonkeyPatch) -> None:
    """Parameter-estimation DoE should report an explicit uncertainty-reduction workflow summary."""
    loader = DataLoader()
    dataset, validation = loader.load(
        DatasetRequest(file_path=runfile.PRESET_DATA1_PATHS[0], dataset_id="DATA1_501_1")
    )
    assert validation[0] is True

    engine = UQEngine(data_loader=loader)
    options = ModelOptions(
        mode=ExperimentMode.DATA,
        parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
        b_form="single",
        nfe=10,
        use_sigma_logit_transform=False,
        process_model_profile=ProcessModelProfile.DATA1,
        paper_profile="DATA1_PAPER",
    )

    def _fake_run_doe_with_pyomo_v24(experiment, *, fd_formula="central", step=0.001, objective_option="determinant", solver_name="ipopt", tee=False):
        return {
            "fim": np.eye(2),
            "d_opt_logdet": 1.23,
            "objective_option": objective_option,
            "fd_formula": fd_formula,
            "step": step,
        }

    monkeypatch.setattr(workflow_module, "run_doe_with_pyomo_v24", _fake_run_doe_with_pyomo_v24)

    result = engine.analyze_doe_parameter_estimation(dataset, options)

    assert result["workflow_type"] == "parameter_estimation"
    assert result["design_goal"] == "reduce_parameter_uncertainty"
    assert result["score_type"] == "fim_d_optimality"
    assert result["candidate_model_count"] == 1
    assert result["recommended_design"]["objective_option"] == "determinant"
    assert result["summary"]["workflow_type"] == "parameter_estimation"


@pytest.mark.regression
def test_uqengine_model_discrimination_doe_reports_explicit_workflow(monkeypatch: pytest.MonkeyPatch) -> None:
    """Model-discrimination DoE should report explicit separation-focused outputs and ranking columns."""
    loader = DataLoader()
    dataset, validation = loader.load(
        DatasetRequest(file_path=runfile.PRESET_DATA2_PATHS[0], dataset_id="DATA2_A")
    )
    assert validation[0] is True

    engine = UQEngine(data_loader=loader)
    candidate_options = [
        ModelOptions(
            mode=ExperimentMode.LAG,
            parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
            b_form="single",
            nfe=5,
            use_sigma_logit_transform=False,
            process_model_profile=ProcessModelProfile.DATA2,
            paper_profile="DATA2_PAPER",
        ),
        ModelOptions(
            mode=ExperimentMode.LAG,
            parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
            b_form="pervial",
            nfe=5,
            use_sigma_logit_transform=False,
            process_model_profile=ProcessModelProfile.DATA2,
            paper_profile="DATA2_PAPER",
        ),
    ]

    def _fake_build_simulation_model(self, *, theta_values=None, specs_override=None, solver_name="ipopt", solver_options=None, tee=False):
        c_d_value = float((specs_override or {}).get("C_D_value", 0.0))
        scale = 1.0 if self.model_options.b_form == "single" else 2.0
        return {"signature": np.asarray([scale * c_d_value, scale], dtype=float)}

    monkeypatch.setattr(
        workflow_module.DiafiltrationExperiment,
        "build_simulation_model",
        _fake_build_simulation_model,
    )
    monkeypatch.setattr(
        workflow_module.UQEngine,
        "_design_signature",
        staticmethod(lambda model: np.asarray(model["signature"], dtype=float)),
    )

    result = engine.analyze_doe_model_discrimination(
        dataset,
        candidate_options,
        candidate_names=["single_b", "per_vial_b"],
        candidate_design_overrides=[{"C_D_value": 0.5}, {"C_D_value": 2.0}],
    )

    assert result["workflow_type"] == "model_discrimination"
    assert result["design_goal"] == "separate_candidate_models"
    assert result["score_type"] == "mean_pairwise_signature_distance"
    assert result["candidate_model_count"] == 2
    assert list(result["candidate_model_names"]) == ["single_b", "per_vial_b"]
    assert {"mean_pairwise_distance", "min_pairwise_distance", "max_pairwise_distance"}.issubset(
        result["ranking"].columns
    )
    assert result["recommended_design"]["design_id"] == "design_2"
    assert result["summary"]["workflow_type"] == "model_discrimination"


@pytest.mark.regression
@pytest.mark.nightly
def test_uqengine_run_full_workflow_reports_distinct_design_workflows(monkeypatch: pytest.MonkeyPatch) -> None:
    """Full-workflow metadata should keep parameter-estimation and model-discrimination design modes separate."""
    engine = UQEngine(data_loader=DataLoader())
    requests = [DatasetRequest(file_path=runfile.PRESET_DATA1_PATHS[0], dataset_id="DATA1_501_1")]
    candidate_options = [
        ModelOptions(
            mode=ExperimentMode.DATA,
            parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
            b_form="single",
            nfe=5,
            use_sigma_logit_transform=False,
            process_model_profile=ProcessModelProfile.DATA1,
            paper_profile="DATA1_PAPER",
        ),
        ModelOptions(
            mode=ExperimentMode.DATA,
            parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
            b_form="pervial",
            nfe=5,
            use_sigma_logit_transform=False,
            process_model_profile=ProcessModelProfile.DATA1,
            paper_profile="DATA1_PAPER",
        ),
    ]

    monkeypatch.setattr(
        workflow_module.UQEngine,
        "analyze_doe_parameter_estimation",
        lambda self, dataset, model_options, **kwargs: {
            "workflow_type": "parameter_estimation",
            "design_goal": "reduce_parameter_uncertainty",
            "score_type": "fim_d_optimality",
            "candidate_model_count": 1,
            "recommended_design": {"objective_option": "determinant", "d_opt_logdet": 1.0},
            "summary": {
                "workflow_type": "parameter_estimation",
                "design_goal": "reduce_parameter_uncertainty",
                "score_type": "fim_d_optimality",
                "candidate_model_count": 1,
                "recommended_design": {"objective_option": "determinant", "d_opt_logdet": 1.0},
            },
        },
    )
    monkeypatch.setattr(
        workflow_module.UQEngine,
        "analyze_doe_model_discrimination",
        lambda self, dataset, candidate_options, **kwargs: {
            "workflow_type": "model_discrimination",
            "design_goal": "separate_candidate_models",
            "score_type": "mean_pairwise_signature_distance",
            "candidate_model_count": len(candidate_options),
            "recommended_design": {"design_id": "design_1", "discrimination_score": 2.0},
            "summary": {
                "workflow_type": "model_discrimination",
                "design_goal": "separate_candidate_models",
                "score_type": "mean_pairwise_signature_distance",
                "candidate_model_count": len(candidate_options),
                "recommended_design": {"design_id": "design_1", "discrimination_score": 2.0},
            },
        },
    )

    workflow = engine.run_full_workflow(
        requests,
        candidate_options,
        run_estimation=False,
        run_parameter_estimation_doe=True,
        run_model_discrimination_doe=True,
    )

    assert workflow.doe_parameter_estimation["workflow_type"] == "parameter_estimation"
    assert workflow.doe_model_discrimination["workflow_type"] == "model_discrimination"
    assert workflow.metadata["design_workflows"]["parameter_estimation"]["score_type"] == "fim_d_optimality"
    assert workflow.metadata["design_workflows"]["model_discrimination"]["score_type"] == "mean_pairwise_signature_distance"
