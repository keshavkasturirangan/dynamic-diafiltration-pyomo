"""Object-oriented unified diafiltration workflow.

This module concentrates the orchestration layer for the unified DATA1/DATA2/DATA3
workflow in one place. It intentionally reuses the existing validated scientific
helpers from ``unified_codebase_library.py`` rather than rewriting the low-level
model equations during the architectural refactor.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from pyomo.contrib.parmest.experiment import Experiment as ParmestExperiment
from pyomo.contrib.parmest.parmest import Estimator

from .output_artifacts import (
    build_data1_stage_a_outputs,
    build_data1_stage_b_outputs,
    build_data2_validation_artifacts,
    legacy_style_objective_terms,
    render_data1_fig3_reference,
    render_data1_main_reference_figures,
)
from ..unified_codebase_library import (
    ExperimentalData,
    ModelOptions,
    ParameterGuess,
    ParameterTreatmentMode,
    SourceType,
    _augment_theta_with_sigma,
    _temporary_ipopt_opt,
    apply_discretization,
    build_guess_from_experiment_v24,
    evaluate_data1_paper_contour_objectives_v24,
    extract_model_trajectories_v24,
    estimate_parameters_with_parmest_v24,
    label_parmest_and_doe_suffixes,
    load_experiment_easy,
    model_construct_for_parmest_v24,
    model_construct_inter_v24,
    plot_experiment,
    run_doe_with_pyomo_v24,
    theta_component_list,
    theta_names,
)


@dataclass(frozen=True)
class DatasetRequest:
    """One dataset load request for the unified object-oriented workflow."""

    file_path: Union[str, Path]
    selector: Optional[Union[str, int]] = None
    specs: Optional[Dict[str, object]] = None
    convert_to_concentration: bool = False
    conductivity_model: str = "msa"
    conductivity_model_params: Optional[Dict[str, object]] = None
    plot: bool = False
    dataset_id: Optional[str] = None

    def resolved_path(self) -> Path:
        """Return the absolute path for the requested dataset file."""
        return Path(self.file_path).expanduser().resolve()


@dataclass(frozen=True)
class ParameterScope:
    """Declare how fitted parameters should be treated across datasets.

    ``shared_parameters`` participate in one joint fit across all datasets.
    ``dataset_specific_parameters`` are identified explicitly so the engine can
    keep their handling visible in result metadata and two-stage refinement.
    """

    shared_parameters: Tuple[str, ...] = ("Lp", "B", "beta_0", "beta_1", "sigma", "sigma_logit")
    dataset_specific_parameters: Tuple[str, ...] = ()

    def shared_name_set(self) -> set[str]:
        return {self._base_name(name) for name in self.shared_parameters}

    def dataset_specific_name_set(self) -> set[str]:
        return {self._base_name(name) for name in self.dataset_specific_parameters}

    @staticmethod
    def _base_name(name: str) -> str:
        if "[" in name:
            return name.split("[", 1)[0]
        return name


@dataclass
class UQResult:
    """Container for the main result payload returned by ``UQEngine``."""

    datasets: List[ExperimentalData]
    model_options: ModelOptions
    experiments: List["DiafiltrationExperiment"]
    load_results: List[Tuple[bool, List[Tuple[str, str]]]]
    estimation: Optional[Dict[str, object]] = None
    uncertainty: Optional[Dict[str, object]] = None
    doe_parameter_estimation: Optional[Dict[str, object]] = None
    doe_model_discrimination: Optional[Dict[str, object]] = None
    comparison: Optional[pd.DataFrame] = None
    plotting: Optional[Dict[str, pd.DataFrame]] = None
    metadata: Dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class RegressionPlan:
    """Explicit plan for one regression workflow across datasets.

    This makes the shared-vs-dataset-specific treatment visible before solving.
    """

    dataset_ids: Tuple[str, ...]
    shared_parameters: Tuple[str, ...]
    dataset_specific_parameters: Tuple[str, ...]
    mode: str
    joint_objective_description: str


@dataclass(frozen=True)
class DatasetRegressionSummary:
    """Per-dataset summary row for a joint or single regression result."""

    dataset_id: str
    n_measurements: int
    objective: float
    shared_theta: Dict[str, float]
    dataset_specific_theta: Dict[str, float]


@dataclass(frozen=True)
class ModelSelectionRow:
    """One candidate-model comparison row."""

    candidate_name: str
    objective: float
    n_measurements: int
    n_parameters: int
    aic: float
    aicc: float
    bic: float


@dataclass(frozen=True)
class DesignWorkflowSummary:
    """Compact, workflow-neutral summary for one experiment-design mode."""

    workflow_type: str
    design_goal: str
    score_type: str
    candidate_model_count: int
    recommended_design: Dict[str, object]


class DataLoader:
    """Load MAT and XLSX experiment files into a common internal representation."""

    def load(self, request: DatasetRequest) -> Tuple[ExperimentalData, Tuple[bool, List[Tuple[str, str]]]]:
        exp, validation = load_experiment_easy(
            str(request.resolved_path()),
            selector=request.selector,
            specs=request.specs,
            convert_to_concentration=bool(request.convert_to_concentration),
            conductivity_model=str(request.conductivity_model),
            conductivity_model_params=request.conductivity_model_params,
            plot=bool(request.plot),
        )
        if request.dataset_id is not None:
            exp.dataset_id = request.dataset_id
        elif exp.dataset_id is None:
            exp.dataset_id = request.resolved_path().stem
        return exp, validation

    def load_many(
        self, requests: Sequence[DatasetRequest]
    ) -> Tuple[List[ExperimentalData], List[Tuple[bool, List[Tuple[str, str]]]]]:
        datasets: List[ExperimentalData] = []
        validations: List[Tuple[bool, List[Tuple[str, str]]]] = []
        for request in requests:
            exp, validation = self.load(request)
            datasets.append(exp)
            validations.append(validation)
        return datasets, validations


class DiafiltrationExperiment(ParmestExperiment):
    """One concrete dataset + model-options combination."""

    def __init__(
        self,
        dataset: ExperimentalData,
        model_options: ModelOptions,
        *,
        guess: Optional[ParameterGuess] = None,
        dataset_label: Optional[str] = None,
        parameter_scope: Optional[ParameterScope] = None,
    ) -> None:
        super().__init__(model=None)
        self.dataset = dataset
        self.model_options = model_options
        self.guess = guess
        self.dataset_label = dataset_label or str(dataset.dataset_id or dataset.filename or "dataset")
        self.parameter_scope = parameter_scope or ParameterScope()

    @property
    def source_type(self) -> Optional[SourceType]:
        return self.dataset.source

    def build_model(self) -> pyo.ConcreteModel:
        """Build and label the model used by Parmest/DoE."""
        return model_construct_for_parmest_v24(self.dataset, options=self.model_options, guess=self.guess)

    def get_labeled_model(self):
        self.model = self.build_model()
        return self.model

    def build_simulation_model(
        self,
        *,
        theta_values: Optional[Dict[str, float]] = None,
        specs_override: Optional[Dict[str, object]] = None,
        solver_name: str = "ipopt",
        solver_options: Optional[Dict[str, object]] = None,
        tee: bool = False,
    ) -> pyo.ConcreteModel:
        """Build a simulation-only model for plotting or discrimination scoring."""
        dataset = deepcopy(self.dataset)
        if specs_override:
            for key, value in specs_override.items():
                setattr(dataset, key, value)
        sim_options = replace(
            self.model_options,
            parameter_treatment_mode=ParameterTreatmentMode.SIMULATION,
        )
        guess = build_guess_from_experiment_v24(dataset, sim_options, override=self.guess)
        if theta_values:
            guess = self._guess_with_theta_overrides(guess, theta_values)
        model = model_construct_inter_v24(dataset, options=sim_options, guess=guess)
        apply_discretization(model, nfe=sim_options.nfe, scheme=sim_options.fd_scheme)
        label_parmest_and_doe_suffixes(model, dataset, sim_options)
        solver = pyo.SolverFactory(solver_name)
        if solver_options:
            for key, value in solver_options.items():
                solver.options[key] = value
        solver.solve(model, tee=tee)
        return model

    def plot_quicklook(self, *, kind: str = "retentate_signal") -> None:
        """Expose the existing quick-look plotting helper through the experiment object."""
        plot_experiment(self.dataset, kind=kind)

    def measurement_frame(self) -> pd.DataFrame:
        """Return measured trajectories in a tidy dataframe for plotting/reporting."""
        records: List[Dict[str, object]] = []
        for vial in self.dataset.vials:
            time_s = np.asarray(vial.time_s if vial.time_s is not None else [], dtype=float)
            mass_g = np.asarray(vial.mass_g if vial.mass_g is not None else np.full(time_s.shape, np.nan), dtype=float)
            ret_signal = np.asarray(
                vial.retentate_signal if vial.retentate_signal is not None else np.full(time_s.shape, np.nan),
                dtype=float,
            )
            perm_signal = np.asarray(
                vial.permeate_signal if vial.permeate_signal is not None else np.full(time_s.shape, np.nan),
                dtype=float,
            )
            for idx, time_value in enumerate(time_s):
                records.append(
                    {
                        "dataset_id": self.dataset_label,
                        "vial": vial.number,
                        "time_s": float(time_value),
                        "mass_g": float(mass_g[idx]) if idx < len(mass_g) else np.nan,
                        "retentate_signal": float(ret_signal[idx]) if idx < len(ret_signal) else np.nan,
                        "permeate_signal": float(perm_signal[idx]) if idx < len(perm_signal) else np.nan,
                    }
                )
        return pd.DataFrame.from_records(records)

    @staticmethod
    def _guess_with_theta_overrides(base_guess: ParameterGuess, theta_values: Dict[str, float]) -> ParameterGuess:
        kwargs = {
            "Lp": float(theta_values.get("Lp", base_guess.Lp)),
            "sigma": float(theta_values.get("sigma", base_guess.sigma)),
            "B": theta_values.get("B", base_guess.B),
            "beta_0": theta_values.get("beta_0", base_guess.beta_0),
            "beta_1": theta_values.get("beta_1", base_guess.beta_1),
            "beta_2": theta_values.get("beta_2", base_guess.beta_2),
            "beta_3": theta_values.get("beta_3", base_guess.beta_3),
            "S0": theta_values.get("S0", base_guess.S0),
            "S": theta_values.get("S", base_guess.S),
        }
        return ParameterGuess(**kwargs)


class _ScopedParmestExperiment(DiafiltrationExperiment):
    """Parmest wrapper that keeps dataset-specific parameters explicit in labels.

    Parmest natively handles shared-theta multi-experiment estimation well. We use
    this scoped variant to keep dataset-specific parameter treatment visible in the
    model metadata and in result reporting. The exact mixed-scope estimation path is
    implemented as a two-stage refinement in ``UQEngine``.
    """

    def get_labeled_model(self):
        model = super().get_labeled_model()
        model._uq_dataset_label = self.dataset_label
        return model


class UQEngine:
    """Unified uncertainty-quantification, comparison, and DoE engine."""

    def __init__(self, *, data_loader: Optional[DataLoader] = None) -> None:
        self.data_loader = data_loader or DataLoader()

    def load_datasets(
        self, requests: Sequence[DatasetRequest]
    ) -> Tuple[List[ExperimentalData], List[Tuple[bool, List[Tuple[str, str]]]]]:
        return self.data_loader.load_many(requests)

    def build_experiments(
        self,
        datasets: Sequence[ExperimentalData],
        model_options: ModelOptions,
        *,
        parameter_scope: Optional[ParameterScope] = None,
        guesses: Optional[Sequence[Optional[ParameterGuess]]] = None,
    ) -> List[DiafiltrationExperiment]:
        scope = parameter_scope or ParameterScope()
        out: List[DiafiltrationExperiment] = []
        guess_list = list(guesses) if guesses is not None else [None] * len(datasets)
        for dataset, guess in zip(datasets, guess_list):
            out.append(
                _ScopedParmestExperiment(
                    dataset=dataset,
                    model_options=model_options,
                    guess=guess,
                    dataset_label=str(dataset.dataset_id or dataset.filename or "dataset"),
                    parameter_scope=scope,
                )
            )
        return out

    def build_regression_plan(
        self,
        datasets: Sequence[ExperimentalData],
        parameter_scope: Optional[ParameterScope] = None,
    ) -> RegressionPlan:
        """Return an explicit regression plan before fitting."""
        scope = parameter_scope or ParameterScope()
        dataset_ids = tuple(str(dataset.dataset_id or dataset.filename or "dataset") for dataset in datasets)
        if len(datasets) <= 1 and not scope.dataset_specific_parameters:
            mode = "single_dataset"
        elif not scope.dataset_specific_parameters:
            mode = "joint_shared"
        else:
            mode = "joint_shared_then_dataset_specific_refinement"
        joint_objective_description = (
            "Weighted SSE summed across datasets with shared-theta joint solve"
            if len(datasets) > 1
            else "Weighted SSE for one dataset"
        )
        return RegressionPlan(
            dataset_ids=dataset_ids,
            shared_parameters=tuple(sorted(scope.shared_name_set())),
            dataset_specific_parameters=tuple(sorted(scope.dataset_specific_name_set())),
            mode=mode,
            joint_objective_description=joint_objective_description,
        )

    def fit_parameters(
        self,
        datasets: Sequence[ExperimentalData],
        model_options: ModelOptions,
        *,
        parameter_scope: Optional[ParameterScope] = None,
        guess: Optional[ParameterGuess] = None,
        calc_cov: bool = True,
        cov_n: Optional[int] = None,
        solver: str = "ipopt",
        solver_options: Optional[Dict[str, object]] = None,
        tee: bool = False,
    ) -> Dict[str, object]:
        """Run the main parameter-estimation workflow.

        Shared-parameter multi-dataset regression uses the existing Parmest path
        directly. When dataset-specific parameters are declared, the engine makes
        the distinction explicit and performs a two-stage refinement:

        1. joint fit of shared parameters across datasets
        2. per-dataset refinement with shared parameters fixed
        """

        scope = parameter_scope or ParameterScope()
        if not datasets:
            raise ValueError("fit_parameters requires at least one dataset.")
        plan = self.build_regression_plan(datasets, scope)

        if len(datasets) == 1 and not scope.dataset_specific_parameters:
            result = estimate_parameters_with_parmest_v24(
                experiments=list(datasets),
                options=model_options,
                guess=guess,
                calc_cov=calc_cov,
                cov_n=cov_n,
                solver=solver,
                solver_options=solver_options,
                tee=tee,
            )
            result["parameter_scope"] = {
                "shared": sorted(scope.shared_name_set()),
                "dataset_specific": sorted(scope.dataset_specific_name_set()),
                "mode": "single_dataset",
            }
            result["regression_plan"] = self._plan_to_dict(plan)
            result["dataset_summaries"] = self._build_dataset_summaries(
                datasets,
                model_options,
                theta_values=self._coerce_theta_dict(result.get("theta")),
            )
            return result

        joint_solver = "ef_ipopt" if solver == "ipopt" else solver

        if not scope.dataset_specific_parameters:
            result = estimate_parameters_with_parmest_v24(
                experiments=list(datasets),
                options=model_options,
                guess=guess,
                calc_cov=calc_cov,
                cov_n=cov_n,
                solver=joint_solver,
                solver_options=solver_options,
                tee=tee,
            )
            result["parameter_scope"] = {
                "shared": sorted(scope.shared_name_set()),
                "dataset_specific": sorted(scope.dataset_specific_name_set()),
                "mode": "joint_shared",
            }
            result["per_dataset_objective"] = self._per_dataset_objective_breakdown(
                datasets,
                model_options,
                theta_values=self._coerce_theta_dict(result.get("theta")),
            )
            result["regression_plan"] = self._plan_to_dict(plan)
            result["dataset_summaries"] = self._build_dataset_summaries(
                datasets,
                model_options,
                theta_values=self._coerce_theta_dict(result.get("theta")),
            )
            return result

        joint_options = replace(model_options)
        joint_result = estimate_parameters_with_parmest_v24(
            experiments=list(datasets),
            options=joint_options,
            guess=guess,
            calc_cov=calc_cov,
            cov_n=cov_n,
            solver=joint_solver,
            solver_options=solver_options,
            tee=tee,
        )
        shared_theta = self._coerce_theta_dict(joint_result.get("theta"))
        stage_two = self._fit_dataset_specific_refinements(
            datasets=datasets,
            model_options=model_options,
            shared_theta=shared_theta,
            dataset_specific_names=scope.dataset_specific_name_set(),
            solver=solver,
            solver_options=solver_options,
            tee=tee,
        )
        joint_result["parameter_scope"] = {
            "shared": sorted(scope.shared_name_set()),
            "dataset_specific": sorted(scope.dataset_specific_name_set()),
            "mode": "joint_shared_then_dataset_specific_refinement",
        }
        joint_result["dataset_specific_refinement"] = stage_two
        joint_result["regression_plan"] = self._plan_to_dict(plan)
        joint_result["dataset_summaries"] = self._build_dataset_summaries(
            datasets,
            model_options,
            theta_values=shared_theta,
            dataset_specific_refinement=stage_two,
        )
        return joint_result

    def summarize_uncertainty(self, estimation_result: Dict[str, object]) -> Dict[str, object]:
        """Summarize covariance-driven uncertainty when available."""
        summary: Dict[str, object] = {
            "covariance_method": estimation_result.get("covariance_method"),
            "covariance_warning": estimation_result.get("covariance_warning"),
        }
        covariance = estimation_result.get("covariance")
        theta = self._coerce_theta_dict(estimation_result.get("theta"))
        if covariance is None:
            summary["parameter_standard_errors"] = {}
            return summary
        cov = np.asarray(covariance, dtype=float)
        diag = np.sqrt(np.clip(np.diag(cov), a_min=0.0, a_max=None))
        param_names = list(theta.keys())[: len(diag)]
        summary["parameter_standard_errors"] = {
            name: float(diag[idx]) for idx, name in enumerate(param_names)
        }
        return summary

    def analyze_doe_parameter_estimation(
        self,
        dataset: ExperimentalData,
        model_options: ModelOptions,
        *,
        guess: Optional[ParameterGuess] = None,
        solver_name: str = "ipopt",
        tee: bool = False,
    ) -> Dict[str, object]:
        """Run the existing Pyomo.DoE path for parameter-estimation-focused design."""
        exp_obj = DiafiltrationExperiment(dataset=dataset, model_options=model_options, guess=guess)
        result = run_doe_with_pyomo_v24(exp_obj, solver_name=solver_name, tee=tee)
        result["workflow_type"] = "parameter_estimation"
        result["design_goal"] = "reduce_parameter_uncertainty"
        result["score_type"] = "fim_d_optimality"
        result["candidate_model_count"] = 1
        result["candidate_model_names"] = [self._model_option_name(model_options)]
        result["recommended_design"] = {
            "objective_option": "determinant",
            "d_opt_logdet": float(result.get("d_opt_logdet", np.nan)),
        }
        result["summary"] = self._design_summary_to_dict(
            DesignWorkflowSummary(
                workflow_type="parameter_estimation",
                design_goal="reduce_parameter_uncertainty",
                score_type="fim_d_optimality",
                candidate_model_count=1,
                recommended_design=dict(result["recommended_design"]),
            )
        )
        return result

    def analyze_doe_model_discrimination(
        self,
        dataset: ExperimentalData,
        candidate_options: Sequence[ModelOptions],
        *,
        fitted_thetas: Optional[Sequence[Optional[Dict[str, float]]]] = None,
        candidate_design_overrides: Optional[Sequence[Dict[str, object]]] = None,
        candidate_names: Optional[Sequence[str]] = None,
    ) -> Dict[str, object]:
        """Score candidate future experiments for model discrimination.

        This is intentionally explicit about being a discrimination workflow. It
        reuses the current process models and ranks candidate design settings by
        predicted separation across model structures rather than by FIM.
        """

        if not candidate_options:
            raise ValueError("analyze_doe_model_discrimination requires at least one candidate model.")
        design_overrides = list(candidate_design_overrides or [{"C_D_value": dataset.C_D_value}])
        theta_list = list(fitted_thetas or [None] * len(candidate_options))
        model_names = list(candidate_names or [self._model_option_name(options) for options in candidate_options])
        rows: List[Dict[str, object]] = []
        for design_idx, override in enumerate(design_overrides, start=1):
            signatures: List[np.ndarray] = []
            pairwise_rows: List[Dict[str, object]] = []
            for options, theta in zip(candidate_options, theta_list):
                experiment = DiafiltrationExperiment(dataset=dataset, model_options=options)
                model = experiment.build_simulation_model(
                    theta_values=theta,
                    specs_override=override,
                )
                signatures.append(self._design_signature(model))
            distances = self._pairwise_signature_distances(signatures, model_names=model_names)
            score = float(np.mean([row["distance"] for row in distances])) if distances else 0.0
            min_distance = float(min((row["distance"] for row in distances), default=0.0))
            max_distance = float(max((row["distance"] for row in distances), default=0.0))
            rows.append(
                {
                    "design_id": f"design_{design_idx}",
                    "design_override": dict(override),
                    "discrimination_score": float(score),
                    "mean_pairwise_distance": float(score),
                    "min_pairwise_distance": min_distance,
                    "max_pairwise_distance": max_distance,
                    "pairwise_distances": distances,
                }
            )
        ranking = pd.DataFrame(rows).sort_values("discrimination_score", ascending=False).reset_index(drop=True)
        recommended_design = ranking.iloc[0].to_dict() if not ranking.empty else None
        return {
            "workflow_type": "model_discrimination",
            "design_goal": "separate_candidate_models",
            "score_type": "mean_pairwise_signature_distance",
            "candidate_model_count": len(candidate_options),
            "candidate_model_names": model_names,
            "ranking": ranking,
            "recommended_design": recommended_design,
            "summary": self._design_summary_to_dict(
                DesignWorkflowSummary(
                    workflow_type="model_discrimination",
                    design_goal="separate_candidate_models",
                    score_type="mean_pairwise_signature_distance",
                    candidate_model_count=len(candidate_options),
                    recommended_design=dict(recommended_design or {}),
                )
            ),
        }

    def compare_models(
        self,
        datasets: Sequence[ExperimentalData],
        candidate_options: Sequence[ModelOptions],
        *,
        candidate_names: Optional[Sequence[str]] = None,
        parameter_scope: Optional[ParameterScope] = None,
        solver: str = "ipopt",
        solver_options: Optional[Dict[str, object]] = None,
        tee: bool = False,
    ) -> pd.DataFrame:
        """Fit and compare candidate model formulations with AIC-style criteria."""
        names = list(candidate_names) if candidate_names is not None else [
            f"candidate_{idx + 1}" for idx in range(len(candidate_options))
        ]
        rows: List[ModelSelectionRow] = []
        for candidate_name, options in zip(names, candidate_options):
            fit = self.fit_parameters(
                datasets,
                options,
                parameter_scope=parameter_scope,
                calc_cov=False,
                solver=solver,
                solver_options=solver_options,
                tee=tee,
            )
            objective = float(fit.get("objective", np.nan))
            theta = self._coerce_theta_dict(fit.get("theta"))
            n_parameters = max(len(theta), 1)
            n_measurements = self._count_measurements(datasets, options)
            aic, aicc, bic = self._information_criteria(
                objective=objective,
                n_measurements=n_measurements,
                n_parameters=n_parameters,
            )
            rows.append(
                ModelSelectionRow(
                    candidate_name=candidate_name,
                    objective=objective,
                    n_measurements=n_measurements,
                    n_parameters=n_parameters,
                    aic=aic,
                    aicc=aicc,
                    bic=bic,
                )
            )
        comparison = pd.DataFrame([row.__dict__ for row in rows]).sort_values(
            ["aic", "objective"], ascending=[True, True]
        ).reset_index(drop=True)
        comparison["delta_aic"] = comparison["aic"] - float(comparison["aic"].min())
        return comparison

    def generate_plot_bundle(self, experiments: Sequence[DiafiltrationExperiment]) -> Dict[str, pd.DataFrame]:
        """Collect per-experiment tidy measurement frames for downstream plotting."""
        return {experiment.dataset_label: experiment.measurement_frame() for experiment in experiments}

    def simulate_trajectories(
        self,
        dataset: ExperimentalData,
        model_options: ModelOptions,
        *,
        theta_values: Optional[Dict[str, float]] = None,
        guess: Optional[ParameterGuess] = None,
        solver_name: str = "ipopt",
        solver_options: Optional[Dict[str, object]] = None,
        tee: bool = False,
    ) -> Dict[int, Dict[str, np.ndarray]]:
        """Build a simulation model and return extracted per-vial trajectories."""
        experiment = DiafiltrationExperiment(dataset=dataset, model_options=model_options, guess=guess)
        model = experiment.build_simulation_model(
            theta_values=theta_values,
            solver_name=solver_name,
            solver_options=solver_options,
            tee=tee,
        )
        return extract_model_trajectories_v24(dataset, model)

    def build_data1_validation_artifacts(
        self,
        *,
        out_dir: Path,
        dataset: ExperimentalData,
        estimation_result: Dict[str, object],
        simulation_options: ModelOptions,
        references: pd.DataFrame,
        tolerances,
        solver_name: str = "ipopt",
        solver_options: Optional[Dict[str, object]] = None,
        tee: bool = False,
    ) -> Dict[str, object]:
        """Build DATA1 Stage A/B validation artifacts through the workflow layer."""
        theta = self._coerce_theta_dict(estimation_result.get("theta"))
        trajectories = self.simulate_trajectories(
            dataset,
            simulation_options,
            theta_values=theta,
            solver_name=solver_name,
            solver_options=solver_options,
            tee=tee,
        )
        artifacts = {
            "stage_a": build_data1_stage_a_outputs(
                out_dir=out_dir,
                exp=dataset,
                est=estimation_result,
                sim=trajectories,
                tolerances=tolerances,
                references=references,
            ),
            "stage_b": build_data1_stage_b_outputs(
                out_dir=out_dir,
                exp=dataset,
                est=estimation_result,
                sim=trajectories,
            ),
        }
        sim_validation_root = (
            Path(__file__).resolve().parents[3]
            / "docs"
            / "validation"
            / "simulation_validation_data_files"
        )
        fig3_refs = render_data1_fig3_reference(
            figures_dir=out_dir / "figures",
            simulation_validation_root=sim_validation_root,
        )
        if fig3_refs:
            artifacts["stage_a"].setdefault("reference_style_figures", []).extend(fig3_refs)
        main_refs = render_data1_main_reference_figures(
            figures_dir=out_dir / "figures",
            simulation_validation_root=sim_validation_root,
        )
        if main_refs:
            artifacts["stage_b"].setdefault("reference_style_figures", []).extend(main_refs)

        # Direct unified sweeps for the fitted DATA1 dataset. These do not yet
        # replace every multi-dataset paper panel, but they ensure the current
        # run owns its own sigma-sensitivity and contour source data.
        tables_dir = out_dir / "tables"
        figures_dir = out_dir / "figures"
        tables_dir.mkdir(parents=True, exist_ok=True)
        figures_dir.mkdir(parents=True, exist_ok=True)

        direct_figures: list[str] = []
        direct_tables: list[str] = []
        t_delay = float(dataset.vials[0].time_s[0])

        # Figure 4-style sigma sensitivity trajectories from the unified model.
        sigma_values = [0.1, 0.5, 0.9]
        sigma_rows: list[dict[str, object]] = []
        for sigma in sigma_values:
            theta_sigma = dict(theta)
            theta_sigma["sigma"] = float(sigma)
            sigma_sim = self.simulate_trajectories(
                dataset,
                simulation_options,
                theta_values=theta_sigma,
                solver_name=solver_name,
                solver_options=solver_options,
                tee=tee,
            )
            for n, vial in enumerate(dataset.vials, start=1):
                for time_s, m_val, cH_val, cF_val in zip(
                    sigma_sim[n]["time_s"],
                    sigma_sim[n]["mV"],
                    sigma_sim[n]["cH"],
                    sigma_sim[n]["cF"],
                ):
                    sigma_rows.append(
                        {
                            "sigma": float(sigma),
                            "vial_number": n,
                            "time_min": float(time_s) / 60.0,
                            "mass_g": float(m_val),
                            "permeate_concentration_mM": float(cH_val),
                            "retentate_concentration_mM": float(cF_val),
                        }
                    )
        sigma_df = pd.DataFrame(sigma_rows)
        sigma_csv = tables_dir / "data1_direct_fig4_sigma_sensitivity.csv"
        sigma_df.to_csv(sigma_csv, index=False)
        direct_tables.append(str(sigma_csv))

        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(12, 4))
        sigma_styles = {
            0.1: {"color": "#FF4D4D", "linestyle": "--"},
            0.5: {"color": "#1F4AFF", "linestyle": "-"},
            0.9: {"color": "#5CB85C", "linestyle": ":"},
        }
        for idx, (value_col, title, ylabel) in enumerate(
            (
                ("mass_g", "Mass", "Mass [g]"),
                ("permeate_concentration_mM", "Permeate", "Permeate Conc. [mM]"),
                ("retentate_concentration_mM", "Retentate", "Retentate Conc. [mM]"),
            ),
            start=1,
        ):
            ax = fig.add_subplot(1, 3, idx)
            for sigma in sigma_values:
                style = sigma_styles[sigma]
                subset = sigma_df[sigma_df["sigma"] == sigma]
                for _, group in subset.groupby("vial_number"):
                    group = group.sort_values("time_min")
                    ax.plot(group["time_min"], group[value_col], color=style["color"], linestyle=style["linestyle"], linewidth=2.0)
            ax.set_title(title, fontsize=16, fontweight="bold")
            ax.set_xlabel("Time [min]", fontsize=13, fontweight="bold")
            ax.set_ylabel(ylabel, fontsize=13, fontweight="bold")
            ax.tick_params(direction="in", top=True, right=True, labelsize=10)
        handles = [plt.Line2D([0], [0], color=sigma_styles[s]["color"], linestyle=sigma_styles[s]["linestyle"], linewidth=2.0) for s in sigma_values]
        fig.axes[-1].legend(handles, [r"$\sigma$ = 0.1", r"$\sigma$ = 0.5", r"$\sigma$ = 0.9"], fontsize=11, loc="best")
        fig.tight_layout()
        fig4_direct = figures_dir / "data1_fig4_direct_sigma_sensitivity.png"
        fig.savefig(fig4_direct, dpi=300, bbox_inches="tight")
        plt.close(fig)
        direct_figures.append(str(fig4_direct))

        # Figure 5/6-style direct contour sweeps from the unified model.
        def _normalize_dataset_label(value: object) -> str:
            text = str(value)
            if "_" in text:
                text = text.split("_")[-1]
            return text.strip()

        def _baseline_grid(fig_key: str, x_name: str) -> tuple[np.ndarray, np.ndarray, str | None]:
            baseline_dir = simulation_validation_root / "data1_main" / fig_key
            dataset_label = _normalize_dataset_label(getattr(dataset, "dataset_id", dataset_id))
            for candidate in sorted(baseline_dir.glob(f"data1_main_{fig_key}_*_mass.csv")):
                baseline = pd.read_csv(candidate)
                baseline_label = _normalize_dataset_label(baseline["dataset"].iloc[0])
                if baseline_label != dataset_label:
                    continue
                x_values = np.sort(baseline[x_name].astype(float).unique())
                lp_values = np.sort(baseline["Lp"].astype(float).unique())
                panel_id = str(baseline["panel_id"].iloc[0]).upper()
                return x_values, lp_values, panel_id

            if x_name == "sigma":
                return np.linspace(0.0, 1.0, 9), np.linspace(0.5, 7.4, 9), None
            return np.linspace(1e-6, 2.0, 9), np.linspace(0.5, 7.4, 9), None

        def _direct_grid(x_name: str, x_values: np.ndarray, lp_values: np.ndarray) -> pd.DataFrame:
            rows: list[dict[str, object]] = []
            for x_val in x_values:
                for lp_val in lp_values:
                    theta_grid = dict(theta)
                    theta_grid["Lp"] = float(lp_val)
                    theta_grid[x_name] = float(x_val)
                    sim_grid = self.simulate_trajectories(
                        dataset,
                        simulation_options,
                        theta_values=theta_grid,
                        solver_name=solver_name,
                        solver_options=solver_options,
                        tee=False,
                    )
                    scores = evaluate_data1_paper_contour_objectives_v24(dataset, sim_grid).as_dict()
                    rows.append(
                        {
                            x_name: float(x_val),
                            "Lp": float(lp_val),
                            "Obj_mass": scores["mass_log10"],
                            "Obj_concentration": scores["permeate_log10"],
                            "Obj_retentate_concentration": scores["retentate_log10"],
                        }
                    )
            return pd.DataFrame(rows)

        def _write_direct_grid_report(
            *,
            direct_df: pd.DataFrame,
            fig_key: str,
            panel_id: str | None,
            x_name: str,
            report_name: str,
        ) -> str | None:
            if not panel_id:
                return None
            objective_map = {
                "mass": "Obj_mass",
                "permeate": "Obj_concentration",
                "retentate": "Obj_retentate_concentration",
            }
            rows: list[dict[str, object]] = []
            for objective_name, direct_col in objective_map.items():
                baseline_path = (
                    simulation_validation_root
                    / "data1_main"
                    / fig_key
                    / f"data1_main_{fig_key}_{panel_id.lower()}_{objective_name}.csv"
                )
                if not baseline_path.exists():
                    continue
                baseline = pd.read_csv(baseline_path)
                merged = baseline.merge(
                    direct_df[[x_name, "Lp", direct_col]],
                    on=[x_name, "Lp"],
                    how="inner",
                )
                if merged.empty:
                    rows.append(
                        {
                            "figure_id": fig_key,
                            "panel_id": panel_id,
                            "objective_name": objective_name,
                            "n_points": 0,
                            "mae": np.nan,
                            "max_abs_err": np.nan,
                            "rmse": np.nan,
                        }
                    )
                    continue
                err = merged[direct_col].to_numpy(dtype=float) - merged["objective_value"].to_numpy(dtype=float)
                rows.append(
                    {
                        "figure_id": fig_key,
                        "panel_id": panel_id,
                        "objective_name": objective_name,
                        "n_points": int(len(merged)),
                        "mae": float(np.mean(np.abs(err))),
                        "max_abs_err": float(np.max(np.abs(err))),
                        "rmse": float(np.sqrt(np.mean(err**2))),
                    }
                )
            if not rows:
                return None
            report_path = tables_dir / report_name
            pd.DataFrame(rows).to_csv(report_path, index=False)
            return str(report_path)

        sigma_values_grid, lp_values, sigma_panel_id = _baseline_grid("fig5", "sigma")
        sigma_grid_df = _direct_grid("sigma", sigma_values_grid, lp_values)
        sigma_grid_csv = tables_dir / "data1_direct_fig5_sigma_lp_grid.csv"
        sigma_grid_df.to_csv(sigma_grid_csv, index=False)
        direct_tables.append(str(sigma_grid_csv))
        sigma_report = _write_direct_grid_report(
            direct_df=sigma_grid_df,
            fig_key="fig5",
            panel_id=sigma_panel_id,
            x_name="sigma",
            report_name="data1_direct_fig5_sigma_lp_comparison.csv",
        )
        if sigma_report:
            direct_tables.append(sigma_report)

        b_values_grid, b_lp_values, b_panel_id = _baseline_grid("fig6", "B")
        b_grid_df = _direct_grid("B", b_values_grid, b_lp_values)
        b_grid_csv = tables_dir / "data1_direct_fig6_b_lp_grid.csv"
        b_grid_df.to_csv(b_grid_csv, index=False)
        direct_tables.append(str(b_grid_csv))
        b_report = _write_direct_grid_report(
            direct_df=b_grid_df,
            fig_key="fig6",
            panel_id=b_panel_id,
            x_name="B",
            report_name="data1_direct_fig6_b_lp_comparison.csv",
        )
        if b_report:
            direct_tables.append(b_report)

        def _render_direct_contour(df: pd.DataFrame, x_name: str, out_name: str) -> str:
            fig = plt.figure(figsize=(12, 4))
            for idx, (col, title) in enumerate(
                (
                    ("Obj_mass", "Mass"),
                    ("Obj_concentration", "Permeate"),
                    ("Obj_retentate_concentration", "Retentate"),
                ),
                start=1,
            ):
                ax = fig.add_subplot(1, 3, idx)
                grid = df[[x_name, "Lp", col]].copy().rename(columns={col: "objective_value"})
                from .output_artifacts import build_paper_contour_grid_v24  # local import to avoid circularity

                contour = build_paper_contour_grid_v24(grid, objective_column="objective_value")
                levels = np.linspace(float(np.nanmin(contour.z_grid)), float(np.nanmax(contour.z_grid)), 10)
                cs = ax.contour(contour.x_grid, contour.y_grid, contour.z_grid, levels=levels, linewidths=1.6)
                ax.clabel(cs, cs.levels[::2], inline=True, fontsize=9, colors="k", fmt="%1.1f")
                ax.plot(contour.optimum_x, contour.optimum_y, "^", markersize=8, markeredgecolor="red", markerfacecolor=(1.0, 0.7, 0.7), clip_on=False)
                ax.set_title(title, fontsize=16, fontweight="bold")
                ax.set_xlabel(contour.x_label, fontsize=12, fontweight="bold")
                ax.set_ylabel(contour.y_label, fontsize=12, fontweight="bold")
                ax.tick_params(direction="in", top=True, right=True, labelsize=10)
                ax.set_ylim(0.5, 7.4)
                if x_name == "sigma":
                    ax.set_xlim(0.0, 1.0)
                else:
                    ax.set_xlim(0.0, 2.0)
            fig.tight_layout()
            out_path = figures_dir / out_name
            fig.savefig(out_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            return str(out_path)

        direct_figures.append(_render_direct_contour(sigma_grid_df, "sigma", "data1_fig5_direct_sigma_lp_contours.png"))
        direct_figures.append(_render_direct_contour(b_grid_df, "B", "data1_fig6_direct_b_lp_contours.png"))

        artifacts["stage_b"]["direct_unified_tables"] = direct_tables
        artifacts["stage_b"]["direct_unified_figures"] = direct_figures
        return artifacts

    def build_data2_validation_artifacts(
        self,
        *,
        out_dir: Path,
        dataset_id: str,
        dataset: ExperimentalData,
        curve_metrics: Dict[str, float],
        simulation_validation_root: Path,
    ) -> Dict[str, object]:
        """Build DATA2 validation artifacts through the workflow layer."""
        return build_data2_validation_artifacts(
            out_dir=out_dir,
            dataset_id=dataset_id,
            exp=dataset,
            curve_metrics=curve_metrics,
            simulation_validation_root=simulation_validation_root,
        )

    def recommend_next_experiment(
        self,
        *,
        dataset: ExperimentalData,
        candidate_options: Sequence[ModelOptions],
        fitted_thetas: Optional[Sequence[Optional[Dict[str, float]]]] = None,
        candidate_design_overrides: Optional[Sequence[Dict[str, object]]] = None,
        workflow: str = "model_discrimination",
    ) -> Dict[str, object]:
        """Return the recommended next experiment for the requested workflow."""
        if workflow == "parameter_estimation":
            if not candidate_options:
                raise ValueError("At least one model option is required for parameter-estimation DoE.")
            doe = self.analyze_doe_parameter_estimation(dataset, candidate_options[0])
            return {
                "workflow_type": workflow,
                "recommended_design": {
                    "objective_option": "determinant",
                    "d_opt_logdet": float(doe.get("d_opt_logdet", np.nan)),
                },
                "details": doe,
            }
        if workflow == "model_discrimination":
            discrimination = self.analyze_doe_model_discrimination(
                dataset,
                candidate_options,
                fitted_thetas=fitted_thetas,
                candidate_design_overrides=candidate_design_overrides,
            )
            return {
                "workflow_type": workflow,
                "recommended_design": discrimination.get("recommended_design"),
                "details": discrimination,
            }
        raise ValueError(f"Unknown workflow '{workflow}'. Expected 'parameter_estimation' or 'model_discrimination'.")

    def run_full_workflow(
        self,
        requests: Sequence[DatasetRequest],
        candidate_options: Sequence[ModelOptions],
        *,
        parameter_scope: Optional[ParameterScope] = None,
        run_estimation: bool = True,
        run_parameter_estimation_doe: bool = False,
        run_model_discrimination_doe: bool = False,
        calc_cov: bool = True,
        solver: str = "ipopt",
        solver_options: Optional[Dict[str, object]] = None,
        tee: bool = False,
    ) -> UQResult:
        """Run the main DataLoader -> UQEngine workflow across datasets/options."""
        datasets, load_results = self.load_datasets(requests)
        if not candidate_options:
            raise ValueError("run_full_workflow requires at least one candidate model option.")

        primary_options = candidate_options[0]
        experiments = self.build_experiments(datasets, primary_options, parameter_scope=parameter_scope)
        plotting = self.generate_plot_bundle(experiments)

        estimation = None
        uncertainty = None
        comparison = None
        doe_estimation = None
        doe_discrimination = None

        if run_estimation:
            estimation = self.fit_parameters(
                datasets,
                primary_options,
                parameter_scope=parameter_scope,
                calc_cov=calc_cov,
                solver=solver,
                solver_options=solver_options,
                tee=tee,
            )
            uncertainty = self.summarize_uncertainty(estimation)
            comparison = self.compare_models(
                datasets,
                candidate_options,
                parameter_scope=parameter_scope,
                solver=solver,
                solver_options=solver_options,
                tee=tee,
            )
        if run_parameter_estimation_doe and datasets:
            doe_estimation = self.analyze_doe_parameter_estimation(datasets[0], primary_options, solver_name=solver, tee=tee)
        if run_model_discrimination_doe and datasets:
            fitted_thetas = None
            if comparison is not None and estimation is not None:
                fitted_thetas = [self._coerce_theta_dict(estimation.get("theta"))] + [None] * (len(candidate_options) - 1)
            doe_discrimination = self.analyze_doe_model_discrimination(
                datasets[0],
                candidate_options,
                fitted_thetas=fitted_thetas,
            )

        return UQResult(
            datasets=datasets,
            model_options=primary_options,
            experiments=experiments,
            load_results=load_results,
            estimation=estimation,
            uncertainty=uncertainty,
            doe_parameter_estimation=doe_estimation,
            doe_model_discrimination=doe_discrimination,
            comparison=comparison,
            plotting=plotting,
            metadata={
                "dataset_count": len(datasets),
                "candidate_model_count": len(candidate_options),
                "regression_plan": self._plan_to_dict(self.build_regression_plan(datasets, parameter_scope)),
                "parameter_scope": {
                    "shared": sorted((parameter_scope or ParameterScope()).shared_name_set()),
                    "dataset_specific": sorted((parameter_scope or ParameterScope()).dataset_specific_name_set()),
                },
                "design_workflows": {
                    "parameter_estimation": self._extract_design_summary(doe_estimation),
                    "model_discrimination": self._extract_design_summary(doe_discrimination),
                },
            },
        )

    def _fit_dataset_specific_refinements(
        self,
        *,
        datasets: Sequence[ExperimentalData],
        model_options: ModelOptions,
        shared_theta: Dict[str, float],
        dataset_specific_names: set[str],
        solver: str,
        solver_options: Optional[Dict[str, object]],
        tee: bool,
    ) -> Dict[str, object]:
        rows: List[Dict[str, object]] = []
        for dataset in datasets:
            fit = self._refine_dataset_specific_with_fixed_shared(
                dataset=dataset,
                model_options=model_options,
                shared_theta=shared_theta,
                dataset_specific_names=dataset_specific_names,
                solver=solver,
                solver_options=solver_options,
                tee=tee,
            )
            theta = self._coerce_theta_dict(fit.get("theta"))
            shared_subset = {name: theta.get(name, shared_theta.get(name)) for name in shared_theta}
            specific_subset = {
                name: value for name, value in theta.items() if ParameterScope._base_name(name) in dataset_specific_names
            }
            rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "objective": float(fit.get("objective", np.nan)),
                    "shared_theta": shared_subset,
                    "dataset_specific_theta": specific_subset,
                }
            )
        return {"rows": rows}

    def _build_dataset_summaries(
        self,
        datasets: Sequence[ExperimentalData],
        model_options: ModelOptions,
        *,
        theta_values: Dict[str, float],
        dataset_specific_refinement: Optional[Dict[str, object]] = None,
    ) -> List[Dict[str, object]]:
        refinement_rows = {
            str(row["dataset_id"]): row for row in (dataset_specific_refinement or {}).get("rows", [])
        }
        shared_names = {ParameterScope._base_name(name) for name in theta_values}
        summaries: List[DatasetRegressionSummary] = []
        objectives = self._per_dataset_objective_breakdown(
            datasets,
            model_options,
            theta_values=theta_values,
        )
        objective_by_dataset = {str(row["dataset_id"]): float(row["objective"]) for row in objectives}
        for dataset in datasets:
            dataset_id = str(dataset.dataset_id or dataset.filename or "dataset")
            refinement = refinement_rows.get(dataset_id, {})
            shared_theta = refinement.get("shared_theta", {name: theta_values.get(name) for name in shared_names})
            dataset_specific_theta = refinement.get("dataset_specific_theta", {})
            summaries.append(
                DatasetRegressionSummary(
                    dataset_id=dataset_id,
                    n_measurements=self._count_measurements([dataset], model_options),
                    objective=objective_by_dataset.get(dataset_id, np.nan),
                    shared_theta=dict(shared_theta),
                    dataset_specific_theta=dict(dataset_specific_theta),
                )
            )
        return [summary.__dict__ for summary in summaries]

    def _refine_dataset_specific_with_fixed_shared(
        self,
        *,
        dataset: ExperimentalData,
        model_options: ModelOptions,
        shared_theta: Dict[str, float],
        dataset_specific_names: set[str],
        solver: str,
        solver_options: Optional[Dict[str, object]],
        tee: bool,
    ) -> Dict[str, object]:
        guess = self._guess_from_theta_dict(dataset, model_options, shared_theta)
        model = model_construct_for_parmest_v24(dataset, options=model_options, guess=guess)
        for component in theta_component_list(model, model_options):
            base_name = ParameterScope._base_name(component.name)
            if base_name in dataset_specific_names:
                continue
            if base_name in shared_theta and hasattr(component, "fix"):
                component.fix(float(shared_theta[base_name]))

        def _weighted_sse_expr():
            terms = []
            for output_var in model.experiment_outputs.keys():
                y_obs = float(model.experiment_outputs[output_var])
                sigma = float(model.measurement_error[output_var])
                sigma_eff = sigma if abs(sigma) > 1e-12 else 1.0
                terms.append(((output_var - y_obs) / sigma_eff) ** 2)
            return sum(terms) if terms else 0.0

        model._dataset_specific_wsse = pyo.Objective(expr=_weighted_sse_expr(), sense=pyo.minimize)
        solver_obj = pyo.SolverFactory(solver)
        if solver_options:
            for key, value in solver_options.items():
                solver_obj.options[key] = value
        raw = solver_obj.solve(model, tee=tee)
        theta = {component.name: float(pyo.value(component)) for component in theta_component_list(model, model_options)}
        theta = _augment_theta_with_sigma(theta, options=model_options)
        return {
            "raw": raw,
            "objective": float(pyo.value(model._dataset_specific_wsse)),
            "theta": pd.Series(theta),
            "model": model,
        }

    def _per_dataset_objective_breakdown(
        self,
        datasets: Sequence[ExperimentalData],
        model_options: ModelOptions,
        *,
        theta_values: Dict[str, float],
    ) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        for dataset in datasets:
            fit = estimate_parameters_with_parmest_v24(
                experiments=[dataset],
                options=model_options,
                guess=self._guess_from_theta_dict(dataset, model_options, theta_values),
                calc_cov=False,
                solver="ipopt",
                solver_options=None,
                tee=False,
            )
            rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "objective": float(fit.get("objective", np.nan)),
                }
            )
        return rows

    @staticmethod
    def _guess_from_theta_dict(
        dataset: ExperimentalData,
        model_options: ModelOptions,
        theta_values: Dict[str, float],
    ) -> ParameterGuess:
        base = build_guess_from_experiment_v24(dataset, model_options)
        return DiafiltrationExperiment._guess_with_theta_overrides(base, theta_values)

    @staticmethod
    def _coerce_theta_dict(theta_obj: object) -> Dict[str, float]:
        if isinstance(theta_obj, pd.Series):
            return {str(key): float(value) for key, value in theta_obj.items()}
        if isinstance(theta_obj, dict):
            out: Dict[str, float] = {}
            for key, value in theta_obj.items():
                try:
                    out[str(key)] = float(value)
                except Exception:
                    continue
            return out
        return {}

    @staticmethod
    def _plan_to_dict(plan: RegressionPlan) -> Dict[str, object]:
        return {
            "dataset_ids": list(plan.dataset_ids),
            "shared_parameters": list(plan.shared_parameters),
            "dataset_specific_parameters": list(plan.dataset_specific_parameters),
            "mode": plan.mode,
            "joint_objective_description": plan.joint_objective_description,
        }

    @staticmethod
    def _design_summary_to_dict(summary: DesignWorkflowSummary) -> Dict[str, object]:
        return {
            "workflow_type": summary.workflow_type,
            "design_goal": summary.design_goal,
            "score_type": summary.score_type,
            "candidate_model_count": summary.candidate_model_count,
            "recommended_design": dict(summary.recommended_design),
        }

    @classmethod
    def _extract_design_summary(cls, result: Optional[Dict[str, object]]) -> Optional[Dict[str, object]]:
        if result is None:
            return None
        summary = result.get("summary")
        if isinstance(summary, dict):
            return dict(summary)
        return {
            "workflow_type": str(result.get("workflow_type", "")),
            "design_goal": str(result.get("design_goal", "")),
            "score_type": str(result.get("score_type", "")),
            "candidate_model_count": int(result.get("candidate_model_count", 0)),
            "recommended_design": dict(result.get("recommended_design") or {}),
        }

    @staticmethod
    def _model_option_name(model_options: ModelOptions) -> str:
        profile = getattr(model_options.process_model_profile, "value", None) or str(model_options.process_model_profile)
        return f"{profile}:{model_options.b_form}"

    @staticmethod
    def _count_measurements(datasets: Sequence[ExperimentalData], model_options: ModelOptions) -> int:
        count = 0
        for dataset in datasets:
            guess = build_guess_from_experiment_v24(dataset, model_options)
            model = model_construct_for_parmest_v24(dataset, options=model_options, guess=guess)
            count += len(model.experiment_outputs)
        return max(count, 1)

    @staticmethod
    def _information_criteria(*, objective: float, n_measurements: int, n_parameters: int) -> Tuple[float, float, float]:
        n = max(int(n_measurements), 1)
        k = max(int(n_parameters), 1)
        sse = max(float(objective), 1e-12)
        aic = float(n * np.log(sse / n) + 2.0 * k)
        if n - k - 1 <= 0:
            aicc = float("inf")
        else:
            aicc = float(aic + (2.0 * k * (k + 1)) / (n - k - 1))
        bic = float(n * np.log(sse / n) + k * np.log(n))
        return aic, aicc, bic

    @staticmethod
    def _design_signature(model: pyo.ConcreteModel) -> np.ndarray:
        tau_points = list(model.tau)
        last_tau = tau_points[-1]
        values = []
        for n in model.n_vial:
            values.extend(
                [
                    float(pyo.value(model.cF[n, last_tau])),
                    float(pyo.value(model.cH[n, last_tau])),
                    float(pyo.value(model.mV[n, last_tau])),
                ]
            )
        return np.asarray(values, dtype=float)

    @staticmethod
    def _pairwise_signature_distance(signatures: Sequence[np.ndarray]) -> float:
        if len(signatures) <= 1:
            return 0.0
        total = 0.0
        comparisons = 0
        for idx, left in enumerate(signatures):
            for right in signatures[idx + 1 :]:
                total += float(np.linalg.norm(left - right))
                comparisons += 1
        return total / max(comparisons, 1)

    @staticmethod
    def _pairwise_signature_distances(
        signatures: Sequence[np.ndarray],
        *,
        model_names: Sequence[str],
    ) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        for idx, left in enumerate(signatures):
            for right_idx, right in enumerate(signatures[idx + 1 :], start=idx + 1):
                rows.append(
                    {
                        "left_model": str(model_names[idx]),
                        "right_model": str(model_names[right_idx]),
                        "distance": float(np.linalg.norm(left - right)),
                    }
                )
        return rows
