#!/usr/bin/env python3
"""Scaffold runner to reproduce DATA1/DATA2 published results with unified v24 code.

This script is intentionally staged:
1) Runs canonical MAT datasets through unified estimation.
2) Writes machine-readable summaries (JSON).
3) Builds side-by-side paper vs unified comparison CSVs (template-driven).
4) Provides paper-plot styling hooks for future figure parity.

Full per-figure/table reproduction logic should be added incrementally by target ID
from UnifiedFramework/DATA3/docs/validation/paper_target_matrix.md.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import asdict, dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

# Ensure repository root is importable when running from scripts/.
REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Ensure matplotlib/font cache paths are writable in restricted environments.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library import (
    ExperimentMode,
    ModelOptions,
    ParameterGuess,
    ProcessModelProfile,
    RunMode,
)
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.workflow.unified_workflow import (
    DataLoader,
    DatasetRequest,
    DiafiltrationExperiment,
    UQEngine,
)
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.workflow.output_artifacts import (
    extract_model_trajectories,
    legacy_style_objective,
)


@dataclass(frozen=True)
class ValidationTolerances:
    """Project-approved tolerances for DATA1/DATA2 reproduction."""

    param_rel_primary: float = 0.05
    param_rel_with_explanation: float = 0.07
    objective_rel: float = 0.05
    curve_nrmse_range: float = 0.05
    table_sig_digits: int = 3


@dataclass(frozen=True)
class DatasetSpec:
    """Defines one dataset run target."""

    dataset_id: str
    stage: str
    paper: str
    data_file: str
    selector: Optional[str] = None
    mode: ExperimentMode = ExperimentMode.DATA
    b_form: str = "single"
    use_sigma_logit_transform: bool = True
    nfe_override: Optional[int] = None
    fix_sigma_in_estimation: bool = False
    beta_0_bounds: Optional[Tuple[float, float]] = None
    beta_1_bounds: Optional[Tuple[float, float]] = None
    estimation_enabled: bool = True
    fixed_guess: Optional[Dict[str, float]] = None
    initial_guess: Optional[Dict[str, float]] = None
    use_multistart_for_xlsx: bool = True
    use_multistart_for_mat_legacy: bool = False
    multistart_iterations: int = 20
    process_model_profile: Optional[ProcessModelProfile] = None
    paper_profile: Optional[str] = None


CANONICAL_DATASETS: Tuple[DatasetSpec, ...] = (
    DatasetSpec(
        dataset_id="DATA1_511.12",
        stage="A",
        paper="DATA1",
        data_file="DATA1_matlab/data_library/data_stru-dataset511.12.mat",
        mode=ExperimentMode.DATA,
        b_form="single",
        use_sigma_logit_transform=False,
        estimation_enabled=True,
        use_multistart_for_mat_legacy=True,
        multistart_iterations=30,
        process_model_profile=ProcessModelProfile.DATA1,
        paper_profile="DATA1_PAPER",
    ),
    DatasetSpec(
        dataset_id="DATA2_270511.123",
        stage="C",
        paper="DATA2",
        data_file="DATA1_matlab/data_library/data_stru-dataset270511.123.mat",
        mode=ExperimentMode.LAG,
        b_form="convection",
        use_sigma_logit_transform=False,
        fix_sigma_in_estimation=True,
        beta_0_bounds=(0.5, 5.0),
        beta_1_bounds=(0.0, 0.1),
        estimation_enabled=True,
        initial_guess={
            "Lp": 11.113241068147595,
            "sigma": 1.0,
            "beta_0": 1.0649294788103598,
            "beta_1": 0.015171421224603474,
        },
        process_model_profile=ProcessModelProfile.DATA2,
        paper_profile="DATA2_PAPER",
    ),
    DatasetSpec(
        dataset_id="DATA2_270611.123",
        stage="C",
        paper="DATA2",
        data_file="DATA1_matlab/data_library/data_stru-dataset270611.123.mat",
        mode=ExperimentMode.LAG,
        b_form="convection",
        use_sigma_logit_transform=False,
        fix_sigma_in_estimation=True,
        beta_0_bounds=(0.5, 5.0),
        beta_1_bounds=(0.0, 0.1),
        estimation_enabled=False,
        fixed_guess={
            "Lp": 11.113241068147595,
            "sigma": 1.0,
            "beta_0": 1.0649294788103598,
            "beta_1": 0.015171421224603474,
        },
        initial_guess={
            "Lp": 11.113241068147595,
            "sigma": 1.0,
            "beta_0": 1.0649294788103598,
            "beta_1": 0.015171421224603474,
        },
        process_model_profile=ProcessModelProfile.DATA2,
        paper_profile="DATA2_PAPER",
    ),
)


DEFAULT_IPOPT_OPTIONS: Dict[str, object] = {
    "max_iter": 5000,
    "tol": 1e-8,
    "acceptable_tol": 1e-6,
}


def configure_paper_plot_style() -> None:
    """Set plotting defaults for paper-style visual parity (first scaffold)."""
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "legend.fontsize": 9,
            "lines.linewidth": 1.8,
            "figure.dpi": 300,
            "savefig.dpi": 600,
        }
    )


def rel_err(unified: float, reference: float) -> float:
    """Relative error with robust zero handling."""
    if reference == 0:
        return abs(unified - reference)
    return abs(unified - reference) / abs(reference)


def nrmse_vs_range(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute normalized RMSE against observed range."""
    if y_true.size == 0 or y_true.size != y_pred.size:
        return math.inf
    rng = float(np.nanmax(y_true) - np.nanmin(y_true))
    if rng == 0:
        return 0.0
    rmse = float(np.sqrt(np.nanmean((y_true - y_pred) ** 2)))
    return rmse / rng


def flatten_theta(theta_obj: object) -> Dict[str, float]:
    """Normalize theta output from parmest into a simple float dict."""
    if isinstance(theta_obj, pd.Series):
        return {str(k): float(v) for k, v in theta_obj.items()}
    if isinstance(theta_obj, dict):
        out: Dict[str, float] = {}
        for k, v in theta_obj.items():
            try:
                out[str(k)] = float(v)
            except Exception:
                continue
        return out
    return {}


def build_model_options(spec: DatasetSpec, *, nfe: int, parameter_mode: RunMode) -> ModelOptions:
    """Create one canonical ModelOptions object for a dataset spec."""
    options = ModelOptions(
        mode=spec.mode,
        parameter_treatment_mode=parameter_mode,
        b_form=spec.b_form,
        use_sigma_logit_transform=bool(spec.use_sigma_logit_transform) if parameter_mode == RunMode.ESTIMATION else False,
        fix_sigma_in_estimation=bool(spec.fix_sigma_in_estimation),
        nfe=int(nfe),
        use_multistart_for_xlsx=bool(spec.use_multistart_for_xlsx),
        use_multistart_for_mat_legacy=bool(spec.use_multistart_for_mat_legacy),
        multistart_iterations=int(spec.multistart_iterations),
        process_model_profile=spec.process_model_profile,
        paper_profile=spec.paper_profile,
    )
    if spec.beta_0_bounds is not None:
        options.beta_0_lb = float(spec.beta_0_bounds[0])
        options.beta_0_ub = float(spec.beta_0_bounds[1])
    if spec.beta_1_bounds is not None:
        options.beta_1_lb = float(spec.beta_1_bounds[0])
        options.beta_1_ub = float(spec.beta_1_bounds[1])
    return options


def run_dataset(
    engine: UQEngine,
    repo_root: Path,
    spec: DatasetSpec,
    *,
    solver: str,
    nfe: int,
    run_doe: bool,
    calc_cov: bool,
    skip_curve_metrics: bool,
    ipopt_max_cpu_time: Optional[float],
) -> Tuple[Dict[str, object], object, ModelOptions, Dict[str, object]]:
    """Run one dataset through unified estimation (+optional DoE)."""
    data_path = repo_root / spec.data_file
    loader_request = DatasetRequest(
        file_path=str(data_path),
        selector=spec.selector,
        specs=None,
        convert_to_concentration=False,
        plot=False,
        dataset_id=spec.dataset_id,
    )
    exp, (ok, issues) = engine.data_loader.load(loader_request)

    nfe_eff = int(spec.nfe_override) if spec.nfe_override is not None else int(nfe)
    options = build_model_options(spec, nfe=nfe_eff, parameter_mode=RunMode.ESTIMATION)

    result: Dict[str, object] = {
        "dataset_id": spec.dataset_id,
        "stage": spec.stage,
        "paper": spec.paper,
        "data_file": str(data_path),
        "load_ok": bool(ok),
        "issues": list(issues),
        "objective": None,
        "theta": {},
        "covariance_method": None,
        "covariance_warning": None,
        "workflow_class": "UQEngine",
        "process_model_profile": spec.process_model_profile.value if spec.process_model_profile is not None else None,
        "paper_profile": spec.paper_profile,
    }

    estimation_payload: Dict[str, object] = {"theta": dict(spec.fixed_guess or {})}
    if spec.estimation_enabled:
        solver_opts = dict(DEFAULT_IPOPT_OPTIONS)
        if ipopt_max_cpu_time is not None:
            solver_opts["max_cpu_time"] = float(ipopt_max_cpu_time)
        guess_obj = _build_guess_from_theta(dict(spec.initial_guess or {}), options=options) if spec.initial_guess else None
        est = engine.fit_parameters(
            [exp],
            model_options=options,
            guess=guess_obj,
            calc_cov=calc_cov,
            solver=solver,
            solver_options=solver_opts,
            tee=False,
        )
        estimation_payload = est
        result["objective"] = est.get("objective")
        result["theta"] = flatten_theta(est.get("theta"))
        result["covariance_method"] = est.get("covariance_method")
        result["covariance_warning"] = est.get("covariance_warning")
    else:
        result["theta"] = dict(spec.fixed_guess or {})

    # Curve-quality metrics are used as cross-validation targets (especially DATA2).
    curve_metrics = (
        {
            "nrmse.mass": float("nan"),
            "nrmse.cF": float("nan"),
            "nrmse.cV": float("nan"),
            "objective_legacy": float("nan"),
        }
        if bool(skip_curve_metrics)
        else _compute_curve_metrics(exp, spec, result["theta"])
    )
    result["curve_metrics"] = curve_metrics
    result["objective_legacy"] = curve_metrics.get("objective_legacy")

    if run_doe:
        doe = engine.analyze_doe_parameter_estimation(exp, options, solver_name="ipopt", tee=False)
        result["doe_d_opt_logdet"] = float(doe["d_opt_logdet"])
        result["doe_fim"] = np.asarray(doe["fim"], dtype=float).tolist()

    return result, exp, options, estimation_payload


def _build_guess_from_theta(theta: Dict[str, float], options: ModelOptions) -> ParameterGuess:
    """Build ParameterGuess from a flattened theta dictionary."""
    guess = ParameterGuess(
        Lp=float(theta.get("Lp", 11.0)),
        sigma=float(theta.get("sigma", 1.0)),
        B=float(theta.get("B", 0.5)) if "B" in theta else None,
        beta_0=float(theta.get("beta_0", 1.0)) if "beta_0" in theta else None,
        beta_1=float(theta.get("beta_1", 0.01)) if "beta_1" in theta else None,
        S0=float(theta.get("S0", 0.0)) if "S0" in theta else None,
        S=float(theta.get("S", 0.0)) if "S" in theta else None,
    )
    return guess


def _compute_curve_metrics(exp, spec: DatasetSpec, theta: Dict[str, float]) -> Dict[str, float]:
    """Compute dataset-level NRMSE metrics by simulating with provided theta."""
    sim_opts = build_model_options(spec, nfe=spec.nfe_override or 30, parameter_mode=RunMode.SIMULATION)
    guess = _build_guess_from_theta(theta, options=sim_opts)
    try:
        m = DiafiltrationExperiment(exp, sim_opts, guess=guess).build_simulation_model(
            theta_values=theta,
            solver_name="ipopt",
            solver_options=DEFAULT_IPOPT_OPTIONS,
            tee=False,
        )
    except Exception:
        return {
            "nrmse.mass": float("nan"),
            "nrmse.cF": float("nan"),
            "nrmse.cV": float("nan"),
            "objective_legacy": float("nan"),
        }
    sim = extract_model_trajectories(exp, m)

    mass_true: List[float] = []
    mass_pred: List[float] = []
    cf_true: List[float] = []
    cf_pred: List[float] = []
    cv_true: List[float] = []
    cv_pred: List[float] = []
    t_delay = float(exp.vials[0].time_s[0])

    for n, v in enumerate(exp.vials, start=1):
        t_meas_s = np.asarray(v.time_s, dtype=float)
        if v.mass_g is not None:
            y = np.asarray(v.mass_g, dtype=float)
            yp = np.interp(t_meas_s, sim[n]["time_s"] + t_delay, sim[n]["mV"])
            mass_true.extend(y.tolist())
            mass_pred.extend(yp.tolist())
        if v.retentate_signal is not None:
            y = np.asarray(v.retentate_signal, dtype=float)
            yp = np.interp(t_meas_s, sim[n]["time_s"] + t_delay, sim[n]["cF"])
            cf_true.extend(y.tolist())
            cf_pred.extend(yp.tolist())
        yv = float(v.cV_avg) if isinstance(v.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(yv):
            yp = float(np.interp(float(v.time_s[-1]), sim[n]["time_s"] + t_delay, sim[n]["cV"]))
            cv_true.append(yv)
            cv_pred.append(yp)

    return {
        "nrmse.mass": nrmse_vs_range(np.asarray(mass_true, dtype=float), np.asarray(mass_pred, dtype=float)),
        "nrmse.cF": nrmse_vs_range(np.asarray(cf_true, dtype=float), np.asarray(cf_pred, dtype=float)),
        "nrmse.cV": nrmse_vs_range(np.asarray(cv_true, dtype=float), np.asarray(cv_pred, dtype=float)),
        "objective_legacy": legacy_style_objective(exp, sim, t_delay=t_delay),
    }


def load_reference_values(path: Path) -> pd.DataFrame:
    """Load paper reference-value table for side-by-side comparison."""
    if not path.exists():
        return pd.DataFrame(
            columns=[
                "dataset_id",
                "metric",
                "paper_value",
                "paper_units",
                "target_id",
                "notes",
            ]
        )
    return pd.read_csv(path)


def build_unified_metric_table(run_rows: Iterable[Dict[str, object]]) -> pd.DataFrame:
    """Create long-format metric table from run results."""
    records: List[Dict[str, object]] = []
    for row in run_rows:
        ds = str(row["dataset_id"])
        obj = row.get("objective")
        if obj is not None:
            records.append({"dataset_id": ds, "metric": "objective", "unified_value": float(obj)})
        obj_legacy = row.get("objective_legacy")
        if obj_legacy is not None and not pd.isna(obj_legacy):
            records.append({"dataset_id": ds, "metric": "objective_legacy", "unified_value": float(obj_legacy)})
        theta = row.get("theta") or {}
        if isinstance(theta, dict):
            for name, value in theta.items():
                records.append({"dataset_id": ds, "metric": f"theta.{name}", "unified_value": float(value)})
        curve_metrics = row.get("curve_metrics") or {}
        if isinstance(curve_metrics, dict):
            for name, value in curve_metrics.items():
                if str(name) == "objective_legacy":
                    continue
                records.append({"dataset_id": ds, "metric": str(name), "unified_value": float(value)})
    return pd.DataFrame.from_records(records)


def build_side_by_side_table(
    references: pd.DataFrame,
    unified_metrics: pd.DataFrame,
    tolerances: ValidationTolerances,
) -> pd.DataFrame:
    """Merge paper and unified values and compute pass/fail columns."""
    merged = references.merge(unified_metrics, on=["dataset_id", "metric"], how="left")
    merged["param_rel_primary_limit"] = tolerances.param_rel_primary
    merged["param_rel_with_expl_limit"] = tolerances.param_rel_with_explanation
    merged["objective_rel_limit"] = tolerances.objective_rel

    def _row_rel_error(row: pd.Series) -> Optional[float]:
        if pd.isna(row.get("paper_value")) or pd.isna(row.get("unified_value")):
            return None
        return rel_err(float(row["unified_value"]), float(row["paper_value"]))

    merged["relative_error"] = merged.apply(_row_rel_error, axis=1)

    def _row_status(row: pd.Series) -> str:
        if pd.isna(row.get("paper_value")) or pd.isna(row.get("unified_value")):
            return "MISSING_VALUE"
        metric = str(row["metric"])
        err = float(row["relative_error"])
        if metric.startswith("nrmse."):
            threshold = float(row.get("paper_value")) if not pd.isna(row.get("paper_value")) else tolerances.curve_nrmse_range
            return "PASS" if float(row.get("unified_value")) <= threshold else "FAIL"
        if metric == "objective":
            return "PASS" if err <= tolerances.objective_rel else "FAIL"
        # treat all other numeric metrics as parameters for now
        if err <= tolerances.param_rel_primary:
            return "PASS"
        if err <= tolerances.param_rel_with_explanation:
            return "PASS_WITH_EXPLANATION"
        return "FAIL"

    merged["status"] = merged.apply(_row_status, axis=1)
    return merged


def parse_args() -> argparse.Namespace:
    """CLI argument parser."""
    parser = argparse.ArgumentParser(description="Reproduce DATA1/DATA2 paper results (scaffold).")
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[4],
        help="Repository root path.",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=Path("UnifiedFramework/DATA3/results/reproduction"),
        help="Output root for run artifacts.",
    )
    parser.add_argument(
        "--reference-csv",
        type=Path,
        default=Path("UnifiedFramework/DATA3/docs/validation/paper_reference_values_template.csv"),
        help="Paper reference values CSV used for side-by-side comparison.",
    )
    parser.add_argument("--solver", default="ipopt", help="Estimator solver.")
    parser.add_argument("--nfe", type=int, default=30, help="Finite elements for DAE discretization.")
    parser.add_argument("--run-doe", action="store_true", help="Run DoE/FIM for each dataset.")
    parser.add_argument(
        "--calc-cov",
        action="store_true",
        help="Attempt covariance estimation (may require additional solver binaries).",
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        default=[d.dataset_id for d in CANONICAL_DATASETS],
        help="Subset of dataset IDs to run.",
    )
    parser.add_argument(
        "--mat-legacy-restarts",
        type=int,
        default=None,
        help="Override MAT legacy multistart iterations for selected datasets.",
    )
    parser.add_argument(
        "--skip-stage-artifacts",
        action="store_true",
        help="Skip Stage A/B artifact generation and run tolerance tables only.",
    )
    parser.add_argument(
        "--skip-curve-metrics",
        action="store_true",
        help="Skip simulation-based curve metric calculation for faster tolerance reruns.",
    )
    parser.add_argument(
        "--ipopt-max-cpu-time",
        type=float,
        default=None,
        help="Optional IPOPT max_cpu_time (seconds) applied to estimation solves.",
    )
    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()
    configure_paper_plot_style()
    engine = UQEngine(data_loader=DataLoader())

    run_id = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = (args.repo_root / args.out_root / run_id).resolve()
    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    summaries_dir = out_dir / "summaries"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    summaries_dir.mkdir(parents=True, exist_ok=True)

    selected = [d for d in CANONICAL_DATASETS if d.dataset_id in set(args.datasets)]
    if args.mat_legacy_restarts is not None:
        selected = [
            replace(d, multistart_iterations=int(args.mat_legacy_restarts))
            if d.data_file.endswith(".mat")
            else d
            for d in selected
        ]
    run_rows: List[Dict[str, object]] = []
    dataset_run_cache: Dict[str, Dict[str, object]] = {}
    for spec in selected:
        row, exp, options, est = run_dataset(
            engine=engine,
            repo_root=args.repo_root.resolve(),
            spec=spec,
            solver=args.solver,
            nfe=args.nfe,
            run_doe=bool(args.run_doe),
            calc_cov=bool(args.calc_cov),
            skip_curve_metrics=bool(args.skip_curve_metrics),
            ipopt_max_cpu_time=args.ipopt_max_cpu_time,
        )
        run_rows.append(row)
        dataset_run_cache[spec.dataset_id] = {
            "exp": exp,
            "options": options,
            "estimation": est,
            "spec": spec,
        }
        with (summaries_dir / f"{spec.dataset_id}.json").open("w") as f:
            json.dump(row, f, indent=2, default=str)

    unified_metrics = build_unified_metric_table(run_rows)
    unified_metrics.to_csv(tables_dir / "unified_metrics.csv", index=False)

    references = load_reference_values(args.repo_root / args.reference_csv)
    side_by_side = build_side_by_side_table(references, unified_metrics, ValidationTolerances())
    side_by_side.to_csv(tables_dir / "paper_vs_unified_side_by_side.csv", index=False)

    # Stage A concrete outputs: DATA1 main (Fig.2 + Table 1 + Fig.3 + Table 2).
    stage_a_artifacts: Dict[str, object] = {}
    stage_b_artifacts: Dict[str, object] = {}
    if not bool(args.skip_stage_artifacts):
        data1_rows = [r for r in run_rows if r.get("dataset_id") == "DATA1_511.12"]
        if data1_rows:
            data1_cache = dataset_run_cache["DATA1_511.12"]
            exp_data1 = data1_cache["exp"]
            est_data1 = data1_cache["estimation"]
            sim_options = build_model_options(
                data1_cache["spec"],
                nfe=30,
                parameter_mode=RunMode.SIMULATION,
            )
            artifacts = engine.build_data1_validation_artifacts(
                out_dir=out_dir,
                dataset=exp_data1,
                estimation_result=est_data1,
                simulation_options=sim_options,
                references=references,
                tolerances=ValidationTolerances(),
                solver_name=args.solver,
                solver_options=DEFAULT_IPOPT_OPTIONS,
                tee=False,
            )
            stage_a_artifacts = artifacts["stage_a"]
            stage_b_artifacts = artifacts["stage_b"]

    run_meta = {
        "run_id": run_id,
        "repo_root": str(args.repo_root.resolve()),
        "solver": args.solver,
        "nfe": args.nfe,
        "run_doe": bool(args.run_doe),
        "calc_cov": bool(args.calc_cov),
        "skip_stage_artifacts": bool(args.skip_stage_artifacts),
        "skip_curve_metrics": bool(args.skip_curve_metrics),
        "ipopt_max_cpu_time": args.ipopt_max_cpu_time,
        "mat_legacy_restarts_override": args.mat_legacy_restarts,
        "datasets": [d.dataset_id for d in selected],
        "tolerances": asdict(ValidationTolerances()),
        "notes": [
            "Figure/table-specific generation functions are scaffold placeholders in this phase.",
            "Canonical DATA1/DATA2 nightly reproduction now routes load/estimation/DoE through UQEngine.",
            "Use UnifiedFramework/DATA3/docs/validation/paper_target_matrix.md target IDs to incrementally implement plotting and table pipelines.",
        ],
        "stage_a_artifacts": stage_a_artifacts,
        "stage_b_artifacts": stage_b_artifacts,
    }
    with (out_dir / "run_metadata.json").open("w") as f:
        json.dump(run_meta, f, indent=2)

    print(f"[reproduce] run_id={run_id}")
    print(f"[reproduce] output={out_dir}")
    print(f"[reproduce] datasets={', '.join(run_meta['datasets'])}")


if __name__ == "__main__":
    main()
