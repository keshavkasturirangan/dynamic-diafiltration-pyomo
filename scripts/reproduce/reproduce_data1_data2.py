#!/usr/bin/env python3
"""Scaffold runner to reproduce DATA1/DATA2 published results with unified v24 code.

This script is intentionally staged:
1) Runs canonical MAT datasets through unified estimation.
2) Writes machine-readable summaries (JSON).
3) Builds side-by-side paper vs unified comparison CSVs (template-driven).
4) Provides paper-plot styling hooks for future figure parity.

Full per-figure/table reproduction logic should be added incrementally by target ID
from docs/validation/paper_target_matrix.md.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

# Ensure repository root is importable when running from scripts/.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Ensure matplotlib/font cache paths are writable in restricted environments.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

from UnifiedFramework.DATA3.ExperimentalDataLoader.UnifiedCode.experiment_dataload_OOP_v24 import (
    DiafiltrationExperimentV24,
    ExperimentMode,
    ModelOptions,
    ParameterGuess,
    RunMode,
    apply_discretization,
    estimate_parameters_with_parmest_v24,
    load_experiment_easy,
    model_construct_inter_v24,
    run_doe_with_pyomo_v24,
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


CANONICAL_DATASETS: Tuple[DatasetSpec, ...] = (
    DatasetSpec(
        dataset_id="DATA1_511.12",
        stage="A",
        paper="DATA1",
        data_file="DATA1_matlab/data_library/data_stru-dataset511.12.mat",
    ),
    DatasetSpec(
        dataset_id="DATA2_270511.123",
        stage="C",
        paper="DATA2",
        data_file="DATA1_matlab/data_library/data_stru-dataset270511.123.mat",
    ),
    DatasetSpec(
        dataset_id="DATA2_270611.123",
        stage="C",
        paper="DATA2",
        data_file="DATA1_matlab/data_library/data_stru-dataset270611.123.mat",
    ),
)


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


def run_dataset(
    repo_root: Path,
    spec: DatasetSpec,
    *,
    solver: str,
    nfe: int,
    run_doe: bool,
    calc_cov: bool,
) -> Dict[str, object]:
    """Run one dataset through unified estimation (+optional DoE)."""
    data_path = repo_root / spec.data_file
    exp, (ok, issues) = load_experiment_easy(
        str(data_path),
        selector=spec.selector,
        specs=None,
        convert_to_concentration=False,
        plot=False,
    )

    options = ModelOptions(
        mode=ExperimentMode.DATA,
        run_mode=RunMode.ESTIMATION,
        b_form="single",
        nfe=nfe,
    )

    est = estimate_parameters_with_parmest_v24(
        experiments=[exp],
        options=options,
        calc_cov=calc_cov,
        solver=solver,
        tee=False,
    )

    result: Dict[str, object] = {
        "dataset_id": spec.dataset_id,
        "stage": spec.stage,
        "paper": spec.paper,
        "data_file": str(data_path),
        "load_ok": bool(ok),
        "issues": list(issues),
        "objective": est.get("objective"),
        "theta": flatten_theta(est.get("theta")),
        "covariance_method": est.get("covariance_method"),
        "covariance_warning": est.get("covariance_warning"),
    }

    if run_doe:
        exp_obj = DiafiltrationExperimentV24(exp=exp, options=options)
        doe = run_doe_with_pyomo_v24(exp_obj, solver_name="ipopt", tee=False)
        result["doe_d_opt_logdet"] = float(doe["d_opt_logdet"])
        result["doe_fim"] = np.asarray(doe["fim"], dtype=float).tolist()

    return result


def _extract_model_trajectories(exp, model) -> Dict[int, Dict[str, np.ndarray]]:
    """Extract per-vial trajectories in physical time from discretized model."""
    import pyomo.environ as pyo

    tau_vals = sorted(float(t) for t in list(model.tau))
    t_delay = float(exp.vials[0].time_s[0])
    out: Dict[int, Dict[str, np.ndarray]] = {}
    for n in range(1, len(exp.vials) + 1):
        vial = exp.vials[n - 1]
        ti = float(vial.time_s[0]) - t_delay
        tf = float(vial.time_s[-1]) - t_delay
        dur = tf - ti
        t_model = np.array([ti + dur * tau for tau in tau_vals], dtype=float)
        mV = np.array([float(pyo.value(model.mV[n, tau])) for tau in tau_vals], dtype=float)
        cF = np.array([float(pyo.value(model.cF[n, tau])) for tau in tau_vals], dtype=float)
        cH = np.array([float(pyo.value(model.cH[n, tau])) for tau in tau_vals], dtype=float)
        cV = np.array([float(pyo.value(model.cV[n, tau])) for tau in tau_vals], dtype=float)
        out[n] = {"time_s": t_model, "mV": mV, "cF": cF, "cH": cH, "cV": cV}
    return out


def _build_data1_stage_a_outputs(
    repo_root: Path,
    out_dir: Path,
    exp,
    est: Dict[str, object],
    tolerances: ValidationTolerances,
    references: pd.DataFrame,
) -> Dict[str, object]:
    """Generate Stage A concrete outputs: DATA1 Fig.2, Fig.3, Table 1, Table 2."""
    import matplotlib.pyplot as plt
    import pyomo.environ as pyo

    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    t_delay = float(exp.vials[0].time_s[0])
    n_vials = len(exp.vials)

    # Build simulation model at estimated theta for overlay curves/metrics.
    theta = flatten_theta(est.get("theta"))
    lp = float(theta.get("Lp", np.nan))
    b = float(theta.get("B", np.nan))
    sigma = float(theta.get("sigma", np.nan))
    options_sim = ModelOptions(
        mode=ExperimentMode.DATA,
        run_mode=RunMode.SIMULATION,
        b_form="single",
        nfe=30,
    )
    guess = ParameterGuess(Lp=lp, B=b, sigma=sigma)
    m = model_construct_inter_v24(exp, options_sim, guess=guess)
    apply_discretization(m, nfe=options_sim.nfe, scheme=options_sim.fd_scheme)
    m.obj = pyo.Objective(expr=0.0)
    pyo.SolverFactory("ipopt").solve(m, tee=False)
    sim = _extract_model_trajectories(exp, m)

    # Assemble global measured/predicted arrays for metrics.
    mass_true: List[float] = []
    mass_pred: List[float] = []
    cf_true: List[float] = []
    cf_pred: List[float] = []
    cv_true: List[float] = []
    cv_pred: List[float] = []

    # =========================
    # Fig. 2 (all four plots)
    # =========================
    # 2a: measured vial mass vs time (per-vial scatter)
    fig2a = plt.figure(figsize=(4, 4))
    for i, v in enumerate(exp.vials, start=1):
        t = (np.asarray(v.time_s, dtype=float) - t_delay) / 60.0
        y = np.asarray(v.mass_g, dtype=float) if v.mass_g is not None else np.array([])
        if t.size and y.size:
            plt.plot(t, y, "r.", markersize=4, label="Measured" if i == 1 else None)
    plt.xlabel("Time [min]")
    plt.ylabel("Mass in Vial [g]")
    plt.ylim(bottom=0)
    plt.tick_params(direction="in")
    plt.legend(loc="best")
    fig2a_path = figures_dir / "data1_fig2_a_mass_measured.png"
    fig2a.savefig(fig2a_path, bbox_inches="tight")
    plt.close(fig2a)

    # 2b: measured retentate concentration
    fig2b = plt.figure(figsize=(4, 4))
    for i, v in enumerate(exp.vials, start=1):
        t = (np.asarray(v.time_s, dtype=float) - t_delay) / 60.0
        y = np.asarray(v.retentate_signal, dtype=float) if v.retentate_signal is not None else np.array([])
        if t.size and y.size:
            plt.plot(t, y, "ms", markersize=4, label="Measured cF" if i == 1 else None)
    plt.xlabel("Time [min]")
    plt.ylabel("Retentate concentration [mM]")
    plt.ylim(bottom=0)
    plt.tick_params(direction="in", top=True, right=True)
    plt.legend(loc="best")
    fig2b_path = figures_dir / "data1_fig2_b_retentate_measured.png"
    fig2b.savefig(fig2b_path, bbox_inches="tight")
    plt.close(fig2b)

    # 2c: measured permeate/vial concentration
    fig2c = plt.figure(figsize=(4, 4))
    for i, v in enumerate(exp.vials, start=1):
        tf = (float(v.time_s[-1]) - t_delay) / 60.0
        yv = float(v.cV_avg) if isinstance(v.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(yv):
            plt.plot([tf], [yv], "cs", markersize=6, label="Measured cV" if i == 1 else None)
    plt.xlabel("Time [min]")
    plt.ylabel("Vial concentration [mM]")
    plt.ylim(bottom=0)
    plt.tick_params(direction="in", top=True, right=True)
    plt.legend(loc="best")
    fig2c_path = figures_dir / "data1_fig2_c_vial_measured.png"
    fig2c.savefig(fig2c_path, bbox_inches="tight")
    plt.close(fig2c)

    # 2d: combined measured overview
    fig2d = plt.figure(figsize=(4, 4))
    for i, v in enumerate(exp.vials, start=1):
        t = (np.asarray(v.time_s, dtype=float) - t_delay) / 60.0
        if v.mass_g is not None:
            plt.plot(t, np.asarray(v.mass_g, dtype=float), "r.", markersize=3, alpha=0.6, label="Mass meas." if i == 1 else None)
        if v.retentate_signal is not None:
            plt.plot(t, np.asarray(v.retentate_signal, dtype=float), "m-", linewidth=1.2, alpha=0.6, label="cF meas." if i == 1 else None)
        yv = float(v.cV_avg) if isinstance(v.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(yv):
            plt.plot([(float(v.time_s[-1]) - t_delay) / 60.0], [yv], "c^", markersize=4, alpha=0.8, label="cV meas." if i == 1 else None)
    plt.xlabel("Time [min]")
    plt.ylabel("Measured responses")
    plt.tick_params(direction="in")
    plt.legend(loc="best")
    fig2d_path = figures_dir / "data1_fig2_d_combined_measured.png"
    fig2d.savefig(fig2d_path, bbox_inches="tight")
    plt.close(fig2d)

    # =========================
    # Fig. 3 (model overlay)
    # =========================
    fig3 = plt.figure(figsize=(4, 4))
    for n in range(1, n_vials + 1):
        v = exp.vials[n - 1]
        t_meas = (np.asarray(v.time_s, dtype=float) - t_delay) / 60.0
        t_sim = sim[n]["time_s"] / 60.0
        if v.mass_g is not None:
            y_mass = np.asarray(v.mass_g, dtype=float)
            plt.plot(t_meas, y_mass, "r.", markersize=3, alpha=0.6, label="Mass meas." if n == 1 else None)
            y_mass_pred = np.interp(t_meas * 60.0, sim[n]["time_s"], sim[n]["mV"])
            mass_true.extend(y_mass.tolist())
            mass_pred.extend(y_mass_pred.tolist())
        if v.retentate_signal is not None:
            y_cf = np.asarray(v.retentate_signal, dtype=float)
            plt.plot(t_meas, y_cf, "ms", markersize=3, alpha=0.6, label="cF meas." if n == 1 else None)
            plt.plot(t_sim, sim[n]["cF"], "g-", linewidth=1.4, alpha=0.8, label="cF pred." if n == 1 else None)
            y_cf_pred = np.interp(t_meas * 60.0, sim[n]["time_s"], sim[n]["cF"])
            cf_true.extend(y_cf.tolist())
            cf_pred.extend(y_cf_pred.tolist())
        y_cv = float(v.cV_avg) if isinstance(v.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(y_cv):
            t_end = float(v.time_s[-1]) - t_delay
            y_cv_pred = float(np.interp(t_end, sim[n]["time_s"], sim[n]["cV"]))
            cv_true.append(y_cv)
            cv_pred.append(y_cv_pred)
            plt.plot([t_end / 60.0], [y_cv], "cs", markersize=5, alpha=0.8, label="cV meas." if n == 1 else None)
            plt.plot([sim[n]["time_s"][-1] / 60.0], [sim[n]["cV"][-1]], "r^", markersize=5, alpha=0.8, label="cV pred." if n == 1 else None)
    plt.xlabel("Time [min]")
    plt.ylabel("Concentration / mass responses")
    plt.tick_params(direction="in", top=True, right=True)
    plt.legend(loc="best")
    fig3_path = figures_dir / "data1_fig3_model_overlay.png"
    fig3.savefig(fig3_path, bbox_inches="tight")
    plt.close(fig3)

    # =========================
    # Table 1 / Table 2
    # =========================
    objective = float(est.get("objective")) if est.get("objective") is not None else np.nan
    cov_method = est.get("covariance_method")
    cov_warn = est.get("covariance_warning")

    t1_metrics = pd.DataFrame(
        [
            {"dataset_id": "DATA1_511.12", "metric": "objective", "unified_value": objective, "target_id": "D1-M-T1"},
            {"dataset_id": "DATA1_511.12", "metric": "theta.Lp", "unified_value": lp, "target_id": "D1-M-T1"},
            {"dataset_id": "DATA1_511.12", "metric": "theta.B", "unified_value": b, "target_id": "D1-M-T1"},
            {"dataset_id": "DATA1_511.12", "metric": "theta.sigma", "unified_value": sigma, "target_id": "D1-M-T1"},
            {"dataset_id": "DATA1_511.12", "metric": "covariance.method", "unified_value": np.nan, "target_id": "D1-M-T1"},
        ]
    )
    t1_metrics.loc[t1_metrics["metric"] == "covariance.method", "unified_text"] = str(cov_method) if cov_method is not None else str(cov_warn)

    def _rmse(y_true: List[float], y_pred: List[float]) -> float:
        yt = np.asarray(y_true, dtype=float)
        yp = np.asarray(y_pred, dtype=float)
        if yt.size == 0 or yt.size != yp.size:
            return np.nan
        return float(np.sqrt(np.mean((yt - yp) ** 2)))

    rmse_mass = _rmse(mass_true, mass_pred)
    rmse_cf = _rmse(cf_true, cf_pred)
    rmse_cv = _rmse(cv_true, cv_pred)
    nrmse_mass = nrmse_vs_range(np.asarray(mass_true, dtype=float), np.asarray(mass_pred, dtype=float))
    nrmse_cf = nrmse_vs_range(np.asarray(cf_true, dtype=float), np.asarray(cf_pred, dtype=float))
    nrmse_cv = nrmse_vs_range(np.asarray(cv_true, dtype=float), np.asarray(cv_pred, dtype=float))

    t2_metrics = pd.DataFrame(
        [
            {"dataset_id": "DATA1_511.12", "metric": "rmse.mass", "unified_value": rmse_mass, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "rmse.cF", "unified_value": rmse_cf, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "rmse.cV", "unified_value": rmse_cv, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "nrmse.mass", "unified_value": nrmse_mass, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "nrmse.cF", "unified_value": nrmse_cf, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "nrmse.cV", "unified_value": nrmse_cv, "target_id": "D1-M-T2"},
            {
                "dataset_id": "DATA1_511.12",
                "metric": "nrmse.pass_threshold",
                "unified_value": tolerances.curve_nrmse_range,
                "target_id": "D1-M-T2",
            },
        ]
    )

    unified_stage_metrics = pd.concat([t1_metrics, t2_metrics], ignore_index=True)
    unified_stage_metrics.to_csv(tables_dir / "data1_stageA_unified_metrics.csv", index=False)

    refs = references.copy()
    refs_t1 = refs[(refs["dataset_id"] == "DATA1_511.12") & (refs["target_id"] == "D1-M-T1")].copy()
    refs_t2 = refs[(refs["dataset_id"] == "DATA1_511.12") & (refs["target_id"] == "D1-M-T2")].copy()
    t1_side = refs_t1.merge(t1_metrics[["dataset_id", "metric", "unified_value"]], on=["dataset_id", "metric"], how="outer")
    t2_side = refs_t2.merge(t2_metrics[["dataset_id", "metric", "unified_value"]], on=["dataset_id", "metric"], how="outer")

    for df in (t1_side, t2_side):
        def _rel(row):
            pv = row.get("paper_value")
            uv = row.get("unified_value")
            if pd.isna(pv) or pd.isna(uv):
                return np.nan
            return rel_err(float(uv), float(pv))
        df["relative_error"] = df.apply(_rel, axis=1)

    t1_side["status"] = t1_side["relative_error"].map(
        lambda e: "MISSING_VALUE"
        if pd.isna(e)
        else ("PASS" if e <= tolerances.objective_rel else "FAIL")
    )
    t2_side["status"] = t2_side["relative_error"].map(
        lambda e: "MISSING_VALUE"
        if pd.isna(e)
        else (
            "PASS"
            if e <= tolerances.param_rel_primary
            else ("PASS_WITH_EXPLANATION" if e <= tolerances.param_rel_with_explanation else "FAIL")
        )
    )
    t1_path = tables_dir / "data1_table1_side_by_side.csv"
    t2_path = tables_dir / "data1_table2_side_by_side.csv"
    t1_side.to_csv(t1_path, index=False)
    t2_side.to_csv(t2_path, index=False)

    return {
        "figures": [
            str(fig2a_path),
            str(fig2b_path),
            str(fig2c_path),
            str(fig2d_path),
            str(fig3_path),
        ],
        "tables": [
            str(t1_path),
            str(t2_path),
            str(tables_dir / "data1_stageA_unified_metrics.csv"),
        ],
        "summary": {
            "objective": objective,
            "theta": {"Lp": lp, "B": b, "sigma": sigma},
            "rmse": {"mass": rmse_mass, "cF": rmse_cf, "cV": rmse_cv},
            "nrmse": {"mass": nrmse_mass, "cF": nrmse_cf, "cV": nrmse_cv},
        },
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
        theta = row.get("theta") or {}
        if isinstance(theta, dict):
            for name, value in theta.items():
                records.append({"dataset_id": ds, "metric": f"theta.{name}", "unified_value": float(value)})
    return pd.DataFrame.from_records(records)


def build_side_by_side_table(
    references: pd.DataFrame,
    unified_metrics: pd.DataFrame,
    tolerances: ValidationTolerances,
) -> pd.DataFrame:
    """Merge paper and unified values and compute pass/fail columns."""
    merged = references.merge(unified_metrics, on=["dataset_id", "metric"], how="outer")
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
        default=Path(__file__).resolve().parents[2],
        help="Repository root path.",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=Path("results/reproduction"),
        help="Output root for run artifacts.",
    )
    parser.add_argument(
        "--reference-csv",
        type=Path,
        default=Path("docs/validation/paper_reference_values_template.csv"),
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
    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()
    configure_paper_plot_style()

    run_id = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = (args.repo_root / args.out_root / run_id).resolve()
    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    summaries_dir = out_dir / "summaries"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    summaries_dir.mkdir(parents=True, exist_ok=True)

    selected = [d for d in CANONICAL_DATASETS if d.dataset_id in set(args.datasets)]
    run_rows: List[Dict[str, object]] = []
    for spec in selected:
        row = run_dataset(
            repo_root=args.repo_root.resolve(),
            spec=spec,
            solver=args.solver,
            nfe=args.nfe,
            run_doe=bool(args.run_doe),
            calc_cov=bool(args.calc_cov),
        )
        run_rows.append(row)
        with (summaries_dir / f"{spec.dataset_id}.json").open("w") as f:
            json.dump(row, f, indent=2, default=str)

    unified_metrics = build_unified_metric_table(run_rows)
    unified_metrics.to_csv(tables_dir / "unified_metrics.csv", index=False)

    references = load_reference_values(args.repo_root / args.reference_csv)
    side_by_side = build_side_by_side_table(references, unified_metrics, ValidationTolerances())
    side_by_side.to_csv(tables_dir / "paper_vs_unified_side_by_side.csv", index=False)

    # Stage A concrete outputs: DATA1 main (Fig.2 + Table 1 + Fig.3 + Table 2).
    stage_a_artifacts: Dict[str, object] = {}
    data1_rows = [r for r in run_rows if r.get("dataset_id") == "DATA1_511.12"]
    if data1_rows:
        p_data1 = args.repo_root.resolve() / "DATA1_matlab/data_library/data_stru-dataset511.12.mat"
        exp_data1, _ = load_experiment_easy(
            str(p_data1),
            selector=None,
            specs=None,
            convert_to_concentration=False,
            plot=False,
        )
        est_data1 = estimate_parameters_with_parmest_v24(
            experiments=[exp_data1],
            options=ModelOptions(mode=ExperimentMode.DATA, run_mode=RunMode.ESTIMATION, b_form="single", nfe=30),
            calc_cov=False,
            solver=args.solver,
            tee=False,
        )
        stage_a_artifacts = _build_data1_stage_a_outputs(
            repo_root=args.repo_root.resolve(),
            out_dir=out_dir,
            exp=exp_data1,
            est=est_data1,
            tolerances=ValidationTolerances(),
            references=references,
        )

    run_meta = {
        "run_id": run_id,
        "repo_root": str(args.repo_root.resolve()),
        "solver": args.solver,
        "nfe": args.nfe,
        "run_doe": bool(args.run_doe),
        "calc_cov": bool(args.calc_cov),
        "datasets": [d.dataset_id for d in selected],
        "tolerances": asdict(ValidationTolerances()),
        "notes": [
            "Figure/table-specific generation functions are scaffold placeholders in this phase.",
            "Use docs/validation/paper_target_matrix.md target IDs to incrementally implement plotting and table pipelines.",
        ],
        "stage_a_artifacts": stage_a_artifacts,
    }
    with (out_dir / "run_metadata.json").open("w") as f:
        json.dump(run_meta, f, indent=2)

    print(f"[reproduce] run_id={run_id}")
    print(f"[reproduce] output={out_dir}")
    print(f"[reproduce] datasets={', '.join(run_meta['datasets'])}")


if __name__ == "__main__":
    main()
