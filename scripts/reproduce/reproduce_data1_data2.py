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
from dataclasses import asdict, dataclass, replace
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

from UnifiedFramework.DATA3.ExperimentalDataLoader.UnifiedCode.unified_codebase_library import (
    DiafiltrationExperimentV24,
    ExperimentMode,
    ModelOptions,
    ParameterGuess,
    RunMode,
    build_guess_from_experiment_v24,
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


def run_dataset(
    repo_root: Path,
    spec: DatasetSpec,
    *,
    solver: str,
    nfe: int,
    run_doe: bool,
    calc_cov: bool,
    skip_curve_metrics: bool,
    ipopt_max_cpu_time: Optional[float],
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

    nfe_eff = int(spec.nfe_override) if spec.nfe_override is not None else int(nfe)
    options = ModelOptions(
        mode=spec.mode,
        run_mode=RunMode.ESTIMATION,
        b_form=spec.b_form,
        use_sigma_logit_transform=bool(spec.use_sigma_logit_transform),
        fix_sigma_in_estimation=bool(spec.fix_sigma_in_estimation),
        nfe=nfe_eff,
        use_multistart_for_xlsx=bool(spec.use_multistart_for_xlsx),
        use_multistart_for_mat_legacy=bool(spec.use_multistart_for_mat_legacy),
        multistart_iterations=int(spec.multistart_iterations),
    )
    if spec.beta_0_bounds is not None:
        options.beta_0_lb = float(spec.beta_0_bounds[0])
        options.beta_0_ub = float(spec.beta_0_bounds[1])
    if spec.beta_1_bounds is not None:
        options.beta_1_lb = float(spec.beta_1_bounds[0])
        options.beta_1_ub = float(spec.beta_1_bounds[1])

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
    }

    if spec.estimation_enabled:
        if spec.initial_guess:
            guess_obj = _build_guess_from_theta(dict(spec.initial_guess or {}), options=options)
        else:
            guess_obj = build_guess_from_experiment_v24(exp, options)
        solver_opts = dict(DEFAULT_IPOPT_OPTIONS)
        if ipopt_max_cpu_time is not None:
            solver_opts["max_cpu_time"] = float(ipopt_max_cpu_time)
        est = estimate_parameters_with_parmest_v24(
            experiments=[exp],
            options=options,
            guess=guess_obj,
            calc_cov=calc_cov,
            solver=solver,
            solver_options=solver_opts,
            tee=False,
        )
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
        else _compute_curve_metrics(exp, options, result["theta"])
    )
    result["curve_metrics"] = curve_metrics
    result["objective_legacy"] = curve_metrics.get("objective_legacy")

    if run_doe:
        exp_obj = DiafiltrationExperimentV24(exp=exp, options=options)
        doe = run_doe_with_pyomo_v24(exp_obj, solver_name="ipopt", tee=False)
        result["doe_d_opt_logdet"] = float(doe["d_opt_logdet"])
        result["doe_fim"] = np.asarray(doe["fim"], dtype=float).tolist()

    return result


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


def _compute_curve_metrics(exp, options: ModelOptions, theta: Dict[str, float]) -> Dict[str, float]:
    """Compute dataset-level NRMSE metrics by simulating with provided theta."""
    import pyomo.environ as pyo

    sim_opts = ModelOptions(
        mode=options.mode,
        run_mode=RunMode.SIMULATION,
        b_form=options.b_form,
        nfe=options.nfe,
        fd_scheme=options.fd_scheme,
        use_sigma_logit_transform=False,
        use_advanced_xlsx_transport_thermo=options.use_advanced_xlsx_transport_thermo,
    )
    guess = _build_guess_from_theta(theta, options=sim_opts)
    m = model_construct_inter_v24(exp, sim_opts, guess=guess)
    apply_discretization(m, nfe=sim_opts.nfe, scheme=sim_opts.fd_scheme)
    m.obj = pyo.Objective(expr=0.0)
    solver = pyo.SolverFactory("ipopt")
    for k, v in DEFAULT_IPOPT_OPTIONS.items():
        solver.options[str(k)] = v
    try:
        res = solver.solve(m, tee=False)
        term = getattr(getattr(res, "solver", None), "termination_condition", None)
        if term is not None and str(term).lower() not in {"optimal", "locallyoptimal"}:
            # Retry with higher iteration budget for difficult simulation NLPs.
            solver.options["max_iter"] = 10000
            solver.solve(m, tee=False)
    except Exception:
        return {
            "nrmse.mass": float("nan"),
            "nrmse.cF": float("nan"),
            "nrmse.cV": float("nan"),
            "objective_legacy": float("nan"),
        }
    sim = _extract_model_trajectories(exp, m)

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
        "objective_legacy": _legacy_style_objective(exp, sim, t_delay=t_delay),
    }


def _legacy_style_objective(exp, sim: Dict[int, Dict[str, np.ndarray]], *, t_delay: float) -> float:
    """Approximate legacy utility.py objective on unified simulated trajectories."""
    n_vials = len(exp.vials)
    n_extra = int(getattr(exp, "n_extra", 0) or 0)
    n_v0 = int(getattr(exp, "n_v0", 1) or 1)
    collect_vial = max(1, n_vials - n_extra)

    obj_m = 0.0
    obj_cp = 0.0
    obj_cf0 = 0.0
    obj_cf = 0.0
    cnt_m = 0
    cnt_cp = 0
    cnt_cf = 0
    cnt_cf0 = 0

    for n, vial in enumerate(exp.vials, start=1):
        t_meas_s = np.asarray(vial.time_s, dtype=float)

        # Mass: 0.01 g absolute error.
        if vial.mass_g is not None:
            y_m = np.asarray(vial.mass_g, dtype=float)
            y_m_pred = np.interp(t_meas_s, sim[n]["time_s"] + t_delay, sim[n]["mV"])
            mask = np.isfinite(y_m)
            if mask.any() and n >= n_v0:
                res = (y_m_pred[mask] - y_m[mask]) / 0.01
                obj_m += float(np.sum(res**2))
                cnt_m += int(np.sum(mask))

        # Permeate/vial concentration: 3% relative error.
        y_cv = float(vial.cV_avg) if isinstance(vial.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(y_cv) and n > n_extra and abs(y_cv) > 1e-12:
            y_cv_pred = float(np.interp(float(vial.time_s[-1]), sim[n]["time_s"] + t_delay, sim[n]["cV"]))
            obj_cp += float(((y_cv_pred - y_cv) / (0.03 * y_cv)) ** 2)
            cnt_cp += 1

        # Retentate concentration: 0.3% relative error.
        if vial.retentate_signal is not None:
            y_cf = np.asarray(vial.retentate_signal, dtype=float)
            y_cf_pred = np.interp(t_meas_s, sim[n]["time_s"] + t_delay, sim[n]["cF"])
            mask = np.isfinite(y_cf) & (np.abs(y_cf) > 1e-12)
            if mask.any():
                res = (y_cf_pred[mask] - y_cf[mask]) / (0.003 * y_cf[mask])
                if n >= n_v0:
                    obj_cf += float(np.sum(res**2))
                    cnt_cf += int(np.sum(mask))
                else:
                    obj_cf0 += float(np.sum(res**2))
                    cnt_cf0 += int(np.sum(mask))

    # Match legacy scaling: 1e4 * (avg mass + avg permeate + avg retentate)
    term_m = obj_m / max(1, cnt_m)
    term_cp = obj_cp / max(1, cnt_cp if cnt_cp > 0 else collect_vial)
    term_cf = (obj_cf0 + obj_cf) / max(1, cnt_cf0 + cnt_cf)
    return float(1e4 * (term_m + term_cp + term_cf))


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
        Jw = np.array([float(pyo.value(model.Jw[n, tau])) for tau in tau_vals], dtype=float)
        Js = np.array([float(pyo.value(model.Js[n, tau])) for tau in tau_vals], dtype=float)
        out[n] = {"time_s": t_model, "mV": mV, "cF": cF, "cH": cH, "cV": cV, "Jw": Jw, "Js": Js}
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


def _build_data1_stage_b_outputs(
    out_dir: Path,
    exp,
    est: Dict[str, object],
) -> Dict[str, object]:
    """Generate Stage B DATA1 artifacts (main Fig.4/5/6 + SI Fig.S2-S7)."""
    import matplotlib.pyplot as plt
    import pyomo.environ as pyo

    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

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

    t_delay = float(exp.vials[0].time_s[0])
    n_vials = len(exp.vials)

    produced: List[Tuple[str, str]] = []

    # DATA1 main Fig.4 (all plots): transport trajectories overview.
    fig4a = plt.figure(figsize=(5, 4))
    for n in range(1, n_vials + 1):
        plt.plot(sim[n]["time_s"] / 60.0, sim[n]["Jw"], linewidth=1.2, alpha=0.7)
    plt.xlabel("Time [min]")
    plt.ylabel("Jw [cm/s]")
    plt.tick_params(direction="in", top=True, right=True)
    p = figures_dir / "data1_fig4_1_jw_trajectories.png"
    fig4a.savefig(p, bbox_inches="tight")
    plt.close(fig4a)
    produced.append(("D1-M-F4", str(p)))

    fig4b = plt.figure(figsize=(5, 4))
    for n in range(1, n_vials + 1):
        plt.plot(sim[n]["time_s"] / 60.0, sim[n]["Js"], linewidth=1.2, alpha=0.7)
    plt.xlabel("Time [min]")
    plt.ylabel("Js [umol/cm2/s]")
    plt.tick_params(direction="in", top=True, right=True)
    p = figures_dir / "data1_fig4_2_js_trajectories.png"
    fig4b.savefig(p, bbox_inches="tight")
    plt.close(fig4b)
    produced.append(("D1-M-F4", str(p)))

    # DATA1 main Fig.5 (all plots): parity comparisons.
    mass_true: List[float] = []
    mass_pred: List[float] = []
    cf_true: List[float] = []
    cf_pred: List[float] = []
    cv_true: List[float] = []
    cv_pred: List[float] = []
    for n, v in enumerate(exp.vials, start=1):
        t_meas = np.asarray(v.time_s, dtype=float) - t_delay
        if v.mass_g is not None:
            y = np.asarray(v.mass_g, dtype=float)
            yp = np.interp(t_meas, sim[n]["time_s"], sim[n]["mV"])
            mass_true.extend(y.tolist())
            mass_pred.extend(yp.tolist())
        if v.retentate_signal is not None:
            y = np.asarray(v.retentate_signal, dtype=float)
            yp = np.interp(t_meas, sim[n]["time_s"], sim[n]["cF"])
            cf_true.extend(y.tolist())
            cf_pred.extend(yp.tolist())
        yv = float(v.cV_avg) if isinstance(v.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(yv):
            yp = float(np.interp(float(v.time_s[-1] - t_delay), sim[n]["time_s"], sim[n]["cV"]))
            cv_true.append(yv)
            cv_pred.append(yp)

    def _parity_plot(y_t: List[float], y_p: List[float], label: str, path: Path) -> None:
        yt = np.asarray(y_t, dtype=float)
        yp = np.asarray(y_p, dtype=float)
        mask = np.isfinite(yt) & np.isfinite(yp)
        yt = yt[mask]
        yp = yp[mask]
        fig = plt.figure(figsize=(4, 4))
        if yt.size:
            plt.plot(yt, yp, "o", markersize=3, alpha=0.75)
            lo = float(min(np.min(yt), np.min(yp)))
            hi = float(max(np.max(yt), np.max(yp)))
            plt.plot([lo, hi], [lo, hi], "k--", linewidth=1.0)
        plt.xlabel(f"{label} measured")
        plt.ylabel(f"{label} predicted")
        plt.tick_params(direction="in", top=True, right=True)
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)

    p = figures_dir / "data1_fig5_1_mass_parity.png"
    _parity_plot(mass_true, mass_pred, "mass", p)
    produced.append(("D1-M-F5", str(p)))
    p = figures_dir / "data1_fig5_2_cF_parity.png"
    _parity_plot(cf_true, cf_pred, "cF", p)
    produced.append(("D1-M-F5", str(p)))
    p = figures_dir / "data1_fig5_3_cV_parity.png"
    _parity_plot(cv_true, cv_pred, "cV", p)
    produced.append(("D1-M-F5", str(p)))

    # DATA1 main Fig.6 (all plots): residual profiles over experiment time.
    fig6 = plt.figure(figsize=(5, 4))
    for n, v in enumerate(exp.vials, start=1):
        if v.retentate_signal is None:
            continue
        t_meas = np.asarray(v.time_s, dtype=float) - t_delay
        y = np.asarray(v.retentate_signal, dtype=float)
        yp = np.interp(t_meas, sim[n]["time_s"], sim[n]["cF"])
        plt.plot(t_meas / 60.0, yp - y, "-", linewidth=1.0, alpha=0.7)
    plt.axhline(0.0, color="k", linestyle="--", linewidth=1.0)
    plt.xlabel("Time [min]")
    plt.ylabel("cF residual (pred - meas)")
    plt.tick_params(direction="in", top=True, right=True)
    p = figures_dir / "data1_fig6_1_cf_residuals.png"
    fig6.savefig(p, bbox_inches="tight")
    plt.close(fig6)
    produced.append(("D1-M-F6", str(p)))

    # DATA1 SI S2-S7: per-vial overlay bundles (6 groups).
    si_targets = ["D1-S-FS2", "D1-S-FS3", "D1-S-FS4", "D1-S-FS5", "D1-S-FS6", "D1-S-FS7"]
    group_size = max(1, int(np.ceil(n_vials / len(si_targets))))
    for i, target in enumerate(si_targets):
        start = i * group_size + 1
        end = min(n_vials, (i + 1) * group_size)
        fig = plt.figure(figsize=(5, 4))
        for n in range(start, end + 1):
            if n < 1 or n > n_vials:
                continue
            v = exp.vials[n - 1]
            t = (np.asarray(v.time_s, dtype=float) - t_delay) / 60.0
            if v.retentate_signal is not None:
                plt.plot(t, np.asarray(v.retentate_signal, dtype=float), ".", markersize=2, alpha=0.5)
                plt.plot(sim[n]["time_s"] / 60.0, sim[n]["cF"], "-", linewidth=1.0, alpha=0.7)
        plt.xlabel("Time [min]")
        plt.ylabel("Retentate concentration [mM]")
        plt.tick_params(direction="in", top=True, right=True)
        p = figures_dir / f"data1_si_figS{i+2}_1_overlay.png"
        fig.savefig(p, bbox_inches="tight")
        plt.close(fig)
        produced.append((target, str(p)))

    manifest = pd.DataFrame(produced, columns=["target_id", "artifact_path"])
    manifest_path = tables_dir / "data1_stageB_artifact_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    return {
        "manifest": str(manifest_path),
        "counts_by_target": manifest["target_id"].value_counts().to_dict(),
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
    for spec in selected:
        row = run_dataset(
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
            stage_b_artifacts = _build_data1_stage_b_outputs(
                out_dir=out_dir,
                exp=exp_data1,
                est=est_data1,
            )

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
            "Use docs/validation/paper_target_matrix.md target IDs to incrementally implement plotting and table pipelines.",
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
