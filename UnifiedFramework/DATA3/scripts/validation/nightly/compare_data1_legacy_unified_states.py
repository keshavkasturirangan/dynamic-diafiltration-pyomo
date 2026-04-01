#!/usr/bin/env python3
"""Compare DATA1 legacy and unified state trajectories at fixed legacy theta."""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyomo.dae import Simulator
import pyomo.environ as pyo
from pyomo.environ import value
from scipy.io import loadmat

REPO_ROOT = Path(__file__).resolve().parents[5]
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from utility import loadmat as legacy_loadmat  # noqa: E402
from utility import solve_model  # noqa: E402
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode import unified_codebase_runfile as runfile  # noqa: E402
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library import (  # noqa: E402
    ExperimentMode,
    ModelOptions,
    ParameterGuess,
    RunMode,
    apply_discretization,
    build_guess_from_experiment_v24,
    enforce_data1_boundary_initialization,
    load_experiment_easy,
    model_construct_inter_v24,
)

OUT_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "pytest_validation"
DATA1_ROOT = REPO_ROOT / "DATA1_matlab" / "data"


def _legacy_fit_params(fit_dir: str) -> dict[str, float]:
    fit_stru = loadmat(DATA1_ROOT / fit_dir / "fit_stru.mat", squeeze_me=True, struct_as_record=False)["fit_stru"]
    return {
        "Lp": float(fit_stru.Lp),
        "B": float(fit_stru.B),
        "sigma": float(fit_stru.sigma),
    }


def _run_legacy(dataset_id: str, theta: dict[str, float]) -> dict[int, dict[str, list[float]]]:
    data_stru = legacy_loadmat(REPO_ROOT / "DATA1_matlab" / "data_library" / f"data_stru-dataset{dataset_id}.mat")["data_stru"]
    data_stru["data_config"].setdefault("M_O", 0.0)
    data_stru["data_config"].setdefault("C_D", 0.0)
    data_stru["data_config"].setdefault("n_v0", 1)
    _fit_stru, sim_stru, _sim_inter = solve_model(
        data_stru,
        mode="DATA",
        theta=theta,
        sim_opt=True,
        B_form="single",
        LOUD=False,
    )
    return sim_stru


def _run_unified(dataset_file: Path, theta: dict[str, float], nfe: int = 300):
    exp, _notes = load_experiment_easy(path=dataset_file, convert_to_concentration=False)
    options = ModelOptions(
        mode=ExperimentMode.DATA,
        run_mode=RunMode.SIMULATION,
        b_form="single",
        nfe=nfe,
        use_sigma_logit_transform=False,
    )
    guess = build_guess_from_experiment_v24(
        exp,
        options,
        override=ParameterGuess(Lp=theta["Lp"], B=theta["B"], sigma=theta["sigma"]),
    )
    model = model_construct_inter_v24(exp=exp, options=options, guess=guess)
    sim = Simulator(model, package="casadi")
    sim.simulate(numpoints=max(300, nfe), integrator="idas")
    apply_discretization(model, nfe=nfe, scheme=options.fd_scheme)
    sim.initialize_model()
    enforce_data1_boundary_initialization(model)
    solver = pyo.SolverFactory("ipopt")
    solver.options["linear_solver"] = "ma97"
    solver.options["max_iter"] = 3000
    result = solver.solve(model, tee=False)
    term = result.solver.termination_condition
    if term != pyo.TerminationCondition.optimal:
        raise RuntimeError(f"Unified fixed-theta solve failed with termination_condition={term}")
    return model


def _extract_unified_states(model) -> dict[int, dict[str, list[float]]]:
    tau_vals = sorted(float(t) for t in list(model.tau))
    ti_vals = list(model.ti)
    tf_vals = list(model.tf)
    out: dict[int, dict[str, list[float]]] = {}
    for n in model.n_vial:
        vial = int(n)
        ti = float(ti_vals[vial - 1])
        tf = float(tf_vals[vial - 1])
        duration = tf - ti
        times = [ti + tau * duration for tau in tau_vals]
        mV = [float(value(model.mV[n, tau])) for tau in tau_vals]
        cV = [float(value(model.cV[n, tau])) for tau in tau_vals]
        out[vial - 1] = {
            "time": times,
            "cF": [float(value(model.cF[n, tau])) for tau in tau_vals],
            "cIn": [float(value(model.cIn[n, tau])) for tau in tau_vals],
            "cH": [float(value(model.cH[n, tau])) for tau in tau_vals],
            "mV": mV,
            "cV": cV,
            "cVmV": [float(value(model.cVmV[n, tau])) for tau in tau_vals],
            "Jw": [float(value(model.Jw[n, tau])) for tau in tau_vals],
            "Js": [float(value(model.Js[n, tau])) for tau in tau_vals],
            "B": [
                float(value(model.B[n, tau])) if hasattr(model.B, "index_set") and model.B.is_indexed() else float(value(model.B))
                for tau in tau_vals
            ],
            "cV_times_mV": [cV[i] * mV[i] for i in range(len(cV))],
        }
    return out


def _interp_compare(
    legacy: dict[int, dict[str, list[float]]],
    unified: dict[int, dict[str, list[float]]],
    dataset: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows: list[dict[str, object]] = []
    trace_rows: list[dict[str, object]] = []
    states = ["cF", "cH", "cV", "cIn", "mV", "Jw", "Js", "B"]

    for vial_key in sorted(legacy):
        legacy_v = legacy[vial_key]
        unified_v = unified[vial_key]
        vial = vial_key + 1

        legacy_time = [float(t) / 60.0 for t in legacy_v["time"]]
        unified_time = [float(t) / 60.0 for t in unified_v["time"]]

        for state in states:
            legacy_vals = pd.Series(legacy_v[state], dtype=float)
            unified_vals = pd.Series(unified_v[state], dtype=float)
            interp = pd.Series(np.interp(legacy_time, unified_time, unified_vals.to_numpy(dtype=float)), dtype=float)
            err = interp - legacy_vals
            metric_rows.append(
                {
                    "dataset": dataset,
                    "vial": vial,
                    "state": state,
                    "n_points": int(len(legacy_vals)),
                    "mae": float(err.abs().mean()),
                    "rmse": float(math.sqrt(float((err**2).mean()))),
                    "max_abs_err": float(err.abs().max()),
                    "legacy_start": float(legacy_vals.iloc[0]),
                    "legacy_end": float(legacy_vals.iloc[-1]),
                    "unified_interp_start": float(interp.iloc[0]),
                    "unified_interp_end": float(interp.iloc[-1]),
                }
            )
            threshold = max(1e-6, 0.05 * max(legacy_vals.abs().max(), 1.0))
            first_bad_idx = next((i for i, e in enumerate(err.abs()) if float(e) > threshold), None)
            for i, tmin in enumerate(legacy_time):
                trace_rows.append(
                    {
                        "dataset": dataset,
                        "vial": vial,
                        "state": state,
                        "time_min": float(tmin),
                        "legacy_value": float(legacy_vals.iloc[i]),
                        "unified_value": float(interp.iloc[i]),
                        "abs_err": float(abs(err.iloc[i])),
                        "first_divergence_time_min": float(legacy_time[first_bad_idx]) if first_bad_idx is not None else "",
                        "divergence_threshold": float(threshold),
                    }
                )

        legacy_cvmv = pd.Series(pd.Series(legacy_v["cV"], dtype=float) * pd.Series(legacy_v["mV"], dtype=float), dtype=float)
        unified_cvmv = pd.Series(unified_v["cVmV"], dtype=float)
        interp_cvmv = pd.Series(np.interp(legacy_time, unified_time, unified_cvmv.to_numpy(dtype=float)), dtype=float)
        err_cvmv = interp_cvmv - legacy_cvmv
        metric_rows.append(
            {
                "dataset": dataset,
                "vial": vial,
                "state": "cVmV",
                "n_points": int(len(legacy_cvmv)),
                "mae": float(err_cvmv.abs().mean()),
                "rmse": float(math.sqrt(float((err_cvmv**2).mean()))),
                "max_abs_err": float(err_cvmv.abs().max()),
                "legacy_start": float(legacy_cvmv.iloc[0]),
                "legacy_end": float(legacy_cvmv.iloc[-1]),
                "unified_interp_start": float(interp_cvmv.iloc[0]),
                "unified_interp_end": float(interp_cvmv.iloc[-1]),
            }
        )

    metrics = pd.DataFrame(metric_rows)
    traces = pd.DataFrame(trace_rows)
    return metrics, traces


def _boundary_compare(
    legacy: dict[int, dict[str, list[float]]],
    unified: dict[int, dict[str, list[float]]],
    dataset: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for vial_key in sorted(legacy):
        next_key = vial_key + 1
        if next_key not in legacy or next_key not in unified:
            continue

        lv = legacy[vial_key]
        ln = legacy[next_key]
        uv = unified[vial_key]
        un = unified[next_key]

        vial = vial_key + 1
        transitions = {
            "cF": (lv["cF"][-1], ln["cF"][0], uv["cF"][-1], un["cF"][0]),
            "cH": (lv["cH"][-1], ln["cH"][0], uv["cH"][-1], un["cH"][0]),
            "mV": (lv["mV"][-1], ln["mV"][0], uv["mV"][-1], un["mV"][0]),
            "cV": (lv["cV"][-1], ln["cV"][0], uv["cV"][-1], un["cV"][0]),
            "cVmV": (
                lv["cV"][-1] * lv["mV"][-1],
                ln["cV"][0] * ln["mV"][0],
                uv["cVmV"][-1],
                un["cVmV"][0],
            ),
        }
        for state, (legacy_end, legacy_next_start, unified_end, unified_next_start) in transitions.items():
            rows.append(
                {
                    "dataset": dataset,
                    "transition": f"{vial}->{vial + 1}",
                    "state": state,
                    "legacy_end": float(legacy_end),
                    "legacy_next_start": float(legacy_next_start),
                    "legacy_boundary_jump": float(legacy_next_start - legacy_end),
                    "unified_end": float(unified_end),
                    "unified_next_start": float(unified_next_start),
                    "unified_boundary_jump": float(unified_next_start - unified_end),
                    "boundary_jump_delta": float((unified_next_start - unified_end) - (legacy_next_start - legacy_end)),
                }
            )
    return pd.DataFrame(rows)


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    dataset_specs = [
        ("501.1", runfile.PRESET_DATA1_PATHS[0], "501.1 concpolar"),
        ("511.12", runfile.PRESET_DATA1_PATHS[3], "511.12 concpolar"),
    ]

    metric_frames: list[pd.DataFrame] = []
    trace_frames: list[pd.DataFrame] = []
    boundary_frames: list[pd.DataFrame] = []

    for dataset, dataset_file, fit_dir in dataset_specs:
        theta = _legacy_fit_params(fit_dir)
        legacy = _run_legacy(dataset, theta)
        unified_model = _run_unified(dataset_file, theta, nfe=300)
        unified = _extract_unified_states(unified_model)
        metrics, traces = _interp_compare(legacy, unified, dataset)
        boundaries = _boundary_compare(legacy, unified, dataset)
        metric_frames.append(metrics)
        trace_frames.append(traces)
        boundary_frames.append(boundaries)

    metrics_df = pd.concat(metric_frames, ignore_index=True)
    traces_df = pd.concat(trace_frames, ignore_index=True)
    boundary_df = pd.concat(boundary_frames, ignore_index=True)
    metrics_path = OUT_ROOT / "data1_legacy_vs_unified_state_metrics.csv"
    traces_path = OUT_ROOT / "data1_legacy_vs_unified_state_traces.csv"
    boundary_path = OUT_ROOT / "data1_legacy_vs_unified_boundary_report.csv"
    metrics_df.to_csv(metrics_path, index=False, quoting=csv.QUOTE_MINIMAL)
    traces_df.to_csv(traces_path, index=False, quoting=csv.QUOTE_MINIMAL)
    boundary_df.to_csv(boundary_path, index=False, quoting=csv.QUOTE_MINIMAL)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
