#!/usr/bin/env python3
"""Build one-time simulation validation CSV files for DATA1 and DATA2 paper figures."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

import numpy as np
import pandas as pd
from pyomo.environ import (
    ConcreteModel,
    Expression,
    Objective,
    Param,
    RangeSet,
    Reals,
    SolverFactory,
    Var,
    exp,
    minimize,
    value,
)

REPO_ROOT = Path(__file__).resolve().parents[5]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from utility import loadmat as loadmat_data2
from utility import solve_model, solve_model_B_fix


OUT_ROOT = (
    REPO_ROOT
    / "UnifiedFramework"
    / "DATA3"
    / "docs"
    / "validation"
    / "simulation_validation_data_files"
)
DATA1_ROOT = REPO_ROOT / "DATA1_matlab" / "data"
DATA2_DATA_ROOT = REPO_ROOT / "data_library"


@dataclass
class ManifestRow:
    paper: str
    figure_id: str
    panel_id: str
    artifact_kind: str
    source_label: str
    relative_csv: str
    notes: str


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--finalize-only",
        action="store_true",
        help="Refresh README/metadata for an already-generated export without rerunning heavy solves.",
    )
    return parser.parse_args()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_csv(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    df.to_csv(path, index=False)


def build_manifest() -> list[ManifestRow]:
    return []


def loadmat_data1(path: Path) -> dict[str, Any]:
    data1_pkg = REPO_ROOT / "DATA1_matlab"
    if str(data1_pkg) not in sys.path:
        sys.path.insert(0, str(data1_pkg))
    from diafiltration_plots import loadmat  # pylint: disable=import-error

    return loadmat(str(path))


def normalise_series(rows: list[dict[str, Any]], value_columns: list[str]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    for col in value_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def export_data1_fig2(manifest: list[ManifestRow]) -> None:
    out_dir = ensure_dir(OUT_ROOT / "data1_main" / "fig2")
    panel_specs = [
        ("a", "501.1", "501.1 concpolar", "mass_g"),
        ("b", "501.1", "501.1 concpolar", "concentration_mM"),
        ("c", "511.12", "511.12 concpolar", "mass_g"),
        ("d", "511.12", "511.12 concpolar", "concentration_mM"),
    ]
    for panel_id, dataset, fit_dir, quantity in panel_specs:
        data_stru = loadmat_data1(DATA1_ROOT / f"data_stru-dataset{dataset}.mat")["data_stru"]
        fit_stru = loadmat_data1(DATA1_ROOT / fit_dir / "fit_stru.mat")["fit_stru"]
        t_delay = float(data_stru["data_raw"][0]["time"][0])
        rows: list[dict[str, Any]] = []

        if quantity == "mass_g":
            for vial_idx in range(data_stru["data_config"]["n"]):
                raw = data_stru["data_raw"][vial_idx]
                for t_val, y_val in zip(np.atleast_1d(raw["time"]), np.atleast_1d(raw["mass"])):
                    rows.append(
                        {
                            "series_id": f"measurement_vial_{vial_idx + 1}",
                            "series_group": "measurement",
                            "dataset": dataset,
                            "time_min": (float(t_val) - t_delay) / 60.0,
                            "mass_g": float(y_val),
                        }
                    )
                sim = fit_stru["sim_stru"][vial_idx]
                for t_val, y_val in zip(np.atleast_1d(sim["time"]), np.atleast_1d(sim["mV"])):
                    rows.append(
                        {
                            "series_id": f"prediction_vial_{vial_idx + 1}",
                            "series_group": "prediction",
                            "dataset": dataset,
                            "time_min": (float(t_val) - t_delay) / 60.0,
                            "mass_g": float(y_val),
                        }
                    )
        else:
            rows.append(
                {
                    "series_id": "retentate_prediction_vial_1",
                    "series_group": "prediction_line",
                    "dataset": dataset,
                    "time_min": (float(data_stru["data_raw"][0]["time"][0]) - t_delay) / 60.0,
                    "concentration_mM": float(fit_stru["sim_stru"][0]["cF"][0]),
                }
            )
            for vial_idx in range(data_stru["data_config"]["n"]):
                raw = data_stru["data_raw"][vial_idx]
                end_time_min = (float(np.atleast_1d(raw["time"])[-1]) - t_delay) / 60.0
                rows.append(
                    {
                        "series_id": f"permeate_measurement_vial_{vial_idx + 1}",
                        "series_group": "measurement_marker",
                        "dataset": dataset,
                        "time_min": end_time_min,
                        "concentration_mM": float(raw["cV_avg"]),
                    }
                )
                cf_exp = raw["cF_exp"]
                if isinstance(cf_exp, float):
                    rows.append(
                        {
                            "series_id": f"retentate_measurement_vial_{vial_idx + 1}",
                            "series_group": "measurement_marker",
                            "dataset": dataset,
                            "time_min": end_time_min,
                            "concentration_mM": float(cf_exp),
                        }
                    )
                else:
                    for t_val, y_val in zip(np.atleast_1d(raw["time"]), np.atleast_1d(cf_exp)):
                        rows.append(
                            {
                                "series_id": f"retentate_measurement_vial_{vial_idx + 1}",
                                "series_group": "measurement_marker",
                                "dataset": dataset,
                                "time_min": (float(t_val) - t_delay) / 60.0,
                                "concentration_mM": float(y_val),
                            }
                        )

                sim = fit_stru["sim_stru"][vial_idx]
                for t_val, y_val in zip(np.atleast_1d(sim["time"]), np.atleast_1d(sim["cF"])):
                    rows.append(
                        {
                            "series_id": f"retentate_prediction_vial_{vial_idx + 1}",
                            "series_group": "prediction_line",
                            "dataset": dataset,
                            "time_min": (float(t_val) - t_delay) / 60.0,
                            "concentration_mM": float(y_val),
                        }
                    )
                for t_val, y_val in zip(np.atleast_1d(sim["time"]), np.atleast_1d(sim["cH"])):
                    rows.append(
                        {
                            "series_id": f"permeate_prediction_vial_{vial_idx + 1}",
                            "series_group": "prediction_line",
                            "dataset": dataset,
                            "time_min": (float(t_val) - t_delay) / 60.0,
                            "concentration_mM": float(y_val),
                        }
                    )
                rows.append(
                    {
                        "series_id": f"vial_prediction_vial_{vial_idx + 1}",
                        "series_group": "prediction_marker",
                        "dataset": dataset,
                        "time_min": (float(np.atleast_1d(sim["time"])[-1]) - t_delay) / 60.0,
                        "concentration_mM": float(np.atleast_1d(sim["cV"])[-1]),
                    }
                )

        df = normalise_series(rows, ["time_min", quantity])
        out_path = out_dir / f"data1_main_fig2_{panel_id}.csv"
        write_csv(df, out_path)
        manifest.append(
            ManifestRow(
                paper="DATA1",
                figure_id="Fig. 2",
                panel_id=panel_id.upper(),
                artifact_kind="timeseries",
                source_label=f"DATA1_matlab data_stru + {fit_dir}/fit_stru",
                relative_csv=str(out_path.relative_to(OUT_ROOT)),
                notes=f"Dataset {dataset}; exported in long format with explicit series IDs.",
            )
        )


def export_data1_fig3(manifest: list[ManifestRow]) -> None:
    out_dir = ensure_dir(OUT_ROOT / "data1_main" / "fig3")
    rows: list[dict[str, Any]] = []

    for regime, filename, family_prefix in [
        ("filtration", "filtration.csv", "F"),
        ("diafiltration", "diafiltration.csv", "D"),
    ]:
        df = pd.read_csv(DATA1_ROOT / "experiment space" / filename, header=2)
        for col_idx in range(0, len(df.columns), 2):
            ret_col = df.columns[col_idx]
            perm_col = df.columns[col_idx + 1]
            family_id = family_prefix + str(col_idx // 2 + 1)
            pair = df[[ret_col, perm_col]].dropna()
            for _, row in pair.iterrows():
                rows.append(
                    {
                        "regime": regime,
                        "series_id": family_id,
                        "retentate_mM": float(row.iloc[0]),
                        "permeate_mM": float(row.iloc[1]),
                    }
                )
    out_path = out_dir / "data1_main_fig3.csv"
    write_csv(normalise_series(rows, ["retentate_mM", "permeate_mM"]), out_path)
    manifest.append(
        ManifestRow(
            paper="DATA1",
            figure_id="Fig. 3",
            panel_id="A",
            artifact_kind="scatter",
            source_label="DATA1_matlab experiment space CSVs",
            relative_csv=str(out_path.relative_to(OUT_ROOT)),
            notes="Combined filtration and diafiltration experiment-space points.",
        )
    )


def export_data1_fig4(manifest: list[ManifestRow]) -> None:
    out_dir = ensure_dir(OUT_ROOT / "data1_main" / "fig4")
    sigma_values = [0.1, 0.5, 0.9]
    panel_specs = [
        ("A", "501.1", "5.2843"),
        ("B", "511.12", "15.2052"),
    ]
    quantity_map = [
        ("mass", "mV", "mass_g"),
        ("retentate", "cF", "retentate_concentration_mM"),
        ("permeate", "cH", "permeate_concentration_mM"),
    ]

    for panel_id, dataset, c_fin in panel_specs:
        sim_bundle = {
            sigma: loadmat_data1(
                DATA1_ROOT / "sigma sensitivity" / f"sim_stru-dat{dataset} C_Fin{c_fin}sig{sigma:.1f}.mat"
            )["sim_stru"]
            for sigma in sigma_values
        }
        for name, field, value_col in quantity_map:
            rows: list[dict[str, Any]] = []
            for sigma, sim_stru in sim_bundle.items():
                t_delay = float(sim_stru[0]["time"][0])
                for vial_idx, vial in enumerate(sim_stru, start=1):
                    for t_val, y_val in zip(np.atleast_1d(vial["time"]), np.atleast_1d(vial[field])):
                        rows.append(
                            {
                                "series_id": f"sigma_{sigma:.1f}_vial_{vial_idx}",
                                "panel_id": panel_id,
                                "dataset": dataset,
                                "sigma": sigma,
                                "time_min": (float(t_val) - t_delay) / 60.0,
                                value_col: float(y_val),
                            }
                        )
            out_path = out_dir / f"data1_main_fig4_{panel_id.lower()}_{name}.csv"
            write_csv(normalise_series(rows, ["sigma", "time_min", value_col]), out_path)
            manifest.append(
                ManifestRow(
                    paper="DATA1",
                    figure_id="Fig. 4",
                    panel_id=f"{panel_id}-{name}",
                    artifact_kind="timeseries",
                    source_label="DATA1_matlab sigma sensitivity sim_stru mats",
                    relative_csv=str(out_path.relative_to(OUT_ROOT)),
                    notes=f"Dataset {dataset}; one row per vial trajectory and sigma value.",
                )
            )


def export_data1_contours(manifest: list[ManifestRow], figure_id: str, axis_key: str, suffix: str) -> None:
    out_dir = ensure_dir(OUT_ROOT / "data1_main" / suffix)
    # Main-paper contour figures use the diafiltration cases only.
    dataset_specs = [
        ("A", "511.12 concpolar"),
        ("B", "511.11 concpolar"),
        ("C", "511.12"),
    ]
    objective_columns = [
        ("mass", "Obj_mass"),
        ("permeate", "Obj_concentration"),
        ("retentate", "Obj_retentate_concentration"),
    ]
    for panel_id, dataset in dataset_specs:
        df = pd.read_csv(DATA1_ROOT / dataset / f"contourdata-x_{axis_key}-y_Lp.csv")
        axis_col = "B" if axis_key == "B" else "sigma"
        min_idx = {obj_name: int(df[obj_col].idxmin()) for obj_name, obj_col in objective_columns}
        for obj_name, obj_col in objective_columns:
            grid = df[[axis_col, "Lp", obj_col]].copy()
            grid.rename(columns={obj_col: "objective_value"}, inplace=True)
            grid["dataset"] = dataset
            grid["panel_id"] = panel_id
            grid["objective_name"] = obj_name
            optimum = grid.iloc[[min_idx[obj_name]]].copy()
            optimum["is_optimum"] = True
            grid["is_optimum"] = False
            out_path = out_dir / f"data1_main_{suffix}_{panel_id.lower()}_{obj_name}.csv"
            write_csv(pd.concat([grid, optimum], ignore_index=True), out_path)
            manifest.append(
                ManifestRow(
                    paper="DATA1",
                    figure_id=figure_id,
                    panel_id=f"{panel_id}-{obj_name}",
                    artifact_kind="contour_grid",
                    source_label=f"DATA1_matlab/{dataset} concpolar contourdata-x_{axis_key}-y_Lp.csv",
                    relative_csv=str(out_path.relative_to(OUT_ROOT)),
                    notes="Objective grid with optimum row flagged via is_optimum.",
                )
            )


def export_data2_pressure_and_calibration(manifest: list[ManifestRow]) -> None:
    out_dir = ensure_dir(OUT_ROOT / "data2_main" / "supporting")
    pressure_specs = [
        ("Fig. 2", "A", "NF270611.1_lag_pressure.csv", 185.0, "pressure_change_lag.csv"),
        ("Fig. 2", "B", "NF270511.4_overflow_pressure.csv", 130.0, "pressure_change_overflow.csv"),
    ]
    for figure_id, panel_id, filename, t0, out_name in pressure_specs:
        df = pd.read_csv(DATA2_DATA_ROOT / filename, skiprows=[1]).rename(columns=str.strip)
        export = pd.DataFrame(
            {
                "time_min": (pd.to_numeric(df["Time"]) - t0) / 60.0,
                "pressure_psi": pd.to_numeric(df["Pressure"]),
                "retentate_mM": pd.to_numeric(df["Concentration"]),
            }
        )
        out_path = out_dir / out_name
        write_csv(export, out_path)
        manifest.append(
            ManifestRow(
                paper="DATA2",
                figure_id=figure_id,
                panel_id=panel_id,
                artifact_kind="timeseries",
                source_label=f"DATA2_visualization.ipynb -> {filename}",
                relative_csv=str(out_path.relative_to(OUT_ROOT)),
                notes="Pressure and retentate concentration share the same x-axis.",
            )
        )

    for panel_id, filename, z_vals in [
        ("S1-A", "conductivity_calibration1.csv", [0.008813, -0.6949]),
        ("S1-B", "conductivity_calibration2.csv", [0.008372, -0.8735]),
    ]:
        df = pd.read_csv(DATA2_DATA_ROOT / filename, header=0, skiprows=[1]).rename(columns=str.strip)
        export = pd.DataFrame(
            {
                "conductivity_uS_per_cm": pd.to_numeric(df["Conductivity"]),
                "concentration_mM": pd.to_numeric(df["Concentration"]),
            }
        )
        export["trendline_concentration_mM"] = z_vals[0] * export["conductivity_uS_per_cm"] + z_vals[1]
        out_path = out_dir / f"data2_figs1_{panel_id.lower().replace('-', '_')}.csv"
        write_csv(export, out_path)
        manifest.append(
            ManifestRow(
                paper="DATA2",
                figure_id="Fig. S1",
                panel_id=panel_id,
                artifact_kind="scatter",
                source_label=f"DATA2_visualization.ipynb -> {filename}",
                relative_csv=str(out_path.relative_to(OUT_ROOT)),
                notes="Includes notebook linear trendline values.",
            )
        )


def export_data2_mass_tc(manifest: list[ManifestRow]) -> None:
    out_dir = ensure_dir(OUT_ROOT / "data2_main" / "fig3")
    data_stru = loadmat_data2(str(DATA2_DATA_ROOT / "data_stru-dataset270611.123.mat"))["data_stru"]
    raw = data_stru["data_raw"][3]
    times_min = np.array([float(t) / 60.0 for t in np.atleast_1d(raw["time"])])
    masses = np.array(np.atleast_1d(raw["mass"]), dtype=float)
    interp_times = np.linspace(140 / 60.0, float(raw["time"][0]) / 60.0, 50)
    interp_values = np.interp(interp_times, times_min, masses)
    rows = []
    for t_val, y_val in zip(times_min, masses):
        rows.append({"series_id": "measurements", "time_min": t_val, "mass_g": y_val})
    for t_val, y_val in zip(interp_times, interp_values):
        rows.append({"series_id": "extrapolation", "time_min": float(t_val), "mass_g": float(y_val)})
    rows.append({"series_id": "origin_marker", "time_min": float(raw["time"][0]) / 60.0, "mass_g": 0.0})
    rows.append({"series_id": "tc_marker", "time_min": 140 / 60.0, "mass_g": float(np.interp(140 / 60.0, times_min, masses))})
    out_path = out_dir / "data2_main_fig3_mass_tc.csv"
    write_csv(normalise_series(rows, ["time_min", "mass_g"]), out_path)
    manifest.append(
        ManifestRow(
            paper="DATA2",
            figure_id="Fig. 3",
            panel_id="A",
            artifact_kind="timeseries",
            source_label="DATA2_visualization.ipynb mass_tc construction",
            relative_csv=str(out_path.relative_to(OUT_ROOT)),
            notes="Measurement points plus extrapolation/marker overlays used in the notebook figure.",
        )
    )


def model_convection(sim_data: pd.DataFrame, pe_fixed_value: float | None = None) -> tuple[ConcreteModel, dict[str, float]]:
    pe_lower_bound = 0.1
    pe_initial_value = 15
    pe_upper_bound = 20

    model = ConcreteModel()
    model.i = RangeSet(sim_data.index[0] + 1, sim_data.index[-1] + 1)
    model.J_w = Param(model.i, initialize=lambda m, i: sim_data["Jw"][i - 1], mutable=True)
    model.J_s = Param(model.i, initialize=lambda m, i: sim_data["Js"][i - 1], mutable=True)
    model.c_in = Param(model.i, initialize=lambda m, i: sim_data["cIn"][i - 1], mutable=True)
    model.c_h = Param(model.i, initialize=lambda m, i: sim_data["cH"][i - 1], mutable=True)
    model.Js = Var(model.i, initialize=lambda m, i: sim_data["Jw"][i - 1])
    model.Kp = Var(model.i, initialize=1.0, bounds=(0.01, 1.2))
    model.Kf = Var(model.i, initialize=1.0, bounds=(0.01, 1.2))
    model.k1 = Var(initialize=1e-2)
    model.k0 = Var(initialize=1, bounds=(1e-4, 2))
    if pe_fixed_value is not None:
        pe_lower = min(pe_lower_bound, pe_fixed_value)
        pe_upper = max(pe_upper_bound, pe_fixed_value)
    else:
        pe_lower = pe_lower_bound
        pe_upper = pe_upper_bound
    model.Pe = Var(initialize=pe_initial_value, within=Reals, bounds=(pe_lower, pe_upper))
    model.expPe = Expression(expr=exp(model.Pe))
    if pe_fixed_value is not None:
        model.Pe.fix(pe_fixed_value)

    @model.Constraint(model.i)
    def partition_f(m, i):
        return m.Kf[i] == m.k1 * m.c_in[i] + m.k0

    @model.Constraint(model.i)
    def partition_p(m, i):
        return m.Kp[i] == m.k1 * m.c_h[i] + m.k0

    @model.Constraint(model.i)
    def convection(m, i):
        return m.Js[i] * (m.expPe - 1) == m.J_w[i] * (m.c_in[i] * m.Kf[i] * m.expPe - m.Kp[i] * m.c_h[i])

    model.FirstStageCost = Expression(rule=0)
    model.SecondStageCost = Expression(
        expr=sum((model.Js[i] / model.J_w[i] - model.J_s[i] / model.J_w[i]) ** 2 for i in model.i)
    )
    model.Total_Cost_Objective = Objective(expr=model.FirstStageCost + model.SecondStageCost, sense=minimize)

    solver = SolverFactory("ipopt")
    solver.options["linear_solver"] = "ma97"
    solver.options["halt_on_ampl_error"] = "yes"
    solver.options["acceptable_tol"] = 1e-8
    solver.solve(model, tee=False)
    theta_fit = {"k0": value(model.k0), "k1": value(model.k1), "Pe": value(model.Pe)}
    return model, theta_fit


def export_data2_fig7_fig8(manifest: list[ManifestRow]) -> None:
    out_dir7 = ensure_dir(OUT_ROOT / "data2_main" / "fig7")
    out_dir8 = ensure_dir(OUT_ROOT / "data2_main" / "fig8")
    sim_data = pd.read_csv(REPO_ROOT / "DATA1_matlab" / "data_library" / "sim_data-dat270611.123.csv")
    no_startup = sim_data.loc[sim_data["time"] > 120].reset_index(drop=True)
    model, theta_fit = model_convection(no_startup, pe_fixed_value=10)
    predicted_js = np.array([value(model.Js[i]) for i in model.i], dtype=float)
    predicted_kp = np.array([value(model.Kp[i]) for i in model.i], dtype=float)
    predicted_kf = np.array([value(model.Kf[i]) for i in model.i], dtype=float)

    fig7_exports = {
        "data2_main_fig7_js_vs_cin.csv": pd.DataFrame(
            {
                "cIn_mM": no_startup["cIn"].to_numpy(dtype=float),
                "empirical_Js_umol_cm2_s": no_startup["Js"].to_numpy(dtype=float),
                "predicted_Js_umol_cm2_s": predicted_js,
            }
        ),
        "data2_main_fig7_js_vs_time.csv": pd.DataFrame(
            {
                "time_min": no_startup["time"].to_numpy(dtype=float) / 60.0,
                "empirical_Js_umol_cm2_s": no_startup["Js"].to_numpy(dtype=float),
                "predicted_Js_umol_cm2_s": predicted_js,
            }
        ),
        "data2_main_fig7_jw_and_conc_vs_time.csv": pd.DataFrame(
            {
                "time_min": no_startup["time"].to_numpy(dtype=float) / 60.0,
                "Jw_um_per_s": no_startup["Jw"].to_numpy(dtype=float) * 1e4,
                "cIn_mM": no_startup["cIn"].to_numpy(dtype=float),
                "cH_mM": no_startup["cH"].to_numpy(dtype=float),
            }
        ),
        "data2_main_fig7_partition_profiles.csv": pd.DataFrame(
            {
                "cIn_mM": no_startup["cIn"].to_numpy(dtype=float),
                "predicted_Kp": predicted_kp,
                "predicted_Kf": predicted_kf,
                "fit_k0": theta_fit["k0"],
                "fit_k1": theta_fit["k1"],
                "fit_Pe": theta_fit["Pe"],
            }
        ),
    }
    for filename, df in fig7_exports.items():
        out_path = out_dir7 / filename
        write_csv(df, out_path)
        manifest.append(
            ManifestRow(
                paper="DATA2",
                figure_id="Fig. 7 / Fig. S8",
                panel_id=filename.replace("data2_main_fig7_", "").replace(".csv", ""),
                artifact_kind="timeseries",
                source_label="DATA2_visualization.ipynb convection-diffusion model",
                relative_csv=str(out_path.relative_to(OUT_ROOT)),
                notes="Uses sim_data-dat270611.123.csv and the notebook convection fit with Pe fixed at 10.",
            )
        )

    def linear_regression(data: pd.DataFrame, pe_value: float) -> tuple[float, float, float, float, float]:
        exp_pe = math.exp(pe_value)
        js = data["Js"].to_numpy(dtype=float)
        jw = data["Jw"].to_numpy(dtype=float)
        c_in = data["cIn"].to_numpy(dtype=float)
        c_h = data["cH"].to_numpy(dtype=float)
        y = js / jw
        x0 = (c_in * exp_pe - c_h) / (exp_pe - 1)
        x1 = (c_in**2 * exp_pe - c_h**2) / (exp_pe - 1)
        x = np.column_stack((x0, x1))
        beta, residuals, _, _ = np.linalg.lstsq(x, y, rcond=None)
        if residuals.size:
            sse = float(residuals[0])
        else:
            sse = float(np.sum((y - x @ beta) ** 2))
        n_obs, n_param = x.shape
        mse = sse / max(n_obs - n_param, 1)
        cov = mse * np.linalg.inv(x.T @ x)
        se = np.sqrt(np.diag(cov))
        return float(beta[0]), float(beta[1]), float(mse), float(se[0]), float(se[1])

    pe_values = np.logspace(-3, 2, 51)
    rows = []
    for pe_value in pe_values:
        k0, k1, mse, k0_se, k1_se = linear_regression(no_startup, float(pe_value))
        rows.append(
            {
                "Pe": float(pe_value),
                "mse": mse,
                "k0": k0,
                "k1": k1,
                "k0_minus_2se": k0 - 2 * k0_se,
                "k0_plus_2se": k0 + 2 * k0_se,
                "k1_minus_2se": k1 - 2 * k1_se,
                "k1_plus_2se": k1 + 2 * k1_se,
            }
        )
    out_path = out_dir8 / "data2_main_fig8_partition_sensitivity.csv"
    write_csv(
        normalise_series(
            rows,
            ["Pe", "mse", "k0", "k1", "k0_minus_2se", "k0_plus_2se", "k1_minus_2se", "k1_plus_2se"],
        ),
        out_path,
    )
    manifest.append(
        ManifestRow(
            paper="DATA2",
            figure_id="Fig. 8",
            panel_id="A-C",
            artifact_kind="parameter_sweep",
            source_label="DATA2_visualization.ipynb linear_regression sensitivity sweep",
            relative_csv=str(out_path.relative_to(OUT_ROOT)),
            notes="Single CSV contains the objective and both partition-coefficient profiles with confidence bands.",
        )
    )


def export_data2_fig9(manifest: list[ManifestRow]) -> None:
    out_dir = ensure_dir(OUT_ROOT / "data2_main" / "fig9")
    startup = pd.DataFrame({"Mode": ["Lag", "Overflow"], "Improvement_percent": [138, -9]})
    startup_path = out_dir / "data2_main_fig9_startup_improvement.csv"
    write_csv(startup, startup_path)
    manifest.append(
        ManifestRow(
            paper="DATA2",
            figure_id="Fig. 9",
            panel_id="A",
            artifact_kind="bar_values",
            source_label="DATA2_model_error_visualization.ipynb startup summary",
            relative_csv=str(startup_path.relative_to(OUT_ROOT)),
            notes="Values are the percentages hard-coded in the published-era notebook.",
        )
    )

    residual_specs = [
        (
            "concentrating",
            "data_stru-dataset270511.123.mat",
            ("Lag", None, "single"),
            ("Lag", {"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01, "sigma": 1.0, "S0": 0}, 1),
        ),
        (
            "diluting",
            "data_stru-dataset270511.821.mat",
            ("DATA", {"Lp": 10.106582659197427, "B": 16.418266360453412, "sigma": 1.0}, "single"),
            ("DATA", {"Lp": 10.167220041683919, "beta_0": 0.994192507828375, "beta_1": 0.013476905272029719, "sigma": 1.0}, 1),
        ),
    ]

    def residual_rows(fit_stru: dict[str, Any], label: str) -> list[dict[str, Any]]:
        rows_local = []
        for metric_name, values in [
            ("Mass", fit_stru["res_std"]["res_m"]),
            ("Permeate", fit_stru["res_std"]["res_cp"]),
            ("Retentate", fit_stru["res_std"]["res_cf"]),
        ]:
            arr = np.array(values, dtype=float)
            weighted = arr / math.sqrt(len(arr))
            for item in weighted:
                rows_local.append(
                    {
                        "solute_transport": label,
                        "residual_type": metric_name,
                        "weighted_residual": float(item),
                    }
                )
        return rows_local

    for regime, data_file, spec_a, spec_b in residual_specs:
        data_stru = loadmat_data2(str(DATA2_DATA_ROOT / data_file))["data_stru"]
        try:
            if regime == "diluting":
                fit_a, _, _ = solve_model_B_fix(data_stru, spec_a[0], theta=spec_a[1], sim_opt=False, B_form=spec_a[2], LOUD=False)
                fit_b, _, _ = solve_model_B_fix(data_stru, spec_b[0], theta=spec_b[1], sim_opt=False, B_form=spec_b[2], LOUD=False)
            else:
                if spec_a[2] == "single":
                    fit_a, _, _ = solve_model(data_stru, spec_a[0], theta=spec_a[1], sim_opt=False, B_form="single", LOUD=False)
                else:
                    fit_a, _, _ = solve_model(data_stru, spec_a[0], theta=spec_a[1], sim_opt=False, B_form=spec_a[2], LOUD=False)
                fit_b, _, _ = solve_model(data_stru, spec_b[0], theta=spec_b[1], sim_opt=False, B_form=spec_b[2], LOUD=False)
        except Exception as exc:  # pragma: no cover
            failure_path = out_dir / f"data2_main_fig9_{regime}_residuals_ERROR.txt"
            failure_path.write_text(str(exc), encoding="utf-8")
            manifest.append(
                ManifestRow(
                    paper="DATA2",
                    figure_id="Fig. 9",
                    panel_id=f"B-{regime}",
                    artifact_kind="error_note",
                    source_label="DATA2_model_error_visualization.ipynb residual boxplot inputs",
                    relative_csv=str(failure_path.relative_to(OUT_ROOT)),
                    notes="Residual export skipped because the local legacy Pyomo solve failed.",
                )
            )
            continue
        rows = residual_rows(fit_a, "Diffusion Only") + residual_rows(fit_b, "Convection-Diffusion")
        out_path = out_dir / f"data2_main_fig9_{regime}_residuals.csv"
        write_csv(normalise_series(rows, ["weighted_residual"]), out_path)
        manifest.append(
            ManifestRow(
                paper="DATA2",
                figure_id="Fig. 9",
                panel_id=f"B-{regime}",
                artifact_kind="residual_distribution",
                source_label="DATA2_model_error_visualization.ipynb residual boxplot inputs",
                relative_csv=str(out_path.relative_to(OUT_ROOT)),
                notes="Weighted residuals exported in long format, grouped by residual type and transport model.",
            )
        )


def flatten_data2_mass_and_concentration(data_stru: dict[str, Any], sim_stru: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    t_delay = float(data_stru["data_raw"][0]["time"][0])
    mass_rows: list[dict[str, Any]] = []
    conc_rows: list[dict[str, Any]] = []
    for vial_idx in range(data_stru["data_config"]["n"]):
        raw = data_stru["data_raw"][vial_idx]
        for t_val, y_val in zip(np.atleast_1d(raw["time"]), np.atleast_1d(raw["mass"])):
            mass_rows.append(
                {
                    "series_id": f"measurement_vial_{vial_idx + 1}",
                    "time_min": (float(t_val) - t_delay) / 60.0,
                    "mass_g": float(y_val),
                }
            )
        sim = sim_stru[vial_idx]
        for t_val, y_val in zip(np.atleast_1d(sim["time"]), np.atleast_1d(sim["mV"])):
            mass_rows.append(
                {
                    "series_id": f"prediction_vial_{vial_idx + 1}",
                    "time_min": float(t_val) / 60.0,
                    "mass_g": float(y_val),
                }
            )

        c_v_avg = raw["cV_avg"]
        if isinstance(c_v_avg, list):
            for t_val, y_val in zip(np.atleast_1d(raw["time"]), np.atleast_1d(c_v_avg)):
                conc_rows.append(
                    {
                        "series_id": f"permeate_measurement_vial_{vial_idx + 1}",
                        "time_min": (float(t_val) - t_delay) / 60.0,
                        "concentration_mM": float(y_val),
                    }
                )
        else:
            conc_rows.append(
                {
                    "series_id": f"permeate_measurement_vial_{vial_idx + 1}",
                    "time_min": (float(np.atleast_1d(raw["time"])[-1]) - t_delay) / 60.0,
                    "concentration_mM": float(c_v_avg),
                }
            )
        for t_val, y_val in zip(np.atleast_1d(raw["time"]), np.atleast_1d(raw["cF_exp"])):
            conc_rows.append(
                {
                    "series_id": f"retentate_measurement_vial_{vial_idx + 1}",
                    "time_min": (float(t_val) - t_delay) / 60.0,
                    "concentration_mM": float(y_val),
                }
            )
        for t_val, y_val in zip(np.atleast_1d(sim["time"]), np.atleast_1d(sim["cF"])):
            conc_rows.append(
                {
                    "series_id": f"retentate_prediction_vial_{vial_idx + 1}",
                    "time_min": float(t_val) / 60.0,
                    "concentration_mM": float(y_val),
                }
            )
        for t_val, y_val in zip(np.atleast_1d(sim["time"]), np.atleast_1d(sim["cH"])):
            conc_rows.append(
                {
                    "series_id": f"permeate_prediction_vial_{vial_idx + 1}",
                    "time_min": float(t_val) / 60.0,
                    "concentration_mM": float(y_val),
                }
            )
    return normalise_series(mass_rows, ["time_min", "mass_g"]), normalise_series(conc_rows, ["time_min", "concentration_mM"])


def export_data2_cross_verification(manifest: list[ManifestRow]) -> None:
    out_dir = ensure_dir(OUT_ROOT / "data2_si" / "cross_verification")
    fit_stru_base = {
        "Lp": 11.113241068147595,
        "beta_0": 1.0649294788103598,
        "beta_1": 0.015171421224603474,
        "sigma": 1.0,
    }
    cases = [
        ("FS2-A1", "data_stru-dataset270611.121.mat", True),
        ("FS2-A2", "data_stru-dataset270711.121.mat", True),
        ("FS3-B1", "data_stru-dataset270511.221.mat", True),
        ("FS3-B2", "data_stru-dataset270511.321.mat", False),
        ("FS4-C1", "data_stru-dataset270511.421.mat", False),
        ("FS4-C2", "data_stru-dataset270511.921.mat", True),
        ("FS5-D1", "data_stru-dataset270511.521.mat", True),
        ("FS5-D2", "data_stru-dataset270511.621.mat", True),
        ("FS6-E1", "data_stru-dataset270511.721.mat", False),
        ("FS6-E2", "data_stru-dataset270511.821.mat", True),
    ]
    for panel_id, data_file, sigma_fixed in cases:
        data_stru = loadmat_data2(str(DATA2_DATA_ROOT / data_file))["data_stru"]
        _, sim_stru, _ = solve_model_B_fix(
            data_stru,
            "DATA",
            theta=fit_stru_base,
            sim_opt=False,
            B_form=1,
            sigma_fixed=sigma_fixed,
            LOUD=False,
        )
        mass_df, conc_df = flatten_data2_mass_and_concentration(data_stru, sim_stru)
        dataset = str(data_stru["dataset"])
        mass_path = out_dir / f"data2_{panel_id.lower()}_mass_dat{dataset}.csv"
        conc_path = out_dir / f"data2_{panel_id.lower()}_concentration_dat{dataset}.csv"
        write_csv(mass_df, mass_path)
        write_csv(conc_df, conc_path)
        manifest.append(
            ManifestRow(
                paper="DATA2",
                figure_id=f"Fig. {panel_id.split('-')[0].replace('FS', 'S')}",
                panel_id=panel_id,
                artifact_kind="timeseries",
                source_label="run_cross_verification.py logic via utility.solve_model_B_fix",
                relative_csv=str(mass_path.relative_to(OUT_ROOT)),
                notes=f"Mass comparison for dataset {dataset}.",
            )
        )
        manifest.append(
            ManifestRow(
                paper="DATA2",
                figure_id=f"Fig. {panel_id.split('-')[0].replace('FS', 'S')}",
                panel_id=panel_id,
                artifact_kind="timeseries",
                source_label="run_cross_verification.py logic via utility.solve_model_B_fix",
                relative_csv=str(conc_path.relative_to(OUT_ROOT)),
                notes=f"Concentration comparison for dataset {dataset}.",
            )
        )


def export_data2_pre_b_dependence(manifest: list[ManifestRow]) -> None:
    out_dir = ensure_dir(OUT_ROOT / "data2_si" / "figs7")
    cases = [
        ("A.0", "data_stru-dataset270511.12.mat"),
        ("A.1", "data_stru-dataset270611.12.mat"),
        ("A.2", "data_stru-dataset270711.12.mat"),
        ("B.1", "data_stru-dataset270511.22.mat"),
        ("B.2", "data_stru-dataset270511.32.mat"),
        ("C.1", "data_stru-dataset270511.42.mat"),
        ("C.2", "data_stru-dataset270511.92.mat"),
        ("D.1", "data_stru-dataset270511.52.mat"),
        ("D.2", "data_stru-dataset270511.62.mat"),
        ("E.1", "data_stru-dataset270511.72.mat"),
        ("E.2", "data_stru-dataset270511.82.mat"),
    ]
    js_rows: list[dict[str, Any]] = []
    b_rows: list[dict[str, Any]] = []
    for label, data_file in cases:
        data_stru = loadmat_data2(str(DATA2_DATA_ROOT / data_file))["data_stru"]
        _, sim_stru, _ = solve_model(data_stru, "DATA", sim_opt=False, B_form="pervial", LOUD=False)
        for vial_idx, vial in sim_stru.items():
            if vial["B"] is None:
                continue
            start_idx = 80 if vial_idx == 0 else 50
            cin = np.array(vial["cIn"], dtype=float)
            js = np.array(vial["Js"], dtype=float)
            jw = np.array(vial["Jw"], dtype=float)
            b = np.array(vial["B"], dtype=float)
            for c_val, ratio in zip(cin[start_idx:], js[start_idx:] / jw[start_idx:]):
                js_rows.append(
                    {
                        "case_id": label,
                        "vial_index": vial_idx + 1,
                        "interface_concentration_mM": float(c_val),
                        "js_over_jw_mM": float(ratio),
                    }
                )
            for c_val, b_val in zip(cin, b):
                b_rows.append(
                    {
                        "case_id": label,
                        "vial_index": vial_idx + 1,
                        "interface_concentration_mM": float(c_val),
                        "B_um_per_s": float(b_val),
                    }
                )
    js_path = out_dir / "data2_figs7_js_over_jw_vs_cin.csv"
    b_path = out_dir / "data2_figs7_B_vs_cin.csv"
    write_csv(normalise_series(js_rows, ["interface_concentration_mM", "js_over_jw_mM"]), js_path)
    write_csv(normalise_series(b_rows, ["interface_concentration_mM", "B_um_per_s"]), b_path)
    manifest.append(
        ManifestRow(
            paper="DATA2",
            figure_id="Fig. S7",
            panel_id="A",
            artifact_kind="timeseries",
            source_label="run_pre_B_dependence.py logic via utility.solve_model(B_form='pervial')",
            relative_csv=str(b_path.relative_to(OUT_ROOT)),
            notes="Discrete-B curves grouped by case and vial.",
        )
    )
    manifest.append(
        ManifestRow(
            paper="DATA2",
            figure_id="Fig. S7",
            panel_id="B",
            artifact_kind="timeseries",
            source_label="run_pre_B_dependence.py logic via utility.solve_model(B_form='pervial')",
            relative_csv=str(js_path.relative_to(OUT_ROOT)),
            notes="Js/Jw versus interface concentration grouped by case and vial.",
        )
    )


def _write_readme(out_root: Path) -> None:
    readme = out_root / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# simulation_validation_data_files",
                "",
                "One-time generated simulation-side CSV files used to validate legacy paper figures.",
                "",
                "Source code paths:",
                "- `utility.py`",
                "- `DATA1_model_demo.ipynb`",
                "- `DATA2_model_demo.ipynb`",
                "- `DATA2_visualization.ipynb`",
                "- `run_DATA2_model_variations.py`",
                "",
                "These CSVs are intended to stay stable after generation so later pytest checks can compare unified-code outputs against the same published-paper baselines.",
                "",
                "This folder is the current source of truth for simulation-side paper-validation CSVs.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_metadata(errors_present: bool) -> None:
    metadata_json = OUT_ROOT / "simulation_validation_metadata.json"
    errors_path = OUT_ROOT / "simulation_validation_errors.json"
    metadata = {
        "generator": str(Path(__file__).relative_to(REPO_ROOT)),
        "output_root": str(OUT_ROOT.relative_to(REPO_ROOT)),
        "manifest_csv": str((OUT_ROOT / "simulation_validation_manifest.csv").relative_to(REPO_ROOT)),
        "manifest_json": str((OUT_ROOT / "simulation_validation_manifest.json").relative_to(REPO_ROOT)),
        "errors_file": (
            str(errors_path.relative_to(REPO_ROOT))
            if errors_present
            else None
        ),
    }
    metadata_json.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    _write_readme(OUT_ROOT)


def _finalize_existing_export() -> None:
    errors_path = OUT_ROOT / "simulation_validation_errors.json"
    _write_metadata(errors_path.exists())


def main() -> int:
    args = _parse_args()
    ensure_dir(OUT_ROOT)

    if args.finalize_only:
        _finalize_existing_export()
        return 0

    manifest: list[ManifestRow] = build_manifest()
    errors: list[dict[str, str]] = []

    def run_step(label: str, func, *func_args) -> None:
        try:
            func(*func_args)
        except Exception as exc:  # pragma: no cover
            errors.append({"step": label, "error": str(exc)})

    run_step("data1_fig2", export_data1_fig2, manifest)
    run_step("data1_fig3", export_data1_fig3, manifest)
    run_step("data1_fig4", export_data1_fig4, manifest)
    run_step("data1_fig5", export_data1_contours, manifest, "Fig. 5", "sigma", "fig5")
    run_step("data1_fig6", export_data1_contours, manifest, "Fig. 6", "B", "fig6")
    run_step("data2_supporting", export_data2_pressure_and_calibration, manifest)
    run_step("data2_fig3", export_data2_mass_tc, manifest)
    run_step("data2_fig7_fig8", export_data2_fig7_fig8, manifest)
    run_step("data2_fig9", export_data2_fig9, manifest)
    run_step("data2_cross_verification", export_data2_cross_verification, manifest)
    run_step("data2_pre_b_dependence", export_data2_pre_b_dependence, manifest)

    manifest_rows = pd.DataFrame([asdict(row) for row in manifest])
    write_csv(manifest_rows, OUT_ROOT / "simulation_validation_manifest.csv")
    (OUT_ROOT / "simulation_validation_manifest.json").write_text(
        json.dumps([asdict(row) for row in manifest], indent=2),
        encoding="utf-8",
    )
    if errors:
        (OUT_ROOT / "simulation_validation_errors.json").write_text(
            json.dumps(errors, indent=2),
            encoding="utf-8",
        )
    _write_metadata(bool(errors))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
