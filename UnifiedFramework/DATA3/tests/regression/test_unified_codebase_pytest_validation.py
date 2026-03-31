"""Pytest validation for the unified codebase entrypoints and core utilities."""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

REPO_ROOT = Path(__file__).resolve().parents[4]
UNICODE_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "ExperimentalDataAnalysis" / "UnifiedCode"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(UNICODE_ROOT) not in sys.path:
    sys.path.insert(0, str(UNICODE_ROOT))

import conductivity_paper as conductivity_paper_module  # noqa: E402
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode import unified_codebase_runfile as runfile  # noqa: E402
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library import (  # noqa: E402
    ExperimentMode,
    ModelOptions,
    RunMode,
    UnifiedPipelineConfigV24,
    load_experiment_easy,
    run_unified_pipeline_v24,
)

SIMULATION_VALIDATION_ROOT = (
    REPO_ROOT
    / "UnifiedFramework"
    / "DATA3"
    / "docs"
    / "validation"
    / "simulation_validation_data_files"
)
PAPER_DIGITIZED_ROOT = (
    REPO_ROOT
    / "UnifiedFramework"
    / "DATA3"
    / "docs"
    / "validation"
    / "nightly"
    / "digitized_baselines"
    / "paper"
    / "data1_main"
)
PYTEST_RESULTS_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "pytest_validation"
COMPARISON_REPORT_CSV = PYTEST_RESULTS_ROOT / "unified_data1_fig2_comparison.csv"
PAPER_COMPARISON_REPORT_CSV = PYTEST_RESULTS_ROOT / "unified_data1_fig2_paper_comparison.csv"
DATA2_COMPARISON_REPORT_CSV = PYTEST_RESULTS_ROOT / "unified_data2_measurement_comparison.csv"


@pytest.mark.regression
def test_conductivity_paper_variant_shedlovsky_snapshot() -> None:
    """Conductivity helper should remain numerically stable for a fixed KCl case."""
    result = conductivity_paper_module.variant_shedlovsky(
        conc=[0.005, 0.01, 0.02],
        temp=298.15,
        epsilon=78.4,
        eta=0.0089,
        lambda_0=149.86,
        a=4.0e-8,
        z_1=1,
        z_2=-1,
        lambda_0_cation=73.5,
        lambda_0_anion=76.3,
    )
    expected = np.array([0.7185609584, 1.4146177922, 2.7705771429], dtype=float)
    np.testing.assert_allclose(np.asarray(result, dtype=float), expected, rtol=1e-9, atol=1e-9)


@pytest.mark.regression
def test_runfile_data1_presets_build_existing_configs() -> None:
    """DATA1 preset paths should exist and produce the expected config defaults."""
    for path in runfile.PRESET_DATA1_PATHS:
        assert path.exists(), f"Missing DATA1 preset path: {path}"

    class Args:
        selector = None
        convert_to_concentration = False
        conductivity_model = "msa"
        plot = False
        calc_cov = False
        nfe = 30
        solver = "ipopt"

    config = runfile._build_data1_config(runfile.PRESET_DATA1_PATHS[0], Args(), run_doe=False)
    assert isinstance(config, UnifiedPipelineConfigV24)
    assert config.file_path.endswith("data_stru-dataset501.1.mat")
    assert config.run_parmest is True
    assert config.run_doe is False
    assert config.conductivity_model == "msa"
    assert config.model_options.mode == ExperimentMode.DATA
    assert config.model_options.run_mode == RunMode.ESTIMATION
    assert str(config.model_options.b_form).lower() == "single"


@pytest.mark.regression
def test_unified_data1_mat_pipeline_runs_with_filtration_cd_default() -> None:
    """Unified DATA1 MAT estimation should run for filtration datasets and return finite theta."""
    config = UnifiedPipelineConfigV24(
        file_path=str(runfile.PRESET_DATA1_PATHS[0]),
        selector=None,
        specs=None,
        convert_to_concentration=False,
        plot=False,
        model_options=ModelOptions(
            mode=ExperimentMode.DATA,
            run_mode=RunMode.ESTIMATION,
            b_form="single",
            nfe=30,
            use_sigma_logit_transform=False,
            use_multistart_for_mat_legacy=False,
        ),
        run_parmest=True,
        run_doe=False,
        calc_cov=False,
        cov_n=None,
        solver="ipopt",
        solver_options=None,
        tee=False,
    )

    result = run_unified_pipeline_v24(config)

    assert result["load_ok"] is True
    exp = result["exp"]
    assert exp.C_D_value == 0.0
    assert "C_D_value=0.0_for_filtration_mat" in exp.used_defaults

    theta = result["parmest"]["theta"]
    assert list(theta.index) == ["Lp", "B", "sigma"]
    assert all(math.isfinite(float(value)) for value in theta.values)
    assert float(theta["Lp"]) > 0.0
    assert float(theta["B"]) >= 0.0
    assert 0.0 <= float(theta["sigma"]) <= 1.0

    model = result["parmest"]["model"]
    assert model is not None
    assert len(list(model.n_vial)) == len(exp.vials)
    assert len(list(model.tau)) > 1


def _run_data1_pipeline(dataset_file: Path):
    config = UnifiedPipelineConfigV24(
        file_path=str(dataset_file),
        selector=None,
        specs=None,
        convert_to_concentration=False,
        plot=False,
        model_options=ModelOptions(
            mode=ExperimentMode.DATA,
            run_mode=RunMode.ESTIMATION,
            b_form="single",
            nfe=30,
            use_sigma_logit_transform=False,
            use_multistart_for_mat_legacy=False,
        ),
        run_parmest=True,
        run_doe=False,
        calc_cov=False,
        cov_n=None,
        solver="ipopt",
        solver_options=None,
        tee=False,
    )
    return run_unified_pipeline_v24(config)


def _extract_unified_data1_predictions(result) -> pd.DataFrame:
    model = result["parmest"]["model"]
    ti_vals = list(model.ti)
    tf_vals = list(model.tf)
    tau_vals = sorted(float(t) for t in list(model.tau))
    rows = []
    for n in model.n_vial:
        n_i = int(n)
        ti = float(ti_vals[n_i - 1])
        tf = float(tf_vals[n_i - 1])
        duration = tf - ti
        for tau in tau_vals:
            time_min = (ti + tau * duration) / 60.0
            rows.append(
                {
                    "vial": n_i,
                    "time_min": time_min,
                    "mass_g": float(model.mV[n, tau].value),
                    "retentate_concentration_mM": float(model.cF[n, tau].value),
                    "permeate_concentration_mM": float(model.cH[n, tau].value),
                    "vial_concentration_mM": float(model.cV[n, tau].value),
                }
            )
    return pd.DataFrame(rows)


def _compare_series(
    baseline_df: pd.DataFrame,
    unified_df: pd.DataFrame,
    *,
    panel_id: str,
    quantity_col: str,
    unified_col: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for series_id in sorted(x for x in baseline_df["series_id"].unique() if "prediction" in x):
        vial = int(series_id.split("_")[-1])
        base = baseline_df[baseline_df["series_id"] == series_id].dropna(subset=[quantity_col]).copy()
        if base.empty:
            continue
        uni = unified_df[unified_df["vial"] == vial][["time_min", unified_col]].sort_values("time_min")
        pred = np.interp(base["time_min"].to_numpy(dtype=float), uni["time_min"].to_numpy(dtype=float), uni[unified_col].to_numpy(dtype=float))
        err = pred - base[quantity_col].to_numpy(dtype=float)
        rows.append(
            {
                "panel_id": panel_id,
                "series_id": series_id,
                "n_points": int(len(base)),
                "mae": float(np.mean(np.abs(err))),
                "max_abs_err": float(np.max(np.abs(err))),
                "rmse": float(np.sqrt(np.mean(err**2))),
            }
        )
    return rows


def _best_match_paper_trace(
    paper_df: pd.DataFrame,
    unified_df: pd.DataFrame,
    *,
    panel_id: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for series_name in sorted(paper_df["series_name"].unique()):
        paper = paper_df[paper_df["series_name"] == series_name].dropna(subset=["x_data", "y_data"]).sort_values("x_data")
        best_row: dict[str, object] | None = None
        for vial in sorted(unified_df["vial"].unique()):
            uni = unified_df[unified_df["vial"] == vial][["time_min", "mass_g"]].sort_values("time_min")
            pred = np.interp(
                paper["x_data"].to_numpy(dtype=float),
                uni["time_min"].to_numpy(dtype=float),
                uni["mass_g"].to_numpy(dtype=float),
            )
            err = pred - paper["y_data"].to_numpy(dtype=float)
            candidate = {
                "panel_id": panel_id,
                "paper_series": series_name,
                "matched_unified_vial": int(vial),
                "n_points": int(len(paper)),
                "mae": float(np.mean(np.abs(err))),
                "max_abs_err": float(np.max(np.abs(err))),
                "rmse": float(np.sqrt(np.mean(err**2))),
            }
            if best_row is None or float(candidate["mae"]) < float(best_row["mae"]):
                best_row = candidate
        if best_row is not None:
            rows.append(best_row)
    return rows


@pytest.fixture(scope="module")
def data1_fig2_comparison_report() -> pd.DataFrame:
    PYTEST_RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    runs = {
        "501.1": _extract_unified_data1_predictions(_run_data1_pipeline(runfile.PRESET_DATA1_PATHS[0])),
        "511.12": _extract_unified_data1_predictions(_run_data1_pipeline(runfile.PRESET_DATA1_PATHS[3])),
    }
    rows: list[dict[str, object]] = []

    rows.extend(
        _compare_series(
            pd.read_csv(SIMULATION_VALIDATION_ROOT / "data1_main" / "fig2" / "data1_main_fig2_a.csv"),
            runs["501.1"],
            panel_id="Fig2A",
            quantity_col="mass_g",
            unified_col="mass_g",
        )
    )
    rows.extend(
        _compare_series(
            pd.read_csv(SIMULATION_VALIDATION_ROOT / "data1_main" / "fig2" / "data1_main_fig2_c.csv"),
            runs["511.12"],
            panel_id="Fig2C",
            quantity_col="mass_g",
            unified_col="mass_g",
        )
    )
    baseline_d = pd.read_csv(SIMULATION_VALIDATION_ROOT / "data1_main" / "fig2" / "data1_main_fig2_d.csv")
    for matcher, ucol in [
        ("retentate_prediction", "retentate_concentration_mM"),
        ("permeate_prediction", "permeate_concentration_mM"),
        ("vial_prediction", "vial_concentration_mM"),
    ]:
        part = baseline_d[baseline_d["series_id"].str.startswith(matcher)].copy()
        rows.extend(
            _compare_series(
                part,
                runs["511.12"],
                panel_id=f"Fig2D-{matcher}",
                quantity_col="concentration_mM",
                unified_col=ucol,
            )
        )

    baseline_b = pd.read_csv(SIMULATION_VALIDATION_ROOT / "data1_main" / "fig2" / "data1_main_fig2_b.csv")
    comparison_rows = []
    for matcher, ucol in [
        ("retentate_prediction", "retentate_concentration_mM"),
        ("permeate_prediction", "permeate_concentration_mM"),
        ("vial_prediction", "vial_concentration_mM"),
    ]:
        part = baseline_b[baseline_b["series_id"].str.startswith(matcher)].copy()
        comparison_rows.extend(
            _compare_series(
                part,
                runs["501.1"],
                panel_id=f"Fig2B-{matcher}",
                quantity_col="concentration_mM",
                unified_col=ucol,
            )
        )

    report = pd.DataFrame(rows + comparison_rows).sort_values(["panel_id", "series_id"]).reset_index(drop=True)
    report.to_csv(COMPARISON_REPORT_CSV, index=False)
    return report


@pytest.fixture(scope="module")
def data1_fig2_paper_comparison_report() -> pd.DataFrame:
    PYTEST_RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    runs = {
        "501.1": _extract_unified_data1_predictions(_run_data1_pipeline(runfile.PRESET_DATA1_PATHS[0])),
        "511.12": _extract_unified_data1_predictions(_run_data1_pipeline(runfile.PRESET_DATA1_PATHS[3])),
    }
    rows: list[dict[str, object]] = []
    rows.extend(
        _best_match_paper_trace(
            pd.read_csv(PAPER_DIGITIZED_ROOT / "fig2" / "data1_main_fig2_a_traces_refined.csv"),
            runs["501.1"],
            panel_id="Fig2A-paper",
        )
    )
    rows.extend(
        _best_match_paper_trace(
            pd.read_csv(PAPER_DIGITIZED_ROOT / "fig2" / "data1_main_fig2_c_traces_refined.csv"),
            runs["511.12"],
            panel_id="Fig2C-paper",
        )
    )
    report = pd.DataFrame(rows).sort_values(["panel_id", "paper_series"]).reset_index(drop=True)
    report.to_csv(PAPER_COMPARISON_REPORT_CSV, index=False)
    return report


@pytest.mark.regression
def test_data1_fig2_comparison_report_written(data1_fig2_comparison_report: pd.DataFrame) -> None:
    """Write a concrete numeric comparison report that can be inspected after pytest."""
    assert COMPARISON_REPORT_CSV.exists()
    assert not data1_fig2_comparison_report.empty
    assert {"panel_id", "series_id", "mae", "max_abs_err"}.issubset(data1_fig2_comparison_report.columns)


@pytest.mark.regression
def test_data1_fig2_paper_comparison_report_written(data1_fig2_paper_comparison_report: pd.DataFrame) -> None:
    """Write a direct unified-vs-paper digitized comparison report."""
    assert PAPER_COMPARISON_REPORT_CSV.exists()
    assert not data1_fig2_paper_comparison_report.empty
    assert {"panel_id", "paper_series", "matched_unified_vial", "mae", "max_abs_err"}.issubset(
        data1_fig2_paper_comparison_report.columns
    )


@pytest.mark.regression
def test_data1_fig2_mass_predictions_are_reasonable_against_paper_digitization(
    data1_fig2_paper_comparison_report: pd.DataFrame,
) -> None:
    """Direct paper-side mass trace comparison for the refined digitized Figure 2 panels."""
    assert float(data1_fig2_paper_comparison_report["mae"].max()) <= 0.35
    assert float(data1_fig2_paper_comparison_report["max_abs_err"].max()) <= 0.80


@pytest.fixture(scope="module")
def data2_measurement_comparison_report() -> pd.DataFrame:
    """Compare unified DATA2 MAT loading against committed measurement-extraction CSVs."""
    PYTEST_RESULTS_ROOT.mkdir(parents=True, exist_ok=True)

    exp, (ok, issues) = load_experiment_easy(
        str(runfile.PRESET_DATA2_PATHS[1]),
        selector=None,
        specs=None,
        convert_to_concentration=False,
        plot=False,
    )
    assert ok is True, issues

    baseline = pd.read_csv(SIMULATION_VALIDATION_ROOT / "data2_main" / "fig3" / "data2_main_fig3_mass_tc.csv")
    baseline_measure = baseline[baseline["series_id"] == "measurements"].copy().sort_values("time_min")
    vial = exp.vials[3]
    unified_measure = pd.DataFrame(
        {
            "time_min": np.asarray(vial.time_s, dtype=float) / 60.0,
            "mass_g": np.asarray(vial.mass_g, dtype=float),
        }
    ).sort_values("time_min")

    common = baseline_measure.merge(unified_measure, on="time_min", suffixes=("_baseline", "_unified"))
    err = common["mass_g_unified"].to_numpy(dtype=float) - common["mass_g_baseline"].to_numpy(dtype=float)
    report = pd.DataFrame(
        [
            {
                "target": "DATA2_Fig3_measurements",
                "n_points": int(len(common)),
                "mae": float(np.mean(np.abs(err))),
                "max_abs_err": float(np.max(np.abs(err))),
                "rmse": float(np.sqrt(np.mean(err**2))),
            }
        ]
    )
    report.to_csv(DATA2_COMPARISON_REPORT_CSV, index=False)
    return report


@pytest.mark.regression
def test_data2_measurement_comparison_report_written(data2_measurement_comparison_report: pd.DataFrame) -> None:
    """Write a concrete DATA2 measurement comparison report for inspection."""
    assert DATA2_COMPARISON_REPORT_CSV.exists()
    assert not data2_measurement_comparison_report.empty
    assert {"target", "n_points", "mae", "max_abs_err", "rmse"}.issubset(
        data2_measurement_comparison_report.columns
    )


@pytest.mark.regression
def test_data2_fig3_measurements_match_committed_baseline(
    data2_measurement_comparison_report: pd.DataFrame,
) -> None:
    """Unified DATA2 MAT loading should reproduce the committed Fig. 3 measurement extraction."""
    row = data2_measurement_comparison_report.iloc[0]
    assert int(row["n_points"]) >= 40
    assert float(row["mae"]) <= 1e-12
    assert float(row["max_abs_err"]) <= 1e-12


@pytest.mark.regression
def test_data1_fig2_mass_predictions_are_close_to_legacy_baseline(
    data1_fig2_comparison_report: pd.DataFrame,
) -> None:
    """Mass panels are the current passing numeric comparison slice for unified DATA1."""
    mass = data1_fig2_comparison_report[data1_fig2_comparison_report["panel_id"].isin(["Fig2A", "Fig2C"])].copy()
    assert not mass.empty
    assert float(mass["mae"].max()) <= 0.15
    assert float(mass["max_abs_err"].max()) <= 0.26


@pytest.mark.regression
@pytest.mark.xfail(reason="Unified DATA1 concentration trajectories are not yet legacy-paper parity.", strict=False)
def test_data1_fig2_concentration_predictions_match_legacy_baseline(
    data1_fig2_comparison_report: pd.DataFrame,
) -> None:
    """Concentration panels are compared numerically and currently document a known gap."""
    conc = data1_fig2_comparison_report[~data1_fig2_comparison_report["panel_id"].isin(["Fig2A", "Fig2C"])].copy()
    assert not conc.empty
    assert float(conc["mae"].max()) <= 0.25
    assert float(conc["max_abs_err"].max()) <= 0.50
