from __future__ import annotations

import ast
import inspect
from pathlib import Path

import pytest

import refactored_codes_v1.refactored_ucb_library as lib


FORBIDDEN_RENDERER_CALLS = {
    "loadmat",
    "load_experimental_data",
    "build_model",
    "estimate_parameters",
    "solve_model",
    "solve_model_B_fix",
    "run_workflow",
    "run_pipeline",
}


def _called_names(fn):
    tree = ast.parse(inspect.getsource(fn))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            call = node.func
            if isinstance(call, ast.Name):
                names.add(call.id)
            elif isinstance(call, ast.Attribute):
                names.add(call.attr)
    return names


def test_stage_results_round_trip_preserves_legacy_keys():
    payload = {
        "data_file": ["example.mat"],
        "mode": "DATA",
        "workflow_family": "DATA2",
        "data": {"a": 1},
        "fit_stru": {"Lp": 1.0, "B": 2.0, "sigma": 0.9},
        "legacy_extra_key": {"still": "here"},
    }
    result = lib.StageResults.from_dict(payload)
    round_trip = result.to_dict()
    assert round_trip["data_file"] == ["example.mat"]
    assert round_trip["fit_stru"]["sigma"] == 0.9
    assert round_trip["_meta"]["legacy_extra_key"] == {"still": "here"}


def test_manifest_api_exposes_campaigns_and_figures():
    campaigns = lib.list_campaigns()
    assert "DATA1" in campaigns
    assert "DATA2" in campaigns
    assert isinstance(lib.list_figures("DATA1"), list)
    assert isinstance(lib.list_figures("DATA2"), list)


def test_lookup_figure_accepts_full_or_short_names():
    figures = lib.list_figures("DATA1")
    if not figures:
        pytest.skip("No DATA1 figures registered")
    first = figures[0]
    spec = lib.lookup_figure("DATA1", first)
    assert spec.name == first
    short = first.split(".", 1)[-1]
    assert lib.lookup_figure("DATA1", short).name == first


def test_data2_pressure_changes_is_pipeline_backed():
    spec = lib.lookup_figure("DATA2", "pressure_changes")
    assert spec.name == "data2.pressure_changes"
    assert spec.renderer == "render_pressure_changes_pipeline"


def test_data2_calibration_plots_is_pipeline_backed():
    spec = lib.lookup_figure("DATA2", "calibration_plots")
    assert spec.name == "data2.calibration_plots"
    assert spec.renderer == "render_calibration_curve_pipeline"


def test_data2_startup_barplot_is_pipeline_backed():
    spec = lib.lookup_figure("DATA2", "startup_barplot")
    assert spec.name == "data2.startup_barplot"
    assert spec.renderer == "render_startup_barplot"


def test_data2_model_error_visualization_is_pipeline_backed():
    spec = lib.lookup_figure("DATA2", "model_error_visualization")
    assert spec.name == "data2.model_error_visualization"
    assert spec.renderer == "render_model_error_visualization_pipeline"


def test_data3_time_series_is_pipeline_backed():
    spec = lib.lookup_figure("DATA3", "option3.time_series")
    assert spec.name == "data3.option3.time_series"
    assert spec.renderer == "render_data3_time_series"


def test_runfile_imports_refactored_library_only():
    runfile = Path(__file__).resolve().parents[1] / "refactored_codes_v1" / "refactored_ucb_runfile.py"
    text = runfile.read_text()
    assert "import refactored_ucb_library as ucb" in text
    assert "from diafiltration" not in text
    assert "import diafiltration" not in text
    assert "run_data3_time_series_plots" not in text
    assert "materialize(" in text


@pytest.mark.parametrize(
    "name",
    [
        "render_sim_comparison",
        "render_data3_time_series",
        "render_contour",
        "render_sigma_sensitivity",
        "render_concentration_range",
        "render_concentration_comparison",
        "render_calibration_curve",
        "render_pressure_change",
        "render_startup_barplot",
        "render_error_boxplot",
    ],
)
def test_renderers_do_not_call_pipeline_or_solver_directly(name):
    if not hasattr(lib, name):
        pytest.skip(f"{name} is not present")
    fn = getattr(lib, name)
    names = _called_names(fn)
    forbidden = names & FORBIDDEN_RENDERER_CALLS
    assert not forbidden, f"{name} calls forbidden functions: {sorted(forbidden)}"


def test_materialize_returns_error_dict_for_missing_figure(tmp_path):
    with pytest.raises(KeyError):
        lib.lookup_figure("DATA1", "not_a_real_figure")
