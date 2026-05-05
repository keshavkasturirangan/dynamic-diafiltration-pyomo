from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from PIL import Image
import pytest

import refactored_codes_v1.refactored_ucb_library as lib
from conftest import data_root


def _has_any_mat(root: Path) -> bool:
    return any(root.rglob("*.mat"))


def _load_png_pixels(path: Path) -> np.ndarray:
    with Image.open(path) as img:
        return np.asarray(img.convert("RGBA"))


def _first_existing_path(*paths: Path) -> Path:
    for path in paths:
        if path.exists():
            return path
    raise FileNotFoundError("None of the candidate baseline images exist.")


def test_data1_materialize_smoke_when_data_available(tmp_path):
    root = data_root("DIAFILTRATION_DATA1_ROOT")
    if root is None or not _has_any_mat(root):
        pytest.skip("Set DIAFILTRATION_DATA1_ROOT to run DATA1 materialization smoke test")
    figures = lib.list_figures("DATA1")
    if not figures:
        pytest.skip("No DATA1 figures registered")
    result = lib.materialize(figures[0], campaign="DATA1", save_dir=tmp_path, data_root=root)
    assert result["status"] in {"ok", "error"}
    if result["status"] == "ok":
        assert "paths" in result


def test_data2_materialize_smoke_when_data_available(tmp_path):
    root = data_root("DIAFILTRATION_DATA2_ROOT")
    if root is None or not _has_any_mat(root):
        pytest.skip("Set DIAFILTRATION_DATA2_ROOT to run DATA2 materialization smoke test")
    figures = lib.list_figures("DATA2")
    if not figures:
        pytest.skip("No DATA2 figures registered")
    result = lib.materialize(figures[0], campaign="DATA2", save_dir=tmp_path, data_root=root)
    assert result["status"] in {"ok", "error"}
    if result["status"] == "ok":
        assert "paths" in result


def test_data2_pressure_changes_pipeline_matches_baseline(tmp_path):
    root = data_root("DIAFILTRATION_DATA2_ROOT")
    if root is None or not _has_any_mat(root):
        pytest.skip("Set DIAFILTRATION_DATA2_ROOT to run DATA2 equivalence test")
    result = lib.materialize("pressure_changes", campaign="DATA2", save_dir=tmp_path, data_root=root)
    assert result["status"] == "ok"
    assert len(result["paths"]) == 2

    new_lag = Path(result["paths"][0])
    new_overflow = Path(result["paths"][1])
    baseline_lag = _first_existing_path(
        Path("UnifiedFramework/DATA3/results/paper_artifacts/data2/notebook_figures/pressure_change_lag.png"),
        Path("artifacts/figures/root_legacy/pressure_change_lag.png"),
        Path("pressure_change_lag.png"),
    )
    baseline_overflow = _first_existing_path(
        Path("UnifiedFramework/DATA3/results/paper_artifacts/data2/notebook_figures/pressure_change_overflow.png"),
        Path("artifacts/figures/root_legacy/pressure_change_overflow.png"),
        Path("pressure_change_overflow.png"),
    )

    np.testing.assert_array_equal(_load_png_pixels(new_lag), _load_png_pixels(baseline_lag))
    np.testing.assert_array_equal(_load_png_pixels(new_overflow), _load_png_pixels(baseline_overflow))


def test_data2_calibration_plots_pipeline_matches_baseline(tmp_path):
    root = data_root("DIAFILTRATION_DATA2_ROOT")
    if root is None or not _has_any_mat(root):
        pytest.skip("Set DIAFILTRATION_DATA2_ROOT to run DATA2 equivalence test")
    result = lib.materialize("calibration_plots", campaign="DATA2", save_dir=tmp_path, data_root=root)
    assert result["status"] == "ok"
    assert len(result["paths"]) >= 4

    new_a = Path(result["paths"][0])
    new_b = Path(result["paths"][1])
    new_s1 = Path(result["paths"][2])
    new_pub = Path(result["paths"][3])
    baseline_a = _first_existing_path(
        Path("UnifiedFramework/DATA3/results/paper_artifacts/data2/notebook_figures/calib_curve_a.png"),
        Path("artifacts/figures/root_legacy/calib_curve_a.png"),
    )
    baseline_b = _first_existing_path(
        Path("UnifiedFramework/DATA3/results/paper_artifacts/data2/notebook_figures/calib_curve_b.png"),
        Path("artifacts/figures/root_legacy/calib_curve_b.png"),
    )
    baseline_s1 = _first_existing_path(
        Path("UnifiedFramework/DATA3/results/paper_artifacts/data2/notebook_figures/figure_s1.png"),
        Path("artifacts/figures/root_legacy/figure_s1.png"),
    )
    baseline_pub = _first_existing_path(
        Path("UnifiedFramework/DATA3/results/paper_artifacts/data2/notebook_figures/calib_curve.png"),
        Path("artifacts/figures/root_legacy/calib_curve.png"),
    )

    np.testing.assert_array_equal(_load_png_pixels(new_a), _load_png_pixels(baseline_a))
    np.testing.assert_array_equal(_load_png_pixels(new_b), _load_png_pixels(baseline_b))
    np.testing.assert_array_equal(_load_png_pixels(new_s1), _load_png_pixels(baseline_s1))
    np.testing.assert_array_equal(_load_png_pixels(new_pub), _load_png_pixels(baseline_pub))


def test_data2_startup_barplot_pipeline_matches_baseline(tmp_path):
    result = lib.materialize(
        "startup_barplot",
        campaign="DATA2",
        save_dir=tmp_path,
        extra_opts={"results": {"meta": {"pipeline": "test"}}},
    )
    assert result["status"] == "ok"
    assert len(result["paths"]) == 1

    new_path = Path(result["paths"][0])
    baseline = _first_existing_path(
        Path("UnifiedFramework/DATA3/results/paper_artifacts/data2/notebook_figures/startup_barplot.png"),
        Path("artifacts/figures/root_legacy/startup_barplot.png"),
        Path("startup_barplot.png"),
    )

    np.testing.assert_array_equal(_load_png_pixels(new_path), _load_png_pixels(baseline))


def test_data2_model_error_visualization_pipeline_matches_baseline(tmp_path):
    root = data_root("DIAFILTRATION_DATA2_ROOT")
    if root is None or not _has_any_mat(root):
        pytest.skip("Set DIAFILTRATION_DATA2_ROOT to run DATA2 equivalence test")
    result = lib.materialize(
        "model_error_visualization",
        campaign="DATA2",
        save_dir=tmp_path,
        data_root=root,
    )
    assert result["status"] == "ok"
    assert len(result["paths"]) == 3

    new_startup = Path(result["paths"][0])
    new_conc = Path(result["paths"][1])
    new_dil = Path(result["paths"][2])
    baseline_startup = _first_existing_path(
        Path("UnifiedFramework/DATA3/results/paper_artifacts/data2/notebook_figures/startup_barplot.png"),
        Path("artifacts/figures/root_legacy/startup_barplot.png"),
    )
    baseline_conc = _first_existing_path(
        Path("UnifiedFramework/DATA3/results/paper_artifacts/data2/notebook_figures/concentrating_residuals_boxplot.png"),
        Path("artifacts/figures/root_legacy/concentrating_residuals_boxplot.png"),
    )
    baseline_dil = _first_existing_path(
        Path("UnifiedFramework/DATA3/results/paper_artifacts/data2/notebook_figures/diluting_residuals_boxplot.png"),
        Path("artifacts/figures/root_legacy/diluting_residuals_boxplot.png"),
    )

    np.testing.assert_array_equal(_load_png_pixels(new_startup), _load_png_pixels(baseline_startup))
    np.testing.assert_array_equal(_load_png_pixels(new_conc), _load_png_pixels(baseline_conc))
    np.testing.assert_array_equal(_load_png_pixels(new_dil), _load_png_pixels(baseline_dil))


@pytest.mark.skipif(
    os.environ.get("DIAFILTRATION_RUN_FULL_EQUIVALENCE", "0") != "1",
    reason="Full output equivalence is opt-in because it runs legacy and refactored workflows.",
)
def test_full_equivalence_placeholder(tmp_path):
    """Fill this test one figure at a time during migration.

    Recommended pattern:
      1. Run legacy function into tmp_path / 'legacy'.
      2. Run lib.materialize(...) into tmp_path / 'new'.
      3. Compare numeric arrays with np.testing.assert_allclose.
      4. Compare images with exact hash or perceptual hash.
    """
    pytest.skip("Implement per-figure equivalence checks as figures are migrated")
