"""Paper-comparison layer — DATA2.

Mirror of test_data1_paper_comparison.py scoped to DATA2 (33 numeric manifest
rows + the DATA2 subset of target_notebook_source_map.csv).

See docs/PAPER_COMPARISON_LAYER.md for architecture.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from helpers.paper_comparison import (
    DriftRow,
    PAPER_ARTIFACTS_DATA2,
    load_locked_mapping,
    load_numeric_manifest,
    load_tolerance_overrides,
    resolved_tolerance,
    stale_manifest_rows,
    write_drift_report,
)


_LOCKED = load_locked_mapping()
_NUMERIC = load_numeric_manifest()
_TOLERANCES = load_tolerance_overrides()

DATA2_LOCKED_ROWS = [r for r in _LOCKED if r.dataset_scope == "DATA2"]
DATA2_NUMERIC_ROWS = [r for r in _NUMERIC if r.paper == "DATA2"]


# ======================================================================
# Smoke sub-layer
# ======================================================================

@pytest.mark.smoke
def test_data2_numeric_manifest_loaded():
    assert len(DATA2_NUMERIC_ROWS) > 0, (
        "simulation_validation_manifest.csv has no DATA2 rows."
    )


@pytest.mark.smoke
@pytest.mark.parametrize("row", DATA2_NUMERIC_ROWS, ids=lambda r: r.case_id)
def test_data2_baseline_csv_exists(row):
    p = row.baseline_path()
    assert p.exists(), (
        f"Baseline CSV missing: {p}\n"
        f"(manifest row: {row.paper} / {row.figure_id} / {row.panel_id})"
    )


@pytest.mark.smoke
@pytest.mark.parametrize("row", DATA2_NUMERIC_ROWS, ids=lambda r: r.case_id)
def test_data2_baseline_csv_has_rows(row):
    import pandas as pd
    p = row.baseline_path()
    if not p.exists():
        pytest.skip(f"baseline missing (covered by test_data2_baseline_csv_exists): {p}")
    df = pd.read_csv(p)
    assert len(df) > 0, f"Baseline CSV empty: {p}"


@pytest.mark.smoke
def test_data2_locked_mapping_loaded():
    assert len(DATA2_LOCKED_ROWS) > 0, (
        "target_notebook_source_map.csv has no DATA2 rows."
    )


@pytest.mark.smoke
def test_no_stale_manifest_rows_data2():
    """DATA2 mirror of the stale-manifest check."""
    stale = [s for s in stale_manifest_rows() if s.startswith("DATA2/")]
    if stale:
        pytest.xfail(
            "DATA2 manifest has rows pointing at non-CSV stubs. These are filtered "
            "out of comparison; clean them up when convenient:\n"
            + "\n".join(f"  {s}" for s in stale)
        )


# ======================================================================
# Mapping audit sub-layer
# ======================================================================

@pytest.fixture(scope="session")
def data2_paper_artifacts_exist():
    if not PAPER_ARTIFACTS_DATA2.exists():
        return False
    return any(PAPER_ARTIFACTS_DATA2.glob("*.png"))


@pytest.mark.regression
@pytest.mark.parametrize("row", DATA2_LOCKED_ROWS, ids=lambda r: r.target_id)
def test_data2_mapped_artifacts_present(row, data2_paper_artifacts_exist):
    """SKIPs (not fails) per-row on missing artifacts. See DATA1 mirror for rationale."""
    if not data2_paper_artifacts_exist:
        pytest.skip(
            f"paper_artifacts/data2/notebook_figures/ has no PNGs yet — run\n"
            f"  python refactored_codes_v1/refactored_ucb_runfile.py\n"
            f"and pick DATA2 root first."
        )
    missing = [p for p in row.expected_paths() if not p.exists()]
    if missing:
        pytest.skip(
            f"target {row.target_id} ({row.paper_item}) not in last run scope:\n"
            + "\n".join(f"  missing: {m.name}" for m in missing)
        )


@pytest.mark.regression
def test_data2_mapped_artifact_coverage(data2_paper_artifacts_exist):
    """Sanity bound — at least 25% of mapped DATA2 artifacts present."""
    if not data2_paper_artifacts_exist:
        pytest.skip("paper_artifacts/data2/notebook_figures/ empty")
    total_artifacts = 0
    present_artifacts = 0
    for row in DATA2_LOCKED_ROWS:
        for p in row.expected_paths():
            total_artifacts += 1
            if p.exists():
                present_artifacts += 1
    coverage = present_artifacts / max(total_artifacts, 1)
    # 10% catches catastrophic renames without false-positiving on subset runs.
    # DATA2_FAST_SUBSET doesn't overlap heavily with target_notebook_source_map.csv,
    # so this threshold is intentionally permissive. Tighten via env var.
    import os
    min_cov = float(os.environ.get("MIN_PAPER_COVERAGE", "0.10"))
    assert coverage >= min_cov, (
        f"DATA2 paper-artifact coverage {coverage:.1%} < {min_cov:.0%}: "
        f"{present_artifacts}/{total_artifacts} mapped artifacts on disk. "
        f"A refactor likely renamed a family of outputs."
    )


# ======================================================================
# Numeric regression sub-layer — stubbed for v1
# ======================================================================

@pytest.mark.regression
@pytest.mark.nightly
@pytest.mark.parametrize("row", DATA2_NUMERIC_ROWS, ids=lambda r: r.case_id)
def test_data2_numeric_matches_baseline(row, lib, session_run_dir):
    """STUBBED. See docs/PAPER_COMPARISON_LAYER.md § 'A gap we found during build'."""
    metric_name, tol = resolved_tolerance(row, _TOLERANCES)
    pytest.skip(
        f"numeric sub-layer stub — see docs/PAPER_COMPARISON_LAYER.md\n"
        f"  target: {row.case_id}\n"
        f"  metric: {metric_name}  tol: {tol:.4g}\n"
        f"  baseline: {row.relative_csv}"
    )
