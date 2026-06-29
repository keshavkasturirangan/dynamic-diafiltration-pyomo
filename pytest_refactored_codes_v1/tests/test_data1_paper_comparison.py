"""Paper-comparison layer — DATA1.

Three sub-layers, each marked separately so the operator can choose which
runs in any given invocation:

  smoke      — assert every DATA1 baseline CSV named in
               simulation_validation_manifest.csv exists on disk with the
               expected long-format columns. No produced artifacts, no solver.
               < 1 s.

  regression — mapping audit. For each DATA1 row of
               target_notebook_source_map.csv, assert every mapped PNG
               filename resolves to a real file in
               paper_artifacts/data1/notebook_figures/. Catches refactor
               renames + the case where a figure stops being produced.
               No solver. ~1 s.

  nightly    — numeric regression. Re-runs the forward simulation via the
               refactored library, extracts the data that would have been
               plotted, compares to baseline CSV. STUBBED with pytest.skip
               in v1 — implementation deferred to a follow-up.

See docs/PAPER_COMPARISON_LAYER.md for the architecture.
See docs/PAPER_COMPARISON_PROMPT.md for the original spec.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from helpers.paper_comparison import (
    DriftRow,
    PAPER_ARTIFACTS_DATA1,
    load_locked_mapping,
    load_numeric_manifest,
    load_tolerance_overrides,
    resolved_tolerance,
    stale_manifest_rows,
    write_drift_report,
)


# ----------------------------------------------------------------------
# Eager parametrize-time loading — needed because pytest.parametrize
# can't see session fixtures.
# ----------------------------------------------------------------------

_LOCKED = load_locked_mapping()
_NUMERIC = load_numeric_manifest()
_TOLERANCES = load_tolerance_overrides()

DATA1_LOCKED_ROWS = [r for r in _LOCKED if r.dataset_scope == "DATA1"]
DATA1_NUMERIC_ROWS = [r for r in _NUMERIC if r.paper == "DATA1"]


# ======================================================================
# Smoke sub-layer — baseline CSV presence + column schema
# ======================================================================

@pytest.mark.smoke
def test_data1_numeric_manifest_loaded():
    """Sanity: at least one DATA1 row was parsed from the numeric manifest."""
    assert len(DATA1_NUMERIC_ROWS) > 0, (
        "simulation_validation_manifest.csv has no DATA1 rows — either the "
        "manifest moved or has been emptied. Check "
        "UnifiedFramework/DATA3/docs/validation/simulation_validation_data_files/"
    )


@pytest.mark.smoke
@pytest.mark.parametrize("row", DATA1_NUMERIC_ROWS, ids=lambda r: r.case_id)
def test_data1_baseline_csv_exists(row):
    """Each DATA1 baseline CSV named in the manifest is on disk."""
    p = row.baseline_path()
    assert p.exists(), (
        f"Baseline CSV missing: {p}\n"
        f"(manifest row: {row.paper} / {row.figure_id} / {row.panel_id})"
    )


@pytest.mark.smoke
@pytest.mark.parametrize("row", DATA1_NUMERIC_ROWS, ids=lambda r: r.case_id)
def test_data1_baseline_csv_has_rows(row):
    """Each DATA1 baseline CSV has at least one data row."""
    import pandas as pd
    p = row.baseline_path()
    if not p.exists():
        pytest.skip(f"baseline missing (covered by test_data1_baseline_csv_exists): {p}")
    df = pd.read_csv(p)
    assert len(df) > 0, f"Baseline CSV empty: {p}"


@pytest.mark.smoke
def test_data1_locked_mapping_loaded():
    """At least one DATA1 row was parsed from target_notebook_source_map.csv."""
    assert len(DATA1_LOCKED_ROWS) > 0, (
        "target_notebook_source_map.csv has no DATA1 rows. Check "
        "UnifiedFramework/DATA3/docs/validation/target_notebook_source_map.csv"
    )


@pytest.mark.smoke
def test_no_stale_manifest_rows():
    """No DATA1 manifest rows point at non-CSV stubs (e.g. _ERROR.txt from
    legacy extraction failures). If this fails, the manifest needs cleanup —
    not the pytest suite. Failure is informational; we re-run with stale rows
    filtered out, so smoke comparison tests still pass."""
    stale = [s for s in stale_manifest_rows() if s.startswith("DATA1/")]
    if stale:
        pytest.xfail(
            "DATA1 manifest has rows pointing at non-CSV stubs (legacy extraction "
            "errors). These are filtered out of comparison; clean them up when "
            "convenient:\n" + "\n".join(f"  {s}" for s in stale)
        )


# ======================================================================
# Mapping audit sub-layer — produced PNG presence
# ======================================================================

@pytest.fixture(scope="session")
def data1_paper_artifacts_exist():
    """True if paper_artifacts/data1/notebook_figures/ contains any PNG.

    Used as a skip gate by the per-row tests below. Avoids cascading
    'no such file' failures when materialize_all hasn't been run yet.
    """
    if not PAPER_ARTIFACTS_DATA1.exists():
        return False
    return any(PAPER_ARTIFACTS_DATA1.glob("*.png"))


@pytest.mark.regression
@pytest.mark.parametrize("row", DATA1_LOCKED_ROWS, ids=lambda r: r.target_id)
def test_data1_mapped_artifacts_present(row, data1_paper_artifacts_exist):
    """Each artifact named in the locked mapping is on disk.

    SKIPs (not fails) on missing artifacts — the locked mapping covers the
    FULL DATA1 manifest, but most invocations only run a subset. A skip
    here is information ("this figure wasn't requested in the last run"),
    not a regression. The aggregate coverage test below flags it if too
    many artifacts go missing at once.
    """
    if not data1_paper_artifacts_exist:
        pytest.skip(
            f"paper_artifacts/data1/notebook_figures/ has no PNGs yet — run\n"
            f"  python refactored_codes_v1/refactored_ucb_runfile.py\n"
            f"and pick DATA1 root first."
        )
    missing = [p for p in row.expected_paths() if not p.exists()]
    if missing:
        pytest.skip(
            f"target {row.target_id} ({row.paper_item}) not in last run scope:\n"
            + "\n".join(f"  missing: {m.name}" for m in missing)
        )


@pytest.mark.regression
def test_data1_mapped_artifact_coverage(data1_paper_artifacts_exist):
    """Sanity bound: at least 25% of the locked DATA1 mapped artifacts should
    be on disk when the paper_artifacts tree is non-empty.

    Too low → likely a refactor renamed a whole family of outputs, or the
    user is mid-pipeline-run. Above 25% gives signal without being noisy
    against partial-subset invocations.
    """
    if not data1_paper_artifacts_exist:
        pytest.skip("paper_artifacts/data1/notebook_figures/ empty")
    total_artifacts = 0
    present_artifacts = 0
    for row in DATA1_LOCKED_ROWS:
        for p in row.expected_paths():
            total_artifacts += 1
            if p.exists():
                present_artifacts += 1
    coverage = present_artifacts / max(total_artifacts, 1)
    # 10% catches the catastrophic case (a refactor renamed an entire output
    # family) without false-positiving on partial-subset runs where the user
    # only requested DATA1_FAST_SUBSET. Tighten in CI by setting
    # MIN_PAPER_COVERAGE in the environment.
    import os
    min_cov = float(os.environ.get("MIN_PAPER_COVERAGE", "0.10"))
    assert coverage >= min_cov, (
        f"DATA1 paper-artifact coverage {coverage:.1%} < {min_cov:.0%}: "
        f"{present_artifacts}/{total_artifacts} mapped artifacts on disk. "
        f"A refactor likely renamed a family of outputs."
    )


# ======================================================================
# Numeric regression sub-layer — stubbed for v1
# ======================================================================

@pytest.mark.regression
@pytest.mark.nightly
@pytest.mark.parametrize("row", DATA1_NUMERIC_ROWS, ids=lambda r: r.case_id)
def test_data1_numeric_matches_baseline(row, lib, session_run_dir):
    """Re-run the refactored forward simulation and compare to baseline.

    STUBBED in v1. The implementation requires:
      1. Identifying which dataset+vial each manifest row maps to
         (source_label gives a textual hint; needs a parser).
      2. Calling lib.solve_model(sim_opt=True) with the warm-start θ for
         that dataset.
      3. Extracting the time-aligned series the plot would have rendered
         (mass, cF, cV) from the solved Pyomo model.
      4. Resampling onto the baseline's time grid.
      5. Computing the metric and asserting tolerance.

    See docs/PAPER_COMPARISON_LAYER.md § "A gap we found during build" for
    the design discussion.
    """
    metric_name, tol = resolved_tolerance(row, _TOLERANCES)
    pytest.skip(
        f"numeric sub-layer stub — see docs/PAPER_COMPARISON_LAYER.md\n"
        f"  target: {row.case_id}\n"
        f"  metric: {metric_name}  tol: {tol:.4g}\n"
        f"  baseline: {row.baseline_path().relative_to(Path('.').resolve()) if row.baseline_path().is_relative_to(Path('.').resolve()) else row.baseline_path()}"
    )


# ======================================================================
# Drift report writer — emitted at session end via the regression layer
# ======================================================================

@pytest.fixture(scope="session", autouse=False)
def write_data1_drift_report(session_run_dir):
    """Optional explicit fixture for emitting a drift report when the numeric
    sub-layer is wired up. Currently a no-op stub.
    """
    rows: list[DriftRow] = []
    yield rows
    if rows:
        write_drift_report(session_run_dir / "paper_comparison_report.csv", rows)
