"""Layer 4 (structural) — DATA1 regression tests.

Migrated from refactored_codes_v1/tests/test_data1_regression.py.
Imports re-pointed at the new helpers location; otherwise byte-equivalent
to the original (same skip semantics, same fixture contract).

Each test runs `solve_model(..., sim_opt=True)`-equivalent structural
fingerprinting: build the Pyomo model via `model_construct_inter` WITHOUT
solving it, capture loader + build fingerprints, compare to the committed
baseline in `baselines/regression_references.json`.

If a test fails:
  1. A library edit silently shifted the model formulation (loader,
     constraints, or DAE discretization) → investigate the change.
  2. Reference values are stale → if the change is intentional, re-run
     `python3 pytest_refactored_codes_v1/_capture_references.py` to refresh.
"""
from __future__ import annotations

import pytest

from helpers.reference_cases import DATA1_CASES
from helpers.trajectory_compare import run_and_compare
from conftest import DATA1_DIR


pytestmark = [pytest.mark.regression]


@pytest.mark.parametrize("case", DATA1_CASES, ids=lambda c: c["case_id"])
def test_data1_structural_matches_reference(case, lib, references):
    """One structural-fingerprint regression test per DATA1 reference case."""
    case_id = case["case_id"]
    if case.get("skip_reason"):
        pytest.skip(f"Test case opted out: {case['skip_reason']}")
    ref = references.get(case_id)
    if ref is None:
        pytest.fail(f"No reference values captured for {case_id!r}")
    if ref.get("skip_reason"):
        pytest.skip(f"Captured as skipped: {ref['skip_reason']}")
    if "error" in ref:
        pytest.skip(f"Reference capture for {case_id} failed: {ref['error']}")

    data_path = DATA1_DIR / case["data_file"]
    if not data_path.exists():
        pytest.skip(f"data file not on disk: {data_path}")

    run_and_compare(lib, case, ref, data_path)
