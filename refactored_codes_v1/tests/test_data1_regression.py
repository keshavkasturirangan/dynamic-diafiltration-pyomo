"""Forward-simulation regression tests for DATA1.

Each test runs `solve_model(..., sim_opt=True)` on a known dataset with a
known theta and asserts the resulting end-of-vial values match the
captured baseline within 1%.

If a test fails:
  1. A library edit silently shifted the model formulation (loader,
     constraints, or DAE discretization) → investigate the change.
  2. Reference values are stale → if the change is intentional, re-run
     `_capture_references.py` to refresh the baseline JSON.
"""
import pytest

from reference_cases import DATA1_CASES
from conftest import DATA1_DIR
from _trajectory_compare import run_and_compare


@pytest.mark.parametrize("case", DATA1_CASES, ids=lambda c: c["case_id"])
def test_data1_trajectory_matches_reference(case, lib, references):
    """One forward-simulation regression test per DATA1 reference case."""
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
