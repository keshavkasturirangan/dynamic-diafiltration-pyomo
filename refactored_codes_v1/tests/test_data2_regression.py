"""Forward-simulation regression tests for DATA2.

Mirrors test_data1_regression.py — same architecture (sim_opt=True,
end-of-vial fingerprint comparison). DATA2 cases use B_form=1 (numeric)
and an explicit theta dict, matching the DATA2 paper demo notebook.
"""
import pytest

from reference_cases import DATA2_CASES
from conftest import DATA2_DIR
from _trajectory_compare import run_and_compare


@pytest.mark.parametrize("case", DATA2_CASES, ids=lambda c: c["case_id"])
def test_data2_trajectory_matches_reference(case, lib, references):
    """One forward-simulation regression test per DATA2 reference case."""
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

    data_path = DATA2_DIR / case["data_file"]
    if not data_path.exists():
        pytest.skip(f"data file not on disk: {data_path}")

    run_and_compare(lib, case, ref, data_path)
