"""Layer 4 (structural) — DATA2 regression tests.

Migrated from refactored_codes_v1/tests/test_data2_regression.py.
Mirrors test_data1_regression.py — same architecture (structural fingerprint
comparison via model_construct_inter, no solver run). DATA2 cases use
B_form=1 (numeric) and an explicit theta dict matching the DATA2 paper
demo notebook.
"""
from __future__ import annotations

import pytest

from helpers.reference_cases import DATA2_CASES
from helpers.trajectory_compare import run_and_compare
from conftest import DATA2_DIR


pytestmark = [pytest.mark.regression]


@pytest.mark.parametrize("case", DATA2_CASES, ids=lambda c: c["case_id"])
def test_data2_structural_matches_reference(case, lib, references):
    """One structural-fingerprint regression test per DATA2 reference case."""
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
