"""pytest fixtures + helpers shared by the DATA1/DATA2 regression tests.

The regression strategy
-----------------------
For each reference case (a specific .mat dataset + a specific fit recipe), we
ran the fit ONCE with the committed code state, captured the fitted
parameters into `regression_references.json`, and committed that JSON.

Each pytest then:
  1. Re-runs the same fit with the current library
  2. Reads the reference parameters from the JSON
  3. Asserts every parameter matches within REL_TOLERANCE

This catches the case where a library edit silently shifts a DATA1 or DATA2
result. If the reference values themselves need to be intentionally updated
(e.g. after a deliberate physics change), re-run `_capture_references.py`.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest


# Tolerance used by every regression test.
# 1e-2 = 1 % — tight enough to catch silent drift, loose enough to absorb
# legitimate solver float-noise variation across Pyomo/IPOPT versions.
REL_TOLERANCE = 1e-2

# Repo root and library import path.
# Path(__file__) → tests/conftest.py
#   parents[0] = tests/
#   parents[1] = refactored_codes_v1/
#   parents[2] = REPO_ROOT
HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
LIB_DIR = REPO_ROOT / "refactored_codes_v1"
# Make the library importable AND make tests/ importable (so test files can
# `from reference_cases import ...` and `from _trajectory_compare import ...`).
for p in (str(LIB_DIR), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

# Working directory matters: the DATA1/DATA2 loaders use relative paths like
# "data_library/data_stru-...mat" so we need to be at the repo root.
os.chdir(str(REPO_ROOT))


# Path constants for the .mat reference data.
DATA1_DIR = REPO_ROOT / "legacy" / "data1_matlab" / "data"
DATA2_DIR = REPO_ROOT / "legacy" / "data1_matlab" / "data_library"

# The reference-values JSON lives next to the test files.
REF_JSON = Path(__file__).resolve().parent / "regression_references.json"


@pytest.fixture(scope="session")
def lib():
    """Import the library once per test session."""
    import refactored_ucb_library as _lib
    return _lib


@pytest.fixture(scope="session")
def references():
    """Load the captured reference parameter values."""
    if not REF_JSON.exists():
        pytest.fail(
            f"Reference values not captured yet. Run:\n"
            f"   python {Path(__file__).resolve().parent}/_capture_references.py"
        )
    with open(REF_JSON) as fh:
        return json.load(fh)


def assert_params_match(fitted: dict, reference: dict, tol: float = REL_TOLERANCE):
    """Assert each numeric parameter in `reference` matches `fitted` within
    relative tolerance `tol`. Compares the keys that exist in BOTH dicts.

    Skips dict-valued entries (e.g. nested theta sub-dicts) and string fields.
    """
    mismatches = []
    for key, ref_val in reference.items():
        if key not in fitted:
            mismatches.append(f"  missing key in fitted: {key!r}")
            continue
        if isinstance(ref_val, dict) or isinstance(ref_val, str):
            continue
        if ref_val is None:
            continue
        try:
            ref_f = float(ref_val)
            fit_f = float(fitted[key])
        except (TypeError, ValueError):
            continue
        # Tolerance is RELATIVE except for values near zero, where we fall
        # back to a small absolute tolerance.
        denom = max(abs(ref_f), 1e-9)
        rel = abs(fit_f - ref_f) / denom
        if rel > tol:
            mismatches.append(
                f"  {key!r}: fitted={fit_f:.6g}, ref={ref_f:.6g}, "
                f"rel_diff={rel:.3%} (tol={tol:.1%})"
            )
    if mismatches:
        pytest.fail(
            "Parameters drifted beyond tolerance:\n" + "\n".join(mismatches)
        )
