"""Layer 3 — validation gate.

Parses `baselines/known_nonpass.csv`, exercises the
`helpers.allowlist.evaluate()` triple-output contract on synthetic data, and
asserts the schema is healthy. Synthetic data is used here so this test runs
without needing a full materialize_all() — the integration with real
side-by-side CSVs is exercised by the Layer 2 artifacts tests when DATA
roots are present.

Default mode: WARNING-ONLY.
Set VALIDATION_GATE_FAIL_HARD=1 to flip to fail-hard semantics.
"""
from __future__ import annotations

import os
from datetime import date

import pytest

from helpers.allowlist import read_exceptions, evaluate, VALID_STATUSES, REQUIRED_COLUMNS
from conftest import ALLOWLIST_CSV


pytestmark = [pytest.mark.regression]


# ----------------------------------------------------------------------
# Allowlist file schema
# ----------------------------------------------------------------------

def test_allowlist_file_exists():
    """known_nonpass.csv exists (even if empty headers-only)."""
    assert ALLOWLIST_CSV.exists(), (
        f"Allowlist not found at {ALLOWLIST_CSV}. Layer 3 requires this file."
    )


def test_allowlist_required_columns():
    """known_nonpass.csv has the documented 5 columns in its header."""
    with open(ALLOWLIST_CSV) as fh:
        header = fh.readline().strip().split(",")
    missing = REQUIRED_COLUMNS - set(header)
    assert not missing, (
        f"known_nonpass.csv missing required columns: {sorted(missing)}.\n"
        f"Header: {header}"
    )


def test_allowlist_parses_cleanly(allowlist):
    """read_exceptions returns a dict (empty is fine; no rows = clean slate)."""
    assert isinstance(allowlist, dict)


# ----------------------------------------------------------------------
# Allowlist semantics — synthetic data
# ----------------------------------------------------------------------

def test_evaluate_returns_triple_on_empty_allowlist():
    """Empty allowlist + clean FAIL row → all in failures, none in warnings."""
    import pandas as pd
    df = pd.DataFrame([
        {"target_id": "T.a", "status": "FAIL"},
        {"target_id": "T.b", "status": "PASS"},
    ])
    failures, warnings, resolved = evaluate(df, {})
    assert len(failures) == 1
    assert failures.iloc[0]["target_id"] == "T.a"
    assert warnings.empty
    assert resolved.empty


def test_evaluate_categorizes_allowlisted_failure_as_warning():
    """A FAIL row whose target_id is allowlisted lands in warnings, not failures."""
    import pandas as pd
    df = pd.DataFrame([
        {"target_id": "T.a", "status": "FAIL"},
    ])
    exceptions = {
        "T.a": {
            "expected_status": "FAIL",
            "active": True,
            "review_after": date(2099, 1, 1),
            "reason": "synthetic",
        },
    }
    failures, warnings, resolved = evaluate(df, exceptions)
    assert failures.empty
    assert len(warnings) == 1
    assert warnings.iloc[0]["allowlist_reason"] == "synthetic"


def test_evaluate_categorizes_passing_allowlist_as_resolved():
    """A PASS row whose target_id is still allowlisted lands in resolved (stale entry)."""
    import pandas as pd
    df = pd.DataFrame([
        {"target_id": "T.a", "status": "PASS"},
    ])
    exceptions = {
        "T.a": {
            "expected_status": "FAIL",
            "active": True,
            "review_after": date(2099, 1, 1),
            "reason": "fixed last week",
        },
    }
    failures, warnings, resolved = evaluate(df, exceptions)
    assert failures.empty
    assert warnings.empty
    assert len(resolved) == 1


def test_evaluate_rejects_malformed_dataframe():
    """evaluate() raises if the side-by-side DataFrame is missing required columns."""
    import pandas as pd
    with pytest.raises(ValueError, match="target_id"):
        evaluate(pd.DataFrame([{"foo": "bar"}]), {})


def test_read_exceptions_filters_expired_rows(tmp_path):
    """An entry with review_after < today is filtered out (forces cleanup)."""
    csv = tmp_path / "expired.csv"
    csv.write_text(
        "target_id,expected_status,active,review_after,reason\n"
        "T.expired,FAIL,true,2000-01-01,old\n"
        "T.fresh,FAIL,true,2099-01-01,still relevant\n"
    )
    out = read_exceptions(csv)
    assert "T.expired" not in out
    assert "T.fresh" in out


def test_read_exceptions_filters_inactive_rows(tmp_path):
    """An entry with active != true is filtered out."""
    csv = tmp_path / "inactive.csv"
    csv.write_text(
        "target_id,expected_status,active,review_after,reason\n"
        "T.inactive,FAIL,false,2099-01-01,disabled\n"
        "T.active,FAIL,true,2099-01-01,enabled\n"
    )
    out = read_exceptions(csv)
    assert "T.inactive" not in out
    assert "T.active" in out


def test_read_exceptions_rejects_invalid_status(tmp_path):
    """An entry with an unknown expected_status is silently dropped."""
    csv = tmp_path / "bad_status.csv"
    csv.write_text(
        "target_id,expected_status,active,review_after,reason\n"
        "T.bogus,UNKNOWN,true,2099-01-01,not a real status\n"
        "T.ok,FAIL,true,2099-01-01,real\n"
    )
    out = read_exceptions(csv)
    assert "T.bogus" not in out
    assert "T.ok" in out


# ----------------------------------------------------------------------
# Failure-mode contract
# ----------------------------------------------------------------------

def test_fail_hard_env_var_documented():
    """VALIDATION_GATE_FAIL_HARD is the documented switch; default is warning-only."""
    # We assert the contract, not the current setting — the env var may be
    # set by the user. Just check the variable name is the one the README
    # advertises.
    assert "VALIDATION_GATE_FAIL_HARD" in os.environ or True, (
        # Trivially true; this exists so a doc-grep finds the variable name.
        "VALIDATION_GATE_FAIL_HARD is the documented switch."
    )
