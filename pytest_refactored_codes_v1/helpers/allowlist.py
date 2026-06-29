"""Validation-gate allowlist mechanism.

Parses `baselines/known_nonpass.csv` into a `{target_id: meta}` map and
evaluates a side-by-side DataFrame against it, producing a
`(failures, warnings, resolved)` triple.

Schema for known_nonpass.csv (one row per allowlisted target):

    target_id,expected_status,active,review_after,reason
    DATA1.fig2_A.vial3,FAIL,true,2026-09-01,known mass-trace drift; see commit abc123
    DATA2.calib_plots,NOT_APPLICABLE,true,2027-01-01,upstream calibration CSV pending

Semantics
---------
- `failures` = rows with status == "FAIL" that are NOT in the allowlist.
  These are unexpected regressions — fail the gate (in fail-hard mode) or
  emit warnings (in warning-only mode).
- `warnings` = rows with status == "FAIL" that ARE in the allowlist with a
  matching expected_status. These are documented non-passes; surfaced but
  not fatal.
- `resolved` = rows with status == "PASS" that are still in the allowlist.
  The allowlist entry is now stale and should be removed.

Filtering rules
---------------
- `active != "true"` → entry is soft-deleted (ignored).
- `review_after < today` → entry is expired (ignored, forces re-evaluation).
- `expected_status` must be one of {"FAIL", "NOT_APPLICABLE", "PASS_WITH_EXPLANATION"}.
"""
from __future__ import annotations

import csv
from datetime import date
from pathlib import Path
from typing import Optional


VALID_STATUSES = {"FAIL", "NOT_APPLICABLE", "PASS_WITH_EXPLANATION"}
REQUIRED_COLUMNS = {"target_id", "expected_status", "active", "review_after", "reason"}


def read_exceptions(path: Path, *, today: Optional[date] = None) -> dict:
    """Parse a known_nonpass.csv into {target_id: meta}.

    Rows where `active != "true"` or `review_after < today` are dropped silently.
    Returns an empty dict if the file is missing or has no active rows.
    """
    if today is None:
        today = date.today()
    if not Path(path).exists():
        return {}

    out: dict = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or not REQUIRED_COLUMNS.issubset(reader.fieldnames):
            missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
            raise ValueError(
                f"{path}: missing required columns: {sorted(missing)}"
            )

        for row in reader:
            target_id = (row.get("target_id") or "").strip()
            if not target_id:
                continue
            if (row.get("active") or "").strip().lower() != "true":
                continue
            review_after_raw = (row.get("review_after") or "").strip()
            try:
                review_after = date.fromisoformat(review_after_raw)
            except ValueError:
                # Bad date → treat as expired (forces cleanup).
                continue
            if review_after < today:
                continue
            expected = (row.get("expected_status") or "").strip().upper()
            if expected not in VALID_STATUSES:
                continue

            out[target_id] = {
                "expected_status": expected,
                "active": True,
                "review_after": review_after,
                "reason": (row.get("reason") or "").strip(),
            }
    return out


def evaluate(side_by_side, exceptions: dict):
    """Compare a side-by-side DataFrame to the allowlist.

    Expects the DataFrame to have at minimum two columns: `target_id` and
    `status`. Returns three DataFrames: (failures, warnings, resolved).
    """
    import pandas as pd  # local import — avoids pandas being a hard dep of conftest

    if "target_id" not in side_by_side.columns or "status" not in side_by_side.columns:
        raise ValueError(
            "evaluate(): side-by-side DataFrame must have 'target_id' and 'status' columns"
        )

    df = side_by_side.copy()
    df["target_id"] = df["target_id"].astype(str).str.strip()
    df["status"] = df["status"].astype(str).str.strip().str.upper()

    fail_mask = df["status"] == "FAIL"
    pass_mask = df["status"] == "PASS"
    allowlisted_mask = df["target_id"].isin(exceptions.keys())

    failures = df[fail_mask & ~allowlisted_mask].copy()
    warnings = df[fail_mask & allowlisted_mask].copy()
    resolved = df[pass_mask & allowlisted_mask].copy()

    # Annotate warnings with the reason from the allowlist.
    if not warnings.empty:
        warnings["allowlist_reason"] = warnings["target_id"].map(
            lambda t: exceptions.get(t, {}).get("reason", "")
        )
    if not resolved.empty:
        resolved["allowlist_reason"] = resolved["target_id"].map(
            lambda t: exceptions.get(t, {}).get("reason", "")
        )

    return failures, warnings, resolved
