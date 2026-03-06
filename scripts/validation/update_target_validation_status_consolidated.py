#!/usr/bin/env python3
"""Update consolidated target validation status with explicit numeric semantics.

Semantics:
- If a target_id has no numeric rows in tolerance CSV -> NOT_APPLICABLE
- Else aggregate row statuses with precedence:
  FAIL > PASS_WITH_EXPLANATION > MISSING_VALUE > PASS
- evidence_run_id is updated for targets that have numeric rows.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Refresh docs/validation/target_validation_status_consolidated.csv")
    parser.add_argument(
        "--consolidated-csv",
        type=Path,
        default=Path("docs/validation/target_validation_status_consolidated.csv"),
        help="Consolidated target status CSV path.",
    )
    parser.add_argument(
        "--tolerance-csv",
        type=Path,
        default=Path("docs/validation/tolerance_eval_refresh_20260221.csv"),
        help="Row-level tolerance evaluation CSV path.",
    )
    parser.add_argument(
        "--evidence-run-id",
        default="20260221-refresh",
        help="Evidence run label written for numeric targets found in tolerance CSV.",
    )
    return parser.parse_args()


def aggregate_target_status(statuses: pd.Series) -> str:
    s = {str(x).strip() for x in statuses.dropna().tolist()}
    if "FAIL" in s:
        return "FAIL"
    if "PASS_WITH_EXPLANATION" in s:
        return "PASS_WITH_EXPLANATION"
    if "MISSING_VALUE" in s:
        return "MISSING_VALUE"
    if "PASS" in s:
        return "PASS"
    return "MISSING_VALUE"


def main() -> None:
    args = parse_args()
    consolidated = pd.read_csv(args.consolidated_csv)
    tolerance = pd.read_csv(args.tolerance_csv)

    required_cons = {
        "target_id",
        "numeric_validation_status",
        "evidence_run_id",
    }
    required_tol = {"target_id", "status"}
    missing_cons = required_cons - set(consolidated.columns)
    missing_tol = required_tol - set(tolerance.columns)
    if missing_cons:
        raise ValueError(f"Missing consolidated columns: {sorted(missing_cons)}")
    if missing_tol:
        raise ValueError(f"Missing tolerance columns: {sorted(missing_tol)}")

    by_target = tolerance.groupby("target_id", dropna=False)["status"].apply(aggregate_target_status)
    target_to_status = by_target.to_dict()

    def _new_numeric_status(target_id: str) -> str:
        if target_id not in target_to_status:
            return "NOT_APPLICABLE"
        return str(target_to_status[target_id])

    consolidated["numeric_validation_status"] = consolidated["target_id"].map(_new_numeric_status)

    has_numeric = consolidated["target_id"].isin(set(target_to_status.keys()))
    consolidated.loc[has_numeric, "evidence_run_id"] = args.evidence_run_id
    consolidated.loc[~has_numeric, "evidence_run_id"] = consolidated.loc[~has_numeric, "evidence_run_id"].fillna("")

    consolidated.to_csv(args.consolidated_csv, index=False)
    print(f"[status-update] updated: {args.consolidated_csv}")
    print(f"[status-update] numeric targets: {int(has_numeric.sum())}")
    print(f"[status-update] non-numeric targets: {int((~has_numeric).sum())}")
    not_applicable_targets = consolidated.loc[
        consolidated["numeric_validation_status"] == "NOT_APPLICABLE", "target_id"
    ].tolist()
    if not_applicable_targets:
        print(
            "[status-update][WARNING] NOT_APPLICABLE numeric targets "
            f"({len(not_applicable_targets)}): {', '.join(not_applicable_targets)}"
        )


if __name__ == "__main__":
    main()
