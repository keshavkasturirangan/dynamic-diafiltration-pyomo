#!/usr/bin/env python3
"""Nightly validation gate for consolidated target statuses.

Policy:
- Ignore NOT_APPLICABLE targets.
- Treat PASS and PASS_WITH_EXPLANATION as passing.
- For FAIL/MISSING_VALUE, allow only documented exceptions.
- Any unexpected FAIL/MISSING_VALUE fails the gate.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple


NON_PASS_STATUSES = {"FAIL", "MISSING_VALUE"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Nightly validation gate against consolidated target status CSV.")
    parser.add_argument(
        "--status-csv",
        type=Path,
        default=Path("docs/validation/target_validation_status_consolidated.csv"),
        help="Consolidated target status CSV.",
    )
    parser.add_argument(
        "--exceptions-csv",
        type=Path,
        default=Path("docs/validation/nightly_known_nonpass.csv"),
        help="Known non-pass exceptions allowlist CSV.",
    )
    return parser.parse_args()


def read_status_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing status CSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    required = {"target_id", "numeric_validation_status"}
    missing = required - set(rows[0].keys() if rows else [])
    if missing:
        raise ValueError(f"Missing required columns in status CSV: {sorted(missing)}")
    return rows


def read_exceptions(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return {}

    required = {"target_id", "expected_status", "active", "reason"}
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Missing required columns in exceptions CSV: {sorted(missing)}")

    exceptions: Dict[str, Dict[str, str]] = {}
    for row in rows:
        target_id = str(row.get("target_id", "")).strip()
        if not target_id:
            continue
        active = str(row.get("active", "true")).strip().lower()
        if active in {"0", "false", "no", "n"}:
            continue
        exceptions[target_id] = row
    return exceptions


def evaluate(
    rows: List[Dict[str, str]], exceptions: Dict[str, Dict[str, str]]
) -> Tuple[List[str], List[str], List[str], Counter]:
    failures: List[str] = []
    warnings: List[str] = []
    resolved: List[str] = []

    counts = Counter()
    seen_targets = set()

    for row in rows:
        target_id = str(row.get("target_id", "")).strip()
        status = str(row.get("numeric_validation_status", "")).strip()
        if not target_id:
            continue
        seen_targets.add(target_id)
        counts[status] += 1

        if status in {"PASS", "PASS_WITH_EXPLANATION", "NOT_APPLICABLE"}:
            continue

        if status in NON_PASS_STATUSES:
            exc = exceptions.get(target_id)
            if exc is None:
                failures.append(f"{target_id}: unexpected {status} (no exception entry)")
                continue

            expected = str(exc.get("expected_status", "ANY")).strip()
            reason = str(exc.get("reason", "")).strip()
            if expected in {"ANY", status}:
                warnings.append(f"{target_id}: allowed {status} ({reason})")
            else:
                failures.append(
                    f"{target_id}: status {status} does not match expected {expected} ({reason})"
                )
            continue

        failures.append(f"{target_id}: unknown status '{status}'")

    for target_id, exc in exceptions.items():
        if target_id not in seen_targets:
            warnings.append(f"{target_id}: listed in exceptions but missing from consolidated CSV")
            continue
        current = next(
            str(r.get("numeric_validation_status", "")).strip() for r in rows if str(r.get("target_id", "")).strip() == target_id
        )
        expected = str(exc.get("expected_status", "ANY")).strip()
        if current not in NON_PASS_STATUSES:
            resolved.append(f"{target_id}: now {current} (exception can likely be removed)")
        elif expected not in {"ANY", current}:
            # Already captured as failure above; leave here for completeness in summary context.
            pass

    return failures, warnings, resolved, counts


def main() -> None:
    args = parse_args()
    rows = read_status_rows(args.status_csv)
    exceptions = read_exceptions(args.exceptions_csv)

    failures, warnings, resolved, counts = evaluate(rows, exceptions)

    print("[nightly-gate] status counts:", dict(counts))
    print(f"[nightly-gate] exceptions active: {len(exceptions)}")

    if warnings:
        print("[nightly-gate] warnings:")
        for w in warnings:
            print("  -", w)

    if resolved:
        print("[nightly-gate] resolved exceptions:")
        for r in resolved:
            print("  -", r)

    if failures:
        print("[nightly-gate] failures:")
        for f in failures:
            print("  -", f)
        sys.exit(1)

    print("[nightly-gate] PASS: no unexpected FAIL/MISSING_VALUE targets.")


if __name__ == "__main__":
    main()
