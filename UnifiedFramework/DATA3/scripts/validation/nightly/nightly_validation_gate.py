#!/usr/bin/env python3
"""Nightly validation gate for consolidated target statuses.

Policy:
- Ignore NOT_APPLICABLE targets.
- Treat PASS and PASS_WITH_EXPLANATION as passing.
- For FAIL/MISSING_VALUE, allow only documented exceptions.
- Any unexpected FAIL/MISSING_VALUE fails the gate.

This gate can also execute nightly pytest regression checks and persist a
lightweight summary CSV so automation runs always publish current test status.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple
from xml.etree import ElementTree as ET


NON_PASS_STATUSES = {"FAIL", "MISSING_VALUE"}


@dataclass
class PytestSummary:
    exit_code: int
    tests: int
    failures: int
    errors: int
    skipped: int
    passed: int
    failed_cases: List[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Nightly validation gate against consolidated target status CSV.")
    parser.add_argument(
        "--status-csv",
        type=Path,
        default=Path("UnifiedFramework/DATA3/docs/validation/target_validation_status_consolidated.csv"),
        help="Consolidated target status CSV.",
    )
    parser.add_argument(
        "--exceptions-csv",
        type=Path,
        default=Path("UnifiedFramework/DATA3/docs/validation/nightly/config/known_nonpass.csv"),
        help="Known non-pass exceptions allowlist CSV.",
    )
    parser.add_argument(
        "--skip-pytest",
        action="store_true",
        help="Skip running nightly pytest checks.",
    )
    parser.add_argument(
        "--pytest-target",
        default="UnifiedFramework/DATA3/tests/regression",
        help="Pytest target path/module for nightly checks.",
    )
    parser.add_argument(
        "--pytest-marker",
        default="nightly",
        help="Pytest marker expression for nightly checks.",
    )
    parser.add_argument(
        "--pytest-junit-xml",
        type=Path,
        default=Path("UnifiedFramework/DATA3/docs/validation/nightly/pytest_latest.xml"),
        help="JUnit XML output path for pytest run.",
    )
    parser.add_argument(
        "--pytest-summary-csv",
        type=Path,
        default=Path("UnifiedFramework/DATA3/docs/validation/nightly/pytest_latest_summary.csv"),
        help="CSV output path for pytest summary metrics.",
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
            pass

    return failures, warnings, resolved, counts


def _sum_attr(node: ET.Element, name: str) -> int:
    total = 0
    for testsuite in node.iter("testsuite"):
        value = testsuite.attrib.get(name, "0")
        try:
            total += int(float(value))
        except (TypeError, ValueError):
            continue
    if total == 0 and node.tag == "testsuite":
        value = node.attrib.get(name, "0")
        try:
            total = int(float(value))
        except (TypeError, ValueError):
            total = 0
    return total


def parse_junit_summary(junit_xml: Path, exit_code: int) -> PytestSummary:
    if not junit_xml.exists():
        return PytestSummary(
            exit_code=exit_code,
            tests=0,
            failures=1,
            errors=0,
            skipped=0,
            passed=0,
            failed_cases=["pytest: junit xml not generated"],
        )

    root = ET.parse(junit_xml).getroot()
    tests = _sum_attr(root, "tests")
    failures = _sum_attr(root, "failures")
    errors = _sum_attr(root, "errors")
    skipped = _sum_attr(root, "skipped")
    passed = max(tests - failures - errors - skipped, 0)

    failed_cases: List[str] = []
    for testcase in root.iter("testcase"):
        if testcase.find("failure") is None and testcase.find("error") is None:
            continue
        classname = testcase.attrib.get("classname", "")
        name = testcase.attrib.get("name", "")
        case_id = "::".join(part for part in [classname, name] if part)
        failed_cases.append(case_id or "unknown_testcase")

    return PytestSummary(
        exit_code=exit_code,
        tests=tests,
        failures=failures,
        errors=errors,
        skipped=skipped,
        passed=passed,
        failed_cases=failed_cases,
    )


def write_pytest_summary_csv(path: Path, summary: PytestSummary) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "timestamp_utc",
        "exit_code",
        "tests",
        "passed",
        "failures",
        "errors",
        "skipped",
        "failed_cases",
    ]
    row = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "exit_code": summary.exit_code,
        "tests": summary.tests,
        "passed": summary.passed,
        "failures": summary.failures,
        "errors": summary.errors,
        "skipped": summary.skipped,
        "failed_cases": " | ".join(summary.failed_cases),
    }
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerow(row)


def run_pytest(pytest_target: str, pytest_marker: str, junit_xml: Path) -> PytestSummary:
    junit_xml.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "-m",
        pytest_marker,
        pytest_target,
        "--junitxml",
        str(junit_xml),
    ]

    print("[nightly-gate] running pytest:", " ".join(cmd))
    completed = subprocess.run(cmd, text=True, capture_output=True)

    if completed.stdout.strip():
        print("[nightly-gate][pytest] stdout:")
        print(completed.stdout.strip())
    if completed.stderr.strip():
        print("[nightly-gate][pytest] stderr:")
        print(completed.stderr.strip())

    return parse_junit_summary(junit_xml, completed.returncode)


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

    if args.skip_pytest:
        print("[nightly-gate] pytest: skipped (--skip-pytest)")
    else:
        pytest_summary = run_pytest(args.pytest_target, args.pytest_marker, args.pytest_junit_xml)
        write_pytest_summary_csv(args.pytest_summary_csv, pytest_summary)
        print(
            "[nightly-gate] pytest summary:",
            {
                "tests": pytest_summary.tests,
                "passed": pytest_summary.passed,
                "failures": pytest_summary.failures,
                "errors": pytest_summary.errors,
                "skipped": pytest_summary.skipped,
                "exit_code": pytest_summary.exit_code,
            },
        )
        print(f"[nightly-gate] pytest summary csv: {args.pytest_summary_csv}")
        if pytest_summary.exit_code != 0:
            failures.append("pytest nightly checks failed; see junit xml/summary csv for details")
            for case in pytest_summary.failed_cases[:10]:
                failures.append(f"pytest failure: {case}")

    if failures:
        print("[nightly-gate] failures:")
        for f in failures:
            print("  -", f)
        sys.exit(1)

    print("[nightly-gate] PASS: no unexpected FAIL/MISSING_VALUE targets.")


if __name__ == "__main__":
    main()
