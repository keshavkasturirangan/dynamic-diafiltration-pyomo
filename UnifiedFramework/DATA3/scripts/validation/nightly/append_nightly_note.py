#!/usr/bin/env python3
"""Append one structured nightly note entry to UnifiedFramework/DATA3/docs/validation/nightly/logs/nightly_test_notes.md."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Append a pytest note entry.")
    p.add_argument("--notes-md", type=Path, default=Path("UnifiedFramework/DATA3/docs/validation/nightly/logs/nightly_test_notes.md"))
    p.add_argument("--run-id", required=True, help="Run identifier, e.g., 20260306-nightly-local")
    p.add_argument("--context", default="Nightly validation run")
    p.add_argument("--gate-status", choices=("PASS", "FAIL", "SKIPPED"), required=True)
    p.add_argument("--gate-summary", default="")
    p.add_argument("--pytest-status", choices=("PASS", "FAIL", "TIMEOUT", "SKIPPED", "UNKNOWN"), default="UNKNOWN")
    p.add_argument("--action", default="")
    p.add_argument("--note", action="append", default=[], help="Additional bullet note; can repeat.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    ts = datetime.now().strftime("%Y-%m-%d %H:%M local")

    lines = [
        f"## {ts} | run_id: {args.run_id}",
        f"- Context: {args.context}",
        f"- Validation gate: {args.gate_status}.",
        f"- Gate summary: {args.gate_summary or 'n/a'}",
        f"- Full nightly pytest: {args.pytest_status}.",
        f"- Action: {args.action or 'n/a'}",
    ]
    for n in args.note:
        lines.append(f"- Note: {n}")
    lines.append("")

    path = args.notes_md
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# Pytest Notes\n\nThis running log captures key outcomes from nightly validation runs.\n\n", encoding="utf-8")

    with path.open("a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(f"[nightly-note] appended: {path}")


if __name__ == "__main__":
    main()
