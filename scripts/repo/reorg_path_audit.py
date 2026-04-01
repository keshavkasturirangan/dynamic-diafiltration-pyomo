#!/usr/bin/env python3
"""Audit hard-coded repository paths to support staged reorganization.

Outputs a markdown report with literal path-hit counts by file.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PATTERNS = [
    "UnifiedFramework/DATA3",
    "DATA1_matlab",
    "DATA2_visualization.ipynb",
    "results/reproduction",
    "docs/validation",
    "run_cross_verification.py",
    "utility.py",
]

TEXT_SUFFIXES = {
    ".py",
    ".md",
    ".txt",
    ".csv",
    ".yml",
    ".yaml",
    ".ini",
    ".toml",
    ".tex",
    ".json",
    ".ipynb",
}


@dataclass
class Hit:
    relpath: str
    pattern: str
    count: int


def iter_files(repo_root: Path) -> Iterable[Path]:
    roots = [
        repo_root / "README.md",
        repo_root / "docs",
        repo_root / "UnifiedFramework" / "DATA3" / "scripts",
        repo_root / "UnifiedFramework" / "DATA3" / "tests",
    ]
    for root in roots:
        if not root.exists():
            continue
        if root.is_file():
            yield root
            continue
        for p in root.rglob("*"):
            if p.is_file() and p.suffix.lower() in TEXT_SUFFIXES:
                yield p

    for p in repo_root.glob("*.py"):
        if p.is_file():
            yield p


def count_hits(repo_root: Path) -> list[Hit]:
    hits: list[Hit] = []
    seen: set[Path] = set()
    for p in iter_files(repo_root):
        rel = p.relative_to(repo_root)
        # Ignore generated audit outputs to avoid self-inflating counts.
        if str(rel).startswith("docs/repo_reorg/PATH_AUDIT_REPORT"):
            continue
        if p in seen:
            continue
        seen.add(p)
        try:
            text = p.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        for pattern in PATTERNS:
            c = text.count(pattern)
            if c:
                hits.append(Hit(relpath=str(p.relative_to(repo_root)), pattern=pattern, count=c))
    return hits


def render_report(repo_root: Path, hits: list[Hit]) -> str:
    by_pattern: dict[str, int] = {k: 0 for k in PATTERNS}
    for h in hits:
        by_pattern[h.pattern] += h.count

    lines: list[str] = []
    lines.append("# Reorg Path Audit Report")
    lines.append("")
    lines.append("Generated for staged repository reorganization.")
    lines.append("")
    lines.append("## Pattern Totals")
    lines.append("")
    for p in PATTERNS:
        lines.append(f"- `{p}`: {by_pattern[p]}")

    lines.append("")
    lines.append("## Top Files By Path-Literal Hits")
    lines.append("")

    file_totals: dict[str, int] = {}
    for h in hits:
        file_totals[h.relpath] = file_totals.get(h.relpath, 0) + h.count

    for relpath, total in sorted(file_totals.items(), key=lambda kv: kv[1], reverse=True)[:40]:
        lines.append(f"- `{relpath}`: {total}")

    lines.append("")
    lines.append("## Detailed Hits")
    lines.append("")
    lines.append("| file | pattern | count |")
    lines.append("|---|---|---:|")
    for h in sorted(hits, key=lambda x: (x.relpath, x.pattern)):
        lines.append(f"| `{h.relpath}` | `{h.pattern}` | {h.count} |")

    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- High hit-count files should be migrated last or updated with adapter constants first.")
    lines.append("- Keep `UnifiedFramework/DATA3` stable until the path-hit count is near zero for moved items.")
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Audit hard-coded paths before repository reorganization")
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    p.add_argument(
        "--out",
        type=Path,
        default=Path("docs/repo_reorg/PATH_AUDIT_REPORT.md"),
        help="Output markdown path relative to repo root unless absolute.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    out = args.out if args.out.is_absolute() else repo_root / args.out
    hits = count_hits(repo_root)
    report = render_report(repo_root, hits)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report, encoding="utf-8")
    print(f"[reorg-audit] wrote: {out}")
    print(f"[reorg-audit] hit_rows: {len(hits)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
