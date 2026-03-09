#!/usr/bin/env python3
"""Storage-aware pruning for UnifiedFramework/DATA3/results/reproduction artifacts.

Policy:
- Keep full artifacts (including per-run `figures/`) for the newest N runs.
- Keep full artifacts for any run referenced by validation breadcrumbs/docs.
- For all other runs, prune per-run `figures/` by either:
  - archiving to .tar.gz then deleting `figures/` (default), or
  - deleting `figures/` directly.
- Always keep run_metadata.json, tables/, and summaries/.

This script is dry-run by default; pass --apply to perform changes.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple


@dataclass(frozen=True)
class RunAction:
    run_id: str
    keep_full: bool
    reason: str
    figures_path: str
    figures_size_mb: float
    action: str
    archive_path: str | None


RUN_DIR_PATTERN = re.compile(r"^20\d{6}(?:-\d{6}|-data1-notebook)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prune reproduction artifacts in a storage-aware way.")
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Repository root path.",
    )
    parser.add_argument(
        "--repro-root",
        type=Path,
        default=Path("UnifiedFramework/DATA3/results/reproduction"),
        help="Reproduction output root (relative to repo-root by default).",
    )
    parser.add_argument(
        "--keep-full",
        type=int,
        default=3,
        help="Number of newest runs to keep fully intact.",
    )
    parser.add_argument(
        "--mode",
        choices=("archive", "remove"),
        default="archive",
        help="How to prune per-run figures/ directories for older runs.",
    )
    parser.add_argument(
        "--archive-root",
        type=Path,
        default=Path("UnifiedFramework/DATA3/results/reproduction/_archives"),
        help="Archive output dir when --mode archive.",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply changes (default is dry-run).",
    )
    return parser.parse_args()


def dir_size_bytes(path: Path) -> int:
    total = 0
    for root, _, files in os.walk(path):
        root_path = Path(root)
        for name in files:
            fp = root_path / name
            try:
                total += fp.stat().st_size
            except OSError:
                continue
    return total


def list_run_dirs(repro_root: Path) -> List[Path]:
    if not repro_root.exists():
        return []
    out: List[Path] = []
    for p in repro_root.iterdir():
        if p.is_dir() and RUN_DIR_PATTERN.match(p.name):
            out.append(p)
    return out


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return ""


def protected_runs_from_docs(repo_root: Path, run_ids: Iterable[str]) -> Set[str]:
    """Mark runs referenced anywhere in docs/validation as protected."""
    docs_dir = repo_root / "UnifiedFramework/DATA3/docs/validation"
    if not docs_dir.exists():
        return set()

    run_ids_set = set(run_ids)
    protected: Set[str] = set()
    for file_path in docs_dir.rglob("*"):
        if not file_path.is_file():
            continue
        text = _read_text(file_path)
        if not text:
            continue
        for run_id in run_ids_set:
            if run_id in text:
                protected.add(run_id)
    return protected


def protected_runs_from_ledger(repo_root: Path, run_ids: Iterable[str]) -> Set[str]:
    """Also parse continuity_ledger.csv columns explicitly."""
    ledger = repo_root / "UnifiedFramework/DATA3/docs/validation/continuity_ledger.csv"
    if not ledger.exists():
        return set()
    run_ids_set = set(run_ids)
    protected: Set[str] = set()
    try:
        with ledger.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                rid = (row.get("run_id_or_doc") or "").strip()
                if rid in run_ids_set:
                    protected.add(rid)
                for field in ("primary_evidence_paths", "notes"):
                    val = row.get(field) or ""
                    for run_id in run_ids_set:
                        if run_id in val:
                            protected.add(run_id)
    except OSError:
        pass
    return protected


def newest_run_ids(run_dirs: List[Path], keep_full: int) -> Set[str]:
    if keep_full <= 0:
        return set()
    sorted_by_mtime = sorted(
        run_dirs,
        key=lambda p: p.stat().st_mtime if p.exists() else 0.0,
        reverse=True,
    )
    return {p.name for p in sorted_by_mtime[:keep_full]}


def plan_actions(
    run_dirs: List[Path],
    protected_ids: Set[str],
    newest_ids: Set[str],
    mode: str,
    archive_root_abs: Path,
) -> Tuple[List[RunAction], float]:
    actions: List[RunAction] = []
    reclaimable_mb = 0.0

    for run_dir in sorted(run_dirs):
        run_id = run_dir.name
        figs = run_dir / "figures"
        fig_size_mb = (dir_size_bytes(figs) / (1024 * 1024)) if figs.exists() else 0.0

        keep_full = False
        reason = ""
        action = "none"
        archive_path: str | None = None

        if run_id in protected_ids:
            keep_full = True
            reason = "referenced_in_docs"
        elif run_id in newest_ids:
            keep_full = True
            reason = "latest_n"
        else:
            keep_full = False
            reason = "prune_eligible"
            if figs.exists():
                reclaimable_mb += fig_size_mb
                action = "archive_and_remove_figures" if mode == "archive" else "remove_figures"
                if mode == "archive":
                    archive_path = str((archive_root_abs / f"{run_id}-figures.tar.gz").resolve())
            else:
                action = "nothing_to_prune"

        actions.append(
            RunAction(
                run_id=run_id,
                keep_full=keep_full,
                reason=reason,
                figures_path=str(figs.resolve()),
                figures_size_mb=round(fig_size_mb, 2),
                action=action,
                archive_path=archive_path,
            )
        )
    return actions, round(reclaimable_mb, 2)


def apply_actions(actions: List[RunAction], repro_root_abs: Path, mode: str, archive_root_abs: Path) -> None:
    if mode == "archive":
        archive_root_abs.mkdir(parents=True, exist_ok=True)

    for a in actions:
        if a.keep_full:
            continue
        run_dir = repro_root_abs / a.run_id
        figs = run_dir / "figures"
        if not figs.exists():
            continue
        if mode == "archive":
            archive_path = archive_root_abs / f"{a.run_id}-figures.tar.gz"
            with tarfile.open(archive_path, "w:gz") as tf:
                tf.add(figs, arcname=f"{a.run_id}/figures")
            shutil.rmtree(figs)
        else:
            shutil.rmtree(figs)


def print_report(
    actions: List[RunAction],
    dry_run: bool,
    mode: str,
    keep_full: int,
    protected_count: int,
    reclaimable_mb: float,
) -> None:
    print(f"[prune] dry_run={dry_run}")
    print(f"[prune] mode={mode}")
    print(f"[prune] keep_full={keep_full}")
    print(f"[prune] protected_runs={protected_count}")
    print(f"[prune] reclaimable_figures_mb={reclaimable_mb}")
    print("[prune] actions:")
    for a in actions:
        print(
            "  - "
            f"{a.run_id}: keep_full={a.keep_full}, reason={a.reason}, "
            f"figures_mb={a.figures_size_mb}, action={a.action}"
        )
        if a.archive_path:
            print(f"    archive={a.archive_path}")


def main() -> None:
    args = parse_args()
    repo_root_abs = args.repo_root.resolve()
    repro_root_abs = (repo_root_abs / args.repro_root).resolve()
    archive_root_abs = (repo_root_abs / args.archive_root).resolve()

    run_dirs = list_run_dirs(repro_root_abs)
    run_ids = [p.name for p in run_dirs]

    protected_ids = set()
    protected_ids.update(protected_runs_from_docs(repo_root_abs, run_ids))
    protected_ids.update(protected_runs_from_ledger(repo_root_abs, run_ids))
    newest_ids = newest_run_ids(run_dirs, keep_full=int(args.keep_full))

    actions, reclaimable_mb = plan_actions(
        run_dirs=run_dirs,
        protected_ids=protected_ids,
        newest_ids=newest_ids,
        mode=args.mode,
        archive_root_abs=archive_root_abs,
    )

    dry_run = not bool(args.apply)
    if not dry_run:
        apply_actions(actions, repro_root_abs=repro_root_abs, mode=args.mode, archive_root_abs=archive_root_abs)

    print_report(
        actions=actions,
        dry_run=dry_run,
        mode=args.mode,
        keep_full=int(args.keep_full),
        protected_count=len(protected_ids),
        reclaimable_mb=reclaimable_mb,
    )


if __name__ == "__main__":
    main()
