#!/usr/bin/env python3
"""Restore pruned reproduction figures from archived tarballs.

Default behavior is dry-run. Use --apply to actually extract.
"""

from __future__ import annotations

import argparse
import shutil
import tarfile
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Restore UnifiedFramework/DATA3/results/reproduction/<run_id>/figures from archive.")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2], help="Repository root path.")
    parser.add_argument("--run-id", required=True, help="Run ID to restore (e.g., 20260221-004354).")
    parser.add_argument(
        "--repro-root",
        type=Path,
        default=Path("UnifiedFramework/DATA3/results/reproduction"),
        help="Reproduction root relative to repo-root.",
    )
    parser.add_argument(
        "--archive-root",
        type=Path,
        default=Path("UnifiedFramework/DATA3/results/reproduction/_archives"),
        help="Archive root relative to repo-root.",
    )
    parser.add_argument(
        "--archive-path",
        type=Path,
        default=None,
        help="Optional explicit archive path; overrides archive-root/run-id default.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace existing run figures/ if present.")
    parser.add_argument("--apply", action="store_true", help="Apply restore (default is dry-run).")
    return parser.parse_args()


def _safe_members(tf: tarfile.TarFile, expected_prefix: str):
    for member in tf.getmembers():
        # Restrict extraction to run_id/figures subtree only.
        if not member.name.startswith(expected_prefix):
            continue
        # Reject absolute or traversal paths.
        p = Path(member.name)
        if p.is_absolute() or ".." in p.parts:
            continue
        yield member


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    repro_root = (repo_root / args.repro_root).resolve()
    archive_root = (repo_root / args.archive_root).resolve()
    run_id = str(args.run_id)

    archive_path = args.archive_path.resolve() if args.archive_path else (archive_root / f"{run_id}-figures.tar.gz").resolve()
    run_dir = repro_root / run_id
    figures_dir = run_dir / "figures"
    dry_run = not bool(args.apply)

    print(f"[restore] dry_run={dry_run}")
    print(f"[restore] run_id={run_id}")
    print(f"[restore] archive={archive_path}")
    print(f"[restore] target={figures_dir}")

    if not archive_path.exists():
        raise SystemExit(f"[restore] archive not found: {archive_path}")
    if not run_dir.exists():
        raise SystemExit(f"[restore] run directory not found: {run_dir}")
    if figures_dir.exists() and not args.overwrite:
        raise SystemExit("[restore] run figures/ already exists. Re-run with --overwrite to replace.")

    expected_prefix = f"{run_id}/figures"
    with tarfile.open(archive_path, "r:gz") as tf:
        members = list(_safe_members(tf, expected_prefix=expected_prefix))
        file_count = sum(1 for m in members if m.isfile())
        print(f"[restore] files_to_extract={file_count}")
        if dry_run:
            return

        if figures_dir.exists() and args.overwrite:
            shutil.rmtree(figures_dir)

        # Archive uses arcname "<run_id>/figures", so extract to repro_root.
        tf.extractall(path=repro_root, members=members)

    print("[restore] completed")


if __name__ == "__main__":
    main()
