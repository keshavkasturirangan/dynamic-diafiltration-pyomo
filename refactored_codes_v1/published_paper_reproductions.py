#!/usr/bin/env python3
"""Post-process workflow outputs into paper-style figures and tables.

This script does not refit the model. It takes the outputs produced by:
- refactored_ucb_runfile.py
- refactored_ucb_library.py
- conductivity_paper.py

and organizes them into the same paper-style bundles used in the DATA1 and
DATA2 publications.

The workflow stays intentionally simple:
- DATA1 = utility-style model results + DATA1 paper plotting recipe
- DATA2 = utility-style model results + DATA2 paper plotting recipe

The script can work from either:
- a full run folder with figures/tables/summaries subfolders, or
- a flat output folder containing PNG/CSV/JSON files
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import shutil
from pathlib import Path
from typing import Iterable

from refactored_ucb_library import (
    build_colored_paper_composites_with_tables,
    get_campaign_root,
    list_supported_campaigns,
    run_paper_reproduction,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts"
DEFAULT_REPRODUCTION_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "reproduction"


def _latest_child_dir(root: Path) -> Path | None:
    """Return the newest child directory under `root`, if any."""
    if not root.exists():
        return None
    candidates = [p for p in root.iterdir() if p.is_dir()]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def _default_source_dir(campaign: str) -> Path | None:
    """Choose the most likely source folder for one campaign."""
    campaign = str(campaign or "DATA1").upper()

    # Prefer the raw input root for end-to-end reproduction.
    campaign_root = Path(get_campaign_root(campaign))
    if campaign_root.exists():
        return campaign_root

    # Fall back to the flat figure output folder used by the runner.
    flat_dir = REPO_ROOT / "UnifiedFramework" / "DATA3" / "figures" / campaign.lower()
    if flat_dir.exists():
        return flat_dir

    # If the runner folder is not present, prefer the latest reproduction run
    # because it usually already has the full figures/tables/summaries layout.
    latest_repro = _latest_child_dir(DEFAULT_REPRODUCTION_ROOT)
    if latest_repro is not None:
        return latest_repro

    return None


def _direct_bucket_files(source_dir: Path) -> bool:
    """Return True if the source folder already has figures/tables/summaries."""
    return any((source_dir / name).is_dir() for name in ("figures", "tables", "summaries"))


def _looks_like_raw_input_root(source_dir: Path) -> bool:
    """Heuristically decide whether a folder is a raw input root."""
    if _direct_bucket_files(source_dir):
        return False
    if (source_dir / "experiment space").exists() or (source_dir / "sigma sensitivity").exists():
        return True
    if any(source_dir.glob("data_stru-*.mat")) or any(source_dir.rglob("data_stru-*.mat")):
        return True
    if any(source_dir.glob("*.mat")):
        return True
    return False


def _iter_source_files(source_dir: Path) -> Iterable[tuple[Path, str, Path]]:
    """Yield (source_path, bucket_name, bucket_relative_path) triples.

    If the source directory already has figures/tables/summaries subfolders,
    we preserve those buckets. Otherwise we bucket files by suffix so a flat
    output directory can still be normalized into the paper layout.
    """
    if _direct_bucket_files(source_dir):
        for bucket in ("figures", "tables", "summaries"):
            bucket_root = source_dir / bucket
            if not bucket_root.exists():
                continue
            for path in bucket_root.rglob("*"):
                if path.is_file():
                    yield path, bucket, path.relative_to(bucket_root)
        return

    for path in source_dir.rglob("*"):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        if suffix in {".png", ".jpg", ".jpeg", ".svg", ".pdf", ".eps"}:
            bucket = "figures"
        elif suffix in {".csv", ".tsv", ".xlsx", ".xls"}:
            bucket = "tables"
        elif suffix in {".json", ".md", ".txt"}:
            bucket = "summaries"
        else:
            continue
        yield path, bucket, path.relative_to(source_dir)


def _stage_source_dir(source_dir: Path, stage_dir: Path) -> list[dict[str, str]]:
    """Copy the source outputs into a normalized paper-artifacts tree."""
    stage_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, str]] = []

    for src, bucket, rel_path in _iter_source_files(source_dir):
        dst = stage_dir / bucket / rel_path
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        manifest.append(
            {
                "source": str(src),
                "bucket": bucket,
                "staged": str(dst),
            }
        )

    return manifest


def _write_manifest(manifest: list[dict[str, str]], out_dir: Path, campaign: str, source_dir: Path) -> Path:
    """Write a small CSV/Markdown manifest for the staged paper bundle."""
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / f"{campaign.lower()}_artifact_manifest.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["source", "bucket", "staged"])
        writer.writeheader()
        writer.writerows(manifest)

    md_path = out_dir / f"{campaign.lower()}_artifact_manifest.md"
    with md_path.open("w") as fh:
        fh.write(f"# {campaign} artifact manifest\n\n")
        fh.write(f"Source directory: `{source_dir}`\n\n")
        fh.write(f"Total staged files: {len(manifest)}\n\n")
        fh.write("| Source | Bucket | Staged |\n")
        fh.write("| --- | --- | --- |\n")
        for row in manifest:
            fh.write(f"| {row['source']} | {row['bucket']} | {row['staged']} |\n")

    return md_path


def build_paper_artifacts(
    campaign: str,
    source_dir: Path | None = None,
    output_root: Path | None = None,
    mode: str = "auto",
) -> list[str]:
    """Build the paper-style bundles either from raw inputs or from staged outputs."""
    campaign = str(campaign or "DATA1").upper()
    if campaign not in list_supported_campaigns():
        raise ValueError(f"Unknown campaign '{campaign}'. Known campaigns: {list_supported_campaigns()}")

    source_dir = Path(source_dir).expanduser().resolve() if source_dir is not None else _default_source_dir(campaign)
    if source_dir is None or not source_dir.exists():
        raise FileNotFoundError(
            f"Could not find a source directory for {campaign}. "
            "Pass --source-dir explicitly or run the workflow first."
        )

    output_root = Path(output_root).expanduser().resolve() if output_root is not None else DEFAULT_OUTPUT_ROOT
    mode = str(mode or "auto").lower()
    if mode not in {"auto", "raw", "stage"}:
        raise ValueError("mode must be one of: auto, raw, stage")

    raw_mode = mode == "raw" or (mode == "auto" and _looks_like_raw_input_root(source_dir))
    stage_mode = mode == "stage" or (mode == "auto" and _direct_bucket_files(source_dir))

    campaign_root = output_root / campaign.lower() / source_dir.name
    run_dir = campaign_root / "workflow_run"
    stage_dir = campaign_root / "staged"

    if raw_mode:
        run_dir.mkdir(parents=True, exist_ok=True)
        produced = run_paper_reproduction(campaign, data_root=source_dir, save_dir=run_dir)
        manifest = _stage_source_dir(run_dir, stage_dir)
        _write_manifest(manifest, campaign_root, campaign, source_dir)
    elif stage_mode:
        manifest = _stage_source_dir(source_dir, stage_dir)
        _write_manifest(manifest, campaign_root, campaign, source_dir)
        produced = [str(p) for p in source_dir.rglob("*") if p.is_file()]
    else:
        raise ValueError(
            "Could not infer whether the source is raw inputs or a finished run. "
            "Pass --mode raw or --mode stage explicitly."
        )

    # Build the colored paper composites from the staged files.
    composites_dir = campaign_root / "composites_colored_tables"
    composite_outputs = build_colored_paper_composites_with_tables(
        stage_dir,
        campaign=campaign,
        out_dir=composites_dir,
    )

    # Return the key files for easy printing in the terminal.
    outputs = [str(p) for p in sorted(campaign_root.rglob("*")) if p.is_file()]
    outputs.extend(produced)
    outputs.extend(composite_outputs)
    return sorted(set(outputs))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize run outputs into paper-style figures and tables for DATA1/DATA2.",
    )
    parser.add_argument(
        "--campaign",
        default="ALL",
        choices=["ALL", *list_supported_campaigns()],
        help="Which campaign to post-process.",
    )
    parser.add_argument(
        "--source-dir",
        default=None,
        help="A run folder to post-process. If omitted, the newest reproduction folder is used.",
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help=f"Where to write the paper-artifacts tree. Defaults to {DEFAULT_OUTPUT_ROOT}.",
    )
    parser.add_argument(
        "--mode",
        default="auto",
        choices=["auto", "raw", "stage"],
        help="Use raw inputs, stage an existing run folder, or auto-detect.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    campaigns = list_supported_campaigns() if args.campaign == "ALL" else [args.campaign]
    all_outputs: list[str] = []
    for campaign in campaigns:
        outputs = build_paper_artifacts(
            campaign,
            source_dir=Path(args.source_dir) if args.source_dir else None,
            output_root=Path(args.output_root) if args.output_root else None,
            mode=args.mode,
        )
        all_outputs.extend(outputs)
        print(f"\n{campaign}: created {len(outputs)} staged or derived files")
        for item in outputs:
            print(f"  - {item}")

    print(f"\nDone. Total files written or collected: {len(set(all_outputs))}")


if __name__ == "__main__":
    main()
