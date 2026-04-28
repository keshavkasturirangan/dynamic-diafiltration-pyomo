#!/usr/bin/env python3
"""Post-process workflow outputs into paper-style figures only.

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
- a full run folder with a figures subfolder, or
- a flat output folder containing image files
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
import sys

import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from refactored_ucb_library import (
    get_campaign_root,
    list_supported_campaigns,
    run_paper_reproduction,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts"
DEFAULT_REPRODUCTION_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "reproduction"
PAPER_CAMPAIGNS = ("DATA1", "DATA2")


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


def _has_figure_bucket(source_dir: Path) -> bool:
    """Return True if the source folder already has a figures subfolder."""
    return (source_dir / "figures").is_dir()


def _looks_like_raw_input_root(source_dir: Path) -> bool:
    """Heuristically decide whether a folder is a raw input root."""
    if _has_figure_bucket(source_dir):
        return False
    if (source_dir / "experiment space").exists() or (source_dir / "sigma sensitivity").exists():
        return True
    if any(source_dir.glob("data_stru-*.mat")) or any(source_dir.rglob("data_stru-*.mat")):
        return True
    if any(source_dir.glob("*.mat")):
        return True
    return False


def _iter_figure_files(source_dir: Path):
    """Yield image files from either a figures bucket or a flat output tree."""
    if _has_figure_bucket(source_dir):
        source_dir = source_dir / "figures"

    for path in source_dir.rglob("*"):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        if suffix not in {".png", ".jpg", ".jpeg", ".svg", ".pdf", ".eps"}:
            continue
        yield path, path.relative_to(source_dir)


def _stage_figure_dir(source_dir: Path, stage_dir: Path) -> list[Path]:
    """Copy the source figures into a normalized paper-artifacts tree."""
    stage_dir.mkdir(parents=True, exist_ok=True)
    staged: list[Path] = []

    for src, rel_path in _iter_figure_files(source_dir):
        dst = stage_dir / "figures" / rel_path
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        staged.append(dst)

    return staged


def build_paper_artifacts(
    campaign: str,
    source_dir: Path | None = None,
    output_root: Path | None = None,
    mode: str = "auto",
    show: bool = True,
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
    stage_mode = mode == "stage" or (mode == "auto" and _has_figure_bucket(source_dir))

    if mode == "auto" and not raw_mode and not stage_mode:
        # Prefer the end-to-end paper workflow when the source points at a
        # campaign root, even if it does not contain a finished run layout yet.
        campaign_root = Path(get_campaign_root(campaign))
        if source_dir.resolve() == campaign_root.resolve() or campaign_root in source_dir.resolve().parents:
            raw_mode = True

    campaign_root = output_root / campaign.lower() / source_dir.name
    run_dir = campaign_root / "workflow_run"
    stage_dir = campaign_root / "staged"

    if raw_mode:
        run_dir.mkdir(parents=True, exist_ok=True)
        run_paper_reproduction(campaign, data_root=source_dir, save_dir=run_dir)
        staged = _stage_figure_dir(run_dir, stage_dir)
    elif stage_mode:
        staged = _stage_figure_dir(source_dir, stage_dir)
    else:
        # If the user did not specify a source folder, default to the raw
        # campaign root and run the workflow end-to-end.
        if mode == "auto":
            run_dir.mkdir(parents=True, exist_ok=True)
            run_paper_reproduction(campaign, data_root=source_dir, save_dir=run_dir)
            staged = _stage_figure_dir(run_dir, stage_dir)
        else:
            raise ValueError(
                "Could not infer whether the source is raw inputs or a finished run. "
                "Pass --mode raw or --mode stage explicitly."
            )

    # Return the key files for easy printing in the terminal.
    outputs = [str(p) for p in sorted((stage_dir / "figures").rglob("*")) if p.is_file()]
    if show:
        _show_figure_files([Path(p) for p in outputs])
    return sorted(set(outputs))


def _show_figure_files(paths: list[Path]) -> None:
    """Open the saved figure files in matplotlib windows."""
    figure_paths = [p for p in paths if p.suffix.lower() in {".png", ".jpg", ".jpeg", ".svg", ".pdf", ".eps"}]
    for idx, fig_path in enumerate(sorted(figure_paths), start=1):
        try:
            image = plt.imread(fig_path)
        except Exception:
            continue
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(image)
        ax.axis("off")
        ax.set_title(f"{idx}: {fig_path.name}", fontsize=10)
    if figure_paths:
        plt.show()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize run outputs into paper-style figures and tables for DATA1/DATA2.",
    )
    parser.add_argument(
        "--campaign",
        default="DATA1",
        choices=["ALL", *list_supported_campaigns()],
        help="Which campaign to post-process. Default is DATA1.",
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
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open the generated figures in matplotlib windows.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    campaigns = list(PAPER_CAMPAIGNS) if args.campaign == "ALL" else [args.campaign]
    all_outputs: list[str] = []
    for campaign in campaigns:
        if args.campaign == "ALL":
            print(f"\nOpening {campaign} figures. Close the figure windows to continue.")
        outputs = build_paper_artifacts(
            campaign,
            source_dir=Path(args.source_dir) if args.source_dir else None,
            output_root=Path(args.output_root) if args.output_root else None,
            mode=args.mode,
            show=not args.no_show,
        )
        all_outputs.extend(outputs)
        print(f"\n{campaign}: created {len(outputs)} staged or derived files")
        for item in outputs:
            print(f"  - {item}")
        if args.campaign == "ALL" and campaign != campaigns[-1]:
            print("\nNext up: DATA2. Close the DATA1 windows to continue.")

    print(f"\nDone. Total files written or collected: {len(set(all_outputs))}")


if __name__ == "__main__":
    main()
