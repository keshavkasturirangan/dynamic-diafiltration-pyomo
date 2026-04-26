#!/usr/bin/env python3
"""Validate simulation-validation figure mappings against extracted published-paper figures."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Iterable

import cv2
import numpy as np
import pandas as pd
from PIL import Image, ImageOps


REPO_ROOT = Path(__file__).resolve().parents[5]
SOURCE_MAP = REPO_ROOT / "UnifiedFramework" / "DATA3" / "docs" / "validation" / "target_notebook_source_map.csv"
PAGE_INDEX = REPO_ROOT / "UnifiedFramework" / "DATA3" / "docs" / "validation" / "target_pdf_page_index.csv"
PAPER_ROOT = (
    REPO_ROOT
    / "UnifiedFramework"
    / "DATA3"
    / "results"
    / "reproduction"
    / "20260306-023356-paper-pdf-extract"
    / "pdf_extract"
)
OUT_ROOT = (
    REPO_ROOT
    / "UnifiedFramework"
    / "DATA3"
    / "results"
    / "simulation_validation"
)
COMPOSITE_ROOT = OUT_ROOT / "composites"
REPORT_CSV = OUT_ROOT / "simulation_validation_figure_report.csv"
REPORT_MD = OUT_ROOT / "simulation_validation_figure_report.md"

PREFERRED_DIRS = [
    REPO_ROOT / "UnifiedFramework" / "DATA3" / "figures",
    REPO_ROOT / "DATA1_matlab",
    REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "reproduction" / "20260221-data1-notebook" / "figures" / "data1_notebook_regen",
    REPO_ROOT,
]

FIGURE_LAYOUT_OVERRIDES = {
    "D1-M-F2": [(2, 2)],
    "D1-M-F4": [(1, 3)],
    "D1-M-F5": [(1, 3)],
    "D1-M-F6": [(1, 3)],
    "D1-S-FS2": [(1, 2), (2, 1)],
    "D1-S-FS3": [(1, 2), (2, 1)],
    "D1-S-FS4": [(1, 3)],
    "D1-S-FS5": [(1, 2), (2, 1)],
    "D1-S-FS6": [(1, 3)],
    "D1-S-FS7": [(3, 4), (4, 3)],
    "D2-M-F7": [(2, 2)],
    "D2-M-F9": [(1, 2), (2, 1)],
    "D2-S-FS2": [(2, 2)],
    "D2-S-FS3": [(2, 2)],
    "D2-S-FS4": [(2, 2)],
    "D2-S-FS5": [(2, 2)],
    "D2-S-FS6": [(2, 2)],
    "D2-S-FS7": [(1, 2), (2, 1)],
    "D2-S-FS8": [(2, 2)],
}

PASS_THRESHOLD = 0.26


def _trim_white(im: Image.Image) -> Image.Image:
    gray = ImageOps.grayscale(im)
    arr = np.asarray(gray, dtype=np.uint8)
    mask = arr < 248
    if not mask.any():
        return gray
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    return gray.crop((int(cols[0]), int(rows[0]), int(cols[-1]) + 1, int(rows[-1]) + 1))


def _resolve_artifact_paths(patterns: Iterable[str]) -> list[Path]:
    resolved: list[Path] = []
    for pattern in patterns:
        pattern = pattern.strip()
        if not pattern:
            continue
        if "*" not in pattern:
            chosen = None
            for base in PREFERRED_DIRS:
                direct = base / pattern
                if direct.exists():
                    chosen = direct
                    break
            if chosen is None:
                matches = sorted(REPO_ROOT.rglob(pattern))
                if matches:
                    chosen = matches[0]
            if chosen is not None:
                resolved.append(chosen)
            continue

        matches: list[Path] = []
        for base in PREFERRED_DIRS:
            base_matches = sorted(base.glob(pattern))
            if base_matches:
                matches = base_matches
                break
        if not matches:
            matches = sorted(REPO_ROOT.rglob(pattern))
        unique = []
        seen = set()
        for match in matches:
            key = str(match.resolve())
            if key in seen:
                continue
            seen.add(key)
            unique.append(match)
        resolved.extend(unique)
    return resolved


def _extract_image_patterns(raw_field: str) -> list[str]:
    patterns = re.findall(r"[A-Za-z0-9_\-.*]+\.png", raw_field)
    patterns += re.findall(r"[A-Za-z0-9_\-.*]+\.jpg", raw_field)
    return patterns


def _candidate_layouts(n_images: int, target_id: str) -> list[tuple[int, int]]:
    if target_id in FIGURE_LAYOUT_OVERRIDES:
        return FIGURE_LAYOUT_OVERRIDES[target_id]
    if n_images <= 1:
        return [(1, 1)]
    layouts = {(1, n_images), (n_images, 1)}
    root = int(math.sqrt(n_images))
    for rows in range(2, root + 2):
        cols = math.ceil(n_images / rows)
        layouts.add((rows, cols))
        layouts.add((cols, rows))
    return sorted(layouts, key=lambda x: abs(x[0] - x[1]))


def _compose_images(images: list[Path], rows: int, cols: int, out_path: Path) -> None:
    pil_images = [_trim_white(Image.open(path).convert("RGB")) for path in images]
    widths = [im.width for im in pil_images]
    heights = [im.height for im in pil_images]
    cell_w = max(widths)
    cell_h = max(heights)
    pad = 24
    canvas = Image.new("RGB", (cols * cell_w + (cols + 1) * pad, rows * cell_h + (rows + 1) * pad), "white")

    for idx, im in enumerate(pil_images):
        row = idx // cols
        col = idx % cols
        x0 = pad + col * (cell_w + pad) + (cell_w - im.width) // 2
        y0 = pad + row * (cell_h + pad) + (cell_h - im.height) // 2
        canvas.paste(im, (x0, y0))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_path)


def _similarity_score(a_path: Path, b_path: Path) -> float:
    a_img = np.asarray(_trim_white(Image.open(a_path)), dtype=np.uint8)
    b_img = np.asarray(_trim_white(Image.open(b_path)), dtype=np.uint8)
    a_img = cv2.resize(a_img, (900, 900), interpolation=cv2.INTER_AREA)
    b_img = cv2.resize(b_img, (900, 900), interpolation=cv2.INTER_AREA)
    a_mask = (255 - a_img) / 255.0
    b_mask = (255 - b_img) / 255.0
    a_mask = (a_mask > 0.08).astype(np.float32)
    b_mask = (b_mask > 0.08).astype(np.float32)
    intersection = float((a_mask * b_mask).sum())
    union = float(np.clip(a_mask + b_mask, 0, 1).sum())
    iou = intersection / union if union else 0.0
    intensity = 1.0 - float(np.mean(np.abs(a_img.astype(np.float32) - b_img.astype(np.float32))) / 255.0)
    return 0.75 * iou + 0.25 * intensity


def _paper_images_for_pdf_key(pdf_key: str) -> list[Path]:
    return sorted(
        path
        for path in (PAPER_ROOT / pdf_key / "images").glob("*")
        if path.suffix.lower() in {".png", ".jpg", ".jpeg"}
    )


def _figure_rows() -> pd.DataFrame:
    source_map = pd.read_csv(SOURCE_MAP)
    page_index = pd.read_csv(PAGE_INDEX)
    merged = source_map.merge(page_index[["target_id", "pdf_key"]], on="target_id", how="inner")
    mask = merged["artifact_type"].astype(str).str.contains("Figure", case=False, na=False)
    return merged.loc[mask].copy()


def run_validation() -> pd.DataFrame:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    COMPOSITE_ROOT.mkdir(parents=True, exist_ok=True)
    rows = []
    for rec in _figure_rows().to_dict(orient="records"):
        target_id = str(rec["target_id"])
        patterns = _extract_image_patterns(str(rec["mapped_artifacts_or_signals"]))
        artifacts = _resolve_artifact_paths(patterns)
        paper_images = _paper_images_for_pdf_key(str(rec["pdf_key"]))
        if not artifacts or not paper_images:
            rows.append(
                {
                    "target_id": target_id,
                    "paper_item": rec["paper_item"],
                    "artifact_count": len(artifacts),
                    "paper_image_count": len(paper_images),
                    "best_score": np.nan,
                    "threshold": PASS_THRESHOLD,
                    "status": "MISSING_VALUE",
                    "best_layout": "",
                    "best_paper_image": "",
                    "composite_path": "",
                }
            )
            continue

        best = None
        for rows_count, cols_count in _candidate_layouts(len(artifacts), target_id):
            if rows_count * cols_count < len(artifacts):
                continue
            composite_path = COMPOSITE_ROOT / f"{target_id}_{rows_count}x{cols_count}.png"
            _compose_images(artifacts, rows_count, cols_count, composite_path)
            for paper_image in paper_images:
                score = _similarity_score(composite_path, paper_image)
                if best is None or score > best["best_score"]:
                    best = {
                        "target_id": target_id,
                        "paper_item": rec["paper_item"],
                        "artifact_count": len(artifacts),
                        "paper_image_count": len(paper_images),
                        "best_score": score,
                        "threshold": PASS_THRESHOLD,
                        "status": "PASS" if score >= PASS_THRESHOLD else "FAIL",
                        "best_layout": f"{rows_count}x{cols_count}",
                        "best_paper_image": str(paper_image.relative_to(REPO_ROOT)),
                        "composite_path": str(composite_path.relative_to(REPO_ROOT)),
                    }
        if best is None:
            rows.append(
                {
                    "target_id": target_id,
                    "paper_item": rec["paper_item"],
                    "artifact_count": len(artifacts),
                    "paper_image_count": len(paper_images),
                    "best_score": np.nan,
                    "threshold": PASS_THRESHOLD,
                    "status": "MISSING_VALUE",
                    "best_layout": "",
                    "best_paper_image": "",
                    "composite_path": "",
                }
            )
        else:
            rows.append(best)

    report = pd.DataFrame(rows).sort_values("target_id").reset_index(drop=True)
    report.to_csv(REPORT_CSV, index=False, quoting=csv.QUOTE_MINIMAL)
    REPORT_MD.write_text(
        report.to_markdown(index=False) + "\n",
        encoding="utf-8",
    )
    return report


def main() -> int:
    run_validation()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
