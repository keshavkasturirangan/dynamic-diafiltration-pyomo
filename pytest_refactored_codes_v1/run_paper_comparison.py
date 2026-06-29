#!/usr/bin/env python3
"""Standalone paper-figure comparison.

For each DATA1 / DATA2 target in target_notebook_source_map.csv, scores the
currently-on-disk refactored output (PNG) against every paper figure crop /
page render in pdf_extract/. The best-scoring paper image is reported per
target.

  - Pure PIL + numpy (no scikit-image dep).
  - Read-only: no edits to refactored_codes_v1/ or to UnifiedFramework/.
  - Writes a per-target report under pytest_refactored_codes_v1/_runs/.

Run:
    cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo
    python3 pytest_refactored_codes_v1/run_paper_comparison.py

Metric: zero-mean normalized cross correlation (NCC) on grayscale, resized
to a common 256x256. NCC range is [-1, 1]; a perfect match is 1.0. For
panel-vs-page comparison (DATA1, no crops) typical good scores are 0.2-0.5
because the panel occupies only ~15% of the page. For panel-vs-crop
(DATA2 main + SI), 0.4-0.8 is typical for a true match.

The threshold of 0.26 used in the (April 17) `simulation_validation_figure_report.csv`
is preserved here for continuity, but the report includes the raw scores so
you can re-threshold without re-running.

See docs/PAPER_COMPARISON_LAYER.md for how this fits the suite's 4-layer
architecture.
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
from PIL import Image


# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------

HERE = Path(__file__).resolve().parent                           # pytest_refactored_codes_v1/
REPO_ROOT = HERE.parent
UF_VALIDATION = REPO_ROOT / "UnifiedFramework" / "DATA3" / "docs" / "validation"
LOCKED_MAPPING = UF_VALIDATION / "target_notebook_source_map.csv"
PAPER_ARTIFACTS = {
    "DATA1": REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "data1" / "notebook_figures",
    "DATA2": REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "data2" / "notebook_figures",
}
PDF_EXTRACT_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "reproduction" / "20260306-023356-paper-pdf-extract" / "pdf_extract"

# Mapping: dataset_scope + paper_item.startswith("Fig. S") tells us which
# subdir of pdf_extract to look in.
PDF_EXTRACT_SUBDIRS = {
    ("DATA1", False): PDF_EXTRACT_ROOT / "data1_main",   # main figures
    ("DATA1", True):  PDF_EXTRACT_ROOT / "data1_si",     # SI figures
    ("DATA2", False): PDF_EXTRACT_ROOT / "data2_main",
    ("DATA2", True):  PDF_EXTRACT_ROOT / "data2_si",
}

OUT_DIR = HERE / "_runs"
OUT_DIR.mkdir(exist_ok=True)
TIMESTAMP = time.strftime("%Y-%m-%d_%H%M%S")
OUT_CSV = OUT_DIR / f"paper_comparison_{TIMESTAMP}.csv"
OUT_MD = OUT_DIR / f"paper_comparison_{TIMESTAMP}.md"
OUT_LATEST_CSV = OUT_DIR / "paper_comparison_LATEST.csv"

DEFAULT_THRESHOLD = 0.26    # matches April 17 baseline
COMMON_SIZE = (256, 256)


# ----------------------------------------------------------------------
# Image loading + similarity
# ----------------------------------------------------------------------

_IMAGE_CACHE: dict = {}


def _load_grayscale(path: Path) -> np.ndarray | None:
    """Load and resize image to COMMON_SIZE, grayscale, float. Cached."""
    if path in _IMAGE_CACHE:
        return _IMAGE_CACHE[path]
    if not path.exists():
        _IMAGE_CACHE[path] = None
        return None
    try:
        img = Image.open(path).convert("L").resize(COMMON_SIZE, Image.LANCZOS)
        arr = np.asarray(img, dtype=float)
        _IMAGE_CACHE[path] = arr
        return arr
    except Exception as exc:
        print(f"  [warn] failed to load {path.name}: {exc}", file=sys.stderr)
        _IMAGE_CACHE[path] = None
        return None


def ncc(a: np.ndarray, b: np.ndarray) -> float:
    """Zero-mean normalized cross correlation. Returns float in [-1, 1]."""
    a = a - a.mean()
    b = b - b.mean()
    denom = float(np.sqrt((a ** 2).sum() * (b ** 2).sum()))
    if denom == 0:
        return 0.0
    return float((a * b).sum() / denom)


# ----------------------------------------------------------------------
# Locked mapping loader
# ----------------------------------------------------------------------

def load_locked_rows() -> list[dict]:
    rows = []
    with open(LOCKED_MAPPING, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            artifacts_raw = (r.get("mapped_artifacts_or_signals") or "").strip()
            artifacts = [a.strip() for a in artifacts_raw.split(";") if a.strip()]
            rows.append({
                "target_id": (r.get("target_id") or "").strip(),
                "paper_item": (r.get("paper_item") or "").strip(),
                "artifact_type": (r.get("artifact_type") or "").strip(),
                "dataset_scope": (r.get("dataset_scope") or "").strip(),
                "mapped_artifacts": artifacts,
            })
    return rows


def pdf_extract_for(row: dict) -> Path:
    """Resolve which pdf_extract subdir holds the paper images for this row."""
    is_si = "S" in row["paper_item"].split(".")[1] if "." in row["paper_item"] else False
    # "Fig. S5" → SI; "Fig. 5" → main. Look for " S" after "Fig."
    is_si = " S" in row["paper_item"] or row["paper_item"].startswith("Fig. S") or "S" in row["target_id"].split("-")[1]
    key = (row["dataset_scope"], is_si)
    return PDF_EXTRACT_SUBDIRS.get(key, PDF_EXTRACT_ROOT / "data1_main")


def list_paper_images(extract_dir: Path) -> list[Path]:
    """Return all paper images for the extract subdir.

    Prefer cropped images if available (data2_main, data2_si, data1_si).
    Fall back to full pages (data1_main has only these).
    """
    images = []
    if (extract_dir / "images").exists():
        images.extend(sorted((extract_dir / "images").glob("*.png")))
        images.extend(sorted((extract_dir / "images").glob("*.jpg")))
    if not images and (extract_dir / "pages").exists():
        images.extend(sorted((extract_dir / "pages").glob("*.png")))
    return images


# ----------------------------------------------------------------------
# Per-target comparison
# ----------------------------------------------------------------------

def score_one_artifact(produced: Path, paper_images: list[Path]) -> tuple[float, Path | None]:
    """Compare one produced PNG against every paper image; return (best, path)."""
    arr_p = _load_grayscale(produced)
    if arr_p is None:
        return float("nan"), None
    best_score = -2.0
    best_path = None
    for paper_img in paper_images:
        arr_b = _load_grayscale(paper_img)
        if arr_b is None:
            continue
        s = ncc(arr_p, arr_b)
        if s > best_score:
            best_score = s
            best_path = paper_img
    return best_score, best_path


def evaluate_target(row: dict) -> list[dict]:
    """Score each mapped artifact in this row. Returns a list of dicts (one per artifact)."""
    paper_dir = pdf_extract_for(row)
    paper_images = list_paper_images(paper_dir)
    paper_kind = "image-crops" if (paper_dir / "images").exists() and any((paper_dir / "images").glob("*")) else "page-renders"

    artifact_base = PAPER_ARTIFACTS.get(row["dataset_scope"])
    if artifact_base is None:
        return []

    out = []
    for name in row["mapped_artifacts"]:
        produced = artifact_base / name
        if not produced.exists():
            out.append({
                "target_id": row["target_id"],
                "paper_item": row["paper_item"],
                "artifact": name,
                "produced_path": str(produced),
                "produced_present": False,
                "paper_kind": paper_kind,
                "score": "",
                "matched_paper_image": "",
                "verdict": "MISSING",
            })
            continue
        score, matched = score_one_artifact(produced, paper_images)
        if not np.isfinite(score):
            verdict = "ERROR"
        elif score >= DEFAULT_THRESHOLD:
            verdict = "MATCHES_PAPER"
        else:
            verdict = "WEAK_MATCH"
        out.append({
            "target_id": row["target_id"],
            "paper_item": row["paper_item"],
            "artifact": name,
            "produced_path": str(produced),
            "produced_present": True,
            "paper_kind": paper_kind,
            "score": f"{score:.4f}" if np.isfinite(score) else "",
            "matched_paper_image": str(matched) if matched else "",
            "verdict": verdict,
        })
    return out


# ----------------------------------------------------------------------
# Main + report writer
# ----------------------------------------------------------------------

def main():
    print(f"=== Paper-comparison run @ {TIMESTAMP} ===")
    print(f"  locked mapping: {LOCKED_MAPPING}")
    print(f"  pdf_extract:    {PDF_EXTRACT_ROOT}")
    print(f"  output:         {OUT_CSV}")
    print()

    locked = load_locked_rows()
    print(f"  loaded {len(locked)} mapped targets ({sum(1 for r in locked if r['dataset_scope']=='DATA1')} DATA1, "
          f"{sum(1 for r in locked if r['dataset_scope']=='DATA2')} DATA2)")
    print()

    all_rows = []
    t0 = time.time()
    for i, row in enumerate(locked, 1):
        if row["dataset_scope"] not in ("DATA1", "DATA2"):
            continue
        results = evaluate_target(row)
        all_rows.extend(results)
        scored = [r for r in results if r["score"]]
        present = [r for r in results if r["produced_present"]]
        if scored:
            best = max(float(r["score"]) for r in scored)
            print(f"  [{i:>2d}/{len(locked)}] {row['target_id']:<10s}  "
                  f"{len(present)}/{len(results)} present, best NCC = {best:.3f}")
        else:
            print(f"  [{i:>2d}/{len(locked)}] {row['target_id']:<10s}  "
                  f"{len(present)}/{len(results)} present  (no scoreable artifacts)")

    print()
    print(f"  elapsed: {time.time()-t0:.1f}s")
    print()

    # Write CSV
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "target_id", "paper_item", "artifact", "produced_path", "produced_present",
            "paper_kind", "score", "matched_paper_image", "verdict",
        ])
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    # Latest symlink-equivalent: copy
    import shutil
    shutil.copy(OUT_CSV, OUT_LATEST_CSV)
    print(f"  wrote CSV → {OUT_CSV.relative_to(REPO_ROOT)}")
    print(f"  latest    → {OUT_LATEST_CSV.relative_to(REPO_ROOT)}")

    # Write Markdown
    md_lines = [f"# Paper Comparison Report — {TIMESTAMP}", ""]
    md_lines.append(f"- Source: {LOCKED_MAPPING.relative_to(REPO_ROOT)}")
    md_lines.append(f"- Metric: Normalized Cross Correlation (NCC), range [-1, 1]")
    md_lines.append(f"- Threshold: {DEFAULT_THRESHOLD} (matches the April 17 baseline)")
    md_lines.append(f"- Total artifacts scored: {sum(1 for r in all_rows if r['score'])}")
    md_lines.append(f"- Matches paper:  {sum(1 for r in all_rows if r['verdict']=='MATCHES_PAPER')}")
    md_lines.append(f"- Weak match:     {sum(1 for r in all_rows if r['verdict']=='WEAK_MATCH')}")
    md_lines.append(f"- Missing on disk: {sum(1 for r in all_rows if r['verdict']=='MISSING')}")
    md_lines.append("")
    md_lines.append("## Per-artifact scores")
    md_lines.append("")
    md_lines.append("| target_id | artifact | produced? | paper kind | NCC | verdict |")
    md_lines.append("|---|---|---|---|---|---|")
    for r in all_rows:
        score_disp = r["score"] if r["score"] else "—"
        present_disp = "✓" if r["produced_present"] else "✗"
        md_lines.append(f"| `{r['target_id']}` | `{r['artifact']}` | {present_disp} | {r['paper_kind']} | {score_disp} | **{r['verdict']}** |")
    OUT_MD.write_text("\n".join(md_lines))
    print(f"  wrote MD  → {OUT_MD.relative_to(REPO_ROOT)}")

    # Console summary
    n_match = sum(1 for r in all_rows if r["verdict"] == "MATCHES_PAPER")
    n_weak = sum(1 for r in all_rows if r["verdict"] == "WEAK_MATCH")
    n_miss = sum(1 for r in all_rows if r["verdict"] == "MISSING")
    n_err  = sum(1 for r in all_rows if r["verdict"] == "ERROR")
    total = len(all_rows)
    print()
    print(f"=== Verdict summary ===")
    print(f"  total artifacts:  {total}")
    print(f"  matches paper:    {n_match} ({100*n_match/max(total,1):.0f}%)")
    print(f"  weak match:       {n_weak}")
    print(f"  missing:          {n_miss}")
    print(f"  error:            {n_err}")


if __name__ == "__main__":
    main()
