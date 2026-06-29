#!/usr/bin/env python3
"""Composite-based paper-figure comparison.

Reconstructs the April 17 methodology (validate_simulation_validation_data_files.py
is no longer on disk) so we can refresh the per-target similarity scores
against the CURRENT state of paper_artifacts/data{1,2}/notebook_figures/.

For each target in target_notebook_source_map.csv (DATA1 + DATA2):
  1. Collect produced PNGs that currently exist on disk
  2. Try every reasonable grid layout (1x1, 1x2, 2x1, 1x3, 3x1, 2x2, 1x4, 4x1, 3x4, 4x3, ...)
  3. Build the composite via PIL (resize each panel to a common cell size, tile)
  4. NCC the composite against every paper figure crop/page in the appropriate
     pdf_extract subdir
  5. Pick the (layout, paper_image) pair with the highest NCC
  6. Save the best composite to _runs/run-<timestamp>/composites/

Pure PIL + numpy. Read-only on refactored_codes_v1/ and on UnifiedFramework/.

Run:
    cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo
    python3 pytest_refactored_codes_v1/run_paper_comparison_composite.py

The output is comparable to UnifiedFramework/DATA3/results/simulation_validation/
simulation_validation_figure_report.csv (April 17 baseline) so we can report
drift target-by-target.

Threshold 0.26 matches the April 17 baseline. A target with NCC >= 0.26 is
considered to MATCH the published paper figure.
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

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
UF_VALIDATION = REPO_ROOT / "UnifiedFramework" / "DATA3" / "docs" / "validation"
LOCKED_MAPPING = UF_VALIDATION / "target_notebook_source_map.csv"
APRIL_REPORT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "simulation_validation" / "simulation_validation_figure_report.csv"

PAPER_ARTIFACTS = {
    "DATA1": REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "data1" / "notebook_figures",
    "DATA2": REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "data2" / "notebook_figures",
}
PDF_EXTRACT_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "reproduction" / "20260306-023356-paper-pdf-extract" / "pdf_extract"
PDF_EXTRACT_SUBDIRS = {
    ("DATA1", False): PDF_EXTRACT_ROOT / "data1_main",
    ("DATA1", True):  PDF_EXTRACT_ROOT / "data1_si",
    ("DATA2", False): PDF_EXTRACT_ROOT / "data2_main",
    ("DATA2", True):  PDF_EXTRACT_ROOT / "data2_si",
}

TIMESTAMP = time.strftime("%Y-%m-%d_%H%M%S")
OUT_RUN_DIR = HERE / "_runs" / f"composite_{TIMESTAMP}"
OUT_RUN_DIR.mkdir(parents=True, exist_ok=True)
OUT_COMPOSITES = OUT_RUN_DIR / "composites"
OUT_COMPOSITES.mkdir(exist_ok=True)
OUT_CSV = OUT_RUN_DIR / "report.csv"
OUT_MD = OUT_RUN_DIR / "report.md"
OUT_LATEST_CSV = HERE / "_runs" / "composite_LATEST.csv"

DEFAULT_THRESHOLD = 0.26
CELL_SIZE = (400, 400)        # per-panel cell in the composite
COMPARE_SIZE = (256, 256)     # NCC target size


# ----------------------------------------------------------------------
# Layout enumeration
# ----------------------------------------------------------------------

def candidate_layouts(n: int) -> list[tuple[int, int]]:
    """Return reasonable (rows, cols) layouts for n artifacts.

    For small n we enumerate; for larger n we pick a few sensible options
    near the square root.
    """
    if n <= 0:
        return []
    if n == 1: return [(1, 1)]
    if n == 2: return [(1, 2), (2, 1)]
    if n == 3: return [(1, 3), (3, 1)]
    if n == 4: return [(2, 2), (1, 4), (4, 1)]
    if n == 5: return [(1, 5), (5, 1)]
    if n == 6: return [(2, 3), (3, 2), (1, 6), (6, 1)]
    if n == 8: return [(2, 4), (4, 2)]
    if n == 9: return [(3, 3)]
    if n == 12: return [(3, 4), (4, 3), (2, 6), (6, 2)]
    if n == 16: return [(4, 4)]
    # Generic fallback: try a few near-square shapes
    out = []
    for r in range(1, n + 1):
        if n % r == 0:
            c = n // r
            out.append((r, c))
    return out[:6] if out else [(1, n)]


# ----------------------------------------------------------------------
# Image loading + cache
# ----------------------------------------------------------------------

_GRAY_CACHE: dict = {}


def load_grayscale(path: Path, size=COMPARE_SIZE) -> np.ndarray | None:
    """Load, grayscale, resize to size — float array. Cached."""
    key = (str(path), size)
    if key in _GRAY_CACHE:
        return _GRAY_CACHE[key]
    if not path.exists():
        _GRAY_CACHE[key] = None
        return None
    try:
        img = Image.open(path).convert("L").resize(size, Image.LANCZOS)
        arr = np.asarray(img, dtype=float)
        _GRAY_CACHE[key] = arr
        return arr
    except Exception as exc:
        print(f"  [warn] failed to load {path.name}: {exc}", file=sys.stderr)
        _GRAY_CACHE[key] = None
        return None


def ncc(a: np.ndarray, b: np.ndarray) -> float:
    a = a - a.mean()
    b = b - b.mean()
    denom = float(np.sqrt((a ** 2).sum() * (b ** 2).sum()))
    if denom == 0:
        return 0.0
    return float((a * b).sum() / denom)


# ----------------------------------------------------------------------
# Composite builder
# ----------------------------------------------------------------------

def build_composite(panels: list[Path], rows: int, cols: int) -> Image.Image | None:
    """Tile panels into a rows×cols grid, each panel at CELL_SIZE.

    Pads with white if fewer panels than rows*cols.
    """
    if rows * cols < len(panels):
        return None  # not enough cells
    cell_w, cell_h = CELL_SIZE
    canvas = Image.new("RGB", (cell_w * cols, cell_h * rows), (255, 255, 255))
    for i, panel in enumerate(panels):
        try:
            img = Image.open(panel).convert("RGB").resize(CELL_SIZE, Image.LANCZOS)
        except Exception:
            continue
        r, c = divmod(i, cols)
        canvas.paste(img, (c * cell_w, r * cell_h))
    return canvas


# ----------------------------------------------------------------------
# Locked mapping + April baseline loaders
# ----------------------------------------------------------------------

def load_locked_rows() -> list[dict]:
    rows = []
    with open(LOCKED_MAPPING, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            artifacts_raw = (r.get("mapped_artifacts_or_signals") or "").strip()
            artifacts = [a.strip() for a in artifacts_raw.split(";") if a.strip()]
            # Skip CSV / Table targets — they don't have image artifacts to composite.
            artifact_type = (r.get("artifact_type") or "").strip()
            if artifact_type == "Table":
                continue
            rows.append({
                "target_id": (r.get("target_id") or "").strip(),
                "paper_item": (r.get("paper_item") or "").strip(),
                "artifact_type": artifact_type,
                "dataset_scope": (r.get("dataset_scope") or "").strip(),
                "mapped_artifacts": artifacts,
            })
    return rows


def load_april_baseline() -> dict[str, dict]:
    """Returns {target_id: {layout, score, status}}."""
    out: dict[str, dict] = {}
    if not APRIL_REPORT.exists():
        return out
    with open(APRIL_REPORT, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            tid = (r.get("target_id") or "").strip()
            if not tid:
                continue
            try:
                score = float(r.get("best_score") or 0)
            except ValueError:
                score = 0.0
            out[tid] = {
                "layout": (r.get("best_layout") or "").strip(),
                "score": score,
                "status": (r.get("status") or "").strip(),
            }
    return out


def is_si(row: dict) -> bool:
    pi = row["paper_item"]
    return " S" in pi or pi.startswith("Fig. S") or "FS" in row["target_id"] or "-S-" in row["target_id"]


def list_paper_images(extract_dir: Path) -> list[Path]:
    images = []
    if (extract_dir / "images").exists():
        images.extend(sorted((extract_dir / "images").glob("*.png")))
        images.extend(sorted((extract_dir / "images").glob("*.jpg")))
    if not images and (extract_dir / "pages").exists():
        images.extend(sorted((extract_dir / "pages").glob("*.png")))
    return images


# ----------------------------------------------------------------------
# Per-target scoring
# ----------------------------------------------------------------------

def evaluate_target(row: dict) -> dict:
    """Score a target by trying every layout × paper image combo. Returns one row."""
    artifact_base = PAPER_ARTIFACTS.get(row["dataset_scope"])
    if artifact_base is None:
        return {
            "target_id": row["target_id"], "paper_item": row["paper_item"],
            "n_mapped": len(row["mapped_artifacts"]), "n_present": 0,
            "best_layout": "", "best_score": 0.0, "best_paper_image": "",
            "composite_path": "", "status": "ERROR_NO_SCOPE",
        }
    produced = [artifact_base / name for name in row["mapped_artifacts"]]
    present = [p for p in produced if p.exists()]
    if not present:
        return {
            "target_id": row["target_id"], "paper_item": row["paper_item"],
            "n_mapped": len(row["mapped_artifacts"]), "n_present": 0,
            "best_layout": "", "best_score": float("nan"), "best_paper_image": "",
            "composite_path": "", "status": "MISSING_ALL_PANELS",
        }

    paper_dir = PDF_EXTRACT_SUBDIRS.get((row["dataset_scope"], is_si(row)))
    if paper_dir is None:
        return {
            "target_id": row["target_id"], "paper_item": row["paper_item"],
            "n_mapped": len(row["mapped_artifacts"]), "n_present": len(present),
            "best_layout": "", "best_score": float("nan"), "best_paper_image": "",
            "composite_path": "", "status": "ERROR_NO_PDF_EXTRACT",
        }
    paper_images = list_paper_images(paper_dir)
    if not paper_images:
        return {
            "target_id": row["target_id"], "paper_item": row["paper_item"],
            "n_mapped": len(row["mapped_artifacts"]), "n_present": len(present),
            "best_layout": "", "best_score": float("nan"), "best_paper_image": "",
            "composite_path": "", "status": "ERROR_NO_PAPER_IMAGES",
        }

    # Try every layout × paper image combo; track the best.
    best_score = -2.0
    best_layout = (1, 1)
    best_paper = None
    best_composite_arr = None
    best_composite_pil = None

    for rows_, cols_ in candidate_layouts(len(present)):
        comp = build_composite(present, rows_, cols_)
        if comp is None:
            continue
        comp_arr = np.asarray(comp.convert("L").resize(COMPARE_SIZE, Image.LANCZOS), dtype=float)
        for paper_img in paper_images:
            paper_arr = load_grayscale(paper_img)
            if paper_arr is None:
                continue
            s = ncc(comp_arr, paper_arr)
            if s > best_score:
                best_score = s
                best_layout = (rows_, cols_)
                best_paper = paper_img
                best_composite_pil = comp

    # Save the best composite for visual inspection
    composite_path = ""
    if best_composite_pil is not None:
        composite_path = OUT_COMPOSITES / f"{row['target_id']}_{best_layout[0]}x{best_layout[1]}.png"
        try:
            best_composite_pil.save(composite_path)
        except Exception:
            composite_path = ""

    status = "PASS" if best_score >= DEFAULT_THRESHOLD else "FAIL"
    return {
        "target_id": row["target_id"], "paper_item": row["paper_item"],
        "n_mapped": len(row["mapped_artifacts"]), "n_present": len(present),
        "best_layout": f"{best_layout[0]}x{best_layout[1]}" if best_composite_pil else "",
        "best_score": best_score if best_composite_pil else float("nan"),
        "best_paper_image": str(best_paper.relative_to(REPO_ROOT)) if best_paper else "",
        "composite_path": str(composite_path.relative_to(REPO_ROOT)) if composite_path else "",
        "status": status,
    }


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    print(f"=== Composite paper-comparison @ {TIMESTAMP} ===")
    print(f"  output: {OUT_RUN_DIR.relative_to(REPO_ROOT)}")
    print()

    rows = load_locked_rows()
    april = load_april_baseline()
    print(f"  locked mapping (non-Table rows): {len(rows)}")
    print(f"  april baseline: {len(april)} targets")
    print()

    results = []
    t0 = time.time()
    for i, row in enumerate(rows, 1):
        if row["dataset_scope"] not in ("DATA1", "DATA2"):
            continue
        result = evaluate_target(row)
        # Append April baseline for drift comparison
        ap = april.get(row["target_id"], {})
        result["april_layout"] = ap.get("layout", "")
        result["april_score"] = ap.get("score", float("nan"))
        result["april_status"] = ap.get("status", "")
        if ap:
            try:
                result["drift"] = result["best_score"] - ap["score"]
            except (TypeError, ValueError):
                result["drift"] = float("nan")
        else:
            result["drift"] = float("nan")
        results.append(result)
        score_disp = f"{result['best_score']:.3f}" if not np.isnan(result['best_score']) else "n/a"
        april_disp = f"{ap.get('score', 0):.3f}" if ap else "—"
        marker = "✓" if result["status"] == "PASS" else "✗"
        print(f"  [{i:>2d}/{len(rows)}] {row['target_id']:<10s}  "
              f"present={result['n_present']}/{result['n_mapped']}  "
              f"NCC={score_disp:>6s}  april={april_disp:>6s}  {marker} {result['status']}")

    print()
    print(f"  elapsed: {time.time()-t0:.1f}s")
    print()

    # Write CSV
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "target_id", "paper_item", "n_mapped", "n_present",
            "best_layout", "best_score", "best_paper_image", "composite_path",
            "status", "april_layout", "april_score", "april_status", "drift",
        ])
        w.writeheader()
        for r in results:
            w.writerow(r)
    import shutil
    shutil.copy(OUT_CSV, OUT_LATEST_CSV)
    print(f"  CSV    → {OUT_CSV.relative_to(REPO_ROOT)}")
    print(f"  latest → {OUT_LATEST_CSV.relative_to(REPO_ROOT)}")

    # Write MD
    md = [f"# Composite Paper-Comparison Report — {TIMESTAMP}", ""]
    md.append(f"- Source: `{LOCKED_MAPPING.relative_to(REPO_ROOT)}`")
    md.append(f"- April baseline: `{APRIL_REPORT.relative_to(REPO_ROOT)}` (2026-04-17)")
    md.append(f"- Metric: Normalized Cross Correlation (NCC) on composite vs paper figure")
    md.append(f"- Threshold: {DEFAULT_THRESHOLD} (April 17 baseline)")
    md.append("")
    n_pass = sum(1 for r in results if r["status"] == "PASS")
    n_fail = sum(1 for r in results if r["status"] == "FAIL")
    n_missing = sum(1 for r in results if r["status"] == "MISSING_ALL_PANELS")
    md.append(f"## Summary")
    md.append("")
    md.append(f"- **PASS (NCC ≥ 0.26):**  {n_pass} / {len(results)}")
    md.append(f"- **FAIL (NCC < 0.26):**  {n_fail} / {len(results)}")
    md.append(f"- **MISSING all panels:**  {n_missing} / {len(results)}")
    md.append("")
    md.append(f"## Per-target results")
    md.append("")
    md.append("| target | paper item | present | layout | NCC (now) | NCC (Apr 17) | drift | status |")
    md.append("|---|---|---|---|---|---|---|---|")
    for r in results:
        score_disp = f"{r['best_score']:.3f}" if not np.isnan(r['best_score']) else "—"
        april_disp = f"{r['april_score']:.3f}" if not np.isnan(r.get('april_score', float('nan'))) else "—"
        drift_disp = f"{r['drift']:+.3f}" if not np.isnan(r.get('drift', float('nan'))) else "—"
        marker = "✓" if r["status"] == "PASS" else "✗"
        md.append(f"| `{r['target_id']}` | {r['paper_item']} | {r['n_present']}/{r['n_mapped']} | "
                  f"{r['best_layout']} | {score_disp} | {april_disp} | {drift_disp} | {marker} **{r['status']}** |")
    md.append("")
    md.append(f"## Composites saved at `{OUT_COMPOSITES.relative_to(REPO_ROOT)}/`")
    OUT_MD.write_text("\n".join(md))
    print(f"  MD     → {OUT_MD.relative_to(REPO_ROOT)}")

    # Console summary
    print()
    print("=== Summary ===")
    print(f"  PASS (NCC ≥ {DEFAULT_THRESHOLD}):  {n_pass} / {len(results)}")
    print(f"  FAIL (NCC < {DEFAULT_THRESHOLD}):  {n_fail} / {len(results)}")
    print(f"  MISSING all panels:    {n_missing} / {len(results)}")
    print()
    drifts = [r for r in results if not np.isnan(r.get('drift', float('nan')))]
    if drifts:
        big_drift = [r for r in drifts if abs(r['drift']) >= 0.05]
        print(f"  Targets with |drift| ≥ 0.05 vs April 17: {len(big_drift)}")
        for r in sorted(big_drift, key=lambda x: x['drift']):
            sign = "+" if r['drift'] >= 0 else ""
            print(f"    {r['target_id']:<10s}  {sign}{r['drift']:.3f}  ({r['april_score']:.3f} → {r['best_score']:.3f})")


if __name__ == "__main__":
    main()
