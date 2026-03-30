"""Lightweight pytest validation checks against committed paper artifacts."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
import pytest
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from UnifiedFramework.DATA3.scripts.validation.nightly.nightly_validation_gate import (  # noqa: E402
    evaluate,
    read_exceptions,
    read_status_rows,
)

VALIDATION_ROOT = REPO_ROOT / "UnifiedFramework/DATA3/docs/validation"
STATUS_CSV = VALIDATION_ROOT / "target_validation_status_consolidated.csv"
EXCEPTIONS_CSV = VALIDATION_ROOT / "nightly/config/known_nonpass.csv"
SOURCE_MAP_CSV = VALIDATION_ROOT / "target_notebook_source_map.csv"
PAGE_INDEX_CSV = VALIDATION_ROOT / "target_pdf_page_index.csv"
PAPER_EXTRACT_ROOT = (
    REPO_ROOT
    / "UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract"
)
ARTIFACT_SEARCH_ROOTS = (
    REPO_ROOT,
    REPO_ROOT / "UnifiedFramework/DATA3/figures",
    REPO_ROOT / "UnifiedFramework/DATA3/results/reproduction/20260221-004009/figures",
    REPO_ROOT / "UnifiedFramework/DATA3/results/reproduction/20260221-data1-notebook/figures/data1_notebook_regen",
)
ARTIFACT_PATTERN_RE = re.compile(r"[\w./*-]+\.(?:png|jpg|jpeg|csv)")


def _figure_target_rows() -> list[dict[str, str]]:
    source_map = pd.read_csv(SOURCE_MAP_CSV)
    page_index = pd.read_csv(PAGE_INDEX_CSV)
    merged = source_map.merge(
        page_index[["target_id", "pdf_key", "matched_pages", "match_status", "evidence_run_id"]],
        on="target_id",
        how="inner",
    )
    figure_rows = merged[merged["artifact_type"].isin(["Figure", "Figure set"])].copy()
    figure_rows = figure_rows[figure_rows["match_status"] == "MATCHED"]
    return figure_rows.to_dict("records")


FIGURE_TARGET_ROWS = _figure_target_rows()


def _parse_pages(raw_pages: str) -> list[int]:
    return [int(token.strip()) for token in str(raw_pages).split(";") if token.strip()]


def _extract_artifact_patterns(raw_value: str) -> list[str]:
    seen: set[str] = set()
    patterns: list[str] = []
    for match in ARTIFACT_PATTERN_RE.findall(str(raw_value)):
        if match not in seen:
            patterns.append(match)
            seen.add(match)
    return patterns


def _resolve_pattern(pattern: str) -> list[Path]:
    raw = Path(pattern)
    matches: list[Path] = []

    if raw.is_absolute():
        return [raw] if raw.exists() else []

    repo_candidate = REPO_ROOT / pattern
    if "*" in pattern:
        matches.extend(path for path in REPO_ROOT.glob(pattern) if path.exists())
        if not matches:
            for root in ARTIFACT_SEARCH_ROOTS:
                matches.extend(path for path in root.glob(raw.name) if path.exists())
    else:
        if repo_candidate.exists():
            matches.append(repo_candidate)
        for root in ARTIFACT_SEARCH_ROOTS:
            candidate = root / raw.name
            if candidate.exists():
                matches.append(candidate)
        if not matches:
            matches.extend(path for path in REPO_ROOT.rglob(raw.name) if path.exists())

    deduped: list[Path] = []
    seen: set[Path] = set()
    for match in matches:
        resolved = match.resolve()
        if resolved not in seen:
            deduped.append(match)
            seen.add(resolved)
    return deduped


def _assert_image_is_readable(path: Path) -> None:
    with Image.open(path) as image:
        image.verify()


@pytest.mark.nightly
@pytest.mark.regression
def test_nightly_status_gate_matches_allowlist() -> None:
    """Nightly gate should fail only on genuinely new regressions."""
    failures, warnings, resolved, counts = evaluate(
        read_status_rows(STATUS_CSV),
        read_exceptions(EXCEPTIONS_CSV),
    )

    assert not failures, (
        "Unexpected pytest validation failures found.\n"
        + "\n".join(f"- {failure}" for failure in failures)
    )
    assert counts, "Expected consolidated validation status rows to be present."
    assert warnings is not None
    assert resolved is not None


@pytest.mark.nightly
@pytest.mark.regression
@pytest.mark.parametrize(
    "row",
    FIGURE_TARGET_ROWS,
    ids=[record["target_id"] for record in FIGURE_TARGET_ROWS],
)
def test_nightly_figure_targets_have_paper_pages_and_artifacts(row: dict[str, str]) -> None:
    """Every mapped figure target should point to readable paper pages and local artifacts."""
    pdf_key = str(row["pdf_key"]).strip()
    assert pdf_key, f"Missing pdf_key for {row['target_id']}"

    pages = _parse_pages(str(row["matched_pages"]))
    assert pages, f"No matched paper pages recorded for {row['target_id']}"
    for page in pages:
        page_path = PAPER_EXTRACT_ROOT / pdf_key / "pages" / f"page-{page:02d}.png"
        assert page_path.exists(), f"Missing extracted paper page for {row['target_id']}: {page_path}"
        _assert_image_is_readable(page_path)

    artifact_patterns = _extract_artifact_patterns(str(row["mapped_artifacts_or_signals"]))
    assert artifact_patterns, f"No artifact filenames parsed for {row['target_id']}"
    for pattern in artifact_patterns:
        matches = _resolve_pattern(pattern)
        assert matches, f"Could not resolve mapped artifact '{pattern}' for {row['target_id']}"
        for match in matches:
            if match.suffix.lower() in {".png", ".jpg", ".jpeg"}:
                _assert_image_is_readable(match)
