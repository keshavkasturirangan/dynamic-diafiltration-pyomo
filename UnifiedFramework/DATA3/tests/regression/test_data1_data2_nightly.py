"""Nightly slow regression scaffold for DATA1/DATA2 paper reproduction."""

from __future__ import annotations

import subprocess
import sys
import json
from pathlib import Path

import pandas as pd
import pytest

from UnifiedFramework.DATA3.scripts.validation.nightly.nightly_validation_gate import read_exceptions


REPO_ROOT = Path(__file__).resolve().parents[4]
SCRIPT = REPO_ROOT / "UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py"
REFERENCE_CSV = REPO_ROOT / "UnifiedFramework/DATA3/docs/validation/paper_reference_values_template.csv"
KNOWN_NONPASS_CSV = REPO_ROOT / "UnifiedFramework/DATA3/docs/validation/nightly/config/known_nonpass.csv"


def _has_reference_values(path: Path) -> bool:
    """Return True if at least one paper_value is populated."""
    if not path.exists():
        return False
    df = pd.read_csv(path)
    if "paper_value" not in df.columns:
        return False
    return bool(df["paper_value"].notna().any())


def _allowed_nonpass_target_ids() -> set[str]:
    """Return target IDs explicitly allowlisted for nightly non-pass statuses."""
    return set(read_exceptions(KNOWN_NONPASS_CSV).keys())


@pytest.mark.slow
@pytest.mark.nightly
@pytest.mark.regression
def test_data1_data2_reproduction_nightly(tmp_path: Path) -> None:
    """Run nightly canonical reproduction and fail on tolerance misses."""
    if not _has_reference_values(REFERENCE_CSV):
        pytest.skip("Reference paper values are not populated yet.")

    out_root = tmp_path / "reproduction"
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--repo-root",
        str(REPO_ROOT),
        "--out-root",
        str(out_root),
        "--reference-csv",
        str(REFERENCE_CSV.relative_to(REPO_ROOT)),
        "--solver",
        "ipopt",
        "--nfe",
        "30",
    ]
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)

    run_dirs = sorted((REPO_ROOT / out_root).glob("*"))
    assert run_dirs, "No reproduction output directory created."
    latest = run_dirs[-1]
    side_by_side = latest / "tables" / "paper_vs_unified_side_by_side.csv"
    assert side_by_side.exists(), "Missing side-by-side comparison CSV."
    meta_path = latest / "run_metadata.json"
    assert meta_path.exists(), "Missing reproduction run metadata JSON."
    meta = json.loads(meta_path.read_text())

    stage_a = meta.get("stage_a_artifacts", {})
    stage_b = meta.get("stage_b_artifacts", {})
    data2_artifacts = meta.get("data2_artifacts", {})
    if stage_a:
        for artifact_path in stage_a.get("figures", []):
            assert Path(artifact_path).exists(), f"Missing DATA1 Stage A figure: {artifact_path}"
        for artifact_path in stage_a.get("tables", []):
            assert Path(artifact_path).exists(), f"Missing DATA1 Stage A table: {artifact_path}"
    if stage_b:
        assert Path(stage_b["manifest"]).exists(), f"Missing DATA1 Stage B manifest: {stage_b['manifest']}"
    for artifact_group in data2_artifacts.values():
        for artifact_path in artifact_group.get("tables", []):
            assert Path(artifact_path).exists(), f"Missing DATA2 artifact table: {artifact_path}"
        for artifact_path in artifact_group.get("figures", []):
            assert Path(artifact_path).exists(), f"Missing DATA2 artifact figure: {artifact_path}"

    df = pd.read_csv(side_by_side)
    failing = df[df["status"] == "FAIL"]
    if failing.empty:
        return

    allowed_targets = _allowed_nonpass_target_ids()
    unexpected = failing[~failing["target_id"].astype(str).isin(allowed_targets)]
    assert unexpected.empty, (
        "Unexpected tolerance failures detected in nightly reproduction run.\n"
        f"{unexpected[['target_id', 'dataset_id', 'metric', 'relative_error', 'status']].to_string(index=False)}"
    )
