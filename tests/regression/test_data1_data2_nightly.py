"""Nightly slow regression scaffold for DATA1/DATA2 paper reproduction."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts/reproduce/reproduce_data1_data2.py"
REFERENCE_CSV = REPO_ROOT / "docs/validation/paper_reference_values_template.csv"


def _has_reference_values(path: Path) -> bool:
    """Return True if at least one paper_value is populated."""
    if not path.exists():
        return False
    df = pd.read_csv(path)
    if "paper_value" not in df.columns:
        return False
    return bool(df["paper_value"].notna().any())


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

    df = pd.read_csv(side_by_side)
    failing = df[df["status"] == "FAIL"]
    assert failing.empty, (
        "Tolerance failures detected in nightly reproduction run.\n"
        f"{failing[['dataset_id', 'metric', 'relative_error', 'status']].to_string(index=False)}"
    )
