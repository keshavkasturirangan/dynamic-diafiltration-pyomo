"""Paper-comparison layer — shared utilities.

Loads the two locked authorities (target_notebook_source_map.csv and
simulation_validation_manifest.csv) and provides per-row comparison helpers.

Design contract:
  - Read-only on the external infrastructure (no edits to
    UnifiedFramework/DATA3/docs/validation/ or to refactored_codes_v1/).
  - Tolerances loaded from pytest_refactored_codes_v1/baselines/paper_comparison_tolerances.csv.
  - Drift report emitted to _runs/run-<id>/paper_comparison_report.csv per session.

For the architectural rationale see docs/PAPER_COMPARISON_LAYER.md.
For the original spec see docs/PAPER_COMPARISON_PROMPT.md.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ----------------------------------------------------------------------
# Path anchors — resolved relative to this file so the helper works no
# matter where pytest is invoked from.
# ----------------------------------------------------------------------

HERE = Path(__file__).resolve().parent             # …/pytest_refactored_codes_v1/helpers/
PYTEST_ROOT = HERE.parent                          # …/pytest_refactored_codes_v1/
REPO_ROOT = PYTEST_ROOT.parent                     # …/ (repo top)
UF_VALIDATION = REPO_ROOT / "UnifiedFramework" / "DATA3" / "docs" / "validation"

LOCKED_MAPPING_CSV = UF_VALIDATION / "target_notebook_source_map.csv"
NUMERIC_MANIFEST_CSV = UF_VALIDATION / "simulation_validation_data_files" / "simulation_validation_manifest.csv"
SIM_VAL_DATA_ROOT = UF_VALIDATION / "simulation_validation_data_files"

PAPER_ARTIFACTS_DATA1 = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "data1" / "notebook_figures"
PAPER_ARTIFACTS_DATA2 = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "data2" / "notebook_figures"

TOLERANCE_OVERRIDES_CSV = PYTEST_ROOT / "baselines" / "paper_comparison_tolerances.csv"


# ----------------------------------------------------------------------
# Default tolerances per artifact_kind (overridden per target_id in tolerances.csv)
# ----------------------------------------------------------------------

DEFAULT_TOLERANCES = {
    "timeseries":    {"metric": "nrmse", "tolerance": 0.01},
    "scatter":       {"metric": "mae",   "tolerance": 0.01},
    "contour_grid":  {"metric": "maxrel", "tolerance": 0.02},
    "table":         {"metric": "relmax", "tolerance": 0.01},
}


# ----------------------------------------------------------------------
# Data classes — one row of each manifest
# ----------------------------------------------------------------------

@dataclass
class LockedMappingRow:
    """One row of target_notebook_source_map.csv."""
    target_id: str                       # e.g. "D1-M-F2"
    paper_item: str                      # e.g. "Fig. 2 (all four plots)"
    artifact_type: str                   # "Figure" | "Figure set" | "Table"
    dataset_scope: str                   # "DATA1" | "DATA2"
    mapped_artifacts: list[str]          # parsed semicolon list of filenames
    mapping_confidence: str              # "HIGH" | "MEDIUM" | "LOW"
    notes: str

    def expected_paths(self) -> list[Path]:
        """Resolve each mapped filename against the right paper_artifacts subdir."""
        base = PAPER_ARTIFACTS_DATA1 if self.dataset_scope == "DATA1" else PAPER_ARTIFACTS_DATA2
        return [base / name for name in self.mapped_artifacts]


@dataclass
class NumericManifestRow:
    """One row of simulation_validation_manifest.csv."""
    paper: str                           # "DATA1" | "DATA2"
    figure_id: str                       # e.g. "Fig. 2"
    panel_id: str                        # e.g. "A"
    artifact_kind: str                   # "timeseries" | "scatter" | "contour_grid" | "table"
    source_label: str
    relative_csv: str                    # e.g. "data1_main/fig2/data1_main_fig2_a.csv"
    notes: str

    @property
    def case_id(self) -> str:
        """Stable identifier for pytest.param ids."""
        fig = self.figure_id.replace(" ", "_").replace(".", "").replace(",", "")
        pan = self.panel_id.replace(" ", "_").replace(",", "_").replace("/", "_")
        return f"{self.paper}__{fig}__{pan}"

    def baseline_path(self) -> Path:
        return SIM_VAL_DATA_ROOT / self.relative_csv


# ----------------------------------------------------------------------
# Loaders — both authorities read once per pytest session
# ----------------------------------------------------------------------

def load_locked_mapping() -> list[LockedMappingRow]:
    """Parse target_notebook_source_map.csv into a list of dataclasses."""
    rows: list[LockedMappingRow] = []
    if not LOCKED_MAPPING_CSV.exists():
        return rows
    with open(LOCKED_MAPPING_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            artifacts_raw = (r.get("mapped_artifacts_or_signals") or "").strip()
            artifacts = [a.strip() for a in artifacts_raw.split(";") if a.strip()]
            rows.append(LockedMappingRow(
                target_id=(r.get("target_id") or "").strip(),
                paper_item=(r.get("paper_item") or "").strip(),
                artifact_type=(r.get("artifact_type") or "").strip(),
                dataset_scope=(r.get("dataset_scope") or "").strip(),
                mapped_artifacts=artifacts,
                mapping_confidence=(r.get("mapping_confidence") or "").strip(),
                notes=(r.get("notes") or "").strip(),
            ))
    return rows


def load_numeric_manifest() -> list[NumericManifestRow]:
    """Parse simulation_validation_manifest.csv into a list of dataclasses.

    Filters out rows whose relative_csv points at a non-CSV file (e.g. _ERROR.txt
    stubs from legacy extraction failures). Such rows are tracked as stale but
    not treated as real comparison targets.
    """
    rows: list[NumericManifestRow] = []
    if not NUMERIC_MANIFEST_CSV.exists():
        return rows
    with open(NUMERIC_MANIFEST_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rel = (r.get("relative_csv") or "").strip()
            # Skip stale rows pointing at non-CSV stubs.
            if not rel.lower().endswith(".csv"):
                continue
            rows.append(NumericManifestRow(
                paper=(r.get("paper") or "").strip(),
                figure_id=(r.get("figure_id") or "").strip(),
                panel_id=(r.get("panel_id") or "").strip(),
                artifact_kind=(r.get("artifact_kind") or "").strip(),
                source_label=(r.get("source_label") or "").strip(),
                relative_csv=rel,
                notes=(r.get("notes") or "").strip(),
            ))
    return rows


def stale_manifest_rows() -> list[str]:
    """List manifest rows whose relative_csv points at non-CSV stubs.

    Used by smoke tests to surface stale extraction failures as warnings
    without failing the suite.
    """
    out: list[str] = []
    if not NUMERIC_MANIFEST_CSV.exists():
        return out
    with open(NUMERIC_MANIFEST_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rel = (r.get("relative_csv") or "").strip()
            if rel and not rel.lower().endswith(".csv"):
                out.append(f"{r.get('paper')}/{r.get('figure_id')}/{r.get('panel_id')} -> {rel}")
    return out


def load_tolerance_overrides() -> dict:
    """Parse paper_comparison_tolerances.csv into {target_id: {metric, tolerance, reason}}."""
    out: dict = {}
    if not TOLERANCE_OVERRIDES_CSV.exists():
        return out
    with open(TOLERANCE_OVERRIDES_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            target_id = (r.get("target_id") or "").strip()
            if not target_id:
                continue
            try:
                tol = float(r.get("tolerance") or 0.01)
            except ValueError:
                tol = 0.01
            out[target_id] = {
                "metric": (r.get("metric") or "").strip().lower() or None,
                "tolerance": tol,
                "reason": (r.get("reason") or "").strip(),
            }
    return out


# ----------------------------------------------------------------------
# Per-row resolution helpers
# ----------------------------------------------------------------------

def resolved_tolerance(row: NumericManifestRow, overrides: dict) -> tuple[str, float]:
    """Return (metric_name, tolerance) for one numeric manifest row.

    Override > default-per-kind. Falls back to ('nrmse', 0.01) if nothing matches.
    """
    case = f"{row.paper}.{row.figure_id}.{row.panel_id}"
    if case in overrides:
        ov = overrides[case]
        if ov.get("metric") and ov.get("tolerance") is not None:
            return ov["metric"], ov["tolerance"]
    defaults = DEFAULT_TOLERANCES.get(row.artifact_kind)
    if defaults:
        return defaults["metric"], defaults["tolerance"]
    return "nrmse", 0.01


# ----------------------------------------------------------------------
# Metric helpers — used by the nightly (solver-driven) sub-layer
# ----------------------------------------------------------------------

def nrmse(predicted, baseline) -> float:
    """Normalized RMSE: RMSE / range(baseline). 0 = perfect, 1 = order-of-magnitude wrong."""
    import numpy as np
    a = np.asarray(predicted, dtype=float)
    b = np.asarray(baseline, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"nrmse: shape mismatch — predicted {a.shape}, baseline {b.shape}")
    rmse_val = float(np.sqrt(np.mean((a - b) ** 2)))
    rng = float(np.ptp(b))
    if rng <= 0:
        # All-baseline-constant → fall back to absolute RMSE.
        return rmse_val
    return rmse_val / rng


def mae(predicted, baseline) -> float:
    """Mean Absolute Error, scaled to mean(|baseline|) so it's dimensionless."""
    import numpy as np
    a = np.asarray(predicted, dtype=float)
    b = np.asarray(baseline, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"mae: shape mismatch — predicted {a.shape}, baseline {b.shape}")
    scale = float(np.mean(np.abs(b))) or 1.0
    return float(np.mean(np.abs(a - b))) / scale


def maxrel(predicted, baseline, *, eps: float = 1e-9) -> float:
    """Max absolute relative diff over the grid. eps floors small-value division."""
    import numpy as np
    a = np.asarray(predicted, dtype=float)
    b = np.asarray(baseline, dtype=float)
    denom = np.maximum(np.abs(b), eps)
    return float(np.max(np.abs(a - b) / denom))


METRIC_FUNCTIONS = {
    "nrmse": nrmse,
    "mae": mae,
    "maxrel": maxrel,
    "relmax": maxrel,  # alias
}


# ----------------------------------------------------------------------
# Baseline CSV loader — small, schema-tolerant
# ----------------------------------------------------------------------

def load_baseline_csv(row: NumericManifestRow):
    """Load the baseline CSV. Returns the pandas DataFrame or None if missing."""
    p = row.baseline_path()
    if not p.exists():
        return None
    import pandas as pd
    return pd.read_csv(p)


# ----------------------------------------------------------------------
# Drift report writer — appended to _runs/run-*/paper_comparison_report.csv
# ----------------------------------------------------------------------

@dataclass
class DriftRow:
    paper: str
    figure_id: str
    panel_id: str
    artifact_kind: str
    target_id: str
    baseline_csv: str
    produced_csv: str
    metric: str
    value: Optional[float]
    tolerance: float
    verdict: str               # PASS | FAIL | SKIP | ERROR
    notes: str = ""


def write_drift_report(report_path: Path, rows: list[DriftRow]) -> None:
    """Write the per-target drift table. Overwrites if exists."""
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "paper", "figure_id", "panel_id", "artifact_kind", "target_id",
            "baseline_csv", "produced_csv", "metric", "value", "tolerance",
            "verdict", "notes",
        ])
        for r in rows:
            w.writerow([
                r.paper, r.figure_id, r.panel_id, r.artifact_kind, r.target_id,
                r.baseline_csv, r.produced_csv, r.metric,
                "" if r.value is None else f"{r.value:.6g}", f"{r.tolerance:.4g}",
                r.verdict, r.notes,
            ])
