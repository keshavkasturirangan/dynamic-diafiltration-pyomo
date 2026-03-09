#!/usr/bin/env python3
"""Compare digitized paper panel curves against unified-generated curves.

Input files:
- manifest CSV: panel/file/column mapping per target.
- thresholds CSV: per-target tolerances (+ DEFAULT fallback).

Status semantics per row:
- PASS: metrics within thresholds.
- FAIL: metrics outside thresholds.
- MISSING_VALUE: missing files/columns, insufficient overlap, or too few points.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Threshold:
    nrmse_max: float
    mae_max: float
    max_abs_err_max: float
    min_points: int


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare digitized paper curves to unified curves.")
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument(
        "--manifest",
        type=Path,
        default=Path("docs/validation/nightly/digitized_baselines/manifest.csv"),
    )
    p.add_argument(
        "--thresholds",
        type=Path,
        default=Path("docs/validation/nightly/digitized_baselines/thresholds.csv"),
    )
    p.add_argument(
        "--out-csv",
        type=Path,
        default=Path("docs/validation/nightly/digitized_baselines/comparison_latest.csv"),
    )
    p.add_argument(
        "--fail-on-fail",
        action="store_true",
        help="Exit non-zero if any row status is FAIL.",
    )
    p.add_argument(
        "--fail-on-missing",
        action="store_true",
        help="Exit non-zero if any row status is MISSING_VALUE.",
    )
    return p.parse_args()


def _resolve(repo_root: Path, p: str) -> Path:
    path = Path(str(p).strip())
    return path if path.is_absolute() else (repo_root / path)


def _to_bool(v: object) -> bool:
    return str(v).strip().lower() in {"1", "true", "yes", "y"}


def _safe_float(v: object) -> float | None:
    s = str(v).strip()
    if s == "" or s.lower() == "nan":
        return None
    return float(s)


def _load_thresholds(path: Path) -> Dict[str, Threshold]:
    df = pd.read_csv(path)
    req = {"target_id", "nrmse_max", "mae_max", "max_abs_err_max", "min_points"}
    missing = req - set(df.columns)
    if missing:
        raise ValueError(f"Missing threshold columns: {sorted(missing)}")

    out: Dict[str, Threshold] = {}
    for _, row in df.iterrows():
        tid = str(row["target_id"]).strip()
        if not tid:
            continue
        out[tid] = Threshold(
            nrmse_max=float(row["nrmse_max"]),
            mae_max=float(row["mae_max"]),
            max_abs_err_max=float(row["max_abs_err_max"]),
            min_points=int(row["min_points"]),
        )
    return out


def _load_xy(df: pd.DataFrame, x_col: str, y_col: str) -> Tuple[np.ndarray, np.ndarray]:
    x = pd.to_numeric(df[x_col], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[y_col], errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size == 0:
        return x, y
    order = np.argsort(x)
    return x[order], y[order]


def _metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[float, float, float]:
    err = y_pred - y_true
    mae = float(np.mean(np.abs(err)))
    rmse = float(np.sqrt(np.mean(err**2)))
    rng = float(np.max(y_true) - np.min(y_true)) if y_true.size else 0.0
    nrmse = rmse / rng if rng > 0 else rmse
    max_abs = float(np.max(np.abs(err))) if err.size else np.nan
    return nrmse, mae, max_abs


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    manifest_path = _resolve(repo_root, str(args.manifest))
    thresholds_path = _resolve(repo_root, str(args.thresholds))
    out_csv = _resolve(repo_root, str(args.out_csv))

    manifest = pd.read_csv(manifest_path)
    thresholds = _load_thresholds(thresholds_path)
    default_thr = thresholds.get("DEFAULT")
    if default_thr is None:
        raise ValueError("thresholds.csv must include a DEFAULT row")

    required_manifest = {
        "target_id",
        "panel_id",
        "active",
        "paper_csv",
        "unified_csv",
        "x_col_paper",
        "y_col_paper",
        "x_col_unified",
        "y_col_unified",
        "x_min",
        "x_max",
    }
    missing_cols = required_manifest - set(manifest.columns)
    if missing_cols:
        raise ValueError(f"Missing manifest columns: {sorted(missing_cols)}")

    rows_out: List[Dict[str, object]] = []

    for _, row in manifest.iterrows():
        if not _to_bool(row.get("active", False)):
            continue

        target_id = str(row["target_id"]).strip()
        panel_id = str(row["panel_id"]).strip()
        thr = thresholds.get(target_id, default_thr)

        paper_path = _resolve(repo_root, str(row["paper_csv"]))
        unified_path = _resolve(repo_root, str(row["unified_csv"]))
        xcp = str(row["x_col_paper"]).strip()
        ycp = str(row["y_col_paper"]).strip()
        xcu = str(row["x_col_unified"]).strip()
        ycu = str(row["y_col_unified"]).strip()
        x_min = _safe_float(row.get("x_min", ""))
        x_max = _safe_float(row.get("x_max", ""))

        out = {
            "target_id": target_id,
            "panel_id": panel_id,
            "paper_csv": str(paper_path),
            "unified_csv": str(unified_path),
            "n_points": 0,
            "nrmse": np.nan,
            "mae": np.nan,
            "max_abs_err": np.nan,
            "nrmse_max": thr.nrmse_max,
            "mae_max": thr.mae_max,
            "max_abs_err_max": thr.max_abs_err_max,
            "min_points": thr.min_points,
            "status": "MISSING_VALUE",
            "notes": "",
        }

        if not paper_path.exists() or not unified_path.exists():
            out["notes"] = "missing input file(s)"
            rows_out.append(out)
            continue

        try:
            df_p = pd.read_csv(paper_path)
            df_u = pd.read_csv(unified_path)
        except Exception as exc:
            out["notes"] = f"csv read error: {exc}"
            rows_out.append(out)
            continue

        if any(c not in df_p.columns for c in (xcp, ycp)) or any(c not in df_u.columns for c in (xcu, ycu)):
            out["notes"] = "missing required x/y column(s)"
            rows_out.append(out)
            continue

        xp, yp = _load_xy(df_p, xcp, ycp)
        xu, yu = _load_xy(df_u, xcu, ycu)
        if xp.size == 0 or xu.size < 2:
            out["notes"] = "insufficient numeric points"
            rows_out.append(out)
            continue

        lo = max(np.min(xp), np.min(xu))
        hi = min(np.max(xp), np.max(xu))
        if x_min is not None:
            lo = max(lo, x_min)
        if x_max is not None:
            hi = min(hi, x_max)

        mask = (xp >= lo) & (xp <= hi)
        x_eval = xp[mask]
        y_true = yp[mask]
        if x_eval.size < thr.min_points:
            out["n_points"] = int(x_eval.size)
            out["notes"] = "insufficient overlap points in x-range"
            rows_out.append(out)
            continue

        y_pred = np.interp(x_eval, xu, yu)
        nrmse, mae, max_abs = _metrics(y_true=y_true, y_pred=y_pred)

        out["n_points"] = int(x_eval.size)
        out["nrmse"] = nrmse
        out["mae"] = mae
        out["max_abs_err"] = max_abs

        if nrmse <= thr.nrmse_max and mae <= thr.mae_max and max_abs <= thr.max_abs_err_max:
            out["status"] = "PASS"
            out["notes"] = "within thresholds"
        else:
            out["status"] = "FAIL"
            out["notes"] = "threshold exceedance"

        rows_out.append(out)

    out_df = pd.DataFrame(rows_out)
    if out_df.empty:
        out_df = pd.DataFrame(
            columns=[
                "target_id",
                "panel_id",
                "paper_csv",
                "unified_csv",
                "n_points",
                "nrmse",
                "mae",
                "max_abs_err",
                "nrmse_max",
                "mae_max",
                "max_abs_err_max",
                "min_points",
                "status",
                "notes",
            ]
        )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_csv, index=False)

    counts = out_df["status"].value_counts(dropna=False).to_dict() if not out_df.empty else {}
    print(f"[digitized-compare] wrote: {out_csv}")
    print(f"[digitized-compare] rows: {len(out_df)}")
    print(f"[digitized-compare] status_counts: {counts}")

    has_fail = bool((out_df["status"] == "FAIL").any()) if not out_df.empty else False
    has_missing = bool((out_df["status"] == "MISSING_VALUE").any()) if not out_df.empty else False

    if args.fail_on_fail and has_fail:
        raise SystemExit(1)
    if args.fail_on_missing and has_missing:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
