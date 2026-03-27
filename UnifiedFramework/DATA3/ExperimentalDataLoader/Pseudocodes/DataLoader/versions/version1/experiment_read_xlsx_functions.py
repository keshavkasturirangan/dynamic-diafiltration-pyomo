#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 11:25:45 2025

@authors: Keshav Kasturi Rangan, Alex Dowling

experiment_read_xlsx_functions.py

One-line purpose:
    Read diafiltration experiments from Excel workbooks where EACH SHEET is ONE experiment.

This module:
    - Reads raw Excel grids to capture header/notes metadata.
    - Detects the measurement-table header row (heuristic).
    - Detects time/mass columns plus optional signals (retentate/permeate) and vial marker.
    - Builds vial ranges + VialData records.
    - Returns file-type agnostic ExperimentData.

Note:
    pandas.read_excel requires an Excel engine (typically openpyxl) for .xlsx.
"""

from __future__ import annotations  # Enable forward references.

from typing import Dict, List, Optional, Tuple  # Explicit types (no Any).

import os  # Basename extraction.
import re  # Pattern matching.

import numpy as np  # Arrays.
import pandas as pd  # Excel parsing.

from experiment_data_types import ExperimentData, VialData, MetaValue, ArrayI  # Shared containers/types.


def list_sheets(path: str) -> List[str]:
    """
    List sheet names in an Excel workbook.

    Arguments:
        path: str, path to an Excel workbook (.xlsx/.xlsm)

    Returns:
        list[str], sheet names in the workbook
    """
    xl = pd.ExcelFile(path)  # Load workbook metadata.
    return list(xl.sheet_names)  # Return sheet names.


def _find_header_row(df_raw: pd.DataFrame, *, max_scan_rows: int = 60) -> int:
    """
    Heuristic: find the measurement-table header row in a sheet read with header=None.

    Arguments:
        df_raw: pd.DataFrame, raw grid
        max_scan_rows: int, how many top rows to scan

    Returns:
        int, header row index (0-based)
    """
    tokens = ["time", "mass", "cond", "conduct", "vial", "swap", "sec", "g"]
    best_row = 0
    best_score = -1

    for r in range(min(max_scan_rows, len(df_raw))):
        row_vals = [str(x).strip().lower() for x in df_raw.iloc[r].tolist()]
        score = sum(any(tok in cell for tok in tokens) for cell in row_vals)
        if score > best_score:
            best_score = score
            best_row = r

    return best_row


def _parse_key_value_metadata(
    df_raw: pd.DataFrame, *, header_row: int
) -> Tuple[Dict[str, MetaValue], Dict[str, MetaValue]]:
    """
    Parse rows above the table into (metadata, config) dicts (best effort).

    Arguments:
        df_raw: pd.DataFrame, raw grid
        header_row: int, row where measurement table header begins

    Returns:
        (metadata, config): tuple of typed dicts
    """
    metadata: Dict[str, MetaValue] = {}
    config: Dict[str, MetaValue] = {}

    patterns: Dict[str, re.Pattern[str]] = {
        "Temp": re.compile(r"\btemp\b|temperature", re.IGNORECASE),
        "delP": re.compile(r"\bdelp\b|delta\s*p|pressure|\bdp\b", re.IGNORECASE),
        "Am": re.compile(r"\bam\b|membrane\s*area|area", re.IGNORECASE),
        "rho": re.compile(r"\brho\b|density", re.IGNORECASE),
        "M_F0": re.compile(r"m[_\s-]*f0|initial\s*feed\s*mass", re.IGNORECASE),
        "C_F0": re.compile(r"c[_\s-]*f0|initial\s*feed\s*conc", re.IGNORECASE),
        "C_D": re.compile(r"c[_\s-]*d|dialysate\s*conc", re.IGNORECASE),
        "M_O": re.compile(r"m[_\s-]*o|overflow\s*mass", re.IGNORECASE),
        "namec": re.compile(r"namec|components?|solute\s*name", re.IGNORECASE),
        "nc": re.compile(r"\bnc\b|num\s*components", re.IGNORECASE),
        "ni": re.compile(r"\bni\b|num\s*species", re.IGNORECASE),
    }

    for r in range(header_row):
        row = df_raw.iloc[r].tolist()
        if len(row) < 2:
            continue

        key_cell = next((x for x in row if str(x).strip() not in ("", "nan", "None")), None)
        if key_cell is None:
            continue

        key_str = str(key_cell).strip()
        key_pos = row.index(key_cell)

        val_cell = next((x for x in row[key_pos + 1 :] if str(x).strip() not in ("", "nan", "None")), None)
        if val_cell is None:
            metadata[f"row_{r}"] = [str(x) for x in row if str(x).strip() not in ("", "nan", "None")]
            continue

        assigned = False
        for cfg_key, pat in patterns.items():
            if pat.search(key_str):
                val_str = str(val_cell).strip()
                try:
                    config[cfg_key] = float(val_str)
                except ValueError:
                    if cfg_key == "namec":
                        config[cfg_key] = [s.strip() for s in val_str.split(",") if s.strip()]
                    else:
                        config[cfg_key] = val_str
                assigned = True
                break

        if not assigned:
            metadata[key_str] = str(val_cell).strip()

    return metadata, config


def _detect_columns(df: pd.DataFrame) -> Dict[str, str]:
    """
    Detect canonical columns in a measurement table.

    Arguments:
        df: pd.DataFrame, measurement table

    Returns:
        dict[str,str], canonical -> actual col names.
        Canonical keys: 'time','mass','vial','ret','perm' (some optional).
    """
    cols = [str(c) for c in df.columns]
    cols_l = [c.lower() for c in cols]

    def find_first(substrings: List[str]) -> Optional[str]:
        for s in substrings:
            for orig, low in zip(cols, cols_l):
                if s in low:
                    return orig
        return None

    out: Dict[str, str] = {}
    tcol = find_first(["time", "seconds", "sec"])
    mcol = find_first(["mass", "grams", " g"])
    vcol = find_first(["vial", "swap", "marker"])
    rcol = find_first(["retent", "ret", "feed cond", "retentate cond", "cf", "c_f"])
    pcol = find_first(["permeat", "perm", "dialysate cond", "permeate cond", "cv", "c_v"])

    if tcol is not None:
        out["time"] = tcol
    if mcol is not None:
        out["mass"] = mcol
    if vcol is not None:
        out["vial"] = vcol
    if rcol is not None:
        out["ret"] = rcol
    if pcol is not None:
        out["perm"] = pcol

    return out


def infer_vial_ranges(vial_marker: ArrayI) -> List[Tuple[int, int]]:
    """
    Convert vial marker labels into contiguous index ranges.

    Arguments:
        vial_marker: np.ndarray[int], label per row

    Returns:
        list[(start,end)] inclusive windows
    """
    ranges: List[Tuple[int, int]] = []
    if vial_marker.size == 0:
        return ranges

    start = 0
    for i in range(1, int(vial_marker.size)):
        if int(vial_marker[i]) != int(vial_marker[i - 1]):
            ranges.append((start, i - 1))
            start = i
    ranges.append((start, int(vial_marker.size) - 1))
    return ranges


def build_vial_marker(vial_ranges: List[Tuple[int, int]], n_total: int) -> ArrayI:
    """
    Build a full-length vial marker array from vial ranges.

    Arguments:
        vial_ranges: list[(start,end)] inclusive
        n_total: int, full length

    Returns:
        np.ndarray[int], labels 1..N
    """
    marker = np.zeros(n_total, dtype=int)
    for k, (s, e) in enumerate(vial_ranges, start=1):
        marker[s : e + 1] = k
    return marker


def read_excel_experiment(path: str, sheet_name: str) -> ExperimentData:
    """
    Read ONE Excel sheet (ONE experiment) and return ExperimentData.

    Arguments:
        path: str, path to Excel workbook
        sheet_name: str, experiment sheet name

    Returns:
        ExperimentData, file-type agnostic container
    """
    df_raw = pd.read_excel(path, sheet_name=sheet_name, header=None)
    header_row = _find_header_row(df_raw)
    metadata, config = _parse_key_value_metadata(df_raw, header_row=header_row)

    df = pd.read_excel(path, sheet_name=sheet_name, header=header_row)
    cols = _detect_columns(df)

    if "time" not in cols or "mass" not in cols:
        raise ValueError(f"Could not detect time/mass columns in '{sheet_name}'. Columns={list(df.columns)}")

    time = pd.to_numeric(df[cols["time"]], errors="coerce").to_numpy(dtype=float)
    mass = pd.to_numeric(df[cols["mass"]], errors="coerce").to_numpy(dtype=float)

    ret = None
    if "ret" in cols:
        ret = pd.to_numeric(df[cols["ret"]], errors="coerce").to_numpy(dtype=float)

    perm = None
    if "perm" in cols:
        perm = pd.to_numeric(df[cols["perm"]], errors="coerce").to_numpy(dtype=float)

    if "vial" in cols:
        vial_marker = (
            pd.to_numeric(df[cols["vial"]], errors="coerce")
            .ffill()
            .fillna(1)
            .to_numpy(dtype=int)
        )
    else:
        vial_marker = np.ones_like(time, dtype=int)

    vial_ranges = infer_vial_ranges(vial_marker)

    vials: List[VialData] = []
    for k, (s, e) in enumerate(vial_ranges, start=1):
        vials.append(
            VialData(
                number=k,
                start_idx=s,
                end_idx=e,
                time=time[s : e + 1],
                mass=mass[s : e + 1],
                retentate_signal=ret[s : e + 1] if ret is not None else None,
                permeate_signal=perm[s : e + 1] if perm is not None else None,
                extras={},
            )
        )

    metadata.setdefault("measurement_df", df)

    return ExperimentData(
        path=path,
        source="xlsx",
        run_id=sheet_name,
        filename=os.path.basename(path),
        metadata=metadata,
        config=config,
        time=time,
        mass=mass,
        retentate_signal=ret,
        permeate_signal=perm,
        vial_marker=vial_marker,
        vial_ranges=vial_ranges,
        vials=vials,
        n_vials=len(vials),
    )


def read_excel_workbook(path: str) -> List[ExperimentData]:
    """
    Read ALL sheets in a workbook, assuming each sheet is one experiment.

    Arguments:
        path: str, path to Excel workbook

    Returns:
        list[ExperimentData], one per sheet (experiment)
    """
    sheets = list_sheets(path)
    experiments: List[ExperimentData] = []

    for sheet in sheets:
        if str(sheet).strip().lower() in {"index", "readme", "notes"}:
            continue
        exp = read_excel_experiment(path, sheet_name=sheet)
        exp.metadata.setdefault("sheet_name", sheet)
        experiments.append(exp)

    return experiments
