#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""experiment_load_v2.py

Compatibility wrapper around the class-field-based, file-agnostic loader:
    experiment_dataload_OOP_v14_simpler.py

Why this file exists
--------------------
Older notebooks/scripts in this project expected a *simple* dictionary output
for quick plotting and sanity checks.

The loader returns an `ExperimentalData` object (plus validation issues).
This wrapper:
    1) calls the class-field loader,
    2) converts `ExperimentalData` -> a small dict of common signals.

This wrapper is intentionally lossy: it does *not* expose all metadata or
assay tables. If you need full fidelity, use the loader module directly.

Author
------
Original: kkasturi (2025)
Refactor: class-field loader wrapper (2026)
"""

from __future__ import annotations

from typing import List, Tuple, Optional, Dict, Union
import os
import sys
from pathlib import Path

import numpy as np

# -----------------------------------------------------------------------------
# Robust local import
# -----------------------------------------------------------------------------
# This repo is often run as "scripts in a folder" (not an installed package).
# To make imports work no matter what your current working directory is,
# we add *this file's directory* to sys.path before importing the loader.
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

# File-agnostic loader (class-field based)
# NOTE: this wrapper intentionally depends only on the *public* entrypoint
#       of the loader module (load_experiment) and the data container class.
from experiment_dataload_OOP_v14_simpler import load_experiment as load_experiment_core
from experiment_dataload_OOP_v14_simpler import ExperimentalData


# =============================================================================
# Public API (kept stable for older notebooks)
# =============================================================================

def load_experiment(
    path: str,
    plot: bool = False,
    save_prefix: Optional[str] = None,
    sheet_name: Optional[Union[str, List[str]]] = None,
    list_sheets: bool = False,
) -> Dict[str, object]:
    """Load one experiment from Excel or MATLAB, optionally plot, return dict.

    Inputs
    ------
    path:
        Path to a .xlsx workbook (multi-sheet) or a .mat file (single experiment).

    plot:
        If True, generate quicklook plots (mass and retentate/permeate signals).

    save_prefix:
        If provided, figures are saved as PNG with this prefix.
        When multiple sheets are requested, the sheet name is appended.

    sheet_name:
        Excel only:
            - None -> auto-select the *first* sheet in the workbook.
            - str  -> load exactly that sheet.
            - list[str] -> load multiple sheets and return a multi-sheet dict.

    list_sheets:
        Excel only: if True, print workbook sheet names and return {}.

    Output
    ------
    dict:
        For a single experiment:
            {
              "times": [np.ndarray, ...],
              "masses": [np.ndarray, ...],
              "retentate_cond": [np.ndarray, ...],
              "permeate_cond": [np.ndarray, ...],
              "windows": [(0, n0-1), (0, n1-1), ...],
              "validation": {"ok_fatal": bool, "issues": [(sev,msg), ...]}
            }

        For multiple Excel sheets:
            {
              "sheets": [...],
              "by_sheet": {"Sheet1": <single-experiment dict>, ...}
            }
    """

    # ---------------------------------------------------------------------
    # 0) Basic checks
    # ---------------------------------------------------------------------
    ext = os.path.splitext(path)[1].lower()
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    # ---------------------------------------------------------------------
    # 1) Sheet listing mode (Excel only)
    # ---------------------------------------------------------------------
    if list_sheets:
        if ext not in {".xlsx", ".xls", ".xlsm", ".xlsb"}:
            print("list_sheets=True is only for Excel workbooks.")
            return {}
        # Import pandas lazily so MAT-only usage doesn't require Excel deps.
        import pandas as pd

        xl = pd.ExcelFile(path)
        print("Available sheets in workbook:")
        for sh in xl.sheet_names:
            print("  -", sh)
        return {}

    # ---------------------------------------------------------------------
    # 2) Excel path: possibly many sheets
    # ---------------------------------------------------------------------
    if ext in {".xlsx", ".xls", ".xlsm", ".xlsb"}:
        # Import pandas lazily so users without Excel support can still load MAT.
        import pandas as pd

        xl = pd.ExcelFile(path)

        # If the caller didn't specify, default to the first sheet.
        if sheet_name is None:
            selected = xl.sheet_names[0]
            # Core loader expects selector = sheet name for XLSX inputs.
            exp, validation = load_experiment_core(path, selector=selected)
            out = _exp_to_simple_dict(exp, validation)
            if plot:
                _plot_quicklook(out, title=f"{os.path.basename(path)} :: {selected}", save_prefix=save_prefix)
            return out

        # One specific sheet.
        if isinstance(sheet_name, str):
            if sheet_name not in xl.sheet_names:
                raise ValueError(f"Sheet '{sheet_name}' not found. Available: {xl.sheet_names}")
            exp, validation = load_experiment_core(path, selector=sheet_name)
            out = _exp_to_simple_dict(exp, validation)
            if plot:
                _plot_quicklook(out, title=f"{os.path.basename(path)} :: {sheet_name}", save_prefix=save_prefix)
            return out

        # Many sheets.
        names: List[str] = list(sheet_name)
        missing = [nm for nm in names if nm not in xl.sheet_names]
        if missing:
            raise ValueError(f"Sheet(s) not found: {missing}. Available: {xl.sheet_names}")

        by_sheet: Dict[str, Dict[str, object]] = {}
        for nm in names:
            exp, validation = load_experiment_core(path, selector=nm)
            out = _exp_to_simple_dict(exp, validation)
            by_sheet[nm] = out
            if plot:
                sp = f"{save_prefix}_{nm}" if save_prefix else None
                _plot_quicklook(out, title=f"{os.path.basename(path)} :: {nm}", save_prefix=sp)

        return {"sheets": names, "by_sheet": by_sheet}

    # ---------------------------------------------------------------------
    # 3) MAT path: single experiment
    # ---------------------------------------------------------------------
    if ext == ".mat":
        # MAT inputs contain a single experiment; selector=None -> default struct key "data_stru".
        exp, validation = load_experiment_core(path, selector=None)
        out = _exp_to_simple_dict(exp, validation)
        if plot:
            _plot_quicklook(out, title=os.path.basename(path), save_prefix=save_prefix)
        return out

    raise ValueError(f"Unsupported file type: {ext}")


# =============================================================================
# Private helpers
# =============================================================================

def _exp_to_simple_dict(
    exp: ExperimentalData,
    validation: Tuple[bool, List[Tuple[str, str]]],
) -> Dict[str, object]:
    """Convert `ExperimentalData` into the legacy quicklook dict.

    Inputs
    ------
    exp:
        ExperimentalData with exp.vials populated.

    validation:
        (ok_fatal, issues) tuple from validate_experiment.

    Output
    ------
    dict with keys: times, masses, retentate_cond, permeate_cond, windows, validation

    Notes
    -----
    - Each vial becomes one entry in each list.
    - `windows` are *per-vial* local indices (0..len(vial)-1), since
      this wrapper does not preserve the original workbook row indices.
    """

    ok_fatal, issues = validation

    times: List[np.ndarray] = []
    masses: List[np.ndarray] = []
    retent: List[np.ndarray] = []
    permea: List[np.ndarray] = []
    windows: List[Tuple[int, int]] = []

    for v in exp.vials:
        # Time is required (validator flags if missing).
        t = np.asarray(v.time_s) if v.time_s is not None else np.array([], dtype=float)
        times.append(t)

        # Mass is optional. If missing, keep a NaN array aligned to time.
        if v.mass_g is None:
            masses.append(np.full_like(t, np.nan, dtype=float))
        else:
            masses.append(np.asarray(v.mass_g, dtype=float))

        # Retentate "signal" is conductivity for XLSX and often conductivity-derived
        # for MAT (we keep it honest and do NOT convert here).
        if v.retentate_signal is None:
            retent.append(np.full_like(t, np.nan, dtype=float))
        else:
            retent.append(np.asarray(v.retentate_signal, dtype=float))

        # Permeate signal is optional.
        if v.permeate_signal is None:
            permea.append(np.full_like(t, np.nan, dtype=float))
        else:
            permea.append(np.asarray(v.permeate_signal, dtype=float))

        # Per-vial "window" in local coordinates.
        windows.append((0, max(len(t) - 1, 0)))

    return {
        "times": times,
        "masses": masses,
        "retentate_cond": retent,
        "permeate_cond": permea,
        "windows": windows,
        "validation": {"ok_fatal": bool(ok_fatal), "issues": list(issues)},
    }


def _plot_quicklook(
    payload: Dict[str, object],
    *,
    title: str,
    save_prefix: Optional[str],
) -> None:
    """Quicklook plots (mass and retentate/permeate signals) for the legacy dict."""

    # Import matplotlib lazily so headless / minimal environments can still load data.
    import matplotlib.pyplot as plt

    times: List[np.ndarray] = payload.get("times", [])
    masses: List[np.ndarray] = payload.get("masses", [])
    retent: List[np.ndarray] = payload.get("retentate_cond", [])
    permea: List[np.ndarray] = payload.get("permeate_cond", [])

    # ---------------------------
    # Mass plot
    # ---------------------------
    plt.figure()
    for t, m in zip(times, masses):
        if len(t) == 0:
            continue
        plt.plot(t, m)
    plt.xlabel("Time (s)")
    plt.ylabel("Mass (g)")
    plt.title(f"Mass :: {title}")
    plt.tight_layout()
    if save_prefix:
        plt.savefig(f"{save_prefix}_mass.png", dpi=200)

    # ---------------------------
    # Conductivity signals plot
    # ---------------------------
    plt.figure()
    for t, r, p in zip(times, retent, permea):
        if len(t) == 0:
            continue
        plt.plot(t, r)
        # Permeate often exists but may be NaN-only.
        if np.isfinite(p).any():
            plt.plot(t, p)
    plt.xlabel("Time (s)")
    plt.ylabel("Signal (uS/cm or source units)")
    plt.title(f"Signals :: {title}")
    plt.tight_layout()
    if save_prefix:
        plt.savefig(f"{save_prefix}_signals.png", dpi=200)

    # If called in scripts, don't block execution.
    plt.close("all")
