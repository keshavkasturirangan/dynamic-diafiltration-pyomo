#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 07:12:33 2025

@author: Keshav Kasturi Rangan, Alex Dowling
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, List, Union, Tuple
from pathlib import Path
import pandas as pd

import numpy as np


class SourceType(str, Enum):
    """
    Enumerates supported experiment data sources.

    Values:
        XLSX: Excel workbook input (.xlsx)
        MAT:  MATLAB file input (.mat)
    """
    XLSX = "xlsx"
    MAT = "mat"


@dataclass
class VialData:
    """
    Stores all data for ONE vial segment (one contiguous interval between vial swaps).

    This class is intentionally "light":
        - It does not load files
        - It does not compute anything
        - It only stores data as explicit, readable fields

    Core idea:
        Each vial segment holds time-aligned arrays (time, signals) plus optional
        per-vial scalar assay values (e.g., ICP concentration).

    Legacy compatibility:
        cV_avg may be stored as:
            - a scalar (float), OR
            - an array (np.ndarray), OR
            - None
        because historical datasets sometimes encode permeate concentration either way.
    """

    # -------------------------
    # Identity / indexing
    # -------------------------
    number: Optional[int] = None  # vial index (1-based); None means "not set yet"

    # -------------------------
    # Core time-aligned arrays
    # -------------------------
    time_s: Optional[np.ndarray] = None  # time stamps [s] for this vial; None until loaded
    mass_g: Optional[np.ndarray] = None  # mass [g] aligned with time_s; optional

    # Honest retentate signal (conductivity unless converted later)
    retentate_signal: Optional[np.ndarray] = None       # e.g., retentate conductivity time series
    retentate_signal_name: Optional[str] = None         # e.g., "retentate_cond_uS_cm"
    retentate_signal_units: Optional[str] = None        # e.g., "uS/cm"

    # -------------------------
    # Optional additional time-series (aligned with time_s if present)
    # -------------------------
    pressure_psi: Optional[np.ndarray] = None
    retentate_temp_C: Optional[np.ndarray] = None
    permeate_temp_C: Optional[np.ndarray] = None

    permeate_signal: Optional[np.ndarray] = None
    permeate_signal_name: Optional[str] = None
    permeate_signal_units: Optional[str] = None

    vial_swap_flag: Optional[np.ndarray] = None

    # -------------------------
    # Optional per-vial assays (scalars)
    # -------------------------
    cV_assay_value: Optional[float] = None
    cV_assay_units: Optional[str] = None

    sample_volume_mL: Optional[float] = None
    nitric_acid_volume_mL: Optional[float] = None

    # -------------------------
    # Legacy-compatible permeate storage (scalar OR array)
    # -------------------------
    cV_avg: Optional[Union[float, np.ndarray]] = None


@dataclass
class ExperimentalData:
    """
    Stores ONE experiment (one Excel sheet OR one MATLAB struct) as explicit fields.

    This class is intentionally "light":
        - It does not load files
        - It does not parse tables
        - It does not validate itself
        - It only stores data in explicit fields

    Missing values:
        All fields default to None so loaders can populate what exists and leave the rest empty.
    """

    # -------------------------
    # Provenance / metadata
    # -------------------------
    source: Optional[SourceType] = None  # use Enum (prevents typos)
    filename: Optional[str] = None
    sheet_name: Optional[str] = None

    # -------------------------
    # Experiment identity
    # -------------------------
    dataset_id: Optional[float] = None
    mode: Optional[str] = None

    # -------------------------
    # Measurement behavior
    # -------------------------
    continuous_retentate: Optional[bool] = None

    # -------------------------
    # Optional operating conditions / config
    # -------------------------
    delP_bar: Optional[float] = None
    Temp_K: Optional[float] = None
    Am_cm2: Optional[float] = None
    rho_g_cm3: Optional[float] = None

    # -------------------------
    # Optional feed/dialysate info
    # -------------------------
    M_F0_g: Optional[float] = None
    C_F0_value: Optional[float] = None
    C_F0_units: Optional[str] = None

    M_O_g: Optional[float] = None
    C_D_value: Optional[float] = None
    C_D_units: Optional[str] = None

    # -------------------------
    # Optional component bookkeeping
    # -------------------------
    component_names: Optional[List[str]] = None
    num_components: Optional[int] = None

    # -------------------------
    # Optional initial guesses / priors (theta0 is numeric)
    # -------------------------
    Lp0: Optional[float] = None
    B0: Optional[float] = None
    sigma0: Optional[float] = None
    theta0: Optional[np.ndarray] = None  # numeric vector only

    # -------------------------
    # Optional calibration storage
    # -------------------------
    salt1_calib_x: Optional[np.ndarray] = None
    salt1_calib_y: Optional[np.ndarray] = None
    salt2_calib_x: Optional[np.ndarray] = None
    salt2_calib_y: Optional[np.ndarray] = None

    # -------------------------
    # Vial list (ordered)
    # -------------------------
    vials: List[VialData] = field(default_factory=list)


def is_missing_scalar(x: object) -> bool:
    """
    Determine whether a scalar-like value should be treated as "missing".

    This helper centralizes missing-value logic so that:
        - Excel inputs (empty cells, NaN)
        - MATLAB inputs (NaN, empty arrays collapsed to scalars)
        - Python-native None / empty strings
    are handled consistently across loaders and validators.

    Inputs:
        x:
            Any scalar-like value (float, int, str, None, etc.).

    Output:
        True  -> value is missing / unusable
        False -> value is present
    """

    # Case 1: Python None always means missing
    if x is None:
        return True

    # Case 2: Numeric NaN (covers float('nan'), numpy.nan, MATLAB NaN)
    # We check type explicitly to avoid raising errors on non-numeric objects
    if isinstance(x, float) and np.isnan(x):
        return True

    # Case 3: Empty or whitespace-only strings (common in Excel metadata cells)
    if isinstance(x, str) and x.strip() == "":
        return True

    # Otherwise, treat the value as present
    return False

def as_float_array(x: object) -> np.ndarray:
    """
    Convert an input into a 1D NumPy array of floats (best-effort).

    Why this exists:
        Excel and MATLAB loaders may hand us values in many shapes/types:
            - Python lists/tuples
            - NumPy arrays
            - pandas Series
            - MATLAB-loaded numeric arrays (often already NumPy arrays)
        This helper standardizes all of those into a single, predictable representation:
            a 1D numpy.ndarray with dtype=float.

    Design rules:
        1) If x is None -> return an empty float array.
        2) Otherwise -> try to convert to np.asarray(..., dtype=float).
        3) Always flatten to 1D so downstream slicing/indexing is consistent.

    Inputs:
        x:
            Array-like object (list, tuple, np.ndarray, pandas Series, etc.)
            or None.

    Output:
        np.ndarray:
            1D float array. If x is None, returns np.array([], dtype=float).
    """

    # If the caller provides nothing, return a standard empty float array.
    # This avoids special-casing None everywhere else in the code.
    if x is None:
        return np.array([], dtype=float)

    # Convert x into a NumPy array of floats.
    # np.asarray is used (instead of np.array) because it avoids unnecessary copying.
    # dtype=float enforces numeric conversion (e.g., strings will raise in real code).
    arr = np.asarray(x, dtype=float)

    # Flatten to 1D so that all loaders and helper functions can assume a simple shape.
    # This also handles MATLAB-style (N,1) and (1,N) arrays.
    return arr.reshape(-1)

def as_int_array(x: object) -> np.ndarray:
    """
    Convert an input into a 1D NumPy array of integers (best-effort).

    Why this exists:
        Certain signals (most importantly vial swap flags) are conceptually
        discrete indicators rather than continuous measurements.

        In Excel and MATLAB, these often appear as:
            - 0 / 1 integers
            - floats like 0.0 / 1.0
            - arrays with NaNs mixed in
            - pandas Series

        This helper standardizes all of those into a clean 1D integer array
        suitable for segmentation logic.

    Design rules:
        1) If x is None -> return an empty integer array.
        2) Convert input to numeric form first (via float).
        3) Treat NaNs as zeros (i.e., "no swap").
        4) Cast to int.
        5) Always flatten to 1D.

    Inputs:
        x:
            Array-like object (list, tuple, np.ndarray, pandas Series, etc.)
            or None.

    Output:
        np.ndarray:
            1D integer array. If x is None, returns np.array([], dtype=int).
    """

    # If nothing was provided, return a standard empty integer array.
    # This mirrors as_float_array and avoids None checks downstream.
    if x is None:
        return np.array([], dtype=int)

    # First convert to a NumPy array of floats.
    # This allows us to safely handle:
    #   - 0 / 1
    #   - 0.0 / 1.0
    #   - NaN values (common in Excel)
    arr = np.asarray(x, dtype=float)

    # Replace NaN values with 0.0.
    # Interpretation:
    #   NaN in a swap flag column means "no swap occurred".
    # This is a deliberate, conservative choice.
    arr = np.nan_to_num(arr, nan=0.0)

    # Convert to integers.
    # After NaN handling, values should be cleanly castable.
    arr = arr.astype(int)

    # Flatten to 1D so segmentation logic can assume simple indexing.
    return arr.reshape(-1)

def segment_indices_by_swap(swap_flags: np.ndarray) -> List[Tuple[int, int]]:
    """
    Convert a 0/1 "vial swap" flag array into a list of (start, end) index segments.

    Why this exists:
        In the Excel workflow, the experiment time-series is one long table.
        "Vial Swap" marks the time index where one vial ends and the next begins.

        Downstream, we want an *ordered* list of vial segments so we can slice:
            time_all[a:b], mass_all[a:b], ret_signal_all[a:b], ...

    Convention (important and explicit):
        If swap_flags[idx] == 1, then the row at idx is included in the *current* vial.
        That means the current vial segment ends at end = idx + 1 (Python slicing).

        Example (indices):       0  1  2  3  4  5
        swap_flags:              0  0  1  0  1  0
        segments returned:      [0:3], [3:5], [5:6]
        Explanation:
            - swap at idx=2 ends vial 1 at end=3
            - swap at idx=4 ends vial 2 at end=5
            - remaining rows form vial 3 [5:6]

    Inputs:
        swap_flags (np.ndarray):
            1D integer array (typically from as_int_array), same length as the time-series arrays.

    Output:
        List[Tuple[int, int]]:
            List of (start, end) pairs, each representing a Python slice [start:end].
            These segments are ordered and non-overlapping.
    """

    # Defensive: ensure we have a 1D array (even if caller passed (N,1) accidentally)
    swap_flags = np.asarray(swap_flags).reshape(-1)

    # Find indices where swap == 1 (i.e., vial boundary markers)
    swap_idx = np.where(swap_flags == 1)[0]

    # This will hold the final list of [start:end] segments
    segments: List[Tuple[int, int]] = []

    # The first vial always starts at index 0
    start = 0

    # For each swap marker, close the current vial and start the next one
    for idx in swap_idx:
        end = idx + 1               # include the swap row in the current vial
        segments.append((start, end))
        start = end                 # next vial begins immediately after

    # After processing all swap markers, add the final segment to the end (if any data remains)
    if start < len(swap_flags):
        segments.append((start, len(swap_flags)))

    return segments

def build_vial_from_slice(
    *,
    number: int,
    a: int,
    b: int,
    time_all: np.ndarray,
    mass_all: Optional[np.ndarray],
    ret_signal_all: np.ndarray,
    ret_signal_name: str,
    ret_signal_units: str,
    swap_all: Optional[np.ndarray] = None,
    pressure_all: Optional[np.ndarray] = None,
    ret_temp_all: Optional[np.ndarray] = None,
    perm_temp_all: Optional[np.ndarray] = None,
    perm_signal_all: Optional[np.ndarray] = None,
    perm_signal_name: Optional[str] = None,
    perm_signal_units: Optional[str] = None,
) -> VialData:
    """
    Create and populate a VialData object by slicing time-aligned experiment arrays [a:b].

    What this function does:
        - Takes full-experiment arrays (time, mass, conductivity, etc.)
        - Takes a segment [a:b] that represents one vial (from segment_indices_by_swap)
        - Creates a VialData container and fills its fields with the sliced arrays
        - Stores signal names/units explicitly to prevent ambiguity later

    What this function does NOT do:
        - It does not read files
        - It does not validate the experiment
        - It does not convert conductivity to concentration
        - It does not drop NaNs in continuous signals (it preserves them)

    Inputs:
        number:
            Vial index (1-based).

        a, b:
            Slice bounds in Python style [a:b] (end is exclusive).

        time_all:
            Full experiment time array (1D).

        mass_all:
            Full experiment mass array (1D) aligned with time_all, or None if unavailable.

        ret_signal_all:
            Full experiment retentate signal array (1D) aligned with time_all.
            For Excel, this is retentate conductivity (uS/cm), stored honestly.

        ret_signal_name / ret_signal_units:
            Labels for retentate_signal (e.g., "retentate_cond_uS_cm", "uS/cm").

        swap_all:
            Full experiment swap flag array (1D) aligned with time_all (optional).
            If provided, a slice is stored in the vial as provenance/debug info.

        pressure_all, ret_temp_all, perm_temp_all:
            Optional full-experiment arrays aligned with time_all.

        perm_signal_all:
            Optional permeate signal array aligned with time_all (e.g., permeate conductivity).
            If provided, perm_signal_name/units should also be provided.

    Output:
        VialData:
            A populated vial container whose array fields have length (b - a).
    """

    # Create an empty vial container (VialData is intentionally "light").
    v = VialData()

    # Store the vial index immediately (helps debugging if something fails later).
    v.number = number

    # -------------------------
    # Slice the required arrays
    # -------------------------

    # Time is always required for a vial (it is the alignment axis).
    v.time_s = time_all[a:b]

    # Retentate signal is required (conductivity for Excel; label honestly via name/units).
    v.retentate_signal = ret_signal_all[a:b]
    v.retentate_signal_name = ret_signal_name
    v.retentate_signal_units = ret_signal_units

    # -------------------------
    # Slice optional arrays only if they exist
    # -------------------------

    # Mass is optional in some datasets (so only slice if present).
    if mass_all is not None:
        v.mass_g = mass_all[a:b]

    # Swap flag slice is optional, but useful for provenance and debugging segmentation.
    if swap_all is not None:
        v.vial_swap_flag = swap_all[a:b]

    # Optional operating-condition signals aligned with time.
    if pressure_all is not None:
        v.pressure_psi = pressure_all[a:b]

    if ret_temp_all is not None:
        v.retentate_temp_C = ret_temp_all[a:b]

    if perm_temp_all is not None:
        v.permeate_temp_C = perm_temp_all[a:b]

    # Optional permeate signal (e.g., permeate conductivity).
    if perm_signal_all is not None:
        v.permeate_signal = perm_signal_all[a:b]
        v.permeate_signal_name = perm_signal_name
        v.permeate_signal_units = perm_signal_units

    # Return the populated vial container.
    return v

def validate_experiment(exp: ExperimentalData) -> Tuple[bool, List[Tuple[str, str]]]:
    """
    Validate basic structural integrity of a loaded experiment and FLAG issues (no exceptions).

    Philosophy:
        - Does NOT modify exp.
        - Does NOT raise exceptions.
        - Reports problems using (severity, message) tuples.

    Severity meanings:
        FATAL:
            The experiment is structurally unusable; downstream code is likely to fail.
        REQUIRED:
            Likely needed for Pyomo / estimation later, but may be missing at ingestion time.

    Inputs:
        exp:
            ExperimentalData object produced by a loader.

    Outputs:
        (ok_fatal, issues)
            ok_fatal:
                True if no FATAL issues were found, else False.
            issues:
                List of (severity, message) tuples.
    """

    # Collect all issues found during validation (we try to report everything, not just first failure).
    issues: List[Tuple[str, str]] = []

    # Tracks whether any FATAL errors were found.
    ok_fatal = True

    # -------------------------------------------------------------------------
    # 1) FATAL: must have at least one vial
    # -------------------------------------------------------------------------
    if exp.vials is None or len(exp.vials) == 0:
        issues.append(("FATAL", "No vials found (exp.vials is empty)."))
        return False, issues  # no further checks possible

    # -------------------------------------------------------------------------
    # 2) FATAL: each vial must have time and retentate signal, aligned in length
    # -------------------------------------------------------------------------
    for v in exp.vials:
        # Vial must have a time vector
        if v.time_s is None or len(v.time_s) == 0:
            issues.append(("FATAL", f"Vial {v.number}: missing/empty time_s."))
            ok_fatal = False

        # Vial must have a retentate signal (conductivity stored honestly)
        if v.retentate_signal is None or len(v.retentate_signal) == 0:
            issues.append(("FATAL", f"Vial {v.number}: missing/empty retentate_signal."))
            ok_fatal = False

        # If both exist, they must match in length (alignment requirement)
        if (v.time_s is not None) and (v.retentate_signal is not None):
            if len(v.time_s) != len(v.retentate_signal):
                issues.append(("FATAL", f"Vial {v.number}: time_s and retentate_signal lengths differ."))
                ok_fatal = False

        # Optional structural check: time should not go backwards
        # This is often required by ODE solvers and interpolation routines.
        if v.time_s is not None and len(v.time_s) > 1:
            if np.any(np.diff(v.time_s) < 0):
                issues.append(("FATAL", f"Vial {v.number}: time_s is not monotonic non-decreasing."))
                ok_fatal = False

    # -------------------------------------------------------------------------
    # 3) REQUIRED (flag only): experiment-level operating conditions
    # -------------------------------------------------------------------------
    if is_missing_scalar(exp.delP_bar):
        issues.append(("REQUIRED", "Missing exp.delP_bar (pressure)."))

    if is_missing_scalar(exp.Temp_K):
        issues.append(("REQUIRED", "Missing exp.Temp_K (temperature)."))

    if is_missing_scalar(exp.Am_cm2):
        issues.append(("REQUIRED", "Missing exp.Am_cm2 (membrane area)."))

    # -------------------------------------------------------------------------
    # 4) REQUIRED (flag only): theta0 must be numeric if provided
    # -------------------------------------------------------------------------
    # Per spec, theta0 is always numeric when present (None means not provided).
    if exp.theta0 is not None:
        try:
            _ = np.asarray(exp.theta0, dtype=float)
        except Exception:
            issues.append(("REQUIRED", "exp.theta0 exists but is not numeric/castable to float array."))

    return ok_fatal, issues

def detect_source_type(file_path: Path) -> Optional[SourceType]:
    """
    Detect the data source type from a file extension.

    Why this exists:
        load_experiment(...) needs a simple, reliable way to decide whether the input
        should be routed to the XLSX loader or the MAT loader.

    Inputs:
        file_path:
            A pathlib.Path object pointing to the input file.

    Output:
        SourceType.XLSX if the extension is .xlsx or .xls
        SourceType.MAT  if the extension is .mat
        None             if the extension is unsupported/unknown
    """

    # Normalize the extension to lowercase so ".XLSX" works the same as ".xlsx".
    suffix = file_path.suffix.lower()

    # Excel inputs (one workbook, one experiment per sheet in your convention).
    if suffix in [".xlsx", ".xls"]:
        return SourceType.XLSX

    # MATLAB structured inputs (legacy .mat experiment struct).
    if suffix == ".mat":
        return SourceType.MAT

    # Any other extension is not supported by this loader entrypoint.
    return None


def load_experiment(path: str, selector: object) -> Tuple[ExperimentalData, Tuple[bool, List[Tuple[str, str]]]]:
    """
    File-agnostic entrypoint: load ONE experiment from either XLSX or MAT.

    This is the function you call from "main code" (or from tests) without caring
    about file format.

    Inputs:
        path:
            Path to the input file. Supported extensions:
                - .xlsx/.xls : Excel workbook
                - .mat       : MATLAB file

        selector:
            Identifies WHICH experiment to load from the file:
                - XLSX:
                    sheet name (str) OR sheet index (int)
                    (because one workbook can contain multiple experiments)
                - MAT:
                    struct key (str) OR None
                    (None means "use default key", e.g., 'data_stru')

    Outputs:
        (exp, validation)

        exp:
            An ExperimentalData object populated by the appropriate loader.

        validation:
            A tuple (ok_fatal, issues) returned by validate_experiment(exp)
            where:
                ok_fatal: False means the experiment is structurally unusable
                issues:   list of (severity, message) flags
    """

    # Convert the input string path into a Path object for robust extension handling.
    file_path = Path(path)

    # Determine whether this is an Excel workbook or a MATLAB .mat file.
    kind = detect_source_type(file_path)

    # If we cannot determine a supported type, fail fast with a clear message.
    # (This is not "data missing"; it is an unsupported file format.)
    if kind is None:
        raise ValueError(
            f"Unsupported file type '{file_path.suffix}'. "
            f"Expected one of: .xlsx, .xls, .mat"
        )

    # Route to the correct loader based on file type.
    # NOTE: load_from_xlsx and load_from_mat are format-specific loaders you will implement next.
    if kind == SourceType.XLSX:
        # For Excel: selector identifies which SHEET corresponds to the experiment.
        exp = load_from_xlsx(file_path, sheet_selector=selector)

    elif kind == SourceType.MAT:
        # For MAT: selector identifies the STRUCT KEY (or None for default).
        exp = load_from_mat(file_path, struct_key=selector)

    else:
        # This should not occur because detect_source_type only returns XLSX/MAT/None,
        # but keeping it makes the logic robust if the enum grows later.
        raise ValueError(f"Unhandled SourceType: {kind}")

    # Run validation as the last step (read-only; does not mutate exp).
    validation = validate_experiment(exp)

    # Return the loaded experiment and its validation report.
    return exp, validation

def load_from_xlsx(xlsx_path: Path, sheet_selector: object) -> ExperimentalData:
    """
    Load ONE experiment from an Excel workbook sheet into ExperimentalData fields.

    Key idea:
        - One workbook can contain multiple experiments.
        - Each SHEET corresponds to ONE experiment (your project convention).
        - This function loads exactly ONE sheet/experiment and returns one ExperimentalData object.

    What this function does:
        1) Reads the requested sheet (I/O stub).
        2) Parses a time-series table from that sheet (parsing stub).
        3) Extracts required arrays (time, mass, swap flags, retentate conductivity).
        4) Segments the experiment into vials using swap flags (shared helper).
        5) Builds one VialData per segment using build_vial_from_slice (shared helper).
        6) Optionally attaches per-vial assay values from a vial table (parsing stub).
        7) Returns an ExperimentalData object with an ordered list of vials.

    What this function does NOT do:
        - It does not validate the experiment (caller runs validate_experiment).
        - It does not convert conductivity to concentration.
        - It does not guess missing metadata; missing fields remain None.

    Inputs:
        xlsx_path:
            Path to the Excel workbook (.xlsx or .xls).

        sheet_selector:
            Sheet name (str) or sheet index (int) identifying the experiment sheet to load.

    Output:
        ExperimentalData:
            Populated experiment container with vials.
    """

    # -------------------------------------------------------------------------
    # 1) Read the sheet (format-specific I/O, implemented later)
    # -------------------------------------------------------------------------
    # This should return something you can pass into parsing helpers below.
    # In real code this may be a pandas DataFrame, an openpyxl worksheet, or a dict of tables.
    sheet = read_excel_sheet(xlsx_path, sheet_selector)

    # -------------------------------------------------------------------------
    # 2) Parse the time-series table from the sheet (format-specific parsing)
    # -------------------------------------------------------------------------
    # The returned object must support column access by name, e.g. ts["Time (s)"].
    # In real code, "ts" will likely be a pandas DataFrame.
    ts = parse_time_series_table(sheet)

    # -------------------------------------------------------------------------
    # 3) Extract REQUIRED arrays from the time-series table
    # -------------------------------------------------------------------------
    # These are required to construct vials and run basic validation.

    # Time [s] (required)
    time_all = as_float_array(ts["Time (s)"])

    # Mass [g] (required in your current Excel format; keep it required here)
    # If you later decide mass is optional, you can change this to a guarded lookup.
    mass_all = as_float_array(ts["Mass (g)"])

    # Vial swap flag (required for segmentation)
    swap_all = as_int_array(ts["Vial Swap"])

    # Retentate conductivity signal (required; stored honestly as conductivity)
    ret_cond_all = as_float_array(ts["Retentate Cond @ Temp (uS/cm)"])
    ret_signal_name = "retentate_cond_uS_cm"
    ret_signal_units = "uS/cm"

    # -------------------------------------------------------------------------
    # 4) Extract OPTIONAL arrays from the time-series table (if columns exist)
    # -------------------------------------------------------------------------
    # We do NOT assume these columns exist in every Excel sheet.

    pressure_all = (
        as_float_array(ts["Pressure (psi)"])
        if table_has_column(ts, "Pressure (psi)")
        else None
    )

    ret_temp_all = (
        as_float_array(ts["Retentate Temp"])
        if table_has_column(ts, "Retentate Temp")
        else None
    )

    perm_temp_all = (
        as_float_array(ts["Permeate Temp"])
        if table_has_column(ts, "Permeate Temp")
        else None
    )

    # Optional permeate conductivity (if present, store honestly as a signal)
    perm_cond_all = (
        as_float_array(ts["Permeate Cond @ Temp (uS/cm)"])
        if table_has_column(ts, "Permeate Cond @ Temp (uS/cm)")
        else None
    )
    perm_signal_name = "permeate_cond_uS_cm"
    perm_signal_units = "uS/cm"

    # -------------------------------------------------------------------------
    # 5) Parse the vial assay table (format-specific parsing; may or may not exist)
    # -------------------------------------------------------------------------
    # This table contains per-vial assay information (e.g., ICP measurements).
    # If the sheet does not have it, parse_vial_data_table may return None.
    vial_table = parse_vial_data_table(sheet)

    # -------------------------------------------------------------------------
    # 6) Segment the time-series into vial slices using swap flags
    # -------------------------------------------------------------------------
    segments = segment_indices_by_swap(swap_all)

    # -------------------------------------------------------------------------
    # 7) Create the experiment container (explicit fields only)
    # -------------------------------------------------------------------------
    exp = ExperimentalData()                # empty experiment container
    exp.source = SourceType.XLSX            # provenance: data came from Excel
    exp.filename = xlsx_path.name           # provenance: file name
    exp.sheet_name = get_sheet_name(sheet)  # provenance: sheet name (if helper can infer it)

    # Excel experiments: retentate is continuously collected (per your spec)
    exp.continuous_retentate = True

    # Note: Excel may not contain dataset_id/mode/operating conditions; leave as None here.
    # These can be filled later from metadata if available.

    # -------------------------------------------------------------------------
    # 8) Build one VialData object per segment
    # -------------------------------------------------------------------------
    for i, (a, b) in enumerate(segments, start=1):
        # Build vial from aligned arrays using the shared helper.
        v = build_vial_from_slice(
            number=i,
            a=a,
            b=b,
            time_all=time_all,
            mass_all=mass_all,
            ret_signal_all=ret_cond_all,
            ret_signal_name=ret_signal_name,
            ret_signal_units=ret_signal_units,
            swap_all=swap_all,
            pressure_all=pressure_all,
            ret_temp_all=ret_temp_all,
            perm_temp_all=perm_temp_all,
            perm_signal_all=perm_cond_all,
            perm_signal_name=perm_signal_name,
            perm_signal_units=perm_signal_units,
        )

        # ---------------------------------------------------------------------
        # 9) Attach per-vial assay information (optional)
        # ---------------------------------------------------------------------
        # If vial_table is present, find the row corresponding to this vial and store values.
        # This keeps ingestion honest: if assays aren't available, these stay None.
        if vial_table is not None:
            assay_row = find_vial_row(vial_table, vial_number=i)

            if assay_row is not None:
                # Example mapping (you will finalize exact column names during implementation):
                # Store assay value + units explicitly as a scalar measurement.
                v.cV_assay_value = assay_row.get("ICP Salt 1 (mg/L)", None)
                v.cV_assay_units = "mg/L"

                # Store sample/dilution metadata if present
                v.sample_volume_mL = assay_row.get("Sample Volume (mL)", None)
                v.nitric_acid_volume_mL = assay_row.get("Nitric Acid Volume (mL)", None)

                # Legacy-friendly: store cV_avg as scalar by default (same as assay)
                v.cV_avg = v.cV_assay_value

        # Add vial to experiment (maintains order)
        exp.vials.append(v)

    # -------------------------------------------------------------------------
    # 10) Return the populated experiment container (validation is done by caller)
    # -------------------------------------------------------------------------
    return exp


# =============================================================================
# Format-specific I/O/parsing stubs (to implement next)
# =============================================================================

def read_excel_sheet(xlsx_path: Path, sheet_selector: object) -> dict:
    """
    Read an Excel workbook and resolve a single sheet into a lightweight "sheet handle".

    Why this exists:
        - One workbook may contain many experiments (one per sheet).
        - We want to resolve the intended sheet robustly (by name or by index).
        - We do NOT parse tables here; we only prepare a handle for parsing helpers.

    Inputs:
        xlsx_path:
            Path to the Excel workbook (.xlsx or .xls).

        sheet_selector:
            Either:
                - sheet name (str), OR
                - sheet index (int), 0-based.

    Output:
        dict:
            A lightweight sheet handle with:
                - "excel_file": pandas.ExcelFile object (holds workbook metadata)
                - "sheet_name": resolved sheet name (str)
                - "path": workbook path (Path)

            This dict is intentionally simple and explicit so downstream parsing code
            can do things like:
                pd.read_excel(sheet["path"], sheet_name=sheet["sheet_name"], ...)
            or:
                pd.read_excel(sheet["excel_file"], sheet_name=sheet["sheet_name"], ...)
    """

    # Create an ExcelFile object (this reads workbook structure once and caches it).
    # This is efficient when you plan to read multiple ranges/tables from the same sheet.
    xls = pd.ExcelFile(xlsx_path)

    # Pull available sheet names (order matters for index selection).
    names = list(xls.sheet_names)

    # Resolve selector -> sheet name
    if isinstance(sheet_selector, int):
        # Interpret as 0-based sheet index.
        if sheet_selector < 0 or sheet_selector >= len(names):
            raise IndexError(
                f"Sheet index {sheet_selector} out of range. "
                f"Workbook has {len(names)} sheets: {names}"
            )
        sheet_name = names[sheet_selector]

    elif isinstance(sheet_selector, str):
        # Interpret as sheet name.
        if sheet_selector not in names:
            raise ValueError(
                f"Sheet name '{sheet_selector}' not found in workbook. "
                f"Available sheets: {names}"
            )
        sheet_name = sheet_selector

    else:
        # Any other type is ambiguous / unsupported.
        raise TypeError(
            f"sheet_selector must be int (sheet index) or str (sheet name), "
            f"got {type(sheet_selector).__name__}."
        )

    # Return a lightweight handle (explicit, readable; not a DataFrame yet).
    return {
        "excel_file": xls,
        "sheet_name": sheet_name,
        "path": xlsx_path,
    }

def parse_time_series_table(sheet: object) -> pd.DataFrame:
    """
    Parse and return the time-series table from an Excel sheet handle.

    What this function does:
        - Reads the entire sheet (or enough of it) as raw cells (no header).
        - Locates the row that contains the true column headers for the time-series table.
        - Constructs a DataFrame whose columns are those headers.
        - Returns only the rows below that header row (i.e., the actual data).

    What this function does NOT do:
        - It does not validate units or values.
        - It does not convert conductivity to concentration.
        - It does not segment into vials.
        - It does not drop NaNs in continuous signals (NaNs are preserved).

    Inputs:
        sheet:
            The object returned by read_excel_sheet(...).
            In our implementation, this is a dict-like "sheet handle" containing:
                - "excel_file": pandas.ExcelFile
                - "sheet_name": str
                - "path": Path

    Output:
        pd.DataFrame:
            A DataFrame representing the time-series table.
            It must contain (at minimum) these required columns:
                - "Time (s)"
                - "Mass (g)"
                - "Vial Swap"
                - "Retentate Cond @ Temp (uS/cm)"

    Raises:
        ValueError:
            If the header row cannot be found or required columns are missing.
        TypeError:
            If the input 'sheet' is not in the expected handle format.
    """

    # -------------------------------------------------------------------------
    # 0) Defensive: ensure sheet is the handle we expect (from read_excel_sheet)
    # -------------------------------------------------------------------------
    if not isinstance(sheet, dict):
        raise TypeError(
            "parse_time_series_table expected 'sheet' to be a dict-like sheet handle "
            "(from read_excel_sheet)."
        )

    # Pull out the information needed to read the sheet again using pandas.
    xls = sheet.get("excel_file", None)      # pandas.ExcelFile object
    sheet_name = sheet.get("sheet_name", None)

    if xls is None or sheet_name is None:
        raise TypeError(
            "Sheet handle is missing required keys. Expected keys: "
            "'excel_file' and 'sheet_name'."
        )

    # -------------------------------------------------------------------------
    # 1) Read the entire sheet as raw cells (no header)
    # -------------------------------------------------------------------------
    # header=None means pandas will label columns as 0,1,2,... and keep all rows raw.
    # dtype=object preserves mixed types (numbers + strings), which is necessary
    # for locating the header row reliably.
    raw = pd.read_excel(xls, sheet_name=sheet_name, header=None, dtype=object)

    # -------------------------------------------------------------------------
    # 2) Define the required columns that identify the time-series table
    # -------------------------------------------------------------------------
    # These are the minimum columns we need to proceed with vial construction.
    required_cols = {
        "Time (s)",
        "Mass (g)",
        "Vial Swap",
        "Retentate Cond @ Temp (uS/cm)",
    }

    # -------------------------------------------------------------------------
    # 3) Find the header row by scanning for the required column names
    # -------------------------------------------------------------------------
    # Excel sheets often have metadata rows above the actual table.
    # So we search for a row where the set of string cells contains our required names.
    header_row_idx: Optional[int] = None

    # Scan the first N rows. 200 is conservative and still fast for typical sheets.
    scan_limit = min(len(raw), 200)

    for i in range(scan_limit):
        # Extract the row as a Python list
        row = raw.iloc[i].tolist()

        # Normalize cells:
        #   - Keep only non-null entries
        #   - Convert to stripped strings
        # This makes matching robust to extra spaces and non-string types.
        normalized = set()
        for cell in row:
            if cell is None or (isinstance(cell, float) and np.isnan(cell)):
                continue
            # Convert to string and strip whitespace.
            s = str(cell).strip()
            if s != "":
                normalized.add(s)

        # If all required column names are present in this row, we found the header.
        if required_cols.issubset(normalized):
            header_row_idx = i
            break

    # If we never found a header row, we cannot parse the time-series table.
    if header_row_idx is None:
        raise ValueError(
            f"Could not locate time-series header row in sheet '{sheet_name}'. "
            f"Expected to find required columns: {sorted(required_cols)}"
        )

    # -------------------------------------------------------------------------
    # 4) Build the time-series DataFrame using that header row
    # -------------------------------------------------------------------------
    # The header row defines the column names.
    header_values = raw.iloc[header_row_idx].tolist()

    # Convert header values to cleaned strings (column names).
    # Empty/NaN headers become "Unnamed: <index>" to keep column positions stable.
    columns: List[str] = []
    for j, h in enumerate(header_values):
        if h is None or (isinstance(h, float) and np.isnan(h)) or str(h).strip() == "":
            columns.append(f"Unnamed: {j}")
        else:
            columns.append(str(h).strip())

    # Data begins immediately after the header row.
    data = raw.iloc[header_row_idx + 1 :].copy()

    # Assign the columns.
    data.columns = columns

    # -------------------------------------------------------------------------
    # 5) Drop rows that are completely empty (common at bottom of Excel sheets)
    # -------------------------------------------------------------------------
    data = data.dropna(how="all")

    # -------------------------------------------------------------------------
    # 6) Final: ensure required columns exist exactly as expected
    # -------------------------------------------------------------------------
    missing = [c for c in required_cols if c not in data.columns]
    if missing:
        raise ValueError(
            f"Time-series table found in sheet '{sheet_name}', but missing required columns: {missing}. "
            f"Available columns: {list(data.columns)}"
        )

    # Return the parsed time-series table (as a DataFrame).
    return data

def parse_time_series_table(sheet: object) -> pd.DataFrame:
    """
    Parse and return the time-series table from an Excel sheet handle.

    What this function does:
        - Reads the entire sheet (or enough of it) as raw cells (no header).
        - Locates the row that contains the true column headers for the time-series table.
        - Constructs a DataFrame whose columns are those headers.
        - Returns only the rows below that header row (i.e., the actual data).

    What this function does NOT do:
        - It does not validate units or values.
        - It does not convert conductivity to concentration.
        - It does not segment into vials.
        - It does not drop NaNs in continuous signals (NaNs are preserved).

    Inputs:
        sheet:
            The object returned by read_excel_sheet(...).
            In our implementation, this is a dict-like "sheet handle" containing:
                - "excel_file": pandas.ExcelFile
                - "sheet_name": str
                - "path": Path

    Output:
        pd.DataFrame:
            A DataFrame representing the time-series table.
            It must contain (at minimum) these required columns:
                - "Time (s)"
                - "Mass (g)"
                - "Vial Swap"
                - "Retentate Cond @ Temp (uS/cm)"

    Raises:
        ValueError:
            If the header row cannot be found or required columns are missing.
        TypeError:
            If the input 'sheet' is not in the expected handle format.
    """

    # -------------------------------------------------------------------------
    # 0) Defensive: ensure sheet is the handle we expect (from read_excel_sheet)
    # -------------------------------------------------------------------------
    if not isinstance(sheet, dict):
        raise TypeError(
            "parse_time_series_table expected 'sheet' to be a dict-like sheet handle "
            "(from read_excel_sheet)."
        )

    # Pull out the information needed to read the sheet again using pandas.
    xls = sheet.get("excel_file", None)      # pandas.ExcelFile object
    sheet_name = sheet.get("sheet_name", None)

    if xls is None or sheet_name is None:
        raise TypeError(
            "Sheet handle is missing required keys. Expected keys: "
            "'excel_file' and 'sheet_name'."
        )

    # -------------------------------------------------------------------------
    # 1) Read the entire sheet as raw cells (no header)
    # -------------------------------------------------------------------------
    # header=None means pandas will label columns as 0,1,2,... and keep all rows raw.
    # dtype=object preserves mixed types (numbers + strings), which is necessary
    # for locating the header row reliably.
    raw = pd.read_excel(xls, sheet_name=sheet_name, header=None, dtype=object)

    # -------------------------------------------------------------------------
    # 2) Define the required columns that identify the time-series table
    # -------------------------------------------------------------------------
    # These are the minimum columns we need to proceed with vial construction.
    required_cols = {
        "Time (s)",
        "Mass (g)",
        "Vial Swap",
        "Retentate Cond @ Temp (uS/cm)",
    }

    # -------------------------------------------------------------------------
    # 3) Find the header row by scanning for the required column names
    # -------------------------------------------------------------------------
    # Excel sheets often have metadata rows above the actual table.
    # So we search for a row where the set of string cells contains our required names.
    header_row_idx: Optional[int] = None

    # Scan the first N rows. 200 is conservative and still fast for typical sheets.
    scan_limit = min(len(raw), 200)

    for i in range(scan_limit):
        # Extract the row as a Python list
        row = raw.iloc[i].tolist()

        # Normalize cells:
        #   - Keep only non-null entries
        #   - Convert to stripped strings
        # This makes matching robust to extra spaces and non-string types.
        normalized = set()
        for cell in row:
            if cell is None or (isinstance(cell, float) and np.isnan(cell)):
                continue
            # Convert to string and strip whitespace.
            s = str(cell).strip()
            if s != "":
                normalized.add(s)

        # If all required column names are present in this row, we found the header.
        if required_cols.issubset(normalized):
            header_row_idx = i
            break

    # If we never found a header row, we cannot parse the time-series table.
    if header_row_idx is None:
        raise ValueError(
            f"Could not locate time-series header row in sheet '{sheet_name}'. "
            f"Expected to find required columns: {sorted(required_cols)}"
        )

    # -------------------------------------------------------------------------
    # 4) Build the time-series DataFrame using that header row
    # -------------------------------------------------------------------------
    # The header row defines the column names.
    header_values = raw.iloc[header_row_idx].tolist()

    # Convert header values to cleaned strings (column names).
    # Empty/NaN headers become "Unnamed: <index>" to keep column positions stable.
    columns: List[str] = []
    for j, h in enumerate(header_values):
        if h is None or (isinstance(h, float) and np.isnan(h)) or str(h).strip() == "":
            columns.append(f"Unnamed: {j}")
        else:
            columns.append(str(h).strip())

    # Data begins immediately after the header row.
    data = raw.iloc[header_row_idx + 1 :].copy()

    # Assign the columns.
    data.columns = columns

    # -------------------------------------------------------------------------
    # 5) Drop rows that are completely empty (common at bottom of Excel sheets)
    # -------------------------------------------------------------------------
    data = data.dropna(how="all")

    # -------------------------------------------------------------------------
    # 6) Final: ensure required columns exist exactly as expected
    # -------------------------------------------------------------------------
    missing = [c for c in required_cols if c not in data.columns]
    if missing:
        raise ValueError(
            f"Time-series table found in sheet '{sheet_name}', but missing required columns: {missing}. "
            f"Available columns: {list(data.columns)}"
        )

    # Return the parsed time-series table (as a DataFrame).
    return data



def find_vial_row(vial_table: object, vial_number: int) -> Optional[object]:
    """
    Return the row corresponding to a given vial number from the vial assay table.

    Inputs:
        vial_table:
            Table returned by parse_vial_data_table(...)

        vial_number:
            1-based vial index.

    Output:
        object or None:
            A dict-like row (recommended) so you can call row.get("colname", None).
            If the vial is not found, return None.
    """
    raise NotImplementedError


def get_sheet_name(sheet: object) -> Optional[str]:
    """
    Extract the sheet name from the sheet handle, if available.

    Why this exists:
        read_excel_sheet(...) returns a lightweight dict-like "sheet handle" that already
        contains the resolved sheet name. This helper keeps load_from_xlsx(...) readable
        and keeps provenance handling in one place.

    Inputs:
        sheet:
            Dict-like sheet handle returned by read_excel_sheet(...)

    Output:
        str or None:
            The resolved Excel sheet name if present, otherwise None.
    """
    if isinstance(sheet, dict):
        return sheet.get("sheet_name", None)
    return None


def table_has_column(table: object, colname: str) -> bool:
    """
    Return True if 'table' has a column named colname.

    Why this exists:
        Different Excel sheets may include different optional columns.
        This helper allows the loader to safely check for a column before accessing it.

    Inputs:
        table:
            Table-like object. In real code, this will typically be a pandas DataFrame,
            but we keep it generic so the loader is not tightly coupled to pandas.

        colname:
            The exact column name to check for (case-sensitive).

    Output:
        bool:
            True if the column exists, otherwise False.
    """

    # Most common case: pandas DataFrame (has a .columns attribute)
    # We check attribute existence so this does not crash on non-DataFrame inputs.
    if hasattr(table, "columns"):
        # For a DataFrame, membership in table.columns checks column existence.
        return colname in table.columns

    # Fallback case: dict-like object (keys represent column names)
    # This supports cases where parsing returns a mapping instead of a DataFrame.
    if isinstance(table, dict):
        return colname in table

    # If we do not recognize the table type, we cannot guarantee column existence.
    # Being conservative here is safer than guessing.
    return False

def parse_vial_data_table(sheet: object) -> Optional[pd.DataFrame]:
    """
    Parse and return the per-vial assay table from an Excel sheet handle (if present).

    What this function is for:
        Many of your Excel experiment sheets contain TWO logical tables:
            (A) a time-series table (time, mass, conductivity, swap flags, ...)
            (B) a per-vial assay table (e.g., ICP Salt 1 (mg/L), sample volumes, dilution info)

        This function extracts table (B), if it exists.

    What this function does:
        - Reads the entire sheet as raw cells (no header).
        - Searches for a row that contains a recognizable assay header signature.
        - Constructs a DataFrame using that row as column headers.
        - Returns rows below the header row as the assay table.

    What this function does NOT do:
        - It does not validate assay values.
        - It does not assign assays to vials (that is done in load_from_xlsx via find_vial_row).
        - It does not require the table to exist (returns None if not found).

    Inputs:
        sheet:
            Dict-like "sheet handle" returned by read_excel_sheet(...):
                - "excel_file": pandas.ExcelFile
                - "sheet_name": str

    Output:
        Optional[pd.DataFrame]:
            If an assay table is found, returns a DataFrame where each row corresponds to one vial.
            If no assay table is found, returns None.

    Notes:
        This function intentionally uses a "signature search" strategy because Excel layouts
        can vary (blank lines, metadata blocks, etc.).

        You may refine the signature keys once you confirm the exact columns in NF270_MC3.xlsx.
    """

    # -------------------------------------------------------------------------
    # 0) Defensive: ensure sheet is the handle we expect (from read_excel_sheet)
    # -------------------------------------------------------------------------
    if not isinstance(sheet, dict):
        raise TypeError(
            "parse_vial_data_table expected 'sheet' to be a dict-like sheet handle "
            "(from read_excel_sheet)."
        )

    xls = sheet.get("excel_file", None)
    sheet_name = sheet.get("sheet_name", None)

    if xls is None or sheet_name is None:
        raise TypeError(
            "Sheet handle is missing required keys. Expected keys: "
            "'excel_file' and 'sheet_name'."
        )

    # -------------------------------------------------------------------------
    # 1) Read the entire sheet as raw cells (no header)
    # -------------------------------------------------------------------------
    raw = pd.read_excel(xls, sheet_name=sheet_name, header=None, dtype=object)

    # -------------------------------------------------------------------------
    # 2) Define a signature for the assay table header row
    # -------------------------------------------------------------------------
    # The assay table typically includes at least one of these columns.
    # We keep it flexible because different sheets may include different assay columns.
    signature_any = {
        "ICP Salt 1 (mg/L)",
        "ICP Salt 2 (mg/L)",
        "Sample Volume (mL)",
        "Nitric Acid Volume (mL)",
        "Vial",           # some sheets may have a "Vial" label column
        "Vial Number",    # alternative naming
    }

    # We also want some clue that rows correspond to vials.
    # Many assay tables have a first column that contains vial labels (e.g., "Vial 1").
    # We'll detect that later in find_vial_row; here we just find the header row.

    header_row_idx: Optional[int] = None

    # Scan the first N rows looking for any assay signature columns.
    scan_limit = min(len(raw), 300)

    for i in range(scan_limit):
        row = raw.iloc[i].tolist()

        # Build a set of normalized string values in this row.
        normalized = set()
        for cell in row:
            if cell is None or (isinstance(cell, float) and np.isnan(cell)):
                continue
            s = str(cell).strip()
            if s != "":
                normalized.add(s)

        # If we see at least one signature column name, treat this as a candidate header.
        # We use "any" rather than "all" because assay tables may have different subsets of columns.
        if len(signature_any.intersection(normalized)) > 0:
            header_row_idx = i
            break

    # If no assay header row was found, there is no assay table in this sheet.
    if header_row_idx is None:
        return None

    # -------------------------------------------------------------------------
    # 3) Build a DataFrame using the detected header row
    # -------------------------------------------------------------------------
    header_values = raw.iloc[header_row_idx].tolist()

    # Convert header row values into cleaned string column names.
    columns: List[str] = []
    for j, h in enumerate(header_values):
        if h is None or (isinstance(h, float) and np.isnan(h)) or str(h).strip() == "":
            columns.append(f"Unnamed: {j}")
        else:
            columns.append(str(h).strip())

    # Data begins below the assay header row.
    data = raw.iloc[header_row_idx + 1 :].copy()
    data.columns = columns

    # Drop fully empty rows (common at bottom or between blocks).
    data = data.dropna(how="all")

    # -------------------------------------------------------------------------
    # 4) If the resulting table is empty, treat as "not found"
    # -------------------------------------------------------------------------
    if len(data) == 0:
        return None

    # -------------------------------------------------------------------------
    # 5) Return the assay table as a DataFrame
    # -------------------------------------------------------------------------
    return data
