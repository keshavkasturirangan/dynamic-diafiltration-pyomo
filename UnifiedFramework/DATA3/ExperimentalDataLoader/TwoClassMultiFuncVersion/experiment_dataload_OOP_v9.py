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
