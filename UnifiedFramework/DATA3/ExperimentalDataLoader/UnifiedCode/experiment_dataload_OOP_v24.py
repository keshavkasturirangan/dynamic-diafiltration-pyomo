#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 07:12:33 2025

@author: Keshav Kasturi Rangan, Alex Dowling
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Optional, List, Sequence, Union, Tuple, Dict, Iterable
from pathlib import Path
from contextlib import contextmanager
import pandas as pd
from scipy.io import loadmat  # MATLAB .mat reader (used for initial-guess inference)

import numpy as np
import re
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.contrib.parmest.experiment import Experiment as ParmestExperiment
from pyomo.contrib.parmest.parmest import Estimator
from pyomo.contrib.doe import DesignOfExperiments


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

    # -------------------------------------------------------------------------
    # Identity / indexing
    # -------------------------------------------------------------------------
    number: Optional[int] = None  # vial index (1-based); None means "not set yet"

    # -------------------------------------------------------------------------
    # Core time-aligned arrays
    # -------------------------------------------------------------------------
    time_s: Optional[np.ndarray] = None  # time stamps [s] for this vial; None until loaded
    mass_g: Optional[np.ndarray] = None  # mass [g] aligned with time_s; optional

    # Honest retentate signal (conductivity unless converted later)
    retentate_signal: Optional[np.ndarray] = None       # e.g., retentate conductivity time series
    retentate_signal_name: Optional[str] = None         # e.g., "retentate_cond_uS_cm"
    retentate_signal_units: Optional[str] = None        # e.g., "uS/cm"

    # -------------------------------------------------------------------------
    # Optional: derived concentrations (computed later from conductivity)
    # -------------------------------------------------------------------------
    # We keep these separate from the raw conductivity signals so:
    #   - loaders only save conductivity data (load what was measured),
    #   - a measurement model can convert conductivity -> concentration when desired.

    retentate_concentration: Optional[np.ndarray] = None        # concentration aligned with time_s
    retentate_concentration_units: Optional[str] = None         # e.g., "mM" or "M"
    retentate_concentration_model: Optional[str] = None         # e.g., "variant_shedlovsky" or "msa"

    # -------------------------
    # Optional additional time-series (aligned with time_s if present)
    # -------------------------
    pressure_psi: Optional[np.ndarray] = None
    retentate_temp_C: Optional[np.ndarray] = None
    permeate_temp_C: Optional[np.ndarray] = None

    permeate_signal: Optional[np.ndarray] = None
    permeate_signal_name: Optional[str] = None
    permeate_signal_units: Optional[str] = None
    
    permeate_concentration: Optional[np.ndarray] = None         # concentration aligned with time_s
    permeate_concentration_units: Optional[str] = None          # e.g., "mM" or "M"
    permeate_concentration_model: Optional[str] = None          # e.g., "variant_shedlovsky" or "msa"
    
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

    def sample_count(self) -> int:
        """Return the number of time samples in this vial (0 if missing)."""
        return 0 if self.time_s is None else int(len(self.time_s))

    def time_bounds(self) -> Optional[Tuple[float, float]]:
        """Return (t_start, t_end) in seconds for this vial, or None if missing."""
        if self.time_s is None or len(self.time_s) == 0:
            return None
        return float(self.time_s[0]), float(self.time_s[-1])

    def validate_structure(self) -> Tuple[bool, List[Tuple[str, str]]]:
        """Validate this vial only; returns (ok_fatal, issues)."""
        issues: List[Tuple[str, str]] = []
        ok_fatal = True

        if self.time_s is None or len(self.time_s) == 0:
            issues.append(("FATAL", f"Vial {self.number}: missing/empty time_s."))
            ok_fatal = False

        if self.retentate_signal is None or len(self.retentate_signal) == 0:
            issues.append(("FATAL", f"Vial {self.number}: missing/empty retentate_signal."))
            ok_fatal = False

        if (self.time_s is not None) and (self.retentate_signal is not None):
            if len(self.time_s) != len(self.retentate_signal):
                issues.append(("FATAL", f"Vial {self.number}: time_s and retentate_signal lengths differ."))
                ok_fatal = False

        if self.time_s is not None and len(self.time_s) > 1:
            if np.any(np.diff(self.time_s) < 0):
                issues.append(("FATAL", f"Vial {self.number}: time_s is not monotonic non-decreasing."))
                ok_fatal = False

        return ok_fatal, issues


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

    # -------------------------
    # Loader provenance flags (honesty about assumptions)
    # -------------------------
    # If a loader fills missing metadata using defaults, it records which fields
    # were default-filled here so validate_experiment() can warn downstream.
    used_defaults: List[str] = field(default_factory=list, repr=False)

    def add_vial(self, vial: VialData) -> None:
        """Append one vial while preserving explicit ordering."""
        self.vials.append(vial)

    def get_vial_switch_times(self) -> Tuple[np.ndarray, np.ndarray, float]:
        """Return shifted per-vial start/end times and common delay shift."""
        if not self.vials:
            raise ValueError("ExperimentalData.vials is empty; cannot compute vial switch times.")
        first = self.vials[0]
        if first.time_s is None or len(first.time_s) == 0:
            raise ValueError("Vial 1 has missing/empty time_s; cannot compute t_delay_s.")

        t_delay_s = float(first.time_s[0])
        ti_s = np.array([float(v.time_s[0]) - t_delay_s for v in self.vials if v.time_s is not None and len(v.time_s) > 0], dtype=float)
        tf_s = np.array([float(v.time_s[-1]) - t_delay_s for v in self.vials if v.time_s is not None and len(v.time_s) > 0], dtype=float)
        return ti_s, tf_s, t_delay_s

    def get_inputs(self) -> Dict[str, Optional[float]]:
        """Return experiment-level operating conditions for model construction."""
        return {
            "delP_bar": self.delP_bar,
            "Temp_K": self.Temp_K,
            "Am_cm2": self.Am_cm2,
            "rho_g_cm3": self.rho_g_cm3,
        }

    def get_measurements(self) -> List[Dict[str, Optional[np.ndarray]]]:
        """Return per-vial measurement arrays in a model-friendly structure."""
        out: List[Dict[str, Optional[np.ndarray]]] = []
        for v in self.vials:
            out.append({
                "time_s": v.time_s,
                "retentate_signal": v.retentate_signal,
                "permeate_signal": v.permeate_signal,
                "mass_g": v.mass_g,
                "retentate_concentration": v.retentate_concentration,
                "permeate_concentration": v.permeate_concentration,
            })
        return out

    def validate(self) -> Tuple[bool, List[Tuple[str, str]]]:
        """Validate structural integrity and required metadata; returns (ok_fatal, issues)."""
        issues: List[Tuple[str, str]] = []
        ok_fatal = True

        if self.vials is None or len(self.vials) == 0:
            issues.append(("FATAL", "No vials found (exp.vials is empty)."))
            return False, issues

        for v in self.vials:
            vial_ok, vial_issues = v.validate_structure()
            issues.extend(vial_issues)
            if not vial_ok:
                ok_fatal = False

        if is_missing_scalar(self.delP_bar):
            issues.append(("REQUIRED", "Missing exp.delP_bar (pressure)."))
        if is_missing_scalar(self.Temp_K):
            issues.append(("REQUIRED", "Missing exp.Temp_K (temperature)."))
        if is_missing_scalar(self.Am_cm2):
            issues.append(("REQUIRED", "Missing exp.Am_cm2 (membrane area)."))

        if self.theta0 is not None:
            try:
                _ = np.asarray(self.theta0, dtype=float)
            except Exception:
                issues.append(("REQUIRED", "exp.theta0 exists but is not numeric/castable to float array."))

        for fname in self.used_defaults or []:
            issues.append(("USED_DEFAULT", f"exp.{fname} was filled from MAT_GLOBAL_DEFAULTS (not provided by source)."))

        return ok_fatal, issues


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

    # Backward-compatible wrapper around the class method.
    return exp.validate()


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


def load_experiment(path: str, selector: object, *, initial_guess_db: Optional[Dict[Tuple[str, ...], Dict[str, object]]] = None) -> Tuple[ExperimentalData, Tuple[bool, List[Tuple[str, str]]]]:
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

        initial_guess_db:
            Optional mapping of initial guesses used to populate exp.Lp0/B0/sigma0/theta0
            when loading from XLSX (because spreadsheets typically do not store these).
            - Keys: tuple(component_names), e.g., ('KCl',) or ('NaCl','MgSO4')
            - Values: dict with keys like 'Lp0','B0','sigma0','theta0'
            If None, XLSX-loaded experiments will leave initial-guess fields as None.

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
        exp = load_from_xlsx(file_path, sheet_selector=selector, initial_guess_db=initial_guess_db)

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

# =============================================================================
# Convenience API (thin wrappers around the core loader)
# =============================================================================

def load_experiment_easy(
    path: str,
    selector: Optional[object] = None,
    *,
    specs: Optional[Dict[str, object]] = None,
    initial_guess_db: Optional[Dict[Tuple[str, ...], Dict[str, object]]] = None,
    # ---------------------------------------------------------------------
    # Conductivity -> concentration conversion (XLSX only)
    # ---------------------------------------------------------------------
    convert_to_concentration: bool = False,
    concentration_units: str = "mM",
    conductivity_model: str = "auto",
    conductivity_model_params: Optional[Dict[str, object]] = None,
    n_cations: Optional[int] = None,
    temp_K: Optional[float] = None,
    # ---------------------------------------------------------------------
    # Optional quick plotting
    # ---------------------------------------------------------------------
    plot: bool = False,
    plot_kind: str = "retentate_signal",
) -> Tuple[ExperimentalData, Tuple[bool, List[Tuple[str, str]]]]:
    """Load ONE experiment with optional user overrides and an optional quick plot.

    Why this exists:
        The core loader (load_experiment + load_from_xlsx/load_from_mat) is written to be
        explicit and debuggable. In day-to-day use, you often want a *one-liner*:
            - point at a file (xlsx/mat),
            - optionally tell it a few experiment-level specs/overrides (pressure, mode, etc.),
            - optionally see a plot right away.

    Inputs:
        path:
            File path to either an Excel workbook (.xlsx/.xls) or a MATLAB file (.mat).

        selector:
            Source-specific experiment selector.
                - XLSX: sheet name (str) or sheet index (int).
                - MAT:  struct key (str). If None, defaults to "data_stru".

            If selector is None:
                - XLSX: use the first sheet (index 0).
                - MAT:  use "data_stru".

        specs:
            Optional overrides applied *after* loading. Keys must match fields on ExperimentalData.
            Example:
                specs={"mode": "Lag", "delP_bar": 10.0, "Temp_K": 298.15}

            NOTE:
                This is intentionally a plain dict (not a new class) to keep the interface light.

        initial_guess_db:
            Optional lookup table used by loaders to populate exp.(Lp0,B0,sigma0,theta0,...) when
            the source does not provide guesses (common for XLSX).

        plot:
            If True, produce a quick matplotlib plot (see plot_kind).

        plot_kind:
            What to plot. Supported:
                - "retentate_signal" (default)
                - "permeate_signal"
                - "mass_g"

    Outputs:
        exp:
            ExperimentalData populated from the file (plus any applied overrides in specs).

        validation:
            (ok_fatal, issues) from validate_experiment(exp).

    """
    # -------------------------------------------------------------------------
    # 0) Auto-default selector for convenience
    # -------------------------------------------------------------------------
    file_path = Path(path)

    kind = detect_source_type(file_path)
    if kind is None:
        raise ValueError(f"Unsupported file type: {file_path.suffix}")

    # Default selection:
    #   - XLSX: first sheet
    #   - MAT:  data_stru
    if selector is None:
        selector = 0 if kind == SourceType.XLSX else "data_stru"

    # -------------------------------------------------------------------------
    # 1) Load experiment using the core, file-agnostic entrypoint
    # -------------------------------------------------------------------------
    exp, validation = load_experiment(path, selector, initial_guess_db=initial_guess_db)

    # -------------------------------------------------------------------------
    # 2) Apply optional user overrides (specs) *by setting class fields directly*
    # -------------------------------------------------------------------------
    # This keeps the overall framework file-type agnostic:
    # regardless of XLSX/MAT source, you can enforce a consistent exp field value.
    if specs:
        for k, v in specs.items():
            if hasattr(exp, k):
                setattr(exp, k, v)
                # If user overrides a value previously default-filled by the XLSX loader,
                # clear the corresponding USED_DEFAULT marker from both provenance and issues.
                if k in getattr(exp, "used_defaults", []):
                    exp.used_defaults = [name for name in exp.used_defaults if name != k]
                    ok_fatal, issues = validation
                    issues = [
                        (sev, msg)
                        for (sev, msg) in issues
                        if not (sev == "USED_DEFAULT" and f"exp.{k} " in msg)
                    ]
                    validation = (ok_fatal, issues)
            else:
                # Don't crash; just flag that the override key was ignored.
                ok_fatal, issues = validation
                issues.append(("WARNING", f"specs override ignored: ExperimentalData has no field '{k}'."))
                validation = (ok_fatal, issues)

    # -------------------------------------------------------------------------
    # 3) Optional conductivity -> concentration conversion (XLSX only)
    # -------------------------------------------------------------------------
    # Important design rule (per our pseudocode + your requirement):
    #   - We store conductivity *honestly* as the measured signal (uS/cm).
    #   - We convert to concentration only when explicitly requested, and only
    #     for XLSX sources (MAT already contains concentration from legacy workflows).
    if convert_to_concentration:
        try:
            if exp.source == SourceType.XLSX:
                # If the user explicitly provides number of cations, honor the rule:
                #   1-2 cations -> variant Shedlovsky
                #   >=3 cations -> MSA
                # This avoids relying on salt-name parsing when the user already knows.
                chosen_model = conductivity_model
                if str(conductivity_model).lower() == "auto" and n_cations is not None:
                    chosen_model = "variant_shedlovsky" if int(n_cations) <= 2 else "msa"

                apply_conductivity_to_concentration(
                    exp,
                    output_units=concentration_units,
                    model=chosen_model,
                    model_params=conductivity_model_params,
                    temp_K=temp_K,
                )
            else:
                # MAT inputs: should already have concentrations; no conversion needed.
                ok_fatal, issues = validation
                issues.append((
                    "INFO",
                    "convert_to_concentration=True ignored for MAT input (MAT is expected to already contain concentrations).",
                ))
                validation = (ok_fatal, issues)
        except Exception as e:
            # Conversion failures should not crash loading; they should surface as warnings.
            ok_fatal, issues = validation
            issues.append(("WARNING", f"Conductivity->concentration conversion failed: {e!r}"))
            validation = (ok_fatal, issues)

    # -------------------------------------------------------------------------
    # 3) Optional quick look plot (for sanity checks during development)
    # -------------------------------------------------------------------------
    if plot:
        try:
            plot_experiment(exp, kind=plot_kind)
        except Exception as e:
            # Plotting should never prevent loading. Flag and continue.
            ok_fatal, issues = validation
            issues.append(("WARNING", f"Plot failed ({plot_kind}): {e!r}"))  # keep it honest + debuggable
            validation = (ok_fatal, issues)

    return exp, validation

# =============================================================================
# Conductivity -> concentration conversion wrapper (uses conductivity_paper.py as-is)
# =============================================================================

def _count_cations_from_salts(component_names: Optional[List[str]]) -> int:
    """
    Best-effort estimate of the number of distinct cations present.

    Inputs:
        component_names:
            List of salts/components, e.g., ["NaCl", "CaCl2", "MgSO4"].

    Output:
        int:
            Count of distinct cation element symbols (heuristic).
            Used ONLY to select between Shedlovsky vs MSA per your rule.

    Notes:
        This is intentionally simple: it extracts the leading element symbol.
        If you later want rigorous parsing, we can swap this with a proper formula parser.
    """
    if not component_names:
        return 0

    cations = set()
    for name in component_names:
        if not isinstance(name, str):
            continue
        m = re.match(r"^([A-Z][a-z]?)", name.strip())
        cations.add(m.group(1) if m else name.strip())

    return len(cations)


# -----------------------------------------------------------------------------
# Conductivity -> concentration: variant Shedlovsky parameter auto-fill helpers
# -----------------------------------------------------------------------------
# Motivation:
#   conductivity_paper.py expects certain electrolyte-specific parameters (e.g., z_1, z_2, a, lambda_0).
#   For XLSX experiments, users often know the salts but not these numeric constants.
#   We therefore provide a small built-in lookup for common ions + a basic formula parser,
#   and we *only* fill in parameters that are missing (user-provided values always win).
#
# IMPORTANT:
#   - This does NOT modify conductivity_paper.py.
#   - This is a convenience layer for XLSX workflows.
#   - If the salt formula is unknown, we raise a ValueError with a clear message.

# Ionic charges (z) used by both Shedlovsky and MSA wrappers.
# Source: common electrochemistry data tables (25 C in water).
ION_CHARGE: Dict[str, int] = {
    "H": +1,
    "Li": +1,
    "Na": +1,
    "K": +1,
    "NH4": +1,
    "Mg": +2,
    "Ca": +2,
    "Co": +2,
    "La": +3,
    "OH": -1,
    "F": -1,
    "Cl": -1,
    "Br": -1,
    "I": -1,
    "NO3": -1,
    "ClO4": -1,
    "SO4": -2,
    "Acetate": -1,
}

# Limiting ionic equivalent conductivities at infinite dilution at 25 C in water.
# Units: S*cm^2/equiv (equivalent to S*cm^2/mol for singly charged ions).
# Values transcribed from your provided Appendix 6.1/6.2 tables.
ION_LAMBDA0_25C_S_CM2_EQUIV: Dict[str, float] = {
    "H": 349.8,
    "Li": 38.6,
    "Na": 50.10,
    "K": 73.50,
    "NH4": 73.5,
    "Mg": 53.06,
    "Ca": 59.50,
    "Co": 55.0,
    "La": 69.7,
    "OH": 199.1,
    "F": 55.4,
    "Cl": 76.35,
    "Br": 78.14,
    "I": 76.84,
    "NO3": 71.46,
    "ClO4": 67.3,
    "SO4": 80.0,
    "Acetate": 40.9,
}

# Backward-compatible alias used elsewhere in this file.
ION_LAMBDA0_S_CM2_MOL = ION_LAMBDA0_25C_S_CM2_EQUIV

def _build_salt_lambda_params_25c(*, cation: str, nu_c: int, anion: str, nu_a: int) -> Dict[str, float]:
    """Build {'lambda_0','lambda_0_cation','lambda_0_anion'} for one binary salt at 25 C."""
    if cation not in ION_LAMBDA0_25C_S_CM2_EQUIV or anion not in ION_LAMBDA0_25C_S_CM2_EQUIV:
        raise ValueError(f"Missing ionic lambda_0 values for '{cation}' or '{anion}'.")
    z_cat = ION_CHARGE.get(cation)
    if z_cat is None or z_cat <= 0:
        raise ValueError(f"Missing/invalid cation charge for '{cation}'.")

    lambda_0_cation = float(ION_LAMBDA0_25C_S_CM2_EQUIV[cation])
    lambda_0_anion = float(ION_LAMBDA0_25C_S_CM2_EQUIV[anion])
    eq_per_mol = float(nu_c * abs(z_cat))
    if eq_per_mol <= 0:
        raise ValueError("Invalid equivalents-per-mole while building salt lambda values.")
    lambda_0 = (nu_c * lambda_0_cation + nu_a * lambda_0_anion) / eq_per_mol
    return {
        "lambda_0": float(lambda_0),
        "lambda_0_cation": float(lambda_0_cation),
        "lambda_0_anion": float(lambda_0_anion),
    }

# Ready-to-use salt dictionary for common electrolytes in this project.
SALT_LAMBDA_PARAMS_25C: Dict[str, Dict[str, float]] = {
    "NaCl": _build_salt_lambda_params_25c(cation="Na", nu_c=1, anion="Cl", nu_a=1),
    "KCl": _build_salt_lambda_params_25c(cation="K", nu_c=1, anion="Cl", nu_a=1),
    "LiCl": _build_salt_lambda_params_25c(cation="Li", nu_c=1, anion="Cl", nu_a=1),
    "NH4Cl": _build_salt_lambda_params_25c(cation="NH4", nu_c=1, anion="Cl", nu_a=1),
    "NaNO3": _build_salt_lambda_params_25c(cation="Na", nu_c=1, anion="NO3", nu_a=1),
    "KNO3": _build_salt_lambda_params_25c(cation="K", nu_c=1, anion="NO3", nu_a=1),
    "NaClO4": _build_salt_lambda_params_25c(cation="Na", nu_c=1, anion="ClO4", nu_a=1),
    "MgCl2": _build_salt_lambda_params_25c(cation="Mg", nu_c=1, anion="Cl", nu_a=2),
    "CaCl2": _build_salt_lambda_params_25c(cation="Ca", nu_c=1, anion="Cl", nu_a=2),
    "CoCl2": _build_salt_lambda_params_25c(cation="Co", nu_c=1, anion="Cl", nu_a=2),
    "Na2SO4": _build_salt_lambda_params_25c(cation="Na", nu_c=2, anion="SO4", nu_a=1),
    "MgSO4": _build_salt_lambda_params_25c(cation="Mg", nu_c=1, anion="SO4", nu_a=1),
    "LaCl3": _build_salt_lambda_params_25c(cation="La", nu_c=1, anion="Cl", nu_a=3),
}

def _get_salt_lambda_params_25c(salt_name: str) -> Dict[str, float]:
    """Return lambda dictionary for a salt using explicit dict first, formula fallback second."""
    s = str(salt_name).strip()
    if s in SALT_LAMBDA_PARAMS_25C:
        return dict(SALT_LAMBDA_PARAMS_25C[s])
    cat, nu_cat, an, nu_an = _parse_simple_salt_formula(s)
    return _build_salt_lambda_params_25c(cation=cat, nu_c=nu_cat, anion=an, nu_a=nu_an)

# # Default closest-approach distance for Shedlovsky-type correlations.
# # Units: cm  (4 Å = 4e-8 cm)
# DEFAULT_A_CM = 4e-8
# Conservative default closest-approach distance used by the conductivity paper model [cm]
DEFAULT_A_CM: float = 1e-8


def _parse_simple_salt_formula(salt: str) -> tuple[str, int, str, int]:
    """Parse a simple binary salt formula like 'NaCl' or 'MgCl2'.

    Inputs:
        salt:
            String chemical formula (no spaces). Supported examples:
                - 'NaCl'   -> (Na, 1, Cl, 1)
                - 'MgCl2'  -> (Mg, 1, Cl, 2)
                - 'Na2SO4' -> (Na, 2, SO4, 1)

    Output:
        (cation, nu_c, anion, nu_a):
            cation: cation symbol (e.g., 'Na')
            nu_c:   cation stoichiometric coefficient (e.g., 2 for Na2SO4)
            anion:  anion symbol (e.g., 'SO4')
            nu_a:   anion stoichiometric coefficient (e.g., 1 for Na2SO4)

    Raises:
        ValueError:
            If the formula cannot be parsed by this simple parser.
    """
    import re

    s = (salt or "").strip()
    if not s:
        raise ValueError("Empty salt formula.")

    # Very small parser:
    #   - Find a leading cation token (letters/numbers), then the rest is anion token.
    #   - Supports multi-character ions like 'SO4' or 'NH4' if provided as-is.

    # Sort known ions by length so we match 'NH4' before 'N', 'SO4' before 'S', etc.
    known_ions = sorted(ION_CHARGE.keys(), key=len, reverse=True)

    cation = None
    for ion in known_ions:
        if s.startswith(ion) and ION_CHARGE[ion] > 0:
            cation = ion
            break
    if cation is None:
        raise ValueError(f"Could not identify cation in salt '{salt}'.")

    rest = s[len(cation):]
    m = re.match(r"^(\d+)?(.*)$", rest)
    if not m:
        raise ValueError(f"Could not parse stoichiometry after cation in salt '{salt}'.")
    nu_c_str, rest2 = m.group(1), m.group(2)
    nu_c = int(nu_c_str) if nu_c_str else 1

    anion = None
    for ion in known_ions:
        if rest2.startswith(ion) and ION_CHARGE[ion] < 0:
            anion = ion
            break
    if anion is None:
        raise ValueError(f"Could not identify anion in salt '{salt}'.")

    rest3 = rest2[len(anion):]
    m2 = re.match(r"^(\d+)?$", rest3)
    if not m2:
        raise ValueError(f"Could not parse anion stoichiometry in salt '{salt}'.")
    nu_a_str = m2.group(1)
    nu_a = int(nu_a_str) if nu_a_str else 1

    return cation, nu_c, anion, nu_a


def _autofill_variant_shedlovsky_params(model_params: dict, *, salt_name: str) -> dict:
    """Fill missing variant-Shedlovsky parameters in-place (returns the same dict)."""
    # Parse the formula (cation/anions + stoichiometry).
    cat, nu_cat, an, nu_an = _parse_simple_salt_formula(salt_name)
    salt_lambda = _get_salt_lambda_params_25c(salt_name)

    # Charges for the electrolyte (z_1, z_2 in conductivity_paper naming).
    z_cat = ION_CHARGE.get(cat, None)
    z_an = ION_CHARGE.get(an, None)
    if z_cat is None or z_an is None:
        raise ValueError(f"Missing ion charge for '{cat}' or '{an}'. Add it to ION_CHARGE.")

    # 1) Fill z_1, z_2 if missing (treat None as missing too)
    if model_params.get("z_1") is None:
        model_params["z_1"] = float(z_cat)  # cation charge
    if model_params.get("z_2") is None:
        model_params["z_2"] = float(z_an)   # anion charge (negative)

    # 2) Fill ionic limiting conductivities if missing (or explicitly None).
    if model_params.get("lambda_0_cation") is None:
        model_params["lambda_0_cation"] = float(salt_lambda["lambda_0_cation"])

    if model_params.get("lambda_0_anion") is None:
        model_params["lambda_0_anion"] = float(salt_lambda["lambda_0_anion"])

    # 3) Fill lambda_0 (electrolyte limiting equivalent conductivity) if missing.
    eq_per_mol = float(nu_cat * abs(z_cat))
    if eq_per_mol <= 0:
        raise ValueError(f"Invalid equivalents per mole for salt '{salt_name}'.")

    lam0 = (nu_cat * float(model_params["lambda_0_cation"]) + nu_an * float(model_params["lambda_0_anion"])) / eq_per_mol
    # 3) Fill overall limiting molar conductivity (lambda_0) only if missing.
    # Use stoichiometry-aware inferred value from ionic limiting conductivities.
    if model_params.get("lambda_0") is None:
        model_params["lambda_0"] = float(salt_lambda.get("lambda_0", lam0))

    # 4) Fill a (closest-approach distance) if missing
    if model_params.get("a") is None:
        model_params["a"] = float(DEFAULT_A_CM)

    return model_params

def _autofill_msa_params(model_params: dict, *, salt_names: Optional[List[str]]) -> dict:
    """
    Fill MSA conductivity parameters from salt names when possible.

    Focus for this auto-fill:
        - lambda_0 (ionic list in S*m^2/mol)
        - valency (if missing)
        - n_salts (if missing)
    """
    mp = dict(model_params)
    if not salt_names:
        return mp

    parsed: List[Tuple[str, int, str, int]] = [_parse_simple_salt_formula(s) for s in salt_names]
    n_salts = int(mp.get("n_salts", len(parsed)))
    mp.setdefault("n_salts", n_salts)

    # This wrapper currently supports mixtures with a shared anion in the underlying MSA call.
    common_anion = parsed[0][2]
    if any(an != common_anion for _, _, an, _ in parsed):
        return mp

    if mp.get("valency") is None:
        valency = [int(ION_CHARGE[cat]) for cat, _, _, _ in parsed] + [int(ION_CHARGE[common_anion])]
        mp["valency"] = valency

    if mp.get("lambda_0") is None:
        # conductivity_paper.msa_transport expects ionic lambda_0 in S*m^2/mol
        # ordered as [cat1, cat2, ..., common_anion].
        lambda_ion_s_m2_mol = [float(ION_LAMBDA0_25C_S_CM2_EQUIV[cat]) * 1e-4 for cat, _, _, _ in parsed]
        lambda_ion_s_m2_mol.append(float(ION_LAMBDA0_25C_S_CM2_EQUIV[common_anion]) * 1e-4)
        mp["lambda_0"] = lambda_ion_s_m2_mol

    return mp


def _invert_monotone_1d(
    *,
    fwd_func,
    y_target: float,
    x_lo: float,
    x_hi: float,
    tol: float = 1e-6,
    max_iter: int = 200,
) -> float:
    """
    Robust 1D inversion helper: find x such that fwd_func(x) ~= y_target.

    Despite the name, this routine does NOT assume strict monotonicity.

    Strategy:
        - Define g(x) = fwd_func(x) - y_target
        - Search for a sign-change bracket [a,b] by scanning a grid over [x_lo, x_hi]
        - If no bracket found, expand x_hi progressively and try again
        - Once bracketed, solve with Brent (if SciPy available), else bisection fallback
    """
    import math

    def g(x: float) -> float:
        """G.

        Args:
            x: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        return float(fwd_func(float(x)) - y_target)

    # Progressive expansions of the upper bound (kept conservative)
    expand_factors = [1.0, 2.0, 5.0, 10.0]
    last_report = None

    for fac in expand_factors:
        lo = float(x_lo)
        hi = float(x_hi) * fac

        # Mixed linear/log grid to find brackets robustly
        eps = max(1e-12, lo + 1e-12)
        log_pts = [eps * (hi / eps) ** (i / 59) for i in range(60)]
        lin_pts = [lo + (hi - lo) * (i / 59) for i in range(60)]
        grid = sorted(set([*log_pts, *lin_pts, hi]))

        vals = []
        for x in grid:
            try:
                gx = g(x)
            except Exception:
                continue
            if math.isnan(gx) or math.isinf(gx):
                continue
            vals.append((x, gx))

        if len(vals) < 2:
            last_report = f"Too few valid forward evaluations on [{lo}, {hi}]."
            continue

        # Find first sign change
        for (x1, g1), (x2, g2) in zip(vals[:-1], vals[1:]):
            if g1 == 0.0:
                return float(x1)
            if g2 == 0.0:
                return float(x2)
            if (g1 < 0 < g2) or (g1 > 0 > g2):
                a, b = float(x1), float(x2)

                # Root solve within bracket
                try:
                    from scipy.optimize import brentq  # type: ignore
                    return float(brentq(g, a, b, xtol=tol, maxiter=max_iter))
                except Exception:
                    # Bisection fallback
                    left, right = a, b
                    gl, gr = g(left), g(right)
                    for _ in range(max_iter):
                        mid = 0.5 * (left + right)
                        gm = g(mid)
                        if abs(right - left) <= tol:
                            return float(mid)
                        if (gl < 0 < gm) or (gl > 0 > gm):
                            right, gr = mid, gm
                        else:
                            left, gl = mid, gm
                    return float(0.5 * (left + right))

        g_values = [gv for _, gv in vals]
        last_report = (
            f"No bracket found on [{lo}, {hi}]. "
            f"g_min={min(g_values):.6g}, g_max={max(g_values):.6g} (target y={y_target})."
        )

    raise ValueError(
        "Target y not bracketed in conductivity inversion. "
        + (last_report or "")
        + " Check units/parameters (e.g., lambda_0, a, z) or widen bounds."
    )


def apply_conductivity_to_concentration(
    exp: ExperimentalData,
    *,
    output_units: str = "mM",
    model: str = "auto",
    model_params: Optional[Dict[str, object]] = None,
    temp_K: Optional[float] = None,
) -> None:
    """
    Single clean wrapper: converts BOTH retentate and permeate conductivity to concentration.

    Inputs:
        exp:
            Loaded ExperimentalData (from XLSX or MAT).
        output_units:
            "mM" or "M" (stored into vial fields).
        model:
            - "auto": choose variant Shedlovsky if 1–2 cations, else MSA.
            - or force "variant_shedlovsky" / "msa".
        model_params:
            Dict of required parameters for the chosen model.
        temp_K:
            Override temperature [K]. If None, uses exp.Temp_K.

    Output:
        None. Mutates exp.vials in-place.

    What gets filled:
        v.retentate_concentration (+ units + model)
        v.permeate_concentration  (+ units + model)
    """
    if model_params is None:
        model_params = {}

    # ---------------------------------------------------------------------
    # Best-effort salt identifier for conductivity models
    # ---------------------------------------------------------------------
    # The XLSX loader attempts to populate exp.component_names from the sheet (e.g., cell E1).
    # Variant Shedlovsky requires salt identity to infer (lambda_0, a, z_1, z_2) when the user
    # does not provide them explicitly.
    primary_salt: Optional[str] = None
    if exp.component_names and len(exp.component_names) > 0:
        primary_salt = str(exp.component_names[0]).strip()  # use the first salt as the "primary" one

    # Choose temperature
    T_K = float(temp_K) if temp_K is not None else (float(exp.Temp_K) if exp.Temp_K is not None else None)
    if T_K is None:
        raise ValueError("Temperature required for conversion (provide temp_K or exp.Temp_K).")

    # Choose model per your rule
    if model.lower() == "auto":
        n_cations = _count_cations_from_salts(exp.component_names)
        chosen = "variant_shedlovsky" if n_cations <= 2 else "msa"
    else:
        chosen = model

    # Prepare model parameters once; per-vial conversion reuses these.
    model_params_prepared = dict(model_params)
    if chosen.lower() in ["msa", "mean_spherical_approximation"]:
        model_params_prepared = _autofill_msa_params(model_params_prepared, salt_names=exp.component_names)

    # One internal helper to avoid permeate/retentate duplication
    def _convert_signal_into_vial(v: VialData, *, which: str) -> None:
        # Select the correct signal + units
        """Convert signal into vial.

        Args:
            v: Parameter description.
            which: Parameter description.

        """
        if which == "retentate":
            sig = v.retentate_signal
            sig_units = v.retentate_signal_units
        elif which == "permeate":
            sig = v.permeate_signal
            sig_units = v.permeate_signal_units
        else:
            raise ValueError("which must be 'retentate' or 'permeate'")

        # No signal => nothing to do
        if sig is None:
            return

        # Safety: this wrapper expects conductivity in uS/cm (per your “store honestly” rule)
        if sig_units is None or "us/cm" not in str(sig_units).lower():
            raise ValueError(
                f"Expected {which} conductivity units of uS/cm, got '{sig_units}'. "
                "Do not pre-convert; store conductivity honestly and convert here."
            )

        conc = _conductivity_to_concentration_series(
            cond_uS_cm=sig,
            temp_K=T_K,
            model=chosen,
            model_params=model_params_prepared,
            output_units=output_units,
            salt_name=primary_salt,  # enables Shedlovsky parameter inference when missing
        )

        # Store into the correct vial fields
        if which == "retentate":
            v.retentate_concentration = conc
            v.retentate_concentration_units = output_units
            v.retentate_concentration_model = chosen
        else:
            v.permeate_concentration = conc
            v.permeate_concentration_units = output_units
            v.permeate_concentration_model = chosen

    # Apply across all vials, both signals
    for v in exp.vials:
        _convert_signal_into_vial(v, which="retentate")
        _convert_signal_into_vial(v, which="permeate")

def plot_experiment(exp: ExperimentalData, *, kind: str = "retentate_signal") -> None:
    """Quick visualization helper for a loaded experiment (debugging / sanity check).

    Inputs:
        exp:
            ExperimentalData instance (already loaded).

        kind:
            What to plot. Accepted values (case-insensitive):
                - "retentate_signal": plots each vial's retentate_signal vs time_s
                - "permeate_signal": plots each vial's permeate_signal vs time_s (if available)
                - "mass_g":          plots each vial's mass_g vs time_s (if available)

            Common aliases (accepted for convenience):
                - "mass", "weight"   -> "mass_g"
                - "retentate"        -> "retentate_signal"
                - "permeate"         -> "permeate_signal"
                - "signals"          -> plot retentate + permeate on the same axes (if permeate exists)

    Output:
        None. Displays a matplotlib figure.

    Notes:
        - This helper is intentionally simple: it does not stitch vials into a single timeline.
        - For publication-quality plotting, build a dedicated plotting module later.
    """
    # Local import keeps the core loader usable without matplotlib installed.
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    if exp.vials is None or len(exp.vials) == 0:
        raise ValueError("Experiment has no vials to plot.")

    # -------------------------------------------------------------------------
    # Normalize 'kind' to avoid brittle user-facing strings.
    # -------------------------------------------------------------------------
    if kind is None:
        kind = "retentate_signal"

    kind_norm = str(kind).strip().lower()

    alias_map = {
        # Mass aliases
        "mass": "mass_g",
        "weight": "mass_g",
        "mass_g": "mass_g",

        # Signal aliases
        "retentate": "retentate_signal",
        "ret": "retentate_signal",
        "retentate_signal": "retentate_signal",
        "permeate": "permeate_signal",
        "perm": "permeate_signal",
        "permeate_signal": "permeate_signal",

        # Convenience composite
        "signals": "signals",
        "signal": "signals",
        "conductivity": "signals",
        "cond": "signals",
    }
    kind_norm = alias_map.get(kind_norm, kind_norm)

    # -------------------------------------------------------------------------
    # Plotting logic
    # -------------------------------------------------------------------------
    plt.figure()

    def _show_or_close() -> None:
        # In headless smoke tests (Agg backend), calling plt.show() emits a warning.
        """Show or close.

        """
        if "agg" in str(plt.get_backend()).strip().lower():
            plt.close()
        else:
            plt.show()

    if kind_norm == "signals":
        # Goal:
        #   - color changes when the vial swaps (one color per vial)
        #   - retentate uses dashed lines, permeate uses solid lines
        #   - TWO legends:
        #       (1) Vial number -> color
        #       (2) Signal type -> line style

        # Matplotlib's default color cycle (keeps the plot readable without us hard-coding colors).
        colors = plt.rcParams.get("axes.prop_cycle", None)
        colors = (colors.by_key().get("color", []) if colors is not None else [])

        # NOTE: We build legends manually (below), so avoid adding per-line labels.

        for i, v in enumerate(exp.vials):
            if v.time_s is None:
                continue

            # Same color for retentate/permeate within a vial segment.
            c = (colors[i % len(colors)] if colors else None)

            # Retentate signal (dashed)
            if v.retentate_signal is not None:
                plt.plot(
                    v.time_s,
                    v.retentate_signal,
                    linestyle="--",
                    color=c,
                    label="_nolegend_",
                )

            # Permeate signal (solid)
            if v.permeate_signal is not None:
                plt.plot(
                    v.time_s,
                    v.permeate_signal,
                    linestyle="-",
                    color=c,
                    label="_nolegend_",
                )

        plt.xlabel("time (s)")
        plt.ylabel("signal (raw units)")
        plt.title(f"{exp.filename or 'experiment'} — signals")

        # -----------------------------
        # Legend 1: signal type (style)
        # -----------------------------
        ax = plt.gca()
        style_handles = [
            Line2D([0], [0], color="black", linestyle="--", label="Retentate"),
            Line2D([0], [0], color="black", linestyle="-", label="Permeate"),
        ]
        style_leg = ax.legend(
            handles=style_handles,
            title="Signal",
            loc="upper left",
        )
        ax.add_artist(style_leg)

        # --------------------------
        # Legend 2: vial id (color)
        # --------------------------
        vial_handles = []
        for i in range(len(exp.vials)):
            c = (colors[i % len(colors)] if colors else "black")
            vial_handles.append(Line2D([0], [0], color=c, linestyle="-", label=f"Vial {i+1}"))
        ax.legend(
            handles=vial_handles,
            title="Vial",
            loc="upper right",
        )

        plt.tight_layout()
        _show_or_close()
        return

    # Single-attribute plots (retentate/permeate/mass)
    attr = {
        "retentate_signal": "retentate_signal",
        "permeate_signal": "permeate_signal",
        "mass_g": "mass_g",
    }.get(kind_norm, None)

    if attr is None:
        raise ValueError(f"Unknown plot kind: {kind!r}.")

    # One line per vial segment (honest segmentation)
    #   - color = vial index
    #   - linestyle = retentate/permeate convention (where applicable)
    colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not colors:
        colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]

    ax = plt.gca()
    vial_handles: list[Line2D] = []

    for i, v in enumerate(exp.vials):
        if v.time_s is None:
            continue

        y = getattr(v, attr)
        if y is None:
            continue

        c = colors[i % len(colors)]

        # Style convention requested:
        #   - retentate = dashed
        #   - permeate = solid
        #   - mass = solid
        if kind_norm == "retentate_signal":
            ls = "--"
        else:
            ls = "-"

        plt.plot(v.time_s, y, linestyle=ls, color=c, label="_nolegend_")
        vial_handles.append(Line2D([0], [0], color=c, linestyle="-", label=f"Vial {v.number}"))

    plt.xlabel("time (s)")
    plt.ylabel(kind_norm)
    plt.title(f"{exp.filename or 'experiment'} — {kind_norm}")

    # Two legends (requested):
    #   1) style legend (what signal type this plot is)
    #   2) vial legend (colors map to vial number)
    if kind_norm in ("retentate_signal", "permeate_signal"):
        style_label = "Retentate" if kind_norm == "retentate_signal" else "Permeate"
        style_ls = "--" if kind_norm == "retentate_signal" else "-"
        style_leg = ax.legend(
            handles=[Line2D([0], [0], color="black", linestyle=style_ls, label=style_label)],
            title="Signal",
            loc="upper left",
        )
        ax.add_artist(style_leg)

    ax.legend(handles=vial_handles, title="Vial", loc="upper right")

    plt.tight_layout()
    _show_or_close()


def load_from_xlsx(xlsx_path: Path, sheet_selector: object, *, initial_guess_db: Optional[Dict[Tuple[str, ...], Dict[str, object]]] = None) -> ExperimentalData:
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
        - It may fill missing operating-condition metadata using MAT_GLOBAL_DEFAULTS (and records this in exp.used_defaults).

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
    # 6) Parse additional Excel blocks (metadata + calibration)
    # -------------------------------------------------------------------------
    # Columns J/K: key/value metadata (experiment name, initial/final weights, salts, etc.)
    metadata_kv = parse_experiment_metadata_kv(sheet)
    # Columns O..V: ICP-OES assay results (kept as raw table for now)
    _icp_assay_df = parse_icp_assay_table(sheet)
    # Columns X..AA: calibration points for conductivity->concentration conversion (raw)
    _calib_df = parse_calibration_block(sheet)


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

    
    # -------------------------------------------------------------------------
    # 7b) Map Excel metadata into canonical ExperimentalData fields (best-effort)
    # -------------------------------------------------------------------------
    # We keep this conservative: only populate fields we can infer reliably.
    # Any workbook-specific keys you adopt later can be enforced in validate_experiment(...).

    def _meta_lookup(*candidates: str):
        # Return the first matching key (case-insensitive) from metadata_kv.
        """Meta lookup.

        Args:
            *candidates: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        if not metadata_kv:
            return None
        lowered = {str(k).strip().lower(): k for k in metadata_kv.keys()}
        for c in candidates:
            k = lowered.get(str(c).strip().lower(), None)
            if k is not None:
                return metadata_kv.get(k, None)
        return None

    # Default all current XLSX experiments to "Lag" unless the sheet explicitly says otherwise.
    mode_val = _meta_lookup("Mode", "Experiment mode", "Run mode")
    if mode_val is None:
        exp.mode = "Lag"
    else:
        exp.mode = str(mode_val).strip()

    # Initial / final solution weights [g] -> used to compute overflow mass [g] if desired.
    w0 = _meta_lookup("Initial weight of solution (g)", "Initial solution weight (g)", "Initial weight (g)")
    wf = _meta_lookup("Final weight of solution (g)", "Final solution weight (g)", "Final weight (g)")
    if w0 is not None:
        try:
            exp.M_F0_g = float(w0)
        except Exception:
            pass
    if wf is not None and w0 is not None:
        try:
            exp.M_O_g = float(w0) - float(wf)  # overflow mass [g] (positive if mass lost)
        except Exception:
            pass

    # Salt identities -> component_names
    # Workbooks vary; we support:
    #   - a single "Salts" cell with comma-separated names
    #   - "Salt 1", "Salt 2", ... entries
    salts_cell = _meta_lookup("Salts", "Salt IDs", "Salt identification", "Salt names")
    salts: list[str] = []
    if salts_cell is not None:
        salts = [s.strip() for s in str(salts_cell).replace(";", ",").split(",") if s.strip()]

    if not salts:
        # Collect any keys like "Salt 1", "Salt 2", ...
        for k, v in (metadata_kv or {}).items():
            ks = str(k).strip().lower()
            if ks.startswith("salt"):
                if v is not None and str(v).strip() != "":
                    salts.append(str(v).strip())


    # If metadata blocks did not specify salts, fall back to the sheet-level convention:
    #   - Cell E1 contains a free-text description that includes salt formulas.
    # This makes the loader robust even if the metadata table (J/K) shifts around.
    if not salts:
        try:
            e1 = get_excel_cell_value(sheet, "E1")          # raw free-text (may be None)
            salts = extract_salt_formulas(e1)              # e.g., ["NaCl"], ["LaCl3","NaCl"]
        except Exception:
            salts = []

    if salts:
        exp.component_names = salts
        exp.num_components = len(salts)

    # ---------------------------------------------------------------------
    # 7b) Optional: populate initial guesses from an inferred MAT library DB
    # ---------------------------------------------------------------------
    # XLSX workbooks usually do not store Lp0/B0/sigma0/theta0. If you provide
    # an initial_guess_db (inferred from historical .mat files), we select the
    # best match by component_names and populate these fields. Otherwise they
    # remain None and can be set later.
    if initial_guess_db is not None:
        apply_initial_guesses_from_db(exp, initial_guess_db)


    # Calibration + assay blocks are parsed but not yet mapped into canonical fields here.
    # They remain available for future conversion logic once you finalize the workbook schema.
    # (Keeping them local avoids expanding ExperimentalData until needed.)

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
        exp.add_vial(v)

    # -------------------------------------------------------------------------
    # 10) Return the populated experiment container (validation is done by caller)
    # -------------------------------------------------------------------------
    # ---------------------------------------------------------------------
    # 10A) Best-effort defaults for operating conditions (it runs now, but honest)
    # ---------------------------------------------------------------------
    # Many XLSX workbooks do not explicitly include operating-condition fields
    # that legacy MAT files always provide in data_config (delP, Temp, Am, rho).
    #
    # To keep the loader file-agnostic for downstream Pyomo code, we fill
    # missing values from MAT_GLOBAL_DEFAULTS, and record each default fill
    # in exp.used_defaults so validate_experiment() can emit a warning.
    if exp.delP_bar is None:
        exp.delP_bar = MAT_GLOBAL_DEFAULTS["delP_bar"]  # [bar]
        exp.used_defaults.append("delP_bar")
    if exp.Temp_K is None:
        exp.Temp_K = MAT_GLOBAL_DEFAULTS["Temp_K"]      # [K]
        exp.used_defaults.append("Temp_K")
    if exp.Am_cm2 is None:
        exp.Am_cm2 = MAT_GLOBAL_DEFAULTS["Am_cm2"]      # [cm^2]
        exp.used_defaults.append("Am_cm2")
    if exp.rho_g_cm3 is None:
        exp.rho_g_cm3 = MAT_GLOBAL_DEFAULTS["rho_g_cm3"]  # [g/cm^3]
        exp.used_defaults.append("rho_g_cm3")

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



# =============================================================================
# XLSX column schema + explicit aliases (prevents silent misreads while tolerating
# harmless header drift across Excel sheets).
# =============================================================================

# Canonical column names used throughout the loader.
CANONICAL_REQUIRED_TS_COLS = {
    "Time (s)",
    "Mass (g)",
    "Vial Swap",
    "Retentate Cond @ Temp (uS/cm)",
}

# -----------------------------------------------------------------------------
# Legacy MAT-derived operating-condition defaults (best-effort).
#
# Policy:
#   - If XLSX provides an explicit value, use it.
#   - Otherwise, fill from these defaults so the loader is "it runs now"
#     while remaining honest: exp.used_defaults records each default fill.
# -----------------------------------------------------------------------------
MAT_GLOBAL_DEFAULTS = {
    "delP_bar": 4.1222251595,  # [bar]
    "Temp_K": 298.0,           # [K]
    "Am_cm2": 4.1,             # [cm^2]
    "rho_g_cm3": 1.0,          # [g/cm^3]
}



# Canonical metadata labels (key/value pairs) expected in the XLSX "metadata" block.
# These labels may appear anywhere on the sheet; the parser searches for them and
# reads the value from the cell immediately to the right.
XLSX_METADATA_KEYS = [
    # Experiment-level metadata labels (stored as row-wise key/value pairs in the sheet)
    "Experiment Name",
    "Initial Solution Weight",
    "Final Solution Weight",
    "Final Vial Weight",
    "Final Vial w/ Solution",
    "ICP Calibration Points (#)",
    "Number of Salts",
    # NOTE: Salt names are parsed dynamically as "Salt 1", "Salt 2", ...
]
# Explicit aliases we accept for each canonical column name.
# This is intentionally NOT fuzzy matching; only whitelisted alternatives are allowed.
TS_COL_ALIASES = {
    "Time (s)": {"Time (s)", "Time(s)", "Time [s]", "Time_sec", "Time seconds"},
    "Mass (g)": {"Mass (g)", "Mass(g)", "Mass [g]", "Mass_g"},
    "Vial Swap": {"Vial Swap", "VialSwap", "Vial swap", "Swap", "Vial Change"},
    "Retentate Cond @ Temp (uS/cm)": {
        "Retentate Cond @ Temp (uS/cm)",
        "Retentate Cond (uS/cm)",
        "Retentate Conductivity (uS/cm)",
        "Retentate Cond @ Temp",
        "Retentate Cond",
    },
}

def canonicalize_ts_header_token(token: str) -> Optional[str]:
    """
    Map a raw header cell token to a canonical time-series column name if recognized.

    Inputs:
        token: header cell value converted to stripped string

    Output:
        canonical column name if token matches an allowed alias, else None
    """
    for canonical, aliases in TS_COL_ALIASES.items():
        if token in aliases:
            return canonical
    return None



def parse_time_series_table(sheet: object) -> pd.DataFrame:
    """
    Parse and return the time-series table from an Excel sheet handle.

    What this function does:
        - Reads the entire sheet as raw cells (no header).
        - Locates the row containing the time-series header using explicit canonical
          column names (with whitelisted aliases).
        - Builds a DataFrame with those headers, then renames recognized aliases
          to canonical names so downstream code can use a stable schema.

    Inputs:
        sheet:
            Dict-like "sheet handle" returned by read_excel_sheet(...)

    Output:
        pd.DataFrame:
            Parsed time-series table with canonical column names.

    Raises:
        ValueError if the header row cannot be found or required columns are missing.
        TypeError if sheet is not the expected handle format.
    """

    # Defensive: ensure sheet is the handle we expect (from read_excel_sheet)
    if not isinstance(sheet, dict):
        raise TypeError(
            "parse_time_series_table expected 'sheet' to be a dict-like sheet handle "
            "(from read_excel_sheet)."
        )

    xls = sheet.get("excel_file", None)
    sheet_name = sheet.get("sheet_name", None)

    if xls is None or sheet_name is None:
        raise TypeError(
            "Sheet handle is missing required keys. Expected keys: "
            "'excel_file' and 'sheet_name'."
        )

    # Read the entire sheet as raw cells (no header)
    raw = pd.read_excel(xls, sheet_name=sheet_name, header=None, dtype=object)

    # Find the header row by scanning for a row that contains all canonical required columns,
    # allowing whitelisted aliases.
    header_row_idx: Optional[int] = None
    scan_limit = min(len(raw), 200)

    for i in range(scan_limit):
        row = raw.iloc[i].tolist()

        # Build a set of recognized CANONICAL column names present in this row.
        recognized: set[str] = set()

        for cell in row:
            if cell is None or (isinstance(cell, float) and np.isnan(cell)):
                continue
            token = str(cell).strip()
            if token == "":
                continue

            canon = canonicalize_ts_header_token(token)
            if canon is not None:
                recognized.add(canon)

        # If we recognized all canonical required columns, this is the header row.
        if CANONICAL_REQUIRED_TS_COLS.issubset(recognized):
            header_row_idx = i
            break

    if header_row_idx is None:
        raise ValueError(
            f"Could not locate time-series header row in sheet '{sheet_name}'. "
            f"Expected to find (canonical or aliased) required columns: "
            f"{sorted(CANONICAL_REQUIRED_TS_COLS)}"
        )

    # Build the DataFrame using the detected header row
    header_values = raw.iloc[header_row_idx].tolist()

    columns: List[str] = []
    for j, h in enumerate(header_values):
        if h is None or (isinstance(h, float) and np.isnan(h)) or str(h).strip() == "":
            columns.append(f"Unnamed: {j}")
        else:
            columns.append(str(h).strip())

    data = raw.iloc[header_row_idx + 1 :].copy()
    data.columns = columns
    data = data.dropna(how="all")

    # Rename any recognized aliases to canonical names (explicit, non-fuzzy).
    rename_map: dict[str, str] = {}
    for col in data.columns:
        canon = canonicalize_ts_header_token(str(col).strip())
        if canon is not None:
            rename_map[col] = canon

    if rename_map:
        data = data.rename(columns=rename_map)

    # Ensure canonical required columns exist in the final table
    missing = [c for c in CANONICAL_REQUIRED_TS_COLS if c not in data.columns]
    if missing:
        raise ValueError(
            f"Time-series table found in sheet '{sheet_name}', but missing canonical required columns: {missing}. "
            f"Available columns: {list(data.columns)}"
        )

    return data


def parse_experiment_metadata_kv(sheet: object) -> dict:
    """
    Parse the Excel *metadata key/value block* into a dict.

    Excel nuance (based on your sheet layout):
        - There is a region where labels like "Experiment Name:" appear in one cell,
          and the corresponding value appears in the *cell immediately to the right*.
        - In your current templates, this is usually columns J (labels) and K (values),
          but we do NOT hard-code those columns: we scan the whole sheet so the block
          can move without breaking the loader.

    Inputs:
        sheet:
            Sheet-handle returned by read_excel_sheet(...).
            In this codebase it is a dict with:
                - "excel_file": pandas.ExcelFile
                - "sheet_name": str
                - "path": Path

    Output:
        metadata: dict
            Dictionary mapping label strings -> raw cell values.
            Example keys (as they appear in the sheet):
                - "Experiment Name:"
                - "Initial Solution Weight:"
                - "Final Solution Weight:"
                - "ICP Calibration Points (#):"
                - "Number of Salts:"
                - "Salt 1:"
                - "Salt 2:"

            Values are kept "as-is" (may be str/float/int).
            Downstream code is responsible for casting numeric fields.

    Design choice:
    We only interpret key/value pairs whose *label* matches one of the expected
    XLSX metadata labels in `XLSX_METADATA_KEYS` (e.g., "Experiment Name",
    "Initial Solution Weight", "Salt 1", ...).

    Matching is whitespace- and colon-tolerant:
        - "Experiment Name", "Experiment Name:", and "  Experiment Name : " all match.

    For each matched label cell, the value is read from the immediate right
    neighbor (same row, col+1).

    This keeps parsing robust to row shifts and to the metadata block moving
    to a different set of columns.
"""
    if not isinstance(sheet, dict):
        raise TypeError("parse_experiment_metadata_kv expected a dict-like sheet handle.")

    xls = sheet.get("excel_file", None)
    sheet_name = sheet.get("sheet_name", None)
    if xls is None or sheet_name is None:
        raise TypeError("Sheet handle missing 'excel_file' or 'sheet_name'.")

    # Read the entire sheet as raw cells (no header), so we can scan for labels anywhere.
    raw = pd.read_excel(xls, sheet_name=sheet_name, header=None, dtype=object)

    metadata: dict = {}

    # Build a normalized lookup for expected metadata keys once.
    #
    # We normalize by:
    #   - stripping whitespace
    #   - removing a trailing ':' if present
    #   - collapsing internal whitespace
    #   - lowercasing
    #
    # This makes matching robust to minor formatting changes in Excel.
    expected = {" ".join(k.strip().split()).lower(): k for k in XLSX_METADATA_KEYS}

    # Scan a bounded region for speed; the metadata block is typically near the top
    # of the sheet. If future templates move it further down, increase this limit.
    max_rows = min(len(raw), 200)
    max_cols = raw.shape[1]

    for i in range(max_rows):
        for j in range(max_cols - 1):  # -1 because we always read value at (j+1)
            label = raw.iat[i, j]
            value = raw.iat[i, j + 1]

            # Only strings can be metadata labels.
            if not isinstance(label, str):
                continue

            # Normalize label for matching.
            # Example: "  Experiment Name : " -> "experiment name"
            label_norm = " ".join(label.strip().rstrip(":").split()).lower()

            canonical_key: Optional[str] = None

            # 1) Direct match against expected (non-salt) keys
            if label_norm in expected:
                canonical_key = expected[label_norm]  # e.g., "Experiment Name"

            # 2) Dynamic salt parsing: accept "Salt 1", "Salt 2", ... up to any N.
            #    We do NOT bake Salt 1 / Salt 2 into XLSX_METADATA_KEYS because the
            #    number of salts varies by experiment.
            if canonical_key is None and label_norm.startswith("salt"):
                # Allow variants like "Salt1", "Salt 1", "salt  12", etc.
                m_salt = re.match(r"^salt\s*(\d+)$", label_norm)
                if m_salt is not None:
                    salt_idx = int(m_salt.group(1))
                    canonical_key = f"Salt {salt_idx}"

            # If this label isn't one we care about, skip.
            if canonical_key is None:
                continue

            # Skip empty/NaN values so we don't populate meaningless entries.
            if value is None:
                continue
            if isinstance(value, float) and np.isnan(value):
                continue
            if isinstance(value, str) and value.strip() == "":
                continue

            # Store using the canonical label (no colon) for stable downstream mapping.
            metadata[canonical_key] = value

    return metadata



def parse_icp_assay_table(sheet: object) -> Optional[pd.DataFrame]:
    """
    Parse the ICP-OES assay table region into a DataFrame.

    Excel nuance (based on your sheet layout):
        - The assay table is a conventional column-oriented table with headers like:
            "Vial Number", "Sample Volume (mL)", "Nitric Acid Volume (mL)", "ICP Salt 1 (mg/L)", ...
        - In your current templates this lives around columns O..V,
          but we do NOT hard-code the column letters:
            *we locate the header row by searching for required header tokens*,
            then take the contiguous block of non-empty headers.

    Inputs:
        sheet:
            Sheet-handle from read_excel_sheet(...)

    Output:
        assay_df:
            DataFrame containing the assay rows (one row per vial).
            Returns None if no assay table can be located.

    Notes:
        - This function does NOT enforce units (it preserves whatever the sheet provides).
        - It does NOT attempt to interpret which salt is "salt 1" vs "salt 2".
          That mapping can be done using metadata keys like "Salt 1:".
    """
    if not isinstance(sheet, dict):
        raise TypeError("parse_icp_assay_table expected a dict-like sheet handle.")

    xls = sheet.get("excel_file", None)
    sheet_name = sheet.get("sheet_name", None)
    if xls is None or sheet_name is None:
        raise TypeError("Sheet handle missing 'excel_file' or 'sheet_name'.")

    raw = pd.read_excel(xls, sheet_name=sheet_name, header=None, dtype=object)

    # Minimal "signature" headers that identify the assay table in your templates.
    required_tokens = {"Vial Number", "Sample Volume (mL)", "Nitric Acid Volume (mL)"}

    header_row_idx: Optional[int] = None
    scan_limit = min(len(raw), 250)

    # Find a row that contains all required assay tokens.
    for i in range(scan_limit):
        row = raw.iloc[i].tolist()
        normalized = set()
        for cell in row:
            if cell is None or (isinstance(cell, float) and np.isnan(cell)):
                continue
            s = str(cell).strip()
            if s:
                normalized.add(s)
        if required_tokens.issubset(normalized):
            header_row_idx = i
            break

    if header_row_idx is None:
        # Assay table not found in this sheet.
        return None

    header_values = raw.iloc[header_row_idx].tolist()

    # Determine the contiguous header block starting from the first required token.
    # This prevents us from accidentally capturing the calibration block that follows
    # after a blank separator column.
    start_col = None
    for j, h in enumerate(header_values):
        if str(h).strip() == "Vial Number":
            start_col = j
            break
    if start_col is None:
        return None

    end_col = start_col
    for j in range(start_col, len(header_values)):
        h = header_values[j]
        if h is None or (isinstance(h, float) and np.isnan(h)) or str(h).strip() == "":
            break  # stop at the first blank column; calibration block comes after a separator
        end_col = j + 1  # slice end is exclusive

    # Build DataFrame for the detected block.
    block_headers = [str(h).strip() for h in header_values[start_col:end_col]]
    data = raw.iloc[header_row_idx + 1 :, start_col:end_col].copy()
    data.columns = block_headers
    data = data.dropna(how="all")

    return data



def parse_calibration_block(sheet: object) -> dict:
    """
    Parse the ICP calibration block(s) into per-salt calibration arrays.

    Excel nuance (based on your sheet layout):
        - Calibration blocks are labeled with headers like:
            "ICP Salt 1 Calibration", "ICP Salt 2 Calibration", ...
        - Under each such header there are typically two columns:
            "Concentration" and "Intensity (cps)"
        - The block may move within the sheet, so we locate these labels by scanning.

    Inputs:
        sheet:
            Sheet-handle from read_excel_sheet(...)

    Output:
        calib: dict
            Mapping from salt index (int) -> dict with keys:
                - "x": 1D numpy float array of concentrations (units as in sheet)
                - "y": 1D numpy float array of intensities (cps)

            Example:
                calib[1]["x"], calib[1]["y"] correspond to "ICP Salt 1 Calibration"

    Notes:
        - This function is intentionally "honest": it does not fit a calibration curve.
          It only returns the raw points so a later measurement model can decide how to use them.
        - If a calibration header is found but no numeric data exists below it, that salt is omitted.
    """
    if not isinstance(sheet, dict):
        raise TypeError("parse_calibration_block expected a dict-like sheet handle.")

    xls = sheet.get("excel_file", None)
    sheet_name = sheet.get("sheet_name", None)
    if xls is None or sheet_name is None:
        raise TypeError("Sheet handle missing 'excel_file' or 'sheet_name'.")

    raw = pd.read_excel(xls, sheet_name=sheet_name, header=None, dtype=object)

    calib: dict = {}

    # We search the top part of the sheet where calibration typically lives.
    # Increase bounds if templates change.
    max_rows = min(len(raw), 300)
    max_cols = raw.shape[1]

    # Match "ICP Salt 1 Calibration" etc (case-insensitive).
    pat = re.compile(r"^ICP\s+Salt\s+(\d+)\s+Calibration", re.IGNORECASE)

    for i in range(max_rows):
        for j in range(max_cols):
            cell = raw.iat[i, j]
            if not isinstance(cell, str):
                continue
            m = pat.match(cell.strip())
            if not m:
                continue

            salt_idx = int(m.group(1))

            # Template pattern:
            #   row i contains the "ICP Salt k Calibration" label
            #   row i+1 contains subheaders ("Concentration", "Intensity (cps)")
            #   data starts at row i+2 and continues until blank
            subheader_row = i + 1
            data_start = i + 2

            if subheader_row >= len(raw):
                continue

            # The two columns we expect are at j (Concentration) and j+1 (Intensity),
            # because the label is typically placed above the first of the two.
            x_col = j
            y_col = j + 1
            if y_col >= max_cols:
                continue

            # Collect numeric pairs until we hit a blank row for x and y.
            xs = []
            ys = []
            for r in range(data_start, len(raw)):
                x = raw.iat[r, x_col]
                y = raw.iat[r, y_col]

                # Stop when both are empty -> end of calibration points.
                x_missing = (x is None) or (isinstance(x, float) and np.isnan(x)) or (isinstance(x, str) and x.strip() == "")
                y_missing = (y is None) or (isinstance(y, float) and np.isnan(y)) or (isinstance(y, str) and y.strip() == "")
                if x_missing and y_missing:
                    break

                # Best-effort numeric cast; skip rows that can't be converted.
                try:
                    xs.append(float(x))
                    ys.append(float(y))
                except Exception:
                    continue

            if len(xs) > 0 and len(ys) > 0:
                calib[salt_idx] = {"x": np.asarray(xs, dtype=float), "y": np.asarray(ys, dtype=float)}

    return calib



def find_vial_row(vial_table: object, vial_number: int) -> Optional[dict]:
    """
    Find and return the assay-table row corresponding to a given vial number.

    Inputs:
        vial_table:
            Table returned by parse_vial_data_table(...). Typically a pandas DataFrame.

        vial_number:
            1-based vial index (1, 2, 3, ...).

    Output:
        dict or None:
            If found, returns a dict {column_name: value} for the matching row.
            If not found, returns None.
    """

    if vial_table is None:
        return None

    # DataFrame-like case (preferred)
    if hasattr(vial_table, "columns") and hasattr(vial_table, "iterrows"):
        df = vial_table

        # Strategy A: explicit numeric vial column
        candidate_cols = ["Vial", "Vial Number", "Vial#", "Vial_ID", "Vial ID"]
        for col in candidate_cols:
            if col in df.columns:
                try:
                    series = pd.to_numeric(df[col], errors="coerce")
                    match_idx = series[series == vial_number].index
                    if len(match_idx) > 0:
                        return df.loc[match_idx[0]].to_dict()
                except Exception:
                    pass

        # Strategy B: "Vial X" label in first column
        if len(df.columns) > 0:
            first_col = df.columns[0]
            target1 = f"vial {vial_number}"
            target2 = f"vial{vial_number}"

            for _, row in df.iterrows():
                cell = row.get(first_col, None)
                if cell is None or (isinstance(cell, float) and np.isnan(cell)):
                    continue
                s = str(cell).strip().lower()
                if s == target1 or s.replace(" ", "") == target2:
                    return row.to_dict()

        # Strategy C: scan whole row for "Vial X"
        token1 = f"vial {vial_number}"
        token2 = f"vial{vial_number}"

        for _, row in df.iterrows():
            for cell in row.values:
                if cell is None or (isinstance(cell, float) and np.isnan(cell)):
                    continue
                s = str(cell).strip().lower()
                if s == token1 or s == token2:
                    return row.to_dict()

        return None

    # Fallback: list-of-dicts
    if isinstance(vial_table, list):
        for row in vial_table:
            if not isinstance(row, dict):
                continue

            for key in ["Vial", "Vial Number", "number", "vial_number"]:
                if key in row:
                    try:
                        if int(row[key]) == vial_number:
                            return row
                    except Exception:
                        pass

            token1 = f"vial {vial_number}"
            token2 = f"vial{vial_number}"
            for val in row.values():
                if isinstance(val, str):
                    s = val.strip().lower()
                    if s == token1 or s == token2:
                        return row

        return None

    return None

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
            The expected column name.

            IMPORTANT:
                Excel headers sometimes include harmless variations (extra spaces,
                different capitalization, or trailing whitespace). This helper
                checks both exact matches and a normalized (case/space-insensitive)
                match to be more robust.

    Output:
        bool:
            True if the column exists, otherwise False.
    """

    def _norm(s: object) -> str:
        """Normalize a potential column label for robust comparison."""
        # Convert to string, strip edges, collapse internal whitespace, and lowercase.
        return " ".join(str(s).strip().split()).casefold()

    # Most common case: pandas DataFrame (has a .columns attribute).
    # We check attribute existence so this does not crash on non-DataFrame inputs.
    if hasattr(table, "columns"):
        # 1) Fast path: exact membership check (preserves the user's explicit name).
        if colname in table.columns:
            return True

        # 2) Robust path: normalized comparison (handles whitespace/case drift).
        target = _norm(colname)
        return any(_norm(c) == target for c in list(table.columns))

    # Fallback case: dict-like object (keys represent column names).
    # This supports cases where parsing returns a mapping instead of a DataFrame.
    if isinstance(table, dict):
        if colname in table:
            return True
        target = _norm(colname)
        return any(_norm(k) == target for k in table.keys())

    # If we do not recognize the table type, we cannot guarantee column existence.
    # Being conservative here is safer than guessing.
    return False

def get_excel_cell_value(sheet: dict, cell_ref: str) -> object:
    """
    Read a single Excel cell value from the resolved sheet.

    Why this exists:
        Your Excel workbooks store some experiment metadata in fixed cells
        (e.g., salts used in the experiment are in cell E1 per your convention).
        We keep this as a tiny helper so the loader stays readable.

    Inputs:
        sheet:
            Dict-like sheet handle returned by read_excel_sheet(...).
            Must contain:
                - "path": Path to workbook
                - "sheet_name": resolved sheet name (str)

        cell_ref:
            Excel A1-style cell reference, e.g. "E1".

    Output:
        object:
            The raw cell value (string/float/etc.) or None if empty.

    Notes:
        Implementation uses pandas with header=None to avoid making any assumptions
        about table headers elsewhere in the sheet.
    """
    if not isinstance(sheet, dict):
        raise TypeError("get_excel_cell_value expects a sheet handle dict from read_excel_sheet().")

    path = sheet.get("path", None)
    sheet_name = sheet.get("sheet_name", None)
    if path is None or sheet_name is None:
        raise TypeError("Sheet handle missing required keys: 'path' and 'sheet_name'.")

    # Parse the cell reference (e.g., 'E1' -> col_letter='E', row_number=1)
    m = re.fullmatch(r"([A-Za-z]+)(\d+)", cell_ref.strip())
    if m is None:
        raise ValueError(f"Invalid cell reference: {cell_ref!r}. Expected like 'E1'.")
    col_letter = m.group(1).upper()
    row_number = int(m.group(2))

    # Read exactly one cell range: the requested column, and the requested row.
    # - usecols selects the Excel column by letter
    # - skiprows jumps to the row-1 (0-based)
    # - nrows=1 reads exactly one row
    df = pd.read_excel(
        path,
        sheet_name=sheet_name,
        header=None,
        usecols=col_letter,
        skiprows=row_number - 1,
        nrows=1,
        dtype=object,
    )

    if df.empty:
        return None

    # df will have one row, one column (index 0, column col_letter or 0 depending on pandas)
    val = df.iloc[0, 0]
    if val is None:
        return None
    if isinstance(val, float) and np.isnan(val):
        return None
    return val


def extract_salt_formulas(text: object) -> List[str]:
    """
    Best-effort extraction of salt formulas from a free-text string.

    Your XLSX sheets store the salt information in cell E1 as phrases like:
        "Feed: 1 mM NaCl, Diafiltrate: 50 mM NaCl, 350RPM"
        "... 0.2mM LaCl3 and 50mM NaCl ..."

    Inputs:
        text:
            Free-text description (often a string). None returns [].

    Output:
        List[str]:
            Unique salt formulas (order preserved) such as ["NaCl"], ["LaCl3", "NaCl"], etc.

    Notes:
        This is intentionally conservative:
            - We only match common electrolyte patterns (cation+anion with optional stoichiometry).
            - If your project later includes new ions, extend the regex or add a mapping table.
    """
    if text is None:
        return []

    s = str(text)

    # Regex: common cations + common anions with optional digits (e.g., CaCl2, Na2SO4, LaCl3)
    # You can extend the cation/anion sets as your dataset grows.
    cations = r"Na|K|Li|Ca|Mg|La|Al|Zn|Fe|Cu|Ni|Co|Mn|Sr|Ba"
    anions = r"Cl|Br|I|NO3|SO4|CO3|HCO3|PO4"
    pattern = re.compile(rf"\b(?:{cations})\d*(?:{anions})\d*\b")

    found = pattern.findall(s)

    # De-duplicate while preserving order
    out: List[str] = []
    for f in found:
        if f not in out:
            out.append(f)
    return out


def _conductivity_to_concentration_series(
    *,
    cond_uS_cm: np.ndarray,
    temp_K: float,
    model: str,
    model_params: Dict[str, object],
    output_units: str = "mM",
    salt_name: Optional[str] = None,
) -> np.ndarray:
    """
    Convert conductivity time-series [uS/cm] -> concentration series (mM or M).

    Inputs:
        cond_uS_cm:
            Conductivity series in uS/cm (from XLSX sensors).
        temp_K:
            Temperature [K].
        model:
            "variant_shedlovsky" or "msa".
        model_params:
            Dict of parameters passed into conductivity_paper.py.
        output_units:
            "mM" or "M".
        salt_name:
            Optional salt identity (e.g. "NaCl").
            Used ONLY for parameter autofill if missing.

    Output:
        np.ndarray:
            Concentration series aligned with input conductivity.
    """

    # ------------------------------------------------------------
    # 1) Import conductivity_paper WITHOUT modifying it
    # ------------------------------------------------------------
    try:
        import conductivity_paper as cp
    except ModuleNotFoundError:
        import importlib.util
        from importlib.machinery import SourceFileLoader

        here = Path(__file__).resolve().parent
        paper_path = here / "conductivity_paper.py"

        spec = importlib.util.spec_from_loader(
            "conductivity_paper",
            SourceFileLoader("conductivity_paper", str(paper_path)),
        )
        cp = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cp)

    # ------------------------------------------------------------
    # 2) Convert units: uS/cm -> mS/cm
    # ------------------------------------------------------------
    cond_uS_cm = np.asarray(cond_uS_cm, dtype=float).reshape(-1)
    cond_mS_cm = cond_uS_cm / 1000.0

    conc_out = np.full_like(cond_mS_cm, fill_value=np.nan, dtype=float)

    model = model.lower()

    # ============================================================
    # VARIANT SHEDLOVSKY
    # ============================================================
    if model in ["variant_shedlovsky", "shedlovsky"]:

        # --------------------------------------------------------
        # Robust + honest parameter priority order
        # --------------------------------------------------------

        # Make local copy so we never mutate caller dict
        mp = dict(model_params)

        # If salt_name available -> attempt autofill
        if salt_name is not None:
            mp = _autofill_variant_shedlovsky_params(mp, salt_name=salt_name)

        # Required parameters (must exist AND not be None)
        required = ["epsilon", "eta", "a", "z_1", "z_2", "lambda_0_cation", "lambda_0_anion"] # We require ionic pieces; lambda_0 itself is optional because we can compute it.

        missing = [k for k in required if mp.get(k) is None]
        if missing:
            raise ValueError(
                f"Missing variant Shedlovsky params: {missing}. "
                f"Provide explicitly or ensure salt_name inference works."
            )

        # Optional ionic lambdas (only needed if conductivity_paper requires them)
        if mp.get("lambda_0_cation") is None or mp.get("lambda_0_anion") is None:
            raise ValueError(
                "lambda_0_cation and lambda_0_anion are required for Shedlovsky model."
            )

        epsilon = float(mp["epsilon"])
        eta = float(mp["eta"])
        # Priority order for lambda_0 (limiting conductivity parameter):
            #   1) If provided explicitly, use it.
            #   2) Else compute lambda_0 = lambda_0_cation + lambda_0_anion
        if mp.get("lambda_0") is not None:
            lambda_0 = float(mp["lambda_0"])
        else:
            lambda_0 = float(mp["lambda_0_cation"]) + float(mp["lambda_0_anion"])

        a = float(mp["a"])
        z_1 = int(mp["z_1"])
        z_2 = int(mp["z_2"])
        lambda_0_cation = float(mp["lambda_0_cation"])
        lambda_0_anion = float(mp["lambda_0_anion"])

        # Forward mapping: conc_M -> conductivity
        def fwd(conc_M: float) -> float:
            """Fwd.

            Args:
                conc_M: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return float(
                cp.variant_shedlovsky(
                    [conc_M], temp_K,
                    epsilon, eta, lambda_0, a, z_1, z_2,
                    lambda_0_cation, lambda_0_anion,
                )[0]
            )

        # Invert each time point
        for i, y in enumerate(cond_mS_cm):
            if np.isnan(y):
                continue
            try:
                conc_M = _invert_monotone_1d(
                    fwd_func=fwd,
                    y_target=float(y),
                    x_lo=0.0,
                    x_hi=6.0,
                )
                conc_out[i] = conc_M
            except Exception:
                # Keep conversion best-effort: if a single point cannot be inverted,
                # mark it missing and continue converting the rest of the series.
                conc_out[i] = np.nan

        # Unit conversion
        if output_units.lower() in ["mm", "mmol/l", "mmolar", "mm"]:
            conc_out *= 1000.0  # M -> mM

        return conc_out

    # ============================================================
    # MSA MODEL
    # ============================================================
    if model in ["msa", "mean_spherical_approximation"]:

        required = ["valency", "diameters", "diff_coeff", "eta", "epsilon", "lambda_0"]
        missing = [k for k in required if model_params.get(k) is None]
        if missing:
            raise ValueError(f"Missing MSA params: {missing}")

        valency = list(model_params["valency"])
        diameters = list(model_params["diameters"])
        diff_coeff = list(model_params["diff_coeff"])
        eta = float(model_params["eta"])
        epsilon = float(model_params["epsilon"])
        lambda_0 = list(model_params["lambda_0"])

        n_salts = int(model_params.get("n_salts", 1))
        if n_salts < 1 or n_salts > 3:
            raise ValueError("MSA wrapper supports 1–3 salts.")

        salt_ratios = model_params.get("salt_ratios", None)
        if salt_ratios is None:
            salt_ratios = [1.0] * n_salts

        salt_ratios = np.asarray(salt_ratios, dtype=float)
        salt_ratios /= np.sum(salt_ratios)

        def fwd(total_mM: float) -> float:
            """Fwd.

            Args:
                total_mM: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            salts = [float(total_mM * r) for r in salt_ratios]

            if n_salts == 1:
                return float(
                    cp.msa_transport(valency, diameters, diff_coeff,
                                     temp_K, eta, epsilon, lambda_0,
                                     salts[0])[0]
                )
            if n_salts == 2:
                return float(
                    cp.msa_transport(valency, diameters, diff_coeff,
                                     temp_K, eta, epsilon, lambda_0,
                                     salts[0], salts[1])[0]
                )

            return float(
                cp.msa_transport(valency, diameters, diff_coeff,
                                 temp_K, eta, epsilon, lambda_0,
                                 salts[0], salts[1], salts[2])[0]
            )

        for i, y in enumerate(cond_mS_cm):
            if np.isnan(y):
                continue
            try:
                total_mM = _invert_monotone_1d(
                    fwd_func=fwd,
                    y_target=float(y),
                    x_lo=0.0,
                    x_hi=6000.0,
                )
                conc_out[i] = total_mM
            except Exception:
                conc_out[i] = np.nan

        if output_units.lower() in ["m", "mol/l", "molar"]:
            conc_out /= 1000.0  # mM -> M

        return conc_out

    raise ValueError(f"Unknown conductivity->concentration model: {model}")



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

# =============================================================================
# MAT schema constants (expected struct keys)
# =============================================================================

# These constants define the expected *shape* of the MATLAB experiment struct.
# They are not used to hard-fail ingestion; they exist to keep MATLAB parsing explicit and auditable.
MAT_REQUIRED_ROOT_KEYS = {
    "dataset",
    "filename",
    "mode",
    "data_config",
    "data_raw",
}

MAT_REQUIRED_CONFIG_KEYS = {
    "delP",
    "Temp",
    "Am",
    "rho",
}

MAT_REQUIRED_VIAL_KEYS = {
    "number",
    "time",
    "mass",
    "cV_avg",
    "cF_exp",
}


def load_from_mat(mat_path: Path, struct_key: Optional[str]) -> ExperimentalData:
    """
    Load ONE experiment from a MATLAB .mat file into ExperimentalData fields.

    Inputs:
        mat_path:
            Path to the .mat file on disk.

        struct_key:
            Name of the MAT variable containing the experiment struct.
            If None, we default to "data_stru" (the legacy convention in this project).

    Output:
        exp:
            ExperimentalData instance with explicit class fields populated wherever possible.

    Key design choices (consistent with pseudocode_v3.py):
        - We do NOT try to "interpret" signals (e.g., convert conductivity->concentration).
          We store what the MAT file provides honestly, with labels.
        - Missing fields remain None; structural problems are handled by validate_experiment(...).
        - Each .mat file contains exactly ONE experiment (legacy assumption for this project).
    """
    # -------------------------------------------------------------------------
    # 1) Load raw MAT content from disk (scipy.io.loadmat wrapper).
    # -------------------------------------------------------------------------
    raw = read_mat_file(mat_path)

    # -------------------------------------------------------------------------
    # 2) Choose which variable in the MAT file holds the experiment struct.
    #    Legacy files typically store this under "data_stru".
    # -------------------------------------------------------------------------
    key = struct_key if struct_key is not None else "data_stru"

    if key not in raw:
        # If the requested key is missing, attempt a conservative fallback:
        #   - ignore MATLAB metadata keys (e.g., "__header__")
        #   - if exactly one "real" key exists, treat it as the experiment struct
        real_keys = [k for k in raw.keys() if not str(k).startswith("__")]
        if len(real_keys) == 1:
            key = real_keys[0]
        else:
            raise KeyError(
                f"MAT file does not contain expected struct key '{key}'. "
                f"Available keys: {sorted(real_keys)}"
            )

    # -------------------------------------------------------------------------
    # 3) Convert MATLAB structs / cell arrays into nested Python types.
    #    This makes downstream parsing deterministic:
    #       MATLAB struct -> dict
    #       MATLAB cell   -> list
    #       numeric arrays stay numpy arrays
    # -------------------------------------------------------------------------
    d = mat_struct_to_python(raw[key])

    # -------------------------------------------------------------------------
    # 4) Initialize an empty ExperimentalData container and populate metadata.
    # -------------------------------------------------------------------------
    exp = ExperimentalData()
    exp.source = SourceType.MAT
    exp.filename = str(d.get("filename", mat_path.name))  # provenance filename
    exp.sheet_name = None  # MAT contains one experiment, no sheet concept

    # Experiment identity (these keys exist in legacy data_stru)
    exp.dataset_id = d.get("dataset", None)
    exp.mode = d.get("mode", None)

    # Retentate measurement behavior toggle (naming varies slightly across datasets)
    # Common legacy field names include "continuous_cF" (boolean).
    exp.continuous_retentate = d.get("continuous_cF", None)

    # -------------------------------------------------------------------------
    # 5) Map data_config fields (experiment-level configuration / operating conditions)
    # -------------------------------------------------------------------------
    cfg = d.get("data_config", {}) if isinstance(d.get("data_config", {}), dict) else {}

    # Operating conditions (units from legacy comments)
    exp.delP_bar = cfg.get("delP", None)    # [bar]
    exp.Temp_K = cfg.get("Temp", None)      # [K]
    exp.Am_cm2 = cfg.get("Am", None)        # [cm^2]
    exp.rho_g_cm3 = cfg.get("rho", None)    # [g/cm^3]

    # Feed/dialysate info (legacy names)
    exp.M_F0_g = cfg.get("M_F0", None)      # [g]
    exp.M_O_g = cfg.get("M_O", None)        # [g] overflow mass
    exp.C_D_value = cfg.get("C_D", None)    # [mMol/L] in legacy files (value only here)
    exp.C_D_units = "mMol/L" if exp.C_D_value is not None else None

    exp.C_F0_value = cfg.get("C_F0", None)  # [mMol/L] in legacy files (value only here)
    exp.C_F0_units = "mMol/L" if exp.C_F0_value is not None else None

    # Component bookkeeping:
    #   - legacy files may store nc (number of components) and namec (names)
    #   - we map into exp.num_components and exp.component_names
    namec = cfg.get("namec", None)
    if isinstance(namec, str):
        exp.component_names = [namec]
    elif isinstance(namec, (list, tuple, np.ndarray)):
        exp.component_names = [str(x) for x in list(namec)]
    else:
        exp.component_names = None

    exp.num_components = cfg.get("nc", None)

    # If 'namec' was not available, try to infer salt formulas from the legacy filename.
    # This is best-effort and only used to make downstream interfaces consistent.
    if exp.component_names is None:
        inferred = extract_salt_formulas(exp.filename)
        if inferred:
            exp.component_names = inferred
            if exp.num_components is None:
                exp.num_components = len(inferred)


    # Initial guesses / priors used in parameter estimation / MBDoE workflows
    exp.Lp0 = cfg.get("Lp0", None)
    exp.B0 = cfg.get("B0", None)
    exp.sigma0 = cfg.get("sigma0", None)

    # theta0 must be numeric (project requirement).
    # Legacy files often include theta0 already; if missing but Lp0/B0/sigma0 exist, assemble it.
    theta0 = cfg.get("theta0", None)
    if theta0 is None and (exp.Lp0 is not None) and (exp.B0 is not None) and (exp.sigma0 is not None):
        exp.theta0 = np.array([exp.Lp0, exp.B0, exp.sigma0], dtype=float)
    else:
        # Best-effort numeric conversion; if it fails, leave as None (validator can flag later).
        try:
            exp.theta0 = None if theta0 is None else np.asarray(theta0, dtype=float).reshape(-1)
        except Exception:
            exp.theta0 = None

    # -------------------------------------------------------------------------
    # 6) Map calibration arrays if present (legacy calibration for conductivity->concentration)
    #    NOTE: We store raw arrays; conversion happens in a measurement model later.
    # -------------------------------------------------------------------------
    # These names may vary across MAT datasets; we populate only if the keys exist.
    # Keep this minimal: if your MAT files use different names, we can extend later.
    for salt_idx in [1, 2]:
        x_key = f"salt{salt_idx}_calib_x"
        y_key = f"salt{salt_idx}_calib_y"
        if x_key in cfg:
            try:
                setattr(exp, x_key, np.asarray(cfg.get(x_key), dtype=float).reshape(-1))
            except Exception:
                setattr(exp, x_key, None)
        if y_key in cfg:
            try:
                setattr(exp, y_key, np.asarray(cfg.get(y_key), dtype=float).reshape(-1))
            except Exception:
                setattr(exp, y_key, None)

    # -------------------------------------------------------------------------
    # 7) Build the list of VialData objects from data_raw.
    #    In the legacy MAT schema, each vial is already a separate struct entry.
    # -------------------------------------------------------------------------
    exp.vials = []  # ensure empty list before appending

    raw_vials = d.get("data_raw", [])
    if isinstance(raw_vials, dict):
        # Defensive: if a single vial was stored as dict instead of list, wrap it.
        raw_vials = [raw_vials]

    if isinstance(raw_vials, (list, tuple)):
        for idx, rv in enumerate(raw_vials, start=1):
            # Each rv should be dict-like after mat_struct_to_python conversion.
            if not isinstance(rv, dict):
                continue

            v = VialData()
            v.number = int(rv.get("number", idx)) if rv.get("number", None) is not None else idx

            # Time series (units: seconds)
            v.time_s = None
            try:
                v.time_s = as_float_array(rv.get("time", None))
            except Exception:
                v.time_s = None

            # Mass series (units: grams) — optional
            if rv.get("mass", None) is not None:
                try:
                    v.mass_g = as_float_array(rv.get("mass", None))
                except Exception:
                    v.mass_g = None

            # Retentate signal series:
            # Legacy key: cF_exp.
            #
            # IMPORTANT legacy detail:
            #   In many published MAT datasets, cF_exp is a *single* retentate assay value
            #   (e.g., the final retentate concentration for that vial), NOT a time-series.
            #
            # Our loader uses ONE field (VialData.retentate_signal) for both MAT and XLSX.
            # To keep the schema consistent and satisfy validate_experiment(), we do:
            #   - If cF_exp is already time-aligned: keep it as-is.
            #   - If cF_exp is scalar / length-1: broadcast it across the vial's time grid.
            # This preserves the information content (it is still the same scalar assay),
            # while providing a time-aligned array for downstream code.
            cF_exp_raw = rv.get("cF_exp", None)

            v.retentate_signal = None
            if cF_exp_raw is not None:
                try:
                    ret_arr = as_float_array(cF_exp_raw)  # best-effort numeric array
                except Exception:
                    ret_arr = None

                if ret_arr is not None:
                    # If ret_arr looks like a scalar assay, broadcast it across time_s
                    if (v.time_s is not None) and (len(v.time_s) > 0) and (len(ret_arr) in (0, 1)):
                        # Empty ret_arr is treated as missing; length-1 is broadcast
                        v.retentate_signal = None if len(ret_arr) == 0 else np.full_like(v.time_s, ret_arr[0], dtype=float)
                    else:
                        # Assume it's already a time-series
                        v.retentate_signal = ret_arr

            # "cF_exp" is historically a concentration-like quantity in the MAT files.
            v.retentate_signal_name = "cF_exp"
            v.retentate_signal_units = "mMol/L"  # legacy data_config uses mMol/L; keep consistent label

            # Permeate concentration or signal (legacy key: cV_avg)
            cV = rv.get("cV_avg", None)
            if cV is None:
                v.cV_avg = None
            else:
                # Preserve scalar vs array when possible
                if isinstance(cV, (int, float, np.number)):
                    v.cV_avg = float(cV)
                else:
                    try:
                        arr = np.asarray(cV, dtype=float)
                        # If it's a 0-d array, convert to scalar float
                        v.cV_avg = float(arr) if arr.shape == () else arr.reshape(-1)
                    except Exception:
                        v.cV_avg = None

            exp.add_vial(v)

    # -------------------------------------------------------------------------
    # 8) Return the populated experiment container.
    #    Structural validity is checked by validate_experiment(...) at the entrypoint.
    # -------------------------------------------------------------------------
    return exp


def read_mat_file(mat_path: Path) -> Dict[str, object]:
    """Read a MATLAB ``.mat`` file and return the raw loaded content.

    Inputs:
        mat_path:
            Path to a MATLAB ``.mat`` file on disk.

    Output:
        raw:
            A dict-like mapping returned by ``scipy.io.loadmat``.
            Keys are variable names stored in the MAT file.
            Notes:
                - MATLAB metadata keys like ``__header__`` / ``__version__`` / ``__globals__``
                  may be present. We keep them here; downstream code should select the
                  experiment struct key (typically ``data_stru``).
                - We load with options that make MATLAB structs easier to work with
                  in Python (squeezed singleton dimensions, no record arrays when possible).

    Raises:
        FileNotFoundError:
            If mat_path does not exist.
        ValueError:
            If mat_path does not look like a .mat file.
    """
    # -------------------------------------------------------------------------
    # 0) Defensive checks
    # -------------------------------------------------------------------------
    if not isinstance(mat_path, Path):
        mat_path = Path(mat_path)

    if not mat_path.exists():
        raise FileNotFoundError(f"MAT file not found: {mat_path}")

    if mat_path.suffix.lower() != ".mat":
        raise ValueError(f"Expected a .mat file, got: {mat_path.name}")

    # -------------------------------------------------------------------------
    # 1) Load the MAT file using SciPy
    # -------------------------------------------------------------------------
    # squeeze_me=True:
    #   - convert MATLAB 1x1 arrays into scalars where possible
    #
    # struct_as_record=False:
    #   - represent MATLAB structs as simple Python objects (mat_struct) instead of
    #     numpy record arrays; this makes recursive conversion easier.
    #
    # simplify_cells=True (SciPy >= 1.7-ish, varies by environment):
    #   - convert MATLAB cell arrays into nested Python lists directly.
    #   - If unavailable, we fall back gracefully.
    try:
        raw = loadmat(
            mat_path,
            squeeze_me=True,
            struct_as_record=False,
            simplify_cells=True,  # may not exist in older SciPy; handled below
        )
    except TypeError:
        # Older SciPy versions don't support simplify_cells
        raw = loadmat(
            mat_path,
            squeeze_me=True,
            struct_as_record=False,
        )

    # -------------------------------------------------------------------------
    # 2) Return raw mapping
    # -------------------------------------------------------------------------
    # Downstream logic will typically do:
    #   key = "data_stru" if present else choose first non-metadata key
    #   d = mat_struct_to_python(raw[key])
    return raw


def mat_struct_to_python(mat_struct: object) -> object:
    """
    Recursively convert MATLAB-loaded objects into Python-native containers.

    Why this exists:
        scipy.io.loadmat returns a nested mix of MATLAB-specific objects:
            - MATLAB structs -> `mat_struct` objects with `_fieldnames`
            - MATLAB cell arrays -> object arrays (dtype=object) or Python lists (if simplify_cells=True)
            - MATLAB numeric arrays -> numpy ndarrays
        For a file-type agnostic loader, we need a predictable representation so
        `load_from_mat(...)` can map fields into ExperimentalData / VialData
        without scattering MATLAB-type checks everywhere.

    Inputs:
        mat_struct:
            Any object returned inside the dict from scipy.io.loadmat(...).
            This can be:
                - a MATLAB struct-like object (has `_fieldnames`)
                - numpy ndarray (numeric, object array, or structured)
                - list / tuple
                - scalar (numpy scalar, float, int, str, etc.)

    Output:
        A nested Python-native object built from:
            - dict for MATLAB structs
            - list for MATLAB cell arrays / object arrays
            - numpy.ndarray for numeric arrays / time-series signals
            - Python scalars for numeric/string scalars

    Important design choices:
        - We keep numeric arrays as numpy arrays (do NOT convert to Python lists),
          because downstream code expects vectorized operations and slicing.
        - We *do* convert object arrays (MATLAB cells) into Python lists of
          recursively converted elements.
        - Strings may appear as bytes; we decode to UTF-8 when possible.
    """

    # -------------------------------------------------------------------------
    # 0) Base cases: None or simple Python scalars
    # -------------------------------------------------------------------------
    if mat_struct is None:
        return None

    # Bytes sometimes appear from MAT files; decode to text if possible.
    if isinstance(mat_struct, (bytes, bytearray)):
        try:
            return mat_struct.decode("utf-8")
        except Exception:
            # If decoding fails, return the raw bytes to avoid data loss.
            return mat_struct

    # Python-native scalar types can pass through unchanged.
    if isinstance(mat_struct, (str, int, float, bool)):
        return mat_struct

    # numpy scalar (e.g., np.float64) -> Python scalar (float/int) where possible
    try:
        import numpy as np

        if isinstance(mat_struct, np.generic):
            return mat_struct.item()
    except Exception:
        # If numpy isn't available or conversion fails, fall through.
        pass

    # -------------------------------------------------------------------------
    # 1) MATLAB struct objects: look for `_fieldnames`
    # -------------------------------------------------------------------------
    # With loadmat(..., struct_as_record=False), MATLAB structs become objects
    # that expose a list of field names on `_fieldnames`, and values as attributes.
    if hasattr(mat_struct, "_fieldnames"):
        out: Dict[str, object] = {}
        for fname in getattr(mat_struct, "_fieldnames", []):
            # Fetch MATLAB field value (may itself be struct/array/cell)
            try:
                raw_val = getattr(mat_struct, fname)
            except Exception:
                raw_val = None
            # Recursively convert
            out[str(fname)] = mat_struct_to_python(raw_val)
        return out

    # -------------------------------------------------------------------------
    # 2) Dictionaries are already Python-native, but recurse on values
    # -------------------------------------------------------------------------
    if isinstance(mat_struct, dict):
        return {str(k): mat_struct_to_python(v) for k, v in mat_struct.items()}

    # -------------------------------------------------------------------------
    # 3) Lists / tuples: recurse elementwise
    # -------------------------------------------------------------------------
    if isinstance(mat_struct, (list, tuple)):
        return [mat_struct_to_python(x) for x in mat_struct]

    # -------------------------------------------------------------------------
    # 4) numpy arrays: handle numeric vs object vs structured arrays
    # -------------------------------------------------------------------------
    try:
        import numpy as np

        if isinstance(mat_struct, np.ndarray):
            # 4a) Structured arrays / record arrays: dtype.names lists field names
            #     Common when MATLAB structs are loaded in certain ways.
            if mat_struct.dtype.names is not None:
                names = list(mat_struct.dtype.names)

                # If this is an array of records, return list[dict] (one dict per record)
                # If it's a single record, return a dict directly.
                if mat_struct.shape == ():  # scalar record (0-d)
                    return {n: mat_struct_to_python(mat_struct[n]) for n in names}

                records: List[Dict[str, object]] = []
                for rec in mat_struct:
                    records.append({n: mat_struct_to_python(rec[n]) for n in names})
                return records

            # 4b) Object arrays: typically MATLAB cell arrays -> convert to list
            if mat_struct.dtype == object:
                # Flatten 0-d or 1-element arrays to avoid awkward nesting
                if mat_struct.shape == ():
                    return mat_struct_to_python(mat_struct.item())

                # Convert to a Python list with recursive conversion.
                # Keep original order (row-major).
                flat = mat_struct.ravel()
                converted = [mat_struct_to_python(x) for x in flat]

                # If it's effectively a vector, return list
                if mat_struct.ndim == 1:
                    return converted

                # If it's multi-dimensional, reshape list-of-values back into nested lists
                # to preserve the original cell array shape.
                # Example: (m,n) -> list of m rows each containing n entries.
                if mat_struct.ndim == 2:
                    m, n = mat_struct.shape
                    return [converted[i * n:(i + 1) * n] for i in range(m)]

                # For ndim > 2, return the flat list (rare in our use case)
                return converted

            # 4c) Numeric arrays: keep as numpy.ndarray (best for time-series)
            #     squeeze_me=True already removes many singleton dimensions.
            return mat_struct

        # numpy record scalar (np.void) can appear when indexing structured arrays
        if isinstance(mat_struct, np.void) and getattr(mat_struct.dtype, "names", None):
            return {n: mat_struct_to_python(mat_struct[n]) for n in mat_struct.dtype.names}

    except Exception:
        # If numpy isn't available or something unexpected happens, fall through.
        pass

    # -------------------------------------------------------------------------
    # 5) Fallback: return object as-is
    # -------------------------------------------------------------------------
    # This ensures we never crash during conversion; unknown types can be handled
    # upstream (or flagged by schema validation) without losing raw data.
    return mat_struct
# =============================================================================
# Initial-guess inference helpers (XLSX does not contain Lp0/B0/sigma0/theta0)
# =============================================================================

def _normalize_component_key(names: object) -> Tuple[str, ...]:
    """Normalize component name(s) into a tuple[str,...] suitable for dict keys.

    Inputs:
        names:
            - None
            - a single string (e.g., 'KCl')
            - a list/tuple/np.ndarray of strings (e.g., ['NaCl','MgSO4'])

    Output:
        key:
            Tuple of cleaned component names in the original order.
            Example: ('NaCl', 'MgSO4')
    """
    if names is None:
        return tuple()
    if isinstance(names, str):
        s = names.strip()
        return (s,) if s else tuple()
    # Best-effort for MATLAB arrays or python sequences
    try:
        seq = list(names)
    except Exception:
        return tuple()
    out: List[str] = []
    for x in seq:
        if x is None:
            continue
        xs = str(x).strip()
        if xs:
            out.append(xs)
    return tuple(out)


def infer_initial_guess_db_from_mat_files(mat_files: Iterable[Path]) -> Dict[Tuple[str, ...], Dict[str, object]]:
    """Infer a component-keyed initial-guess database from legacy MATLAB data_stru files.

    Why this exists:
        Your legacy workflow stores initial guesses inside each MATLAB experiment struct
        (data_stru.data_config.Lp0, B0, sigma0, theta0). Your XLSX workbooks typically do not.
        This function builds a reusable lookup table so XLSX experiments can be assigned
        sensible starting guesses *without* hardcoding them.

    Inputs:
        mat_files:
            Iterable of .mat Paths. Each file is assumed to contain ONE experiment.

    Output:
        db:
            Dict keyed by component tuple, e.g.:
                ('KCl',) -> {'Lp0': 12.3, 'B0': 0.08, 'sigma0': 0.5, 'theta0': np.array([...])}

            Additionally, a global fallback is stored under key ('__GLOBAL__',):
                ('__GLOBAL__',) -> median guesses across all experiments.

    Notes:
        - Files that do not contain 'data_stru' are skipped.
        - If theta0 lengths vary across experiments, we compute the elementwise median
          using the most common theta0 length in that component group.
    """
    # Collect raw samples per component-key
    samples: Dict[Tuple[str, ...], List[Dict[str, object]]] = {}

    for p in mat_files:
        if p.suffix.lower() != '.mat':
            continue

        try:
            blob = loadmat(p, squeeze_me=True, struct_as_record=False)
        except Exception:
            # Skip unreadable files (corrupt or non-MATLAB v5)
            continue

        if 'data_stru' not in blob:
            continue

        ds = blob['data_stru']
        cfg = getattr(ds, 'data_config', None)
        if cfg is None:
            continue

        key = _normalize_component_key(getattr(cfg, 'namec', None))

        rec: Dict[str, object] = {}
        rec['Lp0'] = getattr(cfg, 'Lp0', None)
        rec['B0'] = getattr(cfg, 'B0', None)
        rec['sigma0'] = getattr(cfg, 'sigma0', None)

        # theta0 should be numeric; store as 1D numpy array when possible
        th = getattr(cfg, 'theta0', None)
        if th is not None:
            try:
                rec['theta0'] = np.asarray(th, dtype=float).reshape(-1)
            except Exception:
                rec['theta0'] = None
        else:
            rec['theta0'] = None

        samples.setdefault(key, []).append(rec)

    def _median_ignore_none(vals: List[object]) -> Optional[float]:
        """Median ignore none.

        Args:
            vals: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        nums: List[float] = []
        for v in vals:
            if v is None:
                continue
            try:
                fv = float(v)
            except Exception:
                continue
            if np.isnan(fv):
                continue
            nums.append(fv)
        if not nums:
            return None
        return float(np.median(nums))

    def _median_theta0(arrs: List[Optional[np.ndarray]]) -> Optional[np.ndarray]:
        """Median theta0.

        Args:
            arrs: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        arrs2 = [a for a in arrs if isinstance(a, np.ndarray) and a.size > 0]
        if not arrs2:
            return None
        # Choose the most common length
        lengths = [a.size for a in arrs2]
        L = max(set(lengths), key=lengths.count)
        stack = np.vstack([a[:L] for a in arrs2 if a.size >= L])
        return np.median(stack, axis=0)

    # Aggregate into a compact db with medians
    db: Dict[Tuple[str, ...], Dict[str, object]] = {}

    all_recs: List[Dict[str, object]] = []
    for key, recs in samples.items():
        all_recs.extend(recs)

        db[key] = {
            'Lp0': _median_ignore_none([r.get('Lp0') for r in recs]),
            'B0': _median_ignore_none([r.get('B0') for r in recs]),
            'sigma0': _median_ignore_none([r.get('sigma0') for r in recs]),
            'theta0': _median_theta0([r.get('theta0') for r in recs]),
        }

    # Global fallback across all experiments
    if all_recs:
        db[('__GLOBAL__',)] = {
            'Lp0': _median_ignore_none([r.get('Lp0') for r in all_recs]),
            'B0': _median_ignore_none([r.get('B0') for r in all_recs]),
            'sigma0': _median_ignore_none([r.get('sigma0') for r in all_recs]),
            'theta0': _median_theta0([r.get('theta0') for r in all_recs]),
        }

    return db


# TODO(remove): legacy helper is currently unused in the unified v24 workflow.
# def infer_initial_guess_db_from_mat_folder(mat_root: Path) -> Dict[Tuple[str, ...], Dict[str, object]]:
#     """Convenience wrapper: infer initial guesses from all .mat files under a folder."""
#     mat_files = [p for p in mat_root.rglob('*.mat') if '__MACOSX' not in str(p)]
#     return infer_initial_guess_db_from_mat_files(mat_files)


def apply_initial_guesses_from_db(
    exp: ExperimentalData,
    db: Dict[Tuple[str, ...], Dict[str, object]],
) -> None:
    """Populate exp.Lp0/B0/sigma0/theta0 from an initial-guess DB when possible.

    Inputs:
        exp:
            ExperimentalData whose component_names have ideally already been populated.
        db:
            Output of infer_initial_guess_db_from_mat_files(...) or similar.

    Output:
        None (mutates exp in-place)

    Behavior:
        - Try exact match on tuple(exp.component_names)
        - Else fall back to ('__GLOBAL__',) if available
        - Only sets a field if exp.<field> is currently None
    """
    key = _normalize_component_key(exp.component_names)
    rec = db.get(key, None)

    if rec is None:
        rec = db.get(('__GLOBAL__',), None)

    if rec is None:
        return

    if exp.Lp0 is None and rec.get('Lp0') is not None:
        exp.Lp0 = float(rec['Lp0'])
    if exp.B0 is None and rec.get('B0') is not None:
        exp.B0 = float(rec['B0'])
    if exp.sigma0 is None and rec.get('sigma0') is not None:
        exp.sigma0 = float(rec['sigma0'])
    if exp.theta0 is None and rec.get('theta0') is not None:
        exp.theta0 = np.asarray(rec['theta0'], dtype=float).reshape(-1)



# =============================================================================
# 9) Smoke test (optional): run one MAT and one XLSX experiment through the loader
# =============================================================================

if __name__ == "__main__":
    """
    Minimal smoke test for the file-agnostic loader.

    This block is OPTIONAL and safe:
        - It does not modify model logic
        - It does not write files
        - It simply loads one MAT and one XLSX experiment and prints a short summary

    How to use:
        python experiment_dataload_OOP_v13_step5_smoketest.py
    """
    from pathlib import Path

    # ---- Sample files (edit these paths if needed) ----
    mat_path = Path("/mnt/data/data_stru-dataset501.1.mat")
    xlsx_path = Path("/mnt/data/NF270_MC2.xlsx")

    print("\n==============================")
    print("SMOKE TEST: MAT loader")
    print("==============================")
    try:
        exp_mat, (ok_mat, issues_mat) = load_experiment(str(mat_path), selector=None)
        print(f"Loaded MAT file: {mat_path.name}")
        print(f"ok_fatal: {ok_mat}")
        print(f"n_vials: {len(exp_mat.vials)}")
        if exp_mat.vials:
            v0 = exp_mat.vials[0]
            print(f"First vial number: {v0.number}")
            print(f"len(time_s): {0 if v0.time_s is None else len(v0.time_s)}")
            print(f"len(retentate_signal): {0 if v0.retentate_signal is None else len(v0.retentate_signal)}")
            print(f"retentate_signal_name: {v0.retentate_signal_name}")
            print(f"retentate_signal_units: {v0.retentate_signal_units}")
        if issues_mat:
            print("Issues:")
            for sev, msg in issues_mat:
                print(f"  - [{sev}] {msg}")
    except Exception as e:
        print(f"MAT smoke test failed: {type(e).__name__}: {e}")

    print("\n==============================")
    print("SMOKE TEST: XLSX loader (one sheet)")
    print("==============================")
    try:
        import pandas as pd
        xls = pd.ExcelFile(xlsx_path)
        sheet0 = xls.sheet_names[0]
        exp_xlsx, (ok_xlsx, issues_xlsx) = load_experiment(str(xlsx_path), selector=sheet0)
        print(f"Loaded XLSX file: {xlsx_path.name} | sheet: {sheet0}")
        print(f"ok_fatal: {ok_xlsx}")
        print(f"n_vials: {len(exp_xlsx.vials)}")
        print(f"mode: {exp_xlsx.mode}")
        print(f"component_names: {exp_xlsx.component_names}")
        print(f"M_F0_g: {exp_xlsx.M_F0_g} | M_O_g: {exp_xlsx.M_O_g}")
        if exp_xlsx.vials:
            v0 = exp_xlsx.vials[0]
            print(f"First vial number: {v0.number}")
            print(f"len(time_s): {0 if v0.time_s is None else len(v0.time_s)}")
            print(f"len(retentate_signal): {0 if v0.retentate_signal is None else len(v0.retentate_signal)}")
            print(f"retentate_signal_name: {v0.retentate_signal_name}")
            print(f"retentate_signal_units: {v0.retentate_signal_units}")
            print(f"cV_assay_value: {v0.cV_assay_value} {v0.cV_assay_units}")
        if issues_xlsx:
            print("Issues:")
            for sev, msg in issues_xlsx:
                print(f"  - [{sev}] {msg}")
    except Exception as e:
        print(f"XLSX smoke test failed: {type(e).__name__}: {e}")
class ExperimentMode(str, Enum):
    """ExperimentMode.

    Container class used by the unified loader/model workflow.

    """
    DATA = "DATA"
    LAG = "Lag"
    OVERFLOW = "Overflow"


class RunMode(str, Enum):
    """RunMode.

    Container class used by the unified loader/model workflow.

    """
    SIMULATION = "SIMULATION"
    ESTIMATION = "ESTIMATION"


class BForm(str, Enum):
    """BForm.

    Container class used by the unified loader/model workflow.

    """
    SINGLE = "single"
    PERVIAL = "pervial"
    CONVECTION = "convection"


BFormType = Union[str, float]


@dataclass
class ModelOptions:
    """ModelOptions.

    Container class used by the unified loader/model workflow.

    """
    mode: ExperimentMode = ExperimentMode.DATA
    run_mode: RunMode = RunMode.ESTIMATION
    b_form: BFormType = BForm.SINGLE.value
    nfe: int = 300
    fd_scheme: str = "BACKWARD"
    time_scaled_end: float = 1.0
    use_advanced_xlsx_transport_thermo: bool = True
    membrane_fixed_charge_mM: float = -44.0
    partition_coeff_cation: float = 0.85
    partition_coeff_anion: float = 0.85
    osmotic_factor: float = 1.0
    use_sigma_logit_transform: bool = True
    sigma_eps: float = 1e-3


@dataclass
class ParameterGuess:
    """ParameterGuess.

    Container class used by the unified loader/model workflow.

    """
    Lp: float
    sigma: float
    B: Optional[Union[float, Dict[int, float]]] = None
    beta_0: Optional[float] = None
    beta_1: Optional[float] = None
    beta_2: Optional[float] = None
    beta_3: Optional[float] = None
    S0: Optional[float] = None
    S: Optional[float] = None


R_BAR_CM3_PER_UMOL_K = 8.314e-5
NU_CM2_S = 8.927e-3
CELL_DIAMETER_CM = 2.2860
RPM_DEFAULT = 350.0
MH_ML = 0.25
DIFFUSIVITY_CM2_S = {
    "K": 1.960e-5,
    "Na": 1.334e-5,
    "Li": 1.03e-5,
    "Mg": 0.706e-5,
    "Ca": 0.792e-5,
    "Co": 0.72e-5,
    "La": 0.62e-5,
}


def _first_non_nan(x: Optional[np.ndarray], default: float = 1e-6) -> float:
    """First non nan.

    Args:
        x: Parameter description.
        default: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    if x is None:
        return float(default)
    arr = np.asarray(x, dtype=float).reshape(-1)
    for val in arr:
        if not np.isnan(val):
            return float(val)
    return float(default)


def _salt_key(component_names: Optional[List[str]]) -> Optional[str]:
    """Salt key.

    Args:
        component_names: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    if not component_names:
        return None
    s = str(component_names[0]).strip()
    letters = "".join(c for c in s if c.isalpha())
    for k in DIFFUSIVITY_CM2_S:
        if k in letters:
            return k
    return None


def _infer_primary_ions(exp: ExperimentalData) -> Tuple[str, int, str, int]:
    """Infer one effective cation/anion pair and valencies from component names."""
    cat = "Na"
    an = "Cl"
    z_cat = 1
    z_an = -1

    names = exp.component_names or []
    for name in names:
        try:
            cation, _, anion, _ = _parse_simple_salt_formula(str(name))
            zc = int(ION_CHARGE.get(cation, +1))
            za = int(ION_CHARGE.get(anion, -1))
            cat = cation
            an = anion
            z_cat = zc if zc > 0 else 1
            z_an = za if za < 0 else -1
            break
        except Exception:
            continue
    return cat, z_cat, an, z_an


def compute_mass_transfer_coeff(exp: ExperimentalData) -> float:
    """Compute mass transfer coeff.

    Args:
        exp: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    b = CELL_DIAMETER_CM
    nu = NU_CM2_S
    v = RPM_DEFAULT / 60.0 * np.pi * b

    d_key = _salt_key(exp.component_names)
    if d_key is None:
        raise ValueError(f"Cannot infer diffusivity from component names: {exp.component_names}")
    D = DIFFUSIVITY_CM2_S[d_key]
    return 0.23 * v ** 0.57 * D ** 0.67 / (nu ** 0.24 * b ** 0.43)


def _scaled_times(exp: ExperimentalData) -> Tuple[Dict[int, float], Dict[int, float], float]:
    """Scaled times.

    Args:
        exp: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    t_delay = float(exp.vials[0].time_s[0])
    ti = {}
    tf = {}
    for i, v in enumerate(exp.vials, start=1):
        ti[i] = float(v.time_s[0]) - t_delay
        tf[i] = float(v.time_s[-1]) - t_delay
    return ti, tf, t_delay


def _default_guess(exp: ExperimentalData) -> ParameterGuess:
    """Default guess.

    Args:
        exp: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    return ParameterGuess(
        Lp=float(exp.Lp0 if exp.Lp0 is not None else 5.0),
        sigma=float(exp.sigma0 if exp.sigma0 is not None else 0.9),
        B=float(exp.B0 if exp.B0 is not None else 0.5),
        beta_0=float(exp.B0 if exp.B0 is not None else 0.5),
        beta_1=0.5,
        beta_2=0.0,
        beta_3=0.0,
        S0=0.0,
        S=0.0,
    )


def theta_component_list(m: pyo.ConcreteModel, options: ModelOptions) -> List[pyo.ComponentData]:
    """Theta component list.

    Args:
        m: Parameter description.
        options: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    b_form = str(options.b_form).lower()
    comps: List[pyo.ComponentData] = [m.Lp]

    if b_form == BForm.SINGLE.value:
        comps.append(m.B)
    elif b_form == BForm.PERVIAL.value:
        comps.extend(m.B[n] for n in m.n_vial)
    elif b_form == BForm.CONVECTION.value:
        comps.append(m.beta_0)
        comps.append(m.beta_1)
    else:
        raise ValueError(f"Unsupported b_form '{options.b_form}' for theta list.")

    if hasattr(m, "sigma_logit"):
        comps.append(m.sigma_logit)
    else:
        comps.append(m.sigma)
    return comps


def theta_names(m: pyo.ConcreteModel, options: ModelOptions) -> List[str]:
    """Theta names.

    Args:
        m: Parameter description.
        options: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    return [c.name for c in theta_component_list(m, options)]


def _augment_theta_with_sigma(theta_obj: object, *, options: ModelOptions) -> object:
    """Attach physical sigma to theta output when sigma_logit parametrization is used."""
    if not options.use_sigma_logit_transform:
        return theta_obj
    sigma_eps = float(options.sigma_eps)
    scale = (1.0 - 2.0 * sigma_eps)

    if isinstance(theta_obj, pd.Series) and "sigma_logit" in theta_obj.index:
        logit_val = float(theta_obj["sigma_logit"])
        theta_obj = theta_obj.copy()
        theta_obj["sigma"] = sigma_eps + scale / (1.0 + np.exp(-logit_val))
        return theta_obj

    if isinstance(theta_obj, dict) and "sigma_logit" in theta_obj:
        out = dict(theta_obj)
        logit_val = float(out["sigma_logit"])
        out["sigma"] = sigma_eps + scale / (1.0 + np.exp(-logit_val))
        return out

    return theta_obj


@contextmanager
def _temporary_ipopt_opt(
    *,
    max_iter: int = 12000,
    tol: float = 1e-7,
    acceptable_tol: float = 1e-5,
) -> Iterable[None]:
    """Temporarily create/override local ipopt.opt for covariance solves."""
    opt_path = Path.cwd() / "ipopt.opt"
    backup_text: Optional[str] = None
    had_existing = opt_path.exists()
    if had_existing:
        backup_text = opt_path.read_text()

    opt_path.write_text(
        "\n".join(
            [
                f"max_iter {int(max_iter)}",
                f"tol {float(tol)}",
                f"acceptable_tol {float(acceptable_tol)}",
                "mu_strategy adaptive",
                "bound_push 1e-8",
                "bound_frac 1e-8",
                "print_level 0",
            ]
        )
        + "\n"
    )
    try:
        yield
    finally:
        if had_existing and backup_text is not None:
            opt_path.write_text(backup_text)
        else:
            try:
                opt_path.unlink()
            except FileNotFoundError:
                pass


def model_construct_inter_v23(
    exp: ExperimentalData,
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
) -> pyo.ConcreteModel:
    """Build the dynamic diafiltration model (v23)."""
    if guess is None:
        guess = _default_guess(exp)

    if exp.delP_bar is None or exp.Temp_K is None or exp.Am_cm2 is None or exp.rho_g_cm3 is None:
        raise ValueError("Missing required experiment inputs (delP_bar, Temp_K, Am_cm2, rho_g_cm3).")
    if exp.M_F0_g is None:
        raise ValueError("exp.M_F0_g is required for model construction.")
    if exp.C_D_value is None:
        raise ValueError("exp.C_D_value is required for model construction.")

    delP = float(exp.delP_bar)
    T = float(exp.Temp_K)
    Am = float(exp.Am_cm2)
    rho = float(exp.rho_g_cm3)
    ni = int(exp.num_components if exp.num_components is not None else 1)

    M_F0 = float(exp.M_F0_g)
    M_O = float(exp.M_O_g if exp.M_O_g is not None else 0.0)
    cD = float(exp.C_D_value)
    C_F0 = float(exp.C_F0_value if exp.C_F0_value is not None else _first_non_nan(exp.vials[0].retentate_signal, 1e-6))

    N_VIAL = len(exp.vials)
    TI_dict, TF_dict, _ = _scaled_times(exp)
    k = compute_mass_transfer_coeff(exp)

    m = pyo.ConcreteModel()
    m.n_vial = pyo.RangeSet(1, N_VIAL)
    m.tau = ContinuousSet(bounds=(0.0, float(options.time_scaled_end)))
    m.tf = pyo.Set(initialize=[TF_dict[i] for i in m.n_vial])
    m.ti = pyo.Set(initialize=[TI_dict[i] for i in m.n_vial])

    b_form = str(options.b_form).lower()
    advanced_xlsx_transport = (
        bool(options.use_advanced_xlsx_transport_thermo)
        and exp.source == SourceType.XLSX
    )
    sim_opt = options.run_mode == RunMode.SIMULATION

    m.cD = pyo.Var(domain=pyo.NonNegativeReals, initialize=cD)
    m.cD.fix(cD)
    m.C_H0 = pyo.Param(initialize=1e-6, mutable=True)

    sigma_guess = float(guess.sigma)
    sigma_eps = float(options.sigma_eps)
    sigma_eps = min(max(sigma_eps, 1e-8), 0.49)

    if sim_opt:
        m.Lp = pyo.Param(initialize=float(guess.Lp), mutable=True)
        m.sigma = pyo.Param(initialize=sigma_guess, mutable=True)
    else:
        m.Lp = pyo.Var(bounds=(0.5, 50.0), initialize=float(guess.Lp))
        if options.use_sigma_logit_transform:
            sig0 = min(max(sigma_guess, sigma_eps + 1e-8), 1.0 - sigma_eps - 1e-8)
            frac = (sig0 - sigma_eps) / (1.0 - 2.0 * sigma_eps)
            frac = min(max(frac, 1e-8), 1.0 - 1e-8)
            logit0 = float(np.log(frac / (1.0 - frac)))
            m.sigma_logit = pyo.Var(bounds=(-20.0, 20.0), initialize=logit0)

            def _sigma_rule(mm):
                """Sigma rule.

                Args:
                    mm: Parameter description.

                Returns:
                    object: Computed value or expression.

                """
                return sigma_eps + (1.0 - 2.0 * sigma_eps) / (1.0 + pyo.exp(-mm.sigma_logit))
            m.sigma = pyo.Expression(rule=_sigma_rule)
        else:
            m.sigma = pyo.Var(bounds=(sigma_eps, 1.0 - sigma_eps), initialize=min(max(sigma_guess, sigma_eps), 1.0 - sigma_eps))

    if b_form == BForm.SINGLE.value:
        if sim_opt:
            m.B = pyo.Param(initialize=float(guess.B if guess.B is not None else 0.5), mutable=True)
        else:
            m.B = pyo.Var(bounds=(1e-6, 30.0), initialize=float(guess.B if guess.B is not None else 0.5))
    elif b_form == BForm.PERVIAL.value:
        if isinstance(guess.B, dict):
            default_b = list(guess.B.values())[0] if guess.B else 0.5
            b_init = {i: float(guess.B.get(i, default_b)) for i in range(1, N_VIAL + 1)}
        else:
            b_init = {i: float(guess.B if guess.B is not None else 0.5) for i in range(1, N_VIAL + 1)}
        if sim_opt:
            m.B = pyo.Param(m.n_vial, initialize=b_init, mutable=True)
        else:
            m.B = pyo.Var(m.n_vial, bounds=(1e-6, 30.0), initialize=b_init)
    elif b_form == BForm.CONVECTION.value:
        beta0 = float(guess.beta_0 if guess.beta_0 is not None else 2.0)
        beta1 = float(guess.beta_1 if guess.beta_1 is not None else 0.5)
        if sim_opt:
            m.beta_0 = pyo.Param(initialize=beta0, mutable=True)
            m.beta_1 = pyo.Param(initialize=beta1, mutable=True)
        else:
            m.beta_0 = pyo.Var(bounds=(1 + 1e-6, 50.0), initialize=beta0)
            m.beta_1 = pyo.Var(bounds=(0.0, 1.0), initialize=beta1)
        m.H = pyo.Var(m.n_vial, m.tau, initialize=0.5)
    else:
        raise ValueError(f"Unsupported b_form: {options.b_form}")

    if options.mode in (ExperimentMode.LAG, ExperimentMode.OVERFLOW):
        if sim_opt:
            m.S0 = pyo.Param(initialize=float(guess.S0 if guess.S0 is not None else 0.0), mutable=True)
            m.S = pyo.Param(initialize=float(guess.S if guess.S is not None else 0.0), mutable=True)
        else:
            m.S0 = pyo.Var(initialize=float(guess.S0 if guess.S0 is not None else 0.0))
            m.S = pyo.Var(initialize=float(guess.S if guess.S is not None else 0.0))

    if options.mode != ExperimentMode.DATA:
        m.mF = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=M_F0)
        m.dmF = DerivativeVar(m.mF, wrt=m.tau)

    m.cF = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=C_F0)
    m.cIn = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=C_F0)
    m.cH = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=1e-6)
    m.mV = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=1e-6)
    m.cVmV = pyo.Var(m.n_vial, m.tau, initialize=1e-12)
    m.cV = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=1e-6)

    m.Jw = pyo.Var(m.n_vial, m.tau)
    m.Js = pyo.Var(m.n_vial, m.tau)
    if b_form == BForm.CONVECTION.value:
        m.Js_exp = pyo.Var(m.n_vial, m.tau, bounds=(1 + 1e-6, 1e4))

    if advanced_xlsx_transport:
        _, z_cat, _, z_an = _infer_primary_ions(exp)
        m.z_cat = pyo.Param(initialize=float(z_cat))
        m.z_an = pyo.Param(initialize=float(z_an))
        m.fixed_charge = pyo.Param(initialize=float(options.membrane_fixed_charge_mM))
        m.H_cat = pyo.Param(initialize=float(options.partition_coeff_cation))
        m.H_an = pyo.Param(initialize=float(options.partition_coeff_anion))
        m.osmotic_factor = pyo.Param(initialize=float(options.osmotic_factor))

        m.cM = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=C_F0)
        m.cP = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=max(cD, 1e-6))

    m.dcF = DerivativeVar(m.cF, wrt=m.tau)
    m.dcH = DerivativeVar(m.cH, wrt=m.tau)
    m.dmV = DerivativeVar(m.mV, wrt=m.tau)
    m.dcVmV = DerivativeVar(m.cVmV, wrt=m.tau)

    Tauf = float(options.time_scaled_end)

    def _tf_scale(n: int) -> float:
        """Tf scale.

        Args:
            n: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        return (TF_dict[n] - TI_dict[n]) / Tauf

    if options.mode == ExperimentMode.DATA:
        def ode_cF_rule(mm, n, t):
            """Ode cF rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.dcF[n, t] == Am * rho / M_F0 * (mm.cD * mm.Jw[n, t] - mm.Js[n, t]) * _tf_scale(n)
        m.ode_cF = pyo.Constraint(m.n_vial, m.tau, rule=ode_cF_rule)
    else:
        def ode_mF_rule(mm, n, t):
            """Ode mF rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.dmF[n, t] == (-mm.S0 - Am * rho * mm.Jw[n, t]) * _tf_scale(n)
        m.ode_mF = pyo.Constraint(m.n_vial, m.tau, rule=ode_mF_rule)

        def ode_cF_rule(mm, n, t):
            """Ode cF rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.dcF[n, t] == (1 / mm.mF[n, t]) * ((mm.cF[n, t] - mm.cD) * mm.S0 + Am * rho * (mm.cF[n, t] * mm.Jw[n, t] - mm.Js[n, t])) * _tf_scale(n)
        m.ode_cF = pyo.Constraint(m.n_vial, m.tau, rule=ode_cF_rule)

    def ode_cH_rule(mm, n, t):
        """Ode cH rule.

        Args:
            mm: Parameter description.
            n: Parameter description.
            t: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        return mm.dcH[n, t] == Am * rho / MH_ML * (mm.Js[n, t] - mm.cH[n, t] * mm.Jw[n, t]) * _tf_scale(n)
    m.ode_cH = pyo.Constraint(m.n_vial, m.tau, rule=ode_cH_rule)

    def ode_mV_rule(mm, n, t):
        """Ode mV rule.

        Args:
            mm: Parameter description.
            n: Parameter description.
            t: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        return mm.dmV[n, t] == mm.Jw[n, t] * Am * rho * _tf_scale(n)
    m.ode_mV = pyo.Constraint(m.n_vial, m.tau, rule=ode_mV_rule)

    def ode_cVmV_rule(mm, n, t):
        """Ode cVmV rule.

        Args:
            mm: Parameter description.
            n: Parameter description.
            t: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        return mm.dcVmV[n, t] == mm.Jw[n, t] * mm.cH[n, t] * Am * rho * _tf_scale(n)
    m.ode_cVmV = pyo.Constraint(m.n_vial, m.tau, rule=ode_cVmV_rule)

    if advanced_xlsx_transport:
        def eqn_cM_film_rule(mm, n, t):
            """Eqn cM film rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.cM[n, t] == (mm.cF[n, t] - mm.cH[n, t]) * pyo.exp(mm.Jw[n, t] / k) + mm.cH[n, t]
        m.eqn_cM_film = pyo.Constraint(m.n_vial, m.tau, rule=eqn_cM_film_rule)

        def eqn_cIn_rule(mm, n, t):
            """Eqn cIn rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.cIn[n, t] == mm.cM[n, t]
        m.eqn_cIn = pyo.Constraint(m.n_vial, m.tau, rule=eqn_cIn_rule)
    else:
        def eqn_cIn_rule(mm, n, t):
            """Eqn cIn rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.cIn[n, t] == (mm.cF[n, t] - mm.cH[n, t]) * pyo.exp(mm.Jw[n, t] / k) + mm.cH[n, t]
        m.eqn_cIn = pyo.Constraint(m.n_vial, m.tau, rule=eqn_cIn_rule)

    if advanced_xlsx_transport:
        def delta_pi_rule(mm, n, t):
            """Delta pi rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.osmotic_factor * mm.sigma * R_BAR_CM3_PER_UMOL_K * T * (mm.cM[n, t] - mm.cP[n, t]) * ni
        m.delta_pi = pyo.Expression(m.n_vial, m.tau, rule=delta_pi_rule)

        def eqn_Jw_rule(mm, n, t):
            """Eqn Jw rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.Jw[n, t] * 36000 == mm.Lp * (delP - mm.delta_pi[n, t])
        m.eqn_Jw = pyo.Constraint(m.n_vial, m.tau, rule=eqn_Jw_rule)
    else:
        def eqn_Jw_rule(mm, n, t):
            """Eqn Jw rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.Jw[n, t] * 36000 == mm.Lp * (delP - (mm.cIn[n, t] - mm.cH[n, t]) * ni * mm.sigma * R_BAR_CM3_PER_UMOL_K * T)
        m.eqn_Jw = pyo.Constraint(m.n_vial, m.tau, rule=eqn_Jw_rule)

    def eqn_Js_rule(mm, n, t):
        """Eqn Js rule.

        Args:
            mm: Parameter description.
            n: Parameter description.
            t: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        delta_c = (mm.cM[n, t] - mm.cP[n, t]) if advanced_xlsx_transport else (mm.cIn[n, t] - mm.cH[n, t])
        if b_form == BForm.SINGLE.value:
            return mm.Js[n, t] * 10000 == mm.B * delta_c
        if b_form == BForm.PERVIAL.value:
            return mm.Js[n, t] * 10000 == mm.B[n] * delta_c
        return mm.Js[n, t] == mm.Jw[n, t] * mm.H[n, t] * (mm.cIn[n, t] * mm.Js_exp[n, t] - mm.cH[n, t]) / (mm.Js_exp[n, t] - 1)
    m.eqn_Js = pyo.Constraint(m.n_vial, m.tau, rule=eqn_Js_rule)

    if b_form == BForm.CONVECTION.value:
        def eqn_Js_exp_rule(mm, n, t):
            """Eqn Js exp rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.Js_exp[n, t] == pyo.exp(mm.Jw[n, t] / mm.beta_0 * 10000)
        m.eqn_Js_exp = pyo.Constraint(m.n_vial, m.tau, rule=eqn_Js_exp_rule)

        def eqn_H_rule(mm, n, t):
            """Eqn H rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.H[n, t] == mm.beta_1
        m.eqn_H = pyo.Constraint(m.n_vial, m.tau, rule=eqn_H_rule)

    def eqn_cV_rule(mm, n, t):
        """Eqn cV rule.

        Args:
            mm: Parameter description.
            n: Parameter description.
            t: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        return mm.mV[n, t] * mm.cV[n, t] == mm.cVmV[n, t]
    m.eqn_cV = pyo.Constraint(m.n_vial, m.tau, rule=eqn_cV_rule)

    if advanced_xlsx_transport:
        def eqn_cP_link_rule(mm, n, t):
            """Eqn cP link rule.

            Args:
                mm: Parameter description.
                n: Parameter description.
                t: Parameter description.

            Returns:
                object: Computed value or expression.

            """
            return mm.cP[n, t] == mm.cH[n, t]
        m.eqn_cP_link = pyo.Constraint(m.n_vial, m.tau, rule=eqn_cP_link_rule)

        # Advanced thermodynamic diagnostics are kept as expressions to avoid
        # overconstraining the dynamic estimation model.
        m.cA = pyo.Expression(m.n_vial, m.tau, rule=lambda mm, n, t: (abs(mm.z_cat) / abs(mm.z_an)) * mm.cF[n, t])
        m.cA_m = pyo.Expression(m.n_vial, m.tau, rule=lambda mm, n, t: (abs(mm.z_cat) / abs(mm.z_an)) * mm.cM[n, t] + mm.fixed_charge / abs(mm.z_an))
        m.cA_p = pyo.Expression(m.n_vial, m.tau, rule=lambda mm, n, t: (abs(mm.z_cat) / abs(mm.z_an)) * mm.cP[n, t])
        m.Js_cation = pyo.Expression(m.n_vial, m.tau, rule=lambda mm, n, t: mm.Js[n, t])
        m.Js_anion = pyo.Expression(m.n_vial, m.tau, rule=lambda mm, n, t: (abs(mm.z_cat) / abs(mm.z_an)) * mm.Js[n, t])
        m.partition_cat_residual = pyo.Expression(m.n_vial, m.tau, rule=lambda mm, n, t: mm.cM[n, t] - mm.H_cat * mm.cF[n, t])
        m.partition_an_residual = pyo.Expression(m.n_vial, m.tau, rule=lambda mm, n, t: mm.cA_m[n, t] - mm.H_an * mm.cA[n, t])
        m.electroneutral_ret_residual = pyo.Expression(m.n_vial, m.tau, rule=lambda mm, n, t: mm.z_cat * mm.cF[n, t] + mm.z_an * mm.cA[n, t])
        m.electroneutral_mem_residual = pyo.Expression(m.n_vial, m.tau, rule=lambda mm, n, t: mm.z_cat * mm.cM[n, t] + mm.z_an * mm.cA_m[n, t] + mm.fixed_charge)
        m.electroneutral_perm_residual = pyo.Expression(m.n_vial, m.tau, rule=lambda mm, n, t: mm.z_cat * mm.cP[n, t] + mm.z_an * mm.cA_p[n, t])

    def cF_linking_rule(mm, n):
        """CF linking rule.

        Args:
            mm: Parameter description.
            n: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        if n == mm.n_vial.last():
            return pyo.Constraint.Skip
        return mm.cF[n, Tauf] == mm.cF[n + 1, 0.0]
    m.cF_linking = pyo.Constraint(m.n_vial, rule=cF_linking_rule)

    def cH_linking_rule(mm, n):
        """CH linking rule.

        Args:
            mm: Parameter description.
            n: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        if n == mm.n_vial.last():
            return pyo.Constraint.Skip
        return mm.cH[n, Tauf] == mm.cH[n + 1, 0.0]
    m.cH_linking = pyo.Constraint(m.n_vial, rule=cH_linking_rule)

    m.con_boundary = pyo.ConstraintList()
    m.con_boundary.add(m.cH[1, 0.0] == m.C_H0)
    m.con_boundary.add(m.cF[1, 0.0] == _first_non_nan(exp.vials[0].retentate_signal, default=1e-6))
    m.con_boundary.add(m.cVmV[1, 0.0] == 1e-12)
    m.con_boundary.add(m.mV[1, 0.0] == 1e-6)
    if options.mode != ExperimentMode.DATA:
        m.con_boundary.add(m.mF[1, 0.0] == M_F0)

    if options.mode != ExperimentMode.DATA and options.run_mode == RunMode.ESTIMATION:
        m.eqn_S = pyo.Constraint(expr=m.mF[m.n_vial.last(), Tauf] - M_F0 == M_O)

    return m


def apply_discretization(m: pyo.ConcreteModel, *, nfe: int = 300, scheme: str = "BACKWARD") -> None:
    """Apply discretization.

    Args:
        m: Parameter description.
        nfe: Parameter description.
        scheme: Parameter description.

    """
    pyo.TransformationFactory("dae.finite_difference").apply_to(m, nfe=nfe, scheme=scheme)


def _nearest_tau_value(tau_values: Sequence[float], target: float) -> float:
    """Nearest tau value.

    Args:
        tau_values: Parameter description.
        target: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    arr = np.asarray(tau_values, dtype=float)
    return float(arr[np.argmin(np.abs(arr - target))])


def attach_weighted_least_squares_objective(
    m: pyo.ConcreteModel,
    exp: ExperimentalData,
    *,
    sigma_mass_g: float = 0.01,
    sigma_cV_rel: float = 0.03,
    sigma_cF_rel: float = 0.003,
) -> None:
    """Attach weighted least squares objective.

    Args:
        m: Parameter description.
        exp: Parameter description.
        sigma_mass_g: Parameter description.
        sigma_cV_rel: Parameter description.
        sigma_cF_rel: Parameter description.

    """
    if not list(m.tau):
        raise RuntimeError("Model must be discretized before adding measurement objective.")

    tau_vals = sorted(float(t) for t in list(m.tau))
    TI_dict, TF_dict, t_delay = _scaled_times(exp)

    expr = 0.0
    for n in m.n_vial:
        vial = exp.vials[int(n) - 1]
        t_meas = np.asarray(vial.time_s, dtype=float) - t_delay
        dur = TF_dict[int(n)] - TI_dict[int(n)]
        if dur <= 0:
            continue
        t_scaled = (t_meas - TI_dict[int(n)]) / dur

        if vial.mass_g is not None:
            mass = np.asarray(vial.mass_g, dtype=float)
            for tm, y in zip(t_scaled, mass):
                if np.isnan(y):
                    continue
                tau_m = _nearest_tau_value(tau_vals, float(tm))
                expr += ((m.mV[n, tau_m] - float(y)) / sigma_mass_g) ** 2

        cf = np.asarray(vial.retentate_signal, dtype=float) if vial.retentate_signal is not None else np.array([])
        for tm, y in zip(t_scaled, cf):
            if np.isnan(y) or y == 0:
                continue
            tau_m = _nearest_tau_value(tau_vals, float(tm))
            expr += ((m.cF[n, tau_m] - float(y)) / (sigma_cF_rel * abs(float(y)))) ** 2

        if isinstance(vial.cV_avg, (int, float, np.number)) and not np.isnan(float(vial.cV_avg)) and float(vial.cV_avg) != 0:
            tau_m = _nearest_tau_value(tau_vals, 1.0)
            y = float(vial.cV_avg)
            expr += ((m.cV[n, tau_m] - y) / (sigma_cV_rel * abs(y))) ** 2

    m.FirstStageCost = pyo.Expression(expr=0.0)
    m.SecondStageCost = pyo.Expression(expr=expr)
    m.Total_Cost_Objective = pyo.Objective(expr=m.FirstStageCost + m.SecondStageCost, sense=pyo.minimize)


def label_parmest_and_doe_suffixes(
    m: pyo.ConcreteModel,
    exp: ExperimentalData,
    options: ModelOptions,
    *,
    design_vars: Optional[List[pyo.ComponentData]] = None,
    sigma_mass_g: float = 0.01,
    sigma_cV_rel: float = 0.03,
    sigma_cF_rel: float = 0.003,
) -> None:
    """Attach experiment_outputs / measurement_error / unknown_parameters / experiment_inputs suffixes."""
    if not list(m.tau):
        raise RuntimeError("Model must be discretized before labeling outputs.")

    tau_vals = sorted(float(t) for t in list(m.tau))
    TI_dict, TF_dict, t_delay = _scaled_times(exp)

    m.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    m.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)

    for n in m.n_vial:
        vial = exp.vials[int(n) - 1]
        t_meas = np.asarray(vial.time_s, dtype=float) - t_delay
        dur = TF_dict[int(n)] - TI_dict[int(n)]
        if dur <= 0:
            continue
        t_scaled = (t_meas - TI_dict[int(n)]) / dur

        if vial.mass_g is not None:
            mass = np.asarray(vial.mass_g, dtype=float)
            for tm, y in zip(t_scaled, mass):
                if np.isnan(y):
                    continue
                tau_m = _nearest_tau_value(tau_vals, float(tm))
                key = m.mV[n, tau_m]
                m.experiment_outputs[key] = float(y)
                m.measurement_error[key] = float(sigma_mass_g)

        cf = np.asarray(vial.retentate_signal, dtype=float) if vial.retentate_signal is not None else np.array([])
        for tm, y in zip(t_scaled, cf):
            if np.isnan(y):
                continue
            tau_m = _nearest_tau_value(tau_vals, float(tm))
            key = m.cF[n, tau_m]
            m.experiment_outputs[key] = float(y)
            m.measurement_error[key] = float(max(sigma_cF_rel * abs(float(y)), 1e-8))

        if isinstance(vial.cV_avg, (int, float, np.number)) and not np.isnan(float(vial.cV_avg)):
            tau_m = _nearest_tau_value(tau_vals, 1.0)
            key = m.cV[n, tau_m]
            y = float(vial.cV_avg)
            m.experiment_outputs[key] = y
            m.measurement_error[key] = float(max(sigma_cV_rel * abs(y), 1e-8))

    m.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    for comp in theta_component_list(m, options):
        m.unknown_parameters[comp] = pyo.value(comp)

    m.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    if design_vars is None:
        design_vars = [m.cD]
    for dv in design_vars:
        m.experiment_inputs[dv] = None


def model_construct_for_parmest_v23(
    exp: ExperimentalData,
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
) -> pyo.ConcreteModel:
    """Convenience builder used by parmest/doe Experiment wrapper."""
    m = model_construct_inter_v23(exp, options=options, guess=guess)
    apply_discretization(m, nfe=options.nfe, scheme=options.fd_scheme)
    attach_weighted_least_squares_objective(m, exp)
    label_parmest_and_doe_suffixes(m, exp, options)
    return m


class DiafiltrationExperimentV23(ParmestExperiment):
    """ParmEst/DoE Experiment wrapper for one loaded ExperimentalData object."""

    def __init__(
        self,
        exp: ExperimentalData,
        options: ModelOptions,
        guess: Optional[ParameterGuess] = None,
    ):
        """Init.

        Args:
            exp: Parameter description.
            options: Parameter description.
            guess: Parameter description.

        Returns:
            object: Computed value or expression.

        """
        super().__init__(model=None)
        self.exp = exp
        self.options = options
        self.guess = guess

    def get_labeled_model(self):
        """Get labeled model.

        Returns:
            object: Computed value or expression.

        """
        self.model = model_construct_for_parmest_v23(self.exp, self.options, self.guess)
        return self.model


def build_experiment_list_v23(
    experiments: List[ExperimentalData],
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
) -> List[DiafiltrationExperimentV23]:
    """Build experiment list v23.

    Args:
        experiments: Parameter description.
        options: Parameter description.
        guess: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    return [DiafiltrationExperimentV23(exp=e, options=options, guess=guess) for e in experiments]


# TODO(remove): v23 parmest wrapper is unused after migration to v24 SSE_weighted workflow.
# def estimate_parameters_with_parmest_v23(
#     experiments: List[ExperimentalData],
#     options: ModelOptions,
#     guess: Optional[ParameterGuess] = None,
#     *,
#     calc_cov: bool = True,
#     cov_n: Optional[int] = None,
#     solver: str = "ipopt",
#     solver_options: Optional[Dict[str, object]] = None,
#     tee: bool = False,
# ) -> Dict[str, object]:
#     """Run parameter estimation via pyomo.contrib.parmest.Estimator."""
#     exp_list = build_experiment_list_v23(experiments, options, guess)
#     estimator = Estimator(exp_list, tee=tee, solver_options=solver_options)
#
#     try:
#         out = estimator.theta_est(solver=solver, calc_cov=calc_cov, cov_n=cov_n)
#         result: Dict[str, object] = {"raw": out, "estimator": estimator}
#     except RuntimeError as err:
#         if solver == "ipopt" and "Unknown solver interface" in str(err):
#             result = _estimate_parameters_ipopt_fallback(
#                 experiments=experiments,
#                 options=options,
#                 guess=guess,
#                 solver_options=solver_options,
#                 tee=tee,
#             )
#             result["warning"] = (
#                 "ParmEst solver interface failed; used ipopt fallback on a single labeled experiment model. "
#                 "Covariance is not computed in fallback mode."
#             )
#             return result
#         raise
#
#     if isinstance(out, tuple):
#         if len(out) >= 1:
#             result["objective"] = out[0]
#         if len(out) >= 2:
#             result["theta"] = out[1]
#         if len(out) >= 3:
#             result["returned_values"] = out[2]
#         if len(out) >= 4:
#             result["covariance"] = out[3]
#
#     return result


def _estimate_parameters_ipopt_fallback(
    experiments: List[ExperimentalData],
    options: ModelOptions,
    guess: Optional[ParameterGuess],
    solver_options: Optional[Dict[str, object]],
    tee: bool,
) -> Dict[str, object]:
    """Estimate parameters ipopt fallback.

    Args:
        experiments: Parameter description.
        options: Parameter description.
        guess: Parameter description.
        solver_options: Parameter description.
        tee: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    if len(experiments) != 1:
        raise RuntimeError(
            "ipopt fallback currently supports one experiment at a time. "
            "Use Estimator.theta_est with ipopt for multi-experiment solves."
        )

    exp = experiments[0]
    m = model_construct_for_parmest_v24(exp, options=options, guess=guess)
    solver = pyo.SolverFactory("ipopt")
    if solver_options:
        for k, v in solver_options.items():
            solver.options[k] = v
    # Build weighted SSE objective from ParmEst-style suffix labels.
    # This mirrors obj_function="SSE_weighted" for single-experiment fallback.
    def _weighted_sse_expr():
        terms = []
        for out_v in m.experiment_outputs.keys():
            y_obs = float(m.experiment_outputs[out_v])
            sigma = float(m.measurement_error[out_v])
            sigma_eff = sigma if abs(sigma) > 1e-12 else 1.0
            terms.append(((out_v - y_obs) / sigma_eff) ** 2)
        return sum(terms) if terms else 0.0

    m._ipopt_fallback_wsse = pyo.Objective(expr=_weighted_sse_expr(), sense=pyo.minimize)

    raw = solver.solve(m, tee=tee)

    theta_vals = {c.name: float(pyo.value(c)) for c in theta_component_list(m, options)}
    theta_vals = _augment_theta_with_sigma(theta_vals, options=options)
    objective = float(pyo.value(m._ipopt_fallback_wsse))

    return {
        "raw": raw,
        "objective": objective,
        "theta": pd.Series(theta_vals),
        "model": m,
    }


# TODO(remove): not referenced by current runner (FIM is computed directly by DoE API).
# def covariance_to_fim(covariance: Union[np.ndarray, object]) -> np.ndarray:
#     """Covariance to fim."""
#     cov = np.asarray(covariance, dtype=float)
#     if cov.ndim != 2 or cov.shape[0] != cov.shape[1]:
#         raise ValueError("Covariance must be a square 2D matrix.")
#     return np.linalg.pinv(cov)


def d_optimality(fim: Union[np.ndarray, object], *, log_scale: bool = True) -> float:
    """D optimality.

    Args:
        fim: Parameter description.
        log_scale: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    fim_arr = np.asarray(fim, dtype=float)
    sign, logdet = np.linalg.slogdet(fim_arr)
    if sign <= 0:
        return float("-inf") if log_scale else 0.0
    return float(logdet) if log_scale else float(np.exp(logdet))


def run_doe_with_pyomo_v23(
    experiment: DiafiltrationExperimentV23,
    *,
    fd_formula: str = "central",
    step: float = 0.001,
    objective_option: str = "determinant",
    solver_name: str = "ipopt",
    tee: bool = False,
) -> Dict[str, object]:
    """Run DoE FIM analysis/optimization using pyomo.contrib.doe."""
    solver = pyo.SolverFactory(solver_name)
    doe_obj = DesignOfExperiments(
        experiment=experiment,
        fd_formula=fd_formula,
        step=step,
        objective_option=objective_option,
        solver=solver,
        tee=tee,
    )
    # First attempt a direct FIM computation at nominal settings.
    # If this fails, fall back to full DoE run.
    results = None
    fim = None
    try:
        fim = doe_obj.compute_FIM(method="sequential")
        results = {"method": "compute_FIM_sequential"}
    except Exception:
        results = doe_obj.run_doe()
        fim = doe_obj.get_FIM()
    return {
        "doe": doe_obj,
        "results": results,
        "fim": np.asarray(fim, dtype=float),
        "d_opt_logdet": d_optimality(fim, log_scale=True),
    }


# TODO(remove): legacy callback factory is unused in v24 unified workflow.
# def parmest_callback_factory_v23(
#     experiments: List[ExperimentalData],
#     options: ModelOptions,
#     guess: Optional[ParameterGuess] = None,
# ) -> Callable[[int], pyo.ConcreteModel]:
#     """Legacy-compatible callback factory retained for quick usage."""
#
#     def _callback(idx: int) -> pyo.ConcreteModel:
#         exp = experiments[idx]
#         return model_construct_for_parmest_v23(exp, options=options, guess=guess)
#
#     return _callback


def model_construct_inter_v24(
    exp: ExperimentalData,
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
) -> pyo.ConcreteModel:
    """Model construct inter v24.

    Args:
        exp: Parameter description.
        options: Parameter description.
        guess: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    return model_construct_inter_v23(exp=exp, options=options, guess=guess)


def model_construct_for_parmest_v24(
    exp: ExperimentalData,
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
) -> pyo.ConcreteModel:
    """Builder for ParmEst/DoE in v24.

    For ParmEst built-in objectives (SSE/SSE_weighted), do not attach a custom
    objective here. ParmEst will construct FirstStageCost/SecondStageCost and
    Total_Cost_Objective internally from model suffixes.
    """
    # Build the continuous-time process model (transport + thermo relations).
    m = model_construct_inter_v24(exp=exp, options=options, guess=guess)
    # Discretize the DAE model so ParmEst/DoE can solve an NLP.
    apply_discretization(m, nfe=options.nfe, scheme=options.fd_scheme)
    # Attach measured outputs, errors, unknown parameter labels, and design inputs.
    label_parmest_and_doe_suffixes(m, exp, options)
    # Return labeled/discretized model without a custom objective; ParmEst injects SSE objective.
    return m


class DiafiltrationExperimentV24(DiafiltrationExperimentV23):
    """v24 wrapper with xlsx-only advanced transport/thermo model support."""

    def get_labeled_model(self):
        """Get labeled model.

        Returns:
            object: Computed value or expression.

        """
        self.model = model_construct_for_parmest_v24(self.exp, self.options, self.guess)
        return self.model


def build_experiment_list_v24(
    experiments: List[ExperimentalData],
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
) -> List[DiafiltrationExperimentV24]:
    # Convert raw ExperimentalData list into ParmEst-compatible experiment wrappers.
    """Build experiment list v24.

    Args:
        experiments: Parameter description.
        options: Parameter description.
        guess: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    return [DiafiltrationExperimentV24(exp=e, options=options, guess=guess) for e in experiments]


def estimate_parameters_with_parmest_v24(
    experiments: List[ExperimentalData],
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
    *,
    calc_cov: bool = True,
    cov_n: Optional[int] = None,
    solver: str = "ipopt",
    solver_options: Optional[Dict[str, object]] = None,
    tee: bool = False,
) -> Dict[str, object]:
    """ParmEst run aligned with Pyomo 6.9.5 API (theta_est + cov_est)."""
    # Build experiment wrapper list expected by pyomo.contrib.parmest.Estimator.
    exp_list = build_experiment_list_v24(experiments, options, guess)
    # Use built-in SSE_weighted objective to unlock standard covariance workflows.
    estimator = Estimator(
        exp_list,
        obj_function="SSE_weighted",
        tee=tee,
        solver_options=solver_options,
    )

    warning_msgs: List[str] = []
    try:
        # Solve the parameter estimation problem with chosen solver (ipopt by default).
        out = estimator.theta_est(solver=solver)
        # Store raw estimator return and estimator handle for downstream use.
        result: Dict[str, object] = {"raw": out, "estimator": estimator}
    except RuntimeError as err:
        # Pyomo 6.9.5 ParmEst may still require ef_ipopt token in theta_est.
        if solver == "ipopt" and "Unknown solver in Q_Opt=ipopt" in str(err):
            out = estimator.theta_est(solver="ef_ipopt")
            result = {"raw": out, "estimator": estimator}
            warning_msgs.append(
                "ParmEst theta_est in this environment requires solver token 'ef_ipopt'; "
                "mapped user solver 'ipopt' to ParmEst 'ef_ipopt'."
            )
        # Last-resort path: direct single-experiment solve when ParmEst front-end fails.
        elif (
            solver == "ipopt"
            and len(experiments) == 1
            and ("Unknown solver in Q_Opt=ef_ipopt" in str(err) or "No executable found for solver" in str(err))
        ):
            result = _estimate_parameters_ipopt_fallback(
                experiments=experiments,
                options=options,
                guess=guess,
                solver_options=solver_options,
                tee=tee,
            )
            warning_msgs.append(
                "ParmEst theta_est solver interface failed; used direct ipopt solve on labeled model for one experiment."
            )
            if warning_msgs:
                result["warning"] = " ".join(warning_msgs)
            return result
        else:
            raise

    if isinstance(out, tuple):
        # theta_est typically returns (objective, theta, ...)
        if len(out) >= 1:
            result["objective"] = out[0]
        if len(out) >= 2:
            # Add physical sigma to output when sigma_logit transform is used.
            result["theta"] = _augment_theta_with_sigma(out[1], options=options)
        if len(out) >= 3:
            result["returned_values"] = out[2]

    if calc_cov:
        # Try covariance methods from most direct to more specialized.
        cov_attempts = [
            ("finite_difference", {"solver": "ipopt", "step": 1e-4}),
            ("reduced_hessian", {"solver": "ipopt"}),
            ("automatic_differentiation_kaug", {"solver": "ipopt"}),
        ]
        # Collect method-level errors if all methods fail.
        cov_errs: List[str] = []
        for method, kwargs in cov_attempts:
            try:
                # Use temporary ipopt.opt to raise iteration budget and improve robustness.
                with _temporary_ipopt_opt():
                    cov = estimator.cov_est(method=method, **kwargs)
                # Save first successful covariance estimate and method.
                result["covariance"] = cov
                result["covariance_method"] = method
                break
            except Exception as err:
                cov_errs.append(f"{method}: {type(err).__name__}: {err}")
        if "covariance" not in result:
            # Return actionable diagnostics when covariance fails completely.
            result["covariance_warning"] = (
                "Covariance not computed. Parameter estimates were obtained with "
                "ParmEst built-in SSE_weighted, but covariance solves failed for all methods. "
                "Typical causes are local ill-conditioning, parameters near bounds, or "
                f"failed perturbed NLP solves. Attempts: {' | '.join(cov_errs)}"
            )
        if cov_n is not None:
            # Keep backward-compatibility note for callers that still pass cov_n.
            result["cov_n_note"] = "cov_n is deprecated in Pyomo 6.9.5 and is ignored in v24."

    if warning_msgs:
        result["warning"] = " ".join(warning_msgs)

    return result


def run_doe_with_pyomo_v24(
    experiment: DiafiltrationExperimentV24,
    *,
    fd_formula: str = "central",
    step: float = 0.001,
    objective_option: str = "determinant",
    solver_name: str = "ipopt",
    tee: bool = False,
) -> Dict[str, object]:
    """Run doe with pyomo v24.

    Args:
        experiment: Parameter description.
        fd_formula: Parameter description.
        step: Parameter description.
        objective_option: Parameter description.
        solver_name: Parameter description.
        tee: Parameter description.

    Returns:
        object: Computed value or expression.

    """
    return run_doe_with_pyomo_v23(
        experiment=experiment,
        fd_formula=fd_formula,
        step=step,
        objective_option=objective_option,
        solver_name=solver_name,
        tee=tee,
    )
