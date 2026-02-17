
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Created on Wed Dec 03 07:12:33 2025
@author: Keshav Kasturi Rangan, Alex Dowling
===============================================================================

File-agnostic experimental data loading pseudocode.

Goal:
    Load ONE experiment from either:
        - an Excel workbook (.xlsx) where each sheet is one experiment, OR
        - a MATLAB file (.mat) that stores an experiment struct
    and store the result in explicit, readable dataclass fields 
    (do not dump in dictionary "extras").

Design choices:
    1) Conductivity is stored as a signal and not silently converted to 
       concentration. The Shdlovsly/MSA will be used later to convert 
       conductivity -> concentration.
    2) Two light classes are used: ExperimentalData and VialData.
    3) Validation flags are issued without crashing; missing values remain None.

Notes:
    1) This is pseudocode: functions like read_excel_sheet(), read_mat_file(), 
       parse_time_series_table() represent format-specific I/O/parsing steps 
       you will implement with pandas/openpyxl/scipy.io later.
    2) The loader keeps ExperimentalData/VialData "data-only".
       Process-model code may wrap/extend ExperimentalData with convenience methods
       (e.g., vial time bounds) without changing loader responsibilities.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional, List, Tuple, Union

import numpy as np




# =============================================================================
# 1) Enums (avoid string typos like "xlxs" instead of "xlsx")
# =============================================================================

class SourceType(str, Enum):
    """
    Enumerates supported experiment data sources.

    Values:
        XLSX: Excel workbook input (.xlsx)
        MAT:  MATLAB file input (.mat)
    """
    XLSX = "xlsx"
    MAT = "mat"


# =============================================================================
# 2) Dataclasses: exactly TWO light containers, explicit fields only
# =============================================================================

@dataclass
class VialData:
    """
    Stores all data for one vial segment (one contiguous interval between vial 
                                          swaps).

    Core fields:
        number:
            1-based vial index (1, 2, 3, ...). None until assigned by loader.

        time_s:
            Time stamps [s] for this vial segment. None until loaded.

        mass_g:
            Mass measurements [g], aligned with time_s. None if not available.

        retentate_signal:
            Continuous retentate measurement aligned with time_s.
            For Excel: retentate conductivity (uS/cm) stored honestly.
            For MAT:   whatever cF_exp represents in the struct; label honestly 
                       if known.

        retentate_signal_name / retentate_signal_units:
            Human-readable identifier and units for retentate_signal.

    Optional time-series fields (all aligned with time_s if present):
        pressure_psi, retentate_temp_C, permeate_temp_C, permeate_signal, 
        vial_swap_flag.

    Optional assays (scalars per vial):
        cV_assay_value (+ units), sample_volume_mL, nitric_acid_volume_mL, etc.

    IMPORTANT typing decision:
        cV_avg is a UNION because legacy datasets may store it as:
            - scalar (float)
            - array (np.ndarray) e.g., NaN-filled time-series with one final 
                value
        Using a precise Union keeps schema readable while preserving 
        flexibility.
    """
    number: Optional[int] = None

    # Time-aligned arrays (None until loaded)
    time_s: Optional[np.ndarray] = None
    mass_g: Optional[np.ndarray] = None

    # Honest retentate signal (conductivity unless explicitly converted later)
    retentate_signal: Optional[np.ndarray] = None
    retentate_signal_name: Optional[str] = None
    retentate_signal_units: Optional[str] = None

    # Optional additional time-series signals (aligned with time_s)
    pressure_psi: Optional[np.ndarray] = None
    retentate_temp_C: Optional[np.ndarray] = None
    permeate_temp_C: Optional[np.ndarray] = None

    permeate_signal: Optional[np.ndarray] = None
    permeate_signal_name: Optional[str] = None
    permeate_signal_units: Optional[str] = None

    vial_swap_flag: Optional[np.ndarray] = None

    # Optional vial assays (scalars)
    cV_assay_value: Optional[float] = None
    cV_assay_units: Optional[str] = None

    sample_volume_mL: Optional[float] = None
    nitric_acid_volume_mL: Optional[float] = None

    # Legacy-compatible permeate concentration storage (scalar OR array)
    cV_avg: Optional[Union[float, np.ndarray]] = None


@dataclass
class ExperimentalData:
    """
    Stores one experiment (one Excel sheet OR one MATLAB struct) as explicit 
    fields.

    Metadata/provenance:
        source, filename, sheet_name

    Identity:
        dataset_id, mode

    Measurement behavior:
        continuous_retentate (True for your Excel experiments)

    Optional operating conditions/config:
        delP_bar, Temp_K, Am_cm2, rho_g_cm3, etc.

    Vials:
        vials is an ordered list of VialData objects (empty until loaded).

    IMPORTANT typing decision:
        theta0 is ALWAYS numeric in your workflow, so it is typed as:
            Optional[np.ndarray]
        (None means "not provided", otherwise it is a numeric vector/array.)
    """
    # Provenance
    source: Optional[SourceType] = None
    filename: Optional[str] = None
    sheet_name: Optional[str] = None

    # Experiment identity
    dataset_id: Optional[float] = None
    mode: Optional[str] = None

    # Retentate measurement behavior
    continuous_retentate: Optional[bool] = None

    # Optional operating conditions / config
    delP_bar: Optional[float] = None
    Temp_K: Optional[float] = None
    Am_cm2: Optional[float] = None
    rho_g_cm3: Optional[float] = None

    # Optional feed/dialysate information
    M_F0_g: Optional[float] = None
    C_F0_value: Optional[float] = None
    C_F0_units: Optional[str] = None
    M_O_g: Optional[float] = None
    C_D_value: Optional[float] = None
    C_D_units: Optional[str] = None

    # Optional component bookkeeping
    component_names: Optional[List[str]] = None
    num_components: Optional[int] = None

    # Optional initial guesses / priors
    Lp0: Optional[float] = None
    B0: Optional[float] = None
    sigma0: Optional[float] = None

    # theta0 is always numeric (None if missing)
    theta0: Optional[np.ndarray] = None

    # Optional calibration
    salt1_calib_x: Optional[np.ndarray] = None
    salt1_calib_y: Optional[np.ndarray] = None
    salt2_calib_x: Optional[np.ndarray] = None
    salt2_calib_y: Optional[np.ndarray] = None

    # Vial list
    vials: List[VialData] = field(default_factory=list)


# =============================================================================
# 3) Shared helpers (reused by both MAT and XLSX loaders)
# =============================================================================

def is_missing_scalar(x: object) -> bool:
    """
    Returns True if x should be treated as missing.

    Inputs:
        x: scalar-like value

    Output:
        True if missing (None, NaN, empty string), else False
    """
    if x is None:
        return True
    if isinstance(x, float) and np.isnan(x):
        return True
    if isinstance(x, str) and x.strip() == "":
        return True
    return False


def as_float_array(x: object) -> np.ndarray:
    """
    Converts x into a 1D numpy float array (best-effort).

    Inputs:
        x: list/tuple/np.ndarray/pandas Series/MATLAB array-like

    Output:
        1D float array (empty array if x is None)
    """
    if x is None:
        return np.array([], dtype=float)
    arr = np.asarray(x, dtype=float)
    return arr.reshape(-1)


def as_int_array(x: object) -> np.ndarray:
    """
    Converts x into a 1D numpy int array (best-effort).
    Useful for swap flags (0/1).

    Inputs:
        x: list/tuple/np.ndarray/pandas Series

    Output:
        1D int array (empty array if x is None)
    """
    if x is None:
        return np.array([], dtype=int)
    arr = np.asarray(x, dtype=float)
    arr = np.nan_to_num(arr, nan=0.0)
    return arr.astype(int).reshape(-1)


def segment_indices_by_swap(swap_flags: np.ndarray) -> List[Tuple[int, int]]:
    """
    Converts a 0/1 swap flag array into vial segments (start,end) for Python 
    slicing.

    Convention:
        swap_flags == 1 row is included in the *current* vial (end = idx+1).

    Inputs:
        swap_flags: 1D array of 0/1 markers

    Output:
        List of (start, end) segments representing [start:end]
    """
    swap_idx = np.where(swap_flags == 1)[0]
    segments: List[Tuple[int, int]] = []

    start = 0
    for idx in swap_idx:
        end = idx + 1
        segments.append((start, end))
        start = end

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
    Creates and populates a VialData object by slicing time-aligned arrays 
    [a:b].

    Inputs:
        number: vial number (1-based)
        a, b: slice bounds [a:b]
        time_all: full time array
        mass_all: full mass array (optional)
        ret_signal_all: full retentate signal array (conductivity)
        ret_signal_name/units: honest labels for retentate signal
        swap_all, pressure_all, ret_temp_all, perm_temp_all: optional arrays 
                                                            aligned to time_all
        perm_signal_all: optional permeate signal array aligned to time_all

    Output:
        VialData with fields set for this vial segment
    """
    v = VialData()  # create an empty VialData container
    v.number = number  # store the vial index

    # Slice core aligned arrays
    v.time_s = time_all[a:b]  # time stamps for this vial
    v.retentate_signal = ret_signal_all[a:b]  # retentate conductivity for this
                                              # vial

    # Store honest labels so downstream code knows what "retentate_signal" is
    v.retentate_signal_name = ret_signal_name
    v.retentate_signal_units = ret_signal_units

    # Mass is optional; only slice if it exists
    if mass_all is not None:
        v.mass_g = mass_all[a:b]

    # Store optional time-series signals if present
    if swap_all is not None:
        v.vial_swap_flag = swap_all[a:b]

    if pressure_all is not None:
        v.pressure_psi = pressure_all[a:b]
    if ret_temp_all is not None:
        v.retentate_temp_C = ret_temp_all[a:b]
    if perm_temp_all is not None:
        v.permeate_temp_C = perm_temp_all[a:b]

    # Optional permeate signal (e.g., permeate conductivity)
    if perm_signal_all is not None:
        v.permeate_signal = perm_signal_all[a:b]
        v.permeate_signal_name = perm_signal_name
        v.permeate_signal_units = perm_signal_units

    return v


def validate_experiment(exp: ExperimentalData) -> Tuple[bool, List[Tuple[str, str]]]:
    """
    Validates basic structural integrity and flags missing values without
    raising exceptions.

    Additional checks (still non-throwing):
        - FATAL if a vial's time_s goes backwards (non-monotonic time).
        - REQUIRED if exp.theta0 exists but is not numeric/castable to float.

    Inputs:
        exp: loaded experiment container

    Outputs:
        (ok_fatal, issues)
            ok_fatal: False if experiment is structurally unusable
            issues: list of (severity, message) tuples
    """
    issues: List[Tuple[str, str]] = []
    ok_fatal = True

    # FATAL: must have vials
    if exp.vials is None or len(exp.vials) == 0:
        issues.append(("FATAL", "No vials found (exp.vials is empty)."))
        return False, issues

    # FATAL: each vial must have time and retentate signal, and lengths must 
    # match
    for v in exp.vials:
        if v.time_s is None or len(v.time_s) == 0:
            issues.append(("FATAL", f"Vial {v.number}: missing/empty time_s."))
            ok_fatal = False
        if v.retentate_signal is None or len(v.retentate_signal) == 0:
            issues.append(("FATAL", 
                           f"Vial {v.number}: missing/empty retentate_signal."))
            ok_fatal = False
        if (v.time_s is not None) and (v.retentate_signal is not None) and (len(v.time_s) != len(v.retentate_signal)):
            issues.append(("FATAL", f"Vial {v.number}: time and retentate_signal lengths differ."))
            ok_fatal = False


        # FATAL: time should not go backwards within a vial segment
        # Many solvers / interpolation routines assume time is sorted.
        if v.time_s is not None and len(v.time_s) > 1:
            if np.any(np.diff(v.time_s) < 0):
                issues.append(("FATAL",
                               f"Vial {v.number}: time_s is not monotonic non-decreasing."))
                ok_fatal = False

    # REQUIRED (flag only): common operating conditions (many models need these)
    if is_missing_scalar(exp.delP_bar):
        issues.append(("REQUIRED", "Missing exp.delP_bar (pressure)."))
    if is_missing_scalar(exp.Temp_K):
        issues.append(("REQUIRED", "Missing exp.Temp_K (temperature)."))
    if is_missing_scalar(exp.Am_cm2):
        issues.append(("REQUIRED", "Missing exp.Am_cm2 (membrane area)."))


    # REQUIRED (flag only): theta0 must be numeric if provided
    # By schema: theta0 is either None (not provided) or a numeric vector.
    if exp.theta0 is not None:
        try:
            _ = np.asarray(exp.theta0, dtype=float)
        except Exception:
            issues.append(("REQUIRED",
                           "exp.theta0 exists but is not numeric/castable to float array."))

    return ok_fatal, issues


# =============================================================================
# 4) File-agnostic loader entrypoint
# =============================================================================

def load_experiment(path: str, selector: object) -> Tuple[ExperimentalData, Tuple[bool, List[Tuple[str, str]]]]:
    """
    Loads ONE experiment from either .xlsx or .mat into ExperimentalData fields.

    Inputs:
        path: file path to .xlsx or .mat
        selector:
            - XLSX: sheet name (str) or sheet index (int)
            - MAT:  struct key (str) or None meaning default "data_stru"

    Outputs:
        (exp, validation):
            exp: loaded experiment with explicit fields
            validation: (ok_fatal, issues)
    """
    file_path = Path(path)
    kind = detect_source_type(file_path)

    if kind == SourceType.XLSX:
        exp = load_from_xlsx(file_path, sheet_selector=selector)
    elif kind == SourceType.MAT:
        exp = load_from_mat(file_path, struct_key=selector)
    else:
        raise ValueError(f"Unsupported file type: {file_path.suffix}")

    validation = validate_experiment(exp)
    return exp, validation


def detect_source_type(file_path: Path) -> Optional[SourceType]:
    """
    Detects SourceType from file extension.

    Inputs:
        file_path: input path

    Output:
        SourceType.XLSX, SourceType.MAT, or None if unknown
    """
    suffix = file_path.suffix.lower()
    if suffix in [".xlsx", ".xls"]:
        return SourceType.XLSX
    if suffix == ".mat":
        return SourceType.MAT
    return None



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


# =============================================================================
# 5) XLSX loader (one sheet = one experiment)
# =============================================================================

def load_from_xlsx(xlsx_path: Path, sheet_selector: object) -> ExperimentalData:
    """
    Loads one experiment from an Excel workbook sheet into ExperimentalData fields.

    Inputs:
        xlsx_path: path to the Excel workbook
        sheet_selector: sheet name (str) or sheet index (int)

    Output:
        ExperimentalData populated from the sheet
    """
    # Read the sheet (format-specific I/O)
    sheet = read_excel_sheet(xlsx_path, sheet_selector)

    # Parse the time-series table (format-specific parsing)
    ts = parse_time_series_table(sheet)

    # Extract required arrays using shared helpers
    time_all = as_float_array(ts["Time (s)"])
    mass_all = as_float_array(ts["Mass (g)"])
    swap_all = as_int_array(ts["Vial Swap"])

    # Honest retentate signal: conductivity
    ret_cond_all = as_float_array(ts["Retentate Cond @ Temp (uS/cm)"])
    ret_signal_name = "retentate_cond_uS_cm"
    ret_signal_units = "uS/cm"

    # Optional arrays: check column existence via a helper (avoid pandas-specific assumptions)
    pressure_all = as_float_array(ts["Pressure (psi)"]) if table_has_column(ts, "Pressure (psi)") else None
    ret_temp_all = as_float_array(ts["Retentate Temp"]) if table_has_column(ts, "Retentate Temp") else None
    perm_temp_all = as_float_array(ts["Permeate Temp"]) if table_has_column(ts, "Permeate Temp") else None

    perm_cond_all = as_float_array(ts["Permeate Cond @ Temp (uS/cm)"]) if table_has_column(ts, "Permeate Cond @ Temp (uS/cm)") else None
    perm_signal_name = "permeate_cond_uS_cm"
    perm_signal_units = "uS/cm"

    # Parse vial assay table (format-specific parsing)
    vial_table = parse_vial_data_table(sheet)

    # Segment the experiment into vials using swap flags (shared helper)
    segments = segment_indices_by_swap(swap_all)

    # Create experiment container
    exp = ExperimentalData()
    exp.source = SourceType.XLSX
    exp.filename = xlsx_path.name
    exp.sheet_name = get_sheet_name(sheet)
    exp.dataset_id = None
    exp.mode = None
    exp.continuous_retentate = True

    # Build VialData objects for each segment
    for i, (a, b) in enumerate(segments, start=1):
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

        # Attach vial assays (format-specific row lookup)
        assay_row = find_vial_row(vial_table, vial_number=i)
        if assay_row is not None:
            v.cV_assay_value = assay_row.get("ICP Salt 1 (mg/L)", None)
            v.cV_assay_units = "mg/L"
            v.sample_volume_mL = assay_row.get("Sample Volume (mL)", None)
            v.nitric_acid_volume_mL = assay_row.get("Nitric Acid Volume (mL)", None)

            # Keep legacy-compatible cV_avg as scalar by default
            v.cV_avg = v.cV_assay_value

        exp.vials.append(v)

    return exp



# =============================================================================
# MAT schema constants (expected struct keys)
# =============================================================================

# These constants define the expected *shape* of the MATLAB experiment struct (legacy).
# They are intended for explicit parsing/validation and are not fuzzy.
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


# =============================================================================
# 6) MAT loader (struct already contains vials)
# =============================================================================

def load_from_mat(mat_path: Path, struct_key: Optional[str]) -> ExperimentalData:
    """
    Loads one experiment from a MATLAB .mat file into ExperimentalData fields.

    Inputs:
        mat_path: path to .mat file
        struct_key: key name for experiment struct; if None, default "data_stru"

    Output:
        ExperimentalData populated from the struct
    """
    # Read .mat file (format-specific I/O)
    mat_blob = read_mat_file(mat_path)

    # Choose which struct to use
    key = struct_key if struct_key is not None else "data_stru"

    # Convert MATLAB struct to nested Python types (TRANSIENT only)
    d = mat_struct_to_python(mat_blob[key])

    # Create experiment container
    exp = ExperimentalData()
    exp.source = SourceType.MAT
    exp.filename = str(d.get("filename", mat_path.name))
    exp.sheet_name = None
    exp.dataset_id = d.get("dataset", None)
    exp.mode = d.get("mode", None)
    exp.continuous_retentate = bool(d.get("conductivity_cF", 1))

    # Optional operating conditions/config
    cfg = d.get("data_config", {})
    exp.delP_bar = cfg.get("delP", None)
    exp.Temp_K = cfg.get("Temp", None)
    exp.Am_cm2 = cfg.get("Am", None)
    exp.rho_g_cm3 = cfg.get("rho", None)

    # Optional numeric initial guess vector (if present)
    # Keep as None if not available; if present, ensure it becomes a numeric numpy array.
    theta0_raw = cfg.get("theta0", None)
    if theta0_raw is not None:
        exp.theta0 = np.asarray(theta0_raw, dtype=float).reshape(-1)

    # Vial list
    raw_vials = d.get("data_raw", [])

    # Honest labeling: unless proven otherwise, treat cF_exp as conductivity-derived signal
    ret_signal_name = "retentate_cond_uS_cm"
    ret_signal_units = "uS/cm"

    for raw in raw_vials:
        v = VialData()
        v.number = int(raw.get("number", len(exp.vials) + 1))

        v.time_s = as_float_array(raw.get("time", None))
        v.mass_g = as_float_array(raw.get("mass", None)) if raw.get("mass", None) is not None else None
        v.retentate_signal = as_float_array(raw.get("cF_exp", None))

        v.retentate_signal_name = ret_signal_name
        v.retentate_signal_units = ret_signal_units

        # Preserve cV_avg exactly as stored (scalar or array)
        cV = raw.get("cV_avg", None)
        if isinstance(cV, (float, int)) or cV is None:
            v.cV_avg = None if cV is None else float(cV)
        else:
            v.cV_avg = np.asarray(cV)

        exp.vials.append(v)

    return exp


# =============================================================================
# 7) Format-specific I/O/parsing stubs (to implement in real code)
# =============================================================================

def read_excel_sheet(xlsx_path: Path, sheet_selector: object) -> object:
    """Stub: return a sheet object for the requested selector."""
    raise NotImplementedError



def parse_time_series_table(sheet: object) -> object:
    """
    Stub: return a table-like object with named columns for time-series signals.

    Contract:
        - The returned table MUST support column access by canonical names, e.g.:
              ts["Time (s)"], ts["Mass (g)"], ts["Vial Swap"], ts["Retentate Cond @ Temp (uS/cm)"]

        - When reading Excel, this parser SHOULD:
            1) scan raw sheet rows to locate the header row,
            2) recognize headers using CANONICAL_REQUIRED_TS_COLS with whitelisted TS_COL_ALIASES,
            3) rename any recognized aliases to the canonical names so downstream code
               sees a stable schema, even if spreadsheets drift slightly.

        - This function does NOT:
            - validate units/values,
            - segment into vials,
            - convert conductivity to concentration.

    Inputs:
        sheet: format-specific "sheet handle" returned by read_excel_sheet(...)

    Output:
        table-like object (often a pandas DataFrame in real code)
    """
    raise NotImplementedError



def parse_vial_data_table(sheet: object) -> object:
    """Stub: return a table-like object with per-vial assay rows (e.g., 'Vial 1', ...)."""
    raise NotImplementedError


def find_vial_row(vial_table: object, vial_number: int) -> Optional[object]:
    """Stub: return row corresponding to 'Vial {vial_number}'."""
    raise NotImplementedError


def get_sheet_name(sheet: object) -> Optional[str]:
    """Stub: return the sheet name if available."""
    return None


def read_mat_file(mat_path: Path) -> object:
    """Stub: read MATLAB file and return raw loaded content."""
    raise NotImplementedError


def mat_struct_to_python(mat_struct: object) -> object:
    """Stub: convert MATLAB struct types into nested Python dict/list."""
    raise NotImplementedError


def table_has_column(table: object, colname: str) -> bool:
    """
    Stub: return True if 'table' contains a column named colname.

    This abstracts away whether 'table' is a pandas DataFrame, dict-like, etc.
    """
    raise NotImplementedError
