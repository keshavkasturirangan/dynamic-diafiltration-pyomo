#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# authors: Keshav Kasturi Rangan, Alex Dowling
"""
OOP-structured unified data loader for Excel (.xlsx/.xls/.xlsm/.xlsb) and MATLAB (.mat)
used in diafiltration experiments, with unified plotting that:
- Mass vs time: single solid color for all vials + legend.
- Conductivity vs time: retentate solid (one color), permeate dashed (one color) + legend.
- Top x-axis vial labels and vertical lines at vial boundaries.

Provides:
- DataLoader(path, *, measure_sheet=None, index_sheet=None)
    .load()
    .times, .masses, .retentate_cond, .permeate_cond, .windows
    .plot_quicklook(save_prefix=None)
- list_sheets(path) -> List[str]
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple, Optional, Protocol, runtime_checkable
import os
import re  # <-- added for metadata/config parsing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io as sio


# =============================================================================
# Unified data_stru-style ontology (PyDataStru / PyVial / PyDataConfig)
# =============================================================================

@dataclass
class PyDataConfig:
    """Python mirror of MATLAB data_stru.data_config."""
    M_F0: Optional[float] = None
    C_F0: Optional[np.ndarray] = None
    M_O: Optional[float] = None
    C_D: Optional[np.ndarray] = None
    nc: Optional[int] = None
    namec: Any = None
    ni: Optional[int] = None
    delP: Optional[float] = None
    Temp: Optional[float] = None
    Am: Optional[float] = None
    rho: Optional[float] = None
    Lp0: Optional[float] = None
    B0: Optional[np.ndarray] = None
    sigma0: Optional[float] = None
    theta0: Optional[np.ndarray] = None
    n: Optional[int] = None
    nr: Optional[int] = None
    n_extra: Optional[int] = None
    extras: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PyVial:
    """Python mirror of MATLAB data_stru.data_raw(i)."""
    number: int
    time: np.ndarray
    mass: np.ndarray
    cV_avg: Optional[np.ndarray] = None      # permeate signal
    cF_exp: Optional[np.ndarray] = None      # retentate signal
    nr: Optional[int] = None
    extras: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PyDataStru:
    """Unified experiment structure, mirroring MATLAB data_stru."""
    dataset: Optional[float]
    filename: Optional[str]
    mode: Optional[str]
    conductivity_cF: Optional[bool]
    n_obj: Optional[int]
    data_config: PyDataConfig
    data_raw: List[PyVial]

    # -------------------------------------------------------------------------
    # Constructor from MATLAB data_stru .mat (rich mapping)
    # -------------------------------------------------------------------------
    @staticmethod
    def from_mat(path: str) -> "PyDataStru":
        mat = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
        if "data_stru" not in mat:
            raise ValueError(f"{path} does not contain 'data_stru'")

        ds = mat["data_stru"]

        dataset = float(getattr(ds, "dataset", np.nan))
        filename = str(getattr(ds, "filename", None))
        mode = str(getattr(ds, "mode", None))
        conductivity_cF = bool(getattr(ds, "conductivity_cF", False))
        n_obj = int(getattr(ds, "n_obj", 0))

        dc = getattr(ds, "data_config", None)
        cfg_kwargs: Dict[str, Any] = {}
        extras: Dict[str, Any] = {}

        if dc is not None and hasattr(dc, "_fieldnames"):
            for fn in dc._fieldnames:
                val = getattr(dc, fn)
                if isinstance(val, (list, np.ndarray)):
                    val = np.asarray(val).squeeze()
                if fn in PyDataConfig.__dataclass_fields__:
                    cfg_kwargs[fn] = val
                else:
                    extras[fn] = val

        data_config = PyDataConfig(**cfg_kwargs, extras=extras)

        dr = getattr(ds, "data_raw", [])
        if isinstance(dr, np.ndarray):
            dr = dr.ravel()

        vials: List[PyVial] = []
        for i, v in enumerate(dr):
            t = np.asarray(getattr(v, "time", []), dtype=float).squeeze()
            m = np.asarray(getattr(v, "mass", []), dtype=float).squeeze()
            cV = getattr(v, "cV_avg", None)
            if cV is not None:
                cV = np.asarray(cV).squeeze()
            cF = getattr(v, "cF_exp", None)
            if cF is not None:
                cF = np.asarray(cF).squeeze()
            nr = getattr(v, "nr", None)
            num = getattr(v, "number", i + 1)

            vial_extras: Dict[str, Any] = {}
            if hasattr(v, "_fieldnames"):
                for fn in v._fieldnames:
                    if fn not in ["time", "mass", "cV_avg", "cF_exp", "nr", "number"]:
                        vial_extras[fn] = getattr(v, fn)

            vials.append(PyVial(num, t, m, cV, cF, nr, vial_extras))

        return PyDataStru(dataset, filename, mode, conductivity_cF, n_obj, data_config, vials)

    # -------------------------------------------------------------------------
    # Constructor from Excel (.xlsx) using DataLoader + rich config
    # -------------------------------------------------------------------------
    @staticmethod
    def from_excel(path: str,
                   measure_sheet: Optional[str] = None,
                   index_sheet: Optional[str] = None) -> "PyDataStru":
        """
        Build a PyDataStru from an Excel file using v4-style DataLoader logic
        for time/mass/cond windows, and a richer extraction of experiment
        configuration into PyDataConfig (including a config sheet + header metadata).
        """
        dl = DataLoader(path, measure_sheet=measure_sheet, index_sheet=index_sheet)
        dl.load()

        vials: List[PyVial] = []
        for i, (t, m, r, p) in enumerate(
            zip(dl.times, dl.masses, dl.retentate_cond, dl.permeate_cond),
            start=1
        ):
            t_arr = np.asarray(t, dtype=float).squeeze()
            m_arr = np.asarray(m, dtype=float).squeeze()
            r_arr = np.asarray(r, dtype=float).squeeze() if r is not None else None
            p_arr = np.asarray(p, dtype=float).squeeze() if p is not None else None

            vials.append(
                PyVial(
                    number=i,
                    time=t_arr,
                    mass=m_arr,
                    cF_exp=r_arr,
                    cV_avg=p_arr,
                    nr=int(np.isfinite(m_arr).sum()),
                    extras={}
                )
            )

        # richer data_config from Excel:
        data_config = _build_data_config_from_excel(path, measure_sheet, index_sheet)

        return PyDataStru(
            dataset=None,
            filename=os.path.basename(path),
            mode=None,
            conductivity_cF=False,
            n_obj=len(vials),
            data_config=data_config,
            data_raw=vials,
        )


# =============================================================================
# Excel metadata / config extraction helpers
# =============================================================================

def _detect_header_row_generic(raw: pd.DataFrame, max_rows: int = 25) -> Optional[int]:
    """
    Scan top rows of a sheet to find the header row containing both 'time'
    and ('mass' or 'cond').
    """
    for i in range(min(max_rows, len(raw))):
        row = raw.iloc[i].astype(str).str.lower().tolist()
        has_time = any("time" in s for s in row)
        has_mass_or_cond = any(("mass" in s) or ("cond" in s) for s in row)
        if has_time and has_mass_or_cond:
            return i
    return None


def _extract_excel_metadata(path: str, measure_sheet: Optional[str]) -> Dict[str, Any]:
    """
    Extract human-readable metadata from the rows ABOVE the header row
    in an Excel measurement sheet.

    Returns keys:
        - excel_metadata_rows: list[str]
        - feed_description: str (optional)
        - diafiltrate_description: str (optional)
        - notes: list[str] (optional)
        - stirring_rpm: int (optional)
        - datapoints: int (optional)
    """
    xl = pd.ExcelFile(path)
    sheet = measure_sheet or xl.sheet_names[0]
    raw = xl.parse(sheet, header=None)

    hdr = _detect_header_row_generic(raw, max_rows=25)
    if hdr is None:
        meta_block = raw.iloc[:10]
    else:
        meta_block = raw.iloc[:hdr]

    metadata_rows: List[str] = []
    extras: Dict[str, Any] = {}

    for _, row in meta_block.iterrows():
        # Drop NaNs and 'Unnamed: X'-style junk
        tokens = []
        for v in row:
            if pd.isna(v):
                continue
            s = str(v).strip()
            if not s:
                continue
            if s.lower().startswith("unnamed:"):
                continue
            tokens.append(s)
        if tokens:
            text = " ".join(tokens)
            metadata_rows.append(text)

    extras["excel_metadata_rows"] = metadata_rows

    notes: List[str] = []
    for txt in metadata_rows:
        low = txt.lower()

        # Feed description heuristics
        if "feed" in low:
            extras.setdefault("feed_description", txt)

        # Diafiltrate / dialysate
        if "diafiltrate" in low or "dialysate" in low:
            extras.setdefault("diafiltrate_description", txt)

        # Notes
        if "note" in low:
            notes.append(txt)

        # Stirring RPM, e.g. "350 RPM" or "350RPM"
        m = re.search(r"(\d+)\s*rpm", low)
        if m:
            try:
                extras["stirring_rpm"] = int(m.group(1))
            except ValueError:
                pass

        # Datapoints: 404
        if "datapoints" in low:
            m2 = re.search(r"datapoints\s*[:\-]\s*(\d+)", low)
            if m2:
                try:
                    extras["datapoints"] = int(m2.group(1))
                except ValueError:
                    pass

    if notes:
        extras["notes"] = notes

    return extras


def _extract_excel_config_pairs(path: str,
                                measure_sheet: Optional[str],
                                index_sheet: Optional[str]) -> Dict[str, Any]:
    """
    Look for a 'config-like' sheet (Config, Setup, Parameters, etc.) and
    parse it as key–value pairs. The pattern assumed is:

        Row: [ key, value, ... ]

    where the *first* non-empty cell in the row is interpreted as the key
    and the *second* non-empty cell as the value.

    Returns a dict mapping keys (str) -> values (float, str, or list).
    """
    xl = pd.ExcelFile(path)
    # Candidate sheets: anything that's not the measurement or index sheet
    candidates = []
    for sh in xl.sheet_names:
        if sh == measure_sheet or sh == index_sheet:
            continue
        low = sh.lower()
        if any(k in low for k in ["config", "setup", "parameter", "params"]):
            candidates.append(sh)

    if not candidates:
        return {}

    config_pairs: Dict[str, Any] = {}
    for sh in candidates:
        df = xl.parse(sh, header=None)
        for _, row in df.iterrows():
            # Collect non-empty tokens
            tokens = [str(v).strip() for v in row if pd.notna(v) and str(v).strip()]
            if len(tokens) < 2:
                continue
            key = tokens[0]
            val_raw = tokens[1]
            # Try to numeric-cast, else keep as string
            try:
                val = float(val_raw)
            except ValueError:
                val = val_raw
            # If key already exists, don't overwrite; or store as list if needed
            k_norm = key.strip()
            if k_norm in config_pairs:
                old = config_pairs[k_norm]
                if isinstance(old, list):
                    old.append(val)
                else:
                    config_pairs[k_norm] = [old, val]
            else:
                config_pairs[k_norm] = val

    return config_pairs


def _build_data_config_from_excel(path: str,
                                  measure_sheet: Optional[str],
                                  index_sheet: Optional[str]) -> PyDataConfig:
    """
    Construct a PyDataConfig from:
      - free-text metadata rows (above header),
      - config-sheet key–value pairs.

    Strategy:
      1) Extract metadata rows -> meta_extras dict.
      2) Extract config sheet pairs -> cfg_pairs dict.
      3) For keys in cfg_pairs that match PyDataConfig field names,
         assign them directly (with numeric/array casts where appropriate).
      4) Put everything else into extras, alongside meta_extras.
    """
    meta_extras = _extract_excel_metadata(path, measure_sheet)
    cfg_pairs = _extract_excel_config_pairs(path, measure_sheet, index_sheet)

    # Known numeric fields of PyDataConfig we can safely cast
    numeric_fields = {
        "M_F0", "M_O", "nc", "ni", "delP", "Temp", "Am", "rho",
        "Lp0", "sigma0", "n", "nr", "n_extra"
    }

    # Fields that are arrays / vectors
    array_fields = {"C_F0", "C_D", "B0", "theta0"}

    cfg_kwargs: Dict[str, Any] = {}
    extras: Dict[str, Any] = {}

    # 3) Map config-sheet pairs to PyDataConfig fields where possible
    for key, val in cfg_pairs.items():
        key_norm = key.strip()
        if key_norm in PyDataConfig.__dataclass_fields__:
            if key_norm in numeric_fields:
                try:
                    cfg_kwargs[key_norm] = float(val)
                except Exception:
                    cfg_kwargs[key_norm] = val
            elif key_norm in array_fields:
                if isinstance(val, str):
                    parts = re.split(r"[,\s]+", val.strip())
                    parts = [p for p in parts if p]
                    try:
                        arr = np.array([float(p) for p in parts], dtype=float)
                        cfg_kwargs[key_norm] = arr
                    except Exception:
                        cfg_kwargs[key_norm] = val
                else:
                    cfg_kwargs[key_norm] = np.atleast_1d(val)
            else:
                cfg_kwargs[key_norm] = val
        else:
            extras[key_norm] = val

    # 4) Merge meta_extras (from header rows) into extras
    for k, v in meta_extras.items():
        if k in extras:
            old = extras[k]
            if isinstance(old, list):
                if isinstance(v, list):
                    old.extend(v)
                else:
                    old.append(v)
            else:
                extras[k] = [old, v]
        else:
            extras[k] = v

    return PyDataConfig(**cfg_kwargs, extras=extras)


# -------------------------------- Data container -----------------------------------------------

@dataclass
class DataRun:
    time: np.ndarray
    mass: np.ndarray
    retentate_cond: np.ndarray
    permeate_cond: np.ndarray
    vial_marker: Optional[np.ndarray] = None
    extra_signals: Dict[str, np.ndarray] = None
    metadata: Dict[str, Any] = None

    def __post_init__(self):
        n = len(np.asarray(self.time).squeeze())
        if self.extra_signals is None:
            self.extra_signals = {}
        if self.metadata is None:
            self.metadata = {}
        def _fit(x, fill=np.nan):
            if x is None:
                arr = np.full(n, fill)
            else:
                arr = np.asarray(x).squeeze()
            arr = arr[:n] if arr.size >= n else np.pad(arr, (0, n - arr.size), constant_values=fill)
            return arr
        self.time = _fit(self.time).astype(float)
        self.mass = _fit(self.mass).astype(float)
        self.retentate_cond = _fit(self.retentate_cond).astype(float)
        self.permeate_cond = _fit(self.permeate_cond).astype(float)
        if self.vial_marker is not None:
            vm = np.asarray(self.vial_marker).squeeze()
            self.vial_marker = vm[:n] if vm.size >= n else np.pad(vm, (0, n - vm.size), constant_values=np.nan)

# -------------------------------- Data sources --------------------------------------------------

@runtime_checkable
class DataSource(Protocol):
    def sniff(self) -> bool: ...
    def parse(self) -> DataRun: ...

class ExcelSource:
    """Parses Excel workbooks with offset headers and optional index sheet."""
    def __init__(self, path: str, measure_sheet_override: Optional[str] = None, index_sheet_override: Optional[str] = None):
        self.path = path
        self._measure_sheet_override = measure_sheet_override
        self._index_sheet_override = index_sheet_override
        self.measure_sheet: Optional[str] = None
        self.index_sheet: Optional[str] = None
        self._meas_df: Optional[pd.DataFrame] = None

    def sniff(self) -> bool:
        return os.path.splitext(self.path)[1].lower() in {".xlsx", ".xls", ".xlsm", ".xlsb"}

    # ---------- sheet detection ----------
    def _detect_measurement_sheet(self, xl: pd.ExcelFile) -> Optional[str]:
        for sh in xl.sheet_names:
            df = xl.parse(sh, header=None, nrows=10)
            if self._detect_header_row(df, max_rows=10) is not None:
                return sh
        return None

    def _detect_index_sheet(self, xl: pd.ExcelFile) -> Optional[str]:
        for sh in xl.sheet_names:
            if "index" in sh.lower() or "summary" in sh.lower():
                return sh
        return None

    # ---------- header / column detection ----------
    def _detect_header_row(self, df: pd.DataFrame, max_rows: int = 10) -> Optional[int]:
        # reuse the same logic as _detect_header_row_generic
        for i in range(min(max_rows, len(df))):
            row = df.iloc[i].astype(str).str.lower().tolist()
            has_time = any("time" in s for s in row)
            has_mass_or_cond = any(("mass" in s) or ("cond" in s) for s in row)
            if has_time and has_mass_or_cond:
                return i
        return None

    def _first_col_like(self, df: pd.DataFrame, needle: str) -> Optional[str]:
        n = needle.lower()
        for c in df.columns:
            if n in str(c).lower():
                return c
        return None

    def _first_col_from_list(self, df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
        lowers = {c.lower(): c for c in df.columns}
        for cand in candidates:
            if cand.lower() in lowers:
                return lowers[cand.lower()]
        for cand in candidates:
            for c in df.columns:
                if cand.lower() in str(c).lower():
                    return c
        return None

    # ---------- parse ----------
    def parse(self) -> DataRun:
        xl = pd.ExcelFile(self.path)
        # Resolve measurement sheet
        if self._measure_sheet_override is not None:
            if self._measure_sheet_override not in xl.sheet_names:
                raise ValueError(f"Measurement sheet '{self._measure_sheet_override}' not found.")
            self.measure_sheet = self._measure_sheet_override
        else:
            self.measure_sheet = self._detect_measurement_sheet(xl)
        if self.measure_sheet is None:
            raise ValueError("No measurement-like sheet (Time + Mass/Cond) found.")
        # Resolve index sheet
        if self._index_sheet_override is not None:
            if self._index_sheet_override not in xl.sheet_names:
                raise ValueError(f"Index sheet '{self._index_sheet_override}' not found.")
            self.index_sheet = self._index_sheet_override
        else:
            self.index_sheet = self._detect_index_sheet(xl)

        raw = xl.parse(self.measure_sheet, header=None)
        hdr = self._detect_header_row(raw, max_rows=10)
        df = raw.iloc[hdr + 1:].copy()
        df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
        df = df.reset_index(drop=True)
        self._meas_df = df

        tcol = self._first_col_like(df, "time")
        if tcol is None:
            raise ValueError("Time column not found in measurement sheet.")
        time = pd.to_numeric(df[tcol], errors="coerce").to_numpy()

        mcol = self._first_col_from_list(df, ["Mass", "mass", "M (g)", "Permeate Mass", "Collection mass", "Weight"])
        if mcol is None:
            raise ValueError("Mass column not found in measurement sheet.")
        mass = pd.to_numeric(df[mcol], errors="coerce").to_numpy()

                # Conductivity columns (retentate/permeate) if present
        rcol = self._first_col_from_list(df,
            ["Retentate Cond", "Retentate_Cond", "RetentateCond", "Ret_Cond"],)
        pcol = self._first_col_from_list(df,
            ["Permeate Cond", "Permeate_Cond", "PermeateCond", "Perm_Cond"],)
        rcon = (pd.to_numeric(df[rcol], errors="coerce").to_numpy()
            if rcol is not None
            else None)
        pcon = (pd.to_numeric(df[pcol], errors="coerce").to_numpy()
            if pcol is not None
            else None)

        # Optional vial marker column used by one of the window strategies
        # Use substring matching so columns named "Vial swap", "VialSwap", etc. are picked up.
        vcol = self._first_col_like(df, "vial") or self._first_col_like(df, "vial swap")
        vial = (pd.to_numeric(df[vcol], errors="coerce").to_numpy()
            if vcol is not None
            else None)

        return DataRun(time=time, mass=mass,
                       retentate_cond=rcon if rcon is not None else np.full_like(time, np.nan),
                       permeate_cond=pcon if pcon is not None else np.full_like(time, np.nan),
                       vial_marker=vial)

class MatSource:
    """
    Parses MATLAB .mat files.

    Two modes:

    1) If the file contains a MATLAB data_stru struct (the diafiltration format),
       we:
         - read data_stru.data_raw(i).time, mass, cF_exp, cV_avg
         - concatenate all vials into a single DataRun
         - build a vial_marker array so the windowing logic can
           reconstruct per-vial windows.

    2) Otherwise, we fall back to generic flat arrays:
         Time / Mass / RetentateCond / PermeateCond / Vial{,Swap}
    """
    def __init__(self, path: str):
        self.path = path

    def sniff(self) -> bool:
        return os.path.splitext(self.path)[1].lower() == ".mat"

    def _try_scipy(self):
        try:
            from scipy.io import loadmat
            return loadmat(self.path, squeeze_me=True, struct_as_record=False)
        except Exception:
            return None

    def _try_h5py(self):
        try:
            import h5py
            return h5py.File(self.path, "r")
        except Exception:
            return None

    def parse(self) -> DataRun:
        # First, try the usual MATLAB v5/v7 loader
        mat = self._try_scipy()

        # ------------------------------------------------------------------
        # CASE 1: data_stru-based diafiltration file
        # ------------------------------------------------------------------
        if mat is not None and "data_stru" in mat:
            ds = mat["data_stru"]

            # data_raw can be a struct array or a single struct
            dr = getattr(ds, "data_raw", None)
            if dr is None:
                raise ValueError("MAT file has 'data_stru' but no 'data_raw' field.")

            if isinstance(dr, np.ndarray):
                dr_iter = dr.ravel()
            else:
                dr_iter = [dr]

            times: List[np.ndarray] = []
            masses: List[np.ndarray] = []
            rcons: List[np.ndarray] = []
            pcons: List[np.ndarray] = []
            vmarker: List[np.ndarray] = []

            for i, v in enumerate(dr_iter, start=1):
                # Time and mass for this vial
                t = np.asarray(getattr(v, "time", []), dtype=float).squeeze()
                m = np.asarray(getattr(v, "mass", []), dtype=float).squeeze()

                # Optional conductivity traces
                cF = getattr(v, "cF_exp", None)
                cV = getattr(v, "cV_avg", None)

                # Align lengths to time vector
                def _align(arr, n):
                    if arr is None:
                        return np.full(n, np.nan, dtype=float)
                    a = np.asarray(arr, dtype=float).squeeze()
                    if a.size >= n:
                        return a[:n]
                    return np.pad(a, (0, n - a.size), constant_values=np.nan)

                r = _align(cF, t.size)
                p = _align(cV, t.size)

                times.append(t)
                masses.append(m)
                rcons.append(r)
                pcons.append(p)
                vmarker.append(np.full(t.shape, i, dtype=float))

            if times:
                time = np.concatenate(times)
                mass = np.concatenate(masses)
                rcon = np.concatenate(rcons)
                pcon = np.concatenate(pcons)
                vial = np.concatenate(vmarker)
            else:
                time = np.array([], dtype=float)
                mass = np.array([], dtype=float)
                rcon = np.array([], dtype=float)
                pcon = np.array([], dtype=float)
                vial = None

            return DataRun(
                time=time,
                mass=mass,
                retentate_cond=rcon,
                permeate_cond=pcon,
                vial_marker=vial,
            )

        # ------------------------------------------------------------------
        # CASE 2: fallback – generic MAT file with flat arrays
        # ------------------------------------------------------------------
        if mat is None:
            # Try v7.3 HDF5-style
            h5 = self._try_h5py()
            if h5 is None:
                raise ValueError("Could not read .mat with scipy or h5py.")

            def _h5(name: str):
                return np.array(h5[name]).squeeze() if name in h5 else None

            time = _h5("Time") or _h5("time")
            mass = _h5("Mass") or _h5("mass")
            rcon = _h5("RetentateCond") or _h5("Retentate_Cond") or _h5("retentate_cond")
            pcon = _h5("PermeateCond") or _h5("Permeate_Cond") or _h5("permeate_cond")
            vial = _h5("Vial") or _h5("VialSwap") or _h5("vial") or _h5("vial_swap")
            h5.close()
        else:
            # SciPy MAT, but no data_stru: try to find flat arrays
            def _find_key(d, *cands):
                keys = {k.lower(): k for k in d.keys()}
                for c in cands:
                    if c.lower() in keys:
                        return d[keys[c.lower()]]
                return None

            def _grab(*names):
                v = _find_key(mat, *names)
                if v is not None:
                    return np.array(v).squeeze()
                return None

            time = _grab("Time", "time", "t")
            mass = _grab("Mass", "mass", "m")
            rcon = _grab("RetentateCond", "Retentate_Cond", "retentate_cond", "ret_cond")
            pcon = _grab("PermeateCond", "Permeate_Cond", "permeate_cond", "perm_cond")
            vial = _grab("Vial", "VialSwap", "vial", "vial_swap")

        if time is None or mass is None:
            raise ValueError("MAT file missing required arrays: Time and Mass.")

        return DataRun(
            time=np.asarray(time),
            mass=np.asarray(mass),
            retentate_cond=np.asarray(rcon) if rcon is not None else None,
            permeate_cond=np.asarray(pcon) if pcon is not None else None,
            vial_marker=np.asarray(vial) if vial is not None else None,
        )


# ------------------------------- Windowing strategies -----------------------------------------

class WindowStrategy(Protocol):
    def infer(self, source: ExcelSource, meas_df: Optional[pd.DataFrame], run: DataRun) -> List[Tuple[int, int]]: ...

class IndexSheetStrategy:
    def infer(self, source: ExcelSource, meas_df: Optional[pd.DataFrame], run: DataRun) -> List[Tuple[int, int]]:
        windows: List[Tuple[int, int]] = []
        if not isinstance(source, ExcelSource) or source.index_sheet is None:
            return windows
        xl = pd.ExcelFile(source.path)
        try:
            idx_df = xl.parse(source.index_sheet)
        except Exception:
            return windows
        def first_col_like(df: pd.DataFrame, needle: str):
            n = needle.lower()
            for c in df.columns:
                if n in str(c).lower():
                    return c
            return None
        start_col = first_col_like(idx_df, "start")
        end_col = first_col_like(idx_df, "end")
        if start_col is None or end_col is None:
            return windows
        starts = pd.to_numeric(idx_df[start_col], errors="coerce").dropna().astype(int).tolist()
        ends = pd.to_numeric(idx_df[end_col], errors="coerce").dropna().astype(int).tolist()
        for s, e in zip(starts, ends):
            if 0 <= s < len(run.time) and 0 <= e < len(run.time) and s <= e:
                windows.append((s, e))
        return windows

class VialMarkerStrategy:
    def infer(
        self,
        source: ExcelSource,
        meas_df: Optional[pd.DataFrame],
        run: DataRun,
    ) -> List[Tuple[int, int]]:
        """
        Infer vial windows from a 'vial_marker' signal.

        - Tries to cast vial markers to numeric; if that fails, falls back to
          changes in the string representation (run-length encoding).
        - Splits windows whenever the vial index changes.
        - Enforces a minimum window length (>= 3 points) to avoid spurious splits.
        """
        v = run.vial_marker
        if v is None or np.asarray(v).size == 0:
            return []

        v_arr = np.asarray(v).squeeze()

        # Try to interpret as numeric vial indices
        try:
            vnum = pd.to_numeric(pd.Series(v_arr), errors="coerce").fillna(-1).to_numpy()
        except Exception:
            # Fallback: treat as strings and split when label changes
            v_str = np.array(list(map(str, v_arr)))
            # Change points → cumulative sum → "block id"
            vnum = (v_str != np.roll(v_str, 1)).astype(int).cumsum()

        # Boundaries whenever vial id changes
        boundaries = [0] + [
            i for i in range(1, len(vnum)) if vnum[i] != vnum[i - 1]
        ] + [len(vnum)]

        windows: List[Tuple[int, int]] = []
        for s, e in zip(boundaries[:-1], boundaries[1:]):
            # require at least 3 points in a window to be considered a vial
            if e - s >= 3:
                windows.append((int(s), int(e - 1)))
        return windows

class SingleWindowStrategy:
    def infer(self, source: ExcelSource, meas_df: Optional[pd.DataFrame], run: DataRun) -> List[Tuple[int, int]]:
        if len(run.time) == 0:
            return []
        return [(0, len(run.time) - 1)]

# --------------------------------- DataLoader -------------------------------------------------

class DataLoader:
    """Selects a DataSource, infers windows, exposes per-window arrays, and plots quicklooks."""
    def __init__(self, path: str, *, measure_sheet: Optional[str] = None, index_sheet: Optional[str] = None):
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")
        self.path = path
        self._measure_sheet_override = measure_sheet
        self._index_sheet_override = index_sheet
        # outputs
        self.times: List[np.ndarray] = []
        self.masses: List[np.ndarray] = []
        self.retentate_cond: List[np.ndarray] = []
        self.permeate_cond: List[np.ndarray] = []
        self.windows: List[Tuple[int, int]] = []
        # internals
        self._source: Optional[DataSource] = None
        self._run: Optional[DataRun] = None

    # ---------- source selection ----------
    def _select_source(self) -> None:
        excel = ExcelSource(self.path, self._measure_sheet_override, self._index_sheet_override)
        mats = MatSource(self.path)
        if excel.sniff():
            self._source = excel
        elif mats.sniff():
            self._source = mats
        else:
            raise ValueError(f"Unsupported file type: {self.path}")

    # ---------- window inference ----------
    def _infer_windows(self) -> None:
        assert self._run is not None
        meas_df = self._source._meas_df if isinstance(self._source, ExcelSource) else None
        strategies: List[WindowStrategy] = [
            IndexSheetStrategy(),
            VialMarkerStrategy(),
            SingleWindowStrategy()
        ]
        for strat in strategies:
            w = strat.infer(self._source, meas_df, self._run)
            if w:
                self.windows = w
                return
        self.windows = [(0, len(self._run.time) - 1)] if len(self._run.time) else []

    # ---------- public API ----------
    def load(self) -> None:
        self._select_source()
        assert self._source is not None
        self._run = self._source.parse()
        self._infer_windows()
        self.times.clear()
        self.masses.clear()
        self.retentate_cond.clear()
        self.permeate_cond.clear()
    
        for s, e in self.windows:
            sl = slice(s, e + 1)
            t = self._run.time[sl]
            mraw = self._run.mass[sl]
            # re-baseline per vial
            m = mraw - (mraw[0] if mraw.size and not np.isnan(mraw[0]) else 0.0)
            # spike filter with fixed threshold
            m = np.where(m > 1.2, np.nan, m)
            r = self._run.retentate_cond[sl]
            p = self._run.permeate_cond[sl]
            self.times.append(t)
            self.masses.append(m)
            self.retentate_cond.append(r)
            self.permeate_cond.append(p)

    # ---------- plotting (with vial labels and unified colors) ----------
    def _vial_boundaries_and_centers(self):
        boundary_times: List[float] = []
        centers: List[float] = []
        labels: List[str] = []
        for i, t in enumerate(self.times):
            if len(t):
                centers.append(float(np.nanmean(t)))
                labels.append(f"Vial {i+1}")
        for t in self.times[:-1]:
            if len(t):
                boundary_times.append(float(t[-1]))
        return boundary_times, centers, labels

    def _decorate_with_vials(self, ax, boundary_times, centers, labels):
        for bt in boundary_times:
            ax.axvline(bt, color="gray", linestyle=":", linewidth=1)
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        if centers:
            ax2.set_xticks(centers)
            ax2.set_xticklabels(labels)
        ax2.set_xlabel("Vials")

    def plot_mass(self, save_prefix: Optional[str] = None, *, cumulative: bool = False) -> None:
        """Plot Mass vs Time.
        Args:
            save_prefix: if provided, saves to f"{save_prefix}_mass.png" (per-vial) or
                         f"{save_prefix}_mass_cum.png" (cumulative)
            cumulative: if True, plot a single additive line across vials; otherwise
                        plot per-vial Δmass with a single color and legend.
        """
        boundary_times, centers, labels = self._vial_boundaries_and_centers()
        if not any(len(m) > 0 for m in self.masses):
            return
        fig, ax = plt.subplots()
        color_mass = "blue"
    
        if cumulative:
            cumulative_offset = 0.0
            first = True
            for t, m in zip(self.times, self.masses):
                if len(t) == 0 or len(m) == 0:
                    continue
                mask = ~np.isnan(m)
                if mask.any():
                    y = m[mask] + cumulative_offset
                    ax.plot(
                        t[mask],
                        y,
                        color=color_mass,
                        linewidth=1.5,
                        label="Cumulative mass change" if first else None,
                    )
                    cumulative_offset += float(m[mask][-1])
                    first = False
            ax.set_title("Cumulative Permeate Mass vs Time")
            ax.set_ylabel("Cumulative permeate mass [g]")
        else:
            first = True
            for t, m in zip(self.times, self.masses):
                if len(t) == 0 or len(m) == 0:
                    continue
                mask = ~np.isnan(m)
                if mask.any():
                    ax.plot(
                        t[mask],
                        m[mask],
                        color=color_mass,
                        linewidth=1.5,
                        label="Mass change" if first else None,
                    )
                    first = False
            ax.set_title("Mass vs Time")
            ax.set_ylabel("Mass change [g]")
    
        ax.set_xlabel("Time [s]")
        self._decorate_with_vials(ax, boundary_times, centers, labels)
        if not cumulative or not first:
            ax.legend(loc="best")
        fig.tight_layout()
        if save_prefix:
            suffix = "mass_cum" if cumulative else "mass"
            fig.savefig(f"{save_prefix}_{suffix}.png", dpi=160)
            # optional: plt.close(fig) if you don't want it in Spyder
        else:
            plt.show()

    def plot_conductivity(self, save_prefix: Optional[str] = None) -> None:
        """Plot Conductivity vs Time with unified colors (retentate solid, permeate dashed)."""
    
        boundary_times, centers, labels = self._vial_boundaries_and_centers()
    
        fig, ax = plt.subplots()
    
        color_ret = "black"
        color_perm = "red"
    
        # --- Plot retentate (solid) ---
        first_ret = True
        for t, r in zip(self.times, self.retentate_cond):
            if r is None or len(r) == 0:
                continue
            mask = np.isfinite(r)
            if mask.any():
                ax.plot(
                    t[mask],
                    r[mask],
                    color=color_ret,
                    linewidth=1.5,
                    label="Retentate (solid)" if first_ret else None,
                )
                first_ret = False
    
        # --- Plot permeate (dashed) ---
        first_perm = True
        for t, p in zip(self.times, self.permeate_cond):
            if p is None or len(p) == 0:
                continue
            mask = np.isfinite(p)
            if mask.any():
                ax.plot(
                    t[mask],
                    p[mask],
                    linestyle="--",
                    color=color_perm,
                    linewidth=1.5,
                    label="Permeate (dashed)" if first_perm else None,
                )
                first_perm = False
    
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Conductivity")
        ax.set_title("Conductivity vs Time")
    
        # Vial boundaries & labels
        self._decorate_with_vials(ax, boundary_times, centers, labels)
    
        # --- Only show legend if the signals exist ---
        if not first_ret or not first_perm:
            ax.legend(loc="best")
    
        fig.tight_layout()
    
        if save_prefix:
            fig.savefig(f"{save_prefix}_cond.png", dpi=300)
        else:
            plt.show()


    def plot_quicklook(self, save_prefix: Optional[str] = None, *, cumulative: bool = False) -> None:
        """Two-panel quicklook: Mass & Conductivity vs Time."""
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
        # Mass
        if cumulative:
            for t, m in zip(self.times, self.masses):
                ax1.step(t, m, where="post", alpha=0.7)
        else:
            for t, m in zip(self.times, self.masses):
                ax1.plot(t, m, alpha=0.7)
        ax1.set_ylabel("Mass")
        # Cond
        for t, r in zip(self.times, self.retentate_cond):
            ax2.plot(t, r, label="Retentate", alpha=0.7)
        for t, p in zip(self.times, self.permeate_cond):
            ax2.plot(t, p, "--", label="Permeate", alpha=0.7)
        ax2.set_xlabel("Time")
        ax2.set_ylabel("Conductivity")
        boundary_times, centers, labels = self._vial_boundaries_and_centers()
        self._decorate_with_vials(ax1, boundary_times, centers, labels)
        ax2.legend()
        fig.tight_layout()
        if save_prefix:
            fig.savefig(f"{save_prefix}_mass.png", dpi=300)
            plt.close(fig)          # only close when saving to file
        else:
            plt.show()              # let Spyder manage the figure; don't close


# =============================================================================
# High-level experiment wrappers (.xlsx vs .mat) + unified plotting
# =============================================================================

@dataclass
class ExcelExperimentData:
    """High–level wrapper for experimental data stored in Excel (.xlsx)."""
    path: str
    measure_sheet: Optional[str]
    index_sheet: Optional[str]
    data: PyDataStru

    @classmethod
    def from_file(cls,
                  path: str,
                  measure_sheet: Optional[str] = None,
                  index_sheet: Optional[str] = None) -> "ExcelExperimentData":
        pds = PyDataStru.from_excel(path, measure_sheet=measure_sheet, index_sheet=index_sheet)
        return cls(path=path, measure_sheet=measure_sheet, index_sheet=index_sheet, data=pds)

    @property
    def vials(self) -> List[PyVial]:
        return self.data.data_raw

    @property
    def times(self) -> List[np.ndarray]:
        return [v.time for v in self.vials]

    @property
    def masses(self) -> List[np.ndarray]:
        return [v.mass for v in self.vials]

    @property
    def retentate_cond(self) -> List[Optional[np.ndarray]]:
        return [v.cF_exp for v in self.vials]

    @property
    def permeate_cond(self) -> List[Optional[np.ndarray]]:
        return [v.cV_avg for v in self.vials]


@dataclass
class MatExperimentData:
    """High–level wrapper for experimental data stored as MATLAB data_stru (.mat)."""
    path: str
    data: PyDataStru

    @classmethod
    def from_file(cls, path: str) -> "MatExperimentData":
        pds = PyDataStru.from_mat(path)
        return cls(path=path, data=pds)

    @property
    def vials(self) -> List[PyVial]:
        return self.data.data_raw

    @property
    def times(self) -> List[np.ndarray]:
        return [v.time for v in self.vials]

    @property
    def masses(self) -> List[np.ndarray]:
        return [v.mass for v in self.vials]

    @property
    def retentate_cond(self) -> List[Optional[np.ndarray]]:
        return [v.cF_exp for v in self.vials]

    @property
    def permeate_cond(self) -> List[Optional[np.ndarray]]:
        return [v.cV_avg for v in self.vials]


def plot_mass_experiment(exp: Any, ax=None, label_prefix: str = "Vial"):
    """Plot mass vs time for any ExcelExperimentData or MatExperimentData."""
    if ax is None:
        fig, ax = plt.subplots()

    for i, (t, m) in enumerate(zip(exp.times, exp.masses), start=1):
        if t is None or m is None or len(t) == 0:
            continue
        ax.plot(t, m, label=f"{label_prefix} {i}")

    ax.set_xlabel("Time")
    ax.set_ylabel("Mass")
    ax.legend()
    return ax


def plot_conductivity_experiment(exp: Any, ax=None, label_prefix: str = "Vial"):
    """Plot retentate & permeate conductivity vs time for either experiment type."""
    if ax is None:
        fig, ax = plt.subplots()

    # Retentate (cF_exp)
    for i, (t, cF) in enumerate(zip(exp.times, exp.retentate_cond), start=1):
        if cF is None or len(cF) == 0:
            continue
        ax.plot(t[: len(cF)], cF, label=f"{label_prefix} {i} retentate")

    # Permeate (cV_avg), dashed
    for i, (t, cV) in enumerate(zip(exp.times, exp.permeate_cond), start=1):
        if cV is None or len(cV) == 0:
            continue
        ax.plot(t[: len(cV)], cV, linestyle="--", label=f"{label_prefix} {i} permeate")

    ax.set_xlabel("Time")
    ax.set_ylabel("Conductivity")
    ax.legend()
    return ax


# ------------------------------------------- CLI helpers ---------------------------------------

def list_sheets(path: str) -> List[str]:
    ext = os.path.splitext(path)[1].lower()
    if ext not in {".xlsx", ".xls", ".xlsm", ".xlsb"}:
        return []
    xl = pd.ExcelFile(path)
    return list(xl.sheet_names)


def _cli_list_sheets(path: str) -> None:
    ext = os.path.splitext(path)[1].lower()
    if ext not in {".xlsx", ".xls", ".xlsm", ".xlsb"}:
        print("[info] Not an Excel file; no sheets to list.")
        return
    sheets = list_sheets(path)
    print(f"Sheets in {path}:")
    for s in sheets:
        print("  ", s)


def _cli_dump_windows(loader: DataLoader) -> None:
    print("Windows (start_idx, end_idx):")
    for i, (s, e) in enumerate(loader.windows, start=1):
        print(f"  Vial {i}: ({s}, {e})")


def _cli_export_csv(loader: DataLoader, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    for i, (t, m, r, p) in enumerate(zip(loader.times, loader.masses, loader.retentate_cond, loader.permeate_cond), start=1):
        df = pd.DataFrame({
            "Time": t,
            "Mass": m,
            "RetentateCond": r,
            "PermeateCond": p,
        })
        df.to_csv(os.path.join(outdir, f"vial_{i}.csv"), index=False)
    print(f"[info] Exported {len(loader.times)} vials to {outdir}")


def _main():
    import argparse
    parser = argparse.ArgumentParser(description="OOP-structured unified data loader for Excel/MAT diafiltration experiments.")
    parser.add_argument("path", help="Path to .xlsx/.xls/.xlsm/.xlsb or .mat file")
    parser.add_argument("--quicklook", action="store_true", help="Show two-panel quicklook plot")
    parser.add_argument("--save-prefix", type=str, default=None, help="If provided, save figures with this prefix instead of showing.")
    parser.add_argument("--list-sheets", action="store_true", help="List sheets (for Excel file) and exit.")
    parser.add_argument("--dump-windows", action="store_true", help="Print inferred windows and exit.")
    parser.add_argument("--export-csv", type=str, default=None, help="Directory to export per-vial CSV files.")
    parser.add_argument("--measure-sheet", type=str, default=None, help="Override measurement sheet name (Excel only).")
    parser.add_argument("--index-sheet", type=str, default=None, help="Override index sheet name (Excel only).")
    parser.add_argument("--cumulative", action="store_true", help="Plot cumulative mass instead of raw mass.")
    parser.add_argument("--mass-only", action="store_true", help="Only plot mass, not conductivity.")
    parser.add_argument("--cond-only", action="store_true", help="Only plot conductivity, not mass.")
    args = parser.parse_args()

    if args.list_sheets:
        _cli_list_sheets(args.path)
        return

    loader = DataLoader(args.path, measure_sheet=args.measure_sheet, index_sheet=args.index_sheet)
    loader.load()

    if args.dump_windows:
        _cli_dump_windows(loader)

    if args.quicklook or args.save_prefix:
        if args.mass_only and not args.cond_only:
            loader.plot_mass(save_prefix=args.save_prefix, cumulative=args.cumulative)
        elif args.cond_only and not args.mass_only:
            loader.plot_conductivity(save_prefix=args.save_prefix)
        else:
            loader.plot_quicklook(save_prefix=args.save_prefix, cumulative=args.cumulative)

    if args.export_csv:
        _cli_export_csv(loader, args.export_csv)


if __name__ == '__main__':
    _main()
