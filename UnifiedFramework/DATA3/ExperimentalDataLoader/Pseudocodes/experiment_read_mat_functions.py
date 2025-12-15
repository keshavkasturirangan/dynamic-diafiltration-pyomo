#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 11:27:04 2025

@authors: Keshav Kasturi Rangan, Alex Dowling

experiment_read_mat_functions.py

One-line purpose:
    Read diafiltration experiments from MATLAB .mat files containing a 'data_stru' struct.

This module:
    - Loads .mat via scipy.io.loadmat
    - Extracts dataset/filename, config struct, and per-vial data_raw records
    - Concatenates per-vial arrays to full-run arrays
    - Builds vial_ranges + full-length vial_marker
    - Returns file-type agnostic ExperimentData

No `Any` is used; extra payloads are stored as MetaValue.
"""

from __future__ import annotations  # Enable forward references.

from typing import Dict, List, Tuple  # Explicit typing (no Any).

import os  # Path utilities.

import numpy as np  # Arrays.
import scipy.io as spio  # MATLAB file reading.

from experiment_data_types import ExperimentData, VialData, MetaValue, ArrayF, ArrayI  # Shared containers/types.


def _as_array(x: object) -> ArrayF:
    """
    Convert an object to a float numpy array (best effort).

    Arguments:
        x: object, input value (MAT scalar/list/ndarray)

    Returns:
        np.ndarray, float array (squeezed)
    """
    return np.asarray(x, dtype=float).squeeze()


def _struct_to_dict(s: object) -> Dict[str, MetaValue]:
    """
    Convert a MATLAB struct-like object to a typed dict (best effort).

    Arguments:
        s: object, struct-like (has _fieldnames)

    Returns:
        dict[str, MetaValue], shallow mapping (no Any)
    """
    d: Dict[str, MetaValue] = {}
    for f in getattr(s, "_fieldnames", []):
        val = getattr(s, f)
        if isinstance(val, np.ndarray):
            d[f] = val
        elif isinstance(val, (int, float, bool, str)) or val is None:
            d[f] = val
        else:
            d[f] = str(val)
    return d


def build_vial_marker(vial_ranges: List[Tuple[int, int]], n_total: int) -> ArrayI:
    """
    Build full-length vial labels 1..N from index ranges.

    Arguments:
        vial_ranges: list[(start,end)] inclusive
        n_total: int, total length

    Returns:
        np.ndarray[int], vial label for each row
    """
    marker = np.zeros(n_total, dtype=int)
    for k, (s, e) in enumerate(vial_ranges, start=1):
        marker[s : e + 1] = k
    return marker


def read_mat_experiment(path: str) -> ExperimentData:
    """
    Read a MATLAB .mat experiment file (with 'data_stru') and return ExperimentData.

    Arguments:
        path: str, path to .mat file containing 'data_stru'

    Returns:
        ExperimentData, unified container
    """
    mat = spio.loadmat(path, squeeze_me=True, struct_as_record=False)

    if "data_stru" not in mat:
        raise ValueError(f"MAT file missing 'data_stru'. Keys={list(mat.keys())}")

    ds = mat["data_stru"]

    run_id = str(getattr(ds, "dataset", os.path.basename(path)))
    filename = str(getattr(ds, "filename", os.path.basename(path)))

    config_struct = getattr(ds, "data_config", None)
    config: Dict[str, MetaValue] = _struct_to_dict(config_struct) if config_struct is not None else {}

    metadata: Dict[str, MetaValue] = {}
    if hasattr(ds, "mode"):
        metadata["mode"] = str(getattr(ds, "mode"))
    if hasattr(ds, "continuous_cF"):
        cc = getattr(ds, "continuous_cF")
        metadata["continuous_cF"] = bool(cc) if isinstance(cc, (int, bool)) else str(cc)

    data_raw = getattr(ds, "data_raw", None)
    if data_raw is None:
        raise ValueError("data_stru missing 'data_raw' per-vial data.")

    vstructs = data_raw if isinstance(data_raw, (list, tuple, np.ndarray)) else [data_raw]

    vials: List[VialData] = []
    vial_ranges: List[Tuple[int, int]] = []

    time_parts: List[ArrayF] = []
    mass_parts: List[ArrayF] = []
    ret_parts: List[ArrayF] = []
    perm_parts: List[ArrayF] = []

    cursor = 0

    for k, v in enumerate(vstructs, start=1):
        t = _as_array(getattr(v, "time", np.array([], dtype=float)))
        m = _as_array(getattr(v, "mass", np.array([], dtype=float)))

        # MATLAB naming commonly: cF_exp (retentate) and cV_avg (permeate)
        ret_obj = getattr(v, "cF_exp", None)
        perm_obj = getattr(v, "cV_avg", None)

        ret_arr = _as_array(ret_obj) if ret_obj is not None else None
        perm_arr = _as_array(perm_obj) if perm_obj is not None else None

        start = cursor
        end = cursor + int(t.size) - 1
        vial_ranges.append((start, end))
        cursor = end + 1

        vials.append(
            VialData(
                number=k,
                start_idx=start,
                end_idx=end,
                time=t,
                mass=m,
                retentate_signal=ret_arr,
                permeate_signal=perm_arr,
                extras={},
            )
        )

        time_parts.append(t)
        mass_parts.append(m)
        if ret_arr is not None:
            ret_parts.append(ret_arr)
        if perm_arr is not None:
            perm_parts.append(perm_arr)

    time_full = np.concatenate(time_parts) if time_parts else np.array([], dtype=float)
    mass_full = np.concatenate(mass_parts) if mass_parts else np.array([], dtype=float)

    # Only concatenate signals if present for every vial (strict).
    ret_full = np.concatenate(ret_parts) if len(ret_parts) == len(time_parts) else None
    perm_full = np.concatenate(perm_parts) if len(perm_parts) == len(time_parts) else None

    vial_marker = build_vial_marker(vial_ranges, n_total=int(time_full.size))

    return ExperimentData(
        path=path,
        source="mat",
        run_id=run_id,
        filename=filename,
        metadata=metadata,
        config=config,
        time=time_full,
        mass=mass_full,
        retentate_signal=ret_full,
        permeate_signal=perm_full,
        vial_marker=vial_marker,
        vial_ranges=vial_ranges,
        vials=vials,
        n_vials=len(vials),
    )
