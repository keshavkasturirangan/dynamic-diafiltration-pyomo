#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 11:20:33 2025

@authors: Keshav Kasturi Rangan, Alex Dowling

experiment_data_types.py

One-line purpose:
    Define lightweight, file-type agnostic containers for diafiltration experimental data.

Design:
    - ExperimentData stores ONE experiment run (regardless of .mat or .xlsx source).
    - VialData stores ONE vial segment (contiguous window) within an experiment.
    - No parsing logic lives here; parsing is done by standalone functions in reader modules.

No `Any` is used; instead we use explicit unions (MetaValue) for structured payloads.
"""

from __future__ import annotations  # Enable forward references in type hints.

from dataclasses import dataclass, field  # Lightweight containers.
from typing import Dict, List, Optional, Tuple, Union, Literal  # Explicit typing (no Any).

import numpy as np  # Arrays for signals/markers.
import pandas as pd  # Optional: store raw measurement tables for provenance (still no Any).


# ----------------------------
# Explicit aliases (no Any)
# ----------------------------

ArrayF = np.ndarray  # Numeric arrays (typically float) for signals like time/mass.
ArrayI = np.ndarray  # Integer arrays for markers/labels like vial IDs.
Scalar = Union[int, float, bool, str]  # Basic scalar types.


# JSON-ish recursive type for metadata/config payloads WITHOUT Any.
MetaValue = Union[
    None,
    Scalar,
    ArrayF,
    ArrayI,
    pd.DataFrame,
    List[Scalar],
    List[float],
    Dict[str, "MetaValue"],
]


@dataclass
class VialData:
    """
    Stores ONE vial segment (one contiguous window) of an experiment.

    Arguments:
        number: int
            Vial index (recommend 1-based).
        start_idx: int
            Start index in concatenated arrays (inclusive).
        end_idx: int
            End index in concatenated arrays (inclusive).
        time: np.ndarray
            Time samples [s].
        mass: np.ndarray
            Mass samples [g].
        retentate_signal: np.ndarray | None
            Retentate-side signal (conductivity OR concentration).
        permeate_signal: np.ndarray | None
            Permeate-side signal (conductivity OR concentration).
        extras: dict[str, MetaValue]
            Extra per-vial fields (typed; no Any).
    """

    number: int
    start_idx: int
    end_idx: int

    time: ArrayF
    mass: ArrayF

    retentate_signal: Optional[ArrayF] = None
    permeate_signal: Optional[ArrayF] = None

    extras: Dict[str, MetaValue] = field(default_factory=dict)


@dataclass
class ExperimentData:
    """
    Stores ONE experiment run, regardless of source (.mat or .xlsx).

    Arguments:
        path: str
            Full path to the source file.
        source: Literal['mat','xlsx']
            Input format identifier.
        run_id: str
            Run identifier (MAT dataset id OR Excel sheet name).
        filename: str
            Display filename (usually basename(path)).
        metadata: dict[str, MetaValue]
            Descriptive information (notes, header fields, etc.).
        config: dict[str, MetaValue]
            Configuration/settings/constants (pressure, temp, membrane area, etc.).
        time: np.ndarray
            Full concatenated time vector [s].
        mass: np.ndarray
            Full concatenated mass vector [g].
        retentate_signal: np.ndarray | None
            Full retentate signal series (if available).
        permeate_signal: np.ndarray | None
            Full permeate signal series (if available).
        vial_marker: np.ndarray | None
            Integer label per row indicating which vial the row belongs to.
        vial_ranges: list[tuple[int,int]]
            Inclusive index windows (start,end) for each vial.
        vials: list[VialData]
            Per-vial segmented records.
        n_vials: int | None
            Convenience count of vials.
    """

    path: str
    source: Literal["mat", "xlsx"]
    run_id: str
    filename: str

    metadata: Dict[str, MetaValue] = field(default_factory=dict)
    config: Dict[str, MetaValue] = field(default_factory=dict)

    time: ArrayF = field(default_factory=lambda: np.array([], dtype=float))
    mass: ArrayF = field(default_factory=lambda: np.array([], dtype=float))
    retentate_signal: Optional[ArrayF] = None
    permeate_signal: Optional[ArrayF] = None
    vial_marker: Optional[ArrayI] = None

    vial_ranges: List[Tuple[int, int]] = field(default_factory=list)
    vials: List[VialData] = field(default_factory=list)
    n_vials: Optional[int] = None
