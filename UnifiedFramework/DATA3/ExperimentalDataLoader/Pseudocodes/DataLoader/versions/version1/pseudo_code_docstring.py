#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 11:28:46 2025

@authors: Keshav Kasturi Rangan, Alex Dowling

PSEUDOCODE — file-type agnostic experiment ingestion

Goal:
    Keep classes light (pure containers).
    Do heavy lifting in functions that parse .xlsx or .mat.
    Return a single ExperimentData object per experiment run.
"""

# ============================
# Light containers (no parsing)
# ============================

@dataclass
class VialData:
    """Store one vial segment.

    Arguments:
        number: int
        start_idx: int
        end_idx: int
        time: ndarray[float]
        mass: ndarray[float]
        retentate_signal: ndarray[float] | None
        permeate_signal: ndarray[float] | None
        extras: dict[str, MetaValue]
    """
    ...


@dataclass
class ExperimentData:
    """Store one experiment run (one sheet or one MAT dataset).

    Arguments:
        path: str
        source: 'xlsx'|'mat'
        run_id: str
        filename: str
        metadata: dict[str, MetaValue]
        config: dict[str, MetaValue]
        time: ndarray[float]
        mass: ndarray[float]
        retentate_signal: ndarray[float] | None
        permeate_signal: ndarray[float] | None
        vial_marker: ndarray[int] | None
        vial_ranges: list[(start,end)]
        vials: list[VialData]
        n_vials: int|None
    """
    ...


# ============================
# Excel functional ingestion
# ============================

def read_excel_experiment(path: str, sheet_name: str) -> ExperimentData:
    """Parse one Excel sheet (one experiment) into ExperimentData.

    Arguments:
        path: str
        sheet_name: str

    Returns:
        ExperimentData
    """
    # 1) Read raw sheet grid (header=None) so notes/metadata are preserved
    # 2) Detect where the measurement table header begins
    # 3) Parse top block into (metadata, config)
    # 4) Re-read using header_row to get clean measurement table
    # 5) Detect columns (time, mass, optional ret/perm, optional vial marker)
    # 6) Convert columns to numeric arrays
    # 7) Compute vial_ranges from vial marker (or default single vial)
    # 8) Slice arrays into VialData list
    # 9) Return ExperimentData
    ...


def read_excel_workbook(path: str) -> list[ExperimentData]:
    """Parse every sheet (each sheet is one experiment).

    Arguments:
        path: str

    Returns:
        list[ExperimentData]
    """
    # 1) list_sheets(path)
    # 2) for each sheet: read_excel_experiment(path, sheet)
    # 3) return list of ExperimentData
    ...


# ============================
# MAT functional ingestion
# ============================

def read_mat_experiment(path: str) -> ExperimentData:
    """Parse MATLAB .mat containing 'data_stru' into ExperimentData.

    Arguments:
        path: str

    Returns:
        ExperimentData
    """
    # 1) loadmat(path)
    # 2) ds = mat['data_stru']
    # 3) extract ds.dataset, ds.filename (if present)
    # 4) config = dict(ds.data_config fields)
    # 5) iterate ds.data_raw per-vial:
    #       - read time, mass
    #       - read cF_exp/cV_avg as generic signals
    #       - compute (start,end) indices
    #       - create VialData
    # 6) concatenate full-run arrays
    # 7) build vial_marker from vial_ranges
    # 8) return ExperimentData
    ...
