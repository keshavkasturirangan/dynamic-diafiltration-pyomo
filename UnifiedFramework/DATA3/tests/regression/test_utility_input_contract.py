"""Regression checks for the legacy utility.py input contract."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat


REPO_ROOT = Path(__file__).resolve().parents[4]


def _load_data_stru(path: Path):
    return loadmat(path, squeeze_me=True, struct_as_record=False)["data_stru"]


def _config_keys(data_stru) -> set[str]:
    return {name for name in dir(data_stru.data_config) if not name.startswith("_")}


def _raw_keys(data_stru) -> set[str]:
    return {name for name in dir(data_stru.data_raw[0]) if not name.startswith("_")}


@pytest.mark.regression
def test_data1_501_1_utility_contract_shape() -> None:
    """DATA1 501.1 uses the smaller DATA1-style legacy contract."""
    data_stru = _load_data_stru(REPO_ROOT / "DATA1_matlab" / "data_library" / "data_stru-dataset501.1.mat")

    config = _config_keys(data_stru)
    raw = _raw_keys(data_stru)

    assert {
        "Am",
        "B0",
        "C_F0",
        "Lp0",
        "M_F0",
        "Temp",
        "delP",
        "n",
        "namec",
        "ni",
        "rho",
        "sigma0",
    }.issubset(config)
    assert {"time", "mass", "cV_avg", "cF_exp"}.issubset(raw)

    # This representative DATA1 file does not carry the fuller DATA2-style fields.
    assert "n_v0" not in config
    assert "n_extra" not in config
    assert "n_h" not in config
    assert "n_A" not in config

    # DATA1 501.1 also demonstrates the scalar retentate-assay shape.
    assert isinstance(data_stru.data_raw[0].cF_exp, float)
    assert isinstance(data_stru.data_raw[0].cV_avg, float)
    assert isinstance(data_stru.data_raw[0].mass, np.ndarray)
    assert isinstance(data_stru.data_raw[0].time, np.ndarray)


@pytest.mark.regression
def test_data1_511_12_utility_contract_shape() -> None:
    """DATA1 511.12 is still DATA1-style but includes extra chemistry fields."""
    data_stru = _load_data_stru(REPO_ROOT / "DATA1_matlab" / "data_library" / "data_stru-dataset511.12.mat")

    config = _config_keys(data_stru)
    raw = _raw_keys(data_stru)

    assert {"C_D", "M_O"}.issubset(config)
    assert {"time", "mass", "cV_avg", "cF_exp"}.issubset(raw)

    # Even with C_D/M_O present, this file still lacks the staged gating fields.
    assert "n_v0" not in config
    assert "n_extra" not in config
    assert "n_h" not in config
    assert "n_A" not in config

    assert isinstance(data_stru.data_raw[0].cF_exp, np.ndarray)
    assert isinstance(data_stru.data_raw[0].cV_avg, float)


@pytest.mark.regression
def test_data2_270511_121_utility_contract_shape() -> None:
    """Representative DATA2 input carries the fuller staged contract."""
    data_stru = _load_data_stru(REPO_ROOT / "data_library" / "data_stru-dataset270511.121.mat")

    config = _config_keys(data_stru)
    raw = _raw_keys(data_stru)

    assert {
        "Am",
        "B0",
        "C_D",
        "C_F0",
        "Lp0",
        "M_F0",
        "M_O",
        "Temp",
        "delP",
        "n",
        "n_extra",
        "n_v0",
        "namec",
        "ni",
        "rho",
        "sigma0",
    }.issubset(config)
    assert {"time", "mass", "cV_avg", "cF_exp"}.issubset(raw)

    # DATA2-style files use time-series concentration measurements.
    assert isinstance(data_stru.data_raw[0].cF_exp, np.ndarray)
    assert isinstance(data_stru.data_raw[0].cV_avg, np.ndarray)
    assert isinstance(data_stru.data_raw[0].mass, np.ndarray)
    assert isinstance(data_stru.data_raw[0].time, np.ndarray)
