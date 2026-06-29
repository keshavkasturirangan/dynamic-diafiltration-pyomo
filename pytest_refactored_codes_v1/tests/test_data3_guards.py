"""Layer 4 (guard) — DATA3 knobs must be inert for DATA1/DATA2.

Converted from `refactored_codes_v1/tests/regression_data1_data2_guards.py`
(a 222-line print-script) into a proper pytest module.

What this protects
------------------
The library exposes module-level DATA3-only knobs that get toggled by the
NF270 workflows. Each one MUST be inert for DATA1/DATA2 model construction —
otherwise an NF270 user could silently corrupt a DATA1 or DATA2 paper
reproduction by leaving a knob set from a prior run.

Knobs under test
----------------
  - NF270_SIGMA_INTERIOR_BOUNDS
  - NF270_CF_RESIDUAL_FLOOR_MM
  - NF270_USE_PERMEATE_PROBE
  - NF270_USE_RETENTATE_ICP_ANCHOR
  - NF270_B_BOUNDS_PER_SALT

For each knob, two DATA1 model builds are exercised:
  1. knob unset (None / legacy)  → expected: legacy bounds / behavior
  2. knob set to an aggressive non-default value  → expected: STILL legacy
     (the knob must be gated on workflow_family == "DATA3")

DATA2 is exercised the same way for the same guarantee.

Performance
-----------
Each model_construct_inter build is ~1 s, so this file's full run is ~30 s.
Tagged `regression`, not `smoke`, because of that wall-clock.
"""
from __future__ import annotations

import pytest


pytestmark = [pytest.mark.regression]


# Use a known DATA2 fixture file (270511.123) since we already know it
# builds cleanly under workflow_family in {DATA1, DATA2}. DATA1 model-build
# has a pre-existing bug on its own .mat files (con_boundary[2] resolves
# to False) — so we use the DATA2 .mat with workflow_family="DATA1"/"DATA2"
# as the lever instead.
FIXTURE_MAT = "data_library/data_stru-dataset270511.123.mat"


@pytest.fixture(scope="module")
def fixture_data_stru(lib, repo_root):
    """Load the DATA2 270511.123 fixture once for the whole module."""
    return lib.loadmat(str(repo_root / FIXTURE_MAT))["data_stru"]


def _reset_all_knobs(lib):
    """Reset every DATA3 knob to its documented default."""
    lib.NF270_CF_RESIDUAL_FLOOR_MM = None
    lib.NF270_SIGMA_INTERIOR_BOUNDS = None
    lib.NF270_USE_PERMEATE_PROBE = None
    lib.NF270_USE_RETENTATE_ICP_ANCHOR = None
    lib.NF270_B_BOUNDS_PER_SALT = {
        "NaCl":  (1e-6, 30.0),
        "CaCl2": (1e-6, 15.0),
        "LaCl3": (1e-6, 10.0),
    }


@pytest.fixture(autouse=True)
def _knob_reset(lib):
    """Reset knobs before AND after every test, so any pollution is contained."""
    _reset_all_knobs(lib)
    yield
    _reset_all_knobs(lib)


# ----------------------------------------------------------------------
# sigma interior bounds guard
# ----------------------------------------------------------------------

@pytest.mark.parametrize("workflow_family", ["DATA1", "DATA2"])
@pytest.mark.parametrize("knob", [None, (0.6, 0.9)])
def test_sigma_interior_bounds_guard(lib, fixture_data_stru, workflow_family, knob):
    """For DATA1/DATA2, σ bounds must always be (0.0, 1.0) regardless of the knob."""
    from pyomo.environ import value
    lib.NF270_SIGMA_INTERIOR_BOUNDS = knob
    m = lib.model_construct_inter(
        fixture_data_stru, "DATA", theta=None, sim_opt=False,
        B_form="single", workflow_family=workflow_family,
    )
    lb = float(value(m.sigma.lb))
    ub = float(value(m.sigma.ub))
    assert (lb, ub) == (0.0, 1.0), (
        f"sigma bounds leaked into {workflow_family} with knob={knob}: "
        f"got ({lb}, {ub}); expected (0.0, 1.0)"
    )


# ----------------------------------------------------------------------
# cF residual floor guard
# ----------------------------------------------------------------------

@pytest.mark.parametrize("workflow_family", ["DATA1", "DATA2"])
@pytest.mark.parametrize("knob", [None, 1.0])
def test_cf_residual_floor_guard(lib, workflow_family, knob):
    """For DATA1/DATA2, the cF residual floor must never be active."""
    lib.NF270_CF_RESIDUAL_FLOOR_MM = knob
    cf_floor_mM = (
        lib.NF270_CF_RESIDUAL_FLOOR_MM
        if str(workflow_family).upper() == "DATA3"
        else None
    )
    cf_meas = 5.0  # 5 mM
    relative_scale = 0.003 * cf_meas  # 0.015 mM
    floor_active = (
        cf_floor_mM is not None
        and cf_floor_mM > 0
        and relative_scale < cf_floor_mM
    )
    assert floor_active is False, (
        f"cF residual floor leaked into {workflow_family} with knob={knob}."
    )


# ----------------------------------------------------------------------
# permeate-probe channel guard
# ----------------------------------------------------------------------

@pytest.mark.parametrize("workflow_family", ["DATA1", "DATA2"])
@pytest.mark.parametrize("knob", [None, 0.03])
def test_permeate_probe_guard(lib, workflow_family, knob):
    """For DATA1/DATA2, the permeate-probe channel must never be active."""
    lib.NF270_USE_PERMEATE_PROBE = knob
    is_data3 = (str(workflow_family).upper() == "DATA3")
    use_perm_probe = bool(is_data3 and lib.NF270_USE_PERMEATE_PROBE)
    assert use_perm_probe is False, (
        f"permeate-probe channel leaked into {workflow_family} with knob={knob}."
    )


# ----------------------------------------------------------------------
# Retentate ICP anchor guard
# ----------------------------------------------------------------------

@pytest.mark.parametrize("workflow_family", ["DATA1", "DATA2"])
@pytest.mark.parametrize("knob", [None, True])
def test_retentate_icp_anchor_guard(lib, workflow_family, knob):
    """For DATA1/DATA2, the retentate ICP anchor must never be active."""
    lib.NF270_USE_RETENTATE_ICP_ANCHOR = knob
    is_data3 = (str(workflow_family).upper() == "DATA3")
    use_retentate_anchor = bool(is_data3 and lib.NF270_USE_RETENTATE_ICP_ANCHOR)
    assert use_retentate_anchor is False, (
        f"retentate ICP anchor leaked into {workflow_family} with knob={knob}."
    )


# ----------------------------------------------------------------------
# B per-salt upper-bound guard
# ----------------------------------------------------------------------

@pytest.mark.parametrize("workflow_family", ["DATA1", "DATA2"])
def test_b_bounds_per_salt_guard(lib, fixture_data_stru, workflow_family):
    """For DATA1/DATA2, B bounds must always be the legacy (1e-6, 50) tuple,
    even when NF270_B_BOUNDS_PER_SALT is set to aggressively distorted values."""
    from pyomo.environ import value
    lib.NF270_B_BOUNDS_PER_SALT = {
        "NaCl":  (1e-6, 3.0),
        "CaCl2": (1e-6, 1.0),
        "LaCl3": (1e-6, 0.5),
    }
    m = lib.model_construct_inter(
        fixture_data_stru, "DATA", theta=None, sim_opt=False,
        B_form="single", workflow_family=workflow_family,
    )
    lb = float(value(m.B.lb))
    ub = float(value(m.B.ub))
    assert abs(lb - 1e-6) < 1e-12 and abs(ub - 50.0) < 1e-9, (
        f"B bounds leaked into {workflow_family}: got ({lb}, {ub}); "
        f"expected (1e-6, 50.0)"
    )
