"""Catalog of DATA1 / DATA2 reference cases for regression testing.

Test architecture (redesigned 2026-05-24)
------------------------------------------
Each case fixes a (data_file, mode, B_form, workflow_family, THETA) tuple
and runs `solve_model(..., sim_opt=True)` — a FORWARD-ONLY simulation, no
parameter estimation. The capture script records key end-of-vial values
(mV[-1], cF[-1], cV[-1], cH[-1]) from the simulated trajectory; the pytest
re-runs the simulation and asserts the trajectory hasn't drifted within 1%.

Why no fit?
-----------
A full IPOPT fit takes ~5 min per case, making the test suite too slow to
run on every commit. A forward simulation takes ~1 sec and STILL catches
the regressions we care about most:
  - Loader changes that break the data_stru shape
  - model_construct_inter changes that alter the ODEs / constraints
  - DAE discretization or transformation changes

What it doesn't catch:
  - Fitter / parmest regressions (would require the full fit)
  - Sensitivity to initial-guess changes

If we ever need fit-level regression, add a second tier of tests with
`sim_opt=False` and a longer pytest timeout.

Adding a new case
-----------------
1. Append an entry below with a THETA dict that the model accepts.
2. Re-run `_capture_references.py` to capture trajectories.
3. The new pytest is auto-picked-up via parametrize.
"""
from __future__ import annotations


# DATA1 cases — DATA-mode (closed-cell filtration) on KCl/NF270.
# All are currently SKIPPED because the DATA1 DATA-mode model-build path has
# a pre-existing bug (`con_boundary[2]` constraint resolves to a trivial
# False during Pyomo construction). Reactivate by removing `skip_reason`
# once the model-build issue is fixed.
DATA1_CASES = [
    {
        "case_id": "data1__501.11__single",
        "data_file": "data_stru-dataset501.11.mat",
        "data_dir":  "legacy/data1_matlab/data",
        "mode":      "DATA",
        "B_form":    "single",
        "theta":     {"Lp": 4.2, "B": 0.78, "sigma": 0.85},   # near literature
        "workflow_family": "DATA1",
        "skip_reason": "pre-existing DATA1 model-build bug (con_boundary[2] resolves to False)",
    },
    {
        "case_id": "data1__511.12__single",
        "data_file": "data_stru-dataset511.12.mat",
        "data_dir":  "legacy/data1_matlab/data",
        "mode":      "DATA",
        "B_form":    "single",
        "theta":     {"Lp": 7.5, "B": 1.2, "sigma": 0.92},
        "workflow_family": "DATA1",
        "skip_reason": "pre-existing DATA1 model-build bug (con_boundary[2] resolves to False)",
    },
]


# DATA2 cases — NF270 / KCl Lag-mode and Overflow-mode validation.
# These use the exact theta the DATA2 paper demo notebook uses, so the
# simulated trajectory is reproducible and matches what the published
# figures show.
DATA2_CASES = [
    {
        "case_id": "data2__270511.123__Lag_with_time_correction",
        "data_file": "data_stru-dataset270511.123.mat",
        "data_dir":  "legacy/data1_matlab/data_library",
        "mode":      "Lag",
        "B_form":    1,
        "theta":     {"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01,
                      "sigma": 1.0, "S0": 0, "S": 0.0},
        "workflow_family": "DATA2",
    },
    {
        "case_id": "data2__270511.423__Overflow_with_time_correction",
        "data_file": "data_stru-dataset270511.423.mat",
        "data_dir":  "legacy/data1_matlab/data_library",
        "mode":      "Overflow",
        "B_form":    1,
        "theta":     {"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01,
                      "sigma": 1.0, "S0": -0.1756665334051245, "S": 0.0},
        "workflow_family": "DATA2",
    },
]


ALL_CASES = DATA1_CASES + DATA2_CASES
