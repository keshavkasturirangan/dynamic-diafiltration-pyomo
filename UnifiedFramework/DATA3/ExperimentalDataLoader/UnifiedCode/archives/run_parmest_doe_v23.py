#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unified v23 execution script: dataloader -> process model -> ParmEst/DoE."""

from pathlib import Path

from experiment_dataload_OOP_v23_conductivity_patched import load_experiment_easy
from process_model_builder_v23 import (
    DiafiltrationExperimentV23,
    ExperimentMode,
    ModelOptions,
    RunMode,
    estimate_parameters_with_parmest_v23,
    run_doe_with_pyomo_v23,
    theta_names,
)


# ---------------------------------------------------------------------
# 1) User input file and selector
# ---------------------------------------------------------------------
path = Path(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/"
    "UnifiedFramework/ExperimentalDataFiles/NF270_MC2.xlsx"
)
selector = "05.07.24_NaCl"  # XLSX sheet name (or int index)

# ---------------------------------------------------------------------
# 2) Optional experiment overrides
# ---------------------------------------------------------------------
specs = {
    "mode": "Lag",
    "delP_bar": 10.0,
    "Temp_K": 298.15,
    "Am_cm2": 4.1,
    "rho_g_cm3": 1.0,
    "C_D_value": 50.0,
    "C_D_units": "mM",
    "M_F0_g": 200.0,
}

# ---------------------------------------------------------------------
# 3) Workflow toggles
# ---------------------------------------------------------------------
RUN_PARMEST = True
RUN_DOE = True

# ---------------------------------------------------------------------
# 4) Load experiment
# ---------------------------------------------------------------------
exp, (ok, issues) = load_experiment_easy(
    str(path),
    selector=selector,
    specs=specs,
    convert_to_concentration=False,
    plot=False,
)

print("Load OK:", ok)
for sev, msg in issues:
    print(f"{sev}: {msg}")

# ---------------------------------------------------------------------
# 5) Build labeled experiment for ParmEst/DoE
# ---------------------------------------------------------------------
options = ModelOptions(
    mode=ExperimentMode.DATA,
    run_mode=RunMode.ESTIMATION,
    b_form="single",
    nfe=60,
)

exp_obj = DiafiltrationExperimentV23(exp=exp, options=options)
m = exp_obj.get_labeled_model()

print("Theta names:", theta_names(m, options))
print("Labeled outputs:", len(m.experiment_outputs))
print("Labeled unknowns:", len(m.unknown_parameters))
print("Labeled design inputs:", len(m.experiment_inputs))

# ---------------------------------------------------------------------
# 6) Optional ParmEst run
# ---------------------------------------------------------------------
if RUN_PARMEST:
    est = estimate_parameters_with_parmest_v23(
        experiments=[exp],
        options=options,
        calc_cov=True,
        cov_n=max(1, len(m.experiment_outputs)),
        solver="ipopt",
        tee=False,
    )
    print("ParmEst objective:", est.get("objective"))
    print("ParmEst theta:\n", est.get("theta"))
    if "covariance" in est:
        print("Covariance:\n", est["covariance"])

# ---------------------------------------------------------------------
# 7) Optional DoE run
# ---------------------------------------------------------------------
if RUN_DOE:
    doe = run_doe_with_pyomo_v23(exp_obj, solver_name="ipopt", tee=False)
    print("DoE FIM:\n", doe["fim"])
    print("DoE D-opt (log-det):", doe["d_opt_logdet"])
