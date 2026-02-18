#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simple user-facing runner for the unified v23 workflow."""

from pathlib import Path

from experiment_dataload_OOP_v23 import (
    DiafiltrationExperimentV23,
    ExperimentMode,
    ModelOptions,
    RunMode,
    estimate_parameters_with_parmest_v23,
    load_experiment_easy,
    run_doe_with_pyomo_v23,
    theta_names,
)

# 1) User-selected file and selector
path = Path(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/"
    "UnifiedFramework/ExperimentalDataFiles/NF270_MC2.xlsx"
)
selector = "05.07.24_NaCl"  # XLSX sheet name/index OR None for MAT default key

# 2) Optional overrides (mostly for XLSX)
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

# 3) Workflow switches
RUN_PARMEST = True
RUN_DOE = True

# 4) Load one experiment
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

# 5) Build model options
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

# 6) Parameter estimation
if RUN_PARMEST:
    est = estimate_parameters_with_parmest_v23(
        experiments=[exp],
        options=options,
        calc_cov=True,
        cov_n=max(1, len(m.experiment_outputs)),
        solver="ipopt",
        tee=False,
    )
    if "warning" in est:
        print("WARNING:", est["warning"])
    print("ParmEst objective:", est.get("objective"))
    print("ParmEst theta:\n", est.get("theta"))

# 7) DoE/FIM analysis
if RUN_DOE:
    doe = run_doe_with_pyomo_v23(exp_obj, solver_name="ipopt", tee=False)
    print("DoE FIM:\n", doe["fim"])
    print("DoE D-opt (log-det):", doe["d_opt_logdet"])
