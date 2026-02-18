#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simple user-facing runner for the unified v24 workflow."""

from pathlib import Path

from experiment_dataload_OOP_v24 import (
    DiafiltrationExperimentV24,
    ExperimentMode,
    ModelOptions,
    RunMode,
    estimate_parameters_with_parmest_v24,
    load_experiment_easy,
    run_doe_with_pyomo_v24,
    theta_names,
)

# 1) User-selected file and selector
# Absolute path to the experiment data file to load.
path = Path(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/"
    "UnifiedFramework/ExperimentalDataFiles/NF270_MC2.xlsx"
)
# Sheet selector for XLSX (use None for default MAT struct key behavior).
selector = "05.07.24_NaCl"  # XLSX sheet name/index OR None for MAT default key

# 2) Optional overrides (mostly for XLSX)
# These values are applied after reading metadata from the file.
specs = {
    "mode": "Lag",          # Experiment mode used by the process model.
    "delP_bar": 10.0,       # Applied transmembrane pressure [bar].
    "Temp_K": 298.15,       # Temperature [K].
    "Am_cm2": 4.1,          # Active membrane area [cm^2].
    "rho_g_cm3": 1.0,       # Solution density [g/cm^3].
    "C_D_value": 50.0,      # Diafiltrate concentration value.
    "C_D_units": "mM",      # Diafiltrate concentration units.
    "M_F0_g": 200.0,        # Initial feed mass [g].
}

# 3) Workflow switches
# Toggle parameter estimation stage.
RUN_PARMEST = True
# Toggle DoE/FIM stage.
RUN_DOE = True

# 4) Load one experiment
# Main high-level loader that handles format detection and validation.
exp, (ok, issues) = load_experiment_easy(
    str(path),                      # File path string.
    selector=selector,              # XLSX sheet selector.
    specs=specs,                    # User overrides.
    convert_to_concentration=False, # Keep signals as loaded (no conversion here).
    plot=False,                     # Disable quick diagnostic plots in this run.
)
# Report fatal-structure status.
print("Load OK:", ok)
# Print all loader/validator issues (fatal/required/warnings/default-use).
for sev, msg in issues:
    print(f"{sev}: {msg}")

# 5) Build model options
# Settings used by model construction/discretization/estimation.
options = ModelOptions(
    mode=ExperimentMode.DATA,       # Data-mode dynamic balances.
    run_mode=RunMode.ESTIMATION,    # Treat selected parameters as unknowns.
    b_form="single",                # Single B-parameter transport form.
    nfe=30,                         # Finite elements for DAE discretization.
)

# Wrap loaded data in ParmEst/DoE Experiment adapter.
exp_obj = DiafiltrationExperimentV24(exp=exp, options=options)
# Build and label the Pyomo model once for quick diagnostics.
m = exp_obj.get_labeled_model()
# Report parameter names that are currently estimated.
print("Theta names:", theta_names(m, options))
# Number of measurement points exposed to ParmEst.
print("Labeled outputs:", len(m.experiment_outputs))
# Number of unknown parameters exposed to ParmEst.
print("Labeled unknowns:", len(m.unknown_parameters))
# Number of design inputs exposed to DoE.
print("Labeled design inputs:", len(m.experiment_inputs))

# 6) Parameter estimation
if RUN_PARMEST:
    # Run ParmEst with built-in SSE_weighted objective.
    est = estimate_parameters_with_parmest_v24(
        experiments=[exp],                         # Single-experiment estimation.
        options=options,                           # Model options.
        calc_cov=True,                             # Attempt covariance matrix.
        cov_n=max(1, len(m.experiment_outputs)),   # Deprecated passthrough note.
        solver="ipopt",                            # ParmEst solver.
        tee=False,                                 # Hide solver stream output.
    )
    # Generic warning field from estimation stage.
    if "warning" in est:
        print("WARNING:", est["warning"])
    # Covariance-specific warning with method-level failure diagnostics.
    if "covariance_warning" in est:
        print("WARNING:", est["covariance_warning"])
    # Estimated objective function value.
    print("ParmEst objective:", est.get("objective"))
    # Estimated parameters.
    print("ParmEst theta:\n", est.get("theta"))
    # Covariance result (if a method converged).
    if "covariance" in est:
        print("ParmEst covariance method:", est.get("covariance_method"))
        print("ParmEst covariance:\n", est.get("covariance"))

# 7) DoE/FIM analysis
if RUN_DOE:
    # Run local DoE/FIM analysis on the same experiment wrapper.
    doe = run_doe_with_pyomo_v24(exp_obj, solver_name="ipopt", tee=False)
    # Fisher information matrix at current settings.
    print("DoE FIM:\n", doe["fim"])
    # D-optimality score (log-det of FIM).
    print("DoE D-opt (log-det):", doe["d_opt_logdet"])
