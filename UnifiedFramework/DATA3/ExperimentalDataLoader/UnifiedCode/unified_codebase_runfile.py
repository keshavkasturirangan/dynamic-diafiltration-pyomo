#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simple user-facing runner for the unified codebase pipeline."""

from pathlib import Path

try:
    from .unified_codebase_library import (
        ExperimentMode,
        ModelOptions,
        RunMode,
        UnifiedPipelineConfigV24,
        run_unified_pipeline_v24,
    )
except ImportError:
    from unified_codebase_library import (
        ExperimentMode,
        ModelOptions,
        RunMode,
        UnifiedPipelineConfigV24,
        run_unified_pipeline_v24,
    )

# 1) User-selected file and selector
path = Path(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/"
    "UnifiedFramework/ExperimentalDataFiles/NF270_MC2.xlsx"
)
selector = "05.07.24_NaCl"  # XLSX sheet selector, or None for MAT default.

# 2) Optional metadata overrides (mainly for XLSX files)
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

# 3) Shared model options (used for MAT and XLSX)
model_options = ModelOptions(
    mode=ExperimentMode.DATA,
    run_mode=RunMode.ESTIMATION,
    b_form="single",
    nfe=30,
)

# 4) Unified pipeline execution config
config = UnifiedPipelineConfigV24(
    file_path=str(path),
    selector=selector,
    specs=specs,
    convert_to_concentration=False,
    plot=False,
    model_options=model_options,
    run_parmest=True,
    run_doe=True,
    calc_cov=True,
    cov_n=None,
    solver="ipopt",
    solver_options=None,
    tee=False,
)

# 5) Execute one unified run
result = run_unified_pipeline_v24(config)

print("Load OK:", result["load_ok"])
for sev, msg in result["issues"]:
    print(f"{sev}: {msg}")

print("Theta names:", result["theta_names"])
print("Labeled outputs:", result["labeled_counts"]["outputs"])
print("Labeled unknowns:", result["labeled_counts"]["unknowns"])
print("Labeled design inputs:", result["labeled_counts"]["inputs"])

if "parmest" in result:
    est = result["parmest"]
    if "warning" in est:
        print("WARNING:", est["warning"])
    if "covariance_warning" in est:
        print("WARNING:", est["covariance_warning"])
    print("ParmEst objective:", est.get("objective"))
    print("ParmEst theta:\n", est.get("theta"))
    if "covariance" in est:
        print("ParmEst covariance method:", est.get("covariance_method"))
        print("ParmEst covariance:\n", est.get("covariance"))

if "doe" in result:
    doe = result["doe"]
    print("DoE FIM:\n", doe["fim"])
    print("DoE D-opt (log-det):", doe["d_opt_logdet"])
