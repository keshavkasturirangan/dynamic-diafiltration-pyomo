#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 07:12:33 2025

@author: Keshav Kasturi Rangan, Alex Dowling
"""

from pathlib import Path

from experiment_dataload_OOP_v21_conductivity_patched import load_experiment_easy


# ---------------------------------------------------------------------
# 1) User choices: pick ONE file path and (if XLSX) a sheet selector
# ---------------------------------------------------------------------
path = Path(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/"
    "UnifiedFramework/ExperimentalDataFiles/NF270_MC2.xlsx"
)

# For XLSX: sheet name or index
selector = "05.07.24_NaCl"   # or selector = 0

# For MAT: selector is usually None (defaults to "data_stru")
# selector = None


# ---------------------------------------------------------------------
# 2) Rig/operating specs (ONLY needed for XLSX because those fields are not stored)
# ---------------------------------------------------------------------
xlsx_specs = {
    "mode": "Lag",          # your current XLSX experiments are lag mode
    "delP_bar": 10.0,       # applied pressure [bar]
    "Temp_K": 298.15,       # temperature [K]
    "Am_cm2": 4.1,          # membrane area [cm^2]
    "rho_g_cm3": 1.0,       # density [g/cm^3]
}


# ---------------------------------------------------------------------
# 3) Decide whether to pass specs based on file extension
# ---------------------------------------------------------------------
suffix = path.suffix.lower()

if suffix in [".xlsx", ".xls"]:
    specs = xlsx_specs
elif suffix == ".mat":
    specs = None  # MAT should already contain these in data_config
else:
    raise ValueError(f"Unsupported file type: {suffix}")


# ---------------------------------------------------------------------
# 4) Load (file-agnostic)
# ---------------------------------------------------------------------
exp, (ok, issues) = load_experiment_easy(
    str(path),
    selector=selector,
    specs=specs,

    # -----------------------------------------------------------------
    # Conductivity -> concentration conversion (XLSX only)
    # -----------------------------------------------------------------
    # Set convert_to_concentration=True to fill v.retentate_concentration and
    # v.permeate_concentration for each vial. MAT inputs typically already
    # contain concentration from legacy workflows, so conversion is skipped.
    convert_to_concentration=True,
    concentration_units="mM",
    conductivity_model="auto",   # auto: 1-2 cations -> Shedlovsky, >=3 -> MSA
    n_cations=1,                   # set explicitly when you know it (avoids parsing salt names)
    conductivity_model_params={
        # Pass through to conductivity_paper-based models (values are examples).
        # Keep these consistent with conductivity_paper.py expectations.
        "temp_K": 298.15,
        "eta": 0.00089,
        "epsilon": 78.3,
        "lambda_0_cation": 50.0,
        "lambda_0_anion": 50.0,
    },

    plot=True,
    plot_kind="signals",
)

print("OK:", ok)
for sev, msg in issues:
    print(sev, msg)
