#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Minimal runner script for the unified experiment loader.

This file is intentionally *simple*:
    - choose an input file (.xlsx or .mat)
    - choose selector (sheet name/index for XLSX; struct key for MAT)
    - optionally provide rig/spec overrides (needed for XLSX, usually not for MAT)
    - optionally convert conductivity->concentration (XLSX only)
    - optionally plot a signal for quick sanity-checks
"""

from pathlib import Path

# Import the *unified* loader entrypoint.
# NOTE: this import path assumes this script lives next to the loader file.
from experiment_dataload_OOP_v17_conductivity_workflow import load_experiment_easy


# ---------------------------------------------------------------------
# 1) Choose ONE file to load
# ---------------------------------------------------------------------
# Example A (XLSX): multi-experiment workbook where each sheet is one experiment.
path = Path(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/"
    "UnifiedFramework/ExperimentalDataFiles/NF270_MC2.xlsx"
)
selector = "05.07.24_NaCl"   # sheet name OR selector=0 for first sheet

# Example B (MAT): single-experiment file. Uncomment to use.
# path = Path("/path/to/data_stru-dataset501.1.mat")
# selector = None  # defaults to "data_stru"


# ---------------------------------------------------------------------
# 2) Rig/operating specs (ONLY needed for XLSX; MAT stores these in data_config)
# ---------------------------------------------------------------------
xlsx_specs = {
    # Experiment mode (your current XLSX set is lag mode)
    "mode": "Lag",

    # Operating conditions used in the process model
    "delP_bar": 10.0,       # applied pressure [bar]
    "Temp_K": 298.15,       # temperature [K]
    "Am_cm2": 4.1,          # effective membrane area [cm^2]
    "rho_g_cm3": 1.0,       # density [g/cm^3]
}

# Decide whether to apply specs (based purely on file extension)
suffix = path.suffix.lower()
specs = xlsx_specs if suffix in [".xlsx", ".xls"] else None


# ---------------------------------------------------------------------
# 3) Optional: conductivity -> concentration conversion (XLSX only)
# ---------------------------------------------------------------------
# If you set convert_to_concentration=True, the loader will:
#   - read conductivity time-series (retentate + permeate)
#   - choose a physics model:
#       auto: variant Shedlovsky if <=2 cations, else MSA
#   - write concentrations back onto each vial (fields on VialData)
convert_to_concentration = suffix in [".xlsx", ".xls"]

# Model parameters required by conductivity_paper.py.
# Keep this dict explicit so it's obvious what is being assumed.
conductivity_model_params = {
    # Shared / physical constants
    "epsilon": 78.3,          # [-] relative permittivity (water ~78 at room temp)
    "eta": 0.00089,           # [Pa*s] dynamic viscosity (water ~0.89 mPa*s)
    "lambda_0": 50.0,         # [S cm^2 / mol] limiting molar conductivity (placeholder)
    # Variant Shedlovsky parameters (1–2 cations)
    "a": 1.0,                 # [dimensionless] activity/empirical constant
    "z_1": 1,                 # [-] cation charge
    "z_2": -1,                # [-] anion charge
    "lambda_0_cation": 50.0,  # [S cm^2 / mol]
    "lambda_0_anion": 50.0,   # [S cm^2 / mol]
    # MSA parameters (>=3 cations) would be added here if needed
}


# ---------------------------------------------------------------------
# 4) Load (file-agnostic)
# ---------------------------------------------------------------------
exp, (ok, issues) = load_experiment_easy(
    str(path),
    selector=selector,
    specs=specs,
    convert_to_concentration=convert_to_concentration,
    concentration_units="mM",
    conductivity_model="auto",
    conductivity_model_params=conductivity_model_params,
    plot=True,
    plot_kind="retentate_signal",   # try: "signals", "permeate_signal", "mass_g"
)

print("OK:", ok)
for sev, msg in issues:
    print(f"{sev} {msg}")
