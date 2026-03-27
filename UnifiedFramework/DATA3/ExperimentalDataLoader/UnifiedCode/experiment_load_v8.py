"""Example runner (no CLI): load one experiment, optionally convert conductivity -> concentration,
and then build a diafiltration Pyomo model.

This file is intentionally *simple*:
    - edit the PATH + selector for your use case,
    - tweak a small 'specs' dict if you know experiment-level metadata,
    - call the loader.

Notes:
    - XLSX experiments store raw conductivity signals (uS/cm).
      If you want concentrations for modeling, call apply_conductivity_to_concentration(...).
    - MAT experiments already contain concentration data (per legacy pipeline). We don't reconvert.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# 1) Import the file-agnostic data loader (XLSX + MAT)
# ---------------------------------------------------------------------------
from experiment_dataload_OOP_v16_conductivity import (  # noqa: E402
    load_experiment_easy,
    apply_conductivity_to_concentration,
)

# ---------------------------------------------------------------------------
# 2) Import the process-model builder (kept separate from the loader)
# ---------------------------------------------------------------------------
from process_model_v15 import (  # noqa: E402
    build_diafiltration_model,
    ModelOptions,
    ParameterGuess,
    ExperimentMode,
    RunMode,
    BForm,
)


# =============================================================================
# A) Load a .XLSX file (one sheet = one experiment)
# =============================================================================

# Path to the Excel workbook.
xlsx_path = "NF270_MC2.xlsx"

# Sheet selection:
#   - use a string sheet name (e.g., "Experiment_3"), OR
#   - use an integer index (0 for first sheet).
sheet_selector = "Experiment_3"

# Optional overrides you *know* are correct.
# These are applied after parsing, regardless of input file type.
specs = {
    "mode": "Lag",        # current XLSX experiments are lag mode
    "Temp_K": 298.15,     # temperature [K] (override only if you know it's correct)
    # If your XLSX does NOT contain dialysate concentration, provide it here for modeling.
    # "C_D_value": 0.0,
    # "C_D_units": "mM",
}

# Load the experiment (and show a quick plot of signals).
exp, (ok, issues) = load_experiment_easy(
    xlsx_path,
    selector=sheet_selector,
    specs=specs,
    plot=True,
    plot_kind="signals",  # plots retentate + permeate conductivity if permeate exists
)

# ---------------------------------------------------------------------------
# Optional: Convert conductivity (uS/cm) -> concentration (mM or M)
# ---------------------------------------------------------------------------
# IMPORTANT:
#   - This uses conductivity_paper.py as-is.
#   - The model choice is automatic by your rule:
#       * 1–2 cations -> variant Shedlovsky
#       * 3+ cations  -> MSA
#   - exp.component_names must be populated (XLSX parsing should set this from "Salt 1", ...)

conductivity_model_params = {
    # NOTE: These are placeholders; set them to values consistent with conductivity_paper.py.
    "epsilon": 78.3,
    "eta": 0.0089,               # viscosity units MUST match conductivity_paper.py expectations
    "lambda_0": 1.0,
    "a": 4e-8,
    "z_1": 1,
    "z_2": -1,
    "lambda_0_cation": 50.0,
    "lambda_0_anion": 50.0,

    # MSA-only parameters (only used if auto-select chooses MSA)
    # "valency": [...],
    # "diameters": [...],
    # "diff_coeff": [...],
}

# Convert and store into each VialData:
#   v.retentate_concentration, v.permeate_concentration (arrays aligned with time_s)
apply_conductivity_to_concentration(
    exp,
    output_units="mM",
    model="auto",
    model_params=conductivity_model_params,
)


# =============================================================================
# B) Build the Pyomo model (example)
# =============================================================================

options = ModelOptions(
    mode=ExperimentMode.LAG,
    run_mode=RunMode.ESTIMATION,
    B_form=BForm.SINGLE,
)

guess = ParameterGuess(
    Lp=10.0,
    sigma=0.5,
    theta0=np.array([10.0, 1.0, 0.5], dtype=float),
    B=1.0,
)

# Build model (will raise if required experiment-level fields are missing)
m = build_diafiltration_model(exp, options, guess)


print("OK:", ok)
for sev, msg in issues:
    print(sev, msg)


# =============================================================================
# C) Load a .MAT file (one file = one experiment)
# =============================================================================
# NOTE: MAT files already contain concentrations in legacy datasets.
# Uncomment to test:
#
# mat_path = "data_stru-dataset501.1.mat"
# exp, (ok, issues) = load_experiment_easy(
#     mat_path,
#     selector=None,  # defaults to "data_stru"
#     plot=True,
#     plot_kind="retentate_signal",
# )
#
# print("OK:", ok)
# for sev, msg in issues:
#     print(sev, msg)
