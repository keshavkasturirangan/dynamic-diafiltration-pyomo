# =========================================================================== #
# =========================================================================== #
# =========================================================================== #

from experiment_dataload_OOP_v17_conductivity import (
    load_experiment_easy, build_diafiltration_model,
    ModelOptions, ParameterGuess, ExperimentMode, RunMode, BForm
)
import numpy as np

# =============================================================================
# To load a .XLSX file (one sheet = one experiment)
# =============================================================================

# Optional: experiment-level overrides you *know* are correct.
# These override whatever is parsed from the file.
specs = {
    "mode": "Lag",        # current XLSX experiments are lag mode
    "delP_bar": 10.0,     # applied pressure [bar] (override only if you know it's correct)
    "Temp_K": 298.15,     # temperature [K] (override only if you know it's correct)
}

model_params = {
    # common
    "temp_K": 298.15,

    # Shedlovsky parameters (example placeholders — must be set consistently with the paper code)
    "epsilon": 78.3,
    "eta": 0.0089,             # (watch units; must match conductivity_paper.py)
    "lambda_0": 1.0,
    "a": 4e-8,
    "z_1": 1,
    "z_2": -1,
    "lambda_0_cation": 50.0,
    "lambda_0_anion": 50.0,
}

exp, (ok, issues) = load_experiment_easy(
    "NF270_MC2.xlsx",
    selector="Experiment_3",
    specs={"mode": "Lag", "Temp_K": 298.15},
    convert_to_concentration=True,
    conductivity_model_params=model_params,
    n_cations=1,          # or 2 => Shedlovsky; >=3 => MSA
    conc_units="mM",
    plot=True,
    plot_kind="retentate_signal",
)
options = ModelOptions(mode=ExperimentMode.LAG, run_mode=RunMode.ESTIMATION, B_form=BForm.SINGLE)

guess = ParameterGuess(
    Lp=10.0,
    sigma=0.5,
    theta0=np.array([10.0, 1.0, 0.5], dtype=float),
    B=1.0
)

m = build_diafiltration_model(exp, options, guess)

print("OK:", ok)
for sev, msg in issues:
    print(sev, msg)

# =============================================================================
# To load a .MAT file (one file = one experiment)
# =============================================================================
# Uncomment to run a MAT experiment; selector=None defaults to "data_stru".
#
# exp, (ok, issues) = load_experiment_easy(
#     "data_stru-dataset501.1.mat",
#     selector=None,     # defaults to "data_stru"
#     specs=None,        # usually not needed for MAT because config is inside the file
#     plot=True,
#     plot_kind="retentate_signal",
# )
#
# print("OK:", ok)
# for sev, msg in issues:
#     print(sev, msg)
