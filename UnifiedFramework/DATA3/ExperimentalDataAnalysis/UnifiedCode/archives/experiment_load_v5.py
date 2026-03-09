


# =============================================================================
# =============================================================================
# =============================================================================

from process_model_v15 import (
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

exp, (ok, issues) = load_experiment_easy(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/ExperimentalDataFiles/NF270_MC2.xlsx",
    selector="05.07.24_NaCl",   # sheet name (string) OR selector=0 for first sheet
    specs=specs,
    plot=True,
    # NOTE:
    #   Your current plot_experiment implementation recognizes specific kinds
    #   like "retentate_signal", "permeate_signal", and "mass".
    #   "signals" is NOT currently implemented in the loader, so we call a
    #   supported plot kind here.
    plot_kind="mass",
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
