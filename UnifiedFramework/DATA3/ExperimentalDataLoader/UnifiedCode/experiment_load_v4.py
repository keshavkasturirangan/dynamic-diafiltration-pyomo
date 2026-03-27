# from experiment_dataload_OOP_v14_simpler import load_experiment_easy
from experiment_dataload_OOP_v17_conductivity import load_experiment_easy

# =============================================================================
# To load a .XLSX file
# =============================================================================

specs = {
    "mode": "Lag",        # you said your current XLSX experiments are lag mode
    "delP_bar": 10.0,     # optional override if you know it
    "Temp_K": 298.15,     # optional override if you know it
}

exp, (ok, issues) = load_experiment_easy(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/ExperimentalDataFiles/NF270_MC2.xlsx",
    selector="05.07.24_NaCl",   # sheet name (string) OR selector=0 for first sheet
    specs=specs,
    plot=True,
    plot_kind="signals",        # "signals" plot retentate_signal and (if present) permeate_signal
    )                            # "retentate_signal" plots only retentate
                                # "permeate_signal" plots only permeate
                                # optionally "mass" plots mass

print("OK:", ok)
for sev, msg in issues:
    print(sev, msg)

# =============================================================================
# To load a .MAT file
# =============================================================================
"""
exp, (ok, issues) = load_experiment_easy(
    "data_stru-dataset501.1.mat",
    selector=None,     # defaults to "data_stru"
    specs=None,        # usually not needed for MAT because config is inside the file
    plot=True,
    plot_kind="signals",
)

print("OK:", ok)
for sev, msg in issues:
    print(sev, msg)
"""