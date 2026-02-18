# from experiment_dataload_OOP_v14_simpler import load_experiment_easy

from experiment_dataload_OOP_v17_conductivity import load_experiment_easy

# ---------------------------------------------------------------------
# Specs you KNOW (rig / operating conditions) and that XLSX does not store
# ---------------------------------------------------------------------
specs = {
    "mode": "Lag",          # all your current XLSX experiments are lag mode
    "delP_bar": 10.0,       # applied pressure [bar]
    "Temp_K": 298.15,       # temperature [K]
    "Am_cm2": 4.1,          # membrane area [cm^2]
    "rho_g_cm3": 1.0,       # density [g/cm^3] (water ~1.0)
}

# Load one experiment from Excel
exp, (ok, issues) = load_experiment_easy(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/ExperimentalDataFiles/NF270_MC2.xlsx",
    selector="05.07.24_NaCl",   # or selector=0 for first sheet
    specs=specs,
    plot=True,
    plot_kind="signals",
)

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