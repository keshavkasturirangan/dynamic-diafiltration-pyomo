"""
First-principles conductivity calibration for DATA3 (docs/prompts/
conductivity_model_calibration.md), using the vendored Shedlovsky and MSA
models (analysis/conductivity_models/conductivity.py) from the
Lilonfe/Estrada/Singh/Ouimet/Phillip/Dowling soft-sensor manuscript
(JMS submission MEMSCI-S-26-01544-2, 2026).

Four tasks, run in order:
  1. Generate kappa(c) (MSA primary, Shedlovsky secondary) for NaCl/CaCl2/
     LaCl3 over the DATA3 diafiltration ranges, at the experimental
     retentate temperature and at 298.15 K (manuscript reference T).
  2. Linearity test: fit linear + quadratic kappa(c) over each range; report
     slope/intercept/R^2/curvature/max deviation of the linear fit.
  3. Adjudicate discussion point 9: overlay model kappa(c) against the DATA3
     authors' inline calibration and the Phase-1 vial-ICP fit, per salt.
  4. Reproduce the manuscript's own Table 3 MAPEs (implementation check),
     using the manuscript's experimental data
     (salts-conductivity-modelling/Data/{NaCl,CaCl2,LaCl3}.csv).

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/conductivity_models/calibration_curves.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "preprocessing"))
import conductivity as cd  # noqa: E402  (vendored, unmodified)
import preprocess_nacl as pp  # noqa: E402  (for raw retentate temperatures)
import plot_style  # noqa: E402

DATA_DIR = REPO_ROOT / "data"
FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
# Manuscript's own experimental data (separate sibling repo -- not copied
# into this repo; referenced read-only for Task 4's implementation check).
MANUSCRIPT_DATA_DIR = Path(
    "/Users/adowling/DowlingLab/Membranes/salts-conductivity-modelling/Data"
)

# ---------------------------------------------------------------------------
# Manuscript Table 1: Shedlovsky parameters (T=298.15 K, eps_r=78.43,
# eta=8.9e-3 poise, z_Cl=-1, conc in M)
# ---------------------------------------------------------------------------
T_MANUSCRIPT = 298.15
EPSILON_R = 78.43
ETA_SHED_POISE = 8.9e-3
Z_CL = -1
LAMBDA0_CL = 76.35  # cm^2.S/equiv

SHED_PARAMS = {
    "NaCl": dict(lambda_0=126.45, lambda_0_cation=50.10, a=4.0e-8, z_cation=1),
    "CaCl2": dict(lambda_0=135.85, lambda_0_cation=59.50, a=4.3e-8, z_cation=2),
    "LaCl3": dict(lambda_0=145.90, lambda_0_cation=69.70, a=4.9e-8, z_cation=3),
}
# Shedlovsky validity: ionic strength (equivalent, N) <= 0.1. c_eq = z*c_salt.
SHED_VALID_N = 0.1

# ---------------------------------------------------------------------------
# Manuscript Table 2: MSA parameters (conc in mM, eta=0.89e-3 Pa.s, order
# [cation, Cl])
# ---------------------------------------------------------------------------
ETA_MSA_PAS = 0.89e-3
MSA_IONS = {
    "Na": dict(lambda_0_mol=50.08e-4, sigma=2.04e-10, D=1.33e-9, z=1),
    "Ca": dict(lambda_0_mol=59.47e-4, sigma=2.00e-10, D=0.79e-9, z=2),
    "La": dict(lambda_0_mol=69.70e-4, sigma=2.36e-10, D=0.62e-9, z=3),
    "Cl": dict(lambda_0_mol=76.31e-4, sigma=3.62e-10, D=2.03e-9, z=-1),
}
MSA_SALT_CATION = {"NaCl": "Na", "CaCl2": "Ca", "LaCl3": "La"}

# ---------------------------------------------------------------------------
# DATA3 diafiltration c_int ranges (from the analyzed/reprocessed datasets;
# also sampled near c->0 as requested)
# ---------------------------------------------------------------------------
DATA3_RANGES_MM = {
    "NaCl": (1e-3, 180.0),
    "CaCl2": (1e-3, 108.0),
    "LaCl3": (1e-3, 28.0),
}
N_POINTS = 300

# Experimental retentate temperature source per salt (the run analyzed in
# sec:nacl/sec:cacl2/sec:lacl3 and reprocessed in Phase 1)
T_EXP_SOURCE = {
    "NaCl": ("NF270_MC5.xlsx", "05.27.26_NaCl"),
    "CaCl2": ("NF270_MC3.xlsx", "07.11.24_SCaCl2"),
    "LaCl3": ("NF270_MC2.xlsx", "05.21.24_LaCl3"),
}

# Authors' inline calibration, DATA3 authors' fits, kappa [uS/cm] = a + s*c[mM]
AUTHORS_INLINE = {
    "NaCl": dict(a=-63.706, s=76.685, label="NaCl authors' inline (MC5 05.27.26 embedded)"),
    "CaCl2": dict(a=-132.6, s=152.1, label="CaCl2 authors' inline"),
    "LaCl3": dict(a=-302.6, s=216.3, label="LaCl3 authors' inline"),
}
# Phase-1 vial-ICP fits (analysis/preprocessing/README.md)
VIAL_ICP_FIT = {
    "NaCl": dict(a=111.334, s=93.004, label="NaCl vial-ICP (MC5 07.23.24, representative raw-only run)"),
    "CaCl2": dict(a=96.406, s=187.577, label="CaCl2 vial-ICP (MC3 07.11.24_SCaCl2)"),
    "LaCl3": dict(a=93.837, s=260.084, label="LaCl3 vial-ICP (MC2 05.21.24_LaCl3)"),
}


def get_experimental_temperature_K(salt):
    file, sheet = T_EXP_SOURCE[salt]
    raw = pp.load_raw_sheet(DATA_DIR / file, sheet)
    t_c = raw.ts["ret_temp"].dropna().mean()
    return t_c + 273.15


def kappa_shedlovsky_uS_cm(conc_mM, salt, T):
    """conc_mM: array-like, mM. Returns (kappa_uS_cm, is_extrapolated)."""
    conc_mM = np.asarray(conc_mM, dtype=float)
    conc_M = list(conc_mM * 1e-3)
    p = SHED_PARAMS[salt]
    cond_mS_cm = cd.shedlovsky(
        conc_M, T, EPSILON_R, ETA_SHED_POISE, p["lambda_0"], p["a"],
        p["z_cation"], Z_CL, p["lambda_0_cation"], LAMBDA0_CL,
    )
    kappa_uS_cm = np.asarray(cond_mS_cm) * 1000.0
    c_eq_N = np.abs(p["z_cation"]) * conc_mM * 1e-3  # equivalents/L
    is_extrapolated = c_eq_N > SHED_VALID_N
    return kappa_uS_cm, is_extrapolated


def kappa_msa_uS_cm(conc_mM, salt, T):
    conc_mM = np.asarray(conc_mM, dtype=float)
    cation = MSA_SALT_CATION[salt]
    cat, cl = MSA_IONS[cation], MSA_IONS["Cl"]
    valency = [cat["z"], cl["z"]]
    diameters = [cat["sigma"], cl["sigma"]]
    diff_coeff = [cat["D"], cl["D"]]
    lambda_0 = [cat["lambda_0_mol"], cl["lambda_0_mol"]]
    cond_mS_cm = cd.msa(valency, diameters, diff_coeff, T, ETA_MSA_PAS, EPSILON_R,
                         lambda_0, list(conc_mM))
    return np.asarray(cond_mS_cm) * 1000.0


# ---------------------------------------------------------------------------
# Task 1: generate kappa(c) at experimental T and at 298.15 K
# ---------------------------------------------------------------------------
def task1_generate_curves():
    print("=" * 100)
    print("TASK 1: kappa(c) generation (MSA primary, Shedlovsky secondary)")
    print("=" * 100)
    curves = {}
    for salt, (c_lo, c_hi) in DATA3_RANGES_MM.items():
        c = np.linspace(c_lo, c_hi, N_POINTS)
        T_exp = get_experimental_temperature_K(salt)
        msa_exp = kappa_msa_uS_cm(c, salt, T_exp)
        msa_298 = kappa_msa_uS_cm(c, salt, T_MANUSCRIPT)
        shed_exp, extrap_exp = kappa_shedlovsky_uS_cm(c, salt, T_exp)
        shed_298, extrap_298 = kappa_shedlovsky_uS_cm(c, salt, T_MANUSCRIPT)
        n_extrap = int(np.sum(extrap_298))
        print(f"\n{salt}: range {c_lo:.3g}-{c_hi:.3g} mM, T_exp={T_exp:.2f} K "
              f"({T_exp-273.15:.2f} C)")
        print(f"  MSA(T_exp)  kappa: {msa_exp[0]:.4g} -> {msa_exp[-1]:.4g} uS/cm")
        print(f"  MSA(298.15) kappa: {msa_298[0]:.4g} -> {msa_298[-1]:.4g} uS/cm")
        print(f"  Shedlovsky(298.15) kappa: {shed_298[0]:.4g} -> {shed_298[-1]:.4g} uS/cm "
              f"({n_extrap}/{len(c)} points exceed 0.1 N, extrapolated)")
        curves[salt] = dict(c=c, T_exp=T_exp, msa_exp=msa_exp, msa_298=msa_298,
                             shed_exp=shed_exp, shed_298=shed_298, extrap_298=extrap_298)
    return curves


# ---------------------------------------------------------------------------
# Task 2: linearity test
# ---------------------------------------------------------------------------
def linearity_stats(c, kappa):
    lin = np.polyfit(c, kappa, 1)
    quad = np.polyfit(c, kappa, 2)
    pred_lin = np.polyval(lin, c)
    pred_quad = np.polyval(quad, c)
    ss_tot = np.sum((kappa - kappa.mean()) ** 2)
    r2_lin = 1 - np.sum((kappa - pred_lin) ** 2) / ss_tot
    r2_quad = 1 - np.sum((kappa - pred_quad) ** 2) / ss_tot
    max_dev = np.max(np.abs(kappa - pred_lin))
    max_dev_pct_range = max_dev / (kappa.max() - kappa.min()) * 100
    return dict(slope=lin[0], intercept=lin[1], r2_lin=r2_lin,
                curvature=quad[0], r2_quad=r2_quad,
                max_dev=max_dev, max_dev_pct_range=max_dev_pct_range)


def task2_linearity(curves):
    print("\n" + "=" * 100)
    print("TASK 2: linearity test (298.15 K curves, over the DATA3 range)")
    print("=" * 100)
    rows = []
    for salt, d in curves.items():
        for model_name, kappa in [("MSA", d["msa_298"]), ("Shedlovsky", d["shed_298"])]:
            stats = linearity_stats(d["c"], kappa)
            rows.append((salt, model_name, stats))
            print(f"\n{salt} / {model_name}:")
            print(f"  linear:    kappa = {stats['intercept']:+.4g} + {stats['slope']:.4g}*c  "
                  f"(R^2={stats['r2_lin']:.6f})")
            print(f"  quadratic curvature coeff: {stats['curvature']:.4e}  (R^2={stats['r2_quad']:.6f})")
            print(f"  max |linear fit - model|: {stats['max_dev']:.4g} uS/cm "
                  f"({stats['max_dev_pct_range']:.3f}% of range)")
    return rows


# ---------------------------------------------------------------------------
# Task 3: adjudicate the calibration discrepancy (discussion point 9)
# ---------------------------------------------------------------------------
def task3_adjudicate(curves):
    print("\n" + "=" * 100)
    print("TASK 3: calibration adjudication (discussion point 9)")
    print("=" * 100)
    results = {}
    for salt, d in curves.items():
        c = d["c"]
        model_lin = linearity_stats(c, d["msa_298"])
        model_slope, model_intercept = model_lin["slope"], model_lin["intercept"]

        inline = AUTHORS_INLINE[salt]
        vial = VIAL_ICP_FIT[salt]
        rel_slope_inline = (inline["s"] - model_slope) / model_slope * 100
        rel_slope_vial = (vial["s"] - model_slope) / model_slope * 100

        print(f"\n{salt}:")
        print(f"  MSA linear-fit calibration over DATA3 range: "
              f"kappa = {model_intercept:+.3f} + {model_slope:.3f}*c  [uS/cm, mM]")
        print(f"  Authors' inline: kappa = {inline['a']:+.3f} + {inline['s']:.3f}*c  "
              f"(slope rel. to MSA: {rel_slope_inline:+.1f}%, intercept: {inline['a']:+.2f} uS/cm)")
        print(f"  Vial-ICP fit:    kappa = {vial['a']:+.3f} + {vial['s']:.3f}*c  "
              f"(slope rel. to MSA: {rel_slope_vial:+.1f}%, intercept: {vial['a']:+.2f} uS/cm)")
        print(f"  Physically-admissible intercept from MSA (small-c limit): "
              f"{model_intercept:+.3f} uS/cm ({'small positive' if model_intercept >= 0 else 'small negative'})")
        results[salt] = dict(model_slope=model_slope, model_intercept=model_intercept,
                              rel_slope_inline=rel_slope_inline, rel_slope_vial=rel_slope_vial)
    return results


# ---------------------------------------------------------------------------
# Task 4: reproduce manuscript Table 3 MAPEs (implementation check)
# ---------------------------------------------------------------------------
MANUSCRIPT_DATA_FILES = {
    "NaCl": ("NaCl.csv", "NaCl conc", "NaCl cond cal"),
    "CaCl2": ("CaCl2.csv", "CaCl2_conc", "CaCl2_cond_cal"),
    "LaCl3": ("LaCl3.csv", "LaCl3 conc", "LaCl3 cond cal"),
}


def mape(pred, obs):
    pred = np.asarray(pred, dtype=float)
    obs = np.asarray(obs, dtype=float)
    return float(np.mean(np.abs((pred - obs) / obs)) * 100)


def task4_reproduce_manuscript():
    print("\n" + "=" * 100)
    print("TASK 4: reproduce manuscript Table 3 MAPEs (implementation check)")
    print("=" * 100)
    targets = {"NaCl": (1.0, 1.0), "CaCl2": (2.2, 2.3), "LaCl3": (1.9, 2.2)}
    results = {}
    for salt, (fname, conc_col, cond_col) in MANUSCRIPT_DATA_FILES.items():
        df = pd.read_csv(MANUSCRIPT_DATA_DIR / fname)
        conc_mM = df[conc_col].to_numpy(dtype=float)
        cond_cal_mS_cm = df[cond_col].to_numpy(dtype=float) / 1000.0  # uS/cm -> mS/cm

        shed_uS, _ = kappa_shedlovsky_uS_cm(conc_mM, salt, T_MANUSCRIPT)
        shed_mS = shed_uS / 1000.0
        msa_mS = kappa_msa_uS_cm(conc_mM, salt, T_MANUSCRIPT) / 1000.0

        mape_shed = mape(shed_mS, cond_cal_mS_cm)
        mape_msa = mape(msa_mS, cond_cal_mS_cm)
        tgt_shed, tgt_msa = targets[salt]
        print(f"\n{salt} (n={len(df)}): Shedlovsky MAPE = {mape_shed:.2f}% "
              f"(manuscript: {tgt_shed:.1f}%), MSA MAPE = {mape_msa:.2f}% "
              f"(manuscript: {tgt_msa:.1f}%)")
        results[salt] = dict(mape_shed=mape_shed, mape_msa=mape_msa,
                              target_shed=tgt_shed, target_msa=tgt_msa)
    return results


# ---------------------------------------------------------------------------
# Task 5: figure (plot_style compliant)
# ---------------------------------------------------------------------------
def task5_figure(curves):
    print("\n" + "=" * 100)
    print("TASK 5: figure")
    print("=" * 100)
    plot_style.apply_style()
    salts = ["NaCl", "CaCl2", "LaCl3"]
    colors = plot_style.OKABE_ITO
    fig = plt.figure(figsize=(plot_style.fig_width("full"), plot_style.fig_width("full") * 0.62))
    gs = GridSpec(2, 3, height_ratios=[1, 2.2], hspace=0.08, wspace=0.32, figure=fig)

    for j, salt in enumerate(salts):
        d = curves[salt]
        c = d["c"]
        msa = d["msa_298"]
        shed = d["shed_298"]
        extrap = d["extrap_298"]
        lin = linearity_stats(c, msa)
        pred_lin = np.polyval([lin["slope"], lin["intercept"]], c)

        ax_res = fig.add_subplot(gs[0, j])
        ax_main = fig.add_subplot(gs[1, j], sharex=ax_res)

        # top: curvature panel (model - best linear fit)
        ax_res.axhline(0, color="0.6", lw=1)
        ax_res.plot(c, msa - pred_lin, color=colors["black"], lw=1.5)
        ax_res.set_title(salt)
        ax_res.tick_params(labelbottom=False)
        if j == 0:
            ax_res.set_ylabel(r"MSA $-$ lin.$\,$fit" + "\n[$\\mu$S/cm]", fontsize=8)

        # bottom: kappa(c) curves + calibration overlays
        ax_main.plot(c, msa, "-", color=colors["orange"], lw=2, label="MSA (298.15 K)")
        c_valid = np.where(~extrap, c, np.nan)
        c_extrap_only = np.where(extrap, c, np.nan)
        ax_main.plot(c_valid, np.where(~extrap, shed, np.nan), "--", color=colors["black"],
                     lw=1.5, label="Shedlovsky (valid, $\\leq0.1$N)")
        ax_main.plot(c_extrap_only, np.where(extrap, shed, np.nan), ":", color=colors["black"],
                     lw=1.5, label="Shedlovsky (extrapolated)")
        inline = AUTHORS_INLINE[salt]
        vial = VIAL_ICP_FIT[salt]
        ax_main.plot(c, inline["a"] + inline["s"] * c, "-.", color=colors["vermillion"],
                     lw=1.5, label="Processed-sheet inline")
        ax_main.plot(c, vial["a"] + vial["s"] * c, "-.", color=colors["blue"],
                     lw=1.5, label="Vial-ICP fit")

        ax_main.set_xlabel(r"$c$ [mM]")
        if j == 0:
            ax_main.set_ylabel(r"$\boldsymbol{\kappa}$ [$\mu$S/cm]")
        ax_main.tick_params(direction="in", top=True, right=True)
        if j == 2:
            ax_main.legend(loc="upper left", fontsize=7)

    fig.suptitle("First-principles conductivity models vs. DATA3 calibrations", fontsize=11, fontweight="bold")
    out_path = FIG_DIR / "conductivity_model_calibration.png"
    plot_style.save_fig(fig, out_path)
    print(f"Saved {out_path} (+ .pdf)")
    plt.close(fig)


if __name__ == "__main__":
    curves = task1_generate_curves()
    task2_linearity(curves)
    task3_adjudicate(curves)
    task4_reproduce_manuscript()
    task5_figure(curves)
