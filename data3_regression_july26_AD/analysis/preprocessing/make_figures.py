"""
Illustrative figures for the preprocessing pipeline (docs/reports/main.tex,
sec:preproc): a representative per-experiment calibration fit, and a
reconstructed-vs-reference parity plot against the "(2)" sheet.

Run with: conda activate data3-regression && python analysis/preprocessing/make_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import preprocess_nacl as pp
import plot_style
from validate import load_processed_sheet, GIVEN_CALIB

plot_style.apply_style()

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def plot_calibration_example(outfile):
    """Vial (conductivity, dilution-corrected ICP) pairs and the fitted
    linear calibration, for a representative raw-only run."""
    raw = pp.load_raw_sheet(DATA_DIR / "NF270_MC5.xlsx", "07.23.24_NaCl")
    calib = pp.fit_calibration(raw.vials)

    df = raw.vials.dropna(subset=["cond", "icp_mgL", "sample_vol_mL", "acid_vol_mL"])
    dilution = (df["sample_vol_mL"] + df["acid_vol_mL"]) / df["sample_vol_mL"]
    c_mM = (df["icp_mgL"].to_numpy() * dilution.to_numpy()) / pp.ICP_MOLAR_MASS
    cond = df["cond"].to_numpy()

    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 1.1))
    ax.scatter(c_mM, cond, color=plot_style.OKABE_ITO["blue"], s=40,
               zorder=3, label="vial (ICP, conductivity)")
    c_grid = np.linspace(0, c_mM.max() * 1.05, 50)
    ax.plot(c_grid, calib.intercept + calib.slope * c_grid, "k--", lw=2,
            label=f"fit: cond $=$ {calib.intercept:.1f} $+$ {calib.slope:.1f}$\\cdot c$")
    ax.set_xlabel(r"$\boldsymbol{c}$ [mM] (dilution-corrected ICP)")
    ax.set_ylabel(r"conductivity [$\boldsymbol{\mu}$S/cm]")
    ax.set_title(f"Per-experiment calibration (MC5 07.23.24_NaCl)\n"
                 f"$\\boldsymbol{{R^2}}$={calib.r2:.4f}, n={calib.n_points}", fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1,
              fontsize=8, borderaxespad=0)
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_reconstruction_parity(outfile):
    """Parity plot: our pipeline's c_int (non-CP, embedded calibration) vs.
    the "(2)" sheet's B column, split into the 22 live-formula rows (exact)
    and the 725 hardcoded/pasted rows (approximate)."""
    raw = pp.load_raw_sheet(DATA_DIR / "NF270_MC5.xlsx", "05.27.26_NaCl")
    n = len(raw.ts)
    ref = load_processed_sheet(DATA_DIR / "BoE Analysis.xlsx",
                                "NF270_MC5 05.27.26_NaCl (2)", n)
    result = pp.process_experiment(raw, apply_cp=False, jw_window=51,
                                    calibration=GIVEN_CALIB)
    ours = result.df["c_int_mM"].to_numpy()
    theirs = ref["B"].to_numpy()

    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 1.35))
    ax.scatter(theirs[22:], ours[22:], s=10, alpha=0.5,
               color=plot_style.OKABE_ITO["vermillion"],
               label="rows 24-748 (hardcoded/pasted)")
    ax.scatter(theirs[:22], ours[:22], s=30,
               color=plot_style.OKABE_ITO["bluish_green"],
               label="rows 2-23 (live formula)")
    lims = [0, max(theirs.max(), ours.max()) * 1.05]
    ax.plot(lims, lims, "k:", lw=1.5, label="parity ($y=x$)")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel(r"reference $\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"reconstructed $\boldsymbol{c_{int}}$ [mM]")
    ax.set_title("Reconstruction vs. reference\n(05.27.26_NaCl, non-CP)", fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1,
              fontsize=8, borderaxespad=0)
    ax.set_box_aspect(1)
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    plot_calibration_example(FIG_DIR / "nacl_preproc_calibration.png")
    plot_reconstruction_parity(FIG_DIR / "nacl_preproc_reconstruction.png")
