"""
Step 2: Port of Bill's MATLAB single-salt (NaCl) nonlinear regression to
Python/scipy, reproducing prior_analysis/B_regression_concentration_NAClanal.m.

Model background
-----------------
Donnan equilibrium gives the membrane-phase co-ion (Cl-) concentration as a
function of the adjacent solution-phase concentration c_s:

    c_co^m(c_s) = ( -|chi| + sqrt(|chi|^2 + 4*delta_star*c_s^2) ) / 2

Bill's preliminary MATLAB fits

    f(beta; c_int, c_p) = sqrt(beta1^2 + 4*beta2*c_int^2)
                         - sqrt(beta1^2 + 4*beta2*c_p^2)

to a transformed "response" deltaC_m, where beta1 <-> |chi| and
beta2 <-> delta_star. Note this MATLAB form omits the factor of 1/2 from the
Donnan root above, so it effectively fits 2*(c_co^m(c_int) - c_co^m(c_p))
rather than the physical membrane-phase concentration difference. Bill says
the script was preliminary, so we treat the missing 1/2 as a likely bug.

This script:
  1. Reproduces Bill's regression EXACTLY (same rows, same uncorrected model
     "f") for a clean numerical comparison.
  2. Repeats the same uncorrected-model fit using the full valid data range
     found in the spreadsheet (rather than the truncated MATLAB row window),
     to see how much the truncation (skipping the low-concentration startup
     transient) affects the estimates.
  3. Fits the corrected model f_corrected = 0.5 * f to the SAME deltaC_m data,
     both for the MATLAB row window and the full valid range, to quantify how
     |chi| and delta_star shift once the factor of 1/2 is restored.
  4. Produces the two requested figures and prints/saves a summary table.

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/step2_nacl/regress_nacl.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

# ---------------------------------------------------------------------------
# Paths (all relative to the repo root, so the script can be run from anywhere)
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
import plot_style  # noqa: E402
import cp_authors_bulk  # noqa: E402

plot_style.apply_style()

DATA_FILE = REPO_ROOT / "data" / "Rejection_Analysis.xlsx"
SHEET_NAME = "NF270_MC5 05.27.26_NaCl (2)"
FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# CP Phase 2 (docs/prompts/cp_phase2_regressions.md): apply concentration
# polarization on top of the authors' bulk c_int (the "(2)" sheet's column
# B) as the new standing convention. Default primary run is CP=on;
# APPLY_CP=False reproduces the pre-Phase-2 (bulk == interfacial) numbers
# byte-for-byte (regression-tested).
APPLY_CP = True

# ---------------------------------------------------------------------------
# Physical constants (from the MATLAB script)
# ---------------------------------------------------------------------------
D_NA_M = 1.33e-12  # membrane-phase Na+ diffusivity [m^2/s]
D_CL_M = 2.03e-12  # membrane-phase Cl- diffusivity [m^2/s]
L_MEMBRANE = 80e-9  # effective membrane thickness [m]
# Ambipolar (mixed) NaCl membrane-phase diffusivity [m^2/s]
D_NACL_M = 2 * D_NA_M * D_CL_M / (D_NA_M + D_CL_M)

BETA0 = np.array([22.73, 0.09])  # initial guess: (|chi| [mM], delta_star [-])

# 1-indexed spreadsheet columns used, matching the MATLAB readmatrix ranges
COL_C_INT = "B"  # NaCl Interfacial Concentration (mM)
COL_C_P = "H"    # NaCl permeate interfacial concentration (mM) (smoothed)
COL_JW = "J"     # Jw (m3/m2*s)

# MATLAB used the fixed row window B57:B748 (1-indexed spreadsheet rows,
# header on row 1). This skips the low-concentration startup transient.
MATLAB_FIRST_ROW = 57
MATLAB_LAST_ROW = 748


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_sheet(sheet_name: str = SHEET_NAME, already_cp: bool = False,
               apply_cp: bool | None = None) -> pd.DataFrame:
    """Load a NaCl sheet (defaults to the Step-2 primary sheet), 0-indexed
    DataFrame rows correspond to spreadsheet row (index + 2) since row 1 is
    the header.

    CP Phase 2: if `apply_cp` (default: the module's `APPLY_CP` flag) and
    the sheet is NOT `already_cp` (i.e. its column B is the authors' bulk,
    non-CP convention -- true for the "(2)" sheet, false for the main
    05.27.26/05.07.24 sheets which the authors already CP-corrected), apply
    concentration polarization on top of the loaded bulk c_int
    (cp_authors_bulk.apply_cp, using the Phase-1 k_salt). c_p is returned
    unchanged either way (thin-film model: permeate unpolarized)."""
    if apply_cp is None:
        apply_cp = APPLY_CP
    df = pd.read_excel(DATA_FILE, sheet_name=sheet_name, header=0)
    # pandas column order matches spreadsheet column order (A, B, C, ...).
    # Column A is the leftmost -> index 0; B -> index 1; H -> index 7; J -> index 9.
    col_map = {
        "c_int": df.columns[1],  # column B
        "c_p": df.columns[7],    # column H
        "J_w": df.columns[9],    # column J
    }
    out = df[[col_map["c_int"], col_map["c_p"], col_map["J_w"]]].copy()
    out.columns = ["c_int", "c_p", "J_w"]
    # spreadsheet row number for each DataFrame row (header is row 1, first
    # data row is row 2)
    out.index = np.arange(2, 2 + len(out))
    out.index.name = "sheet_row"
    if apply_cp and not already_cp:
        out["c_int"] = cp_authors_bulk.apply_cp(out["c_int"], out["c_p"], out["J_w"], "NaCl")
    return out


def report_data_extent(df_full: pd.DataFrame) -> pd.DataFrame:
    """Report where valid (non-NaN) data actually lives, and any anomalies."""
    valid_mask = df_full.notna().all(axis=1)
    valid_rows = df_full.index[valid_mask]
    first_valid, last_valid = valid_rows.min(), valid_rows.max()

    n_total = len(df_full)
    n_valid = valid_mask.sum()
    n_nan = n_total - n_valid

    df_valid = df_full.loc[valid_mask]
    n_zero_Jw = int((df_valid["J_w"] == 0).sum())
    n_neg_Jw = int((df_valid["J_w"] < 0).sum())
    n_beyond_748 = int((valid_rows > MATLAB_LAST_ROW).sum())

    print("=" * 70)
    print("DATA EXTENT / ANOMALY CHECK")
    print("=" * 70)
    print(f"Sheet dimension reports rows up to {df_full.index.max()} "
          f"(workbook 'dimension' can overstate actual data extent).")
    print(f"Valid (non-NaN across c_int, c_p, J_w) data spans rows "
          f"{first_valid}-{last_valid}  (n = {n_valid}).")
    print(f"Rows with any NaN in this window: {n_nan}")
    print(f"Rows with J_w == 0: {n_zero_Jw}; J_w < 0: {n_neg_Jw}")
    print(f"Valid rows beyond {MATLAB_LAST_ROW} "
          f"(MATLAB's last row): {n_beyond_748}")
    print(f"MATLAB used rows {MATLAB_FIRST_ROW}-{MATLAB_LAST_ROW} "
          f"(n = {MATLAB_LAST_ROW - MATLAB_FIRST_ROW + 1}), skipping the "
          f"low-concentration startup transient in rows "
          f"{first_valid}-{MATLAB_FIRST_ROW - 1}.")
    print()
    return df_valid


def build_deltaC_m(df: pd.DataFrame) -> pd.DataFrame:
    """Add the transformed response deltaC_m [mM] to the dataframe.

    J_s = c_p * J_w            [mol/m^2/s]  (c_p in mM = mol/m^3, J_w in m/s)
    deltaC_m = J_s * l / D_NaCl_m            [mM]
    """
    df = df.copy()
    J_s = df["c_p"] * df["J_w"]  # mol/m^2/s
    df["J_s"] = J_s
    df["deltaC_m"] = J_s * L_MEMBRANE / D_NACL_M  # mM
    return df


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
def donnan_root(beta, c):
    """sqrt(beta1^2 + 4*beta2*c^2), beta = (|chi|, delta_star)."""
    beta1, beta2 = beta
    return np.sqrt(beta1 ** 2 + 4.0 * beta2 * c ** 2)


def f_matlab(beta, c_int, c_p):
    """Bill's (uncorrected) MATLAB model -- omits the 1/2 factor."""
    return donnan_root(beta, c_int) - donnan_root(beta, c_p)


def f_corrected(beta, c_int, c_p):
    """Physically correct model: 0.5 * f_matlab, i.e. the true Donnan
    co-ion concentration difference c_co^m(c_int) - c_co^m(c_p).

    Note: f_matlab(beta1, beta2; c) is exactly homogeneous under
    (beta1, beta2) -> (2*beta1, 4*beta2): f_matlab(2*b1, 4*b2; c) =
    2*f_matlab(b1, b2; c) for all c (pull a factor of 4 out of the sqrt,
    which becomes a factor of 2). Consequently, on a fixed data set,
    least-squares fits of f_corrected = 0.5*f_matlab and f_matlab are
    related by an EXACT rescaling of the optimum:
        |chi|_corrected      = 2 * |chi|_uncorrected
        delta_star_corrected = 4 * delta_star_uncorrected
    with identical SSE at the optimum (verified numerically below to
    ~1e-10 relative precision). This is a stronger, exact version of the
    qualitative expectation "chi roughly unchanged, delta_star rescales"
    -- for this exactly homogeneous model, chi in fact doubles, and it is
    delta_star that would be "roughly unchanged" only in a loose order-of
    magnitude sense (it actually scales by a clean factor of 4).
    """
    return 0.5 * f_matlab(beta, c_int, c_p)


# ---------------------------------------------------------------------------
# Fitting utilities
# ---------------------------------------------------------------------------
@dataclass
class FitResult:
    label: str
    beta: np.ndarray
    se: np.ndarray
    sse: float
    n: int
    p: int
    cov: np.ndarray


def fit_model(model_fn, c_int, c_p, y, beta0=BETA0, label="fit") -> FitResult:
    """Nonlinear least squares fit using scipy.optimize.least_squares.

    Standard errors are computed from the linearized covariance matrix
    cov = sigma^2 * (J^T J)^-1, with sigma^2 = SSE / (n - p) (as scipy's
    curve_fit does internally), where J is the Jacobian of the residuals
    at the optimum.
    """
    def resid(beta):
        return model_fn(beta, c_int, c_p) - y

    result = least_squares(resid, beta0, method="lm")  # Levenberg-Marquardt,
    # matches MATLAB lsqcurvefit's default trust-region-reflective closely
    # enough for this well-conditioned, unconstrained problem.

    beta_hat = result.x
    residuals = result.fun
    sse = float(np.sum(residuals ** 2))
    n = len(y)
    p = len(beta_hat)
    dof = max(n - p, 1)

    J = result.jac  # Jacobian of residuals w.r.t. beta, shape (n, p)
    JTJ = J.T @ J
    try:
        JTJ_inv = np.linalg.inv(JTJ)
    except np.linalg.LinAlgError:
        JTJ_inv = np.linalg.pinv(JTJ)

    sigma2 = sse / dof
    cov = sigma2 * JTJ_inv
    se = np.sqrt(np.diag(cov))

    return FitResult(label=label, beta=beta_hat, se=se, sse=sse, n=n, p=p, cov=cov)


def print_fit(fr: FitResult):
    print(f"--- {fr.label} ---")
    print(f"  n = {fr.n} data points, p = {fr.p} parameters")
    print(f"  |chi|       = {fr.beta[0]:.6f}  +/- {fr.se[0]:.6f}  mM")
    print(f"  delta_star  = {fr.beta[1]:.6f}  +/- {fr.se[1]:.6f}  (-)")
    print(f"  SSE         = {fr.sse:.6f}")
    print()


# ---------------------------------------------------------------------------
# Residual contour (mirrors MATLAB figure(7): contourf(chi, log10(delta), log10(SSE)))
# ---------------------------------------------------------------------------
def compute_residual_grid(model_fn, c_int, c_p, y,
                           chi_range=None, del_range=None):
    if chi_range is None:
        chi_range = np.linspace(20, 30, 110)
    if del_range is None:
        del_range = np.logspace(-1.6, -1, 110)

    sse_grid = np.zeros((len(chi_range), len(del_range)))
    for i, chi in enumerate(chi_range):
        for j, delta in enumerate(del_range):
            pred = model_fn((chi, delta), c_int, c_p)
            sse_grid[i, j] = np.sum((pred - y) ** 2)
    return chi_range, del_range, sse_grid


def plot_residual_contour(chi_range, del_range, sse_grid, beta_opt, outfile,
                           title):
    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 5.5 / 7))
    # sse_grid indexed [i (chi), j (delta)] -> transpose for (delta, chi) axes
    log_sse = np.log10(sse_grid.T)
    cf = ax.contourf(chi_range, np.log10(del_range), log_sse, levels=30,
                      cmap="viridis")
    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label(r"$\boldsymbol{\log_{10}(\mathrm{SSE})}$")
    ax.plot(beta_opt[0], np.log10(beta_opt[1]), marker="*", color="red",
            markersize=18, markeredgecolor="white", linestyle="none",
            label="fitted optimum")
    ax.set_xlabel(r"$\boldsymbol{|\chi|}$ [mM]")
    ax.set_ylabel(r"$\boldsymbol{\log_{10}(\delta^*)}$")
    ax.set_title(title)
    ax.legend(loc="best")
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_fit(df, fr_matlab, outfile, title):
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()
    y = df["deltaC_m"].to_numpy()

    order = np.argsort(c_int)
    c_int_sorted = c_int[order]
    c_p_sorted = c_p[order]

    y_fit = f_matlab(fr_matlab.beta, c_int_sorted, c_p_sorted)

    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 5.5 / 7))
    ax.plot(c_int, y, "o", color="tab:red", markersize=3, alpha=0.6,
            label=r"data: $\Delta c_m$")
    ax.plot(c_int_sorted, y_fit, "--", color="tab:red", linewidth=2,
            label=r"fitted model $f(\beta;\, c_{int}, c_p)$")
    ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"$\Delta \boldsymbol{c_m}$ [mM]")
    ax.set_title(title)
    ax.legend(loc="best")
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("Step 2: NaCl single-salt Donnan regression (Python port)")
    print("=" * 70)
    print(f"Data file: {DATA_FILE}")
    print(f"Sheet:     {SHEET_NAME}")
    print(f"D_NaCl_m = {D_NACL_M:.6e} m^2/s   l = {L_MEMBRANE:.3e} m")
    print(f"APPLY_CP = {APPLY_CP}  (CP Phase 2: {'CP-corrected c_int (primary)' if APPLY_CP else 'bulk c_int, non-CP (regression-test/comparison mode)'})")
    print()

    df_full = load_sheet()
    df_valid = report_data_extent(df_full)

    # ------------------------------------------------------------------
    # (A) Reproduce MATLAB exactly: rows 57-748, uncorrected model f_matlab
    # ------------------------------------------------------------------
    df_matlab_window = df_full.loc[MATLAB_FIRST_ROW:MATLAB_LAST_ROW].copy()
    assert df_matlab_window.notna().all().all(), (
        "Unexpected NaNs inside the MATLAB row window 57:748"
    )
    df_matlab_window = build_deltaC_m(df_matlab_window)

    print(f"Sanity check units: deltaC_m = J_s * l / D_NaCl_m has units "
          f"[mol/m^2/s * m] / [m^2/s] = [mol/m^3] = [mM] (since c_p is in "
          f"mM = mol/m^3). Example deltaC_m values (mM): "
          f"{df_matlab_window['deltaC_m'].iloc[:3].to_numpy()}")
    print()

    fr_bill = fit_model(
        f_matlab,
        df_matlab_window["c_int"].to_numpy(),
        df_matlab_window["c_p"].to_numpy(),
        df_matlab_window["deltaC_m"].to_numpy(),
        beta0=BETA0,
        label="Bill's MATLAB reproduction (rows 57-748, uncorrected f)",
    )
    print_fit(fr_bill)

    # ------------------------------------------------------------------
    # (B) Same uncorrected model, but using the FULL valid data range
    # ------------------------------------------------------------------
    df_full_valid = build_deltaC_m(df_valid)
    fr_full_uncorrected = fit_model(
        f_matlab,
        df_full_valid["c_int"].to_numpy(),
        df_full_valid["c_p"].to_numpy(),
        df_full_valid["deltaC_m"].to_numpy(),
        beta0=BETA0,
        label=f"Full valid range (rows {df_valid.index.min()}-"
              f"{df_valid.index.max()}, uncorrected f)",
    )
    print_fit(fr_full_uncorrected)

    # ------------------------------------------------------------------
    # (C) Corrected model (factor of 1/2 restored), both row windows
    # ------------------------------------------------------------------
    fr_bill_corrected = fit_model(
        f_corrected,
        df_matlab_window["c_int"].to_numpy(),
        df_matlab_window["c_p"].to_numpy(),
        df_matlab_window["deltaC_m"].to_numpy(),
        beta0=BETA0,
        label="Corrected model (rows 57-748, f_corrected = 0.5*f)",
    )
    print_fit(fr_bill_corrected)

    fr_full_corrected = fit_model(
        f_corrected,
        df_full_valid["c_int"].to_numpy(),
        df_full_valid["c_p"].to_numpy(),
        df_full_valid["deltaC_m"].to_numpy(),
        beta0=BETA0,
        label=f"Corrected model (rows {df_valid.index.min()}-"
              f"{df_valid.index.max()}, f_corrected = 0.5*f)",
    )
    print_fit(fr_full_corrected)

    # ------------------------------------------------------------------
    # Quantify the expectation: |chi| ~ unchanged, delta_star rescales
    # ------------------------------------------------------------------
    chi_ratio = fr_bill_corrected.beta[0] / fr_bill.beta[0]
    delta_ratio = fr_bill_corrected.beta[1] / fr_bill.beta[1]
    print("=" * 70)
    print("Effect of restoring the factor of 1/2 (rows 57-748):")
    print(f"  |chi|_corrected / |chi|_uncorrected       = {chi_ratio:.6f}  (exact identity: 2)")
    print(f"  delta_star_corrected / delta_star_uncorrected = {delta_ratio:.6f}  (exact identity: 4)")
    print(f"  SSE unchanged: {fr_bill.sse:.6f} vs {fr_bill_corrected.sse:.6f}")
    print("  f_matlab(beta1,beta2;c) is exactly homogeneous under "
          "(beta1,beta2) -> (2*beta1,4*beta2), which scales f_matlab by "
          "exactly 2. So fitting 0.5*f_matlab to the same data reproduces "
          "the identical SSE surface with chi doubled and delta_star "
          "quadrupled -- not merely 'roughly unchanged' chi, but an exact "
          "factor of 2.")
    print()

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------
    plot_fit(
        df_matlab_window, fr_bill,
        FIG_DIR / "nacl_fit.png",
        "NaCl Donnan regression: data vs. fitted model\n"
        "(Bill's MATLAB reproduction, rows 57-748)",
    )

    chi_range, del_range, sse_grid = compute_residual_grid(
        f_matlab,
        df_matlab_window["c_int"].to_numpy(),
        df_matlab_window["c_p"].to_numpy(),
        df_matlab_window["deltaC_m"].to_numpy(),
    )
    plot_residual_contour(
        chi_range, del_range, sse_grid, fr_bill.beta,
        FIG_DIR / "nacl_residual_contour.png",
        r"$\boldsymbol{\log_{10}(\mathrm{SSE})}$ contour (rows 57-748, uncorrected $\boldsymbol{f}$)",
    )

    # ------------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------------
    print("=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    header = (f"{'Fit':<55}{'|chi| (mM)':>14}{'SE(chi)':>10}"
              f"{'delta*':>12}{'SE(delta*)':>12}{'SSE':>14}{'n':>6}")
    print(header)
    for fr in [fr_bill, fr_full_uncorrected, fr_bill_corrected, fr_full_corrected]:
        print(f"{fr.label:<55}{fr.beta[0]:>14.4f}{fr.se[0]:>10.4f}"
              f"{fr.beta[1]:>12.6f}{fr.se[1]:>12.6f}{fr.sse:>14.4f}{fr.n:>6}")
    print()

    return {
        "bill": fr_bill,
        "full_uncorrected": fr_full_uncorrected,
        "bill_corrected": fr_bill_corrected,
        "full_corrected": fr_full_corrected,
        "df_valid_extent": (df_valid.index.min(), df_valid.index.max()),
    }


if __name__ == "__main__":
    main()
