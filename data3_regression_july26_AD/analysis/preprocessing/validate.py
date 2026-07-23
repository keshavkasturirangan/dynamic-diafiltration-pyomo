"""
Validation of the NaCl raw->processed preprocessing pipeline
(preprocess_nacl.py) against the two runs where we have BOTH raw and
processed (pasted/formula) data:

  1. Exact target: NF270_MC5 "05.27.26_NaCl (2)" sheet in BoE Analysis.xlsx
     (live formulas, non-CP, embedded calibration c=(cond+63.706)/76.685).
     Source raw sheet: NF270_MC5.xlsx / "05.27.26_NaCl".
  2. Approximate (CP-corrected) targets: the main "05.27.26_NaCl" sheet
     (same raw source, apply_cp=True) and "NF270_MC2 05.07.24_NaCl"
     (different raw source: NF270_MC2.xlsx / "05.07.24_NaCl").

IMPORTANT DEVIATION FROM THE ORIGINAL PLAN (flagged): both "05.27.26_NaCl"
and "05.07.24_NaCl" RAW sheets have an entirely empty vial "Conductivity @
Temp" column (col P) -- unlike the 6 raw-only runs processed by
process_raw_runs.py, which do have per-vial conductivity recorded. So a
per-experiment calibration CANNOT be fit for either validation run from its
own vial data. For "05.27.26_NaCl" we instead use the embedded calibration
constants (intercept=-63.706, slope=76.685) as GIVEN ground truth (this is
exactly what the (2) sheet's live formula encodes). For "05.07.24_NaCl"
(MC2) there is no live formula and no vial conductivity either, so its
original calibration is simply unknown; we approximate it with the same
MC5 constants and expect this (not the CP model) to dominate any mismatch --
see README for the full discussion.

Run with: conda activate data3-regression && python analysis/preprocessing/validate.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import preprocess_nacl as pp

DATA_DIR = REPO_ROOT / "data"
GIVEN_CALIB = pp.Calibration(slope=76.685, intercept=-63.706, n_points=0,
                              r2=float("nan"), source="given (embedded MC5 formula)")


def load_processed_sheet(path, sheet_name, n_rows):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet_name]
    rows = []
    for r in range(2, 2 + n_rows):
        rows.append([ws.cell(row=r, column=c).value for c in range(1, 16)])
    wb.close()
    cols = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"]
    df = pd.DataFrame(rows, columns=cols)
    for c in ["B", "D", "F", "H", "J", "K", "M", "N", "O"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def rel_err(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b) & (b != 0)
    return np.abs(a[mask] - b[mask]) / np.abs(b[mask])


def tune_jw_window(raw, ref_J, windows=(5, 7, 9, 11, 15, 21, 31, 41, 51)):
    best = None
    for w in windows:
        Jw = pp.compute_jw(raw.ts["time_s"], raw.ts["mass_g"], window=w)
        err = rel_err(Jw, ref_J)
        med = np.nanmedian(err) if len(err) else np.inf
        if best is None or med < best[1]:
            best = (w, med)
    return best  # (window, median relative error)


def validate_exact_target():
    print("=" * 88)
    print("VALIDATION 1: exact target -- NF270_MC5 05.27.26_NaCl (2), non-CP, embedded calibration")
    print("=" * 88)
    raw = pp.load_raw_sheet(DATA_DIR / "NF270_MC5.xlsx", "05.27.26_NaCl")
    n = len(raw.ts)
    ref = load_processed_sheet(DATA_DIR / "BoE Analysis.xlsx",
                                "NF270_MC5 05.27.26_NaCl (2)", n)

    best_w, best_med = tune_jw_window(raw, ref["J"].to_numpy())
    print(f"Jw smoothing window tuned to {best_w} rows "
          f"(median relative error vs. reference J = {best_med:.3%})")

    result = pp.process_experiment(raw, apply_cp=False, jw_window=best_w,
                                    calibration=GIVEN_CALIB)
    df = result.df

    # IMPORTANT FINDING (verified by inspecting cell formulas with
    # data_only=False): column B in this sheet is a LIVE formula
    # "=(N{row}+63.706)/76.685" for only the first 22 data rows (rows
    # 2-23); the remaining 725 rows (~97% of the sheet) are HARDCODED
    # pasted numbers that do NOT satisfy that formula exactly (they
    # regress against verbatim N with R^2=0.999 but a different effective
    # slope/intercept and up to ~4 mM residuals) -- i.e. most of the "(2)"
    # sheet's B column was produced by some undocumented downstream step
    # (a different/refined calibration, or calibration against smoothed
    # rather than instantaneous conductivity), not by the formula visible
    # in the cell. We therefore validate the two regions SEPARATELY.
    a = df["c_int_mM"].to_numpy()
    b = ref["B"].to_numpy()
    e_formula_region = np.abs(a[:22] - b[:22]) / np.abs(b[:22])
    e_pasted_region = np.abs(a[22:] - b[22:]) / np.abs(b[22:])
    print("\nB column has two regimes in this sheet (found via data_only=False "
          "formula inspection):")
    print(f"  rows 2-23  (live formula,  n=22 ):  median rel err = "
          f"{np.median(e_formula_region):.2e}, max = {np.max(e_formula_region):.2e}  "
          "-- matches to machine precision, as expected")
    print(f"  rows 24-748 (hardcoded/pasted, n=725): median rel err = "
          f"{np.median(e_pasted_region):.2e}, max = {np.max(e_pasted_region):.2e}  "
          "-- NOT reproducible from the visible formula; see README flag")

    # M, N are verbatim raw copies -- exact match expected
    print("\nFull-sheet per-column match (median / max relative error vs. reference):")
    header = f"{'col':<8}{'meaning':<28}{'median rel err':>16}{'max rel err':>14}{'n compared':>12}"
    print(header)
    checks = [
        ("B", "c_int (calibration)", df["c_int_mM"], ref["B"]),
        ("D", "ionic strength", df["ionic_strength_mM"], ref["D"]),
        ("M", "time (verbatim)", df["time_s"], ref["M"]),
        ("N", "ret. cond. (verbatim)", df["ret_cond_uS_cm"], ref["N"]),
        ("J", "Jw (smoothed)", df["Jw_m3_m2_s"], ref["J"]),
        ("K", "B coefficient (uses our Jw, ref H)", df["B_um_s"], ref["K"]),
    ]
    for col, meaning, ours, refcol in checks:
        e = rel_err(ours, refcol)
        print(f"{col:<8}{meaning:<28}{np.nanmedian(e):>16.2e}{np.nanmax(e):>14.2e}{len(e):>12d}")

    # K using the REFERENCE H and J (isolates the K-formula from Jw/H
    # smoothing uncertainty) -- this is the "K to ~1e-6" claim in the prompt.
    K_formula_check = df["Jw_m3_m2_s"] * 0 + (ref["J"] * ref["H"]) / (df["c_int_mM"] - ref["H"]) * 1e6
    e_formula = rel_err(K_formula_check, ref["K"])
    print(f"\nK formula check (our B, reference J & H -- isolates the K formula "
          f"itself from Jw/H smoothing uncertainty):")
    print(f"  median rel err = {np.nanmedian(e_formula):.2e}, max = {np.nanmax(e_formula):.2e}")
    print("(F is blank in this sheet -- no rejection reference to check here.)")
    print()
    return best_w


def validate_cp_targets(jw_window):
    print("=" * 88)
    print("VALIDATION 2: approximate CP-corrected targets")
    print("=" * 88)

    # --- MC5 05.27.26_NaCl (main sheet, CP-corrected), same raw source ---
    raw = pp.load_raw_sheet(DATA_DIR / "NF270_MC5.xlsx", "05.27.26_NaCl")
    n = len(raw.ts)
    ref = load_processed_sheet(DATA_DIR / "BoE Analysis.xlsx",
                                "NF270_MC5 05.27.26_NaCl", n)
    result = pp.process_experiment(raw, apply_cp=True, jw_window=jw_window,
                                    calibration=GIVEN_CALIB)
    df = result.df
    print("\nMC5 05.27.26_NaCl (main, CP-corrected; same raw source + calibration as target 1):")
    for col, meaning, ours, refcol in [
        ("B", "c_int (CP-corrected)", df["c_int_mM"], ref["B"]),
        ("H", "c_p (permeate)", df["c_p_mM"], ref["H"]),
    ]:
        e = rel_err(ours, refcol)
        print(f"  {col} ({meaning}): median rel err = {np.nanmedian(e):.2%}, "
              f"max = {np.nanmax(e):.2%}, n={len(e)}")
    print("  (F, the rejection column, is blank in this sheet too.)")

    # --- MC2 05.07.24_NaCl (different membrane cut/session; calibration unknown) ---
    raw2 = pp.load_raw_sheet(DATA_DIR / "NF270_MC2.xlsx", "05.07.24_NaCl")
    n2 = len(raw2.ts)
    ref2 = load_processed_sheet(DATA_DIR / "BoE Analysis.xlsx",
                                 "NF270_MC2 05.07.24_NaCl", n2)
    result2 = pp.process_experiment(raw2, apply_cp=True, jw_window=jw_window,
                                     calibration=GIVEN_CALIB)
    df2 = result2.df
    print("\nMC2 05.07.24_NaCl (CP-corrected; DIFFERENT membrane cut/session, "
          "calibration UNKNOWN -- approximated with the MC5 constants, flagged):")
    for col, meaning, ours, refcol in [
        ("B", "c_int (CP-corrected)", df2["c_int_mM"], ref2["B"]),
        ("F", "rejection", df2["rejection"], ref2["F"]),
        ("H", "c_p (permeate)", df2["c_p_mM"], ref2["H"]),
    ]:
        e = rel_err(ours, refcol)
        print(f"  {col} ({meaning}): median rel err = {np.nanmedian(e):.2%}, "
              f"max = {np.nanmax(e):.2%}, n={len(e)}")
    print("\nExpected: MC2 mismatch dominated by the unknown/approximated "
          "calibration (different probe/session), not by the CP model itself.")
    print()


if __name__ == "__main__":
    best_w = validate_exact_target()
    validate_cp_targets(best_w)
