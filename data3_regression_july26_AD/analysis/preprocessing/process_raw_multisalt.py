"""
Phase 1 of the CP redo (docs/prompts/cp_preprocessing_extend.md): apply the
salt-general preprocessing pipeline (preprocess_nacl.py's Salt-parametrized
process_experiment) to every CaCl2/LaCl3 dataset and the two NaCl
"validation twin" runs, with BOTH apply_cp=True and apply_cp=False, and
validate CP-off against the authors' non-CP processed sheets (and CP-on
against the one authors' sheet that is already CP-corrected).

Datasets (raw campaign workbooks in data/, all confirmed present):
  CaCl2  : MC3.xlsx::07.11.24_SCaCl2   (own vial calibration)
           MC5.xlsx::05.27.26_CaCl2    (embedded formula -- see below)
           MC5.xlsx::06.27.26_CaCl2    (NO vial cond, NO embedded formula --
                                        approximated, flagged; see below)
  LaCl3  : MC2.xlsx::05.21.24_LaCl3    (own vial calibration)
  NaCl   : MC5.xlsx::05.27.26_NaCl     (embedded formula, the validation twin)
           MC2.xlsx::05.07.24_NaCl     (calibration unknown, approximated
                                        with the MC5 embedded formula)

IMPORTANT DEVIATION FROM THE TASK PROMPT, FLAGGED: the prompt states "the
CaCl2 and LaCl3 raw sheets DO record [per-vial conductivity]". This is true
for MC3 07.11.24_SCaCl2 and MC2 05.21.24_LaCl3 (both have full 14-point vial
conductivity, used to fit their own calibration below), but FALSE for BOTH
MC5 CaCl2 sheets (05.27.26_CaCl2 and 06.27.26_CaCl2 have an entirely empty
vial "Conductivity @ Temp" column) -- the same situation as the two NaCl
validation-twin sheets. For 05.27.26_CaCl2 we found a live-formula fallback
(see below, analogous to the NaCl "(2)" sheet); for 06.27.26_CaCl2 there is
NO live formula anywhere in its processed sheet either (unlike 05.27.26's),
so its calibration is genuinely unrecoverable from this workbook -- we
approximate it with 05.27.26_CaCl2's calibration (same membrane cut, same
campaign, closest date) and flag the result as PROVISIONAL, exactly
mirroring how MC2 05.07.24_NaCl's calibration was already approximated with
the MC5 embedded formula in the original NaCl validation.

Embedded/given calibrations (found by inspecting cell *formulas*,
data_only=False, mirroring validate.py's finding for the NaCl "(2)" sheet):
  - NF270_MC5 05.27.26_NaCl (2): `=(N{row}+63.706)/76.685` for rows 2-23
    (intercept=-63.706, slope=76.685) -- pre-existing, from validate.py.
  - NF270_MC5 05.27.26_CaCl2:    `=N{row}/149.73` for rows 2-9 (intercept=0,
    slope=149.73) -- NEW finding, this task.
  - NF270_MC5 06.27.26_CaCl2:    no live formula found anywhere in the
    column -- NEW finding, this task; approximated as described above.

Run with: conda activate data3-regression && python analysis/preprocessing/process_raw_multisalt.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import preprocess_nacl as pp  # noqa: E402

DATA_DIR = REPO_ROOT / "data"
OUT_DIR = Path(__file__).resolve().parent / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

JW_WINDOW = 51  # tuned in validate.py against the one reference J column; reused throughout

GIVEN_CALIB_NACL = pp.Calibration(slope=76.685, intercept=-63.706, n_points=0,
                                   r2=float("nan"), source="given (embedded MC5 05.27.26_NaCl formula)")
GIVEN_CALIB_CACL2_MC5 = pp.Calibration(slope=149.73, intercept=0.0, n_points=0,
                                        r2=float("nan"), source="given (embedded MC5 05.27.26_CaCl2 formula)")

# Each entry: raw (file, sheet), salt, calibration (None => fit from own
# vial data), target processed sheet for validation, and whether the target
# sheet is CP-corrected (only the two "main" NaCl sheets are) or non-CP.
DATASETS = [
    dict(key="CaCl2_MC3-0711S", file="NF270_MC3.xlsx", sheet="07.11.24_SCaCl2",
         salt=pp.SALT_CACL2, calibration=None,
         target_sheet="NF270_MC3 07.11.24_SCaCl2", target_is_cp=False,
         out_prefix="CaCl2_MC3_07.11.24_SCaCl2"),
    dict(key="CaCl2_MC5-052726", file="NF270_MC5.xlsx", sheet="05.27.26_CaCl2",
         salt=pp.SALT_CACL2, calibration=GIVEN_CALIB_CACL2_MC5,
         target_sheet="NF270_MC5 05.27.26_CaCl2", target_is_cp=False,
         out_prefix="CaCl2_MC5_05.27.26_CaCl2"),
    dict(key="CaCl2_MC5-062726", file="NF270_MC5.xlsx", sheet="06.27.26_CaCl2",
         salt=pp.SALT_CACL2, calibration=GIVEN_CALIB_CACL2_MC5,
         target_sheet="NF270_MC5 06.27.26_CaCl2", target_is_cp=False,
         out_prefix="CaCl2_MC5_06.27.26_CaCl2", provisional_calib=True),
    dict(key="LaCl3_MC2-052124", file="NF270_MC2.xlsx", sheet="05.21.24_LaCl3",
         salt=pp.SALT_LACL3, calibration=None,
         target_sheet="NF270_MC2 05.21.24_LaCl3", target_is_cp=False,
         out_prefix="LaCl3_MC2_05.21.24_LaCl3"),
    dict(key="NaCl_MC5-052726", file="NF270_MC5.xlsx", sheet="05.27.26_NaCl",
         salt=pp.SALT_NACL, calibration=GIVEN_CALIB_NACL,
         target_sheet="NF270_MC5 05.27.26_NaCl (2)", target_is_cp=False,
         cp_target_sheet="NF270_MC5 05.27.26_NaCl",
         out_prefix="NaCl_MC5_05.27.26_NaCl"),
    dict(key="NaCl_MC2-050724", file="NF270_MC2.xlsx", sheet="05.07.24_NaCl",
         salt=pp.SALT_NACL, calibration=GIVEN_CALIB_NACL,
         target_sheet=None, target_is_cp=False,
         cp_target_sheet="NF270_MC2 05.07.24_NaCl",
         out_prefix="NaCl_MC2_05.07.24_NaCl", provisional_calib=True),
]


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


def pct_shift(nonCP, CP):
    nonCP = np.asarray(nonCP, dtype=float)
    CP = np.asarray(CP, dtype=float)
    mask = np.isfinite(nonCP) & np.isfinite(CP) & (nonCP != 0)
    shift = (CP[mask] - nonCP[mask]) / nonCP[mask] * 100.0
    return np.median(shift), np.max(shift)


def write_csv(out_path, df, salt, calib, apply_cp, source_file, source_sheet, provisional=False):
    k = pp.mass_transfer_coefficient(D=salt.solution_D)
    with open(out_path, "w") as f:
        f.write(f"# source: {source_file} / {source_sheet}\n")
        f.write(f"# salt: {salt.name}, ICP molar mass = {salt.icp_molar_mass} g/mol "
                f"[elemental cation]\n")
        r2_str = f"{calib.r2:.6f}" if calib.r2 == calib.r2 else "n/a (given, not fit)"
        f.write(f"# calibration: cond = {calib.intercept:.6f} + {calib.slope:.6f} * c_mM "
                f"(source: {calib.source}, n={calib.n_points}, R^2={r2_str})\n")
        f.write(f"# ionic strength: I = {salt.i_multiplier} * c_salt (valence-aware, z={salt.z_cation})\n")
        f.write(f"# CP: solution D_s = {salt.solution_D:.6e} m^2/s, k = {k:.6e} m/s\n")
        f.write(f"# apply_cp: {apply_cp}\n")
        if provisional:
            f.write("# WARNING: calibration is APPROXIMATED (not fit from this run's own "
                    "vial data) -- see process_raw_multisalt.py module docstring; PROVISIONAL\n")
        df.to_csv(f, index=False)


def main():
    print("=" * 100)
    print("Phase 1 CP redo: CaCl2/LaCl3 + NaCl-validation-twin preprocessing, CP on & off")
    print("=" * 100)

    validation_rows = []
    shift_rows = []

    for ds in DATASETS:
        path = DATA_DIR / ds["file"]
        raw = pp.load_raw_sheet(path, ds["sheet"])
        salt = ds["salt"]
        provisional = ds.get("provisional_calib", False)

        calib = ds["calibration"]
        if calib is None:
            calib = pp.fit_calibration(raw.vials, icp_molar_mass=salt.icp_molar_mass)

        print(f"\n--- {ds['key']} ({ds['file']}::{ds['sheet']}, salt={salt.name}, n={len(raw.ts)}) ---")
        if calib is None:
            print("  NO CALIBRATION -- skipping (insufficient vial data, no fallback given)")
            continue
        print(f"  calibration: cond = {calib.intercept:.4f} + {calib.slope:.4f}*c  "
              f"(source={calib.source}, n={calib.n_points}, R2={calib.r2})"
              + ("  [PROVISIONAL/APPROXIMATED]" if provisional else ""))
        print(f"  I multiplier = {salt.i_multiplier}, solution D_s = {salt.solution_D:.4e} m^2/s, "
              f"k = {pp.mass_transfer_coefficient(D=salt.solution_D):.4e} m/s")

        result_noncp = pp.process_experiment(raw, apply_cp=False, jw_window=JW_WINDOW,
                                              calibration=calib, salt=salt)
        result_cp = pp.process_experiment(raw, apply_cp=True, jw_window=JW_WINDOW,
                                           calibration=calib, salt=salt)

        out_noncp = OUT_DIR / f"{ds['out_prefix']}_nonCP.csv"
        out_cp = OUT_DIR / f"{ds['out_prefix']}_CP.csv"
        write_csv(out_noncp, result_noncp.df, salt, calib, False, ds["file"], ds["sheet"], provisional)
        write_csv(out_cp, result_cp.df, salt, calib, True, ds["file"], ds["sheet"], provisional)
        print(f"  Saved {out_noncp}")
        print(f"  Saved {out_cp}")

        med_shift, max_shift = pct_shift(result_noncp.df["c_int_mM"], result_cp.df["c_int_mM"])
        print(f"  CP shift in c_int (nonCP -> CP): median {med_shift:+.1f}%, max {max_shift:+.1f}%")
        shift_rows.append((ds["key"], med_shift, max_shift))

        # --- Validation: CP-off vs authors' non-CP sheet ---
        if ds.get("target_sheet"):
            n = len(raw.ts)
            ref = load_processed_sheet(DATA_DIR / "BoE Analysis.xlsx", ds["target_sheet"], n)
            print(f"  Validation (CP-off vs '{ds['target_sheet']}', non-CP authors' sheet):")
            for col, meaning, ours, refcol in [
                ("B", "c_int", result_noncp.df["c_int_mM"], ref["B"]),
                ("H", "c_p", result_noncp.df["c_p_mM"], ref["H"]),
                ("D", "ionic strength", result_noncp.df["ionic_strength_mM"], ref["D"]),
                ("J", "Jw", result_noncp.df["Jw_m3_m2_s"], ref["J"]),
                ("F", "rejection", result_noncp.df["rejection"], ref["F"]),
                ("K", "B coeff", result_noncp.df["B_um_s"], ref["K"]),
            ]:
                e = rel_err(ours, refcol)
                if len(e) == 0:
                    print(f"    {col} ({meaning}): no comparable (finite/nonzero) reference values")
                    continue
                print(f"    {col} ({meaning}): median rel err = {np.median(e):.2%}, "
                      f"max = {np.max(e):.2%}, n={len(e)}")
                validation_rows.append((ds["key"], "CP-off vs non-CP", col, meaning,
                                         np.median(e), np.max(e), len(e)))

        # --- Validation: CP-on vs authors' already-CP sheet (NaCl only) ---
        if ds.get("cp_target_sheet"):
            n = len(raw.ts)
            ref = load_processed_sheet(DATA_DIR / "BoE Analysis.xlsx", ds["cp_target_sheet"], n)
            print(f"  Validation (CP-on vs '{ds['cp_target_sheet']}', authors' CP-corrected sheet):")
            for col, meaning, ours, refcol in [
                ("B", "c_int (CP)", result_cp.df["c_int_mM"], ref["B"]),
                ("H", "c_p", result_cp.df["c_p_mM"], ref["H"]),
            ]:
                e = rel_err(ours, refcol)
                if len(e) == 0:
                    print(f"    {col} ({meaning}): no comparable (finite/nonzero) reference values")
                    continue
                print(f"    {col} ({meaning}): median rel err = {np.median(e):.2%}, "
                      f"max = {np.max(e):.2%}, n={len(e)}")
                validation_rows.append((ds["key"], "CP-on vs CP", col, meaning,
                                         np.median(e), np.max(e), len(e)))

    print("\n" + "=" * 100)
    print("Summary: CP shift per dataset (median/max %% change in c_int, nonCP -> CP)")
    print("=" * 100)
    for key, med, mx in shift_rows:
        print(f"  {key:<24} median {med:+7.1f}%   max {mx:+8.1f}%")

    print("\n" + "=" * 100)
    print("Summary: validation errors")
    print("=" * 100)
    for key, kind, col, meaning, med, mx, n in validation_rows:
        print(f"  {key:<24} {kind:<20} {col} ({meaning:<14}) median {med:.2%}  max {mx:.2%}  n={n}")


if __name__ == "__main__":
    main()
