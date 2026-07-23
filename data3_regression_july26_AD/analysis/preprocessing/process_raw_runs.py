"""
Apply the validated preprocessing pipeline (preprocess_nacl.py) to the 6
raw-only NaCl runs, fitting each run's OWN per-experiment calibration from
its vial (conductivity, ICP) pairs (all 6 have per-vial conductivity
recorded, unlike the two validation-target sheets -- see
analysis/preprocessing/README.md).

Writes TWO CSVs per run to analysis/preprocessing/processed/:
  - `<name>.csv` (unchanged filename/content): apply_cp=False (non-CP, bulk
    == interfacial), matching the convention Step 2/Step 3 and
    analysis/nacl_all_datasets/regress_all_nacl.py already consume (the "(2)"
    sheet) -- kept byte-for-byte identical (docs/prompts/
    cp_preprocessing_extend.md Phase 1 regression test).
  - `<name>_CP.csv` (new): apply_cp=True, same fitted calibration -- the
    CP-corrected twin, for Phase 2's CP-vs-non-CP comparison. Not consumed
    by any current regression driver.

Columns: c_int_mM, ionic_strength_mM, rejection, c_p_mM, Jw_m3_m2_s, B_um_s,
time_s, ret_cond_uS_cm, perm_cond_uS_cm, with the fitted calibration and
apply_cp flag recorded in a header comment.

Run with: conda activate data3-regression && python analysis/preprocessing/process_raw_runs.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import preprocess_nacl as pp

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
OUT_DIR = Path(__file__).resolve().parent / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Jw smoothing window tuned in validate.py against the one reference J
# column we have (05.27.26_NaCl (2)); reused here for consistency.
JW_WINDOW = 51

RUNS = [
    {"file": "NF270_MC3.xlsx", "sheet": "07.09.24_NaCl", "out": "MC3_07.09.24_NaCl.csv"},
    {"file": "NF270_MC3.xlsx", "sheet": "07.22.24_SNaCl", "out": "MC3_07.22.24_SNaCl.csv"},
    {"file": "NF270_MC4.xlsx", "sheet": "07.11.24_SNaCl", "out": "MC4_07.11.24_SNaCl.csv"},
    {"file": "NF270_MC5.xlsx", "sheet": "07.23.24_NaCl", "out": "MC5_07.23.24_NaCl.csv"},
    {"file": "NF270_MC5.xlsx", "sheet": "07.23.24_SNaCl", "out": "MC5_07.23.24_SNaCl.csv"},
    {"file": "NF270_MC5.xlsx", "sheet": "07.23.24_S2NaCl", "out": "MC5_07.23.24_S2NaCl.csv"},
]


def main():
    print("=" * 90)
    print("Processing 6 raw-only NaCl runs with per-experiment fitted calibrations")
    print("=" * 90)
    summary = []
    for run in RUNS:
        path = DATA_DIR / run["file"]
        raw = pp.load_raw_sheet(path, run["sheet"])
        result = pp.process_experiment(raw, apply_cp=False, jw_window=JW_WINDOW)
        calib = result.calibration

        out_path = OUT_DIR / run["out"]
        with open(out_path, "w") as f:
            if calib is not None:
                f.write(f"# source: {run['file']} / {run['sheet']}\n")
                f.write(f"# calibration: cond = {calib.intercept:.6f} + "
                        f"{calib.slope:.6f} * c_mM  (fit from {calib.n_points} "
                        f"vial points, R^2={calib.r2:.6f}, ICP molar mass = "
                        f"{pp.ICP_MOLAR_MASS} g/mol [Na+])\n")
                f.write("# apply_cp: False (non-CP, bulk == interfacial)\n")
            else:
                f.write(f"# source: {run['file']} / {run['sheet']}\n")
                f.write("# WARNING: calibration could not be fit (insufficient "
                        "vial data) -- concentration columns are NaN; PROVISIONAL\n")
            result.df.to_csv(f, index=False)

        # CP-corrected twin, same fitted calibration, same Jw window --
        # written alongside for Phase 2's CP-vs-non-CP comparison (Phase 1
        # scope: reprocess every single-salt dataset from raw with CP
        # applied; docs/prompts/cp_preprocessing_extend.md).
        cp_out_path = OUT_DIR / run["out"].replace(".csv", "_CP.csv")
        if calib is not None:
            cp_result = pp.process_experiment(raw, apply_cp=True, jw_window=JW_WINDOW,
                                               calibration=calib)
            with open(cp_out_path, "w") as f:
                f.write(f"# source: {run['file']} / {run['sheet']}\n")
                f.write(f"# calibration: cond = {calib.intercept:.6f} + "
                        f"{calib.slope:.6f} * c_mM  (fit from {calib.n_points} "
                        f"vial points, R^2={calib.r2:.6f}, ICP molar mass = "
                        f"{pp.ICP_MOLAR_MASS} g/mol [Na+])\n")
                f.write(f"# salt: NaCl, I multiplier = {pp.SALT_NACL.i_multiplier}, "
                        f"CP solution D_s = {pp.SALT_NACL.solution_D:.4e} m^2/s, "
                        f"k = {pp.mass_transfer_coefficient(D=pp.SALT_NACL.solution_D):.4e} m/s\n")
                f.write("# apply_cp: True (CP-corrected retentate-side interfacial concentration)\n")
                cp_result.df.to_csv(f, index=False)

        n_flags = len(result.flags)
        cal_str = (f"cond={calib.intercept:.3f}+{calib.slope:.3f}*c  "
                   f"(n={calib.n_points}, R2={calib.r2:.5f})"
                   if calib is not None else "NO CALIBRATION (provisional)")
        print(f"{run['sheet']:<20} n={len(raw.ts):>5}  {cal_str}")
        for flag in result.flags:
            print(f"  FLAG: {flag}")
        summary.append((run["sheet"], len(raw.ts), calib, result.flags))
        print(f"  Saved {out_path}")
        if calib is not None:
            print(f"  Saved {cp_out_path}")

    print()
    print("=" * 90)
    print("Calibration summary (compare to embedded MC5 05.27.26 formula: "
          "cond = -63.706 + 76.685*c)")
    print("=" * 90)
    for sheet, n, calib, flags in summary:
        if calib is not None:
            print(f"{sheet:<20} intercept={calib.intercept:>10.3f}  "
                  f"slope={calib.slope:>10.3f}  n_vial_pts={calib.n_points}  "
                  f"R2={calib.r2:.5f}  {'PROVISIONAL' if flags else ''}")
        else:
            print(f"{sheet:<20} NO CALIBRATION FIT (provisional/skip)")


if __name__ == "__main__":
    main()
