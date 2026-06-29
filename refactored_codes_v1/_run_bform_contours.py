#!/usr/bin/env python3
"""DATA3 NF270 B-level identifiability contours — the three pairwise WSSE maps
the seeding workflow reads:  (B, Lp), (sigma, Lp), (sigma, B), with B as a single
LEVEL (shape held fixed) via the lumped B_form='single' engine path.  Forward-sim
only (no fitting) in the canonical figure_s5 line style; CSV sidecars saved.

Cross-sheet parallel (each sheet's grid is sequential internally), IPOPT pinned to
one thread per worker so N workers use N cores cleanly.

Writes ONLY to paper_artifacts/nf270/bform_study/contours/<run_id>/.
Usage: python3 _run_bform_contours.py [grid=30] [workers=6] [sheets=all|<rid,..>]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, traceback
from pathlib import Path
from multiprocessing import Pool
import matplotlib
matplotlib.use("Agg")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib

ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
STUDY = ART / "bform_study"
CONT = STUDY / "contours"
CONT.mkdir(parents=True, exist_ok=True)

SHEET_IDS = [
    "MC3.07.22.24_SNaCl", "MC5.07.23.24_S2NaCl", "MC4.07.11.24_SNaCl",
    "MC2.05.07.24_NaCl", "MC5.07.23.24_NaCl", "MC5.07.23.24_SNaCl",
    "MC2.05.07.24_CaCl2", "MC3.07.11.24_SCaCl2", "MC3.07.12.24_S2CaCl2",
    "MC2.05.21.24_LaCl3", "MC4.07.11.24_SLaCl3",
]
GRID = int(sys.argv[1]) if len(sys.argv) > 1 else 30
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 6


def one_sheet(rid):
    lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
    try:
        fam = lib.NF270_RUN_REGISTRY[rid]
        ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
        mode = ds["mode"]
        # reuse the campaign's constant-B fit to center the grid if available
        theta_fit = None
        rp = STUDY / rid / "result_single.json"
        if rp.exists():
            try:
                tj = json.loads(rp.read_text()).get("parameters")
                if tj and "Lp" in tj and "B" in tj and "sigma" in tj:
                    theta_fit = tj
            except Exception:
                theta_fit = None
        out = CONT / rid
        lib.run_nf270_contour_for_sheet(
            rid, save_dir=out, grid_density=GRID, theta_fit=theta_fit,
            B_form="single", nfe=80, mode=mode)
        return (rid, "ok")
    except Exception as e:
        traceback.print_exc()
        return (rid, f"FAIL: {e}")


def main():
    sheets = SHEET_IDS
    if len(sys.argv) > 3 and sys.argv[3] != "all":
        sheets = sys.argv[3].split(",")
    print(f"[contours] {len(sheets)} sheets, grid={GRID}, workers={WORKERS}", flush=True)
    with Pool(WORKERS) as pool:
        for rid, status in pool.imap_unordered(one_sheet, sheets):
            print(f"[contours] {rid:22s} {status}", flush=True)
    print("CONTOURS DONE", flush=True)


if __name__ == "__main__":
    main()
