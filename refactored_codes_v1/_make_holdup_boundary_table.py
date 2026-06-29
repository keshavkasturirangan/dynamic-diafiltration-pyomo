#!/usr/bin/env python3
"""Generate the curated holdup-boundary table + per-sheet diagnostic plots.

DATA2-style curated holdup boundaries for the NF270 / DATA3 sheets. This script
loads each single-salt sheet, records the AUTO-detected holdup boundary (the
mass-rise split the loader uses by default), writes a pre-filled CSV, and saves
a vial-1 mass-vs-time diagnostic per sheet with the boundary marked.

Workflow to curate:
  1. python3 _make_holdup_boundary_table.py
  2. Open nf270_holdup_boundaries.csv, inspect the diagnostic PNGs, and fill the
     `holdup_end_s` column with the TRUE holdup-fill end time (seconds, sheet
     time axis) for any sheet you want to override. Leave blank to keep auto.
  3. Re-run the loader / fits — _nf270_curated_holdup_end_s() honors the column.

Usage: python3 _make_holdup_boundary_table.py
Output: nf270_holdup_boundaries.csv, holdup_diagnostics/<run_id>.png
"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib

OUT_CSV = HERE / "nf270_holdup_boundaries.csv"
DIAG_DIR = HERE / "holdup_diagnostics"
DIAG_DIR.mkdir(exist_ok=True)

# The 11 single-salt campaign sheets (same source of truth as the runfile).
RUN_IDS = [
    "MC2.05.07.24_NaCl", "MC3.07.22.24_SNaCl", "MC4.07.11.24_SNaCl",
    "MC5.07.23.24_NaCl", "MC5.07.23.24_SNaCl", "MC5.07.23.24_S2NaCl",
    "MC2.05.07.24_CaCl2", "MC3.07.11.24_SCaCl2", "MC3.07.12.24_S2CaCl2",
    "MC2.05.21.24_LaCl3", "MC4.07.11.24_SLaCl3",
]


def _load(rid):
    fam = lib.NF270_RUN_REGISTRY[rid]
    ds = lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]
    return fam, ds


def main():
    rows = []
    for rid in RUN_IDS:
        fam, ds = _load(rid)
        cfg = ds["data_config"]
        n_hold = int(cfg.get("n_holdup_samples", 0) or 0)
        t_start = cfg.get("t_perm_start_s")
        source = cfg.get("holdup_boundary_source")
        v1 = ds["data_raw"][0]
        t1 = np.asarray(v1["time"], dtype=float).reshape(-1)
        m1 = np.asarray(v1["mass"], dtype=float).reshape(-1)
        # When the loader split vial 1, data_raw[0] is the holdup vial; rejoin it
        # with the real vial (data_raw[1]) to draw the full pre-split vial 1.
        if cfg.get("vial1_split_status") == "split" and len(ds["data_raw"]) > 1:
            v2 = ds["data_raw"][1]
            t2 = np.asarray(v2["time"], dtype=float).reshape(-1)
            m2 = np.asarray(v2["mass"], dtype=float).reshape(-1)
            t_full = np.concatenate([t1, t2])
            # de-reset the real-vial mass back onto the holdup tail for a clean line
            offset = (m1[-1] if m1.size else 0.0)
            m_full = np.concatenate([m1, m2 + offset])
        else:
            t_full, m_full = t1, m1

        rows.append({
            "run_id": rid,
            "workbook": fam["workbook"],
            "sheet": fam["sheet"],
            "mode": ds.get("mode"),
            "auto_split_status": cfg.get("vial1_split_status"),
            "auto_holdup_end_s": round(float(t_start), 1) if t_start is not None else "",
            "auto_n_holdup": n_hold,
            "holdup_end_s": "",   # <-- experimenter fills this to OVERRIDE
            "notes": "",
        })

        # Diagnostic plot: vial-1 mass vs time with the auto boundary marked.
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(t_full, m_full, ".", ms=4, color="tab:blue", label="vial-1 mass")
        if t_start is not None:
            ax.axvline(float(t_start), color="red", ls="--", lw=1.5,
                       label=f"auto boundary = {float(t_start):.0f} s  (n_holdup={n_hold})")
        ax.axhline(0.05, color="gray", ls=":", lw=1, label="0.05 g rise threshold")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Cumulative mass [g]")
        ax.set_title(f"{rid}  ({ds.get('mode')}, source={source})")
        ax.legend(fontsize=8, loc="best")
        fig.tight_layout()
        out = DIAG_DIR / f"{rid}.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"  {rid:24s} auto_end={str(rows[-1]['auto_holdup_end_s']):>7} s  n_holdup={n_hold:3d}  -> {out.name}")

    df = pd.DataFrame(rows, columns=[
        "run_id", "workbook", "sheet", "mode", "auto_split_status",
        "auto_holdup_end_s", "auto_n_holdup", "holdup_end_s", "notes",
    ])
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWROTE {OUT_CSV}  ({len(df)} sheets)")
    print(f"Diagnostics in {DIAG_DIR}/  — fill `holdup_end_s` to override; blank keeps auto.")


if __name__ == "__main__":
    main()
