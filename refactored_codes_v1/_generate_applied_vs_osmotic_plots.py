"""Generate applied-pressure-vs-osmotic-pressure plots for every DATA3
single-salt sheet, without re-running fits.

The plot is built entirely from measured cF (Shedlovsky-inverted from the
retentate conductivity probe) — no model prediction is needed.  Its purpose
is to explain the collaborator-flagged per-vial mass pattern:
  - Concentrating regime: cF rises → Δπ rises → ΔP-σ·Δπ shrinks → mass per
    vial DECREASES over the run.
  - Diluting regime: cF falls → Δπ falls → ΔP-σ·Δπ grows → mass per vial
    INCREASES over the run.

Output goes into the standard campaign_figures/ directory alongside the
existing mass-/concentration-/pressure-/osmotic-/conductivity-*.png files.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import refactored_ucb_library as lib  # noqa: E402

DATA_ROOT = ROOT.parent / "UnifiedFramework" / "ExperimentalDataFiles"
SAVE_DIR = (
    ROOT.parent / "UnifiedFramework" / "DATA3" / "results"
    / "paper_artifacts" / "nf270" / "campaign_figures"
)
SAVE_DIR.mkdir(parents=True, exist_ok=True)

SHEETS = [
    ("NF270_MC2", "05.07.24_NaCl"),
    ("NF270_MC2", "05.07.24_CaCl2"),
    ("NF270_MC2", "05.21.24_LaCl3"),
    ("NF270_MC3", "07.09.24_NaCl"),
    ("NF270_MC3", "07.22.24_SNaCl"),
    ("NF270_MC3", "07.11.24_CaCl2"),
    ("NF270_MC3", "07.11.24_SCaCl2"),
    ("NF270_MC3", "07.12.24_S2CaCl2"),
    ("NF270_MC4", "07.11.24_SNaCl"),
    ("NF270_MC4", "07.11.24_SLaCl3"),
    ("NF270_MC5", "07.23.24_NaCl"),
    ("NF270_MC5", "07.23.24_SNaCl"),
    ("NF270_MC5", "07.23.24_S2NaCl"),
]


def main():
    print(f"[applied-vs-osmotic] generating plots for {len(SHEETS)} single-salt sheets")
    print(f"[applied-vs-osmotic] save_dir = {SAVE_DIR}")
    print()
    all_paths = []
    n_ok = 0
    n_err = 0
    for wb_stem, sheet in SHEETS:
        run_id = f"{wb_stem.replace('NF270_', '')}.{sheet}"
        wb_path = DATA_ROOT / f"{wb_stem}.xlsx"
        try:
            ds = lib.loadxlsx(str(wb_path), sheet=sheet)["data_stru"]
        except Exception as exc:
            print(f"  [{run_id}] LOAD ERROR: {type(exc).__name__}: {exc}")
            n_err += 1
            continue
        results = {
            "data": ds,
            "sim_stru": [],   # not needed — plot uses measured cF only
            "model_settings": {"workflow_family": "DATA3"},
            "fit_meta": {},
        }
        try:
            paths = lib.run_data3_applied_vs_osmotic_plots(results, save_dir=SAVE_DIR, show=False)
        except Exception as exc:
            print(f"  [{run_id}] PLOT ERROR: {type(exc).__name__}: {exc}")
            n_err += 1
            continue
        if not paths:
            print(f"  [{run_id}] (no output)")
            n_err += 1
            continue
        all_paths.extend(paths)
        n_ok += 1
        print(f"  [{run_id}] OK  →  {Path(paths[0]).name}")
    print()
    print(f"[applied-vs-osmotic] {n_ok} ok / {n_err} err out of {len(SHEETS)}")
    print(f"[applied-vs-osmotic] generated files:")
    for p in all_paths:
        print(f"   {p}")


if __name__ == "__main__":
    main()
