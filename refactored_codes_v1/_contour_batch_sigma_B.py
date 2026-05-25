"""Run the 11-sheet contour batch, filling in the missing sigma-B slice.

Per-slice checkpointing in run_nf270_contour_for_sheet means the existing
B-Lp and sigma-Lp CSVs on disk get reused; only the new sigma-B slice gets
freshly computed (one ~10-15 min sweep per sheet). All three slices are
re-rendered in the new DATA1-style format as a bonus.

Total expected wall time: ~2-3 hours for 11 sheets.
"""
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import refactored_ucb_library as lib  # noqa: E402

PANELS_DIR = (
    ROOT.parent / "UnifiedFramework" / "DATA3" / "results"
    / "paper_artifacts" / "nf270" / "contour_panels"
)
PANELS_DIR.mkdir(parents=True, exist_ok=True)

SINGLE_SALT_RUN_IDS = [
    "MC3.07.22.24_SNaCl",   # dilution baseline — fast smoke test on the well-fit sheet
    "MC2.05.07.24_NaCl",
    "MC4.07.11.24_SNaCl",
    "MC5.07.23.24_NaCl",
    "MC5.07.23.24_SNaCl",
    "MC5.07.23.24_S2NaCl",
    "MC2.05.07.24_CaCl2",
    "MC3.07.11.24_SCaCl2",
    "MC3.07.12.24_S2CaCl2",
    "MC2.05.21.24_LaCl3",
    "MC4.07.11.24_SLaCl3",
]

t0 = time.time()
print(f"[contour-batch] starting σ-B fill + DATA1-style re-render for {len(SINGLE_SALT_RUN_IDS)} sheets")
print(f"[contour-batch] panels dir = {PANELS_DIR}")

lib.run_nf270_contour_branch(
    SINGLE_SALT_RUN_IDS,
    save_dir=PANELS_DIR,
    grid_density=10,   # matches the previous batch's resolution
    nfe=80,
)

print(f"[contour-batch] done in {(time.time()-t0)/60:.1f} min")
