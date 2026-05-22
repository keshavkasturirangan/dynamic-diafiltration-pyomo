"""Roll out the paper-style concentration plot to all NF270 single-salt sheets.

Re-renders each sheet's mass + concentration figures using the published-style
markers/colors/legend. Output goes into the existing campaign_figures/ folder,
overwriting the older plots in place.
"""
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import refactored_ucb_library as lib  # noqa: E402

SAVE_DIR = (
    ROOT.parent / "UnifiedFramework" / "DATA3" / "results"
    / "paper_artifacts" / "nf270" / "campaign_figures"
)
SAVE_DIR.mkdir(parents=True, exist_ok=True)

SINGLE_SALT_RUN_IDS = [
    "MC2.05.07.24_NaCl",
    "MC3.07.22.24_SNaCl",
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

only = [f"nf270.{r}" for r in SINGLE_SALT_RUN_IDS]
t0 = time.time()
print(f"[rollout] starting paper-style rerender for {len(only)} sheets")
print(f"[rollout] save_dir = {SAVE_DIR}")
results = lib.materialize_all(
    campaign="NF270",
    save_dir=str(SAVE_DIR),
    only=only,
)
for r in results:
    name = r.get("name", "?")
    status = r.get("status", "?")
    arts = r.get("artifacts") or r.get("error") or ""
    print(f"  {name:40s}  {status:8s}  {arts}")
print(f"[rollout] done in {(time.time()-t0)/60:.1f} min")
