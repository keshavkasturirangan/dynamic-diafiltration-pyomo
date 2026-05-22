"""One-shot preview of paper-styled concentration plot for MC3.07.22.24_SNaCl.

Disposable. Renders the single sheet, then writes the PNG into a sibling
preview directory so it does not collide with the existing campaign_figures.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import refactored_ucb_library as lib  # noqa: E402

PREVIEW_DIR = ROOT.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270" / "paper_styling_preview"
PREVIEW_DIR.mkdir(parents=True, exist_ok=True)

results = lib.materialize_all(
    campaign="NF270",
    save_dir=str(PREVIEW_DIR),
    only=["nf270.MC3.07.22.24_SNaCl"],
)

for r in results:
    print(r.get("name"), "->", r.get("status"), r.get("artifacts") or r.get("error"))
