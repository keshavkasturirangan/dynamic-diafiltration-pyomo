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

# ─────────────────────────────────────────────────────────────────────────────
# DATA3-specific add-ons (RE-ENABLED — 2026-05-24 evening decision)
# ─────────────────────────────────────────────────────────────────────────────
# Earlier today we tried the DATA2 recipe verbatim and it failed to
# generalize: the DATA2 paper-demo theta + B_form=1 only converged on the
# clean dilution baseline (MC3.07.22.24_SNaCl). Concentration-regime
# sheets thrashed in IPOPT regardless of seed strategy. Per the pivot:
# re-enable the two DATA3-specific knobs that ARE validated to produce
# sensible fits, and accept the deviation from pure DATA2-recipe parity.
# Comment out either knob to revert to the DATA2-recipe path.
#
# (a) Contour-seeded multistart — uses per-sheet contour-grid minima as
#     IPOPT warm starts. Each sheet's "best basin" is on disk already
#     from the May 22 contour batch (paper_artifacts/nf270/contour_panels).
#     Avoids degenerate B=1e-6 / σ=wall starts.
lib.NF270_MULTISTART_USE_CONTOUR_SEEDS = True

# (b) cF residual floor — caps the cF residual scale at 1 mM (the default
#     0.3 %-relative weight is tighter than the measurement uncertainty
#     above cF ≈ 100 mM, pulling σ to the wall). Validated on MC4 SNaCl:
#     σ moved 0.0 → 0.461 (interior), obj_cr improved 200×, all three
#     channels agreed on the same θ.
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0

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
print(f"[rollout] DATA2-recipe mode: NF270_CF_RESIDUAL_FLOOR_MM={lib.NF270_CF_RESIDUAL_FLOOR_MM}, "
      f"NF270_MULTISTART_USE_CONTOUR_SEEDS={lib.NF270_MULTISTART_USE_CONTOUR_SEEDS}")
# multistart_iterations=1: v11 showed all 12 manifest LHS restarts landing
# on the same θ for sheet 1 (Lp=9.11, σ=0.66) — the deterministic seeds
# (3 contour + 1 cross-salt) are doing all the basin-finding work.  One
# LHS perturbation is enough as a safety net.
results = lib.materialize_all(
    campaign="NF270",
    save_dir=str(SAVE_DIR),
    only=only,
    extra_opts={"request_overrides": {"multistart_iterations": 1}},
)
for r in results:
    name = r.get("name", "?")
    status = r.get("status", "?")
    arts = r.get("artifacts") or r.get("error") or ""
    print(f"  {name:40s}  {status:8s}  {arts}")
print(f"[rollout] done in {(time.time()-t0)/60:.1f} min")
