"""Phase 1 rollout — per-sheet-seeded fits on all 11 NF270 single-salt sheets.

Architecture
------------
    mode    = autodetected (Lag for all 11 currently; detector in loader)
    theta   = None  → solve_model fills in per-sheet defaults from
                       data_stru["data_config"]["Lp0", "B0", "sigma0"]
    B_form  = 1     (the DATA2-paper-demo power-law form)
    workflow_family = "DATA3" (Excel loader + EC25 + time correction)

Why theta=None instead of the DATA2 demo's explicit dict?
---------------------------------------------------------
v7 attempt: passed DATA2's hardcoded theta={"Lp": 11, "beta_0": 1,
"sigma": 1.0, ...} verbatim. Worked on the dilution baseline
(MC3.07.22.24_SNaCl, 19 sec) but the first concentration-regime sheet
(MC2.05.07.24_NaCl, cF≈0.9 mM start) couldn't converge from a
σ=1.0 starting point — at very low cF the osmotic resistance is near
zero, IPOPT starts in a degenerate region and thrashes.

Solution: let solve_model pull seeds from the per-sheet data_config:
    Lp0    = 5.0       (literature midpoint for NF270)
    B0     = 0.5
    sigma0 = 0.9       (interior, safer than the σ=1.0 wall)
    beta_0 = computed from B0, Lp0, ΔP
The loader already populates these. theta=None tells solve_model to
use them instead of overriding.

The previous explicit-theta path is preserved as DATA2_RECIPE_THETA
below, commented out, in case we want to revisit it.
"""
from __future__ import annotations

import signal
import sys
import time
import traceback
from pathlib import Path


# Wall-clock timeout enforced from OUTSIDE Pyomo/IPOPT. The solver's
# internal max_cpu_time only checks at iteration boundaries, so if IPOPT
# gets stuck inside a single iteration (e.g., MUMPS linear-solver
# factorization on a poorly-conditioned step), the cap never fires.
# signal.alarm is a brute-force interrupt that kills the call no matter
# where it is. Unix-only; fine for our macOS workflow.
PER_FIT_WALL_TIMEOUT_S = 120  # 2 min hard cap

class _FitTimeout(Exception):
    pass

def _alarm_handler(signum, frame):
    raise _FitTimeout("wall-clock timeout reached")

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import refactored_ucb_library as lib  # noqa: E402

# DATA3-specific add-ons stay INERT for this rollout — we use DATA2 recipe.
assert lib.NF270_CF_RESIDUAL_FLOOR_MM is None, "cf-floor should be disabled"
assert lib.NF270_MULTISTART_USE_CONTOUR_SEEDS is False, "contour seeds should be disabled"

SAVE_DIR = (
    ROOT.parent / "UnifiedFramework" / "DATA3" / "results"
    / "paper_artifacts" / "nf270" / "campaign_figures"
)
SAVE_DIR.mkdir(parents=True, exist_ok=True)

NF270_ROOT = ROOT.parent / "UnifiedFramework" / "ExperimentalDataFiles"

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

# DATA2 paper-demo theta verbatim — KEPT for reference; v7 used this
# explicitly and the σ=1.0 starting point thrashed on low-cF concentration
# sheets. We now pass theta=None and let the loader's per-sheet defaults
# (data_config["Lp0", "B0", "sigma0"]) take over.
# DATA2_RECIPE_THETA = {
#     "Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01,
#     "sigma": 1.0, "S0": 0,
# }
USE_EXPLICIT_THETA = False   # set True + uncomment above to re-test the v7 path

t0 = time.time()
print(f"[rollout-D2] starting per-sheet-seeded rollout for {len(SINGLE_SALT_RUN_IDS)} sheets")
print(f"[rollout-D2] save_dir = {SAVE_DIR}")
print(f"[rollout-D2] theta    = None (use data_config['Lp0','B0','sigma0'] per sheet)")
print(f"[rollout-D2] B_form   = 1, workflow_family = DATA3")
print()

summary = []  # (run_id, status, dt_s, fitted_params)
for rid in SINGLE_SALT_RUN_IDS:
    fam = lib.NF270_RUN_REGISTRY.get(rid)
    if fam is None:
        print(f"  [{rid}] SKIP — not in registry")
        summary.append((rid, "skip-not-in-registry", 0.0, {}))
        continue

    wb_path = NF270_ROOT / fam["workbook"]
    if not wb_path.exists():
        print(f"  [{rid}] SKIP — workbook missing: {wb_path}")
        summary.append((rid, "skip-workbook-missing", 0.0, {}))
        continue

    print(f"  [{rid}]  loading...")
    t_case = time.time()
    try:
        ds = lib.loadxlsx(wb_path, sheet=fam["sheet"])["data_stru"]
    except Exception as exc:
        print(f"    LOAD ERROR ({time.time()-t_case:.1f}s): {type(exc).__name__}: {exc}")
        summary.append((rid, f"load-error: {type(exc).__name__}", time.time()-t_case, {}))
        continue

    print(f"    mode={ds['mode']}, n_vials={ds['data_config']['n']}, "
          f"M_F0={ds['data_config']['M_F0']:.2f}, M_O={ds['data_config']['M_O']:.3f}")

    # Pass theta=None unless USE_EXPLICIT_THETA is True (kept for diagnostics).
    # When theta=None, solve_model reads Lp0/B0/sigma0/etc. from data_config
    # — the per-sheet defaults the loader already computed.
    #
    # Two-layer timeout protection:
    #   1. solver_max_cpu_time=90 — IPOPT internal cap (checked at iter
    #      boundaries; fails to fire if IPOPT is stuck inside an iteration).
    #   2. signal.alarm(120) — hard wall-clock cap enforced from outside
    #      Pyomo. Wraps the solve_model call so we can guarantee no single
    #      fit hangs the rollout, even when IPOPT's internal cap doesn't fire.
    signal.signal(signal.SIGALRM, _alarm_handler)
    signal.alarm(PER_FIT_WALL_TIMEOUT_S)
    try:
        fit_stru, sim_stru, sim_inter = lib.solve_model(
            ds, ds["mode"], None,
            sim_opt=False,
            B_form=1,
            LOUD=False,
            workflow_family="DATA3",
            solver_max_cpu_time=90,
        )
        signal.alarm(0)         # cancel alarm on success
    except _FitTimeout:
        signal.alarm(0)
        dt = time.time() - t_case
        print(f"    FIT WALL-TIMEOUT after {dt:.1f}s (cap was {PER_FIT_WALL_TIMEOUT_S}s)")
        summary.append((rid, f"wall-timeout-{PER_FIT_WALL_TIMEOUT_S}s", dt, {}))
        continue
    except Exception as exc:
        dt = time.time() - t_case
        print(f"    FIT ERROR ({dt:.1f}s): {type(exc).__name__}: {exc}")
        summary.append((rid, f"fit-error: {type(exc).__name__}", dt, {}))
        continue

    dt = time.time() - t_case
    params = fit_stru.get("parameters", {})
    # Flatten any dict-valued entries (e.g. B per vial) for the summary print.
    p_flat = {}
    for k, v in params.items():
        if isinstance(v, (int, float)):
            p_flat[k] = float(v)
    print(f"    FIT OK ({dt:.1f}s)  "
          f"Lp={p_flat.get('Lp', float('nan')):.3f}, "
          f"sigma={p_flat.get('sigma', float('nan')):.3f}, "
          f"beta_0={p_flat.get('beta_0', float('nan')):.3g}, "
          f"obj_m={fit_stru.get('obj_m', 0):.2g}/obj_cv={fit_stru.get('obj_cv', 0):.2g}/obj_cr={fit_stru.get('obj_cr', 0):.2g}")

    # Render the paper-style plots via the existing pipeline.
    try:
        results_dict = {
            "data": ds,
            "sim_stru": sim_stru,
            "model_settings": {"workflow_family": "DATA3"},
        }
        plot_paths = lib.run_data3_time_series_plots(
            results_dict, save_dir=SAVE_DIR, show=False,
        )
        print(f"    PLOT OK ({len(plot_paths)} files)")
    except Exception as exc:
        print(f"    PLOT ERROR: {type(exc).__name__}: {exc}")
        traceback.print_exc()

    summary.append((rid, "ok", dt, p_flat))

print()
print(f"[rollout-D2] === SUMMARY ({(time.time()-t0)/60:.1f} min total) ===")
n_ok = sum(1 for s in summary if s[1] == "ok")
print(f"[rollout-D2] {n_ok}/{len(summary)} sheets fit successfully")
print()
print(f"{'Run ID':<24}  {'status':<22}  {'dt(s)':>6}  {'Lp':>6}  {'sigma':>6}  {'beta_0':>8}  {'beta_1':>8}")
print("-" * 96)
for rid, status, dt, p in summary:
    Lp = p.get("Lp", float("nan"))
    sig = p.get("sigma", float("nan"))
    b0 = p.get("beta_0", float("nan"))
    b1 = p.get("beta_1", float("nan"))
    print(f"{rid:<24}  {status:<22}  {dt:>6.1f}  "
          f"{Lp:>6.2f}  {sig:>6.3f}  {b0:>8.3g}  {b1:>8.3g}")
