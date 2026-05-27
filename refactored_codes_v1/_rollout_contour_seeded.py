"""Phase 1 (v3): Contour-seeded rollout on ALL NaCl single-salt sheets.

Recipe = v5 (the run that produced Lp=9.112, σ=0.665 for MC2.NaCl):
  B_form         = 'single'   (scalar B per sheet, NOT convection)
  cF-floor       = 1.0
  perm-probe     = OFF
  contour seeds  = ON  (per-sheet panel minima as initial guesses)
  cross-salt     = ON  (NaCl ref injected as one of the seeds)

Initial-guess source: per-sheet contour maps.
  - For sheets WITH a cF-floor contour panel
    (MC3.SNaCl, MC4.SNaCl, MC5.S2NaCl): use cF-floor panel
  - For sheets WITHOUT a cF-floor panel (MC2.NaCl, MC5.NaCl, MC5.SNaCl):
    fall back to legacy contour_panels/ — same as v5's behavior
  This is wired in refactored_ucb_library._resolve_nf270_panel_dir.

Sheets in scope (6 NaCl sheets, ordered for diagnostic value):
  1. MC2.05.07.24_NaCl   — SMOKE TEST: must reproduce v5's σ=0.665, else recipe is still wrong
  2. MC3.07.22.24_SNaCl  — dilution, user-marked "already done"; cf-floor panel
  3. MC4.07.11.24_SNaCl  — concentration, cf-floor panel
  4. MC5.07.23.24_NaCl   — concentration, legacy panel only
  5. MC5.07.23.24_SNaCl  — concentration, legacy panel only
  6. MC5.07.23.24_S2NaCl — dilution, cf-floor panel

CaCl2 / LaCl3 sheets are deliberately out of scope — separate strategy planned.
"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

NF270_ROOT = ROOT.parent / "UnifiedFramework" / "ExperimentalDataFiles"
SAVE_DIR = (
    ROOT.parent / "UnifiedFramework" / "DATA3" / "results"
    / "paper_artifacts" / "nf270" / "campaign_figures"
)
SAVE_DIR.mkdir(parents=True, exist_ok=True)

SCRATCH = Path("/tmp/data3_contour_seeded")
SCRATCH.mkdir(parents=True, exist_ok=True)

WORKER_SCRIPT = ROOT / "_contour_seeded_worker.py"

# All 6 NaCl single-salt sheets, ordered smoke-test → user's TODO list.
# Library auto-resolves per-sheet contour panel (cF-floor preferred when
# available, legacy fallback). CaCl2 / LaCl3 are out of scope.
CF_FLOOR_SHEETS = [
    "MC2.05.07.24_NaCl",    # SMOKE TEST: must hit Lp≈9.11, σ≈0.66
    "MC3.07.22.24_SNaCl",   # cf-floor panel; dilution
    "MC4.07.11.24_SNaCl",   # cf-floor panel; concentration
    "MC5.07.23.24_NaCl",    # legacy panel only; concentration
    "MC5.07.23.24_SNaCl",   # legacy panel only; concentration
    "MC5.07.23.24_S2NaCl",  # cf-floor panel; dilution
]

PER_SHEET_TIMEOUT_S = 900   # 11 starts × ~30s each + overhead
MULTISTART_ITERATIONS = 8   # v7 slide 29 recipe: 2 contour + 1 cross-salt + 8 LHS = 11 starts/sheet


def run_one(run_id, fam, timeout_s, multistart_iter):
    """Spawn worker in its own process group; SIGKILL the whole group on timeout."""
    wb_path = NF270_ROOT / fam["workbook"]
    if not wb_path.exists():
        return {"run_id": run_id, "error": f"workbook missing: {wb_path}"}
    out_json = SCRATCH / f"{run_id.replace('/', '_')}.json"
    if out_json.exists():
        out_json.unlink()

    cmd = [
        sys.executable, "-u", str(WORKER_SCRIPT),
        run_id, str(wb_path), fam["sheet"], str(out_json),
        str(SAVE_DIR), str(multistart_iter),
    ]
    t0 = time.time()
    proc = subprocess.Popen(
        cmd,
        start_new_session=True,        # new process group → killpg can reach IPOPT
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    try:
        stdout, _ = proc.communicate(timeout=timeout_s)
        elapsed = time.time() - t0
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            pass
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            pass
        return {
            "run_id": run_id,
            "error": "wall-timeout",
            "fit_time_s": round(elapsed, 1),
            "timeout_s": timeout_s,
        }

    if not out_json.exists():
        return {
            "run_id": run_id,
            "error": f"worker died (rc={proc.returncode})",
            "stdout_tail": (stdout or "")[-1000:],
            "fit_time_s": round(elapsed, 1),
        }
    try:
        return json.loads(out_json.read_text())
    except Exception as exc:
        return {"run_id": run_id, "error": f"could not parse worker JSON: {exc!r}"}


def main():
    import refactored_ucb_library as lib  # for the run registry

    print(f"[contour-seeded] Phase 1 v5 — slide-29 recipe + Phase-A corrections — ALL 6 NaCl sheets")
    print(f"[contour-seeded] B_form=2 (power-law), cF-floor 3 mM (slide-26 noise), perm-probe OFF")
    print(f"[contour-seeded] contour seeds: top-1 from σ-Lp + top-1 from B-Lp (per-slice)")
    print(f"[contour-seeded] cross-salt anchor: MC3.SNaCl gold standard "
          f"(Lp=8.22, B=14.7, σ=0.45)")
    print(f"[contour-seeded] multistart_iterations={MULTISTART_ITERATIONS} "
          f"(2 contour + 1 cross-salt + {MULTISTART_ITERATIONS} LHS = {2 + 1 + MULTISTART_ITERATIONS} starts/sheet)")
    print(f"[contour-seeded] Per-sheet hard cap: {PER_SHEET_TIMEOUT_S}s (killpg-enforced)")
    print(f"[contour-seeded] save_dir = {SAVE_DIR}")
    print()

    t0 = time.time()
    summary = []
    for rid in CF_FLOOR_SHEETS:
        fam = lib.NF270_RUN_REGISTRY.get(rid)
        if fam is None:
            print(f"  [{rid}] SKIP — not in registry")
            summary.append({"run_id": rid, "status": "skip-not-in-registry"})
            continue
        print(f"  [{rid}]  starting subprocess (cap {PER_SHEET_TIMEOUT_S}s, process group)...")
        out = run_one(rid, fam, PER_SHEET_TIMEOUT_S, MULTISTART_ITERATIONS)
        if "error" in out:
            tail = ""
            if "stdout_tail" in out:
                tail = "\n      stdout_tail: " + out["stdout_tail"].replace("\n", "\n      ")
            print(f"    {out['error']}  ({out.get('fit_time_s','?')}s){tail}")
            summary.append({"run_id": rid, "status": out["error"], "out": out})
            continue
        p = out.get("parameters", {})
        plot_n = len(out.get("plot_paths", []))
        plot_err = f" [plot:{out['plot_error']}]" if "plot_error" in out else ""
        print(f"    FIT OK ({out['fit_time_s']}s)  "
              f"Lp={p.get('Lp', float('nan')):.3f}, "
              f"sigma={p.get('sigma', float('nan')):.3f}, "
              f"beta_0={p.get('beta_0', float('nan')):.3g}, "
              f"obj_m={out['obj_m']:.2g}/cv={out['obj_cv']:.2g}/cr={out['obj_cr']:.2g}, "
              f"{plot_n} plots{plot_err}")
        summary.append({"run_id": rid, "status": "ok", "out": out})

    print()
    print(f"[contour-seeded] === SUMMARY ({(time.time()-t0)/60:.1f} min) ===")
    n_ok = sum(1 for s in summary if s["status"] == "ok")
    print(f"[contour-seeded] {n_ok}/{len(summary)} sheets fit successfully\n")
    print(f"{'Run ID':<24}  {'status':<22}  {'Lp':>6}  {'sigma':>6}  {'beta_0':>10}  {'obj_cr':>10}")
    print("-" * 90)
    for s in summary:
        if s["status"] == "ok":
            p = s["out"]["parameters"]
            Lp = p.get("Lp", float("nan"))
            sig = p.get("sigma", float("nan"))
            b0 = p.get("beta_0", float("nan"))
            cr = s["out"]["obj_cr"]
            print(f"{s['run_id']:<24}  {'ok':<22}  {Lp:>6.2f}  {sig:>6.3f}  {b0:>10.3g}  {cr:>10.3g}")
        else:
            print(f"{s['run_id']:<24}  {s['status']:<22}  {'-':>6}  {'-':>6}  {'-':>10}  {'-':>10}")


if __name__ == "__main__":
    main()
