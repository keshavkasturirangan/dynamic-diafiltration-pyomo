"""Phase 1.5: Warm-anchor rollout for NaCl sheets WITHOUT cF-floor contour panels.

Scope:
  - MC5.07.23.24_NaCl   (concentration regime; same as MC2.NaCl)
  - MC5.07.23.24_SNaCl  (concentration regime; replicate of MC4.SNaCl)

Both sheets lack cF-floor contour panels (only wall-geometry ones exist
on disk), so we cannot use the contour-seeded strategy. Instead we anchor
on the MC2.05.07.24_NaCl converged θ (Lp=9.11, B=8.31, σ=0.66) and run
K=4 trials per sheet: 1 deterministic seed + 3 LHS-jitter samples in a
±15-20 % box around the seed. Each trial is a standalone IPOPT call with
multistart=False.

Per-sheet wall-cap: 600s (covers 4 trials × ~90s/trial + overhead).
Trials run sequentially inside the worker (one IPOPT child at a time).
killpg-enforced timeout reaps the worker + its IPOPT child if anything
hangs.
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

SCRATCH = Path("/tmp/data3_warm_anchor")
SCRATCH.mkdir(parents=True, exist_ok=True)

WORKER_SCRIPT = ROOT / "_warm_anchor_worker.py"

WARM_ANCHOR_SHEETS = [
    "MC5.07.23.24_NaCl",
    "MC5.07.23.24_SNaCl",
]

PER_SHEET_TIMEOUT_S = 600   # 4 trials × 90s cap + overhead
K_TRIALS_PER_SHEET = 4      # 1 seed + 3 LHS jitter


def run_one(run_id, fam, timeout_s, k_trials):
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
        str(SAVE_DIR), str(k_trials),
    ]
    t0 = time.time()
    proc = subprocess.Popen(
        cmd,
        start_new_session=True,
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
    import refactored_ucb_library as lib

    print(f"[warm-anchor] Phase 1.5 rollout — NaCl sheets WITHOUT cF-floor contour panels")
    print(f"[warm-anchor] Seed source: MC2.05.07.24_NaCl converged θ "
          f"(Lp=9.11, B=8.31, σ=0.66)")
    print(f"[warm-anchor] cF-floor ON (1 mM), contour seeds OFF, multistart in-worker")
    print(f"[warm-anchor] K trials per sheet: {K_TRIALS_PER_SHEET} "
          f"(1 deterministic + {K_TRIALS_PER_SHEET - 1} LHS-jitter)")
    print(f"[warm-anchor] Per-sheet hard cap: {PER_SHEET_TIMEOUT_S}s (killpg-enforced)")
    print(f"[warm-anchor] save_dir = {SAVE_DIR}")
    print()

    t0 = time.time()
    summary = []
    for rid in WARM_ANCHOR_SHEETS:
        fam = lib.NF270_RUN_REGISTRY.get(rid)
        if fam is None:
            print(f"  [{rid}] SKIP — not in registry")
            summary.append({"run_id": rid, "status": "skip-not-in-registry"})
            continue
        print(f"  [{rid}]  starting subprocess (cap {PER_SHEET_TIMEOUT_S}s, process group)...")
        out = run_one(rid, fam, PER_SHEET_TIMEOUT_S, K_TRIALS_PER_SHEET)
        if "error" in out and out.get("status") != "ok":
            tail = ""
            if "stdout_tail" in out:
                tail = "\n      stdout_tail: " + out["stdout_tail"].replace("\n", "\n      ")
            print(f"    {out['error']}  ({out.get('fit_time_s', out.get('total_time_s', '?'))}s){tail}")
            summary.append({"run_id": rid, "status": out["error"], "out": out})
            continue
        p = out.get("parameters", {})
        plot_n = len(out.get("plot_paths", []))
        plot_err = f" [plot:{out['plot_error']}]" if "plot_error" in out else ""
        n_ok_trials = sum(1 for t in out.get("trials", []) if "error" not in t)
        print(f"    BEST trial #{out['best_idx']} of {out['k_trials']} "
              f"({n_ok_trials} ok, {out['total_time_s']}s total)")
        print(f"      Lp={p.get('Lp', float('nan')):.3f}, "
              f"sigma={p.get('sigma', float('nan')):.3f}, "
              f"beta_0={p.get('beta_0', float('nan')):.3g}, "
              f"obj_m={out['obj_m']:.2g}/cv={out['obj_cv']:.2g}/cr={out['obj_cr']:.2g}, "
              f"{plot_n} plots{plot_err}")
        # Show all trials' σ landings to diagnose multimodality
        sigmas = []
        for t in out.get("trials", []):
            if "parameters" in t:
                sigmas.append(f"{t['parameters'].get('sigma', float('nan')):.3f}")
            else:
                sigmas.append("err")
        print(f"      σ landings across trials: [{', '.join(sigmas)}]")
        summary.append({"run_id": rid, "status": "ok", "out": out})

    print()
    print(f"[warm-anchor] === SUMMARY ({(time.time()-t0)/60:.1f} min) ===")
    n_ok = sum(1 for s in summary if s["status"] == "ok")
    print(f"[warm-anchor] {n_ok}/{len(summary)} sheets fit successfully\n")
    print(f"{'Run ID':<24}  {'status':<22}  {'Lp':>6}  {'sigma':>6}  {'beta_0':>10}  {'obj_cr':>10}  {'trials':>8}")
    print("-" * 100)
    for s in summary:
        if s["status"] == "ok":
            p = s["out"]["parameters"]
            Lp = p.get("Lp", float("nan"))
            sig = p.get("sigma", float("nan"))
            b0 = p.get("beta_0", float("nan"))
            cr = s["out"]["obj_cr"]
            n_ok_trials = sum(1 for t in s["out"].get("trials", []) if "error" not in t)
            print(f"{s['run_id']:<24}  {'ok':<22}  {Lp:>6.2f}  {sig:>6.3f}  {b0:>10.3g}  {cr:>10.3g}  {n_ok_trials:>3}/{s['out']['k_trials']:>3}")
        else:
            print(f"{s['run_id']:<24}  {s['status']:<22}  {'-':>6}  {'-':>6}  {'-':>10}  {'-':>10}  {'-':>8}")

    print()
    print(f"[warm-anchor] σ ≈ 0.66 across trials → MC2 anchor is the right basin for these replicates")
    print(f"[warm-anchor] σ scattered or pinned at 1.0 → multimodal / data genuinely prefers different basin")


if __name__ == "__main__":
    main()
