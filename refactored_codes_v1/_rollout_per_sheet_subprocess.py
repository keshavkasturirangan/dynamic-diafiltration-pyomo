"""Per-sheet rollout with REAL OS-level timeouts (subprocess-based).

Why subprocess instead of signal.alarm
--------------------------------------
signal.alarm queues a SIGALRM to the Python process. When Python is
blocked inside a compiled-C extension (Pyomo → CasADi → IPOPT → MUMPS),
Python cannot dispatch the signal handler until control returns to
Python bytecode. IPOPT never returns if MUMPS is stuck inside a single
linear-solver factorization. The alarm never fires.

This rollout works around that by running each sheet's fit as a
SEPARATE OS PROCESS via subprocess.run(timeout=N). When the timeout
expires, subprocess.run sends SIGTERM (then SIGKILL) to the entire
child process tree — including IPOPT — regardless of what state any
of them are in.

Per-sheet outcome is captured in a JSON written by the worker, then
read back by this parent. The plots get rendered in the parent from
the JSON (or skipped if the fit failed).
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import refactored_ucb_library as lib  # noqa: E402

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

PER_SHEET_TIMEOUT_S = 120
WORKER_SCRIPT = ROOT / "_b_form_worker.py"     # generated below if missing
SCRATCH = Path("/tmp/data3_rollout_scratch")
SCRATCH.mkdir(parents=True, exist_ok=True)


WORKER_SOURCE = """
\"\"\"One-sheet fit worker. Invoked as a subprocess so the parent can
hard-kill it via OS timeout regardless of what IPOPT is doing.\"\"\"
import json, os, sys, time, traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
os.chdir(str(ROOT.parent))

import refactored_ucb_library as lib

def main():
    run_id  = sys.argv[1]
    wb_path = sys.argv[2]
    sheet   = sys.argv[3]
    out_json = sys.argv[4]

    out = {"run_id": run_id}
    t0 = time.time()
    try:
        ds = lib.loadxlsx(wb_path, sheet=sheet)["data_stru"]
    except Exception as exc:
        out["error"] = f"load: {type(exc).__name__}: {exc}"
        Path(out_json).write_text(json.dumps(out))
        return

    out["mode"]    = ds.get("mode")
    out["n_vials"] = ds["data_config"]["n"]
    out["M_F0"]    = float(ds["data_config"]["M_F0"])
    out["M_O"]     = float(ds["data_config"]["M_O"])

    try:
        fit_stru, sim_stru, _ = lib.solve_model(
            ds, ds["mode"], None,        # theta=None → per-sheet seeds
            sim_opt=False,
            B_form=1,                    # DATA2 recipe
            LOUD=False,
            workflow_family="DATA3",
            solver_max_cpu_time=90,
        )
        out["fit_time_s"] = round(time.time() - t0, 1)
        params = fit_stru.get("parameters", {})
        out["parameters"] = {
            k: float(v) for k, v in params.items()
            if isinstance(v, (int, float))
        }
        out["obj_m"]  = float(fit_stru.get("obj_m",  0))
        out["obj_cv"] = float(fit_stru.get("obj_cv", 0))
        out["obj_cr"] = float(fit_stru.get("obj_cr", 0))
        # Also dump sim_stru so the parent can render plots without re-fitting
        out["sim_stru"] = [
            {
                "time": list(map(float, v.get("time", []))),
                "mF":   list(map(float, v.get("mF",   []))) if v.get("mF") is not None else None,
                "mV":   list(map(float, v.get("mV",   []))),
                "cF":   list(map(float, v.get("cF",   []))),
                "cH":   list(map(float, v.get("cH",   []))),
                "cV":   list(map(float, v.get("cV",   []))),
            }
            for v in sim_stru
        ]
    except Exception as exc:
        out["error"] = f"fit: {type(exc).__name__}: {exc}"
        out["fit_time_s"] = round(time.time() - t0, 1)
        out["traceback"] = traceback.format_exc()

    Path(out_json).write_text(json.dumps(out))

if __name__ == "__main__":
    main()
"""

def ensure_worker():
    if not WORKER_SCRIPT.exists():
        WORKER_SCRIPT.write_text(WORKER_SOURCE)
    return WORKER_SCRIPT


def run_one(run_id, fam, timeout_s):
    """Run one sheet's fit in a subprocess. Returns the worker's JSON or an
    error dict."""
    wb_path = NF270_ROOT / fam["workbook"]
    if not wb_path.exists():
        return {"run_id": run_id, "error": f"workbook missing: {wb_path}"}
    out_json = SCRATCH / f"{run_id.replace('/', '_')}.json"
    if out_json.exists():
        out_json.unlink()
    cmd = [
        sys.executable, "-u", str(WORKER_SCRIPT),
        run_id, str(wb_path), fam["sheet"], str(out_json),
    ]
    t0 = time.time()
    try:
        completed = subprocess.run(
            cmd,
            timeout=timeout_s,
            capture_output=True,
            text=True,
            check=False,
        )
        elapsed = time.time() - t0
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        return {
            "run_id": run_id,
            "error": "wall-timeout",
            "fit_time_s": round(elapsed, 1),
            "timeout_s": timeout_s,
        }
    if not out_json.exists():
        # Worker died before writing
        return {
            "run_id": run_id,
            "error": f"worker died (rc={completed.returncode})",
            "stderr_tail": (completed.stderr or "")[-500:],
            "fit_time_s": round(elapsed, 1),
        }
    try:
        return json.loads(out_json.read_text())
    except Exception as exc:
        return {
            "run_id": run_id,
            "error": f"could not parse worker JSON: {exc!r}",
        }


def render_plots_from_json(out):
    """If the worker reported a fit, build a data_stru dict and render
    paper-style plots via run_data3_time_series_plots."""
    if "error" in out or "sim_stru" not in out:
        return []
    fam = lib.NF270_RUN_REGISTRY.get(out["run_id"])
    if fam is None:
        return []
    try:
        wrapped = lib.loadxlsx(NF270_ROOT / fam["workbook"], sheet=fam["sheet"])
        ds = wrapped["data_stru"]
    except Exception:
        return []
    # Reconstruct sim_stru shape that run_data3_time_series_plots expects.
    sim_stru = out["sim_stru"]
    results_dict = {
        "data": ds,
        "sim_stru": sim_stru,
        "model_settings": {"workflow_family": "DATA3"},
    }
    try:
        return lib.run_data3_time_series_plots(results_dict, save_dir=SAVE_DIR, show=False)
    except Exception as exc:
        print(f"      PLOT ERROR: {type(exc).__name__}: {exc}")
        return []


def main():
    ensure_worker()
    t0 = time.time()
    print(f"[rollout-sub] starting subprocess-based rollout for {len(SINGLE_SALT_RUN_IDS)} sheets")
    print(f"[rollout-sub] timeout: {PER_SHEET_TIMEOUT_S}s per sheet (OS-level — subprocess.terminate)")
    print()

    summary = []
    for rid in SINGLE_SALT_RUN_IDS:
        fam = lib.NF270_RUN_REGISTRY.get(rid)
        if fam is None:
            print(f"  [{rid}] SKIP — not in registry")
            summary.append({"run_id": rid, "status": "skip-not-in-registry"})
            continue
        print(f"  [{rid}]  starting subprocess (cap {PER_SHEET_TIMEOUT_S}s)...")
        out = run_one(rid, fam, PER_SHEET_TIMEOUT_S)
        if "error" in out:
            print(f"    {out['error']}  ({out.get('fit_time_s','?')}s)")
            summary.append({"run_id": rid, "status": out["error"], "out": out})
            continue
        p = out.get("parameters", {})
        print(f"    FIT OK ({out['fit_time_s']}s)  "
              f"Lp={p.get('Lp', float('nan')):.3f}, "
              f"sigma={p.get('sigma', float('nan')):.3f}, "
              f"beta_0={p.get('beta_0', float('nan')):.3g}, "
              f"obj_m={out['obj_m']:.2g}/cv={out['obj_cv']:.2g}/cr={out['obj_cr']:.2g}")
        # Render plots in the parent (cheap; just plotting, no IPOPT)
        plot_paths = render_plots_from_json(out)
        if plot_paths:
            print(f"    PLOT OK ({len(plot_paths)} files)")
        summary.append({"run_id": rid, "status": "ok", "out": out})

    print()
    print(f"[rollout-sub] === SUMMARY ({(time.time()-t0)/60:.1f} min) ===")
    n_ok = sum(1 for s in summary if s["status"] == "ok")
    print(f"[rollout-sub] {n_ok}/{len(summary)} sheets fit successfully\n")
    print(f"{'Run ID':<24}  {'status':<22}  {'Lp':>6}  {'sigma':>6}  {'beta_0':>8}  {'beta_1':>8}")
    print("-" * 80)
    for s in summary:
        if s["status"] == "ok":
            p = s["out"]["parameters"]
            Lp = p.get("Lp", float("nan"))
            sig = p.get("sigma", float("nan"))
            b0 = p.get("beta_0", float("nan"))
            b1 = p.get("beta_1", float("nan"))
            print(f"{s['run_id']:<24}  {'ok':<22}  {Lp:>6.2f}  {sig:>6.3f}  {b0:>8.3g}  {b1:>8.3g}")
        else:
            print(f"{s['run_id']:<24}  {s['status']:<22}  {'-':>6}  {'-':>6}  {'-':>8}  {'-':>8}")


if __name__ == "__main__":
    main()
