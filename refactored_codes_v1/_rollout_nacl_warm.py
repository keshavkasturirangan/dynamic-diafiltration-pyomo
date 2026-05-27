"""Focused NaCl rollout — seed every NaCl sheet from MC2 NaCl's optimum.

Strategy
--------
MC2.05.07.24_NaCl converged tonight under cf-floor + contour-seeded
multistart to an INTERIOR θ:
    Lp = 9.11   B = 8.31   sigma = 0.66
(The first NaCl fit all session to land off the σ wall.)

This rollout uses that θ as the explicit seed for the OTHER 4 NaCl
sheets that haven't converged yet. Rationale: same salt, same membrane
family (NF270), same operator/lab, similar regime for 3 of 4 → MC2's
basin should be the right neighborhood.

Knob settings
-------------
* cf-floor:        ON (1.0 mM) — without it σ pins at the wall
* contour seeds:   OFF — the May 22 contours are σ=wall geometry,
                   incompatible with cf-floor's interior basin
* multistart:      0 LHS — we have a known-good interior seed, no
                   perturbation needed
* B_form:          1 (matches the converged fit)

Process control
---------------
Each sheet's fit runs in a SEPARATE PROCESS GROUP (start_new_session=True)
so we can kill the entire group (Python worker + IPOPT child) with
os.killpg(SIGKILL) on timeout. This solves the v10 orphan problem:
subprocess.terminate() only sent SIGTERM to the worker, leaving IPOPT
reparented to launchd. killpg propagates to every descendant.

Worst-case wall time: 4 sheets × 300s cap = 20 min.
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

import refactored_ucb_library as lib  # noqa: E402

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
SAVE_DIR = (
    ROOT.parent / "UnifiedFramework" / "DATA3" / "results"
    / "paper_artifacts" / "nf270" / "campaign_figures"
)
SAVE_DIR.mkdir(parents=True, exist_ok=True)

NF270_ROOT = ROOT.parent / "UnifiedFramework" / "ExperimentalDataFiles"
SCRATCH = Path("/tmp/data3_nacl_warm")
SCRATCH.mkdir(parents=True, exist_ok=True)

WORKER_SCRIPT = ROOT / "_nacl_warm_worker.py"

# MC2.05.07.24_NaCl converged optimum — the seed for every other NaCl sheet.
# For B_form=1 the model uses (Lp, beta_0, beta_1, sigma); B is computed
# from beta_0 internally. beta_0 default = B; beta_1 = 1.0 (DATA2 default).
NACL_SEED_THETA = {
    "Lp":     9.11,
    "B":      8.31,
    "beta_0": 8.31,
    "beta_1": 1.0,
    "sigma":  0.66,
    "S0":     0,
}

# 6 NaCl sheets total in the single-salt subset; 2 already converged this
# session (MC2.NaCl in v11, MC3.SNaCl in earlier smoke test). The 4 below
# are the ones still needing a converged fit.
NACL_TARGETS = [
    "MC4.07.11.24_SNaCl",      # concentration regime
    "MC5.07.23.24_NaCl",       # concentration regime
    "MC5.07.23.24_SNaCl",      # concentration regime
    "MC5.07.23.24_S2NaCl",     # dilution regime
]

PER_SHEET_TIMEOUT_S = 300   # 5 min hard cap; converged MC2.NaCl took ~25 min
                            # with 13 multistart restarts. With one good seed
                            # and no multistart, each fit should be <90 sec.

WORKER_SOURCE = """\\
\"\"\"Per-sheet NaCl warm-start fit worker.

Runs as a subprocess so the parent can kill the whole process group
(worker + IPOPT child) via os.killpg on timeout. Reads its config from
argv, writes result to a JSON file.
\"\"\"
import json
import os
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
os.chdir(str(ROOT.parent))

import refactored_ucb_library as lib

# Enable cf-floor; disable contour seeds (May 22 wall-geometry incompatible
# with cf-floor's interior basin).
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0
lib.NF270_MULTISTART_USE_CONTOUR_SEEDS = False

def main():
    run_id    = sys.argv[1]
    wb_path   = sys.argv[2]
    sheet     = sys.argv[3]
    theta_json = sys.argv[4]
    out_json  = sys.argv[5]

    theta = json.loads(theta_json)
    out = {"run_id": run_id, "seed_theta": theta}
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
            ds, ds["mode"], theta,
            sim_opt=False,
            B_form=1,
            LOUD=False,
            workflow_family="DATA3",
            multistart=False,       # known-good seed, no LHS perturbation
            solver_max_cpu_time=240,
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

# ----------------------------------------------------------------------------
# Process-group helpers (the actually-reliable timeout strategy)
# ----------------------------------------------------------------------------

def ensure_worker():
    if not WORKER_SCRIPT.exists() or WORKER_SCRIPT.read_text() != WORKER_SOURCE:
        WORKER_SCRIPT.write_text(WORKER_SOURCE)
    return WORKER_SCRIPT


def run_one_sheet(run_id, fam, theta, timeout_s):
    """Run one sheet's fit in a separate process group with hard OS-level cap."""
    wb_path = NF270_ROOT / fam["workbook"]
    if not wb_path.exists():
        return {"run_id": run_id, "error": f"workbook missing: {wb_path}"}
    out_json = SCRATCH / f"{run_id.replace('/', '_')}.json"
    if out_json.exists():
        out_json.unlink()

    cmd = [
        sys.executable, "-u", str(WORKER_SCRIPT),
        run_id, str(wb_path), fam["sheet"],
        json.dumps(theta), str(out_json),
    ]
    t0 = time.time()
    # start_new_session=True makes the child the LEADER of a new process group.
    # On timeout we can os.killpg(pgid, SIGKILL) to nuke the whole group —
    # IPOPT included — which subprocess.terminate cannot do.
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        _stdout, stderr = proc.communicate(timeout=timeout_s)
        elapsed = time.time() - t0
        rc = proc.returncode
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        # KILL THE WHOLE PROCESS GROUP — Python + IPOPT + anything they spawned
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            pass
        # Reap the zombie
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
            "error": f"worker died (rc={rc})",
            "stderr_tail": (stderr or "")[-400:],
            "fit_time_s": round(elapsed, 1),
        }
    return json.loads(out_json.read_text())


def render_plots_from_json(out):
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
    results_dict = {
        "data": ds,
        "sim_stru": out["sim_stru"],
        "model_settings": {"workflow_family": "DATA3"},
    }
    try:
        return lib.run_data3_time_series_plots(results_dict, save_dir=SAVE_DIR, show=False)
    except Exception as exc:
        print(f"      PLOT ERROR: {type(exc).__name__}: {exc}")
        return []


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    ensure_worker()
    t0 = time.time()
    print(f"[nacl-warm] Focused NaCl rollout — seed every sheet from MC2 NaCl optimum")
    print(f"[nacl-warm] Seed theta: {NACL_SEED_THETA}")
    print(f"[nacl-warm] cf-floor ON, contour seeds OFF, multistart=False")
    print(f"[nacl-warm] Per-sheet hard cap: {PER_SHEET_TIMEOUT_S}s (killpg-enforced)")
    print(f"[nacl-warm] Targets: {len(NACL_TARGETS)} NaCl sheets ({', '.join(NACL_TARGETS)})")
    print()

    summary = []
    for rid in NACL_TARGETS:
        fam = lib.NF270_RUN_REGISTRY.get(rid)
        if fam is None:
            print(f"  [{rid}] SKIP — not in registry")
            summary.append({"run_id": rid, "status": "skip-not-in-registry"})
            continue

        print(f"  [{rid}]  starting subprocess (cap {PER_SHEET_TIMEOUT_S}s, process group)...")
        out = run_one_sheet(rid, fam, NACL_SEED_THETA, PER_SHEET_TIMEOUT_S)

        if "error" in out:
            print(f"    {out['error']}  ({out.get('fit_time_s','?')}s)")
            if "stderr_tail" in out:
                print(f"    stderr: {out['stderr_tail'][:200]}")
            summary.append({"run_id": rid, "status": out["error"], "out": out})
            continue

        p = out.get("parameters", {})
        print(f"    FIT OK ({out['fit_time_s']}s)  "
              f"Lp={p.get('Lp', float('nan')):.3f}, "
              f"sigma={p.get('sigma', float('nan')):.3f}, "
              f"beta_0={p.get('beta_0', float('nan')):.3g}, "
              f"obj_m={out['obj_m']:.2g}/cv={out['obj_cv']:.2g}/cr={out['obj_cr']:.2g}")

        plot_paths = render_plots_from_json(out)
        if plot_paths:
            print(f"    PLOT OK ({len(plot_paths)} files)")
        summary.append({"run_id": rid, "status": "ok", "out": out})

    print()
    print(f"[nacl-warm] === SUMMARY ({(time.time()-t0)/60:.1f} min) ===")
    n_ok = sum(1 for s in summary if s["status"] == "ok")
    print(f"[nacl-warm] {n_ok}/{len(summary)} NaCl sheets converged with MC2-NaCl seed")
    print()
    print(f"{'Run ID':<24}  {'status':<22}  {'Lp':>7}  {'B':>7}  {'sigma':>7}  {'beta_0':>8}  {'beta_1':>8}")
    print("-" * 90)
    for s in summary:
        if s["status"] == "ok":
            p = s["out"]["parameters"]
            Lp = p.get("Lp", float("nan"))
            B  = p.get("B",  float("nan"))
            sg = p.get("sigma", float("nan"))
            b0 = p.get("beta_0", float("nan"))
            b1 = p.get("beta_1", float("nan"))
            print(f"{s['run_id']:<24}  {'ok':<22}  {Lp:>7.3f}  {B:>7.3f}  {sg:>7.3f}  {b0:>8.3g}  {b1:>8.3g}")
        else:
            print(f"{s['run_id']:<24}  {s['status']:<22}  {'-':>7}  {'-':>7}  {'-':>7}  {'-':>8}  {'-':>8}")
    print()
    print(f"For comparison — MC2.05.07.24_NaCl (seed source):")
    p = NACL_SEED_THETA
    print(f"{'MC2.05.07.24_NaCl':<24}  {'ok (v11)':<22}  {p['Lp']:>7.3f}  {p['B']:>7.3f}  {p['sigma']:>7.3f}  {p['beta_0']:>8.3g}  {p['beta_1']:>8.3g}")


if __name__ == "__main__":
    main()
