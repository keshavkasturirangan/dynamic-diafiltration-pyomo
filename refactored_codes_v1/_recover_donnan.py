#!/usr/bin/env python3
"""Recover Donnan DAE fits that railed in the campaign, using fixed-k retries each
run in an isolated subprocess with a hard wall-clock timeout (kills the whole
process group, including any orphaned ipopt, so a slow-failing k can't hang).

For each sheet whose result_donnan.json is missing or errored, tries k in
{0.6, 0.5, 0.45, 0.7, 0.35}; the first feasible fit within the budget wins.
Usage: python3 _recover_donnan.py [timeout_s=150]
"""
import sys, os, json, signal, subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _run_bform_campaign as camp

STUDY = camp.OUT
TIMEOUT = int(sys.argv[1]) if len(sys.argv) > 1 else 150
K_LADDER = [0.6, 0.5, 0.45, 0.7, 0.35]


def run_with_timeout(args, timeout):
    p = subprocess.Popen(args, start_new_session=True, cwd=str(HERE),
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    try:
        out, _ = p.communicate(timeout=timeout)
        return p.returncode == 0, out
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(p.pid), signal.SIGKILL)
        except Exception:
            pass
        p.wait()
        return False, "TIMEOUT"


def needs_recovery(rid):
    rp = STUDY / rid / "result_donnan.json"
    if not rp.exists():
        return True
    try:
        r = json.loads(rp.read_text())
        return ("error" in r) or ("parameters" not in r)
    except Exception:
        return True


def main():
    for rid, salt, stage in camp.SHEETS:
        if not needs_recovery(rid):
            print(f"[ok]    {rid} donnan already fit", flush=True)
            continue
        done = False
        for kfix in K_LADDER:
            ok, out = run_with_timeout(
                [sys.executable, "_donnan_fit_worker.py", rid, str(kfix), salt], TIMEOUT)
            # success iff the worker wrote a clean result with this k
            rp = STUDY / rid / "result_donnan.json"
            if rp.exists():
                try:
                    rr = json.loads(rp.read_text())
                    if "error" not in rr and rr.get("k_fixed") == kfix:
                        print(f"[recovered] {rid} k={kfix} WSSE3={rr.get('WSSE3'):.2f}", flush=True)
                        done = True
                        break
                except Exception:
                    pass
            tag = "timeout" if out == "TIMEOUT" else "infeasible"
            print(f"[try]   {rid} k={kfix} -> {tag}", flush=True)
        if not done:
            print(f"[FAIL]  {rid} donnan unrecovered (offline fit backs the curve)", flush=True)
    print("DONNAN RECOVERY DONE", flush=True)


if __name__ == "__main__":
    main()
