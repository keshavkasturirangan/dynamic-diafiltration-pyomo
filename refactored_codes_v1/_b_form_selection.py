"""B_form model selection per salt for the NF270 single-salt rollout.

For each of 4 representative sheets (one per salt regime), fits 4 B_form
variants and picks the winner by AIC.

Run me from the repo root:
    python refactored_codes_v1/_b_form_selection.py

Outputs:
    refactored_codes_v1/B_FORM_SELECTION.md       — human-readable report
    refactored_codes_v1/B_FORM_SELECTION.json     — machine-readable results

Total wall time: ~45-80 min (16 fits, each capped at 5 min of IPOPT time).
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

os.chdir(str(REPO_ROOT))

import numpy as np  # noqa: E402
import refactored_ucb_library as lib  # noqa: E402

# Match the rollout's runtime flags: contour-seeded multistart + cf-floor.
lib.NF270_MULTISTART_USE_CONTOUR_SEEDS = True
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0

OUT_MD = HERE / "B_FORM_SELECTION.md"
OUT_JSON = HERE / "B_FORM_SELECTION.json"

# 4 representative sheets — one per salt × regime cell that matters.
REP_SHEETS = [
    {"run_id": "MC3.07.22.24_SNaCl",  "salt": "NaCl",  "regime": "dilution"},
    {"run_id": "MC5.07.23.24_NaCl",   "salt": "NaCl",  "regime": "concentration"},
    {"run_id": "MC2.05.07.24_CaCl2",  "salt": "CaCl2", "regime": "concentration"},
    {"run_id": "MC4.07.11.24_SLaCl3", "salt": "LaCl3", "regime": "concentration"},
]

# 4 B_form variants — each requires a different theta shape.
# We provide the smallest sensible seed; the contour-seed pipeline (and
# multistart with box_frac perturbation) will explore from there.
B_FORMS = [
    # NOTE: every seed theta must include S0 because Lag mode reads it as a
    # Param (fixed during the fit). When theta=None the loader supplies
    # S0=0 by default, but we're passing explicit theta dicts here so we
    # must provide it ourselves. S0=0 = "no continuous inflow during the
    # active leg" — works for both DATA3 NaCl dilution and concentration
    # sheets, where M_O captures the net mass change instead.
    {
        "label":  "single",
        "B_form": "single",
        "seed_theta": {"Lp": 5.0, "B": 1.0, "sigma": 0.5, "S0": 0},
        # 3 fitted scalars: Lp, B, sigma. S is also fit in Lag mode → +1.
    },
    {
        "label":  "pervial",
        "B_form": "pervial",
        # For pervial, B is a vector — solve_model treats B0 from data_stru
        # as the initial guess; seed_theta supplies the scalars.
        "seed_theta": {"Lp": 5.0, "B": 1.0, "sigma": 0.5, "S0": 0},
        # k = 2 + n_vials
    },
    {
        "label":  "convection",
        "B_form": "convection",
        "seed_theta": {"Lp": 5.0, "beta_0": 2.0, "beta_1": 0.5, "sigma": 0.5, "S0": 0},
        # 4 fitted scalars: Lp, beta_0, beta_1, sigma.
    },
    {
        "label":  "power-law (B_form=1)",
        "B_form": 1,
        "seed_theta": {"Lp": 5.0, "beta_0": 1.0, "beta_1": 0.5, "sigma": 0.5, "S0": 0},
        # 4 fitted scalars: Lp, beta_0, beta_1, sigma (B is reconstructed).
    },
]


def _load_sheet(run_id):
    family = lib.NF270_RUN_REGISTRY.get(run_id)
    if family is None:
        return None, f"run_id not in registry"
    wb = REPO_ROOT / "UnifiedFramework" / "ExperimentalDataFiles" / family["workbook"]
    if not wb.exists():
        return None, f"workbook missing: {wb}"
    try:
        return lib.loadxlsx(wb, sheet=family["sheet"])["data_stru"], None
    except Exception as exc:
        return None, f"load failed: {exc!r}"


def _count_parameters(b_form_label, fit_stru, data_stru):
    """How many scalar parameters did this fit actually estimate?"""
    params = fit_stru.get("parameters", {}) if isinstance(fit_stru, dict) else {}
    if not params:
        return 0
    if b_form_label == "single":
        return 3                              # Lp, B, sigma
    if b_form_label == "pervial":
        n_vials = data_stru.get("data_config", {}).get("n", 0)
        return 2 + int(n_vials)               # Lp, sigma, B[1..N]
    if b_form_label in ("convection", "power-law (B_form=1)"):
        return 4                              # Lp, beta_0, beta_1, sigma
    # Default: count scalar entries in the parameters dict.
    return sum(1 for v in params.values() if isinstance(v, (int, float)))


def _wsse_total(fit_stru):
    """Sum the three channel objectives — the WSSE the optimizer minimized."""
    if not isinstance(fit_stru, dict):
        return float("nan"), 0
    obj_m  = float(fit_stru.get("obj_m",  0.0))
    obj_cv = float(fit_stru.get("obj_cv", 0.0))
    obj_cr = float(fit_stru.get("obj_cr", 0.0))
    total = obj_m + obj_cv + obj_cr
    # N is the total number of residuals contributing. Use the counts the
    # fitter reports if available; otherwise estimate from data shape.
    count_m  = int(fit_stru.get("count_m",  0))
    count_cv = int(fit_stru.get("count_cv", 0))
    count_cr = int(fit_stru.get("count_cr", 0))
    n = count_m + count_cv + count_cr
    return total, n


def _aic(wsse, n, k):
    """AIC for least-squares: 2k + n·ln(WSSE/n)."""
    if n <= 0 or not math.isfinite(wsse) or wsse <= 0:
        return float("nan")
    return 2 * k + n * math.log(wsse / n)


def main():
    t0 = time.time()
    print(f"[b_form_select] starting; {len(REP_SHEETS)} rep sheets × {len(B_FORMS)} B_forms = {len(REP_SHEETS)*len(B_FORMS)} fits")
    print(f"[b_form_select] cf-floor: {lib.NF270_CF_RESIDUAL_FLOOR_MM} mM")
    print(f"[b_form_select] contour-seeded multistart: {lib.NF270_MULTISTART_USE_CONTOUR_SEEDS}")
    print()

    results = {}  # results[run_id][b_form_label] = {...}
    for sheet in REP_SHEETS:
        rid = sheet["run_id"]
        results[rid] = {"salt": sheet["salt"], "regime": sheet["regime"], "fits": {}}

        print(f"=== {rid} ({sheet['salt']} {sheet['regime']}) ===")
        data_stru, err = _load_sheet(rid)
        if err:
            results[rid]["error"] = err
            print(f"  LOAD ERROR: {err}")
            continue

        for bf in B_FORMS:
            label = bf["label"]
            print(f"  --- B_form = {label} ---")
            t_fit = time.time()
            try:
                fit_stru, _sim, _inter = lib.solve_model(
                    data_stru,
                    "Lag",
                    bf["seed_theta"],
                    sim_opt=False,
                    B_form=bf["B_form"],
                    LOUD=False,
                    workflow_family="DATA3",
                    solver_max_cpu_time=300,       # 5-min IPOPT cap per fit
                )
            except Exception as exc:
                dt = time.time() - t_fit
                print(f"    FIT ERROR ({dt:.1f}s): {type(exc).__name__}: {exc}")
                results[rid]["fits"][label] = {
                    "error": f"{type(exc).__name__}: {exc}",
                    "fit_time_s": round(dt, 1),
                }
                continue

            dt = time.time() - t_fit
            if not isinstance(fit_stru, dict) or "parameters" not in fit_stru:
                print(f"    FAILED ({dt:.1f}s): no parameters returned")
                results[rid]["fits"][label] = {
                    "error": "no parameters in fit_stru",
                    "fit_time_s": round(dt, 1),
                }
                continue

            wsse, n_data = _wsse_total(fit_stru)
            k = _count_parameters(label, fit_stru, data_stru)
            aic = _aic(wsse, n_data, k)

            # Flatten the parameters for the report.
            params_flat = {}
            for pk, pv in fit_stru.get("parameters", {}).items():
                try:
                    params_flat[pk] = float(pv)
                except (TypeError, ValueError):
                    if isinstance(pv, dict) and len(pv) <= 5:
                        params_flat[pk] = {kk: float(vv) for kk, vv in pv.items()
                                            if isinstance(vv, (int, float))}

            results[rid]["fits"][label] = {
                "wsse": float(wsse),
                "n_data": int(n_data),
                "k_params": int(k),
                "aic": float(aic),
                "fit_time_s": round(dt, 1),
                "parameters": params_flat,
            }
            print(f"    OK ({dt:.1f}s)  WSSE={wsse:.3g}, n={n_data}, k={k}, AIC={aic:.2f}")

        # Pick winner by lowest AIC.
        valid_fits = {
            label: f for label, f in results[rid]["fits"].items()
            if "aic" in f and math.isfinite(f["aic"])
        }
        if valid_fits:
            winner = min(valid_fits.items(), key=lambda kv: kv[1]["aic"])
            results[rid]["winner"] = winner[0]
            results[rid]["winner_aic"] = float(winner[1]["aic"])
            print(f"  → WINNER: {winner[0]} (AIC={winner[1]['aic']:.2f})")
        else:
            results[rid]["winner"] = None
            print(f"  → NO VALID WINNER (all fits failed)")
        print()

    OUT_JSON.write_text(json.dumps(results, indent=2, sort_keys=True))

    # Write the markdown report.
    md = ["# B_form selection results\n",
          f"_Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}_\n",
          f"Total wall time: {(time.time()-t0)/60:.1f} min\n",
          "## Summary table — winner per salt\n",
          "| Salt | Regime | Sheet | Winning B_form | AIC | k_params |",
          "|---|---|---|---|---|---|"]
    for sheet in REP_SHEETS:
        rid = sheet["run_id"]
        r = results.get(rid, {})
        if r.get("winner"):
            fit = r["fits"][r["winner"]]
            md.append(
                f"| {sheet['salt']} | {sheet['regime']} | `{rid}` | "
                f"**{r['winner']}** | {fit['aic']:.2f} | {fit['k_params']} |"
            )
        else:
            md.append(f"| {sheet['salt']} | {sheet['regime']} | `{rid}` | "
                      f"(no winner) | — | — |")

    md.append("\n## Per-sheet detail (all B_form fits)\n")
    for sheet in REP_SHEETS:
        rid = sheet["run_id"]
        r = results.get(rid, {})
        md.append(f"### `{rid}` — {sheet['salt']} {sheet['regime']}\n")
        md.append("| B_form | WSSE | n | k | AIC | fit time (s) | parameters |")
        md.append("|---|---|---|---|---|---|---|")
        for label in (bf["label"] for bf in B_FORMS):
            fit = r.get("fits", {}).get(label, {})
            if "error" in fit:
                md.append(f"| {label} | (error) | — | — | — | {fit.get('fit_time_s','—')} | {fit['error'][:60]} |")
            elif "aic" in fit:
                params_str = ", ".join(
                    f"{k}={v:.3g}" if isinstance(v, float) else f"{k}=…"
                    for k, v in list(fit.get("parameters", {}).items())[:5]
                )
                md.append(
                    f"| {label} | {fit['wsse']:.3g} | {fit['n_data']} | {fit['k_params']} | "
                    f"{fit['aic']:.2f} | {fit['fit_time_s']} | {params_str} |"
                )
            else:
                md.append(f"| {label} | (not run) | — | — | — | — | — |")
        if r.get("winner"):
            md.append(f"\n**Winner:** `{r['winner']}` (lowest AIC)\n")
        md.append("")

    OUT_MD.write_text("\n".join(md))

    print(f"[b_form_select] wrote {OUT_MD}")
    print(f"[b_form_select] wrote {OUT_JSON}")
    print(f"[b_form_select] total wall time: {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
