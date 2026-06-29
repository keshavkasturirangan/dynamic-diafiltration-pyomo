"""Capture structural fingerprints of the Pyomo models for DATA1/DATA2 cases.

Migrated from refactored_codes_v1/tests/_capture_references.py, with paths
re-pointed at the new home in pytest_refactored_codes_v1/.

For each reference case, we BUILD the Pyomo model via `model_construct_inter`
but DON'T solve it. We record a structural fingerprint:
   - data_config keys + values (loader regression)
   - per-vial array lengths (loader regression)
   - Pyomo model variable count + names
   - Pyomo model constraint count + names

The pytest then rebuilds the model with the same inputs and asserts every
captured field matches.

What this catches: Loader regressions (data_stru shape changes), model-build
regressions (new constraint added, variable bounds drift, DAE mesh changes),
any silent shape/structure drift in the Pyomo build path.

What this does NOT catch: Fit-result drift (would need IPOPT solve), ODE
solve correctness (would need a feasible forward simulation).

Usage
-----
    cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo
    python3 pytest_refactored_codes_v1/_capture_references.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent                  # pytest_refactored_codes_v1/
REPO_ROOT = HERE.parent
LIB_DIR = REPO_ROOT / "refactored_codes_v1"
HELPERS_DIR = HERE / "helpers"
BASELINES_DIR = HERE / "baselines"

for p in (str(LIB_DIR), str(HELPERS_DIR), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

os.chdir(str(REPO_ROOT))

import numpy as np  # noqa: E402
import refactored_ucb_library as lib  # noqa: E402
from reference_cases import ALL_CASES  # noqa: E402

OUT_JSON = BASELINES_DIR / "regression_references.json"


def _fingerprint_data_stru(data_stru):
    """Loader-side fingerprint: shape of what the loader produced."""
    cfg = data_stru.get("data_config", {})
    raw = data_stru.get("data_raw", [])

    config_scalars = {}
    for k, v in cfg.items():
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            try:
                config_scalars[k] = float(v)
            except (TypeError, ValueError):
                pass

    per_vial_shape = []
    for i, vial in enumerate(raw):
        rec = {"vial": i + 1}
        for key in ("time", "mass", "cF_exp", "cV_avg"):
            arr = vial.get(key)
            if arr is None:
                rec[f"{key}_len"] = 0
            else:
                try:
                    rec[f"{key}_len"] = int(np.asarray(arr).size)
                except Exception:
                    rec[f"{key}_len"] = -1
        per_vial_shape.append(rec)

    return {
        "n_vials": len(raw),
        "config_scalars": config_scalars,
        "per_vial_shape": per_vial_shape,
    }


def _fingerprint_pyomo_model(model):
    """Build-side fingerprint: Pyomo model variable + constraint counts."""
    from pyomo.environ import Var, Constraint
    var_names = sorted(c.name for c in model.component_objects(Var, descend_into=True))
    con_names = sorted(c.name for c in model.component_objects(Constraint, descend_into=True))
    n_var_components = len(var_names)
    n_con_components = len(con_names)
    n_atomic_vars = sum(1 for _ in model.component_data_objects(Var, descend_into=True))
    n_atomic_cons = sum(1 for _ in model.component_data_objects(Constraint, descend_into=True))
    return {
        "n_var_components": n_var_components,
        "n_con_components": n_con_components,
        "n_atomic_vars": n_atomic_vars,
        "n_atomic_cons": n_atomic_cons,
        "var_component_names": var_names,
        "con_component_names": con_names,
    }


def main():
    t0 = time.time()
    print(f"[capture] starting structural-fingerprint capture for {len(ALL_CASES)} cases")
    print(f"[capture] output → {OUT_JSON}")
    print(f"[capture] working dir: {os.getcwd()}")
    print()

    captured = {}
    for case in ALL_CASES:
        case_id = case["case_id"]
        data_file = REPO_ROOT / case["data_dir"] / case["data_file"]

        print(f"  [{case_id}]")

        skip_reason = case.get("skip_reason")
        if skip_reason:
            print(f"    SKIPPED (opt-out): {skip_reason}")
            captured[case_id] = {"skip_reason": skip_reason}
            continue

        if not data_file.exists():
            print(f"    SKIPPED — data file not found: {data_file}")
            captured[case_id] = {"error": "data file not found"}
            continue

        t_case = time.time()

        try:
            data_stru = lib.loadmat(str(data_file))["data_stru"]
        except Exception as exc:
            print(f"    LOAD ERROR: {exc!r}")
            captured[case_id] = {"error": f"load failed: {exc!r}"}
            continue
        loader_fp = _fingerprint_data_stru(data_stru)

        try:
            model = lib.model_construct_inter(
                data_stru,
                case["mode"],
                case["theta"],
                False,
                case["B_form"],
                workflow_family=case["workflow_family"],
            )
        except Exception as exc:
            print(f"    BUILD ERROR: {type(exc).__name__}: {exc}")
            captured[case_id] = {
                "error": f"model build failed: {type(exc).__name__}: {exc}",
                "loader_fingerprint": loader_fp,
            }
            continue
        build_fp = _fingerprint_pyomo_model(model)

        captured[case_id] = {
            "mode": case["mode"],
            "B_form": str(case["B_form"]),
            "workflow_family": case["workflow_family"],
            "seed_theta": case["theta"],
            "capture_time_s": round(time.time() - t_case, 2),
            "captured_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "loader_fingerprint": loader_fp,
            "build_fingerprint": build_fp,
        }
        print(
            f"    OK ({time.time() - t_case:.1f}s)  "
            f"{loader_fp['n_vials']} vials, "
            f"{build_fp['n_atomic_vars']} vars, "
            f"{build_fp['n_atomic_cons']} cons"
        )

    OUT_JSON.write_text(json.dumps(captured, indent=2, sort_keys=True))
    print()
    print(f"[capture] wrote {OUT_JSON} in {time.time()-t0:.1f}s total")
    n_ok = sum(1 for v in captured.values() if "build_fingerprint" in v)
    n_skip = sum(1 for v in captured.values() if "skip_reason" in v)
    n_err = sum(1 for v in captured.values() if "error" in v)
    print(f"[capture] captured={n_ok}, skipped={n_skip}, errors={n_err}")


if __name__ == "__main__":
    main()
