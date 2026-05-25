"""Shared helper: rebuild the Pyomo model and compare structural fingerprints.

Used by both test_data1_regression.py and test_data2_regression.py.

Compared fields (all EXACT match, no tolerance):
  - data_config scalar values  (allowed within 1% — these are real numbers)
  - per-vial array lengths     (must match exactly)
  - Pyomo Var/Constraint counts (must match exactly)
  - Pyomo Var/Constraint names  (must match exactly as a set)

Any drift in these signals a real structural change to the loader or the
model build path — exactly the regressions we want to catch.
"""
from __future__ import annotations

import pytest

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from _capture_references import _fingerprint_data_stru, _fingerprint_pyomo_model
from conftest import REL_TOLERANCE


def run_and_compare(lib, case, ref, data_path):
    """Rebuild model for `case`, compare loader + build fingerprints to `ref`."""
    # Build fresh fingerprints from current code state.
    data_stru = lib.loadmat(str(data_path))["data_stru"]
    cur_loader_fp = _fingerprint_data_stru(data_stru)

    try:
        model = lib.model_construct_inter(
            data_stru,
            case["mode"],
            case["theta"],
            False,                  # sim_opt=False
            case["B_form"],
            workflow_family=case["workflow_family"],
        )
    except Exception as exc:
        pytest.fail(
            f"model_construct_inter raised {type(exc).__name__}: {exc}\n"
            f"(Reference capture succeeded; current code build failed — "
            f"a regression has been introduced.)"
        )
    cur_build_fp = _fingerprint_pyomo_model(model)

    ref_loader_fp = ref["loader_fingerprint"]
    ref_build_fp = ref["build_fingerprint"]

    mismatches = []

    # --- Loader side: vial count + per-vial array shapes ---
    if cur_loader_fp["n_vials"] != ref_loader_fp["n_vials"]:
        mismatches.append(
            f"  n_vials: now={cur_loader_fp['n_vials']}, "
            f"ref={ref_loader_fp['n_vials']}"
        )

    # data_config scalars within 1% (these are physical values that may
    # legitimately fluctuate at the last decimal due to floating-point).
    for k, ref_v in ref_loader_fp["config_scalars"].items():
        cur_v = cur_loader_fp["config_scalars"].get(k)
        if cur_v is None:
            mismatches.append(f"  config_scalars.{k}: now=<missing>, ref={ref_v}")
            continue
        denom = max(abs(ref_v), 1e-9)
        rel = abs(cur_v - ref_v) / denom
        if rel > REL_TOLERANCE:
            mismatches.append(
                f"  config_scalars.{k}: now={cur_v:.6g}, ref={ref_v:.6g}, "
                f"rel_diff={rel:.3%} (tol={REL_TOLERANCE:.1%})"
            )

    # Per-vial array lengths (exact match).
    for cur_row, ref_row in zip(cur_loader_fp["per_vial_shape"],
                                 ref_loader_fp["per_vial_shape"]):
        for key in ("time_len", "mass_len", "cF_exp_len", "cV_avg_len"):
            if cur_row.get(key) != ref_row.get(key):
                mismatches.append(
                    f"  vial {ref_row.get('vial')} {key}: "
                    f"now={cur_row.get(key)}, ref={ref_row.get(key)}"
                )

    # --- Build side: variable + constraint counts (exact match) ---
    for key in ("n_var_components", "n_con_components",
                "n_atomic_vars", "n_atomic_cons"):
        if cur_build_fp[key] != ref_build_fp[key]:
            mismatches.append(
                f"  {key}: now={cur_build_fp[key]}, ref={ref_build_fp[key]}"
            )

    # Variable & constraint names (exact set match).
    cur_vars = set(cur_build_fp["var_component_names"])
    ref_vars = set(ref_build_fp["var_component_names"])
    added_v = sorted(cur_vars - ref_vars)
    removed_v = sorted(ref_vars - cur_vars)
    if added_v:
        mismatches.append(f"  Vars added: {added_v}")
    if removed_v:
        mismatches.append(f"  Vars removed: {removed_v}")

    cur_cons = set(cur_build_fp["con_component_names"])
    ref_cons = set(ref_build_fp["con_component_names"])
    added_c = sorted(cur_cons - ref_cons)
    removed_c = sorted(ref_cons - cur_cons)
    if added_c:
        mismatches.append(f"  Constraints added: {added_c}")
    if removed_c:
        mismatches.append(f"  Constraints removed: {removed_c}")

    if mismatches:
        pytest.fail(
            "Structural regression detected:\n" + "\n".join(mismatches[:30])
            + (f"\n... ({len(mismatches)} total)" if len(mismatches) > 30 else "")
        )
