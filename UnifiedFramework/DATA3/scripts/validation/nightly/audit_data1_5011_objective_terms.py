#!/usr/bin/env python3
"""Audit legacy vs unified DATA1 objective terms for dataset 501.1."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from scipy.io import loadmat

REPO_ROOT = Path(__file__).resolve().parents[5]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import utility  # noqa: E402
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode import unified_codebase_runfile as runfile  # noqa: E402
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library import (  # noqa: E402
    ExperimentMode,
    ModelOptions,
    ParameterGuess,
    RunMode,
    _build_measurement_records_by_type,
    _legacy_vial_gate_config,
    build_guess_from_experiment_v24,
    load_experiment_easy,
    model_construct_for_legacy_fit_v24,
)
from UnifiedFramework.DATA3.tests.regression.test_unified_codebase_pytest_validation import _run_data1_pipeline  # noqa: E402

DATASET = "501.1"
DATASET_PATH = runfile.PRESET_DATA1_PATHS[0]
FIT_PATH = REPO_ROOT / "DATA1_matlab" / "data" / "501.1 concpolar" / "fit_stru.mat"
OUT_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "pytest_validation"


def _load_legacy_theta() -> dict[str, float]:
    fit_stru = loadmat(FIT_PATH, squeeze_me=True, struct_as_record=False)["fit_stru"]
    return {
        "Lp": float(fit_stru.Lp),
        "B": float(fit_stru.B),
        "sigma": float(fit_stru.sigma),
    }


def _load_unified_theta() -> dict[str, float]:
    result = _run_data1_pipeline(DATASET_PATH)
    return {k: float(v) for k, v in result["parmest"]["theta"].to_dict().items()}


def _legacy_eval(theta: dict[str, float]) -> dict[str, float | int | str]:
    data_stru = utility.loadmat(str(DATASET_PATH))["data_stru"]
    fit_stru, _sim_stru, _sim_inter = utility.solve_model(
        data_stru,
        "DATA",
        theta=theta,
        sim_opt=True,
        B_form="single",
        LOUD=False,
    )
    return {
        "implementation": "legacy_utility",
        "nfe": 300,
        "sigma_requested": float(theta["sigma"]),
        "sigma_applied": float(theta["sigma"]),
        "objective_total": float(fit_stru["Obj"]),
        "objective_mass": float(fit_stru["obj_m"]) * 1e4,
        "objective_cv": float(fit_stru["obj_cv"]) * 1e4,
        "objective_cf": float(fit_stru["obj_cr"]) * 1e4,
        "n_mass_terms": int(len(fit_stru["res_std"]["res_m"])),
        "n_cv_terms": int(len(fit_stru["res_std"]["res_cp"])),
        "n_cf_terms": int(len(fit_stru["res_std"]["res_cf"])),
    }


def _unified_eval(theta: dict[str, float], *, nfe: int) -> dict[str, float | int | str]:
    exp, _validation = load_experiment_easy(
        str(DATASET_PATH),
        selector=None,
        convert_to_concentration=False,
        conductivity_model="msa",
        plot=False,
    )
    options = ModelOptions(
        mode=ExperimentMode.DATA,
        run_mode=RunMode.ESTIMATION,
        b_form="single",
        nfe=nfe,
        use_sigma_logit_transform=False,
        use_multistart_for_mat_legacy=False,
    )
    guess = build_guess_from_experiment_v24(
        exp,
        options,
        override=ParameterGuess(Lp=theta["Lp"], B=theta["B"], sigma=theta["sigma"]),
    )
    m = model_construct_for_legacy_fit_v24(exp, options=options, guess=guess)
    sigma_lb, sigma_ub = m.sigma.bounds
    sigma_applied = float(min(max(theta["sigma"], sigma_lb), sigma_ub))
    m.Lp.fix(theta["Lp"])
    m.B.fix(theta["B"])
    m.sigma.fix(sigma_applied)
    solver = pyo.SolverFactory("ipopt")
    raw = solver.solve(m, tee=False)

    terms = _build_measurement_records_by_type(
        m,
        exp,
        sigma_mass_g=0.01,
        sigma_cV_rel=0.03,
        sigma_cF_rel=0.003,
    )
    mass_vals = [float(pyo.value(t)) for t in terms["mass"]]
    cv_vals = [float(pyo.value(t)) for t in terms["cV"]]
    cf_vals = [float(pyo.value(t)) for t in terms["cF"]]
    _n_v0, n_extra = _legacy_vial_gate_config(exp)
    collect_vial = max(1, len(exp.vials) - n_extra)
    expr_mass = sum(mass_vals) / max(1, len(mass_vals))
    expr_cv = sum(cv_vals) / collect_vial
    expr_cf = sum(cf_vals) / max(1, len(cf_vals))

    return {
        "implementation": "unified_v24",
        "nfe": int(nfe),
        "solver_status": str(raw.solver.termination_condition),
        "sigma_requested": float(theta["sigma"]),
        "sigma_applied": sigma_applied,
        "objective_total": float(pyo.value(m.Total_Cost_Objective)),
        "objective_mass": float(1e4 * expr_mass),
        "objective_cv": float(1e4 * expr_cv),
        "objective_cf": float(1e4 * expr_cf),
        "n_mass_terms": int(len(mass_vals)),
        "n_cv_terms": int(len(cv_vals)),
        "n_cf_terms": int(len(cf_vals)),
    }


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    legacy_theta = _load_legacy_theta()
    unified_theta = _load_unified_theta()

    rows: list[dict[str, object]] = []
    for theta_name, theta in [("legacy_theta", legacy_theta), ("unified_theta", unified_theta)]:
        legacy_row = _legacy_eval(theta)
        legacy_row["theta_source"] = theta_name
        legacy_row["Lp"] = theta["Lp"]
        legacy_row["B"] = theta["B"]
        legacy_row["sigma"] = theta["sigma"]
        rows.append(legacy_row)

        for nfe in (30, 300):
            unified_row = _unified_eval(theta, nfe=nfe)
            unified_row["theta_source"] = theta_name
            unified_row["Lp"] = theta["Lp"]
            unified_row["B"] = theta["B"]
            unified_row["sigma"] = theta["sigma"]
            rows.append(unified_row)

    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_ROOT / "data1_5011_objective_audit.csv", index=False)

    pivot = summary.pivot_table(
        index=["theta_source", "Lp", "B", "sigma"],
        columns=["implementation", "nfe"],
        values=["objective_total", "objective_mass", "objective_cv", "objective_cf"],
        aggfunc="first",
    )
    pivot.to_csv(OUT_ROOT / "data1_5011_objective_audit_pivot.csv")

    note = {
        "dataset": DATASET,
        "legacy_theta": legacy_theta,
        "unified_theta": unified_theta,
        "artifacts": [
            "UnifiedFramework/DATA3/results/pytest_validation/data1_5011_objective_audit.csv",
            "UnifiedFramework/DATA3/results/pytest_validation/data1_5011_objective_audit_pivot.csv",
        ],
    }
    (OUT_ROOT / "data1_5011_objective_audit_note.json").write_text(
        json.dumps(note, indent=2),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
