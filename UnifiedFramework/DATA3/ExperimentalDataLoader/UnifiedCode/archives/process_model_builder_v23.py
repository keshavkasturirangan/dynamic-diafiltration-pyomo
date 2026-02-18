#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Process-model builder (v23): unified Pyomo + ParmEst + DoE workflow.

This module keeps data classes light and moves heavy lifting into functions.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.contrib.parmest.experiment import Experiment as ParmestExperiment
from pyomo.contrib.parmest.parmest import Estimator
from pyomo.contrib.doe import DesignOfExperiments

from experiment_dataload_OOP_v23_conductivity_patched import ExperimentalData


class ExperimentMode(str, Enum):
    DATA = "DATA"
    LAG = "Lag"
    OVERFLOW = "Overflow"


class RunMode(str, Enum):
    SIMULATION = "SIMULATION"
    ESTIMATION = "ESTIMATION"


class BForm(str, Enum):
    SINGLE = "single"
    PERVIAL = "pervial"
    CONVECTION = "convection"


BFormType = Union[str, float]


@dataclass
class ModelOptions:
    mode: ExperimentMode = ExperimentMode.DATA
    run_mode: RunMode = RunMode.ESTIMATION
    b_form: BFormType = BForm.SINGLE.value
    nfe: int = 300
    fd_scheme: str = "BACKWARD"
    time_scaled_end: float = 1.0


@dataclass
class ParameterGuess:
    Lp: float
    sigma: float
    B: Optional[Union[float, Dict[int, float]]] = None
    beta_0: Optional[float] = None
    beta_1: Optional[float] = None
    beta_2: Optional[float] = None
    beta_3: Optional[float] = None
    S0: Optional[float] = None
    S: Optional[float] = None


R_BAR_CM3_PER_UMOL_K = 8.314e-5
NU_CM2_S = 8.927e-3
CELL_DIAMETER_CM = 2.2860
RPM_DEFAULT = 350.0
MH_ML = 0.25
DIFFUSIVITY_CM2_S = {
    "K": 1.960e-5,
    "Na": 1.334e-5,
    "Li": 1.03e-5,
    "Mg": 0.706e-5,
    "Ca": 0.792e-5,
    "La": 0.62e-5,
}


def _first_non_nan(x: Optional[np.ndarray], default: float = 1e-6) -> float:
    if x is None:
        return float(default)
    arr = np.asarray(x, dtype=float).reshape(-1)
    for val in arr:
        if not np.isnan(val):
            return float(val)
    return float(default)


def _salt_key(component_names: Optional[List[str]]) -> Optional[str]:
    if not component_names:
        return None
    s = str(component_names[0]).strip()
    letters = "".join(c for c in s if c.isalpha())
    for k in DIFFUSIVITY_CM2_S:
        if k in letters:
            return k
    return None


def compute_mass_transfer_coeff(exp: ExperimentalData) -> float:
    b = CELL_DIAMETER_CM
    nu = NU_CM2_S
    v = RPM_DEFAULT / 60.0 * np.pi * b

    d_key = _salt_key(exp.component_names)
    if d_key is None:
        raise ValueError(f"Cannot infer diffusivity from component names: {exp.component_names}")
    D = DIFFUSIVITY_CM2_S[d_key]
    return 0.23 * v ** 0.57 * D ** 0.67 / (nu ** 0.24 * b ** 0.43)


def _scaled_times(exp: ExperimentalData) -> Tuple[Dict[int, float], Dict[int, float], float]:
    t_delay = float(exp.vials[0].time_s[0])
    ti = {}
    tf = {}
    for i, v in enumerate(exp.vials, start=1):
        ti[i] = float(v.time_s[0]) - t_delay
        tf[i] = float(v.time_s[-1]) - t_delay
    return ti, tf, t_delay


def _default_guess(exp: ExperimentalData) -> ParameterGuess:
    return ParameterGuess(
        Lp=float(exp.Lp0 if exp.Lp0 is not None else 5.0),
        sigma=float(exp.sigma0 if exp.sigma0 is not None else 0.9),
        B=float(exp.B0 if exp.B0 is not None else 0.5),
        beta_0=float(exp.B0 if exp.B0 is not None else 0.5),
        beta_1=0.5,
        beta_2=0.0,
        beta_3=0.0,
        S0=0.0,
        S=0.0,
    )


def theta_component_list(m: pyo.ConcreteModel, options: ModelOptions) -> List[pyo.ComponentData]:
    b_form = str(options.b_form).lower()
    comps: List[pyo.ComponentData] = [m.Lp]

    if b_form == BForm.SINGLE.value:
        comps.append(m.B)
    elif b_form == BForm.PERVIAL.value:
        comps.extend(m.B[n] for n in m.n_vial)
    elif b_form == BForm.CONVECTION.value:
        comps.append(m.beta_0)
        comps.append(m.beta_1)
    else:
        raise ValueError(f"Unsupported b_form '{options.b_form}' for theta list.")

    comps.append(m.sigma)
    return comps


def theta_names(m: pyo.ConcreteModel, options: ModelOptions) -> List[str]:
    return [c.name for c in theta_component_list(m, options)]


def model_construct_inter_v23(
    exp: ExperimentalData,
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
) -> pyo.ConcreteModel:
    """Build the dynamic diafiltration model (v23)."""
    if guess is None:
        guess = _default_guess(exp)

    if exp.delP_bar is None or exp.Temp_K is None or exp.Am_cm2 is None or exp.rho_g_cm3 is None:
        raise ValueError("Missing required experiment inputs (delP_bar, Temp_K, Am_cm2, rho_g_cm3).")
    if exp.M_F0_g is None:
        raise ValueError("exp.M_F0_g is required for model construction.")
    if exp.C_D_value is None:
        raise ValueError("exp.C_D_value is required for model construction.")

    delP = float(exp.delP_bar)
    T = float(exp.Temp_K)
    Am = float(exp.Am_cm2)
    rho = float(exp.rho_g_cm3)
    ni = int(exp.num_components if exp.num_components is not None else 1)

    M_F0 = float(exp.M_F0_g)
    M_O = float(exp.M_O_g if exp.M_O_g is not None else 0.0)
    cD = float(exp.C_D_value)
    C_F0 = float(exp.C_F0_value if exp.C_F0_value is not None else _first_non_nan(exp.vials[0].retentate_signal, 1e-6))

    N_VIAL = len(exp.vials)
    TI_dict, TF_dict, _ = _scaled_times(exp)
    k = compute_mass_transfer_coeff(exp)

    m = pyo.ConcreteModel()
    m.n_vial = pyo.RangeSet(1, N_VIAL)
    m.tau = ContinuousSet(bounds=(0.0, float(options.time_scaled_end)))
    m.tf = pyo.Set(initialize=[TF_dict[i] for i in m.n_vial])
    m.ti = pyo.Set(initialize=[TI_dict[i] for i in m.n_vial])

    b_form = str(options.b_form).lower()
    sim_opt = options.run_mode == RunMode.SIMULATION

    m.cD = pyo.Var(domain=pyo.NonNegativeReals, initialize=cD)
    m.cD.fix(cD)
    m.C_H0 = pyo.Param(initialize=1e-6, mutable=True)

    if sim_opt:
        m.Lp = pyo.Param(initialize=float(guess.Lp), mutable=True)
        m.sigma = pyo.Param(initialize=float(guess.sigma), mutable=True)
    else:
        m.Lp = pyo.Var(bounds=(0.5, 50.0), initialize=float(guess.Lp))
        m.sigma = pyo.Var(bounds=(0.0, 1.0), initialize=float(guess.sigma))

    if b_form == BForm.SINGLE.value:
        if sim_opt:
            m.B = pyo.Param(initialize=float(guess.B if guess.B is not None else 0.5), mutable=True)
        else:
            m.B = pyo.Var(bounds=(1e-6, 30.0), initialize=float(guess.B if guess.B is not None else 0.5))
    elif b_form == BForm.PERVIAL.value:
        if isinstance(guess.B, dict):
            default_b = list(guess.B.values())[0] if guess.B else 0.5
            b_init = {i: float(guess.B.get(i, default_b)) for i in range(1, N_VIAL + 1)}
        else:
            b_init = {i: float(guess.B if guess.B is not None else 0.5) for i in range(1, N_VIAL + 1)}
        if sim_opt:
            m.B = pyo.Param(m.n_vial, initialize=b_init, mutable=True)
        else:
            m.B = pyo.Var(m.n_vial, bounds=(1e-6, 30.0), initialize=b_init)
    elif b_form == BForm.CONVECTION.value:
        beta0 = float(guess.beta_0 if guess.beta_0 is not None else 2.0)
        beta1 = float(guess.beta_1 if guess.beta_1 is not None else 0.5)
        if sim_opt:
            m.beta_0 = pyo.Param(initialize=beta0, mutable=True)
            m.beta_1 = pyo.Param(initialize=beta1, mutable=True)
        else:
            m.beta_0 = pyo.Var(bounds=(1 + 1e-6, 50.0), initialize=beta0)
            m.beta_1 = pyo.Var(bounds=(0.0, 1.0), initialize=beta1)
        m.H = pyo.Var(m.n_vial, m.tau, initialize=0.5)
    else:
        raise ValueError(f"Unsupported b_form: {options.b_form}")

    if options.mode in (ExperimentMode.LAG, ExperimentMode.OVERFLOW):
        if sim_opt:
            m.S0 = pyo.Param(initialize=float(guess.S0 if guess.S0 is not None else 0.0), mutable=True)
            m.S = pyo.Param(initialize=float(guess.S if guess.S is not None else 0.0), mutable=True)
        else:
            m.S0 = pyo.Var(initialize=float(guess.S0 if guess.S0 is not None else 0.0))
            m.S = pyo.Var(initialize=float(guess.S if guess.S is not None else 0.0))

    if options.mode != ExperimentMode.DATA:
        m.mF = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=M_F0)
        m.dmF = DerivativeVar(m.mF, wrt=m.tau)

    m.cF = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=C_F0)
    m.cIn = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=C_F0)
    m.cH = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=1e-6)
    m.mV = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=1e-6)
    m.cVmV = pyo.Var(m.n_vial, m.tau, initialize=1e-12)
    m.cV = pyo.Var(m.n_vial, m.tau, domain=pyo.NonNegativeReals, initialize=1e-6)

    m.Jw = pyo.Var(m.n_vial, m.tau)
    m.Js = pyo.Var(m.n_vial, m.tau)
    if b_form == BForm.CONVECTION.value:
        m.Js_exp = pyo.Var(m.n_vial, m.tau, bounds=(1 + 1e-6, 1e4))

    m.dcF = DerivativeVar(m.cF, wrt=m.tau)
    m.dcH = DerivativeVar(m.cH, wrt=m.tau)
    m.dmV = DerivativeVar(m.mV, wrt=m.tau)
    m.dcVmV = DerivativeVar(m.cVmV, wrt=m.tau)

    Tauf = float(options.time_scaled_end)

    def _tf_scale(n: int) -> float:
        return (TF_dict[n] - TI_dict[n]) / Tauf

    if options.mode == ExperimentMode.DATA:
        def ode_cF_rule(mm, n, t):
            return mm.dcF[n, t] == Am * rho / M_F0 * (mm.cD * mm.Jw[n, t] - mm.Js[n, t]) * _tf_scale(n)
        m.ode_cF = pyo.Constraint(m.n_vial, m.tau, rule=ode_cF_rule)
    else:
        def ode_mF_rule(mm, n, t):
            return mm.dmF[n, t] == (-mm.S0 - Am * rho * mm.Jw[n, t]) * _tf_scale(n)
        m.ode_mF = pyo.Constraint(m.n_vial, m.tau, rule=ode_mF_rule)

        def ode_cF_rule(mm, n, t):
            return mm.dcF[n, t] == (1 / mm.mF[n, t]) * ((mm.cF[n, t] - mm.cD) * mm.S0 + Am * rho * (mm.cF[n, t] * mm.Jw[n, t] - mm.Js[n, t])) * _tf_scale(n)
        m.ode_cF = pyo.Constraint(m.n_vial, m.tau, rule=ode_cF_rule)

    def ode_cH_rule(mm, n, t):
        return mm.dcH[n, t] == Am * rho / MH_ML * (mm.Js[n, t] - mm.cH[n, t] * mm.Jw[n, t]) * _tf_scale(n)
    m.ode_cH = pyo.Constraint(m.n_vial, m.tau, rule=ode_cH_rule)

    def ode_mV_rule(mm, n, t):
        return mm.dmV[n, t] == mm.Jw[n, t] * Am * rho * _tf_scale(n)
    m.ode_mV = pyo.Constraint(m.n_vial, m.tau, rule=ode_mV_rule)

    def ode_cVmV_rule(mm, n, t):
        return mm.dcVmV[n, t] == mm.Jw[n, t] * mm.cH[n, t] * Am * rho * _tf_scale(n)
    m.ode_cVmV = pyo.Constraint(m.n_vial, m.tau, rule=ode_cVmV_rule)

    def eqn_cIn_rule(mm, n, t):
        return mm.cIn[n, t] == (mm.cF[n, t] - mm.cH[n, t]) * pyo.exp(mm.Jw[n, t] / k) + mm.cH[n, t]
    m.eqn_cIn = pyo.Constraint(m.n_vial, m.tau, rule=eqn_cIn_rule)

    def eqn_Jw_rule(mm, n, t):
        return mm.Jw[n, t] * 36000 == mm.Lp * (delP - (mm.cIn[n, t] - mm.cH[n, t]) * ni * mm.sigma * R_BAR_CM3_PER_UMOL_K * T)
    m.eqn_Jw = pyo.Constraint(m.n_vial, m.tau, rule=eqn_Jw_rule)

    def eqn_Js_rule(mm, n, t):
        if b_form == BForm.SINGLE.value:
            return mm.Js[n, t] * 10000 == mm.B * (mm.cIn[n, t] - mm.cH[n, t])
        if b_form == BForm.PERVIAL.value:
            return mm.Js[n, t] * 10000 == mm.B[n] * (mm.cIn[n, t] - mm.cH[n, t])
        return mm.Js[n, t] == mm.Jw[n, t] * mm.H[n, t] * (mm.cIn[n, t] * mm.Js_exp[n, t] - mm.cH[n, t]) / (mm.Js_exp[n, t] - 1)
    m.eqn_Js = pyo.Constraint(m.n_vial, m.tau, rule=eqn_Js_rule)

    if b_form == BForm.CONVECTION.value:
        def eqn_Js_exp_rule(mm, n, t):
            return mm.Js_exp[n, t] == pyo.exp(mm.Jw[n, t] / mm.beta_0 * 10000)
        m.eqn_Js_exp = pyo.Constraint(m.n_vial, m.tau, rule=eqn_Js_exp_rule)

        def eqn_H_rule(mm, n, t):
            return mm.H[n, t] == mm.beta_1
        m.eqn_H = pyo.Constraint(m.n_vial, m.tau, rule=eqn_H_rule)

    def eqn_cV_rule(mm, n, t):
        return mm.mV[n, t] * mm.cV[n, t] == mm.cVmV[n, t]
    m.eqn_cV = pyo.Constraint(m.n_vial, m.tau, rule=eqn_cV_rule)

    def cF_linking_rule(mm, n):
        if n == mm.n_vial.last():
            return pyo.Constraint.Skip
        return mm.cF[n, Tauf] == mm.cF[n + 1, 0.0]
    m.cF_linking = pyo.Constraint(m.n_vial, rule=cF_linking_rule)

    def cH_linking_rule(mm, n):
        if n == mm.n_vial.last():
            return pyo.Constraint.Skip
        return mm.cH[n, Tauf] == mm.cH[n + 1, 0.0]
    m.cH_linking = pyo.Constraint(m.n_vial, rule=cH_linking_rule)

    m.con_boundary = pyo.ConstraintList()
    m.con_boundary.add(m.cH[1, 0.0] == m.C_H0)
    m.con_boundary.add(m.cF[1, 0.0] == _first_non_nan(exp.vials[0].retentate_signal, default=1e-6))
    m.con_boundary.add(m.cVmV[1, 0.0] == 1e-12)
    m.con_boundary.add(m.mV[1, 0.0] == 1e-6)
    if options.mode != ExperimentMode.DATA:
        m.con_boundary.add(m.mF[1, 0.0] == M_F0)

    if options.mode != ExperimentMode.DATA and options.run_mode == RunMode.ESTIMATION:
        m.eqn_S = pyo.Constraint(expr=m.mF[m.n_vial.last(), Tauf] - M_F0 == M_O)

    return m


def apply_discretization(m: pyo.ConcreteModel, *, nfe: int = 300, scheme: str = "BACKWARD") -> None:
    pyo.TransformationFactory("dae.finite_difference").apply_to(m, nfe=nfe, scheme=scheme)


def _nearest_tau_value(tau_values: Sequence[float], target: float) -> float:
    arr = np.asarray(tau_values, dtype=float)
    return float(arr[np.argmin(np.abs(arr - target))])


def attach_weighted_least_squares_objective(
    m: pyo.ConcreteModel,
    exp: ExperimentalData,
    *,
    sigma_mass_g: float = 0.01,
    sigma_cV_rel: float = 0.03,
    sigma_cF_rel: float = 0.003,
) -> None:
    if not list(m.tau):
        raise RuntimeError("Model must be discretized before adding measurement objective.")

    tau_vals = sorted(float(t) for t in list(m.tau))
    TI_dict, TF_dict, t_delay = _scaled_times(exp)

    expr = 0.0
    for n in m.n_vial:
        vial = exp.vials[int(n) - 1]
        t_meas = np.asarray(vial.time_s, dtype=float) - t_delay
        dur = TF_dict[int(n)] - TI_dict[int(n)]
        if dur <= 0:
            continue
        t_scaled = (t_meas - TI_dict[int(n)]) / dur

        if vial.mass_g is not None:
            mass = np.asarray(vial.mass_g, dtype=float)
            for tm, y in zip(t_scaled, mass):
                if np.isnan(y):
                    continue
                tau_m = _nearest_tau_value(tau_vals, float(tm))
                expr += ((m.mV[n, tau_m] - float(y)) / sigma_mass_g) ** 2

        cf = np.asarray(vial.retentate_signal, dtype=float) if vial.retentate_signal is not None else np.array([])
        for tm, y in zip(t_scaled, cf):
            if np.isnan(y) or y == 0:
                continue
            tau_m = _nearest_tau_value(tau_vals, float(tm))
            expr += ((m.cF[n, tau_m] - float(y)) / (sigma_cF_rel * abs(float(y)))) ** 2

        if isinstance(vial.cV_avg, (int, float, np.number)) and not np.isnan(float(vial.cV_avg)) and float(vial.cV_avg) != 0:
            tau_m = _nearest_tau_value(tau_vals, 1.0)
            y = float(vial.cV_avg)
            expr += ((m.cV[n, tau_m] - y) / (sigma_cV_rel * abs(y))) ** 2

    m.FirstStageCost = pyo.Expression(expr=0.0)
    m.SecondStageCost = pyo.Expression(expr=expr)
    m.Total_Cost_Objective = pyo.Objective(expr=m.FirstStageCost + m.SecondStageCost, sense=pyo.minimize)


def label_parmest_and_doe_suffixes(
    m: pyo.ConcreteModel,
    exp: ExperimentalData,
    options: ModelOptions,
    *,
    design_vars: Optional[List[pyo.ComponentData]] = None,
    sigma_mass_g: float = 0.01,
    sigma_cV_rel: float = 0.03,
    sigma_cF_rel: float = 0.003,
) -> None:
    """Attach experiment_outputs / measurement_error / unknown_parameters / experiment_inputs suffixes."""
    if not list(m.tau):
        raise RuntimeError("Model must be discretized before labeling outputs.")

    tau_vals = sorted(float(t) for t in list(m.tau))
    TI_dict, TF_dict, t_delay = _scaled_times(exp)

    m.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    m.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)

    for n in m.n_vial:
        vial = exp.vials[int(n) - 1]
        t_meas = np.asarray(vial.time_s, dtype=float) - t_delay
        dur = TF_dict[int(n)] - TI_dict[int(n)]
        if dur <= 0:
            continue
        t_scaled = (t_meas - TI_dict[int(n)]) / dur

        if vial.mass_g is not None:
            mass = np.asarray(vial.mass_g, dtype=float)
            for tm, y in zip(t_scaled, mass):
                if np.isnan(y):
                    continue
                tau_m = _nearest_tau_value(tau_vals, float(tm))
                key = m.mV[n, tau_m]
                m.experiment_outputs[key] = float(y)
                m.measurement_error[key] = float(sigma_mass_g)

        cf = np.asarray(vial.retentate_signal, dtype=float) if vial.retentate_signal is not None else np.array([])
        for tm, y in zip(t_scaled, cf):
            if np.isnan(y):
                continue
            tau_m = _nearest_tau_value(tau_vals, float(tm))
            key = m.cF[n, tau_m]
            m.experiment_outputs[key] = float(y)
            m.measurement_error[key] = float(max(sigma_cF_rel * abs(float(y)), 1e-8))

        if isinstance(vial.cV_avg, (int, float, np.number)) and not np.isnan(float(vial.cV_avg)):
            tau_m = _nearest_tau_value(tau_vals, 1.0)
            key = m.cV[n, tau_m]
            y = float(vial.cV_avg)
            m.experiment_outputs[key] = y
            m.measurement_error[key] = float(max(sigma_cV_rel * abs(y), 1e-8))

    m.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    for comp in theta_component_list(m, options):
        m.unknown_parameters[comp] = pyo.value(comp)

    m.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    if design_vars is None:
        design_vars = [m.cD]
    for dv in design_vars:
        m.experiment_inputs[dv] = None


def model_construct_for_parmest_v23(
    exp: ExperimentalData,
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
) -> pyo.ConcreteModel:
    """Convenience builder used by parmest/doe Experiment wrapper."""
    m = model_construct_inter_v23(exp, options=options, guess=guess)
    apply_discretization(m, nfe=options.nfe, scheme=options.fd_scheme)
    attach_weighted_least_squares_objective(m, exp)
    label_parmest_and_doe_suffixes(m, exp, options)
    return m


class DiafiltrationExperimentV23(ParmestExperiment):
    """ParmEst/DoE Experiment wrapper for one loaded ExperimentalData object."""

    def __init__(
        self,
        exp: ExperimentalData,
        options: ModelOptions,
        guess: Optional[ParameterGuess] = None,
    ):
        super().__init__(model=None)
        self.exp = exp
        self.options = options
        self.guess = guess

    def get_labeled_model(self):
        self.model = model_construct_for_parmest_v23(self.exp, self.options, self.guess)
        return self.model


def build_experiment_list_v23(
    experiments: List[ExperimentalData],
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
) -> List[DiafiltrationExperimentV23]:
    return [DiafiltrationExperimentV23(exp=e, options=options, guess=guess) for e in experiments]


def estimate_parameters_with_parmest_v23(
    experiments: List[ExperimentalData],
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
    *,
    calc_cov: bool = True,
    cov_n: Optional[int] = None,
    solver: str = "ef_ipopt",
    solver_options: Optional[Dict[str, object]] = None,
    tee: bool = False,
) -> Dict[str, object]:
    """Run parameter estimation via pyomo.contrib.parmest.Estimator."""
    exp_list = build_experiment_list_v23(experiments, options, guess)
    estimator = Estimator(exp_list, tee=tee, solver_options=solver_options)

    try:
        out = estimator.theta_est(solver=solver, calc_cov=calc_cov, cov_n=cov_n)
        result: Dict[str, object] = {"raw": out, "estimator": estimator}
    except RuntimeError as err:
        # Parmest in this release only accepts ef_ipopt. Provide an ipopt
        # fallback to keep workflows moving in environments without ef_ipopt.
        if solver == "ipopt" and "Unknown solver in Q_Opt=ipopt" in str(err):
            result = _estimate_parameters_ipopt_fallback(
                experiments=experiments,
                options=options,
                guess=guess,
                solver_options=solver_options,
                tee=tee,
            )
            result["warning"] = (
                "ParmEst ef_ipopt is unavailable; used ipopt fallback on a single labeled experiment model. "
                "Covariance is not computed in fallback mode."
            )
            return result
        raise

    if isinstance(out, tuple):
        if len(out) >= 1:
            result["objective"] = out[0]
        if len(out) >= 2:
            result["theta"] = out[1]
        if len(out) >= 3:
            result["returned_values"] = out[2]
        if len(out) >= 4:
            result["covariance"] = out[3]

    return result


def _estimate_parameters_ipopt_fallback(
    experiments: List[ExperimentalData],
    options: ModelOptions,
    guess: Optional[ParameterGuess],
    solver_options: Optional[Dict[str, object]],
    tee: bool,
) -> Dict[str, object]:
    if len(experiments) != 1:
        raise RuntimeError(
            "ipopt fallback currently supports one experiment at a time. "
            "Install ef_ipopt for multi-experiment parmest solve."
        )

    exp = experiments[0]
    m = model_construct_for_parmest_v23(exp, options=options, guess=guess)
    solver = pyo.SolverFactory("ipopt")
    if solver_options:
        for k, v in solver_options.items():
            solver.options[k] = v
    raw = solver.solve(m, tee=tee)

    theta_vals = {c.name: float(pyo.value(c)) for c in theta_component_list(m, options)}
    objective = float(pyo.value(m.Total_Cost_Objective))

    return {
        "raw": raw,
        "objective": objective,
        "theta": pd.Series(theta_vals),
        "model": m,
    }


def covariance_to_fim(covariance: Union[np.ndarray, object]) -> np.ndarray:
    cov = np.asarray(covariance, dtype=float)
    if cov.ndim != 2 or cov.shape[0] != cov.shape[1]:
        raise ValueError("Covariance must be a square 2D matrix.")
    return np.linalg.pinv(cov)


def d_optimality(fim: Union[np.ndarray, object], *, log_scale: bool = True) -> float:
    fim_arr = np.asarray(fim, dtype=float)
    sign, logdet = np.linalg.slogdet(fim_arr)
    if sign <= 0:
        return float("-inf") if log_scale else 0.0
    return float(logdet) if log_scale else float(np.exp(logdet))


def run_doe_with_pyomo_v23(
    experiment: DiafiltrationExperimentV23,
    *,
    fd_formula: str = "central",
    step: float = 0.001,
    objective_option: str = "determinant",
    solver_name: str = "ipopt",
    tee: bool = False,
) -> Dict[str, object]:
    """Run DoE FIM analysis/optimization using pyomo.contrib.doe."""
    solver = pyo.SolverFactory(solver_name)
    doe_obj = DesignOfExperiments(
        experiment=experiment,
        fd_formula=fd_formula,
        step=step,
        objective_option=objective_option,
        solver=solver,
        tee=tee,
    )
    # First attempt a direct FIM computation at nominal settings.
    # If this fails, fall back to full DoE run.
    results = None
    fim = None
    try:
        fim = doe_obj.compute_FIM(method="sequential")
        results = {"method": "compute_FIM_sequential"}
    except Exception:
        results = doe_obj.run_doe()
        fim = doe_obj.get_FIM()
    return {
        "doe": doe_obj,
        "results": results,
        "fim": np.asarray(fim, dtype=float),
        "d_opt_logdet": d_optimality(fim, log_scale=True),
    }


def parmest_callback_factory_v23(
    experiments: List[ExperimentalData],
    options: ModelOptions,
    guess: Optional[ParameterGuess] = None,
) -> Callable[[int], pyo.ConcreteModel]:
    """Legacy-compatible callback factory retained for quick usage."""

    def _callback(idx: int) -> pyo.ConcreteModel:
        exp = experiments[idx]
        return model_construct_for_parmest_v23(exp, options=options, guess=guess)

    return _callback
