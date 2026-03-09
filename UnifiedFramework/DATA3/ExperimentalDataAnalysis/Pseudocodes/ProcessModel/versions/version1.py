#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 11:15:52 2026

@author: Keshav Kasturi Rangan, Alex Dowling
"""

"""
Diafiltration process model construction (Pyomo) — OOP-style pseudocode.

Goal:
    Replace the legacy dict-based `model_construct_inter(data_stru, ...)` with a
    clear, class-field-based workflow that consumes `ExperimentalData` + options
    and returns a Pyomo ConcreteModel.

Design rules (matching your loader work):
    1) "Light" classes: they only store fields; no heavy logic inside them.
    2) Supporting functions do the work (compute constants, attach variables,
       attach constraints, discretize, etc.).
    3) Keep signals/measurements separate from the mechanistic model:
       a later measurement model maps model states -> sensor signals.

Notes:
    - This is pseudocode: it uses Pyomo names, but many details are stubs.
    - The structure mirrors your legacy function but makes responsibilities explicit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Dict, Any, Tuple

import numpy as np

# Pyomo pseudocode imports (real code will import from pyomo.environ / pyomo.dae)
# from pyomo.environ import ConcreteModel, Var, Param, Constraint, Set, RangeSet, NonNegativeReals, ConstraintList
# from pyomo.dae import ContinuousSet, DerivativeVar
# from pyomo.environ import exp  # Pyomo exp()

# -----------------------------------------------------------------------------#
# 1) Light option / spec classes (fields only)
# -----------------------------------------------------------------------------#

@dataclass
class ModelOptions:
    """
    Stores modeling choices that change the Pyomo structure.

    Inputs:
        mode:
            Experiment mode: {"DATA", "Lag", "Overflow"} (match legacy semantics).
        sim_opt:
            True  -> simulation with fixed Params (no estimation).
            False -> parameter estimation (Vars for parameters).
        B_form:
            Solute permeability parameterization:
                "single"      -> constant B
                "pervial"     -> discrete B per vial
                "convection"  -> convection-diffusion style
                number (0, -0.5, 0.5, 1, 2, 3, ...) -> polynomial/exponent dependence forms
        time_scaled_end:
            End of scaled time domain (tau in [0, time_scaled_end]). Legacy uses 1.
    """
    mode: str = "DATA"
    sim_opt: bool = False
    B_form: Any = "single"            # keep flexible; validated elsewhere
    time_scaled_end: float = 1.0


@dataclass
class ParameterGuess:
    """
    Stores initial guesses / fixed values for parameters used by the model.

    For sim_opt=True:
        These become mutable Params (or Params) in the Pyomo model.

    For sim_opt=False:
        These become initial values for Vars.

    NOTE:
        theta0 should be numeric (vector/array). Keep it numeric here.
    """
    Lp: float
    B: Any
    sigma: float
    theta0: np.ndarray

    # Optional coefficients for B-forms (legacy beta_0..beta_3)
    beta_0: Optional[float] = None
    beta_1: Optional[float] = None
    beta_2: Optional[float] = None
    beta_3: Optional[float] = None

    # Optional lag/overflow parameters
    S0: Optional[float] = None
    S: Optional[float] = None


@dataclass
class KnownConstants:
    """
    Stores known physical constants and experiment-level constants.

    Keep them centralized so the "model builder" functions do not
    re-derive or hardcode scattered values.
    """
    # Known physical constants
    R_bar_cm3_umol_K: float = 8.314e-5  # [cm^3*bar / umol / K]
    nu_cm2_s: float = 8.927e-3          # [cm^2/s] kinematic viscosity (legacy)
    b_cm: float = 2.2860                # [cm] stirred cell diameter (legacy)
    rpm: float = 350.0                  # [rev/min] (legacy)

    # A small holdup volume used in legacy (mH = 0.25 mL)
    mH_mL: float = 0.25


# -----------------------------------------------------------------------------#
# 2) Supporting helper functions (work happens here, not inside classes)
# -----------------------------------------------------------------------------#

def compute_avg_velocity(constants: KnownConstants) -> float:
    """
    Computes an average velocity-like quantity used in the legacy mass transfer correlation.

    Returns:
        v_cm_s: float
            A velocity scale [cm/s] computed from rpm and geometry.

    Note:
        This is legacy logic. If you later replace with a better correlation, this is
        the single function to swap.
    """
    # Legacy: v = (rpm/60) * pi * b
    v_cm_s = (constants.rpm / 60.0) * np.pi * constants.b_cm
    return v_cm_s


def select_diffusivity_from_components(component_names: Any) -> float:
    """
    Selects a diffusion coefficient D based on solute identity.

    Inputs:
        component_names:
            Could be a string or list of strings. Legacy checks 'K' substring.

    Returns:
        D_cm2_s: float

    Raises:
        NotImplementedError if unknown solute identity.
    """
    # Legacy behavior: if name contains 'K' -> use D for K+.
    if isinstance(component_names, str) and ("K" in component_names):
        return 1.960e-5  # [cm^2/s] - legacy value for K+
    raise NotImplementedError("Diffusivity selection not implemented for these components.")


def compute_mass_transfer_k(v_cm_s: float, D_cm2_s: float, constants: KnownConstants) -> float:
    """
    Computes mass transfer coefficient k using the legacy empirical correlation.

    Returns:
        k: float
            Mass transfer coefficient (units consistent with the legacy model).
    """
    # Legacy: k = 0.23 * v^0.57 * D^0.67 / (nu^0.24 * b^0.43)
    k = 0.23 * (v_cm_s ** 0.57) * (D_cm2_s ** 0.67) / ((constants.nu_cm2_s ** 0.24) * (constants.b_cm ** 0.43))
    return k


def compute_vial_time_bounds(exp: "ExperimentalData") -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Computes per-vial unscaled initial/final times and a common delay shift.

    Inputs:
        exp:
            ExperimentalData containing ordered vials with time_s arrays.

    Returns:
        ti_s: np.ndarray
            Initial times per vial (shifted by delay), shape (N_vial,)
        tf_s: np.ndarray
            Final times per vial (shifted by delay), shape (N_vial,)
        t_delay_s: float
            The time delay shift used (legacy uses first timestamp of vial 1).
    """
    # Legacy: t_delay = first time in first vial
    t_delay_s = float(exp.vials[0].time_s[0])

    ti_s = np.array([float(v.time_s[0]) - t_delay_s for v in exp.vials], dtype=float)
    tf_s = np.array([float(v.time_s[-1]) - t_delay_s for v in exp.vials], dtype=float)

    return ti_s, tf_s, t_delay_s


# -----------------------------------------------------------------------------#
# 3) Model construction pipeline (each step mutates/attaches to the Pyomo model)
# -----------------------------------------------------------------------------#

def build_diafiltration_model(
    exp: "ExperimentalData",
    options: ModelOptions,
    guess: ParameterGuess,
    constants: Optional[KnownConstants] = None,
) -> "ConcreteModel":
    """
    Main entrypoint: constructs and returns a Pyomo model for diafiltration.

    Inputs:
        exp:
            ExperimentalData (class-field-based) loaded from XLSX/MAT.
        options:
            ModelOptions controlling structure (mode, sim_opt, B_form).
        guess:
            ParameterGuess providing initialization / fixed values.
        constants:
            KnownConstants bundle (defaults if None).

    Output:
        m:
            Pyomo ConcreteModel with sets, variables, constraints, and ICs attached.
    """
    if constants is None:
        constants = KnownConstants()

    # -------------------------------------------------------------------------
    # A) Extract experiment-level values (these should come from exp fields)
    # -------------------------------------------------------------------------
    delP_bar = exp.delP_bar   # applied pressure [bar]
    T_K = exp.Temp_K          # temperature [K]
    Am_cm2 = exp.Am_cm2       # membrane area [cm^2]
    rho_g_cm3 = exp.rho_g_cm3 # density [g/cm^3]
    ni = exp.num_components   # legacy "ni" = number of dissolved species (be careful: naming mismatch possible)

    # Dialysate concentration (legacy uses cD as constant input)
    cD_value = exp.C_D_value

    # Initial feed mass and concentration (legacy uses M_F0 and C_F0)
    M_F0_g = exp.M_F0_g
    C_F0 = exp.C_F0_value

    # Overflow mass (legacy uses M_O for constraints in Lag/Overflow modes)
    M_O_g = exp.M_O_g

    # -------------------------------------------------------------------------
    # B) Compute physical correlation terms used inside constraints
    # -------------------------------------------------------------------------
    v_cm_s = compute_avg_velocity(constants)
    D_cm2_s = select_diffusivity_from_components(exp.component_names if exp.component_names is not None else "")
    k_mass = compute_mass_transfer_k(v_cm_s, D_cm2_s, constants)

    # -------------------------------------------------------------------------
    # C) Time bookkeeping (per-vial ti/tf) and scaled time domain tau ∈ [0, 1]
    # -------------------------------------------------------------------------
    ti_s, tf_s, t_delay_s = compute_vial_time_bounds(exp)

    # Scaled time is always [0, Tauf] where legacy uses Tauf=1.
    Tauf = float(options.time_scaled_end)

    # Useful mapping per vial:
    #   (tf - ti)/Tauf is the scale factor used in ODE constraints (legacy pattern).
    # scale_factor[n] = (tf_s[n] - ti_s[n]) / Tauf

    # -------------------------------------------------------------------------
    # D) Create Pyomo model and define sets
    # -------------------------------------------------------------------------
    # m = ConcreteModel()
    m = ConcreteModel()  # pseudocode: in real code import ConcreteModel

    # Vial index set
    N_VIAL = len(exp.vials)
    m.n_vial = RangeSet(N_VIAL)

    # Scaled continuous time set
    m.tau = ContinuousSet(bounds=(0.0, Tauf))

    # Store ti/tf as Params so they are available inside constraints
    # (In legacy, these were also stored in TI_dict / TF_dict.)
    # m.ti = Param(m.n_vial, initialize={n: ti_s[n-1] for n in m.n_vial}, mutable=False)
    # m.tf = Param(m.n_vial, initialize={n: tf_s[n-1] for n in m.n_vial}, mutable=False)
    m.ti = Param(m.n_vial, initialize={n: float(ti_s[n - 1]) for n in m.n_vial})
    m.tf = Param(m.n_vial, initialize={n: float(tf_s[n - 1]) for n in m.n_vial})

    # Convenience expression for the scale factor used in each vial’s ODE
    # m.time_scale[n] = (m.tf[n] - m.ti[n]) / Tauf
    m.time_scale = Param(
        m.n_vial,
        initialize={n: float(tf_s[n - 1] - ti_s[n - 1]) / Tauf for n in m.n_vial},
    )

    # -------------------------------------------------------------------------
    # E) Define known inputs / fixed constants as Params on the model
    # -------------------------------------------------------------------------
    m.delP = Param(initialize=float(delP_bar))
    m.T = Param(initialize=float(T_K))
    m.Am = Param(initialize=float(Am_cm2))
    m.rho = Param(initialize=float(rho_g_cm3))
    m.R = Param(initialize=float(constants.R_bar_cm3_umol_K))
    m.k_mass = Param(initialize=float(k_mass))

    m.cD = Param(initialize=float(cD_value))  # dialysate concentration

    # Small baseline concentration / volume (legacy uses 1e-6)
    m.C_H0 = Param(initialize=1e-6)
    m.C_V0 = Param(initialize=1e-6)

    # -------------------------------------------------------------------------
    # F) Define parameters-to-estimate or fixed parameters depending on sim_opt
    # -------------------------------------------------------------------------
    attach_parameter_block(m, options=options, guess=guess)

    # -------------------------------------------------------------------------
    # G) Define state variables and intermediate variables
    # -------------------------------------------------------------------------
    attach_state_variables(m, options=options, M_F0_g=M_F0_g, C_F0=C_F0)
    attach_intermediate_variables(m, options=options)

    # -------------------------------------------------------------------------
    # H) Define derivatives (w.r.t scaled time tau)
    # -------------------------------------------------------------------------
    attach_derivatives(m, options=options)

    # -------------------------------------------------------------------------
    # I) Define ODE constraints (mass + concentrations) and algebraic relations
    # -------------------------------------------------------------------------
    attach_ode_constraints(
        m,
        options=options,
        M_F0_g=M_F0_g,
        M_O_g=M_O_g,
        mH_mL=constants.mH_mL,
        ni=ni,
    )
    attach_algebraic_constraints(m, options=options)

    # -------------------------------------------------------------------------
    # J) Link end-of-vial to start-of-next-vial constraints
    # -------------------------------------------------------------------------
    attach_vial_linking_constraints(m, options=options)

    # -------------------------------------------------------------------------
    # K) Initial conditions (ConstraintList generator style, like legacy)
    # -------------------------------------------------------------------------
    attach_initial_conditions(m, exp=exp, options=options, M_F0_g=M_F0_g)

    # -------------------------------------------------------------------------
    # L) Discretize the DAE system (required before solving)
    # -------------------------------------------------------------------------
    discretize_model(m)

    return m


# -----------------------------------------------------------------------------#
# 4) “Attach” functions: keep the main builder readable and modular
# -----------------------------------------------------------------------------#

def attach_parameter_block(m: "ConcreteModel", options: ModelOptions, guess: ParameterGuess) -> None:
    """
    Attaches Lp, B, sigma, and any B_form-specific coefficients to the model.

    Behavior:
        - sim_opt=True: use Params (often mutable for scenario tests)
        - sim_opt=False: use Vars for estimation, with bounds + initial values
    """
    if options.sim_opt:
        # Simulation mode: parameters are fixed (Params)
        m.Lp = Param(initialize=float(guess.Lp), mutable=True)

        # B_form determines whether B is scalar Param, indexed Param, or something else
        if options.B_form == "single":
            m.B = Param(initialize=float(guess.B), mutable=True)
        elif options.B_form == "pervial":
            # Initialize per-vial values (guess.B should be dict-like or list-like)
            # m.B = Param(m.n_vial, initialize=..., mutable=True)
            m.B = Param(m.n_vial, initialize=guess.B, mutable=True)
        elif options.B_form == "convection":
            # Legacy uses beta_0 + H[n,tau] with Js_exp
            m.beta_0 = Param(initialize=float(guess.beta_0), mutable=True)
            m.beta_1 = Param(initialize=float(guess.beta_1), mutable=True)
            m.H = Var(m.n_vial, m.tau, initialize=0.5)
        else:
            # Numeric B_form: allow B to vary over (n,tau) with polynomial coefficients
            m.beta_0 = Param(initialize=float(guess.beta_0), mutable=True)
            # Optional higher-order betas depending on exponent form
            if isinstance(options.B_form, (int, float)) and options.B_form != 0:
                m.beta_1 = Param(initialize=float(guess.beta_1), mutable=True)
            if isinstance(options.B_form, (int, float)) and options.B_form > 1:
                m.beta_2 = Param(initialize=float(guess.beta_2), mutable=True)
            if isinstance(options.B_form, (int, float)) and options.B_form > 2:
                m.beta_3 = Param(initialize=float(guess.beta_3), mutable=True)

            # Let B be a Var over time (legacy bounds shown)
            m.B = Var(m.n_vial, m.tau, bounds=(1e-6, 50.0))

        m.sigma = Param(initialize=float(guess.sigma), mutable=True)

    else:
        # Estimation mode: parameters are decision variables (Vars)
        m.Lp = Var(bounds=(0.5, 50.0), initialize=float(guess.Lp))

        if options.B_form == "single":
            m.B = Var(bounds=(1e-6, 30.0), initialize=float(guess.B))
        elif options.B_form == "pervial":
            m.B = Var(m.n_vial, bounds=(1e-6, 30.0), initialize=guess.B)
        elif isinstance(options.B_form, str) and ("convection" in options.B_form):
            m.beta_0 = Var(bounds=(1.0 + 1e-6, 50.0), initialize=float(guess.beta_0))
            m.beta_1 = Var(bounds=(0.0, 1.0), initialize=float(guess.beta_1) if guess.beta_1 is not None else 0.5)
            m.H = Var(m.n_vial, m.tau, initialize=0.5)
        else:
            m.beta_0 = Var(initialize=float(guess.beta_0) if guess.beta_0 is not None else 1.0)
            m.B = Var(m.n_vial, m.tau)
            if isinstance(options.B_form, (int, float)) and options.B_form != 0:
                m.beta_1 = Var(bounds=(-20.0, 20.0), initialize=float(guess.beta_1))
            if isinstance(options.B_form, (int, float)) and options.B_form > 1:
                m.beta_2 = Var(bounds=(-20.0, 20.0), initialize=float(guess.beta_2))
            if isinstance(options.B_form, (int, float)) and options.B_form > 2:
                m.beta_3 = Var(bounds=(-20.0, 20.0), initialize=float(guess.beta_3))

        m.sigma = Var(bounds=(0.0, 1.0), initialize=float(guess.sigma))


def attach_state_variables(m: "ConcreteModel", options: ModelOptions, M_F0_g: float, C_F0: float) -> None:
    """
    Attaches the main state variables (mF, cF, cIn, cH, mV, cVmV, cV) to the model.

    Notes:
        - mode != DATA: mF is a state (Var); in DATA mode legacy treats feed mass as fixed M_F0.
        - All are indexed by (n_vial, tau).
    """
    if options.mode != "DATA":
        m.mF = Var(m.n_vial, m.tau, domain=NonNegativeReals, initialize=float(M_F0_g))

    m.cF = Var(m.n_vial, m.tau, domain=NonNegativeReals, initialize=float(C_F0))
    m.cIn = Var(m.n_vial, m.tau, domain=NonNegativeReals, initialize=float(C_F0))
    m.cH = Var(m.n_vial, m.tau, domain=NonNegativeReals, initialize=1e-6)

    m.mV = Var(m.n_vial, m.tau, domain=NonNegativeReals, initialize=1e-6)
    m.cVmV = Var(m.n_vial, m.tau, initialize=float(m.C_V0) * 1e-6)
    m.cV = Var(m.n_vial, m.tau, domain=NonNegativeReals, initialize=1e-6)


def attach_intermediate_variables(m: "ConcreteModel", options: ModelOptions) -> None:
    """
    Attaches intermediate flux variables Jw, Js, and any special intermediate variables.

    Notes:
        - Some B_forms introduce Js_exp (legacy) and/or other helper vars.
    """
    m.Jw = Var(m.n_vial, m.tau)
    m.Js = Var(m.n_vial, m.tau)

    if isinstance(options.B_form, str) and ("convection" in options.B_form):
        m.Js_exp = Var(m.n_vial, m.tau, bounds=(1.0 + 1e-6, 1e4))
    elif options.B_form == "K":
        m.Js_exp = Var(m.n_vial, m.tau)


def attach_derivatives(m: "ConcreteModel", options: ModelOptions) -> None:
    """
    Attaches derivative variables with respect to scaled time tau.

    Important:
        In pyomo.dae, you use DerivativeVar(state, wrt=m.tau).
        This only becomes a solvable NLP after discretization.
    """
    if options.mode != "DATA":
        m.dmF = DerivativeVar(m.mF, wrt=m.tau)

    m.dcF = DerivativeVar(m.cF, wrt=m.tau)
    m.dcH = DerivativeVar(m.cH, wrt=m.tau)
    m.dmV = DerivativeVar(m.mV, wrt=m.tau)
    m.dcVmV = DerivativeVar(m.cVmV, wrt=m.tau)


def attach_ode_constraints(
    m: "ConcreteModel",
    options: ModelOptions,
    M_F0_g: float,
    M_O_g: float,
    mH_mL: float,
    ni: int,
) -> None:
    """
    Attaches the ODE constraints for the diafiltration system.

    Core idea:
        Each vial uses a scaled time domain tau ∈ [0,1], but represents a real
        time interval [ti, tf]. The mapping introduces a factor:

            d(state)/d(tau) = d(state)/dt * (tf-ti)/Tauf

        The legacy code multiplies RHS by (tf-ti)/Tauf; we do the same via m.time_scale[n].
    """
    # -------------------------------------------------------------------------
    # Example: ODE for cF (legacy structure, simplified pseudocode)
    # -------------------------------------------------------------------------
    def ode_cF_rule(m, n, tau):
        # Shorthand factor for this vial
        s = m.time_scale[n]

        if options.mode == "DATA":
            # Legacy DATA: dcF = Am*rho/M_F0 * (cD*Jw - Js) * s
            return m.dcF[n, tau] == (m.Am * m.rho / float(M_F0_g)) * (m.cD * m.Jw[n, tau] - m.Js[n, tau]) * s

        # For Lag/Overflow, you’ll include the S0/S terms (not expanded here)
        # return m.dcF[n,tau] == ... * s
        return Constraint.Skip

    m.ode_cF = Constraint(m.n_vial, m.tau, rule=ode_cF_rule)

    # Other ODEs (mF, cH, mV, cVmV) are attached similarly:
    #   - ode_mF_rule depends on mode and S0/S
    #   - ode_cH_rule depends on whether this vial is holdup-limited (legacy N_H logic)
    #   - ode_mV_rule and ode_cVmV_rule depend on Jw/Js

    # In the real implementation, you’ll port each legacy branch into its own rule
    # function so each remains readable.


def attach_algebraic_constraints(m: "ConcreteModel", options: ModelOptions) -> None:
    """
    Attaches algebraic relations:
        cIn relation, Jw equation, Js equation, cV relation, plus B_form relations.

    Notes:
        - Keep each equation in a separate rule for readability.
        - Prefer writing constraints in dimensional-consistent form and only
          apply unit conversions in one place (legacy uses factors like 36000, 10000).
    """
    # Example: cIn equation (legacy: cIn = (cF - cH)*exp(Jw/k) + cH)
    def eqn_cIn_rule(m, n, tau):
        return m.cIn[n, tau] == (m.cF[n, tau] - m.cH[n, tau]) * exp(m.Jw[n, tau] / m.k_mass) + m.cH[n, tau]

    m.eqn_cIn = Constraint(m.n_vial, m.tau, rule=eqn_cIn_rule)

    # Example: Jw equation (legacy: Jw*36000 == Lp*(delP - (cIn-cH)*ni*sigma*R*T))
    def eqn_Jw_rule(m, n, tau):
        return m.Jw[n, tau] * 36000.0 == m.Lp * (m.delP - (m.cIn[n, tau] - m.cH[n, tau]) * ni * m.sigma * m.R * m.T)

    m.eqn_Jw = Constraint(m.n_vial, m.tau, rule=eqn_Jw_rule)

    # Example: Js equation depends on B_form (split into helper function in real code)
    def eqn_Js_rule(m, n, tau):
        if options.B_form == "single":
            return m.Js[n, tau] * 10000.0 == m.B * (m.cIn[n, tau] - m.cH[n, tau])
        if options.B_form == "pervial":
            return m.Js[n, tau] * 10000.0 == m.B[n] * (m.cIn[n, tau] - m.cH[n, tau])
        if options.B_form == "convection":
            # legacy: Js = Jw * H * (cIn*Js_exp - cH)/(Js_exp-1)
            return m.Js[n, tau] == m.Jw[n, tau] * m.H[n, tau] * (m.cIn[n, tau] * m.Js_exp[n, tau] - m.cH[n, tau]) / (m.Js_exp[n, tau] - 1.0)

        # numeric exponent forms: Js uses m.B[n,tau] possibly tied to cIn
        return m.Js[n, tau] * 10000.0 == (m.Jw[n, tau] * m.B[n, tau] * 10000.0) * (m.cIn[n, tau] - m.cH[n, tau])

    m.eqn_Js = Constraint(m.n_vial, m.tau, rule=eqn_Js_rule)

    # Example: cV relation (legacy: mV*cV == cVmV)
    def eqn_cV_rule(m, n, tau):
        return m.mV[n, tau] * m.cV[n, tau] == m.cVmV[n, tau]

    m.eqn_cV = Constraint(m.n_vial, m.tau, rule=eqn_cV_rule)

    # B_form tying constraints (beta-polynomials, etc.) should be attached in a dedicated function:
    # attach_B_form_constraints(m, options)


def attach_vial_linking_constraints(m: "ConcreteModel", options: ModelOptions) -> None:
    """
    Links the terminal state of vial n to the initial state of vial n+1 at tau=0.

    Notes:
        - This implements the “multiple shooting” / piecewise time segmentation logic.
        - Legacy links mF, cF, cH, mV, cVmV with many special cases.
        - Keep the special-case logic but isolate it into small readable rules.
    """
    def cF_link_rule(m, n):
        if n == m.n_vial.last():
            return Constraint.Skip
        return m.cF[n, options.time_scaled_end] == m.cF[n + 1, 0.0]

    m.cF_link = Constraint(m.n_vial, rule=cF_link_rule)

    # Similarly: cH, mV, cVmV, and mF (if mode != DATA)
    # Each can be its own rule to keep complexity localized.


def attach_initial_conditions(m: "ConcreteModel", exp: "ExperimentalData", options: ModelOptions, M_F0_g: float) -> None:
    """
    Adds initial conditions using a ConstraintList-like pattern.

    Key idea:
        Initial conditions depend on mode and which signals exist in exp.vials.
        In your new architecture, any “firstNonNan” logic should operate on:
            exp.vials[0].retentate_signal  (or on a measurement model output)
        rather than on dict lookups.

    Note:
        The generator/yield pattern matches the legacy style.
    """
    def first_non_nan(arr: np.ndarray) -> float:
        for x in arr:
            if not np.isnan(x):
                return float(x)
        return float("nan")

    def _init_rule(m):
        # Example: mF initial condition if not DATA mode
        if options.mode != "DATA":
            yield m.mF[1, 0.0] == float(M_F0_g)

        # Example: cF initial condition from first vial’s retentate signal (honest signal)
        # In practice, you may *not* want to set cF from conductivity directly;
        # you would set it from an inferred/assayed concentration (measurement model).
        cF0 = first_non_nan(exp.vials[0].retentate_signal)
        yield m.cF[1, 0.0] == cF0

        # Holdup initial concentration
        yield m.cH[1, 0.0] == float(m.C_H0)

        # Small initial permeate/volume states (legacy uses 1e-6)
        yield m.cVmV[1, 0.0] == 1e-6 * 1e-6
        yield m.mV[1, 0.0] == 1e-6

        # End marker (real Pyomo uses ConstraintList.End)
        yield ConstraintList.End

    m.con_boundary = ConstraintList(rule=_init_rule)


def discretize_model(m: "ConcreteModel") -> None:
    """
    Discretizes the DAE system so the model becomes a finite-dimensional NLP.

    Typical approaches:
        - Finite difference
        - Collocation on finite elements

    In real code, use TransformationFactory('dae.finite_difference') or 'dae.collocation'.
    """
    # Example (finite difference):
    # disc = TransformationFactory("dae.finite_difference")
    # disc.apply_to(m, nfe=50, wrt=m.tau, scheme="BACKWARD")

    # Example (collocation):
    # disc = TransformationFactory("dae.collocation")
    # disc.apply_to(m, nfe=20, ncp=3, wrt=m.tau)

    pass


# -----------------------------------------------------------------------------#
# 5) Placeholder type for your loader’s ExperimentalData (already defined by you)
# -----------------------------------------------------------------------------#

class ExperimentalData:
    """
    Placeholder: your real ExperimentalData dataclass already exists.

    Required for this model pseudocode:
        exp.vials: list of vial objects with time_s and retentate_signal arrays
        exp.delP_bar, exp.Temp_K, exp.Am_cm2, exp.rho_g_cm3, exp.num_components
        exp.M_F0_g, exp.C_F0_value, exp.M_O_g, exp.C_D_value
        exp.component_names
    """
    ...


class ConcreteModel:
    """Placeholder for Pyomo ConcreteModel."""
    ...


class RangeSet:
    """Placeholder for Pyomo RangeSet."""
    def __init__(self, n): ...


class ContinuousSet:
    """Placeholder for Pyomo ContinuousSet."""
    def __init__(self, bounds): ...


class Param:
    """Placeholder for Pyomo Param."""
    def __init__(self, *args, **kwargs): ...


class Var:
    """Placeholder for Pyomo Var."""
    def __init__(self, *args, **kwargs): ...


class Constraint:
    """Placeholder for Pyomo Constraint."""
    Skip = "SKIP"
    def __init__(self, *args, **kwargs): ...


class DerivativeVar:
    """Placeholder for Pyomo DerivativeVar."""
    def __init__(self, *args, **kwargs): ...


class ConstraintList:
    """Placeholder for Pyomo ConstraintList."""
    End = "END"
    def __init__(self, *args, **kwargs): ...


NonNegativeReals = object()  # placeholder


def exp(x):  # placeholder for pyomo.environ.exp
    return np.exp(x)
