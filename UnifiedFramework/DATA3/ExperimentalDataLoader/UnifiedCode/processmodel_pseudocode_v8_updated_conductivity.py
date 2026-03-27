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
    class-based workflow that consumes `ExperimentalData`
    and returns a Pyomo ConcreteModel.

Design rules (matching your loader work):
    1) Classes: they only store fields; no heavy logic inside them.
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
from enum import Enum
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

# -----------------------------------------------------------------------------
# Typing helpers (pseudocode-level)
# -----------------------------------------------------------------------------
# NOTE:
#   In the legacy code, the meaning/shape of "B" depends on `B_form`.
#   We keep this explicit via a small Union rather than leaving it as `Any`.
BInit = Union[
    float,
    np.ndarray,          # e.g., time-series initialization for B[n,t]
    Sequence[float],     # e.g., per-vial list/tuple of B values
    Dict[int, float],    # e.g., {vial_index: B_value}
]

# Pyomo pseudocode imports (real code will import from pyomo.environ / pyomo.dae)
# from pyomo.environ import ConcreteModel, Var, Param, Constraint, Set, RangeSet, NonNegativeReals, ConstraintList
# from pyomo.dae import ContinuousSet, DerivativeVar
# from pyomo.environ import exp  # Pyomo exp()

# -----------------------------------------------------------------------------#
# 1) Classes (fields only)
# -----------------------------------------------------------------------------#

class ExperimentMode(str, Enum):
    """Supported experiment modes (legacy semantics)."""
    DATA = "DATA"
    LAG = "Lag"
    OVERFLOW = "Overflow"


class RunMode(str, Enum):
    """Workflow modes: simulation vs estimation vs (future) OED."""
    SIMULATION = "SIMULATION"
    ESTIMATION = "ESTIMATION"
    OED = "OED"


class BForm(str, Enum):
    """Supported named solute-permeability parameterizations."""
    SINGLE = "single"
    PERVIAL = "pervial"
    CONVECTION = "convection"


@dataclass
class ModelOptions:
    """
    Stores modeling choices that change the Pyomo structure.

    Key design note (per review):
        Prefer an explicit RunMode enum over a boolean like sim_opt, so the
        workflow can grow naturally to include simulation, estimation, and OED.
    """
    mode: ExperimentMode = ExperimentMode.DATA
    run_mode: RunMode = RunMode.ESTIMATION

    # B_form accepts either:
    #   - a named enum case (BForm), OR
    #   - a numeric exponent for legacy polynomial/exponent dependence.
    B_form: Union[BForm, float] = BForm.SINGLE


@dataclass
class ParameterGuess:
    """
    Stores initial guesses / fixed values for parameters used by the model.

    Design choice:
        Build parameter objects as Vars by default.
        In SIMULATION mode, we fix those Vars to supplied values after creation.

    Coupling to B_form:
        The legacy model supports multiple parameterizations for solute transport.
        Which "B-related" fields are used depends on `options.B_form`:

        - BForm.SINGLE:
            Use `B` as a scalar initial guess.

        - BForm.PERVIAL:
            Use `B` as a per-vial initializer (dict keyed by vial index, a list/array,
            or a scalar if you want the same initial guess for every vial).

        - BForm.CONVECTION:
            Use `beta_0` (and optionally `beta_1`), plus an initial guess `H0` for
            the holdup-like state H (if you expose it as a Var).

        - Numeric exponent (legacy polynomial/exponent dependence):
            Use `beta_0...beta_3` as available; `B` may be omitted or used only as a
            convenience initializer for derived terms.

    NOTE:
        theta0 should be numeric (vector/array). Keep it numeric here.
    """
    Lp: float
    B: Optional[BInit] = None
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

    # Optional convection-related initial guess (only meaningful for BForm.CONVECTION)
    H0: Optional[float] = None




# -----------------------------------------------------------------------------
# Optional convenience wrapper (future): bundle options + guesses coherently
# -----------------------------------------------------------------------------

@dataclass
class ModelConfig:
    """Bundles model-building configuration.

    Why this exists (design note for future refactors):
        - ModelOptions influences which parameters exist and how they are used
          (e.g., the chosen B_form changes whether B is scalar / per-vial / time-varying).
        - ParameterGuess must therefore be interpreted *relative to* ModelOptions.

    For now (minimal diffs / clarity):
        - We keep ModelOptions and ParameterGuess as separate inputs to
          build_diafiltration_model(...).
        - If/when the parameterization stabilizes (especially for B_form variants and OED),
          we can pass a single ModelConfig instead.

    Fields:
        options:
            Structural modeling choices (mode, run_mode, B_form, scaling).
        guess:
            Initial values / fixed values used to initialize (and possibly fix) Vars.
    """

    options: ModelOptions
    guess: ParameterGuess


# =============================================================================
# Specified constants (global): physical constants + legacy geometry
# =============================================================================

# Gas constant in legacy units: [cm^3 * bar / (umol * K)]
R_BAR_CM3_BAR_PER_UMOL_K: float = 8.314e-5

# Kinematic viscosity in legacy units: [cm^2 / s]
NU_CM2_PER_S: float = 8.927e-3

# Stirred cell diameter in legacy units: [cm]
STIRRED_CELL_DIAMETER_CM: float = 2.2860

# Legacy stirring speed: [rev / min]
RPM_DEFAULT: float = 350.0

# Legacy holdup volume used in the model: [mL]
MH_ML: float = 0.25

# Optional: diffusivities by component name (extensible; legacy has only K+).
DIFFUSIVITY_CM2_S: Dict[str, float] = {"K+": 1.960e-5}  # [cm^2/s]


# -----------------------------------------------------------------------------#
# 2) Supporting helper functions (work happens here, not inside classes)
# -----------------------------------------------------------------------------#

def compute_avg_velocity_cm_s(*, rpm: float = RPM_DEFAULT, diameter_cm: float = STIRRED_CELL_DIAMETER_CM) -> float:
    """
    Compute the legacy average velocity scale used in the mass-transfer correlation.

    Inputs:
        rpm: stirring speed [rev/min]
        diameter_cm: stirred cell diameter [cm]

    Output:
        v_cm_s: velocity scale [cm/s] used by the legacy correlation.
    """
    # Legacy: v = (rpm/60) * pi * b
    v_cm_s = (rpm / 60.0) * np.pi * diameter_cm
    return v_cm_s


def select_diffusivity_from_components(
    *,
    component_names: Optional[Sequence[str]],
    diffusivity_dict: Dict[str, float] = DIFFUSIVITY_CM2_S,
) -> float:
    """
    Select a diffusion coefficient D [cm^2/s] using component names and a diffusivity mapping.

    Inputs:
        component_names:
            Sequence of component/species names (e.g., ["K+"]).
            If None/empty, the lookup cannot proceed.
        diffusivity_dict:
            Mapping: component name -> diffusion coefficient [cm^2/s].

    Output:
        D_cm2_s:
            Diffusion coefficient [cm^2/s] for the first matching component in component_names.

    Raises:
        ValueError:
            If component_names is None/empty.
        KeyError:
            If no component name matches the keys in diffusivity_dict.
    """
    if not component_names:
        raise ValueError("component_names is missing/empty; cannot select diffusivity D.")

    # Return the first match in the mapping (supports multi-component lists).
    for name in component_names:
        if name in diffusivity_dict:
            return float(diffusivity_dict[name])

    # No matches: raise a helpful error that lists the available keys.
    raise KeyError(
        f"No matching diffusivity for components={list(component_names)}. "
        f"Known diffusivities: {sorted(diffusivity_dict.keys())}"
    )



def compute_mass_transfer_k(*, v_cm_s: float, D_cm2_s: float, nu_cm2_s: float = NU_CM2_PER_S, diameter_cm: float = STIRRED_CELL_DIAMETER_CM) -> float:
    """
    Compute the legacy mass-transfer coefficient k using the empirical correlation.

    Inputs:
        v_cm_s: velocity scale [cm/s]
        D_cm2_s: diffusion coefficient [cm^2/s]
        nu_cm2_s: kinematic viscosity [cm^2/s]
        diameter_cm: stirred cell diameter [cm]

    Output:
        k: mass transfer coefficient (legacy units consistent with the model).
    """
    # Legacy: k = 0.23 * v^0.57 * D^0.67 / (nu^0.24 * b^0.43)
    k = 0.23 * (v_cm_s ** 0.57) * (D_cm2_s ** 0.67) / ((nu_cm2_s ** 0.24) * (diameter_cm ** 0.43))
    return k


def compute_vial_time_bounds(exp: "ExperimentalData") -> Tuple[np.ndarray, np.ndarray, float]:
    """Backward-compatible wrapper around exp.get_vial_switch_times()."""
    return exp.get_vial_switch_times()


# -----------------------------------------------------------------------------#
# 3) Model construction pipeline (each step mutates/attaches to the Pyomo model)
# -----------------------------------------------------------------------------#

def build_diafiltration_model(
    exp: "ExperimentalData",
    options: ModelOptions,
    guess: ParameterGuess,
    *,
    time_scaled_end: float = 1.0,
) -> "ConcreteModel":
    """
    Main entrypoint: constructs and returns a Pyomo model for diafiltration.

    Inputs:
        exp:
            ExperimentalData (class-field-based) loaded from XLSX/MAT.
        options:
            ModelOptions controlling structure (mode, run_mode, B_form).
        guess:
            ParameterGuess providing initialization / fixed values.
        time_scaled_end:
            End of the scaled time domain tau ∈ [0, time_scaled_end].
            This is a scaling/discretization convenience knob (default 1.0) and is
            intentionally not part of ModelOptions.

    Output:
        m:
            Pyomo ConcreteModel with sets, variables, constraints, and ICs attached.
    """


    # -------------------------------------------------------------------------
    # 0) Measurement pre-processing (ONLY if needed)
    # -------------------------------------------------------------------------
    # The process model equations are written in terms of concentrations.
    #
    # For MAT-based legacy data:
    #   - concentration measurements are already stored in the MAT struct (via an
    #     ICP-linear conductivity correlation used in prior publications).
    #
    # For XLSX-based new experiments:
    #   - sensors record conductivity, so you must either:
    #       (A) run the conductivity->concentration conversion step during loading, OR
    #       (B) run it here before model construction.
    #
    # In this pseudocode, we assume you did (A) via load_experiment_easy(..., convert_to_concentration=True).
    # If not, call: apply_conductivity_to_concentration(exp, ...) before continuing.

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


    # Quick sanity checks: the mechanistic model cannot be built without these.
    # The loader/validator may allow None, but model construction should fail fast.
    if delP_bar is None or T_K is None or Am_cm2 is None or rho_g_cm3 is None:
        raise ValueError(
            "Missing required experiment fields for model construction: "
            "delP_bar, Temp_K, Am_cm2, rho_g_cm3 must all be provided."
        )

    # -------------------------------------------------------------------------
    # B) Compute physical correlation terms used inside constraints
    # -------------------------------------------------------------------------
    # Compute the legacy velocity scale used in the mass-transfer correlation.
    # We keep rpm/diameter explicit (with defaults) so the dependency is obvious.
    v_cm_s = compute_avg_velocity_cm_s(
        rpm=RPM_DEFAULT,
        diameter_cm=STIRRED_CELL_DIAMETER_CM,
    )
    D_cm2_s = select_diffusivity_from_components(
        component_names=exp.component_names,
    )
    # Compute the mass-transfer coefficient used in cIn (concentration polarization) relation.
    k_mass = compute_mass_transfer_k(
        v_cm_s=v_cm_s,
        D_cm2_s=D_cm2_s,
        nu_cm2_s=NU_CM2_PER_S,
        diameter_cm=STIRRED_CELL_DIAMETER_CM,
    )

    # -------------------------------------------------------------------------
    # C) Time bookkeeping (per-vial ti/tf) and scaled time domain tau ∈ [0, Tauf]
    # -------------------------------------------------------------------------
    ti_s, tf_s, t_delay_s = exp.get_vial_switch_times()  # per-vial ti/tf and common delay shift

    # Scaled time is always [0, Tauf] where legacy uses Tauf=1.
    Tauf = float(time_scaled_end)

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
    m.R = Param(initialize=float(R_BAR_CM3_BAR_PER_UMOL_K))
    m.k_mass = Param(initialize=float(k_mass))

    m.cD = Param(initialize=float(cD_value))  # dialysate concentration

    # Small baseline concentration / volume (legacy uses 1e-6)
    m.C_H0 = Param(initialize=1e-6)
    m.C_V0 = Param(initialize=1e-6)

    # -------------------------------------------------------------------------
    # F) Declare *all* model parameters as Vars; simulation vs estimation is handled by fixing/unfixing
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
        mH_mL=MH_ML,
        ni=ni,
    )
    attach_algebraic_constraints(m, options=options)

    # -------------------------------------------------------------------------
    # J) Link end-of-vial to start-of-next-vial constraints
    # -------------------------------------------------------------------------
    attach_vial_linking_constraints(m, options=options, tau_end=Tauf)

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

    Design choice (review-driven):
        We *always* create parameters as Vars so the model structure is identical for:
            - simulation (Vars are fixed)
            - estimation (Vars are free)
            - OED (future: Vars are free + design variables are added)

        This avoids having two different model “shapes” (Params vs Vars) depending on run mode.

    Inputs:
        m:
            Pyomo ConcreteModel (mutated in-place).
        options:
            Structural switches (mode/run_mode/B_form).
        guess:
            Numeric initial guesses and/or fixed values.

    Output:
        None (mutates m by attaching Vars and potentially fixing them).
    """

    # -------------------------------------------------------------------------
    # 1) Core transport parameters (always Vars)
    # -------------------------------------------------------------------------
    # Lp: water permeability-like parameter
    m.Lp = Var(bounds=(0.5, 50.0), initialize=float(guess.Lp))
    # NOTE: We always use Vars, even in simulation mode. In simulation we simply fix these Vars.

    # sigma: reflection coefficient (dimensionless)
    m.sigma = Var(bounds=(0.0, 1.0), initialize=float(guess.sigma))
    # In SIMULATION run_mode, we will fix sigma (and any other parameter Vars) to the provided guess.

    # -------------------------------------------------------------------------
    # 2) Solute permeability parameterization depends on B_form
    # -------------------------------------------------------------------------
    if options.B_form == BForm.SINGLE:
        # One constant B for the entire experiment
        m.B = Var(bounds=(1e-6, 30.0), initialize=float(guess.B))

    elif options.B_form == BForm.PERVIAL:
        # One B per vial (discrete, vial-indexed)
        # guess.B is expected to be dict-like {n: value} or list-like length N_VIAL.
        m.B = Var(m.n_vial, bounds=(1e-6, 30.0), initialize=guess.B)

    elif options.B_form == BForm.CONVECTION:
        # Convection-style model: introduces beta_0 (scale) + H[n,tau] + Js_exp.
        # NOTE: Js_exp and H are attached in attach_intermediate_variables().
        if guess.beta_0 is None:
            raise ValueError("B_form=CONVECTION requires guess.beta_0.")
        m.beta_0 = Var(bounds=(1.0 + 1e-6, 50.0), initialize=float(guess.beta_0))

        # beta_1 is used as a constant H in the legacy code; keep it as a Var so it can be estimated.
        m.beta_1 = Var(bounds=(0.0, 1.0), initialize=float(guess.beta_1) if guess.beta_1 is not None else 0.5)

        # H is a time-varying helper state (dimensionless) in this formulation.
        m.H = Var(m.n_vial, m.tau, initialize=0.5)

    else:
        # Numeric exponent forms (legacy): B can vary with time and concentration via beta coefficients.
        # We store beta_0..beta_3 and define constraints later that tie B[n,tau] to cIn or cF.
        if guess.beta_0 is None:
            raise ValueError("Numeric B_form requires guess.beta_0 (baseline).")

        m.beta_0 = Var(initialize=float(guess.beta_0))

        # Higher-order coefficients are optional depending on exponent order.
        if isinstance(options.B_form, (int, float)) and options.B_form != 0:
            m.beta_1 = Var(bounds=(-20.0, 20.0), initialize=float(guess.beta_1) if guess.beta_1 is not None else 0.0)
        if isinstance(options.B_form, (int, float)) and options.B_form > 1:
            m.beta_2 = Var(bounds=(-20.0, 20.0), initialize=float(guess.beta_2) if guess.beta_2 is not None else 0.0)
        if isinstance(options.B_form, (int, float)) and options.B_form > 2:
            m.beta_3 = Var(bounds=(-20.0, 20.0), initialize=float(guess.beta_3) if guess.beta_3 is not None else 0.0)

        # Time-varying B (later tied to concentrations by constraints)
        m.B = Var(m.n_vial, m.tau, bounds=(1e-6, 50.0))

    # -------------------------------------------------------------------------
    # 3) Simulation mode: fix parameter Vars to supplied values
    # -------------------------------------------------------------------------
    if options.run_mode == RunMode.SIMULATION:
        # Simulation: treat 'guess' as fixed parameter values by fixing the parameter Vars.
        # Estimation: leave parameter Vars free so the optimizer can fit them.
        # Fix core parameters
        m.Lp.fix(float(guess.Lp))
        m.sigma.fix(float(guess.sigma))

        # Fix B-related parameters depending on form
        if options.B_form == BForm.SINGLE:
            m.B.fix(float(guess.B))
        elif options.B_form == BForm.PERVIAL:
            # Fix each vial’s B
            for n in m.n_vial:
                m.B[n].fix(float(guess.B[n]) if isinstance(guess.B, dict) else float(guess.B))
        elif options.B_form == BForm.CONVECTION:
            m.beta_0.fix(float(guess.beta_0))
            if guess.beta_1 is not None:
                m.beta_1.fix(float(guess.beta_1))
        else:
            # Fix polynomial coefficients if provided
            m.beta_0.fix(float(guess.beta_0))
            if hasattr(m, "beta_1") and guess.beta_1 is not None:
                m.beta_1.fix(float(guess.beta_1))
            if hasattr(m, "beta_2") and guess.beta_2 is not None:
                m.beta_2.fix(float(guess.beta_2))
            if hasattr(m, "beta_3") and guess.beta_3 is not None:
                m.beta_3.fix(float(guess.beta_3))


def attach_state_variables(m: "ConcreteModel", options: ModelOptions, M_F0_g: float, C_F0: float) -> None:
    """
    Attaches the main state variables (mF, cF, cIn, cH, mV, cVmV, cV) to the model.

    Notes:
        - mode != DATA: mF is a state (Var); in DATA mode legacy treats feed mass as fixed M_F0.
        - All are indexed by (n_vial, tau).
    """
    if options.mode != ExperimentMode.DATA:
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

    if options.B_form == BForm.CONVECTION:
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
    if options.mode != ExperimentMode.DATA:
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

        if options.mode == ExperimentMode.DATA:
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
        if options.B_form == BForm.SINGLE:
            return m.Js[n, tau] * 10000.0 == m.B * (m.cIn[n, tau] - m.cH[n, tau])
        if options.B_form == BForm.PERVIAL:
            return m.Js[n, tau] * 10000.0 == m.B[n] * (m.cIn[n, tau] - m.cH[n, tau])
        if options.B_form == BForm.CONVECTION:
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


def attach_vial_linking_constraints(m: "ConcreteModel", options: ModelOptions, *, tau_end: float) -> None:
    """
    Links the terminal state of vial n to the initial state of vial n+1 at tau=0.

    Inputs:
        m:
            Pyomo ConcreteModel.
        options:
            Structural switches (mainly mode).
        tau_end:
            The end-point of the scaled time domain (typically 1.0). We pass this
            explicitly so we do not rely on ModelOptions for scaling parameters.

    Notes:
        - This implements the “multiple shooting” / piecewise time segmentation logic.
        - The legacy model links several states (mF, cF, cH, mV, cVmV) with special cases.
        - Start with cF linking as a template; add others similarly.
    """

    # Link cF at the end of vial n to the start of vial n+1.
    def cF_link_rule(m, n):
        if n == m.n_vial.last():
            return Constraint.Skip
        return m.cF[n, tau_end] == m.cF[n + 1, 0.0]

    m.cF_link = Constraint(m.n_vial, rule=cF_link_rule)

    # Similarly: cH, mV, cVmV, and mF (if mode != DATA)
    # Each should be its own rule to keep branching logic localized.


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
        if options.mode != ExperimentMode.DATA:
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
    

    Future architecture note:
        This pseudocode keeps ExperimentalData *interface-compatible* with common
        experiment containers used by Pyomo parameter estimation / design-of-experiments
        tooling (e.g., ParmEst / Pyomo.DoE). We are not implementing inheritance here,
        but we leave stubs for methods like get_inputs() / get_measurements() so a
        later subclass can slot in cleanly.

    """

    vials: Any  # expected: ordered list where each vial has time_s array

    # ---------------------------------------------------------------------
    # Optional interface hooks (future): align with ParmEst / Pyomo.DoE patterns
    # ---------------------------------------------------------------------
    def get_inputs(self) -> Any:
        """Return manipulated inputs / operating conditions for this experiment.

        Stub in pseudocode:
            In a real implementation, this could return time-series inputs or per-vial
            operating conditions (pressure schedules, dialysate settings, etc.).
        """
        raise NotImplementedError

    def get_measurements(self) -> Any:
        """Return measured outputs associated with this experiment.

        Stub in pseudocode:
            In a real implementation, this could return mappings from measurement names
            to arrays (or DataFrames) aligned with the model's time grid.
        """
        raise NotImplementedError

    def get_vial_switch_times(self) -> Tuple[np.ndarray, np.ndarray, float]:
        """Compute per-vial (ti, tf) arrays and a shared delay shift.

        Why this is a method:
            This computation uses *only* fields already stored on ExperimentalData
            (namely the ordered list of vials and their time_s arrays), so keeping
            it here makes call sites read naturally:
                ti_s, tf_s, t_delay_s = exp.get_vial_switch_times()

        Legacy-compatible definition:
            - t_delay_s is defined as the first timestamp in vial 1.
            - ti_s[i] = vials[i].time_s[0]  - t_delay_s
            - tf_s[i] = vials[i].time_s[-1] - t_delay_s

        Returns:
            ti_s: shifted vial start times [s], shape (N_vial,)
            tf_s: shifted vial end times [s], shape (N_vial,)
            t_delay_s: delay shift applied to all vials [s]
        """
        # Defensive: this method is only meaningful if vials exist and each vial has time.
        if not self.vials:
            raise ValueError("ExperimentalData.vials is empty; cannot compute vial switch times.")

        if self.vials[0].time_s is None or len(self.vials[0].time_s) == 0:
            raise ValueError("Vial 1 has missing/empty time_s; cannot compute t_delay_s.")

        # Legacy delay shift: align all vials relative to the first timestamp in vial 1.
        t_delay_s = float(self.vials[0].time_s[0])

        # Per-vial initial/final times shifted by t_delay_s.
        ti_s = np.array([float(v.time_s[0]) - t_delay_s for v in self.vials], dtype=float)
        tf_s = np.array([float(v.time_s[-1]) - t_delay_s for v in self.vials], dtype=float)

        return ti_s, tf_s, t_delay_s

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