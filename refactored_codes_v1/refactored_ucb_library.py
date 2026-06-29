"""
Unified Codebase Library for memebrane separation - Modeling, ParmEst, UQ, DoE
Keshav Kasturi Rangan
University of Notre Dame
"""

import numpy as np
import pandas as pd
import scipy.io as spio
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.patheffects as patheffects
from matplotlib.lines import Line2D
from pathlib import Path
import idaes
import time
import copy
import os
import json
import math
import ctypes.util
import re
from sklearn.metrics import r2_score
from dataclasses import dataclass, field, replace

try:
    import seaborn as sns
except ImportError:  # pragma: no cover - optional plotting dependency
    sns = None

from pyomo.environ import *
from pyomo.dae import *
import idaes.core.util.scaling as iscale
from pyomo.core.base.componentuid import ComponentUID
from pyomo.contrib.parmest.experiment import Experiment as ParmestExperiment
from pyomo.contrib.parmest.parmest import Estimator, SSE, SSE_weighted, compute_covariance_matrix
from pyomo.contrib.doe.doe import DesignOfExperiments


FIGURES_DIR = os.path.join("UnifiedFramework", "DATA3", "figures")


def _has_casadi():
    """Return True when CasADi is importable in the current environment."""
    try:
        import casadi  # noqa: F401
    except ImportError:
        return False
    return True


def _preferred_ipopt_linear_solver():
    """Prefer HSL when present, otherwise fall back to the portable MUMPS build."""
    # Keep this overridable from the environment so a machine-specific solver
    # can be selected without touching the code.
    override = os.environ.get("DIAFILTRATION_IPOPT_LINEAR_SOLVER", "").strip().lower()
    if override in {"ma97", "ma57", "mumps"}:
        return override
    # The packaged Ipopt build in this environment advertises HSL support, but
    # ctypes lookup is not always reliable. Prefer the sparse HSL solvers when
    # a user has not explicitly overridden the choice; fall back to MUMPS if
    # the solver binary does not accept the requested linear solver at runtime.
    if ctypes.util.find_library("hsl"):
        return "ma97"
    return "mumps"


def _build_ipopt_solver():
    """Create an Ipopt solver configured for the current machine."""
    solver = SolverFactory("ipopt")
    linear_solver = _preferred_ipopt_linear_solver()
    solver.options["linear_solver"] = linear_solver
    return solver


def _resolve_ipopt_cpu_time_limit(*, solver_max_cpu_time=None, multistart=False):
    """Resolve an optional Ipopt CPU-time cap.

    Multistart fits are exploratory by nature, so they get a modest default
    cap unless the caller or environment overrides it.
    """
    if solver_max_cpu_time is not None:
        try:
            value = float(solver_max_cpu_time)
            return value if value > 0 else None
        except Exception:
            return None

    env_value = os.environ.get("DIAFILTRATION_IPOPT_MAX_CPU_TIME", "").strip()
    if env_value:
        try:
            value = float(env_value)
            return value if value > 0 else None
        except Exception:
            pass

    if multistart:
        return 600.0
    return None


def _canonical_mode(mode):
    """Normalize legacy mode spellings used by the interactive CLI."""
    text = str(mode or "DATA").strip()
    lower = text.lower()
    if lower == "data":
        return "DATA"
    if lower == "overflow":
        return "Overflow"
    if lower == "lag":
        return "Lag"
    return text


@dataclass(frozen=True)
class ExperimentalRun:
    """
    One experimental run family with multiple file-backed variants.

    The same physical run can appear in different notebook/paper forms, such as
    the base fit and the concentration-polarization variant. Each variant is a
    separate instance of the same run family.
    """

    run_id: str
    root: Path
    data_file: str
    variant: str = "base"
    subdir: str = ""
    fit_file: str = "fit_stru.mat"
    contour_files: tuple[str, ...] = ("contourdata-x_B-y_Lp.csv", "contourdata-x_sigma-y_Lp.csv")
    extra_files: tuple[str, ...] = field(default_factory=tuple)

    @property
    def data_path(self) -> Path:
        return self.root / self.data_file

    @property
    def variant_root(self) -> Path:
        return self.root / self.subdir if self.subdir else self.root

    @property
    def fit_path(self) -> Path:
        return self.variant_root / self.fit_file

    def contour_path(self, name: str) -> Path:
        return self.variant_root / name

    def load_data(self):
        """Load and normalize the main MAT file for this run instance."""
        return _normalize_conductivity_measurements(loadmat(str(self.data_path)).get("data_stru"))

    def load_fit(self):
        """Load the saved fit bundle for this run instance."""
        return loadmat(str(self.fit_path)).get("fit_stru")

    def load_contour(self, name: str):
        """Load a contour CSV for this run instance."""
        return pd.read_csv(self.contour_path(name))


def _figure_output_base(name):
    # Make sure the shared figures folder exists before saving anything there.
    os.makedirs(FIGURES_DIR, exist_ok=True)
    return os.path.join(FIGURES_DIR, name)


def _data1_time_origin(data_stru):
    """Return the DATA1 notebook time origin used for all shifted plots."""
    return float(np.asarray(data_stru["data_raw"][0]["time"], dtype=float).reshape(-1)[0])


def _data1_shifted_time(data_stru, vial_index):
    """Return DATA1 vial time shifted by the shared notebook time origin."""
    return np.asarray(data_stru["data_raw"][vial_index]["time"], dtype=float) - _data1_time_origin(data_stru)


def _data1_shifted_time_minutes(data_stru, vial_index):
    """Return DATA1 vial time shifted and converted to minutes."""
    return _data1_shifted_time(data_stru, vial_index) / 60.0


def _parse_simple_salt_formula(salt):
    """Parse a simple binary salt formula into cation/anion labels.

    This helper is only used by the local DATA3 option-3 branch to derive a
    reasonable diffusivity from the selected workbook sheet.
    """
    salt = str(salt or "").strip()
    if not salt:
        return "", 1, "", 1
    match = re.fullmatch(r"([A-Z][a-z]?)(\d*)([A-Z][a-zA-Z0-9]*)(\d*)", salt)
    if not match:
        return salt, 1, "", 1
    cation, cation_n, anion, anion_n = match.groups()
    return cation, int(cation_n or 1), anion, int(anion_n or 1)


def _data3_diffusivity_cm2_s(namec):
    """Return a salt diffusivity for the local DATA3/custom workflow."""
    name = str(namec or "").strip()
    if not name:
        return 1.960e-5
    salt_map = {
        "K": 1.960e-5,
        "Na": 1.334e-5,
        "Li": 1.03e-5,
        "Mg": 0.706e-5,
        "Ca": 0.792e-5,
        "Co": 0.72e-5,
        "La": 0.62e-5,
    }
    letters = "".join(c for c in name if c.isalpha())
    for key, value in salt_map.items():
        if letters.startswith(key) or key in letters:
            return value
    cation, _, _, _ = _parse_simple_salt_formula(name)
    return salt_map.get(cation, 1.960e-5)


def _load_conductivity_paper():
    """Load the original paper conductivity model used by the DATA3 workflow."""
    import importlib.util
    from importlib.machinery import SourceFileLoader

    paper_path = (
        Path(__file__).resolve().parents[1]
        / "UnifiedFramework"
        / "DATA3"
        / "ExperimentalDataAnalysis"
        / "UnifiedCode"
        / "conductivity_paper.py"
    )
    if not paper_path.exists():
        raise FileNotFoundError(f"Conductivity paper source not found: {paper_path}")
    spec = importlib.util.spec_from_loader(
        "conductivity_paper",
        SourceFileLoader("conductivity_paper", str(paper_path)),
    )
    cp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cp)
    return cp


# Module-level toggle for the mass-litmus test. When True, solve_model pins
# sigma=0 and replaces the WSSE objective with a mass-only term, even if
# mass_litmus_test=False is passed explicitly. The runfile sets this for
# the diagnostic menu option; otherwise it stays False.
#
# What the test does
# ------------------
# Collapses the parameter estimation to its simplest possible form:
#   * sigma = 0 (no rejection) → flux equation reduces to Jw = Lp × ΔP
#   * objective = mass residuals only (drop cV and cF terms)
# This isolates the volumetric flux model from the concentration model and
# leaves a one-parameter convex fit (Lp from vial mass slopes). If the
# mass plots overlay cleanly under this test, the model + data + units
# pipeline is sound. If they don't, the issue is upstream of the optimizer.
NF270_MASS_LITMUS_TEST_ACTIVE = False


# Optional absolute floor on the cF residual scale (in mM). When set, the
# obj_rule replaces the relative 0.3 % cF weight `0.003 * cF_meas[i]` with
# `max(0.003 * cF_meas[i], NF270_CF_RESIDUAL_FLOOR_MM)`. This is a
# diagnostic for cases where the cF channel over-weights at high
# concentration (a 0.3 % weight at cF ~ 100 mM is 0.3 mM, tighter than
# the actual measurement uncertainty). Default None preserves legacy
# behavior exactly.
NF270_CF_RESIDUAL_FLOOR_MM = None


# DATA3-only: retentate (cF) and permeate (cV) concentration MEASUREMENT ERROR,
# expressed as a fraction of the measured concentration, applied CONSISTENTLY in
# the WSSE objective (obj_rule), the Fisher-information covariance (calc_FIM), and
# the Pyomo.parmest measurement_error suffix.  Keeping all three on one source
# avoids an inconsistent reported covariance.
#
# Basis — Lilonfe, Estrada, Singh, Ouimet, Phillip, Dowling, "Soft Sensors Enable
# Real-Time Ion Concentration Measurements", ChemRxiv 2026
# (doi:10.26434/chemrxiv.15002541/v1):
#   * cF (retentate) is conductivity / soft-sensor-derived.  End-to-end error =
#     soft-sensor MAPE (NaCl 1.0 %, CaCl2 2.2 %, LaCl3 1.9 %; Table 3) in
#     quadrature with the regressed cell-constant uncertainty (~1.38 %, Eq. S5)
#     -> ~2 %.  The legacy 0.3 % had NO physical basis and over-weighted cF at
#     high concentration (0.3 % of 100 mM = 0.3 mM, tighter than the probe
#     itself), which pushed sigma to its bound.  DATA3 default = 0.02 (2 %).
#   * cV (permeate) is an offline ICP-OES scalar — the REFERENCE.  ICP-OES
#     accuracy is ~1-3 %, so the legacy 3 % is already physical.  The paper's
#     "within 5 %" is the soft-sensor-vs-ICP AGREEMENT, not ICP's own error, so
#     it does NOT justify loosening the ICP channel.  DATA3 default = None (3 %).
#     (The 5 % soft-sensor figure applies to the conductivity permeate PROBE
#     channel cH, configured separately via NF270_USE_PERMEATE_PROBE.)
#
# GUARDED to workflow_family == "DATA3": when a knob is None, or for DATA1/DATA2,
# the legacy 0.003 / 0.03 weights are used unchanged (byte-identical published
# behavior).  See Architecture.md.
NF270_CF_RESIDUAL_SCALE_FRACTION = 0.02
NF270_CP_RESIDUAL_SCALE_FRACTION = None

# DATA3-workflow toggle. When True, the DATA3 measurement-error spec (cF 2% /
# cV 3%) is used regardless of the `workflow_family` string threaded through a
# given call — i.e. it can only PROMOTE to DATA3, never demote. Default False, so
# DATA1 and DATA2 are completely unaffected. This guards the DATA3 path against
# the fragility the port audit flagged (e.g. calc_FIM's default workflow_family
# is 'DATA1', so a DATA3 call that forgets to pass workflow_family='DATA3' would
# silently revert the retentate weight to the legacy 0.3%). The DATA3 workflow
# entry points flip this on FOR THE DURATION OF A RUN and restore it in a
# finally, so it never bleeds into a subsequent DATA1/DATA2 run. DATA1/DATA2
# never set it. Set it directly (with try/finally) or via run_with_data3_spec().
NF270_FORCE_DATA3_SPEC = False


def run_with_data3_spec(func, *args, **kwargs):
    """Call ``func(*args, **kwargs)`` with NF270_FORCE_DATA3_SPEC forced on, then
    restore the prior value — a safe, exception-proof way for a DATA3 workflow to
    guarantee the DATA3 error spec without leaking the toggle into later runs."""
    global NF270_FORCE_DATA3_SPEC
    _prior = NF270_FORCE_DATA3_SPEC
    NF270_FORCE_DATA3_SPEC = True
    try:
        return func(*args, **kwargs)
    finally:
        NF270_FORCE_DATA3_SPEC = _prior


def _nf270_conc_scales(workflow_family):
    """Return (cf_frac, cp_frac): the retentate/permeate concentration relative
    measurement-error fractions for the WLS weight 1/(frac*c).  DATA3 uses the
    paper-justified NF270_*_RESIDUAL_SCALE_FRACTION when set; DATA1/DATA2 (and
    unset DATA3 knobs) fall back to the legacy 0.003 / 0.03.  The DATA3 path is
    selected by workflow_family=='DATA3' OR the NF270_FORCE_DATA3_SPEC toggle
    (promote-only; never demotes DATA1/DATA2)."""
    is_data3 = (str(workflow_family).upper() == "DATA3") or NF270_FORCE_DATA3_SPEC
    cf = (float(NF270_CF_RESIDUAL_SCALE_FRACTION)
          if is_data3 and NF270_CF_RESIDUAL_SCALE_FRACTION is not None else 0.003)
    cp = (float(NF270_CP_RESIDUAL_SCALE_FRACTION)
          if is_data3 and NF270_CP_RESIDUAL_SCALE_FRACTION is not None else 0.03)
    return cf, cp


# Ionic strength  I = 1/2 * sum_j c_j z_j^2.  For a chloride salt M^{z+}Cl_z,
# electroneutrality (anion conc = z * cation conc, z_anion = 1) gives
#     I = 1/2 * c * (z^2 + z) = 1/2 * z * (z + 1) * c,
# i.e. I = k_I * c with k_I = {NaCl 1, CaCl2 3, LaCl3 6}.  k_I converts the
# (single-salt) cation/salt concentration to ionic strength and is the basis for
# re-parameterizing the solute permeability as B = f(I) instead of B = f(c) —
# because I = k_I*c is LINEAR, a single-salt B(I) fit is just a rescaling of the
# B(c) fit (beta_k -> beta_k / k_I^k), but the resulting beta become transferable
# ACROSS salts, which is what multisalt prediction needs.  See Architecture.md.
NF270_IONIC_STRENGTH_FACTOR = {"NaCl": 1.0, "KCl": 1.0, "CaCl2": 3.0, "LaCl3": 6.0}


def nf270_ionic_strength_factor(salt, cation_valence=None):
    """Return k_I such that ionic strength I = k_I * c for a 1:z chloride salt.

    Looks the salt up by name; otherwise falls back to 0.5*z*(z+1) computed from
    the cation valence z (which can be read from data_config as ni-1, since the
    van't Hoff count ni = 1 cation + z anions = z + 1 for MCl_z)."""
    if salt in NF270_IONIC_STRENGTH_FACTOR:
        return NF270_IONIC_STRENGTH_FACTOR[salt]
    if cation_valence is not None:
        z = float(cation_valence)
        return 0.5 * z * (z + 1.0)
    return 1.0


# Re-parameterize the concentration-dependent solute permeability as B = f(I)
# instead of B = f(c): inside the numeric-B_form polynomial, the membrane-
# interface concentration cIn that drives B (Js = B*(cIn - cH)) is replaced by
# the INTERFACIAL ionic strength I_in = k_I * cIn, with k_I from the salt
# valence (NaCl 1, CaCl2 3, LaCl3 6).  Because I = k_I*c is LINEAR, a single-
# salt B(I) fit is numerically identical to the B(c) fit after rescaling
# beta_k -> beta_k / k_I^(power of cIn that beta_k multiplies) — but the beta
# then become TRANSFERABLE across salts, which is the point for multisalt
# prediction.  GUARDED to workflow_family == "DATA3"; when None (default) or for
# DATA1/DATA2 the legacy B(c) constraint is used byte-identically regardless of
# value.  Set to True (or any truthy) to activate.  See Architecture.md §18.3.
NF270_B_USE_IONIC_STRENGTH = None


# Optional fix of the Donnan-dielectric steric*Born factor k.  The donnan B-form
# has three correlated parameters (P0, X, k); fixing k at a physical value turns
# the DAE fit into a well-conditioned 2-parameter (P0, X) problem.  When set to a
# float (DATA3 only), m.k_dd is a fixed Param instead of a free Var; when None,
# all three are estimated (legacy).  Used as a fallback retry by the campaign.
NF270_DONNAN_FIX_K = None


# Optional fix of the Donnan fixed-charge scale X [mM].  X carries the entire
# c-dependence of the Donnan B(c): X->0 collapses B to the constant Jw*P0*k
# (flat curve).  Fitting it freely on single-salt DATA3 rails X to its lower
# bound (1e-3), so this toggle pins X at a chosen value (DATA3 only) to PROFILE
# the objective cost of forcing the curve to bend.  When None, X is a free Var
# (legacy).  Mirrors NF270_DONNAN_FIX_K.
NF270_DONNAN_FIX_X = None


# Optional interior bound on sigma to prevent the optimizer from pinning at
# the {0, 1} corners on sheets where the sigma surface is genuinely flat
# (multivalent salts especially). When set to (lo, hi) with 0 <= lo < hi <= 1,
# the Var bounds change from (0, 1) to (lo, hi). When None, legacy bounds
# (0, 1) are preserved exactly.
NF270_SIGMA_INTERIOR_BOUNDS = None


# Per-salt upper bounds on the solute-permeability coefficient B, applied
# ONLY to DATA3 NF270 single-salt fits (workflow_family == "DATA3" and
# B_form == 'single').  Replaces the legacy single 50 upper with a
# valence-graded triple derived from dielectric exclusion + Stokes-Einstein
# scaling on a negatively-charged NF270 polyamide layer.  Salts not present
# in the dict fall back to NF270_B_BOUNDS_DEFAULT, so adding a new chloride
# salt is a one-line addition.  GUARDED to DATA3 — DATA1/DATA2 fits always
# see (1e-6, 50) regardless of this dict's contents, preserving byte-
# identical legacy behavior.  See Architecture.md §17.
NF270_B_BOUNDS_PER_SALT = {
    "NaCl":  (1e-6, 30.0),
    "CaCl2": (1e-6, 15.0),
    "LaCl3": (1e-6, 10.0),
}
NF270_B_BOUNDS_DEFAULT = (1e-6, 50.0)


# Adds the continuous permeate-probe conductivity (Shedlovsky-inverted to mM)
# as a new residual channel in the WSSE objective.  NF270 sheets carry per-vial
# cV_perm_cond arrays (5-Hz probe) that today only feed visualization; turning
# this on lets ~500 samples per vial constrain B and sigma instead of the lone
# ICP scalar.  GUARDED to workflow_family == "DATA3" so DATA1/DATA2 are byte-
# identical regardless of the value.  Set to a positive float to enable;
# the value is the relative residual scale (e.g. 0.03 for 3 %).  None = off.
NF270_USE_PERMEATE_PROBE = None
NF270_PERMEATE_PROBE_FLOOR_MM = 0.1     # absolute floor on permeate residual scale (mM)


# When truthy and workflow_family == "DATA3", add a SECOND end-of-experiment
# residual term to obj_cr using `data_config['cF_retentate_icp_mM']` as the
# target — alongside the existing `cF_final_meas` (Final Tube ICP) anchor.
# The two physical samples are not equivalent (Retentate = bulk stirred-cell
# at experiment end; Final Tube = connecting-tubing dead-volume sample), so
# enabling this constraint forces the fit to reconcile both.  Default None
# (off) preserves byte-equivalent behavior with the legacy code path and the
# DATA1/DATA2 paths regardless of value.  Set to True (or any truthy) to
# activate.  See PREFLIGHT_AUDIT.md row "Retentate row" — historically
# documented as a "cross-check only" measurement; this toggle promotes it to
# a fit constraint.
NF270_USE_RETENTATE_ICP_ANCHOR = None


# Seeded multistart: instead of LHS over the full (Lp, B, sigma) bounds, cluster
# starts around explicit seed thetas — typically the top-K grid points from the
# sheet's own contour panel, plus one optimal theta from a different salt's
# fit (cross-salt seeding).  Default off.  When True, NF270/DATA3 fits with
# multistart=True will use the seeded strategy instead of pure LHS.
NF270_MULTISTART_USE_CONTOUR_SEEDS = False
NF270_MULTISTART_CONTOUR_DIR = None      # if None, library auto-detects
# When seeded multistart is on, decide whether to include the cross-salt
# reference theta as one of the seed points.  Default True preserves the
# May-24 behavior.  Set False to use ONLY the per-sheet contour minima as
# seeds (pure contour-seeded mode, no cross-salt anchor).
NF270_MULTISTART_INCLUDE_CROSS_SALT = True
NF270_MULTISTART_CROSS_SALT_REFERENCES = {
    # Reference (Lp, B, sigma) used as the cross-salt seed in seeded
    # multistart.  Restored 2026-05-25 to match the recipe documented in
    # DATA3_single_salt_analysis_v7.pptx slide 29:
    #
    #   "When fitting CaCl₂ or LaCl₃: use MC3 SNaCl's interior-σ fit.
    #    (Lp = 8.22, B = 14.7, σ = 0.45) — the campaign's gold-standard θ."
    #
    # MC3.07.22.24_SNaCl is the ONLY sheet in the campaign with interior σ,
    # and slide 12 documents that all three channels (mass, retentate cF,
    # vial-ICP) agree on the same Lp ≈ 8 there.  Anchoring every cross-salt
    # seed on this θ gives the optimizer a known-physical starting point.
    #
    # Previous (rejected) anchor was MC2.NaCl's Lp=9.11, B=8.31, σ=0.66 —
    # but that fit's σ was tied to objective formulations that produced
    # σ=1.0 under today's recipe.  The MC3.SNaCl anchor is robust because
    # the dilution-regime data themselves DEMAND an interior σ.
    "NaCl":  {"Lp": 8.22, "B": 14.7, "sigma": 0.45},  # MC3.SNaCl gold standard
    "CaCl2": {"Lp": 8.22, "B": 14.7, "sigma": 0.45},  # MC3.SNaCl gold standard
    "LaCl3": {"Lp": 8.22, "B": 14.7, "sigma": 0.45},  # MC3.SNaCl gold standard
}


CONDUCTIVITY_SALT_PARAMS_25C = {
    # Simple built-in defaults for the paper data files.
    "KCl": {
        "epsilon": 78.4,
        "eta": 0.0089,
        "a": 1e-8,
        "z_1": 1,
        "z_2": -1,
        "lambda_0_cation": 73.50,
        "lambda_0_anion": 76.35,
        "lambda_0": 149.85,
    },
    "NaCl": {
        "epsilon": 78.4,
        "eta": 0.0089,
        "a": 1e-8,
        "z_1": 1,
        "z_2": -1,
        "lambda_0_cation": 50.11,
        "lambda_0_anion": 76.35,
        "lambda_0": 126.46,
    },
    "CaCl2": {
        "epsilon": 78.4,
        "eta": 0.0089,
        "a": 1e-8,
        "z_1": 2,
        "z_2": -1,
        "lambda_0_cation": 119.0,
        "lambda_0_anion": 76.35,
        "lambda_0": 195.35,
    },
    "LaCl3": {
        "epsilon": 78.4,
        "eta": 0.0089,
        "a": 1e-8,
        "z_1": 3,
        "z_2": -1,
        "lambda_0_cation": 208.8,
        "lambda_0_anion": 76.35,
        "lambda_0": 285.15,
    },
}

# Linear 25 C compensation factors used to normalize measured conductivity
# before it enters the 25 C-anchored inversion tables.
CONDUCTIVITY_TEMP_COEFF_PER_C = {
    "KCl": 0.019,
    "NaCl": 0.021,
    "CaCl2": 0.023,
    "LaCl3": 0.025,
}
CONDUCTIVITY_TEMP_COEFF_DEFAULT = 0.024

# DEPRECATED / UNUSED — permeate tube-transit shift (dead-volume model).
# The loader no longer applies this: per collaborator feedback we follow the
# DATA2 convention and anchor the single permeate ICP value at vial close
# (compared against the model cV at tau.last()), keeping mass/retentate/permeate
# on one shared time axis. This constant + `_nf270_corrected_cv_index` are kept
# only for the retired-model diagnostic in `_preflight_audit.py`.
# Original note: single apparatus constant for the MC2-MC5 UCB diafiltration rig,
# reverse-engineered from the DATA2 paper's time-corrected dataset 270511.123
# (cV_avg ~67-80 s before vial midpoint, avg flow ~5 mg/s → V_tube ≈ 0.3 g).
# Formula: t_corr = (t_open + t_close)/2 - V_tube / (dm/dt)_avg
NF270_TUBE_VOLUME_G = 0.3


def _ec25_compensate(cond_at_T_uS_per_cm, T_celsius, salt_name):
    """Rescale a conductivity reading taken at temperature T to what
    it WOULD have read at 25 °C.

    Why this matters
    ----------------
    A conductivity meter reads HIGHER when water is warm (ions move
    faster) and LOWER when cold — even though the salt concentration
    is identical. Our concentration look-up tables (the Shedlovsky /
    MSA models in conductivity_paper.py) are anchored at 25 °C.
    If the lab water is at 21 °C, the meter reads ~8 % low vs. the
    25 °C table, so the concentration we infer would be ~8 % wrong
    unless we correct first.

    The equation
    ------------
        sigma_25  =  sigma_T  *  ( 1  +  alpha * (25 - T) )
        ^^^^^^^^     ^^^^^^^      ^^^^^^^^^^^^^^^^^^^^^^^^^
        what we      what the     small bump-up factor:
        WANT (a      meter        alpha ≈ 2 %/°C, so at T = 21 °C
        25 °C        ACTUALLY     we bump up by ~8 %
        equivalent)  READ

    Sign of the correction:
        T < 25 °C → (25 − T) > 0 → bump UP    (meter read too low)
        T > 25 °C → (25 − T) < 0 → bump DOWN  (meter read too high)

    Reference
    ---------
    APHA Standard Methods 2510 B, "Conductivity: Laboratory Method".
    This is the small-ΔT linearization of  sigma_T = sigma_25 / [1 + α(T − 25)],
    accurate to <0.2 % within ±10 °C of 25 °C for dilute electrolytes.

    Parameters
    ----------
    cond_at_T_uS_per_cm : float or array-like
        Raw conductivity reading(s), in microsiemens per cm (μS/cm).
    T_celsius : float or array-like, same shape as cond
        Temperature(s) of the solution when the meter took the reading,
        in degrees Celsius.
    salt_name : str
        Salt identifier ("NaCl", "CaCl2", "LaCl3", "KCl"). Used to pick
        the salt-specific value of alpha from the lookup table at the top
        of this file.

    Returns
    -------
    Same shape as input — the 25 °C-equivalent conductivity, also μS/cm.
    """
    # Step 1: Convert inputs to NumPy arrays so we can do the math on
    # a whole column of measurements at once (vectorized arithmetic) instead
    # of a slow Python for-loop over thousands of timesteps.
    cond = np.asarray(cond_at_T_uS_per_cm, dtype=float)
    temp = np.asarray(T_celsius, dtype=float)

    # Step 2: Look up the per-degree correction "alpha" for our specific
    # salt. alpha ≈ 2 %/°C for most dilute salts, but it varies a little
    # because different ions move at slightly different rates:
    #   KCl ≈ 0.019/°C,  NaCl ≈ 0.021/°C,
    #   CaCl₂ ≈ 0.023/°C, LaCl₃ ≈ 0.025/°C.
    # If the salt isn't in our table, we fall back to a generic 0.024/°C.
    alpha = float(CONDUCTIVITY_TEMP_COEFF_PER_C.get(
        str(salt_name), CONDUCTIVITY_TEMP_COEFF_DEFAULT
    ))

    # Step 3: Apply the formula
    #     sigma_25  =  sigma_T  *  (1 + alpha * (25 - T))
    # NumPy broadcasts the scalar `alpha` across every entry, so this
    # one line does the bump-up correction for every measurement at once.
    # NaN entries propagate naturally — useful for sensor-dropout rows.
    compensated = cond * (1.0 + alpha * (25.0 - temp))

    # Step 4: Preserve the caller's input shape — scalar in → scalar out;
    # array in → array out — so downstream code doesn't have to unwrap.
    if np.isscalar(cond_at_T_uS_per_cm):
        return float(compensated)
    return compensated.astype(cond.dtype, copy=False)


def _nf270_corrected_cv_index(vial_time, vial_mass, V_tube_g=NF270_TUBE_VOLUME_G):
    """DEPRECATED / UNUSED by the loader. Find the timestamp inside a vial where
    the ICP-OES concentration measurement *would* belong under the retired
    tube-transit (dead-volume) model.

    The loader now follows the DATA2 convention and anchors the single permeate
    ICP value at vial close instead (see the per-vial build loop). This helper is
    retained only for the `_preflight_audit.py` retired-model diagnostic.

    Why this matters — the picture
    -------------------------------
    The setup looks like this:

         membrane  ─────tubing─────  vial
         (where               (where the drop
         permeate              lands and the ICP
         is generated)         later measures it)

    Permeate generated by the membrane at time t doesn't reach the vial
    until time (t + tau), where tau is the time it took to traverse the
    tube. Said the other way: a drop that LANDS in the vial at time t
    was actually generated by the membrane back at time (t − tau).

    The single ICP-OES reading represents the AVERAGE concentration of
    every drop that landed in that vial. If we plotted that ICP value
    at "vial close time" (t_close), we'd be ~tau seconds LATE relative
    to when the membrane physics actually produced the average.

    So we shift the measurement backward in time to anchor it at the
    membrane-event time, not the vial-arrival time.

    The equation
    ------------
        tau    =  V_tube  /  (dm/dt)_avg
              =  (mass of fluid trapped in the tube)  /  (avg mass flow rate)

        t_corr =  t_mid  −  tau
              =  geometric midpoint of the vial − transit delay

    where
        V_tube     ≈ 0.3 g (apparatus-specific dead volume, reverse-engineered
                            from the DATA2 paper experiments; see
                            NF270_TUBE_VOLUME_G at module top)
        (dm/dt)_avg = (total mass collected) / (vial duration)
        t_mid       = (t_open + t_close) / 2

    Example
    -------
    Vial opens at 270 s, closes at 455 s, collects 0.93 g.
    Avg flow = 0.93 / 185 = 5.03 mg/s  →  tau = 0.3 / 0.00503 = 60 s.
    Vial midpoint = (270 + 455) / 2 = 362.5 s.
    t_corr = 362.5 − 60 = 302.5 s (≈ 32 s into the vial, not 92 s into it).

    Returns
    -------
    int  : the index into vial_time CLOSEST to t_corr — that's where the
           single ICP value gets planted in the NaN-padded cV_avg array.
    None : if the vial is too short, too empty, or otherwise unusable
           (e.g. a startup vial that collected only droplets).

    Reference
    ---------
    Liu et al., "Time Correction for Permeate" section (DATA2 paper).
    """
    # Step 1: Make sure we got NumPy arrays, not Python lists.
    # .reshape(-1) flattens any (N,1) shape to (N,) so .size etc. behave.
    t = np.asarray(vial_time, dtype=float).reshape(-1)
    m = np.asarray(vial_mass, dtype=float).reshape(-1)

    # Step 2: Reject vials too short to do anything sensible with.
    # Need at least 2 timestamps to define a duration.
    if t.size < 2:
        return None
    finite = np.isfinite(m)
    if finite.sum() < 2:
        return None

    # Step 3: Total mass collected = last finite mass reading − first.
    # (We use the "finite" mask in case the balance had a NaN dropout at
    # the very start or end.)
    m_clean = m[finite]
    m_total = float(m_clean[-1] - m_clean[0])

    # Step 4: Reject the "startup vial" case. If the vial collected
    # less than 1 mg, it's basically empty (the pump was warming up,
    # or the vial was swapped before any real permeate dripped in).
    # The correction would divide by a near-zero flow rate and blow up.
    if m_total <= 1e-3:                  # < 1 mg collected — startup / empty vial
        return None

    # Step 5: Vial duration in seconds and average mass flow rate.
    duration = float(t[-1] - t[0])
    if duration <= 0:
        return None
    dmdt_avg = m_total / duration         # g/s — the (dm/dt)_avg in the formula
    if dmdt_avg <= 0:
        return None

    # Step 6: The actual physics: tau = V_tube / (dm/dt).
    # Units check: g / (g/s) = s. Good.
    tau = float(V_tube_g) / dmdt_avg

    # Step 7: t_corr = vial midpoint − transit delay.
    # Then clamp inside [t_open, t_close] so we never plant the
    # measurement BEFORE the vial actually opened (defensive — for
    # the rare case when tau exceeds vial-half-width).
    t_corr = (t[0] + t[-1]) / 2.0 - tau
    t_corr = max(float(t[0]), min(float(t[-1]), t_corr))   # clamp to vial bounds

    # Step 8: Return the INDEX of the timestamp closest to t_corr.
    # The caller uses this index to plant the single ICP value inside
    # a NaN-padded length-N array, so the fit/plot can interpolate the
    # simulated cV at the corrected timestamp.
    return int(np.argmin(np.abs(t - t_corr)))


def _row_temperature_series(row, *candidate_keys):
    """Return the first usable temperature series stored on a data row."""
    if not isinstance(row, dict):
        return None
    for key in candidate_keys:
        if key not in row:
            continue
        value = row.get(key)
        if value is None:
            continue
        temp = np.asarray(value, dtype=float).reshape(-1)
        if temp.size:
            return temp
    return None


def _next_temperature_slice(series, offset, signal_len):
    """Slice a flattened temperature trajectory for the next row, if present."""
    if not isinstance(series, np.ndarray) or series.size == 0 or signal_len <= 0:
        return None, offset
    next_offset = offset + signal_len
    if next_offset <= series.size:
        return series[offset:next_offset], next_offset
    return None, next_offset


def _apply_ec25_compensation_to_data_stru(data_stru):
    """Normalize conductivity measurements to 25 C when temperatures are present.

    This is the shared load-time hook for any conductivity-backed experimental
    file. If a dataset carries per-row temperature traces for the retentate or
    permeate stream, the corresponding conductivity series is pre-compensated to
    a 25 C-equivalent value before the existing 25 C-anchored inversion runs.
    """
    if not isinstance(data_stru, dict):
        return data_stru
    if not data_stru.get("conductivity_cF") or data_stru.get("conductivity_cF_converted"):
        return data_stru
    if data_stru.get("conductivity_temp_compensated"):
        return data_stru

    cfg = data_stru.get("data_config", {})
    salt_name = str(cfg.get("namec") or "").strip()
    if not salt_name:
        return data_stru

    rows = data_stru.get("data_raw") or []
    if not rows:
        return data_stru

    ret_temp_trajectory = np.asarray(data_stru.get("measured_ret_temp_C", []), dtype=float).reshape(-1)
    perm_temp_trajectory = np.asarray(data_stru.get("measured_perm_temp_C", []), dtype=float).reshape(-1)
    ret_offset = 0
    perm_offset = 0
    planned_rows = []
    saw_conductivity = False
    missing_temperature = False

    for row in rows:
        row_plan = {"row": row}
        ret_signal = row.get("cF_exp")
        if ret_signal is not None:
            saw_conductivity = True
            ret_signal_arr = np.asarray(ret_signal, dtype=float).reshape(-1)
            sig_len = ret_signal_arr.size
            ret_temp = _row_temperature_series(
                row,
                "retentate_temp",
                "retentate_temp_C",
                "measured_ret_temp_C",
            )
            if ret_temp is None:
                ret_temp, ret_offset = _next_temperature_slice(
                    ret_temp_trajectory,
                    ret_offset,
                    sig_len,
                )
            ret_offset += sig_len
            if ret_temp is None:
                missing_temperature = True
            row_plan["ret_signal"] = ret_signal
            row_plan["ret_temp"] = ret_temp
            row_plan["ret_signal_shape"] = np.asarray(ret_signal).shape
            row_plan["ret_temp_shape"] = np.asarray(ret_temp).shape if ret_temp is not None else None

        perm_signal = row.get("cV_perm_cond")
        if perm_signal is not None:
            saw_conductivity = True
            perm_signal_arr = np.asarray(perm_signal, dtype=float).reshape(-1)
            sig_len = perm_signal_arr.size
            perm_temp = _row_temperature_series(
                row,
                "permeate_temp",
                "permeate_temp_C",
                "measured_perm_temp_C",
            )
            if perm_temp is None:
                perm_temp, perm_offset = _next_temperature_slice(
                    perm_temp_trajectory,
                    perm_offset,
                    sig_len,
                )
            perm_offset += sig_len
            if perm_temp is None:
                missing_temperature = True
            row_plan["perm_signal"] = perm_signal
            row_plan["perm_temp"] = perm_temp
            row_plan["perm_signal_shape"] = np.asarray(perm_signal).shape
            row_plan["perm_temp_shape"] = np.asarray(perm_temp).shape if perm_temp is not None else None

        planned_rows.append(row_plan)

    if not saw_conductivity or missing_temperature:
        return data_stru

    compensated_any = False
    for row_plan in planned_rows:
        row = row_plan["row"]
        ret_signal = row_plan.get("ret_signal")
        ret_temp = row_plan.get("ret_temp")
        if ret_signal is not None and ret_temp is not None:
            if "cF_exp_conductivity" not in row:
                row["cF_exp_conductivity"] = copy.deepcopy(ret_signal)
            if "cF_exp_conductivity_temp_C" not in row:
                row["cF_exp_conductivity_temp_C"] = copy.deepcopy(ret_temp)
            compensated_ret = _ec25_compensate(ret_signal, ret_temp, salt_name)
            ret_shape = row_plan.get("ret_signal_shape")
            if ret_shape is not None and ret_shape == ():
                row["cF_exp"] = float(np.asarray(compensated_ret).reshape(-1)[0])
            elif isinstance(ret_signal, list):
                row["cF_exp"] = np.asarray(compensated_ret, dtype=float).tolist()
            else:
                row["cF_exp"] = np.asarray(compensated_ret, dtype=float).reshape(ret_shape)
            compensated_any = True

        perm_signal = row_plan.get("perm_signal")
        perm_temp = row_plan.get("perm_temp")
        if perm_signal is not None and perm_temp is not None:
            if "cV_perm_cond_conductivity" not in row:
                row["cV_perm_cond_conductivity"] = copy.deepcopy(perm_signal)
            if "cV_perm_cond_conductivity_temp_C" not in row:
                row["cV_perm_cond_conductivity_temp_C"] = copy.deepcopy(perm_temp)
            compensated_perm = _ec25_compensate(perm_signal, perm_temp, salt_name)
            perm_shape = row_plan.get("perm_signal_shape")
            if perm_shape is not None and perm_shape == ():
                row["cV_perm_cond"] = float(np.asarray(compensated_perm).reshape(-1)[0])
            elif isinstance(perm_signal, list):
                row["cV_perm_cond"] = np.asarray(compensated_perm, dtype=float).tolist()
            else:
                row["cV_perm_cond"] = np.asarray(compensated_perm, dtype=float).reshape(perm_shape)
            compensated_any = True

    if compensated_any:
        cfg["Temp"] = 298.15
        data_stru["conductivity_temp_compensated"] = True
        data_stru["conductivity_temp_compensation"] = {
            "salt_name": salt_name,
            "alpha_per_C": float(CONDUCTIVITY_TEMP_COEFF_PER_C.get(
                salt_name, CONDUCTIVITY_TEMP_COEFF_DEFAULT
            )),
            "reference_C": 25.0,
            "method": "sigma_25 = sigma_T * (1 + alpha * (25 - T_celsius))",
        }

    return data_stru


def loadmat(filename):# for fun!
    '''
    Read in nested structure(mat file) generated from MATLAB and output dictionaries.

    This function is called instead of using spio.loadmat directly as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries which are still mat-objects.
    Adapted from https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
    
    Arguments:
        filename: str, location + filename
    
    Returns:
        dictionary contains structured data
    '''

    print("\nLoading data file =",filename,"\n")

    def _coerce_scalar(value):
        """Convert MATLAB numeric scalars to plain Python numbers.

        MATLAB MAT files frequently encode small integers as numpy scalar types
        such as uint8. Those types behave badly later in the workflow when we
        subtract times or build negative offsets, so we normalize them here.
        """
        if isinstance(value, np.generic):
            if np.issubdtype(type(value), np.number):
                return float(value)
            return value.item()
        return value

    def _check_keys(d):
        '''
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        '''
        for key in d:
            if isinstance(d[key], spio.matlab.mat_struct):
                d[key] = _todict(d[key])
            elif isinstance(d[key], np.ndarray):
                d[key] = _tolist(d[key])

        return d

    def _todict(matobj):
        '''
        A recursive function which constructs from matobjects nested dictionaries
        '''
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, spio.matlab.mat_struct):
                d[strg] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                d[strg] = _tolist(elem)
            else:
                d[strg] = _coerce_scalar(elem)
        return d

    def _tolist(ndarray):
        '''
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        '''
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, spio.matlab.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(_tolist(sub_elem))
            elif isinstance(sub_elem, (int, float, np.integer, np.floating, np.bool_)):
                elem_list.append(_coerce_scalar(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)

    return _check_keys(data)


def loadxlsx(filename, *, sheet=None):
    '''
    Read in an NF270 / DATA3 experimental Excel workbook into the same
    nested-dict shape that loadmat() produces for DATA1 / DATA2 .mat files.

    The DATA3 workflow is now syntactically identical to DATA1 / DATA2:

        from refactored_ucb_library import *

        # DATA1 / DATA2 (.mat):
        data_stru = loadmat('data_stru-dataset270511.123.mat')['data_stru']

        # DATA3 / NF270 (.xlsx):
        data_stru = loadxlsx('NF270_MC2.xlsx', sheet='05.07.24_NaCl')['data_stru']

        sim_stru = []
        plot_sim_comparison(data_stru, sim_stru, plot_pred=False, lg=True)

        mode = 'Lag'
        fit_stru, sim_stru, sim_inter = solve_model(
            data_stru, mode, sim_opt=False, B_form='single')
        plot_sim_comparison(data_stru, sim_stru, plot_pred=True)

    The returned data_stru carries every field DATA2 .mat data_stru carries:
    data_config (n, n_v0, n_extra, n_h, n_A, M_F0, M_O, C_F0, C_D, namec,
    ni, nc, nr, delP, Temp, Am, rho, Lp0, B0, sigma0, theta0) and data_raw
    (per-vial dicts with time, mass [reset to start at 0g], cF_exp [mM,
    conductivity-derived], cV_avg [mM, ICP-OES scalar]). The vial-1
    holdup-fill period is automatically split into a leading "extra"
    startup vial (n_v0=2, n_extra=1) so the model integrates Jw across
    the holdup but the objective skips its residuals — same convention
    DATA1 / DATA2 .mat already use.

    Arguments
    ---------
    filename : str or Path
        Path to the .xlsx file (e.g. "NF270_MC2.xlsx").
    sheet : str, optional
        Sheet name to load (e.g. "05.07.24_NaCl"). If None, loads the
        first sheet in the workbook.

    Returns
    -------
    dict
        {'data_stru': data_stru} — wrapped to match loadmat()'s call
        pattern so calling code reads identically across DATA1, DATA2,
        and DATA3 workflows.
    '''
    print("\nLoading XLSX file =", filename, "\n")
    if sheet is not None:
        print("  sheet =", sheet)
    data_stru = _load_legacy_data_stru_from_excel(filename, selector=sheet)
    data_stru["source_format"] = "xlsx"
    return {"data_stru": data_stru}


def plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=False,cond=True,preface=False):
    '''
    Plot simulation results comparing with measurements
    '''

    t_delay = _data1_time_origin(data_stru)
    figures = []

    # plot mass data/prediction comparison
    fig = plt.figure(figsize=(4,4))
    n_v0 = int(data_stru['data_config'].get('n_v0', 1))
    for i in range(data_stru['data_config']['n']):
        time = np.asarray(data_stru['data_raw'][i]['time'], dtype=float)
        mass = np.asarray(data_stru['data_raw'][i]['mass'], dtype=float)
        plt.plot((time-t_delay)/60, mass, 'r.', markersize=4)
        if plot_pred and i >= n_v0 - 1:
            plt.plot((np.asarray(sim_stru[i]['time'], dtype=float)-t_delay)/60,
                     np.asarray(sim_stru[i]['mV'], dtype=float),
                     'b', linewidth=3, alpha=.6)

    plt.plot([],[],'r.',markersize=4,label='Measurements')
    if plot_pred:
        plt.plot([],[],'b',linewidth=3,alpha=.6,label='Predictions')

    if preface:
        plt.xlabel('Time',fontsize=24,fontweight='bold')
        plt.ylabel('Mass',fontsize=24,fontweight='bold')
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    else:
        plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
        plt.ylabel('Mass [g]',fontsize=16,fontweight='bold')

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    if xticks:
        xticks[0].label1.set_visible(False)
    if lg:
        plt.legend(fontsize=10,loc='best')
    if LOUD:
        fname = _figure_output_base('mass-dat'+str(data_stru['dataset']))
        fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
    figures.append(fig)

    # plot concentration data/prediction comparison
    fig = plt.figure(figsize=(4,4))
    first_time = _data1_time_origin(data_stru)
    try:
        is_data2_workflow = float(data_stru.get("dataset", 0)) > 1_000
    except Exception:
        is_data2_workflow = False
    plt.plot((first_time-t_delay)/60,
             np.asarray(sim_stru[0]['cF'], dtype=float)[0],
             'ms',markersize=8,clip_on=False)

    plt.plot([],[],'ms',markersize=8,clip_on=False,label='Retentate (Measurement)')
    if plot_pred:
        plt.plot([],[],'g',linewidth=3,label='Retentate (Prediction)')
    plt.plot([],[],'cs',markersize=8,label='Vial (Measurement)')
    if plot_pred:
        plt.plot([],[],'r^',markersize=8,label='Vial (Prediction)')
        plt.plot([],[],'r-',linewidth=3,alpha=.6,label='Permeate (Prediction)')

    n_v0 = int(data_stru['data_config'].get('n_v0', 1))
    for i in range(data_stru['data_config']['n']):
        time = np.asarray(data_stru['data_raw'][i]['time'], dtype=float)
        c_v_avg = np.asarray(data_stru['data_raw'][i]['cV_avg'], dtype=float)
        if (not is_data2_workflow) or c_v_avg.ndim == 0 or c_v_avg.size == 1:
            plt.plot((time[-1]-t_delay)/60,
                     float(c_v_avg.reshape(-1)[0]),
                     'cs',markersize=8)
        else:
            n = min(len(time), len(c_v_avg))
            plt.plot((time[:n]-t_delay)/60,
                     c_v_avg[:n],
                     'cs',markersize=8,clip_on=False)
        if cond:
            cF_exp = data_stru['data_raw'][i]['cF_exp']
            if is_data2_workflow:
                if isinstance(cF_exp, float):
                    plt.plot((time[-1]-t_delay)/60, cF_exp,'ms',markersize=8)
                elif len(cF_exp)>1:
                    plt.plot((time-t_delay)/60,
                             np.asarray(cF_exp, dtype=float),'ms',markersize=8)
            else:
                if isinstance(cF_exp, float):
                    plt.plot((time[-1]-t_delay)/60, cF_exp,'ms',markersize=8)
                elif len(cF_exp)>1:
                    plt.plot((time-t_delay)/60,
                             np.asarray(cF_exp, dtype=float),'ms',markersize=8)
        if plot_pred and i >= n_v0 - 1:
            sim_time = np.asarray(sim_stru[i]['time'], dtype=float)
            plt.plot((sim_time-t_delay)/60,
                     np.asarray(sim_stru[i]['cF'], dtype=float),
                     'g',linewidth=3,alpha=.6)
            plt.plot((sim_time-t_delay)/60,
                     np.asarray(sim_stru[i]['cH'], dtype=float),
                     'r-',linewidth=3,alpha=.6)
            if not preface:
                sim_c_v = np.asarray(sim_stru[i]['cV'], dtype=float)
                if (not is_data2_workflow) or c_v_avg.ndim == 0 or c_v_avg.size == 1:
                    plt.plot((sim_time[-1]-t_delay)/60,
                             sim_c_v[-1],
                             'r^',markersize=8,alpha=.6)
                else:
                    valid = ~np.isnan(c_v_avg[: min(len(time), len(c_v_avg))])
                    if np.any(valid):
                        time_shifted = time[: min(len(time), len(c_v_avg))] - t_delay
                        interp_cv = interpolate.interp1d(sim_time, sim_c_v, fill_value='extrapolate')
                        plt.plot(time_shifted[valid]/60,
                                 interp_cv(time_shifted[valid]),
                                 'r^',markersize=8,alpha=.6)

    if not cond:
        plt.plot((np.asarray(data_stru['data_raw'][-1]['time'], dtype=float)[-1]-t_delay)/60,
                 np.asarray(data_stru['data_raw'][-1]['cF_exp'], dtype=float),
                 'ms',markersize=8)

    if preface:
        plt.xlabel('Time',fontsize=24,fontweight='bold')
        plt.ylabel('Concentration',fontsize=24,fontweight='bold')
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    else:
        plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
        plt.ylabel('Concentration [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in",top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    if xticks:
        xticks[0].label1.set_visible(False)
    if lg:
        plt.legend(fontsize=10,loc='best')
    if LOUD:
        fname = _figure_output_base('concentration-dat'+str(data_stru['dataset']))
        fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
    figures.append(fig)

    # plot mass of stirred cell
    if stirc_mass:
        fig = plt.figure(figsize=(4,4))
        for i in range(data_stru['data_config']['n']):
            if plot_pred:
                plt.plot((np.asarray(sim_stru[i]['time'], dtype=float)-t_delay)/60,
                         np.asarray(sim_stru[i]['mF'], dtype=float),
                         'b',linewidth=3,alpha=.6)

        plt.plot([],[],'r.',markersize=4,label='Mass (Measurements)')
        if plot_pred:
            plt.plot([],[],'b',linewidth=3,alpha=.6,label='Mass (Predictions)')

        plt.plot([],[],'ms',markersize=8,clip_on=False,label='Retentate (Measurement)')
        if plot_pred:
            plt.plot([],[],'g',linewidth=3,label='Retentate (Prediction)')
        plt.plot([],[],'cs',markersize=8,label='Vial (Measurement)')
        if plot_pred:
            plt.plot([],[],'r^',markersize=8,label='Vial (Prediction)')
            plt.plot([],[],'r-',linewidth=3,alpha=.6,label='Permeate (Prediction)')

        plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
        plt.ylabel('Mass in Stirred Cell [g]',fontsize=16,fontweight='bold')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tick_params(direction="in")
        ytop = max([max(np.asarray(sim_stru[i]['mF'], dtype=float)) for i in range(data_stru['data_config']['n'])])
        plt.ylim(bottom=0,top=ytop+0.5)
        if lg:
            plt.legend(fontsize=12.5,loc='best')
        if LOUD:
            fname = _figure_output_base('stirc_mass-dat'+str(data_stru['dataset']))
            fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
    figures.append(fig)

    return figures


def _plot_sim_comparison_data1_legacy(data_stru, sim_stru, plot_pred=True, lg=False, LOUD=False):
    """Mirror the notebook Figure 2 plotting behavior from diafiltration_plots.py."""
    t_delay = float(data_stru['data_raw'][0]['time'][0])
    figures = []

    fig = plt.figure(figsize=(4, 4))
    for i in range(data_stru['data_config']['n']):
        time = np.asarray(data_stru['data_raw'][i]['time'], dtype=float)
        mass = np.asarray(data_stru['data_raw'][i]['mass'], dtype=float)
        plt.plot((time - t_delay) / 60, mass, 'r.', markersize=4)
        if plot_pred:
            sim_time = np.asarray(sim_stru[i]['time'], dtype=float)
            plt.plot((sim_time - t_delay) / 60,
                     np.asarray(sim_stru[i]['mV'], dtype=float),
                     'b', linewidth=3, alpha=.6)

    plt.xlabel('Time [min]', fontsize=16, fontweight='bold')
    plt.ylabel('Mass [g]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    if xticks:
        xticks[0].label1.set_visible(False)
    if lg:
        plt.legend(fontsize=10, loc='best')
    if LOUD:
        fname = _figure_output_base('mass-dat' + str(data_stru['dataset']))
        fig.savefig(fname + '.png', dpi=300, bbox_inches='tight')
    figures.append(fig)

    fig = plt.figure(figsize=(4, 4))
    plt.plot((float(np.asarray(data_stru['data_raw'][0]['time'], dtype=float)[0]) - t_delay) / 60,
             float(np.asarray(sim_stru[0]['cF'], dtype=float)[0]),
             'ms', markersize=8, clip_on=False)
    for i in range(data_stru['data_config']['n']):
        time = np.asarray(data_stru['data_raw'][i]['time'], dtype=float)
        c_v_avg = data_stru['data_raw'][i]['cV_avg']
        if not isinstance(c_v_avg, list):
            plt.plot((float(time[-1]) - t_delay) / 60,
                     float(c_v_avg),
                     'cs', markersize=8)
        else:
            plt.plot((time - t_delay) / 60,
                     np.asarray(c_v_avg, dtype=float),
                     'cs', markersize=8, clip_on=False)

        c_f_exp = data_stru['data_raw'][i]['cF_exp']
        c_f_exp_arr = np.asarray(c_f_exp, dtype=float)
        if c_f_exp_arr.ndim == 0 or c_f_exp_arr.size == 1:
            plt.plot((float(time[-1]) - t_delay) / 60,
                     float(c_f_exp_arr.reshape(-1)[0]),
                     'ms', markersize=8, clip_on=False)
        else:
            plt.plot((time - t_delay) / 60,
                     c_f_exp_arr,
                     'ms', markersize=8, clip_on=False)

        if plot_pred:
            sim_time = np.asarray(sim_stru[i]['time'], dtype=float)
            sim_c_f = np.asarray(sim_stru[i]['cF'], dtype=float)
            sim_c_h = np.asarray(sim_stru[i]['cH'], dtype=float)
            sim_c_v = np.asarray(sim_stru[i]['cV'], dtype=float)
            plt.plot((sim_time - t_delay) / 60, sim_c_f, 'g', linewidth=3, alpha=.6)
            plt.plot((sim_time - t_delay) / 60, sim_c_h, 'r-', linewidth=3, alpha=.6)
            if not isinstance(c_v_avg, list):
                plt.plot((float(sim_time[-1]) - t_delay) / 60,
                         float(sim_c_v[-1]),
                         'r^', markersize=8, alpha=.6)
            else:
                c_v_avg_arr = np.asarray(c_v_avg, dtype=float)
                n_v = min(len(time), len(c_v_avg_arr))
                valid_v = ~np.isnan(c_v_avg_arr[:n_v])
                interp_cv = interpolate.interp1d(sim_time, sim_c_v, fill_value='extrapolate')
                shifted_time = time[:n_v] - t_delay
                plt.plot(shifted_time[valid_v] / 60,
                         interp_cv(shifted_time[valid_v]),
                         'r^', markersize=8, alpha=.6)

    plt.xlabel('Time [min]', fontsize=16, fontweight='bold')
    plt.ylabel('Concentration [mM]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in", top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    if xticks:
        xticks[0].label1.set_visible(False)
    if lg:
        plt.legend(fontsize=10, loc='best')
    if LOUD:
        fname = _figure_output_base('concentration-dat' + str(data_stru['dataset']))
        fig.savefig(fname + '.png', dpi=300, bbox_inches='tight')
    figures.append(fig)
    return figures


def _resolve_data1_figure2_root(data_root):
    """Prefer the notebook-style DATA1 data_library when present for Figure 2."""
    root = _resolve_data_root(data_root)
    if root.name == "data":
        sibling = root.parent / "data_library"
        if sibling.exists():
            return sibling
    return root


def run_data1_figure2_workflow(data_root=None, save_dir=None):
    """Recreate DATA1 Figure 2 using the notebook's legacy plot_sim_comparison path."""
    root = _resolve_data1_figure2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root
    save_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        ("501.1", root / "data_stru-dataset501.1.mat", root / "501.1 concpolar" / "fit_stru.mat"),
        ("511.12", root / "data_stru-dataset511.12.mat", root / "511.12 concpolar" / "fit_stru.mat"),
    ]

    outputs = []
    panel_paths = []
    for dat, data_path, fit_path in cases:
        if not data_path.exists() or not fit_path.exists():
            continue
        data_stru = loadmat(str(data_path)).get("data_stru")
        fit_stru = loadmat(str(fit_path)).get("fit_stru")
        sim_stru = fit_stru.get("sim_stru") if isinstance(fit_stru, dict) else None
        if sim_stru is None:
            continue
        figs = _plot_sim_comparison_data1_legacy(data_stru, sim_stru, plot_pred=True, lg=False, LOUD=False)
        mass_out = save_dir / f"mass-dat{dat}.png"
        conc_out = save_dir / f"concentration-dat{dat}.png"
        figs[0].savefig(mass_out, dpi=300, bbox_inches="tight")
        figs[1].savefig(conc_out, dpi=300, bbox_inches="tight")
        outputs.extend([str(mass_out), str(conc_out)])
        panel_paths.extend([mass_out, conc_out])
        plt.close(figs[0])
        plt.close(figs[1])

    try:
        from PIL import Image, ImageDraw, ImageFont

        ordered = [
            save_dir / "mass-dat501.1.png",
            save_dir / "concentration-dat501.1.png",
            save_dir / "mass-dat511.12.png",
            save_dir / "concentration-dat511.12.png",
        ]
        if all(path.exists() for path in ordered):
            imgs = [Image.open(path).convert("RGB") for path in ordered]
            target_w = max(im.width for im in imgs)
            target_h = max(im.height for im in imgs)
            resized = [im.resize((target_w, target_h), Image.Resampling.LANCZOS) for im in imgs]
            gap = 30
            label_pad = 35
            canvas = Image.new("RGB", (2 * target_w + 3 * gap, 2 * target_h + 3 * gap + 2 * label_pad), "white")
            draw = ImageDraw.Draw(canvas)
            try:
                font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 28)
            except Exception:
                font = None
            positions = [
                (gap, gap + label_pad),
                (2 * gap + target_w, gap + label_pad),
                (gap, 2 * gap + target_h + 2 * label_pad),
                (2 * gap + target_w, 2 * gap + target_h + 2 * label_pad),
            ]
            for label, img, pos in zip(["A", "B", "C", "D"], resized, positions):
                draw.text((pos[0], pos[1] - label_pad), label, fill="black", font=font)
                canvas.paste(img, pos)
            composite = save_dir / "figure_2.png"
            canvas.save(composite)
            outputs.append(str(composite))
    except Exception:
        pass

    return outputs


def save_model(m, LO=True, B_form='single', LOUD=False):
    '''
    Save model results
    
    Arguments:
        m: pyomo model instance
        LO: boolean, if mode is lag/overflow
        B_form: float/str, different solute permeability coefficient B formula 
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration              
        LOUD: boolean, if print out parameter results
    
    Returns:
        fit_stru: dict, parameter fit results
        sim_stru: dict, model predictions
    '''
    fit_stru = dict()
    fit_stru['parameters'] = dict()
    fit_stru['sideparameters'] = dict()
    
    fit_stru['parameters']['Lp'] = value(m.Lp) 
    
    if B_form=='single':
        fit_stru['parameters']['B'] = value(m.B)
    elif B_form=='pervial':
        fit_stru['parameters']['B'] = dict()
        for i in m.n_vial:
            fit_stru['parameters']['B'][i] = value(m.B[i])
    elif B_form=='convection':   
        fit_stru['parameters']['beta_0'] = value(m.beta_0)
        fit_stru['parameters']['beta_1'] = value(m.beta_1)
    elif isinstance(B_form, str) and B_form == 'sat':
        fit_stru['parameters']['B_inf'] = value(m.B_inf)
        fit_stru['parameters']['c_star'] = value(m.c_star)
    elif isinstance(B_form, str) and B_form == 'donnan':
        fit_stru['parameters']['P0'] = value(m.P0)
        fit_stru['parameters']['X'] = value(m.X)
        fit_stru['parameters']['k_dd'] = value(m.k_dd)
    else:
        fit_stru['parameters']['beta_0'] = value(m.beta_0)
        if not isinstance(B_form, str) and B_form != 0:
            fit_stru['parameters']['beta_1'] = value(m.beta_1)
        if not isinstance(B_form, str) and B_form > 1:
            fit_stru['parameters']['beta_2'] = value(m.beta_2)
        if not isinstance(B_form, str) and B_form > 2:
            fit_stru['parameters']['beta_3'] = value(m.beta_3)
            
    fit_stru['parameters']['sigma'] = value(m.sigma)
    if LO:
        fit_stru['parameters']['S0'] = value(m.S0)
        fit_stru['parameters']['S'] = value(m.S)   
    fit_stru['Obj'] = value(m.Obj)
    fit_stru['obj_m'] = value(m.obj_m)
    fit_stru['obj_cv'] = value(m.obj_cv)
    fit_stru['obj_cr'] = value(m.obj_cr)
    # New optional permeate-probe channel — defaults to 0 when toggle is off.
    if hasattr(m, "obj_cv_perm"):
        try:
            fit_stru['obj_cv_perm'] = float(value(m.obj_cv_perm))
        except Exception:
            fit_stru['obj_cv_perm'] = 0.0
    else:
        fit_stru['obj_cv_perm'] = 0.0
    fit_stru['Obj_tru'] = value(m.obj_tru)
    fit_stru['obj_cr_tru'] = value(m.obj_cr_tru)
    fit_stru['res_std'] = {'res_m': [value(res) for res in m.res_m],
                 'res_cp': [value(res) for res in m.res_cp],
                 'res_cf': [value(res) for res in m.res_cf]}   
        
    if LOUD:
        # print parameters
        print('Lp = ',value(m.Lp),' L / m / m / hr / bar')
        if isinstance(B_form, str):
            if B_form=='single':
                print('B = ',value(m.B),' micrometers / s')
            elif B_form=='pervial':    
                print('B = ',[value(m.B[i]) for i in m.n_vial],' micrometers / s')
            elif B_form=='convection':   
                print('D/l = ',value(m.beta_0))
                print('H = ',value(m.beta_1))
            elif B_form=='sat':
                print('B = Jw*B_inf*(1 - exp(-c_in/c_star)); B_inf =',value(m.B_inf),' c_star =',value(m.c_star),' mM')
            elif B_form=='donnan':
                print('B = Jw*P0*(-X+sqrt(X^2+4k^2 c_in^2))/(2 c_in); P0 =',value(m.P0),' X =',value(m.X),' mM  k =',value(m.k_dd))
        else:
            if B_form > 2:
                print('B = Jw*[',value(m.beta_0),'+',value(m.beta_1),'* c_in ^',B_form-2,'+',value(m.beta_2),'* c_in ^',B_form-1,'+',value(m.beta_3),'* c_in ^',B_form, ']')
            elif B_form > 1:
                print('B = Jw*[',value(m.beta_0),'+',value(m.beta_1),'* c_in ^',B_form-1,'+',value(m.beta_2),'* c_in ^',B_form,']')
            elif B_form == 0:
                print('B = Jw*[',value(m.beta_0), ']')
            else:
                print('B = Jw*[',value(m.beta_0),'+',value(m.beta_1),'* c_in ^',B_form,']')

        print('sigma = ',value(m.sigma),' dimensionless')
        print('cD = ',value(m.cD),' mM')

        if LO:
            print('S0 = ',value(m.S0),' g / s')
            print('S = ',value(m.S),' g / hr')
        print('Obj = ' ,value(m.Obj))
        print('Obj(m) = ' ,value(m.obj_m))
        print('Obj(cv) = ' ,value(m.obj_cv))
        print('Obj(cr) = ' ,value(m.obj_cr))
        if LO:
            print('Obj(truncate) = ' ,value(m.obj_tru))
            print('Obj(cr_truncate) = ' ,value(m.obj_cr_tru))
            print('count(cr0) = ' ,value(m.count_cr0))
        print('count = ' ,value(m.count))
        print('count(m) = ' ,value(m.count_m))
        print('count(cv) = ' ,value(m.count_cv))
        print('count(cr) = ' ,value(m.count_cr))
        print('likelihood1 = ' ,value(m.llh1))
        print('likelihood2 = ' ,value(m.llh2)) 
    
    sim_stru = dict()
    
    for n_vial in m.n_vial:
        time = [t * (m.tf[n_vial]-m.ti[n_vial]) + m.ti[n_vial] for t in m.tau]
        cF = [value(m.cF[n_vial,t]) for t in m.tau]
        cIn = [value(m.cIn[n_vial,t]) for t in m.tau]
        cH = [value(m.cH[n_vial,t]) for t in m.tau]
        mV = [value(m.mV[n_vial,t]) for t in m.tau]
        cV = [value(m.cV[n_vial,t]) for t in m.tau]
        Jw = [value(m.Jw[n_vial,t]) for t in m.tau]
        Js = [value(m.Js[n_vial,t]) for t in m.tau]
        
        sim_stru[n_vial-1] = {'time': time,
                'cF': cF,
                'cIn': cIn,
                'cH': cH,
                'mV': mV,
                'cV': cV,
                'Jw': Jw,
                'Js': Js}
        if LO:
            mF = [value(m.mF[n_vial,t]) for t in m.tau]
            sim_stru[n_vial-1]['mF'] = mF
        
        if isinstance(B_form, str) and 'convection' in B_form:
            B = None       
        elif B_form=='single':
            B = [value(m.B) for t in m.tau]
        elif B_form=='pervial':
            B = [value(m.B[n_vial]) for t in m.tau]
        else:
            B = [value(m.B[n_vial,t]) for t in m.tau]        
        sim_stru[n_vial-1]['B'] = B

    return fit_stru, sim_stru


def inter_model(data_stru, sim_stru, fit_stru):
    '''
    Interpolate model predictions at same time stamps of experimental measurements
    
    Arguments:
        data_stru: dict, experimental data dictionary
        sim_stru: dict, model predictions
        fit_stru: dict, parameter fit results
    
    Returns:
        sim_inter: dict, model predictions for experimental measurements 
    '''
    sim_inter = dict()
    t_delay = _data1_time_origin(data_stru)
    for i in range(data_stru['data_config']['n']):
        t_meas = _data1_shifted_time(data_stru, i)
        f_m = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['mV'], kind='linear', fill_value='extrapolate')
        if i >= data_stru['data_config']['n_v0']-1:            
            ind_m = ~np.isnan(data_stru['data_raw'][i]['mass'])
            t_m = [t_meas[i] for i, val in enumerate(ind_m) if val]
        else:
            t_m=[]
        if np.ndim(data_stru['data_raw'][i]['cV_avg']) > 0:
            f_cv = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['cV'], kind='linear', fill_value='extrapolate')
            ind_cv = ~np.isnan(data_stru['data_raw'][i]['cV_avg'])
            t_cv = [t_meas[i] for i, val in enumerate(ind_cv) if val]
            cV = f_cv(t_cv)
        else:
            cV = sim_stru[i]['cV'][-1]
        f_cf = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['cF'], kind='linear', fill_value='extrapolate')
        ind_cf = ~np.isnan(data_stru['data_raw'][i]['cF_exp'])
        t_cf = [t_meas[i] for i, val in enumerate(ind_cf) if val]
        sim_inter[i] = {'time': t_meas,
                        'mV': f_m(t_m),
                        'cF': f_cf(t_cf),
                        'cV': cV}       

    return sim_inter


def model_construct_inter(data_stru, mode, theta=None, sim_opt=False, B_form='single', workflow_family='DATA1'):
    '''
    Build the diafiltration ODE system as a Pyomo nonlinear program.

    The physical setup
    ------------------
    A stirred cell (volume mF, concentration cF) is pressurized to ΔP.
    Solution leaks through a membrane area Am into a vial (mass mV,
    concentration cV). Solute can be rejected by the membrane (σ ≠ 0)
    and/or driven through it by diffusion (B). For DIAFILTRATION specifically
    we *also* feed fresh diafiltrate of concentration cD into the cell so
    the cell volume stays approximately constant.

          ┌──────────────────────┐
          │   Stirred cell       │
          │   mass = mF          │       ΔP applied across membrane
          │   concentration = cF │       ┌─────┐
          │                      │── Am ─┤ mem │── Jw, Js ──► drips into vial
          │   diafiltrate cD     │       └─────┘                   (mV, cV)
          │   feeds in (Lag mode)│
          └──────────────────────┘

    The unknowns the fit estimates
    -------------------------------
    Lp    : water permeability     [L · m⁻² · h⁻¹ · bar⁻¹]
    B     : solute permeability    [μm · s⁻¹]   (membrane mass-transfer coefficient)
    σ     : reflection coefficient [dimensionless, 0–1]
            σ = 1 → perfectly rejecting; σ = 0 → osmotically transparent
    S0/S  : feed flow rate         [g · s⁻¹]    (only Lag / Overflow modes)

    The state variables Pyomo solves for vs (vial n, time τ)
    --------------------------------------------------------
    mF[n,τ]   : cell mass                                [g]
    cF[n,τ]   : cell concentration                       [mM]
    cIn[n,τ] : membrane-interface concentration (incl. concentration polarization) [mM]
    cH[n,τ]   : permeate-side membrane-interface conc.   [mM]
    mV[n,τ]   : cumulative vial mass                     [g]
    cV[n,τ]   : vial-averaged permeate concentration     [mM]
    cVmV[n,τ]: cV × mV   (we integrate the product, then divide back out for cV)
    Jw[n,τ]   : water flux                               [μm · s⁻¹]
    Js[n,τ]   : solute flux                              [μmol · m⁻² · s⁻¹]

    The physical laws encoded as 11 Pyomo Constraint() blocks
    ----------------------------------------------------------
    Each Pyomo Constraint corresponds to one piece of physics:

       ode_mF_rule    →  mass balance on the stirred cell
       ode_cF_rule    →  solute balance in the cell (where Lp, B, σ couple in)
       ode_cH_rule    →  permeate-side back-diffusion
       ode_mV_rule    →  mass accumulating in the vial: dmV/dt = Am·ρ·Jw
       ode_cVmV_rule  →  solute accumulating in the vial: d(cV·mV)/dt = Am·Js
       eqn_S_rule     →  feed flow constraint (Lag / Overflow)
       eqn_cIn_rule   →  concentration polarization: cIn = cF · exp(Jw/k)
       eqn_Jw_rule    →  WATER FLUX:  Jw = Lp · (ΔP − σ · Δπ)
       eqn_Js_rule    →  SOLUTE FLUX from Spiegler–Kedem
       eqn_Js_exp_rule→  the auxiliary exp[Jw·(1−σ)/B] used inside Js
       eqn_H_rule     →  permeate–retentate concentration linkage
       eqn_cV_rule    →  cV · mV identity:  cVmV = cV * mV

    Plus "linking constraints" (mF_linking, cF_linking, ...) that join
    consecutive vials so the cell state is continuous across vial swaps.

    Three operating modes, three boundary cases
    --------------------------------------------
       Lag      : feed inflow = − S₀ − Am·ρ·Jw   (cell loses mass; Lag mode)
       Overflow : feed inflow = + S/3600         (cell stays constant volume)
       DATA     : no inflow (pure dead-end filtration; Am·ρ / M_F0 form)

    These appear as if/elif/elif branches inside ode_mF_rule and ode_cF_rule.

    Arguments
    ---------
    data_stru : dict, experimental data dictionary
    mode      : str, experiment mode, {DATA, Lag, Overflow}
    theta     : dict, preset parameter values
    sim_opt   : boolean, if True run a forward simulation with theta fixed
                         (no parameter estimation, just solve the ODEs)
    B_form    : float/str, different solute permeability formulations
                'single'     - constant B (scalar decision variable)
                'pervial'    - one B per vial
                'convection' - convection–diffusion formulation
                0, -0.5, ... - power-law dependence B(c)
    workflow_family : 'DATA1' / 'DATA2' / 'DATA3' — picks the right
                      defaults and bounds.

    Returns
    -------
    m : pyomo.ConcreteModel — fully-built model, ready for solver or ParmEst.
    '''
    # known parameters
    # applied pressure
    delP = data_stru['data_config']['delP'] #[bar]
    # gas constant
    R = 8.314e-5 #[cm^3 bar / micromol / K]
    # temperature
    T = data_stru['data_config']['Temp'] #[K]
    # membrane area
    Am = data_stru['data_config']['Am'] #[cm^2]
    # density
    rho = data_stru['data_config']['rho'] #[g/cm^3]
    # number of dissolved species
    ni = data_stru['data_config']['ni']
    # kinematic viscosity
    nu = 8.927e-3 #[cm^2/s]
    # diameter of stirred cell
    b = 2.2860 #[cm]
    # Average velocity within the system. DATA3 workbooks can supply the
    # stir speed in the notes; otherwise we keep the historical 350 RPM default.
    rpm = data_stru['data_config'].get('stir_rpm', data_stru['data_config'].get('rpm', 350.0))
    try:
        rpm = float(rpm)
    except Exception:
        rpm = 350.0
    v = rpm / 60 * np.pi * b #[cm/s]
    # diffusion coefficient
    if str(workflow_family).upper() == "DATA3":
        D = _data3_diffusivity_cm2_s(data_stru['data_config'].get('namec'))
    elif isinstance(data_stru['data_config']['namec'], str) and 'K' in data_stru['data_config']['namec']:
        D = 1.960e-5 #[cm^2/s]  - K+
    else:
        raise NotImplementedError
    # mass transfer coefficien
    k = 0.23 * v**0.57 * D**0.67 / (nu**0.24 * b**0.43)

    # known constants
    t_delay = _data1_time_origin(data_stru)
    M_F0 = data_stru['data_config']['M_F0']
    # M_O (mass overflow) is only physically meaningful for Lag/Overflow
    # modes. DATA1 DATA-mode .mat files (closed-cell filtration) sometimes
    # omit it; default to 0 so the load doesn't raise and the eqn_S_rule
    # closure constraint (only enforced for mode != 'DATA') still works.
    M_O = data_stru['data_config'].get('M_O', 0.0)
    # C_D (diafiltrate concentration) only applies when diafiltrate is being
    # fed in. DATA-mode (closed-cell filtration) has no diafiltrate, so the
    # underlying physical default is 0 mM. The DATA-mode ODE multiplies cD
    # by Jw, so cD=0 cleanly zeroes the solute-inflow term.
    cD = data_stru['data_config'].get('C_D', 0.0)
    mH = 0.25 #ml
    C_F0 = data_stru['data_config']['C_F0']
    C_V0 = 1e-6

    N_VIAL = data_stru['data_config']['n']
    # n_v0 = first vial whose residuals enter the objective; defaults to 1
    # (start from vial 1) for DATA1 .mat files that don't set it explicitly.
    N_V0 = data_stru['data_config'].get('n_v0', 1)
    if mode !='DATA':
        N_extra = data_stru['data_config']['n_extra']
        N_H = data_stru['data_config']['n_h']
        N_A = data_stru['data_config']['n_A']       
        if mode == 'Overflow':
            N_A0 = 1
            
    Tauf = 1 #scaled ending time

    # create a model object
    m = ConcreteModel()

    # define the independent variable
    m.n_vial = RangeSet(N_VIAL) # number of vials
    m.tau = ContinuousSet(bounds=(0, Tauf))#scaled time
    
    TF_list = [data_stru['data_raw'][i]['time'][-1]-t_delay for i in range(N_VIAL)]
    TF_dict = dict(zip(m.n_vial,TF_list)) # unscaled time elapse for each vial 
    m.tf = Set(initialize=TF_list)
    
    TI_list = [data_stru['data_raw'][i]['time'][0]-t_delay for i in range(N_VIAL)]
    TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time for each vial
    m.ti = Set(initialize=TI_list)

    # Defensive guard: if a caller explicitly passes workflow_family=None or
    # an unrecognized string, fall back to DATA3. (Note: the function
    # signature default is still 'DATA1' for backward compatibility with
    # legacy DATA1/DATA2 paper-reproduction callers that omit the argument.
    # This guard only fires for None or invalid strings.)
    workflow_family = str(workflow_family or "DATA3").upper()
    if workflow_family not in {"DATA1", "DATA2", "DATA3"}:
        workflow_family = "DATA3"

    # parameter initialization
    param_in = dict()
    mode = _canonical_mode(mode)
    if theta is None:
        param_in['Lp'] = data_stru['data_config']['Lp0']
        param_in['B'] = data_stru['data_config']['B0']
        param_in['sigma'] = data_stru['data_config']['sigma0']
        param_in['beta_0'] = param_in['B']*36000/param_in['Lp']/delP #param_in['B']
        param_in['beta_1'] = 1
        param_in['beta_2'] = 0
        param_in['beta_3'] = 0
        # Seeds for the DATA3 mechanistic partition B-forms (Js = Jw*B*Δc, so the
        # plateau amplitude is seeded from the same B->partition conversion as beta_0).
        _B_amp0 = param_in['B']*36000/param_in['Lp']/delP
        param_in.setdefault('B_inf', _B_amp0)   # saturating plateau  B_inf
        param_in.setdefault('c_star', 20.0)     # saturating knee scale  c*  [mM]
        param_in.setdefault('P0', _B_amp0)       # Donnan bare permeability  P0
        param_in.setdefault('X', 8.0)            # Donnan fixed-charge scale  X  [mM]
        param_in.setdefault('k_dd', 0.5)         # Donnan steric x Born factor  k
        if mode == 'Lag':
            param_in['S0'] = 0      
        elif mode == 'Overflow':
            param_in['S0'] = -M_O/10 *3600
    else:
        param_in = theta
    # define model parameters
    if sim_opt:        
        m.cD = cD
        m.C_H0 = 1e-6
        #simulation parameters
        m.Lp = Param(initialize=param_in['Lp'],mutable=True)
        
        if B_form=='single':
            m.B = Param(initialize=param_in['B'],mutable=True)
        elif B_form=='pervial':    
            m.B = Param(m.n_vial,initialize=param_in['B'],mutable=True)
        elif B_form=='convection':    
            m.beta_0 = Param(initialize=param_in['beta_0'])
            m.H = Var(m.n_vial, m.tau,initialize=0.5)    
        elif isinstance(B_form, str) and B_form in ('sat', 'donnan'):
            # DATA3 NF270 mechanistic partition B-forms (Js = Jw*B*Δc).
            if B_form == 'sat':
                m.B_inf  = Param(initialize=param_in.get('B_inf', 10.0), mutable=True)
                m.c_star = Param(initialize=param_in.get('c_star', 20.0), mutable=True)
            else:
                m.P0   = Param(initialize=param_in.get('P0', 10.0), mutable=True)
                m.X    = Param(initialize=param_in.get('X', 8.0), mutable=True)
                m.k_dd = Param(initialize=param_in.get('k_dd', 0.5), mutable=True)
            # unbounded B[n,t] (matches the polynomial forms): positivity is
            # enforced by the defining constraint, and a hard box makes the DAE
            # infeasible on transient IPOPT steps.
            m.B = Var(m.n_vial, m.tau)
        else:
            m.beta_0 = Param(initialize=param_in['beta_0'],mutable=True)
            if not isinstance(B_form, str):
                if B_form != 0:
                    m.beta_1 = Param(initialize=param_in['beta_1'],mutable=True)
                if B_form > 1:
                    m.beta_2 = Param(initialize=param_in['beta_2'],mutable=True)
                if B_form > 2:
                    m.beta_3 = Param(initialize=param_in['beta_3'],mutable=True) 
                m.B = Var(m.n_vial, m.tau, bounds=(1e-6,50))
        m.sigma = Param(initialize=param_in['sigma'],mutable=True)

    else: 
        #parameter estimation
        m.C_H0 = 1e-6
        m.cD = cD
        m.Lp = Var(bounds=(0.5,50),initialize=param_in['Lp'])

        if B_form=='single':
            # B upper bound bumped 30 → 50 (2026-05-24): the 3 CaCl₂ sheets in
            # the existing centering-fit batch were all pinned EXACTLY at B=30,
            # indicating the bound was clipping. Bumping to 50 gives the
            # optimizer headroom; if the fits still want > 50, that's a
            # separate physics question (the cf-floor change should help too).
            #
            # Per-salt upper override (2026-06-03): for DATA3 NF270 single-salt
            # fits, replace the legacy 50 upper with a valence-graded triple
            # NaCl=30, CaCl2=15, LaCl3=10 (see Architecture.md §17 and the
            # NF270_B_BOUNDS_PER_SALT constant). GUARDED to workflow_family ==
            # "DATA3"; DATA1/DATA2 always see (1e-6, 50), exactly matching the
            # legacy call site above. Salt is read from data_config['namec'],
            # the canonical salt-name field populated by both the DATA3 Excel
            # ingest and DATA1/DATA2 .mat loaders.
            if workflow_family == "DATA3":
                _salt_key = str(data_stru['data_config'].get('namec') or "").strip()
                _b_bounds = NF270_B_BOUNDS_PER_SALT.get(_salt_key, NF270_B_BOUNDS_DEFAULT)
            else:
                _b_bounds = (1e-6, 50)
            m.B = Var(bounds=_b_bounds, initialize=param_in['B'])
        elif B_form=='pervial':
            m.B = Var(m.n_vial, bounds=(1e-6,50),initialize=param_in['B'])
        elif isinstance(B_form, str) and 'convection' in B_form:
            m.beta_0 = Var(bounds=(1+1e-6,50),initialize=param_in['beta_0'])
            m.H = Var(m.n_vial, m.tau,initialize=0.5)  
            m.beta_1 = Var(bounds=(0.0,1.0),initialize=0.5)         
        elif isinstance(B_form, str) and B_form in ('sat', 'donnan'):
            # DATA3 NF270 mechanistic partition B-forms.
            #   saturating:  B = B_inf*(1 - exp(-cc/c_star))
            #   donnan:      B*2cc = P0*(-X + sqrt(X^2 + 4 k^2 cc^2))   (cc = cIn or k_I*cIn)
            if B_form == 'sat':
                m.B_inf  = Var(bounds=(1e-6, 50.0), initialize=param_in.get('B_inf', 10.0))
                m.c_star = Var(bounds=(1e-3, 500.0), initialize=param_in.get('c_star', 20.0))
            else:
                m.P0   = Var(bounds=(1e-6, 100.0), initialize=param_in.get('P0', 10.0))
                if workflow_family == "DATA3" and NF270_DONNAN_FIX_X is not None:
                    m.X = Param(initialize=float(NF270_DONNAN_FIX_X), mutable=True)  # X-profile probe
                else:
                    m.X    = Var(bounds=(1e-3, 200.0), initialize=param_in.get('X', 8.0))
                if workflow_family == "DATA3" and NF270_DONNAN_FIX_K is not None:
                    m.k_dd = Param(initialize=float(NF270_DONNAN_FIX_K))  # fixed-k fallback
                else:
                    m.k_dd = Var(bounds=(1e-3, 5.0), initialize=param_in.get('k_dd', 0.5))
            # unbounded B[n,t] (matches the polynomial forms): positivity comes
            # from the defining constraint; a hard box makes the DAE infeasible
            # on transient IPOPT steps.
            m.B = Var(m.n_vial, m.tau)
        else:
            m.beta_0 = Var(initialize=param_in['beta_0'])
            m.B = Var(m.n_vial, m.tau)
            # Coefficient bounds.  DATA1/DATA2 keep the legacy +/-20 exactly.  For
            # DATA3 the term beta_k * c_in^k is badly scaled at high feed
            # concentration (c_in ~ 95 mM lets beta_2*c^2 reach ~2e5 under the
            # legacy bound, so the quadratic/cubic DAE fit is numerically
            # infeasible).  Scale each bound by the in-run feed range so the term
            # stays physically bounded (|beta_k| * c_max^k <= 30).  Gated on
            # workflow_family => DATA1/DATA2 reproductions are byte-identical.
            def _beta_bnd(_k):
                if workflow_family != "DATA3":
                    return 20.0
                try:
                    _cm = np.nanmax([np.nanmax(np.asarray(r.get('cF_exp', [np.nan]), float))
                                     for r in data_stru['data_raw']])
                    _cm = float(_cm) if (np.isfinite(_cm) and _cm > 1.0) else 100.0
                except Exception:
                    _cm = 100.0
                return 30.0 / (_cm ** _k)
            if B_form != 0:
                _bb1 = _beta_bnd(1)
                m.beta_1 = Var(bounds=(-_bb1, _bb1), initialize=param_in['beta_1'])
            if B_form > 1:
                _bb2 = _beta_bnd(2)
                m.beta_2 = Var(bounds=(-_bb2, _bb2), initialize=param_in['beta_2'])
            if B_form > 2:
                _bb3 = _beta_bnd(3)
                m.beta_3 = Var(bounds=(-_bb3, _bb3), initialize=param_in['beta_3'])
                
        # Honor module-level interior-bound override for NF270 corner-pin
        # diagnostics.  GUARDED to DATA3 only — DATA1/DATA2 always use (0, 1)
        # regardless of the constant's value, so byte-identical legacy behavior
        # is preserved even if a sibling script flips the toggle globally.
        if workflow_family == "DATA3" and NF270_SIGMA_INTERIOR_BOUNDS is not None:
            _sig_bounds = NF270_SIGMA_INTERIOR_BOUNDS
            _sig_init   = max(_sig_bounds[0] + 1e-3, min(_sig_bounds[1] - 1e-3, param_in['sigma']))
        else:
            _sig_bounds = (0.0, 1.0)
            _sig_init   = param_in['sigma']
        m.sigma = Var(bounds=_sig_bounds, initialize=_sig_init)
    
    if mode == 'Lag':
        if sim_opt:
            m.S0 = Param(initialize=param_in['S0'],mutable=True) 
            m.S = Param(initialize=param_in['S'],mutable=True)
        else:
            # g/s
            m.S0 = Param(initialize=param_in['S0'])
            # g/hr
            m.S = Var(initialize=M_O/sum(m.tf-m.ti)*3600)
    elif mode == 'Overflow':
        if sim_opt:
            m.S0 = Param(initialize=param_in['S0'],mutable=True)
            m.S = Param(initialize=param_in['S'],mutable=True)
        else:
            # g/s
            m.S = Var(initialize=0)
            # g/hr
            m.S0 = Var(initialize=-M_O/10)
            
    # define dependent variables 
    if mode !='DATA':
        m.mF = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=M_F0)
    m.cF = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=C_F0)
    m.cIn = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=C_F0)
    m.cH = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=1e-6)
    m.mV = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=1e-6)
    m.cVmV = Var(m.n_vial, m.tau,initialize=C_V0 * 1e-6)
    m.cV = Var(m.n_vial,  m.tau, domain=NonNegativeReals,initialize=1e-6)

    # define intermediate variables
    m.Jw = Var(m.n_vial, m.tau)
    m.Js = Var(m.n_vial, m.tau)
    if isinstance(B_form, str) and 'convection' in B_form:
        m.Js_exp = Var(m.n_vial, m.tau, bounds=(1+1e-6,1e4))
    elif B_form=='K':   
        m.Js_exp = Var(m.n_vial, m.tau)

    # define derivatives
    if mode !='DATA':
        m.dmF = DerivativeVar(m.mF,wrt=m.tau)
    m.dcF = DerivativeVar(m.cF,wrt=m.tau)
    m.dcH = DerivativeVar(m.cH,wrt=m.tau)
    m.dmV = DerivativeVar(m.mV,wrt=m.tau)
    m.dcVmV = DerivativeVar(m.cVmV,wrt=m.tau)
    
    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 1 of 11:  CELL MASS BALANCE                             ║
    # ║                                                                    ║
    # ║   dmF/dt = (feed in)  −  Am · ρ · Jw                              ║
    # ║                          └─────┬─────┘                            ║
    # ║                          mass leaving through the membrane         ║
    # ║                          per unit time:                            ║
    # ║                            Am [m²] × ρ [g/cm³] × Jw [cm/s]         ║
    # ║                                                                    ║
    # ║ "feed in" depends on what mode we're in:                          ║
    # ║   Lag      : feed inflow = − S₀  (S₀ < 0 means cell loses mass)   ║
    # ║   Overflow : feed inflow = − S₀  while concentrating               ║
    # ║              feed inflow = + S/3600  while diafiltrating           ║
    # ║   DATA     : no inflow; cell is closed (no constraint added)      ║
    # ║                                                                    ║
    # ║ The trailing (TF_dict[n] − TI_dict[n]) / Tauf factor converts      ║
    # ║ the equation onto the "normalized time" axis Pyomo.DAE uses.       ║
    # ║ Each vial has its own physical time window [TI, TF]; Tauf is the   ║
    # ║ normalized end (always 1.0). The factor rescales d/dt accordingly. ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def ode_mF_rule(m, n, t):
        if mode == 'Overflow':
            if n < N_A0:
                # Initial concentration leg (before diafiltrate starts): pure
                # filtration, no feed in. dmF/dt = − Am·ρ·Jw
                return m.dmF[n,t] == (0 - Am * rho * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
            elif N_A0-1 < n <= N_A:
                # Active diafiltration leg (Overflow specific): outflow S₀
                # dominates. dmF/dt = − S₀ − Am·ρ·Jw
                return m.dmF[n,t] == (- m.S0  - Am * rho * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
            else:
                # Final concentration leg: inflow S exactly balances outflow,
                # so dmF/dt = S/3600 (units: g/hr → g/s).
                return m.dmF[n,t] == m.S / 3600  * (TF_dict[n]-TI_dict[n])/Tauf
        elif mode == 'Lag':
            if n <= N_A:
                # Active Lag leg: dmF/dt = − S₀ − Am·ρ·Jw
                return m.dmF[n,t] == (- m.S0 - Am * rho * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
            else:
                # Final leg (closed): dmF/dt = S/3600
                return m.dmF[n,t] == m.S / 3600  * (TF_dict[n]-TI_dict[n])/Tauf
    if mode !='DATA':
        m.ode_mF = Constraint(m.n_vial, m.tau, rule=ode_mF_rule)
     
    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 2 of 11:  CELL SOLUTE (CONCENTRATION) BALANCE           ║
    # ║                                                                    ║
    # ║ Apply mass balance to the SOLUTE in the stirred cell:              ║
    # ║                                                                    ║
    # ║   d(mF · cF)/dt  =  (solute in from diafiltrate)                  ║
    # ║                   − (solute out through membrane)                  ║
    # ║                                                                    ║
    # ║ Using the product rule and substituting dmF/dt from Eq. 1:         ║
    # ║                                                                    ║
    # ║   dcF/dt  =  1/mF · [ (diafiltrate−cell) · feed_rate              ║
    # ║                       + Am·ρ · (cell_or_diafiltrate · Jw − Js) ]  ║
    # ║                                                                    ║
    # ║ where:                                                             ║
    # ║   cF  · Jw  = solute carried OUT with the water (convection)       ║
    # ║   Js        = solute moving by diffusion (Spiegler–Kedem, Eq. 8)   ║
    # ║   cD · Jw  = solute carried IN with the diafiltrate (Overflow)    ║
    # ║                                                                    ║
    # ║ DATA mode is simpler: closed cell, fixed total mass M_F0, so:      ║
    # ║   dcF/dt = Am·ρ / M_F0 · (cD·Jw − Js)                              ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def ode_cF_rule(m, n, t):
        if mode == 'Overflow':
            if n < N_A0:
                # Initial leg, no diafiltrate yet: solute only leaves the cell.
                return m.dcF[n,t] == 1 / m.mF[n,t] * (Am * rho * (m.cF[n,t] * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf
            elif N_A0-1 < n <= N_A:
                # Active leg: (cF − cD) outflow term + membrane transport.
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cF[n,t] - m.cD ) * m.S0 + Am * rho * (m.cF[n,t] * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf
            else:
                # Final leg: (cD − cF) inflow drives the cell toward cD.
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cD - m.cF[n,t]) * m.S / 3600 + Am * rho * (m.cD * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf
        elif mode == 'Lag':
            if n <= N_A:
                # Active Lag leg: same form as Overflow active, different S₀ sign.
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cF[n,t] - m.cD ) * m.S0 + Am * rho * (m.cF[n,t] * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf
            else:
                # Final Lag leg: cell relaxes toward cD via inflow.
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cD - m.cF[n,t]) * m.S / 3600 + Am * rho * (m.cD * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf
        elif mode == 'DATA':
            # DATA mode = closed cell, no inflow; one constant M_F0 in the denominator.
            return m.dcF[n,t] == Am * rho / M_F0 * (m.cD * m.Jw[n,t] - m.Js[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf

    m.ode_cF = Constraint(m.n_vial, m.tau, rule=ode_cF_rule)

    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 3 of 11:  PERMEATE-SIDE (cH) BACK-DIFFUSION             ║
    # ║                                                                    ║
    # ║ cH is the solute concentration on the PERMEATE side of the         ║
    # ║ membrane (between membrane and vial).                              ║
    # ║                                                                    ║
    # ║   dcH/dt = Am·ρ / mV · ( Js  −  cH · Jw )                          ║
    # ║                          ^^^^^   ^^^^^^^^^                          ║
    # ║                          solute    solute carried                   ║
    # ║                          arriving  away with the                    ║
    # ║                          via the   water as it                      ║
    # ║                          membrane  drains into the vial             ║
    # ║                                                                    ║
    # ║ Net rate of change of cH is the inflow from the membrane minus      ║
    # ║ the outflow via water carrying solute out. Same form as the cell    ║
    # ║ balance but for the small permeate-side volume.                     ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def ode_cH_rule(m, n, t):
        if mode !='DATA' and n <= N_H:
            # Holdup-vial branch: use the time-varying mV in the denominator.
            return m.dcH[n,t] == Am * rho / m.mV[n,t] * (m.Js[n,t] - m.cH[n,t] * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
        else:
            # Steady-volume branch: use the fixed holdup mass mH.
            return m.dcH[n,t] == Am * rho / mH * (m.Js[n,t] - m.cH[n,t] * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_cH = Constraint(m.n_vial, m.tau, rule=ode_cH_rule)

    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 4 of 11:  VIAL MASS ACCUMULATION                         ║
    # ║                                                                    ║
    # ║   dmV/dt  =  Am · ρ · Jw                                           ║
    # ║                                                                    ║
    # ║ Mass accumulates in the collection vial at exactly the rate water  ║
    # ║ is fluxing through the membrane. This is what the lab balance      ║
    # ║ reads in real time.                                                ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def ode_mV_rule(m, n, t):
        return m.dmV[n,t] == m.Jw[n,t] * Am * rho * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_mV = Constraint(m.n_vial, m.tau, rule=ode_mV_rule)

    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 5 of 11:  VIAL SOLUTE ACCUMULATION                       ║
    # ║                                                                    ║
    # ║   d(mV · cV)/dt  =  Am · ρ · cH · Jw                               ║
    # ║                                                                    ║
    # ║ Solute mass in the vial grows at the rate (concentration of fluid  ║
    # ║ leaving the membrane) × (mass flow of water). We track the PRODUCT ║
    # ║ mV·cV instead of cV directly because it integrates cleanly even    ║
    # ║ when mV starts near zero (where cV is numerically ill-defined).    ║
    # ║ Then Equation 11 (eqn_cV_rule) recovers cV = mV·cV / mV.           ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def ode_cVmV_rule(m, n, t):
        if mode !='DATA' and n <= N_H:
            # Holdup-vial form: incoming solute is Js directly (matches ode_cH form).
            return m.dcVmV[n,t] == Am * rho * m.Js[n,t] * (TF_dict[n]-TI_dict[n])/Tauf
        else:
            # Standard form: water (Jw) carries solute (cH) into the vial.
            return m.dcVmV[n,t] == m.Jw[n,t] * m.cH[n,t] * Am * rho * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_cVmV = Constraint(m.n_vial, m.tau, rule=ode_cVmV_rule)

    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 6 of 11:  GLOBAL MASS CLOSURE (Lag / Overflow only)     ║
    # ║                                                                    ║
    # ║   mF(end of last vial)  −  M_F0  =  M_O                            ║
    # ║                                                                    ║
    # ║ Total cell mass change across the WHOLE experiment must equal      ║
    # ║ M_O (the measured signed mass loss, = Final − Initial weight).     ║
    # ║ This pins the otherwise-free S (final-leg inflow rate) so the      ║
    # ║ bookkeeping balances at the end. Skipped when sim_opt=True because ║
    # ║ then S is a fixed parameter, not a decision variable.              ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def eqn_S_rule(m):
        return m.mF[m.n_vial.last(),Tauf] - M_F0 == M_O
    if mode !='DATA' and not sim_opt:
        m.eqn_S = Constraint(rule=eqn_S_rule)

    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 7 of 11:  CONCENTRATION POLARIZATION (cIn vs cF)         ║
    # ║                                                                    ║
    # ║   cIn  =  (cF − cH) · exp(Jw / k)  +  cH                           ║
    # ║                                                                    ║
    # ║ Right next to the membrane there is a stagnant film of solution    ║
    # ║ that the stirrer can't reach. As water rushes toward the membrane  ║
    # ║ at rate Jw, it carries solute into this film faster than the       ║
    # ║ stirrer can carry it back out. So the concentration RIGHT AT the   ║
    # ║ membrane surface (cIn) is HIGHER than the bulk cell concentration  ║
    # ║ (cF). This is "concentration polarization."                        ║
    # ║                                                                    ║
    # ║ The classic film-model expression for it is:                       ║
    # ║   (cIn − cH) / (cF − cH)  =  exp(Jw / k)                           ║
    # ║ where k is the back-diffusion mass-transfer coefficient (set       ║
    # ║ earlier from membrane and Reynolds-number correlations).           ║
    # ║                                                                    ║
    # ║ Rearranged: cIn = (cF − cH) · exp(Jw / k) + cH.                    ║
    # ║                                                                    ║
    # ║ Note: ALL the other equations use cIn (the membrane-surface conc.) ║
    # ║ as the driving force for flux, not the bulk cF.                    ║
    # ║                                                                    ║
    # ║ Reference: Mulder, "Basic Principles of Membrane Technology"       ║
    # ║ Ch. 7 (concentration polarization in pressure-driven processes).   ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def eqn_cIn_rule(m, n, t):
        return m.cIn[n,t] == (m.cF[n,t] - m.cH[n,t]) * exp(m.Jw[n,t]/ k) + m.cH[n,t]
    m.eqn_cIn = Constraint(m.n_vial, m.tau, rule=eqn_cIn_rule)

    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 8 of 11:  WATER FLUX (the heart of the Lp/σ fit)        ║
    # ║                                                                    ║
    # ║   Jw  =  Lp · ( ΔP  −  σ · Δπ )                                    ║
    # ║                  ^^      ^^^^^^                                     ║
    # ║                  applied  osmotic                                   ║
    # ║                  pressure pressure                                  ║
    # ║                  pushing  pushing                                   ║
    # ║                  water    water                                     ║
    # ║                  THROUGH  BACK                                      ║
    # ║                                                                    ║
    # ║ where the osmotic pressure difference is given by van 't Hoff:    ║
    # ║                                                                    ║
    # ║   Δπ  =  (cIn − cH) · ni · R · T                                   ║
    # ║                                                                    ║
    # ║   Lp  : water permeability                                         ║
    # ║   σ   : reflection coefficient (0 = transparent; 1 = perfect       ║
    # ║         rejection, full osmotic pressure pushes back)              ║
    # ║   ni  : van 't Hoff factor (number of dissolved species per salt   ║
    # ║         formula unit, e.g. 2 for NaCl, 3 for CaCl₂)                ║
    # ║   R   : gas constant in [cm³·bar / μmol / K]                       ║
    # ║   T   : temperature [K]                                            ║
    # ║                                                                    ║
    # ║ Implementation note: the LEFT side carries a factor of 36000       ║
    # ║ because Jw is stored in [cm/s × 10⁻⁴] (i.e. μm/s) while Lp uses    ║
    # ║ [L/m²/h/bar] — the 36000 reconciles the two unit conventions.      ║
    # ║                                                                    ║
    # ║ Reference: Kedem & Katchalsky (1958), "Thermodynamic analysis of   ║
    # ║ the permeability of biological membranes to non-electrolytes."     ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def eqn_Jw_rule(m, n, t):
        return m.Jw[n,t]*36000 == m.Lp *(delP - (m.cIn[n,t] - m.cH[n,t])*ni*m.sigma*R*T)
    m.eqn_Jw = Constraint(m.n_vial, m.tau, rule=eqn_Jw_rule)

    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 9 of 11:  SOLUTE FLUX                                    ║
    # ║                                                                    ║
    # ║ The simple form (B_form='single' or 'pervial') is Fickian:         ║
    # ║                                                                    ║
    # ║   Js  =  B · (cIn − cH)                                            ║
    # ║                                                                    ║
    # ║   B   : solute permeability [μm/s]                                 ║
    # ║   The driving force is the SURFACE concentration difference, NOT   ║
    # ║   the bulk-vs-permeate difference (concentration polarization is   ║
    # ║   already baked in via cIn — see Eq. 7).                           ║
    # ║                                                                    ║
    # ║ The convection form (B_form='convection') uses Spiegler–Kedem:    ║
    # ║                                                                    ║
    # ║   Js  =  Jw · H · (cIn · F − cH) / (F − 1)                         ║
    # ║                                                                    ║
    # ║ with F = exp(Jw·(1−σ)/Pm) — this is the "F factor" everyone in     ║
    # ║ the RO/NF literature refers to. It captures the coupling between   ║
    # ║ water transport and solute back-diffusion through the membrane.    ║
    # ║                                                                    ║
    # ║ The 10000 multiplier on the left rescales Js's storage units to    ║
    # ║ keep the residuals well-conditioned numerically (no physics change). ║
    # ║                                                                    ║
    # ║ Reference: Spiegler & Kedem (1966), "Thermodynamics of hyper-      ║
    # ║ filtration (reverse osmosis): criteria for efficient membranes."   ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def eqn_Js_rule(m, n, t):
        if B_form=='single':
            # Most common case: single scalar B. Fickian diffusion.
            return m.Js[n,t]*10000 == m.B * (m.cIn[n,t] - m.cH[n,t])
        elif B_form=='pervial':
            # Same form but B is a separate decision variable per vial.
            return m.Js[n,t]*10000 == m.B[n] * (m.cIn[n,t] - m.cH[n,t])
        elif B_form=='convection':
            # Spiegler–Kedem convection–diffusion form (eq. above).
            return m.Js[n,t] == m.Jw[n,t] * m.H[n,t] * (m.cIn[n,t] * m.Js_exp[n,t] - m.cH[n,t]) / (m.Js_exp[n,t]-1)
        else:
            #return m.Js[n,t]*10000 ==  m.B[n,t] * (m.cIn[n,t] - m.cH[n,t]) # No Jw in Js
            return m.Js[n,t]*10000 ==  (m.Jw[n,t] * m.B[n,t] * 10000) * (m.cIn[n,t] - m.cH[n,t]) #Jw in Js

    m.eqn_Js = Constraint(m.n_vial, m.tau, rule=eqn_Js_rule)

    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 10a of 11:  SPIEGLER–KEDEM F-FACTOR (auxiliary)         ║
    # ║                                                                    ║
    # ║   F  =  exp( Jw · (1 − σ) / Pm )      (Pm ≡ beta_0 in the code)   ║
    # ║                                                                    ║
    # ║ This is the dimensionless exponential that appears throughout the  ║
    # ║ Spiegler–Kedem formulation. It captures the coupling: at high Jw   ║
    # ║ or low membrane permeability Pm, F → large, and rejection → σ.    ║
    # ║ Only active when B_form='convection'; otherwise skipped.            ║
    # ║                                                                    ║
    # ║ The extra factor of 10000 inside the exp() rescales Jw and beta_0  ║
    # ║ into a consistent unit system for the exponent's argument.         ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def eqn_Js_exp_rule(m, n, t):
        if B_form=='convection':
            return m.Js_exp[n,t]  == exp(m.Jw[n,t]/m.beta_0*10000)
        else:
            return Constraint.Skip
    m.eqn_Js_exp = Constraint(m.n_vial, m.tau, rule=eqn_Js_exp_rule)

    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 10b of 11:  H factor (convection mode only)              ║
    # ║                                                                    ║
    # ║   H  =  beta_1                                                     ║
    # ║                                                                    ║
    # ║ H is a piece-wise constant (per vial × time) that lets H vary if   ║
    # ║ you turn it into a Var later. Currently pinned to the scalar beta_1║
    # ║ for the convection-mode fit.                                       ║
    # ╚══════════════════════════════════════════════════════════════════╝
    def eqn_H_rule(m,n,t):
        if B_form=='convection':
            return m.H[n,t]  == m.beta_1
        else:
            return Constraint.Skip
    m.eqn_H = Constraint(m.n_vial, m.tau, rule=eqn_H_rule)

    # ╔══════════════════════════════════════════════════════════════════╗
    # ║ EQUATION 11 of 11:  VIAL-CONCENTRATION IDENTITY                   ║
    # ║                                                                    ║
    # ║   cVmV  =  mV · cV                                                 ║
    # ║                                                                    ║
    # ║ Recall Eq. 5 integrated d(mV·cV)/dt directly to get cVmV. This     ║
    # ║ algebraic constraint recovers cV from cVmV by dividing by mV.      ║
    # ║                                                                    ║
    # ║ We write it as MULTIPLICATION  (cVmV = mV·cV)  instead of          ║
    # ║ DIVISION  (cV = cVmV/mV)  because IPOPT handles the multiplicative ║
    # ║ form much better numerically — division by a small mV near vial    ║
    # ║ open would create huge Jacobian entries that destabilize the solve.║
    # ╚══════════════════════════════════════════════════════════════════╝
    def eqn_cV_rule(m, n, t):
        return m.mV[n,t] * m.cV[n,t] == m.cVmV[n,t]
    m.eqn_cV = Constraint(m.n_vial, m.tau, rule=eqn_cV_rule)# numercally better without division
    
    # link the variables from different vials 
    def mF_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        else:
            return m.mF[n,Tauf]==m.mF[n+1,0]  
    if mode !='DATA':
        m.mF_linking = Constraint(m.n_vial,rule=mF_linking_rule) 
    
    def cF_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        else:
            return m.cF[n,Tauf]==m.cF[n+1,0]                                                                  
    m.cF_linking = Constraint(m.n_vial,rule=cF_linking_rule)   

    def cH_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        else:
            return m.cH[n,Tauf]==m.cH[n+1,0]                                                                  
    m.cH_linking = Constraint(m.n_vial,rule=cH_linking_rule)      

    def mV_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        elif mode !='DATA'and N_H < N_V0 and n < N_H:
            return m.mV[n,Tauf]==m.mV[n+1,0]
        elif mode !='DATA'and N_H < N_V0 and N_H < n < N_V0-1:
            return m.mV[n,Tauf]==m.mV[n+1,0]
        elif mode !='DATA'and N_V0-1 < N_extra and N_V0 <= n < N_extra+1:
            return m.mV[n,Tauf]==m.mV[n+1,0]
        else:
            return 1e-6==m.mV[n+1,0]                                                                  
    m.mV_linking = Constraint(m.n_vial,rule=mV_linking_rule) 
    
    def cVmV_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        elif mode !='DATA'and N_H < N_V0 and n < N_H:
            return m.cVmV[n,Tauf]==m.cVmV[n+1,0]
        elif mode !='DATA'and N_H < N_V0 and N_H < n < N_V0-1:
            return m.cVmV[n,Tauf]==m.cVmV[n+1,0]
        elif mode !='DATA'and N_V0-1 < N_extra and N_V0 <= n < N_extra+1:
            return m.cVmV[n,Tauf]==m.cVmV[n+1,0]
        else:
            return m.cV[n,Tauf]*1e-6==m.cVmV[n+1,0]                                                                  
    m.cVmV_linking = Constraint(m.n_vial,rule=cVmV_linking_rule)

    # Optional ionic-strength reparameterization of the concentration-dependent
    # B polynomial (DATA3-only, gated on NF270_B_USE_IONIC_STRENGTH).  When
    # enabled, the concentration argument that drives B is scaled by k_I so that
    # B is a function of ionic strength I = k_I*c rather than bare concentration
    # c.  When the toggle is off (or workflow_family != DATA3), _b_kI == 1.0 and
    # the local `cc` binding below IS the same Pyomo component as the bare
    # concentration, so each constraint expression is byte-identical to the
    # legacy path — DATA1/DATA2 reproduction is unaffected regardless of value.
    _b_use_I = bool(workflow_family == "DATA3" and NF270_B_USE_IONIC_STRENGTH)
    if _b_use_I:
        _b_salt = str(data_stru['data_config'].get('namec') or "").strip()
        _b_ni = data_stru['data_config'].get('ni')
        _b_zc = (float(_b_ni) - 1.0) if (_b_ni and float(_b_ni) > 1.0) else None
        _b_kI = nf270_ionic_strength_factor(_b_salt, cation_valence=_b_zc)
    else:
        _b_kI = 1.0

    # B continuity across vials in the per-vial form.
    def B_form_rule1(m,n):
        if n < N_V0-1:
            return m.B[n] == m.B[n+1]
        else:
            return Constraint.Skip

    # B(cF) polynomial form used in the paper variants.
    # Example: B = beta_0 + beta_1*cF^p + beta_2*cF^(p+1) + ...
    def B_form_rule2(m,n,t):
        if not isinstance(B_form, str):
            # cc == m.cF[n,t] when the ionic-strength toggle is off (byte-identical);
            # k_I*m.cF[n,t] when on (B driven by bulk ionic strength).
            cc = (_b_kI * m.cF[n,t]) if _b_use_I else m.cF[n,t]
            if B_form > 2:
                return m.B[n,t] == m.beta_0 + m.beta_1 * cc**(B_form - 2) + m.beta_2 * cc**(B_form - 1) + m.beta_3 * cc**B_form #+ m.beta_4 * m.cH[n,t]
            elif B_form > 1:
                return m.B[n,t] == m.beta_0 + m.beta_1 * cc**(B_form - 1) + m.beta_2 * cc**B_form #+ m.beta_4 * m.cH[n,t]
            else:
                return m.B[n,t] == m.beta_0 + m.beta_1 * cc**B_form
        else:
            return Constraint.Skip

    # B(cIn) polynomial form used by the DATA1/DATA2 model variants.
    # The exponent p is controlled by B_form.
    def B_form_rule3(m,n,t):
        # DATA3 mechanistic partition forms (string B_form): defined on the same
        # interfacial concentration cc (= cIn, or k_I*cIn under the ionic-strength
        # toggle) and the same Jw-scaled Js as the polynomial forms.
        if isinstance(B_form, str) and B_form in ('sat', 'donnan'):
            cc = (_b_kI * m.cIn[n,t]) if _b_use_I else m.cIn[n,t]
            if B_form == 'sat':
                return m.B[n,t] == m.B_inf * (1 - exp(-cc / m.c_star))
            return m.B[n,t] * (2.0*cc) == m.P0 * (-m.X + sqrt(m.X**2 + 4.0*(m.k_dd**2)*(cc**2)))
        if not isinstance(B_form, str):
            # cc == m.cIn[n,t] when the ionic-strength toggle is off (byte-identical);
            # k_I*m.cIn[n,t] when on, so B is driven by the INTERFACIAL ionic
            # strength I_in = k_I*cIn (the membrane-wall concentration the solute
            # flux Js = B*(cIn - cH) actually sees).  See Architecture.md §18.3.
            cc = (_b_kI * m.cIn[n,t]) if _b_use_I else m.cIn[n,t]
            if B_form > 2:
                return m.B[n,t] == m.beta_0 + m.beta_1 * cc**(B_form - 2) + m.beta_2 * cc**(B_form - 1) + m.beta_3 * cc**B_form# + m.beta_4 * m.cH[n,t]
            elif B_form > 1:
                return m.B[n,t] == m.beta_0 + m.beta_1 * cc**(B_form - 1) + m.beta_2 * cc**B_form# + m.beta_4 * m.cH[n,t]
            elif B_form == 0:
                return m.B[n,t] == m.beta_0
            elif B_form < 0:
                return m.B[n,t] * cc**(-B_form) == m.beta_0 * cc**(-B_form) + m.beta_1
            else:
                return m.B[n,t] == m.beta_0 + m.beta_1 * cc**B_form
        else:
            return Constraint.Skip
    
    if B_form=='pervial'and not sim_opt: 
        m.b_form = Constraint(m.n_vial,rule=B_form_rule1)
        pass
    else:
        m.b_form = Constraint(m.n_vial,m.tau,rule=B_form_rule3) 

    # Initial conditions for the ODE system.
    def _init(m): 
        def firstNonNan(listfloats):
            for item in listfloats:
                if np.isnan(item) == False:
                    return item
        if mode !='DATA':
            yield m.mF[1,0]==data_stru['data_config']['M_F0']
            yield m.cH[1,0]==m.C_H0 #1e-6
            yield m.cF[1,0]==firstNonNan(data_stru['data_raw'][0]['cF_exp'])
        else:
            # DATA mode initializes the permeate concentration from the first available measurement.
            # If the saved file has only NaN values, fall back to a tiny positive default.
            if np.ndim(data_stru['data_raw'][0]['cV_avg']) > 0:
                C_H0 = firstNonNan(data_stru['data_raw'][0]['cV_avg'])
            else:
                C_H0 = data_stru['data_raw'][0]['cV_avg']
            if C_H0 is None or (isinstance(C_H0, float) and np.isnan(C_H0)):
                C_H0 = 1e-6
            yield m.cH[1,0]==C_H0 * 0.8
            yield m.cF[1,0]==firstNonNan(data_stru['data_raw'][0]['cF_exp'])
               
        yield m.cVmV[1,0]==1e-6*1e-6
        yield m.mV[1,0]==1e-6
        yield ConstraintList.End
        
    m.con_boundary = ConstraintList(rule=_init)
        
    return m


def solve_model(
    data_stru,
    mode,
    theta=None,
    sim_opt=False,
    B_form='single',
    LOUD=False,
    workflow_family='DATA1',
    nfe=300,
    solver_max_iter=3000,
    solver_retry_max_iter=5000,
    solver_max_cpu_time=None,
    mass_litmus_test=False,
    multistart=False,
    multistart_iterations=10,
    multistart_seed=13,
    band_vials=None,
    skip_sim_init=False,
):
    """
    Solve pyomo model

    Arguments:
        band_vials: optional iterable of 1-indexed vial numbers. When given,
            the WSSE objective is MASKED to only those vials (out-of-band vials
            are skipped in the per-vial residual loop), so one theta is fit to a
            single concentration band. None (default) fits all vials and is
            byte-identical to the legacy objective. Used by the per-concentration
            -band diagnostic (Architecture.md §18.4); single-shot only (not
            wired through multistart).
        data_stru: dict, experimental data dictionary
        mode: str, experiment mode, {DATA, lag, overflow}
        theta: dict, preset parameter values
        sim_opt: boolean, if run simulation with fixed parameter values
        B_form: float/str, different solute permeability coefficient formula
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration
        LOUD: boolean, if print out parameter results and store model predictions
        mass_litmus_test: boolean, if True pin sigma=0 and use a mass-only
            WSSE objective. Collapses the fit to its simplest form — a
            one-parameter Lp fit against vial mass slopes, no concentration
            coupling. If even THIS doesn't match, the issue is upstream of
            the optimizer (model, data loading, or units).

    Returns:
        fit_stru: dict, parameter fit results
        sim_stru: dict, model predictions
        sim_inter: dict, model predictions for experimental measurements
    """

    # Honor the module-level diagnostic toggle so the runfile can flip it.
    mass_litmus_test = bool(mass_litmus_test) or bool(NF270_MASS_LITMUS_TEST_ACTIVE)

    # When the litmus test is active, redirect to the closed-form
    # mass-balance-only solver. That bypasses Pyomo / IPOPT entirely and
    # eliminates the implementation slop where B was a free Var with no
    # objective constraint. The fit_stru / sim_stru shape is identical.
    if mass_litmus_test and not sim_opt:
        print("###################################################################")
        print("[mass-litmus test] redirecting to solve_mass_balance_only")
        print("                   (sigma=0, B dropped, closed-form np regression)")
        fit_stru, sim_stru = solve_mass_balance_only(
            data_stru, B_form=B_form, LOUD=LOUD,
        )
        # solve_model historically returned (fit_stru, sim_stru, sim_inter)
        # but many call sites pack it as (fit_stru, sim_stru) and discard
        # the third. Synthesize a sim_inter with the same shape as sim_stru
        # so neither pattern breaks.
        sim_inter = [{"time": s["time"], "mV": s["mV"]} for s in sim_stru]
        return fit_stru, sim_stru, sim_inter

    # Per-concentration-band masking is single-shot (the obj_rule closure below
    # honors band_vials); it is not threaded through the multistart helper.
    if band_vials is not None:
        band_vials = set(int(v) for v in band_vials)
        if bool(multistart):
            print("[band] band_vials set -> disabling multistart (single-shot band fit).")
            multistart = False

    if bool(multistart) and not sim_opt:
        return _solve_model_multistart(
            data_stru,
            mode,
            theta=theta,
            B_form=B_form,
            LOUD=LOUD,
            workflow_family=workflow_family,
            nfe=nfe,
            solver_max_iter=solver_max_iter,
            solver_retry_max_iter=solver_retry_max_iter,
            solver_max_cpu_time=solver_max_cpu_time,
            mass_litmus_test=mass_litmus_test,
            multistart_iterations=multistart_iterations,
            multistart_seed=multistart_seed,
        )

    print("###################################################################")
    print("Creating and solving the Pyomo model with the following settings: ")
    print("mode =", mode)
    print("theta =",theta)
    print("sim_opt =", sim_opt)
    print("B_form =", B_form)
    print(" ")

    data_stru = _normalize_conductivity_measurements(data_stru)

    # interpolation function
    def interpolation(m,var,n_vial,t):
        """
        Interpolation function for pyomo model
        """
        tp = list(m.tau)
        for j in range(0,len(tp)-1):
            if tp[j+1]>=t and tp[j]<=t:
                var_inter = (var[n_vial,tp[j+1]]-var[n_vial,tp[j]]) * (t - tp[j])/(tp[j+1]-tp[j]) + var[n_vial,tp[j]]
        return var_inter
    
    # Objective function:
    # minimize a weighted sum of squared residuals for mass, permeate concentration,
    # and retentate concentration.
    def obj_rule(m):
        """
        pyomo Objective function (residual calculated using interpolation from model)
        """
        # Keep the DATA2-style 0.01 g mass weight unless a caller has
        # explicitly injected a different value into data_config.
        mass_scale = float(data_stru["data_config"].get("mass_scale_g", 0.01))
        if not np.isfinite(mass_scale) or mass_scale <= 0:
            mass_scale = 0.01
        # DATA3-only paper-justified retentate/permeate concentration sigmas
        # (legacy 0.003 / 0.03 for DATA1/DATA2 and unset knobs).
        _cf_frac, _cp_frac = _nf270_conc_scales(workflow_family)
        obj_m = 0
        obj_cp = 0
        obj_cf0 = 0
        obj_cf = 0
        ob_m = 0
        ob_cp = 0
        ob_cf0 = 0
        ob_cf = 0
        res_m_assemble = []
        res_cp_assemble = []
        res_cf_assemble = []
        # Optional permeate-probe residual accumulator (DATA3-only).  Stays
        # 0 unless NF270_USE_PERMEATE_PROBE is on AND cV_perm_exp is loaded.
        obj_cp_perm = 0
        Count_cp_perm = 0
        res_cp_perm_assemble = []
        _is_data3 = (str(workflow_family).upper() == "DATA3")
        _use_perm_probe = bool(_is_data3 and NF270_USE_PERMEATE_PROBE)
        _perm_rel_scale = (
            float(NF270_USE_PERMEATE_PROBE)
            if isinstance(NF270_USE_PERMEATE_PROBE, (int, float)) and NF270_USE_PERMEATE_PROBE
            else 0.03
        )
        _perm_floor_mM = (
            float(NF270_PERMEATE_PROBE_FLOOR_MM)
            if NF270_PERMEATE_PROBE_FLOOR_MM is not None else 0.0
        )

        # Number of vials that contribute to the fitted objective.
        collect_vial = data_stru['data_config']['n'] - data_stru['data_config']['n_extra']
        Count_m = 0
        Count_cp = 0
        Count_cf = 0
        count_cf0 = 0
        t_delay = _data1_time_origin(data_stru)
        TF_list = [data_stru['data_raw'][i]['time'][-1]-t_delay for i in range(data_stru['data_config']['n'])]
        TF_dict = dict(zip(m.n_vial,TF_list)) # unscaled time elapse for each vial 
        TI_list = [data_stru['data_raw'][i]['time'][0]-t_delay for i in range(data_stru['data_config']['n'])]
        TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time for each vial
    
        for n_vial in m.n_vial:
            # Per-concentration-band mask: skip vials outside the requested band
            # so theta is fit to one band only.  No-op when band_vials is None.
            if band_vials is not None and n_vial not in band_vials:
                continue
            t_meas = _data1_shifted_time(data_stru, n_vial-1)
            mv_meas = np.asarray(data_stru['data_raw'][n_vial-1]['mass'], dtype=float)
            cp_meas = np.asarray(data_stru['data_raw'][n_vial-1]['cV_avg'], dtype=float)
            cf_meas = np.asarray(data_stru['data_raw'][n_vial-1]['cF_exp'], dtype=float)

            t_meas_scaled = [(t-TI_dict[n_vial])/(TF_dict[n_vial]-TI_dict[n_vial]) for t in t_meas]

            mv_pred=[]
            res_m_=[]
            obj_mi = 0
            ob_mi = 0
            count_m = 0
            for i in range(0,len(t_meas_scaled)):
                mv_inter = interpolation(m,m.mV,n_vial,t_meas_scaled[i])
                mv_pred.append(mv_inter)
                
                if not np.isnan(mv_meas[i]):
                    res_m = mv_pred[i]-mv_meas[i]            
                    count_m += 1
                    res_m_.append(res_m/mass_scale)
                    obj_mi += (res_m/mass_scale)**2
                    ob_mi += res_m**2
            if n_vial >= data_stru['data_config']['n_v0']: 
                Count_m += count_m
                res_m_assemble.extend(res_m_)
                obj_m += obj_mi#/count_m/collect_vial
                ob_m += ob_mi
    
            # permeate residual squared / uncertainty / # of measurements
            if np.asarray(data_stru['data_raw'][n_vial-1]['cV_avg']).ndim == 0:
                if n_vial > data_stru['data_config']['n_extra']:
                    res_cp = m.cV[n_vial,m.tau.last()]-data_stru['data_raw'][n_vial-1]['cV_avg']
                    res_cp_assemble.append(res_cp/(_cp_frac*data_stru['data_raw'][n_vial-1]['cV_avg']))
                    obj_cp += (res_cp/(_cp_frac*data_stru['data_raw'][n_vial-1]['cV_avg']))**2 # permeate cV weight (DATA3: NF270_CP_RESIDUAL_SCALE_FRACTION; else 3%)
                    ob_cp += (res_cp/(data_stru['data_raw'][n_vial-1]['cV_avg']))**2
                    Count_cp = collect_vial
            else:
                cp_pred=[]
                obj_cpi = 0
                ob_cpi = 0
                count_cp = 0
                if not np.all(np.isnan(cp_meas)):
                    for i in range(0,len(t_meas_scaled)):
                        cp_inter = interpolation(m,m.cV,n_vial,t_meas_scaled[i])
                        cp_pred.append(cp_inter)
                        if not np.isnan(cp_meas[i]):   
                            res_cp = cp_pred[i]-cp_meas[i]
                            count_cp += 1
                            res_cp_assemble.append(res_cp/(_cp_frac*cp_meas[i]))
                            obj_cpi += (res_cp/(_cp_frac*cp_meas[i]))**2 # permeate cV weight (DATA3: NF270_CP_RESIDUAL_SCALE_FRACTION; else 3%)
                            ob_cpi += (res_cp/(cp_meas[i]))**2
                    Count_cp = collect_vial
                    obj_cp += obj_cpi#/count_cp/collect_vial
                    ob_cp += ob_cpi
         
            # retentate residual squared / uncertainty / # of measurements
            cf_pred=[]
            obj_cfi = 0
            ob_cfi = 0
            count_cf = 0
            if not np.all(np.isnan(cf_meas)):
                # Optional absolute floor on the cF residual scale (mM).
                # When NF270_CF_RESIDUAL_FLOOR_MM is set AND workflow_family
                # is DATA3, use:
                #     scale = max(0.003 * cf_meas[i], floor)
                # GUARDED to DATA3 only — DATA1/DATA2 always use the legacy
                # 0.3 % relative-only weight regardless of the constant's
                # value, so byte-equivalent published-paper behavior is
                # preserved even if a sibling script flips the toggle.
                cf_floor_mM = (
                    NF270_CF_RESIDUAL_FLOOR_MM
                    if str(workflow_family).upper() == "DATA3"
                    else None
                )
                for i in range(0,len(t_meas_scaled)):
                    cf_inter = interpolation(m,m.cF,n_vial,t_meas_scaled[i])
                    cf_pred.append(cf_inter)
                    if not np.isnan(cf_meas[i]):
                        res_cf = cf_pred[i]-cf_meas[i]
                        count_cf += 1
                        relative_scale = _cf_frac * cf_meas[i]  # retentate cF weight (DATA3: NF270_CF_RESIDUAL_SCALE_FRACTION=2%; else 0.3%)
                        if cf_floor_mM is not None and cf_floor_mM > 0 and relative_scale < cf_floor_mM:
                            cf_scale = float(cf_floor_mM)
                        else:
                            cf_scale = relative_scale
                        res_cf_assemble.append(res_cf/cf_scale)
                        obj_cfi += (res_cf/cf_scale)**2 # cF error: 0.3% relative or floored absolute
                        ob_cfi += (res_cf/(cf_meas[i]))**2
                if n_vial >= data_stru['data_config']['n_v0']:
                    Count_cf += count_cf
                    obj_cf += obj_cfi#/count_cf/collect_vial
                    ob_cf += ob_cfi
                else:
                    count_cf0 += count_cf
                    obj_cf0 += obj_cfi#/count_cf/collect_vial
                    ob_cf0 += ob_cfi

            # ---- Optional permeate-probe channel (cV_perm_exp continuous) ----
            # GUARDED:  only runs when (a) workflow_family == "DATA3", (b) the
            # module-level toggle NF270_USE_PERMEATE_PROBE is truthy, and
            # (c) the data carries cV_perm_exp.  Otherwise the block is a
            # complete no-op and the objective is byte-equivalent to the
            # legacy three-channel form.
            #
            # PHYSICS:  the in-line permeate probe measures the INSTANTANEOUS
            # permeate-stream concentration at the membrane outlet — this
            # corresponds to the model's `cH` (wall concentration on the
            # permeate side), NOT to `cV` (vial-integrated average).  Compare
            # cV_perm_exp (probe) to m.cH at the matching time, not to m.cV.
            if _use_perm_probe and n_vial >= data_stru['data_config']['n_v0']:
                cp_perm_meas = np.asarray(
                    data_stru['data_raw'][n_vial-1].get('cV_perm_exp', []),
                    dtype=float,
                )
                if cp_perm_meas.size == len(t_meas_scaled) and cp_perm_meas.size > 0:
                    obj_cp_perm_i = 0.0
                    count_cp_perm = 0
                    for i in range(0, len(t_meas_scaled)):
                        if not np.isfinite(cp_perm_meas[i]) or cp_perm_meas[i] <= 0:
                            continue
                        # cH = instantaneous permeate-side wall concentration
                        cH_inter = interpolation(m, m.cH, n_vial, t_meas_scaled[i])
                        res_cp_perm = cH_inter - cp_perm_meas[i]
                        relative_scale = _perm_rel_scale * cp_perm_meas[i]
                        if _perm_floor_mM > 0 and relative_scale < _perm_floor_mM:
                            scale = _perm_floor_mM
                        else:
                            scale = relative_scale
                        res_cp_perm_assemble.append(res_cp_perm / scale)
                        obj_cp_perm_i += (res_cp_perm / scale) ** 2
                        count_cp_perm += 1
                    Count_cp_perm += count_cp_perm
                    obj_cp_perm += obj_cp_perm_i

        final_tube_meas = data_stru["data_config"].get("cF_final_meas", np.nan)
        if not np.isfinite(final_tube_meas):
            final_tube_meas = data_stru["data_config"].get("icp_final_tube_mM", np.nan)
        if np.isfinite(final_tube_meas):
            final_scale = _abs_scale(final_tube_meas, _cf_frac, 0.003)  # cF anchor: DATA3 measurement error
            res_cf_final = m.cF[m.n_vial.last(), m.tau.last()] - final_tube_meas
            res_cf_assemble.append(res_cf_final / final_scale)
            obj_cf += (res_cf_final / final_scale) ** 2
            final_denom = final_tube_meas if np.isfinite(final_tube_meas) and abs(final_tube_meas) > 0 else final_scale
            ob_cf += (res_cf_final / final_denom) ** 2
            Count_cf += 1

        # ── Optional SECOND end-of-experiment anchor (Retentate ICP) ──
        # Added 2026-05-27 per user audit finding: cF_retentate_icp_mM
        # (bulk stirred-cell concentration at experiment end) was loaded
        # but never used as a fit constraint — only Final Tube ICP was.
        # When NF270_USE_RETENTATE_ICP_ANCHOR is truthy AND workflow_family
        # is DATA3, add a parallel residual term to obj_cr using the
        # Retentate ICP measurement.  GUARDED so DATA1/DATA2 remain
        # byte-equivalent to the legacy code path regardless of the
        # constant's value.
        _is_data3 = (str(workflow_family).upper() == "DATA3")
        if _is_data3 and NF270_USE_RETENTATE_ICP_ANCHOR:
            retentate_icp_meas = data_stru["data_config"].get("cF_retentate_icp_mM", np.nan)
            try:
                retentate_icp_meas = float(retentate_icp_meas)
            except (TypeError, ValueError):
                retentate_icp_meas = float("nan")
            if np.isfinite(retentate_icp_meas):
                # Use the same 0.3 % relative scaling as the Final Tube anchor.
                retentate_scale = _abs_scale(retentate_icp_meas, _cf_frac, 0.003)  # retentate cF anchor: DATA3 measurement error
                res_cf_retentate = m.cF[m.n_vial.last(), m.tau.last()] - retentate_icp_meas
                res_cf_assemble.append(res_cf_retentate / retentate_scale)
                obj_cf += (res_cf_retentate / retentate_scale) ** 2
                ret_denom = retentate_icp_meas if abs(retentate_icp_meas) > 0 else retentate_scale
                ob_cf += (res_cf_retentate / ret_denom) ** 2
                Count_cf += 1

        m.count_m = Count_m
        m.count_cv = Count_cp
        m.count_cr = Count_cf
        m.count_cr0 = count_cf0
        m.count = Count_m+Count_cp+Count_cf+count_cf0
        
        m.res_m = res_m_assemble
        m.res_cp = res_cp_assemble
        m.res_cf = res_cf_assemble

        # Guarded denominators for per-band masked fits: identical to the raw
        # counts when band_vials is None (so the legacy objective is byte-
        # identical), but floored at 1 when a band excludes a whole channel so
        # the per-channel normalization never divides by zero.
        _den_m  = Count_m if band_vials is None else max(Count_m, 1)
        _den_cp = Count_cp if band_vials is None else max(Count_cp, 1)
        _den_cr = (count_cf0 + Count_cf) if band_vials is None else max(count_cf0 + Count_cf, 1)
        _den_cf = Count_cf if band_vials is None else max(Count_cf, 1)

        m.obj_m = obj_m/_den_m
        m.obj_cv = obj_cp/_den_cp
        m.obj_cr = (obj_cf0+obj_cf)/_den_cr
        m.obj_cr_tru = obj_cf/_den_cf
        # Optional permeate-probe channel: stored on m for reporting even
        # when inactive (will be 0/None then).
        if _use_perm_probe and Count_cp_perm > 0:
            m.obj_cv_perm = obj_cp_perm / Count_cp_perm
        else:
            m.obj_cv_perm = 0.0
        # Overall normalized objective used by the solver.
        m.obj_tru = 1e4*(obj_m/_den_m+obj_cp/_den_cp+obj_cf/_den_cf)
        m.llh1 = m.count_m*log(m.obj_m) + m.count_cv*log(m.obj_cv) + (m.count_cr0+m.count_cr)*log(m.obj_cr)#m.count * log((obj_m+obj_cp+obj_cf0+obj_cf)/m.count)
        m.llh2 = Count_m*log(ob_m/_den_m) + Count_cp*log(ob_cp/_den_cp) + (count_cf0+Count_cf)*log((ob_cf0+ob_cf)/_den_cr)

        if mass_litmus_test:
            # Mass-litmus test: collapse the WSSE to mass only. Drop the cV
            # and cF residual terms; ignore Count_cp / Count_cf to avoid
            # divide-by-zero. One-parameter Lp fit against vial mass slopes.
            return 1e4 * (obj_m / max(Count_m, 1))
        # Legacy three-channel WSSE.  The optional fourth-channel
        # (continuous permeate-probe) is added ONLY when the DATA3 toggle
        # is active AND samples were collected; otherwise the objective is
        # byte-identical to the published form.
        legacy_obj = 1e4*(obj_m/_den_m + obj_cp/_den_cp + (obj_cf0+obj_cf)/_den_cr)
        if _use_perm_probe and Count_cp_perm > 0:
            legacy_obj = legacy_obj + 1e4 * (obj_cp_perm / Count_cp_perm)
        return legacy_obj

    # pyomo model instance
    instance = model_construct_inter(data_stru, mode, theta, sim_opt, B_form, workflow_family=workflow_family)
    #instance.pprint()
    if mass_litmus_test and not sim_opt and hasattr(instance, "sigma"):
        # Pin sigma = 0 so the flux equation collapses to Jw = Lp × ΔP
        # (no osmotic correction). The Var still exists in the constraint
        # graph; we just fix it.
        try:
            instance.sigma.fix(0.0)
            instance.sigma.setlb(0.0)
            instance.sigma.setub(0.0)
        except AttributeError:
            pass
    if sim_opt:
        instance.Obj_1 = Objective(expr = 1)
        instance.Obj = Expression(rule=obj_rule)
    else:
        instance.Obj = Objective(rule=obj_rule, sense=minimize)

    # Try initialize
    if skip_sim_init:
        # DATA3 opt-in (default False -> DATA1/DATA2 byte-identical): skip the uncapped
        # casadi/idas Simulator init, which grinds on stiff/fragile cases (e.g. LaCl3)
        # regardless of seed.  Discretize and let CPU-capped IPOPT solve from the seed.
        TransformationFactory('dae.finite_difference').apply_to(instance, nfe=nfe, scheme='BACKWARD')
    else:
        try:
            # Simulate the model using scipy
            sim = Simulator(instance, package='casadi')
            tsim, profiles = sim.simulate(numpoints=300, integrator='idas')
            # Discretize the model using finite_difference
            TransformationFactory('dae.finite_difference').apply_to(instance, nfe=nfe, scheme='BACKWARD')
            # Initialize the discretized model using the simulator profiles
            sim.initialize_model()
        except Exception as e:
            print(f"Initialization failed: {e}. Applying discretization without initialization.")
            TransformationFactory('dae.finite_difference').apply_to(instance, nfe=nfe, scheme='BACKWARD')

    solver = _build_ipopt_solver()
    cpu_time_limit = _resolve_ipopt_cpu_time_limit(
        solver_max_cpu_time=solver_max_cpu_time,
        multistart=multistart,
    )
    if cpu_time_limit is not None:
        solver.options["max_cpu_time"] = cpu_time_limit
        print(f"max_cpu_time={cpu_time_limit}")
    solver.options["max_iter"] = solver_max_iter
    #solver.options["halt_on_ampl_error"] = "yes" # option 1
    #solver.options["print_level"] = 1 # option 2
    try:
        results = solver.solve(instance,tee=True)
        term = results.solver.termination_condition
        if term == TerminationCondition.optimal:
            results.write()
        elif _is_infeasible_termination(term):
            print(f"Solver termination: {term}. Skipping repeated solve.")
            return None, None, None
        else:
            print(f"Solver termination: {term}, resolving...")
            try:
                solver.options["linear_solver"] = "mumps"
                solver.options["max_iter"] = solver_retry_max_iter
                results = solver.solve(instance,tee=True)
                assert results.solver.termination_condition == TerminationCondition.optimal, (
                        f"Solver failed: Non-optimal termination condition "
                        f"{results.solver.termination_condition}. ")
            except Exception as exc:
                print(f"Second solve failed: {exc}")
                return None, None, None
    except Exception as exc:
        print("Retrying with adjusted solver settings...")
        solver.options["linear_solver"] = "mumps"
        solver.options["max_iter"] = solver_retry_max_iter
        try:
            results = solver.solve(instance,tee=True)
            assert results.solver.termination_condition == TerminationCondition.optimal, (
                    f"Solver failed again: Non-optimal termination condition "
                    f"{results.solver.termination_condition}. "
                    "Try alternative initialization or further debugging.")
        except Exception as exc2:
            print(f"Adjusted solve failed: {exc2}")
            return None, None, None

    # Process results if optimal
    if results.solver.termination_condition == TerminationCondition.optimal:
        if mode !='DATA':
            fit_stru, sim_stru=save_model(instance, LO=True, B_form=B_form, LOUD=True)
        else:
            fit_stru, sim_stru=save_model(instance, LO=False, B_form=B_form, LOUD=True)
        sim_inter = inter_model(data_stru, sim_stru, fit_stru)
    else:
        print("Non-optimal solution after retries. Check initialization or parameter settings.")
        fit_stru, sim_stru, sim_inter = None, None, None
        return fit_stru, sim_stru, sim_inter

    if LOUD:
        time = []
        cIn = []
        cH = []
        Jw = []
        Js = []       
        for i in range(data_stru['data_config']['n']):   
            time.extend(sim_stru[i]['time'])
            cIn.extend(sim_stru[i]['cIn'])
            cH.extend(sim_stru[i]['cH'])
            Jw.extend(sim_stru[i]['Jw'])
            Js.extend(sim_stru[i]['Js'])
        
        sim_data = {'time': time,
                'cIn': cIn,
                'cH': cH,
                'Jw': Jw,
                'Js': Js}
        print(sim_data)
        fname = 'sim_data-dat'+str(data_stru['dataset'])       
        #create data frame from dictionary
        sim_datapd = pd.DataFrame(sim_data)
        
        #save dataframe to csv file
        sim_datapd.to_csv(fname+".csv", index=False)
        
        #validate the csv file by importing it
        #print(pd.read_csv(fname+".csv"))

    # Finish print statement
    print("###################################################################")

    return fit_stru, sim_stru, sim_inter


def nested_dict_values(dict_nest):
    ''' 
    This function accepts a nested dictionary as argument
        and iterates over all values of nested dictionaries
    '''
    dict_value=[]
    # Iterate over all values of given dictionary
    for value in dict_nest.values():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over all its values
            for v in  nested_dict_values(value):
                dict_value.append(v)
        else:
            # If value is not dict type then append the value
            dict_value.append(value)
    return dict_value

def nested_dict_keys(dict_nest):
    ''' 
    This function accepts a nested dictionary as argument
        and iterates over all keys of nested dictionaries
    '''   
    dict_key=[]
    # Iterate over all values of given dictionary
    for key in dict_nest.keys():
        # Check if value is of dict type
        if isinstance(dict_nest[key], dict):           
            # If value is dict then iterate over all its keys
            for v in nested_dict_keys(dict_nest[key]):
                dict_key.append(key+'-'+str(v))
        else:
            # If value is not dict type then yield the value
            dict_key.append(key)
    return dict_key

def nested_dict_update(dict_nest,new_value,i=0):
    ''' 
    This function accepts a nested dictionary and a new list of values as arguments
        and updates all values of nested dictionaries
    '''    
    # Iterate over all keys of given dictionary
    for key in dict_nest.keys():
        # Check if value is of dict type
        if isinstance(dict_nest[key], dict):
            # If value is dict then iterate over all its keys
            for v in nested_dict_update(dict_nest[key],new_value,i=i):
                v = new_value[i]
                i += 1
        else:
            # If value is not dict type then update the value
            dict_nest[key] = new_value[i]
            i += 1
    return dict_nest


def calc_FIM(data_stru, mode, theta=None, step=1e-8, formula='backward', B_form='single', workflow_family='DATA1', nfe=300, band_vials=None):
    """
    Calculate FIM

    Arguments:
        data_stru: dict, experimental data dictionary
        mode: str, experiment mode, {DATA, lag, overflow}
        theta: dict, preset parameter values
        step: float, relative step change for parameters
        formula: str, finite difference scheme option, {backward, forward, central}
        B_form: float/str, different solute permeability coefficient formula
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration
        band_vials: optional iterable of 1-indexed vial numbers.  When given, the
                FIM is built from ONLY those vials' predictions (both the
                measurement covariance and the Jacobian are masked identically),
                giving a per-concentration-band FIM whose non-singularity tells
                you whether that band alone identifies theta.  None (default) uses
                all vials and is byte-identical to the legacy FIM.

    Returns:
        doe_stru: dict, FIM results
    """
    formula = str(formula).lower()
    data_stru = _normalize_conductivity_measurements(data_stru)
    workflow_tag = str(workflow_family).upper()
    sim_opt = True
    # Per-band FIM mask (1-indexed vials -> 0-indexed loop index n_vial+1).
    _band = set(int(v) for v in band_vials) if band_vials is not None else None

    def _solve_fim_case(label, theta_guess, *, sim_opt_flag):
        solve_kwargs = dict(
            theta=theta_guess,
            sim_opt=sim_opt_flag,
            B_form=B_form,
            workflow_family=workflow_family,
            nfe=nfe,
        )
        if workflow_tag == "DATA2":
            solve_kwargs["solver_max_iter"] = 3000
            solve_kwargs["solver_retry_max_iter"] = 5000
            return _safe_solve_case_multistart(
                label,
                solve_model,
                data_stru,
                mode,
                **solve_kwargs,
            )
        return solve_model(
            data_stru,
            mode,
            **solve_kwargs,
        )

    if theta is None:
        sim_opt = False
        fit_stru_p, sim_stru_p, sim_inter_p = _solve_fim_case(
            f"{workflow_tag or 'DATA'} FIM base fit",
            theta,
            sim_opt_flag=sim_opt,
        )
        if fit_stru_p is None or sim_inter_p is None:
            if workflow_tag == "DATA2":
                print("Skipping DATA2 FIM because the base fit did not converge.")
                return None
            raise RuntimeError("Base solve failed while building the FIM.")
        theta = fit_stru_p['parameters']
        sim_opt = True
    theta_p = copy.deepcopy(theta)
    theta_p_v=nested_dict_values(theta_p)
    theta_p1_v = []
    theta_p2_v = []

    # Build finite-difference perturbations for each parameter.
    for i in theta_p_v: 
        
        if formula == 'central':
            if i!=0:
                theta_p1_v.append(i * (1-step))
                theta_p2_v.append(i * (1+step))
            else:
                theta_p1_v.append(-step)
                theta_p2_v.append(step)
        elif formula == 'backward':
            if i!=0:
                theta_p1_v.append(i * (1-step))
            else:
                theta_p1_v.append(-step)
            theta_p2_v.append(i)
        else: #forward
            theta_p1_v.append(i)
            if i!=0:
                theta_p2_v.append(i * (1+step))
            else:
                theta_p2_v.append(step)
    
    # Prediction covariance is built from the per-channel measurement errors:
    # 0.01 g for mass; permeate cV and retentate cF use _nf270_conc_scales() so
    # the FIM covariance MATCHES the objective weights (DATA3: cF 2%, cV 3%;
    # DATA1/DATA2: cF 0.3%, cV 3%).  Must stay in sync with obj_rule.
    _fim_cf_frac, _fim_cp_frac = _nf270_conc_scales(workflow_family)
    if sim_opt == True:
        fit_stru_p, sim_stru_p, sim_inter_p = _solve_fim_case(
            f"{workflow_tag or 'DATA'} FIM reference solve",
            theta_p,
            sim_opt_flag=sim_opt,
        )
        if fit_stru_p is None or sim_inter_p is None:
            if workflow_tag == "DATA2":
                print("Skipping DATA2 FIM because the reference simulation did not converge.")
                return None
            raise RuntimeError("Reference solve failed while building the FIM.")
    var_pred=[]
    for n_vial in range(data_stru['data_config']['n']):
        if _band is not None and (n_vial + 1) not in _band:
            continue
        var_pred = np.append(var_pred,0.01 ** 2 * np.ones(len(sim_inter_p[n_vial]['mV'])))
        var_pred = np.append(var_pred,(_fim_cp_frac * sim_inter_p[n_vial]['cV'])**2)
        var_pred = np.append(var_pred,(_fim_cf_frac * sim_inter_p[n_vial]['cF'])**2)
    cov_pred = np.diag(var_pred)
    
    # The Jacobian is the sensitivity of all predictions with respect to all parameters.
    Jac = []
    for i in range(len(theta_p_v)):
        theta_pk1_v = theta_p_v.copy()
        theta_pk2_v = theta_p_v.copy()
        theta_pk1_v[i] = theta_p1_v[i]    
        theta_pk2_v[i] = theta_p2_v[i]       
        theta_pk1 = copy.deepcopy(theta_p)
        theta_pk2 = copy.deepcopy(theta_p)
        theta_pk1 = nested_dict_update(theta_pk1,theta_pk1_v)
        theta_pk2 = nested_dict_update(theta_pk2,theta_pk2_v)
        print(theta_pk1)
        print(theta_pk2)
        
        sim_opt = True
        fit_stru_pk1, sim_stru_pk1, sim_inter_pk1 = _solve_fim_case(
            f"{workflow_tag or 'DATA'} FIM perturbation {i+1}a",
            theta_pk1,
            sim_opt_flag=sim_opt,
        )
        fit_stru_pk2, sim_stru_pk2, sim_inter_pk2 = _solve_fim_case(
            f"{workflow_tag or 'DATA'} FIM perturbation {i+1}b",
            theta_pk2,
            sim_opt_flag=sim_opt,
        )
        if sim_inter_pk1 is None or sim_inter_pk2 is None:
            if workflow_tag == "DATA2":
                print(f"Skipping DATA2 FIM because perturbation solve {i+1} did not converge.")
                return None
            raise RuntimeError(f"Perturbation solve {i+1} failed while building the FIM.")
        jac=[]
        for n_vial in range(data_stru['data_config']['n']):
            if _band is not None and (n_vial + 1) not in _band:
                continue
            jac = np.append(jac,sim_inter_pk2[n_vial]['mV']-sim_inter_pk1[n_vial]['mV'])
            jac = np.append(jac,sim_inter_pk2[n_vial]['cV']-sim_inter_pk1[n_vial]['cV'])
            jac = np.append(jac,sim_inter_pk2[n_vial]['cF']-sim_inter_pk1[n_vial]['cF'])
        # Predictions in row
        if theta_p_v[i]!=0:
            delta = step*theta_p_v[i]
        else:
            delta = step
        if len(Jac) == 0:
            Jac = jac/delta
        else:
            Jac = np.vstack([Jac, jac/delta])

    doe_stru = dict()
    if formula == 'central':
        doe_stru['Jac'] = Jac/2
    else:
        doe_stru['Jac'] = Jac
    # Fisher Information Matrix:
    # FIM = J * Cov(y)^(-1) * J^T
    FIM = doe_stru['Jac'] @ np.linalg.inv(cov_pred) @ doe_stru['Jac'].T
    doe_stru['Jac'] = doe_stru['Jac'].tolist()
    doe_stru['FIM'] = FIM.tolist()

    # Eigenvalues and eigenvectors summarize identifiability and uncertainty directions.
    w, v = np.linalg.eigh(FIM)
    doe_stru['eig_val'] = w.tolist()
    doe_stru['eig_vec'] = v.tolist() # in columns
    theta_p_name=nested_dict_keys(theta_p)
    maxind = np.argmax(abs(v), axis=0)
    doe_stru['eig_dir']=[theta_p_name[i] for i in (maxind)]
    doe_stru['trace'] = np.trace(FIM).tolist()
    doe_stru['det'] = np.linalg.det(FIM).tolist()
    doe_stru['min_eig'] = min(w)
    doe_stru['cond'] = max(w) / min(w)
    try:
        doe_stru['V'] = np.linalg.inv(FIM).tolist()
        doe_stru['std'] = np.sqrt(np.diag(doe_stru['V'])).tolist()
    except:
        print('No inverse of FIM')

    print(doe_stru)
    return doe_stru


def _vial_terminal_cf(data_stru, n_vial):
    """Last finite cF_exp sample for a 1-indexed vial (its terminal retentate
    concentration, used to assign the vial to a concentration band)."""
    a = np.asarray(data_stru['data_raw'][n_vial - 1]['cF_exp'], dtype=float)
    a = a[np.isfinite(a)]
    return float(a[-1]) if a.size else float('nan')


def partition_vials_by_terminal_cf(data_stru, n_bands=2):
    """Partition the FITTED vials (n_vial >= n_v0) into ``n_bands`` contiguous
    concentration bands by terminal cF, low->high, each band a sorted list of
    1-indexed vial numbers.  A median split (n_bands=2) yields balanced low/high
    bands; vials with no finite cF are dropped."""
    n = data_stru['data_config']['n']
    n_v0 = data_stru['data_config'].get('n_v0', 1)
    fitted = [i for i in range(1, n + 1) if i >= n_v0]
    pairs = sorted((_vial_terminal_cf(data_stru, i), i)
                   for i in fitted if np.isfinite(_vial_terminal_cf(data_stru, i)))
    if not pairs:
        return []
    order = [i for _, i in pairs]
    groups = np.array_split(np.asarray(order), max(1, int(n_bands)))
    return [sorted(int(x) for x in g) for g in groups if len(g) > 0]


def solve_model_per_concentration_band(
    data_stru,
    mode,
    *,
    bands=None,
    n_bands=2,
    B_form='single',
    workflow_family='DATA3',
    nfe=80,
    seed_theta=None,
    compute_fim=True,
    min_vials=6,
    max_bands=2,
    solver_max_cpu_time=200,
    LOUD=False,
):
    """Fit one theta per concentration band by MASKING the WSSE objective to each
    band's vials (Architecture.md §18.4).

    Bands partition the fitted vials by terminal cF (median/quantile split via
    ``partition_vials_by_terminal_cf``) unless an explicit ``bands`` list of
    1-indexed vial-number lists is supplied.  Each band is fit with
    ``solve_model(..., band_vials=band)`` warm-started from a full-sheet seed
    fit, and (when ``compute_fim``) its per-band FIM is checked for
    non-singularity — the statistical guard from §18.4.

    Restricted to >= ``min_vials`` fitted vials and <= ``max_bands`` bands; a
    sheet that violates either still runs but the returned dict carries a
    ``warnings`` list so callers can refuse to trust the per-band UQ.

    Returns
    -------
    dict with keys: run_mode, n_fitted, bands, warnings, full (seed theta +
    objectives), per_band (list of {band_vials, cf_range, theta, obj_m, obj_cv,
    obj_cr, fim}).
    """
    n = data_stru['data_config']['n']
    n_v0 = data_stru['data_config'].get('n_v0', 1)
    fitted = [i for i in range(1, n + 1) if i >= n_v0]

    if bands is None:
        bands = partition_vials_by_terminal_cf(data_stru, n_bands=n_bands)
    bands = [sorted(int(v) for v in b) for b in bands if len(b) > 0]

    warnings = []
    if len(fitted) < min_vials:
        warnings.append(
            f"only {len(fitted)} fitted vials (< min_vials={min_vials}); per-band "
            f"theta is weakly determined — treat per-band UQ with caution")
    if len(bands) > max_bands:
        warnings.append(f"{len(bands)} bands > max_bands={max_bands}; capping to {max_bands}")
        bands = bands[:max_bands]
    for b in bands:
        if len(b) < 3:
            warnings.append(f"band {b} has < 3 vials; theta under-determined for that band")

    def _summ(fit):
        p = dict(fit.get('parameters', {}))
        return {
            'theta': p,
            'obj_m': fit.get('obj_m'), 'obj_cv': fit.get('obj_cv'), 'obj_cr': fit.get('obj_cr'),
        }

    # Full-sheet seed fit (also the comparison baseline).
    if seed_theta is None:
        full_fit, _, _ = solve_model(
            data_stru, mode, sim_opt=False, B_form=B_form,
            workflow_family=workflow_family, nfe=nfe,
            solver_max_cpu_time=solver_max_cpu_time, LOUD=LOUD)
        seed_theta = dict(full_fit['parameters'])
        full = _summ(full_fit)
    else:
        full = {'theta': dict(seed_theta), 'obj_m': None, 'obj_cv': None, 'obj_cr': None}

    per_band = []
    for bi, band in enumerate(bands):
        cfs = [_vial_terminal_cf(data_stru, i) for i in band]
        cfs = [c for c in cfs if np.isfinite(c)]
        entry = {
            'band_index': bi,
            'band_vials': band,
            'cf_range': [min(cfs), max(cfs)] if cfs else [None, None],
        }
        fit, _, _ = solve_model(
            data_stru, mode, theta=seed_theta, sim_opt=False, B_form=B_form,
            workflow_family=workflow_family, nfe=nfe, band_vials=band,
            solver_max_cpu_time=solver_max_cpu_time, LOUD=LOUD)
        entry.update(_summ(fit))
        if compute_fim:
            try:
                doe = calc_FIM(
                    data_stru, mode, theta=entry['theta'], B_form=B_form,
                    workflow_family=workflow_family, nfe=nfe, band_vials=band)
                if doe is None:
                    entry['fim'] = {'error': 'calc_FIM returned None'}
                else:
                    min_eig = float(doe.get('min_eig', float('nan')))
                    cond = float(doe.get('cond', float('nan')))
                    entry['fim'] = {
                        'min_eig': min_eig,
                        'cond': cond,
                        'det': float(doe.get('det', float('nan'))),
                        'singular': (not np.isfinite(min_eig)) or min_eig <= 0
                                    or (not np.isfinite(cond)) or cond > 1e12,
                    }
            except Exception as exc:
                entry['fim'] = {'error': repr(exc)}
        per_band.append(entry)

    return {
        'run_mode': mode,
        'B_form': B_form,
        'n_fitted': len(fitted),
        'bands': bands,
        'warnings': warnings,
        'full': full,
        'per_band': per_band,
    }


def correlation_from_covariance(covariance):
    ''' 
    This function calculate correlation matrix from covariance matrix.
    '''    
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    correlation[covariance == 0] = 0
    print(correlation)
    return correlation

def store_json(file_name,structure):
    with open(file_name, "w") as json_file:
        json.dump(structure, json_file)
        
        
def solve_model_B_fix(
    data_stru,
    mode,
    theta=None,
    sim_opt=False,
    B_form=1,
    sigma_fixed=True,
    LOUD=False,
    workflow_family='DATA1',
    nfe=300,
    solver_max_iter=3000,
    solver_retry_max_iter=5000,
    solver_max_cpu_time=None,
    fix_vars=None,
    skip_sim_init=False,
):
    """
    Solve pyomo model
    
    Arguments:
        data_stru: dict, experimental data dictionary
        mode: str, experiment mode, {DATA, lag, overflow}
        theta: dict, preset parameter values
        sim_opt: boolean, if run simulation with fixed parameter values
        B_form: float/str, different solute permeability coefficient formula
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration
        sigma_fixed: boolean, if fix sigma value
        LOUD: boolean, if print out parameter results and store model predictions
    
    Returns:
        fit_stru: dict, parameter fit results
        sim_stru: dict, model predictions
        sim_inter: dict, model predictions for experimental measurements 
    """

    print("\n\n###################################################################")
    print("Creating and solving the Pyomo model with the following settings: ")
    print("mode =", mode)
    print("theta =",theta)
    print("sim_opt =", sim_opt)
    print("B_form =", B_form)
    print("sigma_fixed =",sigma_fixed)
    print(" ")

    data_stru = _normalize_conductivity_measurements(data_stru)
    if str(workflow_family).upper() == "DATA2" and solver_max_iter == 3000:
        # DATA2 contains a few intentionally hard / locally infeasible cases.
        # Fail fast on those so the paper workflow can keep moving.
        solver_max_iter = 300
        solver_retry_max_iter = min(solver_retry_max_iter, 600)

    # interpolation function
    def interpolation(m,var,n_vial,t):
        """
        Interpolation function for pyomo model
        """
        tp = list(m.tau)
        for j in range(0,len(tp)-1):
            if tp[j+1]>=t and tp[j]<=t:
                var_inter = (var[n_vial,tp[j+1]]-var[n_vial,tp[j]]) * (t - tp[j])/(tp[j+1]-tp[j]) + var[n_vial,tp[j]]
        return var_inter
    
    # objective function
    def obj_rule(m):
        """
        pyomo Objective function (residual calculated using interpolation from model)
        """
        # Salt-aware mass normalization.  The current loader keeps the
        # DATA2-style 0.01 g weight in data_config, so the fit stays aligned
        # with the published convention unless a caller deliberately overrides
        # it for an experiment.
        mass_scale = float(data_stru["data_config"].get("mass_scale_g", 0.01))
        if not np.isfinite(mass_scale) or mass_scale <= 0:
            mass_scale = 0.01
        # DATA3-only paper-justified retentate/permeate concentration sigmas
        # (legacy 0.003 / 0.03 for DATA1/DATA2 and unset knobs).
        _cf_frac, _cp_frac = _nf270_conc_scales(workflow_family)
        obj_m = 0
        obj_cp = 0
        obj_cf0 = 0
        obj_cf = 0
        ob_m = 0
        ob_cp = 0
        ob_cf0 = 0
        ob_cf = 0
        res_m_assemble = []
        res_cp_assemble = []
        res_cf_assemble = []

        collect_vial = data_stru['data_config']['n'] - data_stru['data_config']['n_extra']
        Count_m = 0
        Count_cp = 0
        Count_cf = 0
        count_cf0 = 0
        t_delay = _data1_time_origin(data_stru)
        TF_list = [data_stru['data_raw'][i]['time'][-1]-t_delay for i in range(data_stru['data_config']['n'])]
        TF_dict = dict(zip(m.n_vial,TF_list)) # unscaled time elapse for each vial
        TI_list = [data_stru['data_raw'][i]['time'][0]-t_delay for i in range(data_stru['data_config']['n'])]
        TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time for each vial

        for n_vial in m.n_vial:
            t_meas = _data1_shifted_time(data_stru, n_vial-1)
            mv_meas = np.asarray(data_stru['data_raw'][n_vial-1]['mass'], dtype=float)
            cp_meas = np.asarray(data_stru['data_raw'][n_vial-1]['cV_avg'], dtype=float)
            cf_meas = np.asarray(data_stru['data_raw'][n_vial-1]['cF_exp'], dtype=float)

            t_meas_scaled = [(t-TI_dict[n_vial])/(TF_dict[n_vial]-TI_dict[n_vial]) for t in t_meas]

            mv_pred=[]
            res_m_=[]
            obj_mi = 0
            ob_mi = 0
            count_m = 0
            for i in range(0,len(t_meas_scaled)):
                mv_inter = interpolation(m,m.mV,n_vial,t_meas_scaled[i])
                mv_pred.append(mv_inter)

                if not np.isnan(mv_meas[i]):
                    res_m = mv_pred[i]-mv_meas[i]
                    count_m += 1
                    res_m_.append(res_m/mass_scale)
                    obj_mi += (res_m/mass_scale)**2  # salt-aware mass weight
                    ob_mi += res_m**2
            if n_vial >= data_stru['data_config']['n_v0']: 
                Count_m += count_m
                res_m_assemble.extend(res_m_)
                obj_m += obj_mi#/count_m/collect_vial
                ob_m += ob_mi
    
            # permeate residual squared / uncertainty / # of measurements
            if np.asarray(data_stru['data_raw'][n_vial-1]['cV_avg']).ndim == 0:
                if n_vial > data_stru['data_config']['n_extra']:
                    res_cp = m.cV[n_vial,m.tau.last()]-data_stru['data_raw'][n_vial-1]['cV_avg']
                    res_cp_assemble.append(res_cp/(_cp_frac*data_stru['data_raw'][n_vial-1]['cV_avg']))
                    obj_cp += (res_cp/(_cp_frac*data_stru['data_raw'][n_vial-1]['cV_avg']))**2 # permeate cV weight (DATA3: NF270_CP_RESIDUAL_SCALE_FRACTION; else 3%)
                    ob_cp += (res_cp/(data_stru['data_raw'][n_vial-1]['cV_avg']))**2
                    Count_cp = collect_vial
            else:
                cp_pred=[]
                obj_cpi = 0
                ob_cpi = 0
                count_cp = 0
                if not np.all(np.isnan(cp_meas)):
                    for i in range(0,len(t_meas_scaled)):
                        cp_inter = interpolation(m,m.cV,n_vial,t_meas_scaled[i])
                        cp_pred.append(cp_inter)
                        if not np.isnan(cp_meas[i]):   
                            res_cp = cp_pred[i]-cp_meas[i]
                            count_cp += 1
                            res_cp_assemble.append(res_cp/(_cp_frac*cp_meas[i]))
                            obj_cpi += (res_cp/(_cp_frac*cp_meas[i]))**2 # permeate cV weight (DATA3: NF270_CP_RESIDUAL_SCALE_FRACTION; else 3%)
                            ob_cpi += (res_cp/(cp_meas[i]))**2
                    Count_cp = collect_vial
                    obj_cp += obj_cpi#/count_cp/collect_vial
                    ob_cp += ob_cpi
         
            # retentate residual squared / uncertainty / # of measurements
            cf_pred=[]
            obj_cfi = 0
            ob_cfi = 0
            count_cf = 0
            if not np.all(np.isnan(cf_meas)):
                # Optional absolute floor on the cF residual scale (mM).
                # GUARDED to DATA3 only — see model_construct_inter for the
                # full rationale.  DATA1/DATA2 always use the legacy 0.3 %
                # relative weight regardless of the constant's value.
                cf_floor_mM = (
                    NF270_CF_RESIDUAL_FLOOR_MM
                    if str(workflow_family).upper() == "DATA3"
                    else None
                )
                for i in range(0,len(t_meas_scaled)):
                    cf_inter = interpolation(m,m.cF,n_vial,t_meas_scaled[i])
                    cf_pred.append(cf_inter)
                    if not np.isnan(cf_meas[i]):
                        res_cf = cf_pred[i]-cf_meas[i]
                        count_cf += 1
                        relative_scale = _cf_frac * cf_meas[i]  # retentate cF weight (DATA3: NF270_CF_RESIDUAL_SCALE_FRACTION=2%; else 0.3%)
                        if cf_floor_mM is not None and cf_floor_mM > 0 and relative_scale < cf_floor_mM:
                            cf_scale = float(cf_floor_mM)
                        else:
                            cf_scale = relative_scale
                        res_cf_assemble.append(res_cf/cf_scale)
                        obj_cfi += (res_cf/cf_scale)**2 # cF error: 0.3% relative or floored absolute
                        ob_cfi += (res_cf/(cf_meas[i]))**2
                if n_vial >= data_stru['data_config']['n_v0']:
                    Count_cf += count_cf
                    obj_cf += obj_cfi#/count_cf/collect_vial
                    ob_cf += ob_cfi
                else: 
                    count_cf0 += count_cf
                    obj_cf0 += obj_cfi#/count_cf/collect_vial
                    ob_cf0 += ob_cfi
        
        m.count_m = Count_m
        m.count_cv = Count_cp
        m.count_cr = Count_cf
        m.count_cr0 = count_cf0
        m.count = Count_m+Count_cp+Count_cf+count_cf0
        
        m.res_m = res_m_assemble
        m.res_cp = res_cp_assemble
        m.res_cf = res_cf_assemble
        
        m.obj_m = obj_m/Count_m
        m.obj_cv = obj_cp/Count_cp
        m.obj_cr = (obj_cf0+obj_cf)/(count_cf0+Count_cf)
        m.obj_cr_tru = obj_cf/Count_cf
        m.obj_tru = 1e4*(obj_m/Count_m+obj_cp/Count_cp+obj_cf/Count_cf)
        m.llh1 = m.count_m*log(m.obj_m) + m.count_cv*log(m.obj_cv) + (m.count_cr0+m.count_cr)*log(m.obj_cr)#m.count * log((obj_m+obj_cp+obj_cf0+obj_cf)/m.count)
        m.llh2 = Count_m*log(ob_m/Count_m) + Count_cp*log(ob_cp/Count_cp) + (count_cf0+Count_cf)*log((ob_cf0+ob_cf)/(count_cf0+Count_cf))
        
        return 1e4*(obj_m/Count_m + obj_cp/Count_cp + (obj_cf0+obj_cf)/(count_cf0+Count_cf))
    
    # pyomo model instance
    instance = model_construct_inter(data_stru, mode, theta, sim_opt, B_form, workflow_family=workflow_family)
    if fix_vars:
        # Profile-contour mode: fix exactly the named Vars at the given grid values
        # and leave every other parameter (incl. the B(c) coefficients) FREE to be
        # solved.  This is the inverse of the default (fix-B) behaviour and is used
        # to draw WSSE landscapes where B genuinely varies with concentration.
        for _name, _val in fix_vars.items():
            if hasattr(instance, _name):
                getattr(instance, _name).fix(float(_val))
    else:
        if B_form=='single':
            instance.B.fixed=True
        elif isinstance(B_form, str) and B_form in ('sat', 'donnan'):
            for _p in ('B_inf', 'c_star', 'P0', 'X', 'k_dd'):
                if hasattr(instance, _p):
                    getattr(instance, _p).fixed = True
        else:
            instance.beta_0.fixed=True
            instance.beta_1.fixed=True
            if not isinstance(B_form, str) and B_form > 1:
                instance.beta_2.fixed=True
        if sigma_fixed:
            instance.sigma.fixed=True
    if sim_opt:
        instance.Obj_1 = Objective(expr = 1)
        instance.Obj = Expression(rule=obj_rule)
    else:
        instance.Obj = Objective(rule=obj_rule, sense=minimize)

    # Try initialize
    if skip_sim_init:
        # DATA3 profile-contour path (opt-in only -> DATA1/DATA2 byte-identical): skip the
        # UNCAPPED casadi/idas Simulator init, which grinds for minutes on intrinsically
        # stiff (high-sigma) grid nodes regardless of the seed.  Discretize directly and
        # let the warm-started, CPU-capped IPOPT solve from the seed Var values; a stiff
        # node then fails FAST under the cap instead of hanging the integrator.
        TransformationFactory('dae.finite_difference').apply_to(instance, nfe=nfe, scheme='BACKWARD')
    else:
        try:
            # Simulate the model using scipy
            sim = Simulator(instance, package='casadi')
            tsim, profiles = sim.simulate(numpoints=300, integrator='idas')
            # Discretize the model using finite difference
            TransformationFactory('dae.finite_difference').apply_to(instance, nfe=nfe, scheme='BACKWARD')
            # Initialize the discretized model using the simulator profiles
            sim.initialize_model()
        except:
            TransformationFactory('dae.finite_difference').apply_to(instance, nfe=nfe, scheme='BACKWARD')

    solver = _build_ipopt_solver()
    cpu_time_limit = _resolve_ipopt_cpu_time_limit(
        solver_max_cpu_time=solver_max_cpu_time,
        multistart=False,
    )
    if cpu_time_limit is not None:
        solver.options["max_cpu_time"] = cpu_time_limit
        print(f"max_cpu_time={cpu_time_limit}")
    solver.options["max_iter"] = solver_max_iter

    results = solver.solve(instance,tee=True)
    term = results.solver.termination_condition
    if _is_infeasible_termination(term):
        print(f"Solver termination: {term}. Skipping case.")
        return None, None, None
    assert term == TerminationCondition.optimal, "Optimization fail"
    results.write()
    #solver.options["halt_on_ampl_error"] = "yes" # option 1
    #solver.options["print_level"] = 1 # option 2
    #solver.solve(instance,tee=True).write()
    #instance.display()

    if mode !='DATA':
        fit_stru, sim_stru=save_model(instance, LO=True, B_form=B_form, LOUD=True)
    else:
        fit_stru, sim_stru=save_model(instance, LO=False, B_form=B_form, LOUD=True)
    sim_inter = inter_model(data_stru, sim_stru, fit_stru)

    if LOUD:
        time = []
        cIn = []
        cH = []
        Jw = []
        Js = []       
        for i in range(data_stru['data_config']['n']):   
            time.extend(sim_stru[i]['time'])
            cIn.extend(sim_stru[i]['cIn'])
            cH.extend(sim_stru[i]['cH'])
            Jw.extend(sim_stru[i]['Jw'])
            Js.extend(sim_stru[i]['Js'])
        
        sim_data = {'time': time,
                'cIn': cIn,
                'cH': cH,
                'Jw': Jw,
                'Js': Js}
        print(sim_data)
        fname = 'sim_data-dat'+str(data_stru['dataset'])       
        #create data frame from dictionary
        sim_datapd = pd.DataFrame(sim_data)
        
        #save dataframe to csv file
        sim_datapd.to_csv(fname+".csv", index=False)
        
        #validate the csv file by importing it
        #print(pd.read_csv(fname+".csv"))

    print("###################################################################")

    return fit_stru, sim_stru, sim_inter


# -----------------------------------------------------------------------------
# =============================================================================
# Mass-balance-only solver (mass-litmus test)
# -----------------------------------------------------------------------------
#
# Closed-form, numpy-only solver for the simplest possible reduction of the
# diafiltration model:
#
#     sigma = 0      (no rejection)
#     B     = N/A    (concentration physics removed entirely)
#     Jw    = Lp × ΔP / 36000                              [algebraic]
#     dmV/dt = Jw × Am × ρ                                 [vial mass]
#     dmF/dt = −S₀ − Am × ρ × Jw                           [feed mass]
#     mV(t) − mV(t_start) = (Lp × ΔP × Am × ρ / 36000) × (t − t_start)
#
# The concentration ODEs (dcF/dt, dcH/dt, dcVmV/dt) are dropped from the
# model entirely — they couple to Jw only through sigma, which is zero here.
# What remains is a one-parameter (Lp) linear regression problem, convex
# with a single global minimum, solvable in closed form.
#
# No Pyomo, no IPOPT — np.polyfit-equivalent. The fit_stru / sim_stru output
# shape matches solve_model so render_nf270_fit_v2 (and any DATA1/DATA2
# renderer) consumes it without changes.
# =============================================================================

def solve_mass_balance_only(data_stru, *, B_form='single', LOUD=False, **_unused):
    """Closed-form mass-balance-only solver.

    Builds the same fit_stru / sim_stru shape that solve_model returns,
    but the underlying model is the stripped-down mass balance:
        mV(t) = (Lp × ΔP × Am × ρ / 36000) × (t − t_vial_start)
    Lp is fit by linear least-squares (zero-intercept) against the
    aggregated (Δt, mass) pairs from every "real" vial (n_vial >= n_v0;
    extras are skipped). The prediction trace for extra vials is set to
    zero so the renderer leaves a gap there.

    No concentration variables, no B, no sigma — those fields are present
    in the returned fit_stru only as None / 0.0 placeholders so the
    downstream consumer code doesn't KeyError.

    Returns
    -------
    fit_stru : dict
        Parameter fit results with diagnostic_mode='mass_balance_only'.
    sim_stru : list of dicts
        One dict per vial with 'time', 'mV', 'cF', 'cH', 'cV' arrays.
        cF / cH / cV are zero — concentrations weren't fit.
    """
    cfg     = data_stru['data_config']
    delP    = float(cfg['delP'])
    Am      = float(cfg['Am'])
    rho     = float(cfg['rho'])
    n_v0    = int(cfg.get('n_v0', 1))
    data_raw = data_stru.get('data_raw', [])
    if not data_raw:
        raise ValueError("solve_mass_balance_only: data_raw is empty.")

    # ---- Aggregate (Δt, mass) pairs across real vials -------------------
    dts, ms = [], []
    for i, row in enumerate(data_raw):
        if (i + 1) < n_v0:
            continue
        t = np.asarray(row.get('time', []), dtype=float).reshape(-1)
        m = np.asarray(row.get('mass', []), dtype=float).reshape(-1)
        if t.size == 0 or m.size == 0:
            continue
        n = min(t.size, m.size)
        t = t[:n]; m = m[:n]
        valid = np.isfinite(t) & np.isfinite(m)
        if not valid.any():
            continue
        t_start = float(t[valid][0])
        dts.extend((t[valid] - t_start).tolist())
        ms.extend(m[valid].tolist())

    dts = np.asarray(dts, dtype=float)
    ms  = np.asarray(ms,  dtype=float)
    if dts.size < 2:
        raise ValueError(
            "solve_mass_balance_only: need at least 2 valid (t, mass) "
            "points across real vials. Got {}.".format(dts.size))

    # Zero-intercept linear least-squares: slope = Σ(Δt·m) / Σ(Δt²)
    denom = float(np.sum(dts * dts))
    if denom <= 0:
        raise ValueError("solve_mass_balance_only: degenerate Δt window.")
    slope_g_per_s = float(np.sum(dts * ms) / denom)

    # Lp [L/(m²·hr·bar)] = slope × 36000 / (Am · ρ · ΔP)
    if Am <= 0 or rho <= 0 or delP <= 0:
        raise ValueError(
            f"solve_mass_balance_only: bad physical constants "
            f"(Am={Am}, rho={rho}, delP={delP}).")
    Lp = slope_g_per_s * 36000.0 / (Am * rho * delP)

    if LOUD:
        print("###################################################################")
        print("Mass-balance-only solver (closed-form, numpy)")
        print(f"  N points   = {dts.size}")
        print(f"  slope      = {slope_g_per_s*1e3:.4f} mg/s")
        print(f"  Lp         = {Lp:.4f} L/m²/hr/bar")
        print(f"  ΔP / Am / ρ = {delP:.3f} bar / {Am} cm² / {rho} g/cm³")

    # ---- Build sim_stru: predicted mV per vial --------------------------
    Jw           = Lp * delP / 36000.0      # cm/s
    flux_g_per_s = Jw * Am * rho            # g/s; equiv to fitted slope

    sim_stru = []
    for i, row in enumerate(data_raw):
        t = np.asarray(row.get('time', []), dtype=float).reshape(-1)
        if t.size == 0:
            sim_stru.append({"time": np.array([]), "mV": np.array([]),
                              "cF": np.array([]), "cH": np.array([]),
                              "cV": np.array([])})
            continue
        if (i + 1) < n_v0:
            # Extra startup vial: emit zeros so the renderer skips the trace.
            mV = np.zeros_like(t)
        else:
            t_start = float(t[0])
            mV = flux_g_per_s * (t - t_start)
        sim_stru.append({
            "time": t.copy(),
            "mV":   mV,
            "cF":   np.zeros_like(t),
            "cH":   np.zeros_like(t),
            "cV":   np.array([0.0]),
        })

    # ---- Residual statistics --------------------------------------------
    pred = flux_g_per_s * dts
    resid = pred - ms
    rmse  = float(np.sqrt(np.mean(resid ** 2)))
    sse   = float(np.sum(resid ** 2))

    fit_stru = {
        "parameters": {
            "Lp":    float(Lp),
            "B":     None,    # not estimated in mass-balance-only mode
            "sigma": 0.0,     # pinned by construction
        },
        "sideparameters":  {},
        "diagnostic_mode": "mass_balance_only",
        "solver_status":   "closed_form",
        "Obj":             float(1e4 * sse / max(dts.size, 1)),
        "obj_m":           float(sse / max(dts.size, 1)),
        "obj_cv":          0.0,
        "obj_cr":          0.0,
        "rmse_g":          rmse,
        "n_points":        int(dts.size),
        "slope_g_per_s":   slope_g_per_s,
    }
    return fit_stru, sim_stru


def _last_valid_value(values, default=np.nan):
    """Return the last non-NaN value from a list-like input."""
    if isinstance(values, (list, tuple, np.ndarray)):
        for item in reversed(values):
            if not (isinstance(item, float) and np.isnan(item)):
                return item
        return default
    if values is None:
        return default
    return values


def _abs_scale(value, fraction, floor):
    """Make a simple measurement error that is never zero."""
    try:
        scale = abs(float(value)) * fraction
    except Exception:
        scale = floor
    return max(scale, floor)


def _invert_monotone_1d(*, fwd_func, y_target, x_lo, x_hi, tol=1e-6, max_iter=200):
    """Invert a scalar conductivity model by bracketing and bisection."""
    import math

    def g(x):
        return float(fwd_func(float(x)) - y_target)

    lo = float(x_lo)
    hi = float(x_hi)
    for fac in [1.0, 2.0, 5.0, 10.0]:
        a = lo
        b = hi * fac
        grid = np.linspace(a, b, 80)
        vals = []
        for x in grid:
            try:
                gx = g(x)
            except Exception:
                continue
            if math.isnan(gx) or math.isinf(gx):
                continue
            vals.append((x, gx))
        for (x0, g0), (x1, g1) in zip(vals[:-1], vals[1:]):
            if g0 == 0:
                return float(x0)
            if g0 * g1 <= 0:
                left, right = x0, x1
                gl, gr = g0, g1
                for _ in range(max_iter):
                    mid = 0.5 * (left + right)
                    gm = g(mid)
                    if abs(gm) <= tol:
                        return float(mid)
                    if (gl < 0 < gm) or (gl > 0 > gm):
                        right, gr = mid, gm
                    else:
                        left, gl = mid, gm
                return float(0.5 * (left + right))
    raise ValueError("Target conductivity value could not be inverted.")


def _conductivity_to_concentration_series(
    *,
    cond_signal,
    temp_K,
    salt_name,
    model="variant_shedlovsky",
    model_params=None,
    output_units="mM",
):
    """Convert conductivity measurements to concentration using conductivity_paper.py.

    Paper equation path:
    - Shedlovsky equivalent conductivity model: \u03bb = f(c)
    - Variant Shedlovsky specific conductivity model: \u03ba = c \u00d7 \u03bb
    - Numerical inversion step: find c such that \u03ba(c) matches the measured conductivity

    The code uses the paper model forward and then inverts it point by point so
    the rest of the diafiltration workflow can continue to treat cF_exp as a
    concentration field.
    """
    # Load the conductivity model from the standalone paper file.
    cp = _load_conductivity_paper()
    # Start from a small built-in table and let the caller override it if needed.
    params = dict(CONDUCTIVITY_SALT_PARAMS_25C.get(str(salt_name), {}))
    if model_params:
        params.update(model_params)
    if not params:
        raise ValueError(
            f"No built-in conductivity parameters for salt '{salt_name}'. "
            "Pass model_params explicitly for this salt."
        )

    # Work with a lowercase model name so the caller can use either case.
    model = str(model).lower()
    cond_signal = np.asarray(cond_signal, dtype=float).reshape(-1)
    conc_out = np.full_like(cond_signal, np.nan, dtype=float)

    if model not in {"variant_shedlovsky", "shedlovsky"}:
        raise ValueError(
            "This refactor currently supports the paper's Shedlovsky / variant Shedlovsky path only."
        )

    epsilon = float(params["epsilon"])
    eta = float(params["eta"])
    a = float(params["a"])
    z_1 = int(params["z_1"])
    z_2 = int(params["z_2"])
    lambda_0_cation = float(params["lambda_0_cation"])
    lambda_0_anion = float(params["lambda_0_anion"])
    lambda_0 = float(params.get("lambda_0", lambda_0_cation + lambda_0_anion))

    def fwd(conc_M):
        # Forward model: concentration -> specific conductivity (mS/cm).
        return float(
            cp.variant_shedlovsky(
                [conc_M],
                temp_K,
                epsilon,
                eta,
                lambda_0,
                a,
                z_1,
                z_2,
                lambda_0_cation,
                lambda_0_anion,
            )[0]
        )

    for i, y in enumerate(cond_signal):
        if np.isnan(y):
            continue
        # Invert the paper conductivity curve to recover concentration from conductivity.
        conc_out[i] = _invert_monotone_1d(
            fwd_func=fwd,
            y_target=float(y / 1000.0),
            x_lo=0.0,
            x_hi=6.0,
        )

    if str(output_units).lower() in {"mm", "mmol/l", "mmolar"}:
        conc_out *= 1000.0
    return conc_out


def convert_experimental_conductivity_to_concentration(
    data_stru,
    *,
    output_units="mM",
    model="auto",
    model_params=None,
):
    """
    Convert retentate conductivity measurements to concentration in one place.

    This is the only data-conversion function that reaches into
    conductivity_paper.py. Everything else should work with the converted
    diafiltration dataset and stay inside the refactored code path.
    """
    # If the data structure is not the expected kind, leave it alone.
    if not isinstance(data_stru, dict) or not data_stru.get("conductivity_cF"):
        return data_stru

    # Read the salt name and temperature from the loaded experimental file.
    cfg = data_stru.get("data_config", {})
    salt_name = cfg.get("namec")
    temp_K = cfg.get("Temp", cfg.get("Temp_K"))
    if temp_K is None:
        raise ValueError("Temperature is required to convert conductivity to concentration.")
    if not salt_name:
        raise ValueError("Salt name is required to convert conductivity to concentration.")

    if model_params is None:
        model_params = {}

    # Default to the paper model unless the caller asks for something else.
    chosen_model = "variant_shedlovsky" if str(model).lower() == "auto" else str(model)
    for row in data_stru.get("data_raw", []):
        signal = row.get("cF_exp")
        if signal is None:
            continue
        # Keep the original conductivity values before replacing them.
        if "cF_exp_conductivity" not in row:
            row["cF_exp_conductivity"] = copy.deepcopy(signal)
        # Convert the conductivity series into concentration using the paper model.
        conc = _conductivity_to_concentration_series(
            cond_signal=signal,
            temp_K=float(temp_K),
            salt_name=str(salt_name),
            model=chosen_model,
            model_params=dict(model_params),
            output_units=output_units,
        )
        if np.asarray(signal).ndim == 0 or not isinstance(signal, (list, tuple, np.ndarray)):
            row["cF_exp"] = float(np.asarray(conc).reshape(-1)[0])
        else:
            row["cF_exp"] = np.asarray(conc, dtype=float).tolist()

    # Mark the file as converted so the rest of the workflow can tell what happened.
    data_stru["conductivity_cF_converted"] = True
    data_stru["conductivity_cF"] = False
    return data_stru


def convert_experimental_permeate_conductivity_to_concentration(
    data_stru, *, output_units="mM", model="auto", model_params=None,
):
    """Mirror of :func:`convert_experimental_conductivity_to_concentration`
    for the **permeate** probe.  Reads ``cV_perm_cond`` (uS/cm) per vial,
    inverts via the same Shedlovsky path, stores ``cV_perm_exp`` (mM) per
    vial alongside the raw conductivity field.

    Idempotent — running a second time is a no-op.  Bad / non-physical
    samples are set to NaN so the objective can safely skip them.

    GUARDED:  this is only meaningful for DATA3/NF270 sheets, which are the
    only ones that carry ``cV_perm_cond``.  DATA1/DATA2 .mat files do not
    have this field, so calling this on them is a silent no-op.
    """
    if not isinstance(data_stru, dict):
        return data_stru
    if data_stru.get("conductivity_cV_perm_converted"):
        return data_stru
    cfg = data_stru.get("data_config", {})
    salt_name = cfg.get("namec")
    temp_K = cfg.get("Temp", cfg.get("Temp_K"))
    if temp_K is None or not salt_name:
        return data_stru   # no salt or no temperature → can't invert; leave alone
    chosen_model = "variant_shedlovsky" if str(model).lower() == "auto" else str(model)
    if model_params is None:
        model_params = {}

    any_converted = False
    for row in data_stru.get("data_raw", []):
        signal = row.get("cV_perm_cond")
        if signal is None:
            continue
        signal_arr = np.asarray(signal, dtype=float).reshape(-1)
        if signal_arr.size == 0:
            continue
        if "cV_perm_cond_uS_per_cm" not in row:
            row["cV_perm_cond_uS_per_cm"] = copy.deepcopy(signal)
        # Loader stores conductivity in uS/cm; paper model wants mS/cm.
        signal_mS = signal_arr / 1000.0
        try:
            conc = _conductivity_to_concentration_series(
                cond_signal=signal_mS, temp_K=float(temp_K),
                salt_name=str(salt_name), model=chosen_model,
                model_params=dict(model_params),
                output_units=output_units,
            )
            conc = np.asarray(conc, dtype=float)
            # Filter unreliable samples: negative concentrations, or NaN
            # from a failed bisection inversion.  Downstream code already
            # masks NaN cells, so this is a clean no-op for those rows.
            conc[~np.isfinite(conc)] = np.nan
            conc[conc < 0] = np.nan
            row["cV_perm_exp"] = conc.tolist()
            any_converted = True
        except Exception:
            row["cV_perm_exp"] = [float("nan")] * int(signal_arr.size)
    if any_converted:
        data_stru["conductivity_cV_perm_converted"] = True
    return data_stru


def _normalize_conductivity_measurements(data_stru, *, output_units="mM", model="auto", model_params=None):
    """Convert retentate conductivity measurements to concentration when needed.

    If the file says conductivity is present, we first pre-compensate any
    temperature-tagged conductivity traces to a 25 C-equivalent value, then
    replace the fitting field with concentration so the rest of the
    diafiltration equations stay unchanged.

    This keeps the data-preprocessing step aligned with the paper workflow:
    - preserve the original conductivity values for traceability
    - compensate conductivity to 25 C when per-row temperatures are available
    - convert cF_exp using the paper conductivity equation
    - feed the converted concentration into the existing first-principles model

    Permeate-probe conversion (cV_perm_cond -> cV_perm_exp) is ALSO run here
    when the module-level toggle ``NF270_USE_PERMEATE_PROBE`` is truthy and
    the data carries a permeate-conductivity field.  Guarded so DATA1/DATA2
    sheets (no cV_perm_cond) silently no-op.
    """
    if not isinstance(data_stru, dict):
        return data_stru

    if not data_stru.get("conductivity_cF"):
        return data_stru

    data_stru = _apply_ec25_compensation_to_data_stru(data_stru)

    data_stru = convert_experimental_conductivity_to_concentration(
        data_stru,
        output_units=output_units,
        model=model,
        model_params=model_params,
    )

    # Optional permeate-probe channel — only when toggle is on AND the
    # data carries cV_perm_cond (DATA3 / NF270 sheets).
    if NF270_USE_PERMEATE_PROBE:
        data_stru = convert_experimental_permeate_conductivity_to_concentration(
            data_stru,
            output_units=output_units,
            model=model,
            model_params=model_params,
        )
    return data_stru


def _theta_components(model):
    """Pick the model components that should be estimated."""
    theta_components = []

    if hasattr(model, "Lp"):
        theta_components.append(model.Lp)

    if hasattr(model, "B"):
        try:
            theta_components.extend([comp for _, comp in model.B.items()])
        except Exception:
            theta_components.append(model.B)

    for name in ("beta_0", "beta_1", "beta_2", "beta_3", "sigma", "S0", "S"):
        if hasattr(model, name):
            theta_components.append(getattr(model, name))

    return theta_components


def _label_parmest_model(model, data_stru, mode="DATA", workflow_family="DATA1"):
    """Attach the labels ParmEst needs."""
    data_stru = _normalize_conductivity_measurements(data_stru)
    # Keep parmest's measurement_error consistent with obj_rule + calc_FIM.
    _pe_cf_frac, _pe_cp_frac = _nf270_conc_scales(workflow_family)
    model.experiment_outputs = Suffix(direction=Suffix.LOCAL)
    model.measurement_error = Suffix(direction=Suffix.LOCAL)
    model.unknown_parameters = Suffix(direction=Suffix.LOCAL)
    model.experiment_inputs = Suffix(direction=Suffix.LOCAL)

    theta_components = _theta_components(model)
    model.unknown_parameters.update(
        (comp, ComponentUID(comp)) for comp in theta_components
    )

    final_t = model.tau.last()
    output_items = []
    n_v0 = int(data_stru.get("data_config", {}).get("n_v0", 1))

    for i in model.n_vial:
        if i < n_v0:
            continue
        row = data_stru["data_raw"][i - 1]
        cF_obs = _last_valid_value(row.get("cF_exp"))
        cV_obs = _last_valid_value(row.get("cV_avg"))
        mass_obs = _last_valid_value(row.get("mass"))

        output_items.append((model.cF[i, final_t], cF_obs, _abs_scale(cF_obs, _pe_cf_frac, 0.003)))
        output_items.append((model.cV[i, final_t], cV_obs, _abs_scale(cV_obs, _pe_cp_frac, 0.03)))

        if hasattr(model, "mF"):
            output_items.append((model.mF[i, final_t], mass_obs, _abs_scale(mass_obs, 0.01, 0.01)))

    final_tube_meas = _last_valid_value(data_stru.get("data_config", {}).get("cF_final_meas"))
    if final_tube_meas is not None and np.isfinite(final_tube_meas):
        model.cF_terminal_obs = Expression(expr=model.cF[model.n_vial.last(), final_t])
        output_items.append(
            (model.cF_terminal_obs, final_tube_meas, _abs_scale(final_tube_meas, _pe_cf_frac, 0.003))
        )

    for expr, obs, err in output_items:
        model.experiment_outputs[expr] = obs
        model.measurement_error[expr] = err

    return model


class DiafiltrationParmestExperiment(ParmestExperiment):
    """Small ParmEst wrapper around the diafiltration model."""

    def __init__(self, data_stru, mode="DATA", theta=None, B_form="single", nfe=300, workflow_family="DATA1"):
        super().__init__()
        self.data_stru = data_stru
        self.mode = mode
        self.theta = theta
        self.B_form = B_form
        self.nfe = nfe
        self.workflow_family = workflow_family

    def get_labeled_model(self):
        """Build, discretize, and label the model for ParmEst."""
        model = model_construct_inter(
            self.data_stru,
            self.mode,
            theta=self.theta,
            sim_opt=False,
            B_form=self.B_form,
            workflow_family=self.workflow_family,
        )
        TransformationFactory("dae.finite_difference").apply_to(
            model, nfe=self.nfe, scheme="BACKWARD"
        )
        return _label_parmest_model(model, self.data_stru, mode=self.mode, workflow_family=self.workflow_family)


def build_parmest_experiments(data_structures, mode="DATA", theta=None, B_form="single", nfe=300, workflow_family="DATA1"):
    """Turn one data set or many data sets into a ParmEst experiment list."""
    if not isinstance(data_structures, list):
        data_structures = [data_structures]
    return [
        DiafiltrationParmestExperiment(data_stru, mode=mode, theta=theta, B_form=B_form, nfe=nfe, workflow_family=workflow_family)
        for data_stru in data_structures
    ]


def _normalize_experimental_source(source, *, campaign=None, variant="base", data_root=None, selector=None):
    """Turn one experimental source into a loaded run instance and data dict."""
    if isinstance(source, ExperimentalRun):
        run = source
        return run, run.load_data()

    if isinstance(source, dict) and {"run_id", "root", "data_file"}.issubset(source.keys()):
        run = ExperimentalRun(
            run_id=str(source["run_id"]),
            root=Path(source["root"]),
            data_file=str(source["data_file"]),
            variant=str(source.get("variant", variant)),
            subdir=str(source.get("subdir", "")),
            fit_file=str(source.get("fit_file", "fit_stru.mat")),
        )
        return run, run.load_data()

    path = Path(source).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Experimental source not found: {path}")
    if path.suffix.lower() == ".csv":
        skiprows = [1] if any(token in path.name.lower() for token in ("pressure", "calibration")) else None
        return None, _load_csv_artifact(path, skiprows=skiprows)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return None, _load_legacy_data_stru_from_excel(path, selector=selector)
    raw = loadmat(str(path))
    data_stru = raw.get("data_stru")
    if data_stru is None:
        raise ValueError(f"The file {path} does not contain a 'data_stru' entry.")
    return None, _normalize_conductivity_measurements(data_stru)


def _load_csv_artifact(path: Path, skiprows=None) -> dict:
    """Load a CSV artifact into the pipeline as data-only StageResults input."""
    df = pd.read_csv(path, skiprows=skiprows)
    return {
        "artifact_kind": "csv",
        "artifact_path": str(path),
        "data": df,
        "data_file": [str(path)],
        "is_batch": False,
    }


def _parse_excel_metadata_and_table(sheet_df: pd.DataFrame) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Split the workbook into metadata, time-series, and calibration tables."""
    metadata = {}
    if sheet_df.shape[0] >= 1 and sheet_df.shape[1] >= 11:
        for row_idx in range(min(sheet_df.shape[0], 20)):
            key = sheet_df.iat[row_idx, 9] if 9 < sheet_df.shape[1] else None
            value = sheet_df.iat[row_idx, 10] if 10 < sheet_df.shape[1] else None
            if isinstance(key, str) and key.strip():
                metadata[key.strip()] = value
    time_cols = [0, 1, 2, 3, 4, 5, 6, 7]
    ts = sheet_df.iloc[3:, time_cols].copy()
    ts.columns = ["Time (s)", "Mass (g)", "Pressure (psi)", "Retentate Temp", "Retentate Cond @ Temp (uS/cm)", "Permeate Temp", "Permeate Cond @ Temp (uS/cm)", "Vial Swap"]
    ts = ts.dropna(how="all")
    calib = sheet_df.iloc[:, 23:27].copy() if sheet_df.shape[1] >= 27 else pd.DataFrame()
    return metadata, ts, calib


def _parse_excel_salt_name(sheet_name: str, metadata: dict, note_text: str = "") -> str:
    """Infer the salt identity from the sheet name or notes."""
    candidates = ["NaCl", "CaCl2", "LaCl3", "KCl"]
    blob = " ".join([str(sheet_name or ""), note_text, " ".join(str(v) for v in metadata.values())])
    for cand in candidates:
        if cand.lower() in blob.lower():
            return cand
    return "NaCl"


def _detect_nf270_mode(note_text: str, mass_all_cum=None) -> tuple[str, str]:
    """Detect Lag vs Overflow mode for an NF270 .xlsx sheet.

    Returns
    -------
    (mode, source) : tuple of str
        mode   — 'Lag' or 'Overflow' (the value going into data_stru["mode"])
        source — 'notes_text', 'mass_trajectory', or 'default' (audit trail)

    Detection order (first match wins):
      1. Notes text explicitly mentions 'overflow' or 'lag'
      2. Notes text mentions 'wash with' (a Lag operational signal) or
         'continuous diafiltrate' (an Overflow signal)
      3. Mass-trajectory heuristic: if the cell mass dips significantly
         and then RECOVERS toward the initial value, it's Lag (the wash
         restored M_F). If mass stays approximately flat, it's Overflow.
      4. Default to 'Lag' if nothing matches.

    Same architectural contract as DATA2: the user passes a known mode at
    call-time; this helper just supplies a sensible default detected from
    the workbook so the user doesn't have to spell it out per sheet.
    """
    text = str(note_text or "").lower()
    # 1. Explicit mode keyword
    if "overflow" in text:
        return ("Overflow", "notes_text")
    if "lag" in text and "lag" not in {"lagging", "lagged"}:  # avoid false-positive substrings
        return ("Lag", "notes_text")
    # 2. Operational signal
    if "continuous diafiltrate" in text or "continuous diafilt" in text:
        return ("Overflow", "notes_text")
    if "wash with" in text or "batch wash" in text:
        return ("Lag", "notes_text")
    # 3. Mass trajectory: Lag = significant dip + recovery; Overflow = flat
    if mass_all_cum is not None:
        try:
            mass = np.asarray(mass_all_cum, dtype=float).reshape(-1)
            finite = mass[np.isfinite(mass)]
            if finite.size > 10:
                m_start = float(finite[0])
                m_end = float(finite[-1])
                m_max = float(np.max(finite))
                m_min = float(np.min(finite))
                span = max(abs(m_max - m_min), 1e-9)
                start_to_end = abs(m_end - m_start)
                # Heuristic threshold: if start ≈ end (within 20% of span)
                # AND there was significant variation in between (>3× the
                # start-to-end delta), the cell mass dipped and recovered →
                # Lag. Otherwise default Lag (safe).
                if start_to_end < 0.2 * span and span > 3 * start_to_end:
                    return ("Lag", "mass_trajectory")
        except Exception:
            pass
    # 4. Default
    return ("Lag", "default")


def _segment_indices_by_swap(swap_flags: np.ndarray) -> list[tuple[int, int]]:
    """Split a time-series into vial segments using a 0/1 swap flag column."""
    swap_flags = np.asarray(swap_flags).reshape(-1)
    swap_idx = np.where(swap_flags == 1)[0]
    segments: list[tuple[int, int]] = []
    start = 0
    for idx in swap_idx:
        end = idx + 1
        segments.append((start, end))
        start = end
    if start < len(swap_flags):
        segments.append((start, len(swap_flags)))
    return segments


# Sidecar table of CURATED holdup-fill boundaries (DATA2-style), keyed by
# (workbook filename, sheet name). The experimenter records the true end of the
# tube-fill / holdup period in seconds (sheet time axis) in the `holdup_end_s`
# column; the loader honors it instead of the mass-rise auto-detect. Blank →
# fall back to auto-detect. See _make_holdup_boundary_table.py to (re)generate
# the pre-filled table + the per-sheet diagnostic plots used to curate it.
_NF270_HOLDUP_TABLE_PATH = Path(__file__).resolve().parent / "nf270_holdup_boundaries.csv"
_NF270_HOLDUP_TABLE_CACHE = None  # {(workbook, sheet): holdup_end_s}


def _nf270_curated_holdup_end_s(workbook_name, sheet_name):
    """Return the curated holdup-end time [s] for a sheet, or None.

    Reads the sidecar CSV once and caches it. A row contributes an override only
    when its `holdup_end_s` cell is a finite number; otherwise the loader falls
    back to the mass-rise auto-detect. Any read/parse problem degrades silently
    to None so the loader never breaks on a missing or malformed table."""
    global _NF270_HOLDUP_TABLE_CACHE
    if _NF270_HOLDUP_TABLE_CACHE is None:
        table = {}
        try:
            if _NF270_HOLDUP_TABLE_PATH.exists():
                df = pd.read_csv(_NF270_HOLDUP_TABLE_PATH)
                for _, row in df.iterrows():
                    wb = str(row.get("workbook", "")).strip()
                    sh = str(row.get("sheet", "")).strip()
                    val = pd.to_numeric(row.get("holdup_end_s"), errors="coerce")
                    if wb and sh and np.isfinite(val):
                        table[(wb, sh)] = float(val)
        except Exception:
            table = {}
        _NF270_HOLDUP_TABLE_CACHE = table
    return _NF270_HOLDUP_TABLE_CACHE.get((str(workbook_name).strip(), str(sheet_name).strip()))


_CATION_MW_G_PER_MOL = {
    "NaCl":  22.99,    # Na
    "KCl":   39.10,    # K
    "CaCl2": 40.08,    # Ca
    "LaCl3": 138.91,   # La
}


# ---------------------------------------------------------------------------
# ICP-OES calibration audit (read the calibration table embedded in each
# workbook and cross-check the precomputed mg/L values).
# ---------------------------------------------------------------------------
# Each NF270 workbook sheet carries an ICP calibration block at cols 23-26:
#     col 23/25: Concentration (mg/L)        ← calibration standards
#     col 24/26: Intensity (cps)             ← raw ICP reading at each standard
# The ICP machine fits a line through these and uses it to convert each
# vial's intensity (col 18) to mg/L (col 19). The cross-check below
# re-fits the line ourselves and verifies the stored mg/L matches what
# the calibration curve would predict from the same intensity.
#
# This is a pure DATA AUDIT — does NOT change the cV_avg values; it just
# adds a residual field per ICP row so we can flag any disagreement > 5%.

def _parse_icp_calibration(sheet_df, salt_slot=1):
    """Read the ICP calibration table and fit a linear curve.

    Salt 1 calibration lives at cols 23-24, Salt 2 at cols 25-26.
    Header rows are 1-2; data starts at row 3 and runs until NaN.

    Returns a dict with the fitted curve + the raw points, or None if the
    calibration table is missing or too short to fit.
    """
    base_col = 23 if salt_slot == 1 else 25
    if base_col + 1 >= sheet_df.shape[1]:
        return None

    mg_L_list = []
    cps_list = []
    for r in range(3, min(sheet_df.shape[0], 60)):
        c = sheet_df.iat[r, base_col]
        i = sheet_df.iat[r, base_col + 1]
        c_f = pd.to_numeric(c, errors="coerce")
        i_f = pd.to_numeric(i, errors="coerce")
        if pd.notna(c_f) and pd.notna(i_f):
            mg_L_list.append(float(c_f))
            cps_list.append(float(i_f))

    if len(mg_L_list) < 2:
        return None

    cps_arr = np.asarray(cps_list, dtype=float)
    mg_L_arr = np.asarray(mg_L_list, dtype=float)

    # Fit  mg_L = slope * cps + intercept  (intensity is the independent
    # variable because we'll invert "stored mg/L vs stored cps" downstream).
    slope, intercept = np.polyfit(cps_arr, mg_L_arr, 1)
    mg_L_pred = slope * cps_arr + intercept
    ss_res = float(np.sum((mg_L_arr - mg_L_pred) ** 2))
    ss_tot = float(np.sum((mg_L_arr - np.mean(mg_L_arr)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return {
        "slope_mg_per_L_per_cps": float(slope),
        "intercept_mg_per_L": float(intercept),
        "r_squared": float(r_squared),
        "n_points": int(len(mg_L_arr)),
        "mg_L_range": [float(np.min(mg_L_arr)), float(np.max(mg_L_arr))],
        "cps_range": [float(np.min(cps_arr)), float(np.max(cps_arr))],
        "points_mg_L": [float(x) for x in mg_L_arr.tolist()],
        "points_cps":  [float(x) for x in cps_arr.tolist()],
    }


def _icp_calibration_residual(stored_mg_L, intensity_cps, calibration):
    """Cross-check one stored ICP mg/L value against the fitted calibration.

    Returns a small dict with the predicted mg/L from the calibration, the
    absolute residual, the relative residual (%), and a 5%-tolerance flag.
    Returns None if the calibration or the intensity is unusable.
    """
    if calibration is None:
        return None
    if not (np.isfinite(intensity_cps) and intensity_cps > 0):
        return None
    if not np.isfinite(stored_mg_L):
        return None

    slope = calibration["slope_mg_per_L_per_cps"]
    intercept = calibration["intercept_mg_per_L"]
    predicted = slope * intensity_cps + intercept

    abs_res = float(stored_mg_L - predicted)
    rel_pct = float(abs_res / predicted * 100.0) if abs(predicted) > 1e-9 else float("nan")

    return {
        "predicted_mg_L": float(predicted),
        "stored_mg_L": float(stored_mg_L),
        "absolute_residual_mg_L": float(abs_res),
        "relative_residual_pct": float(rel_pct),
        "within_5pct_tolerance": bool(abs(rel_pct) < 5.0) if np.isfinite(rel_pct) else False,
    }

# ---------------------------------------------------------------------------
# Salt-aware mass-residual uncertainty (used by legacy diagnostics)
# ---------------------------------------------------------------------------
# The WSSE objective now keeps the DATA2-style 0.01 g default in the loader.
# This helper remains available for older notebooks / diagnostics that want
# to inspect the balance-noise + model-misfit decomposition explicitly.
MASS_BALANCE_SIGMA_G = 0.01         # analytical balance precision, g/reading
MASS_MODEL_FIT_TOLERANCE = 0.05     # expected model misfit, fraction of |signal|


def _data3_mass_scale_g(M_F0, final_solution_weight):
    """Per-sheet mass uncertainty for legacy diagnostics.

    Combines balance precision (instrument noise floor, constant across salts)
    with a model-misfit budget proportional to the per-sheet mass-loss signal
    |Final Solution Weight - M_F0|, in quadrature.  The loader does not use
    this value as the default WSSE weight anymore; it is retained for
    compatibility with old plots and reports.

    Returns
    -------
    mass_scale_g : float
        sigma on a single mass measurement, in grams.
    components : dict
        Breakdown of the two contributions, useful for plot annotation /
        diagnostics:
            {"signal_g": float, "sigma_balance_g": float,
             "sigma_model_g": float, "sigma_total_g": float,
             "model_fraction": float}
    """
    try:
        M_F0_f = float(M_F0)
    except Exception:
        M_F0_f = 0.0
    try:
        final_f = float(final_solution_weight)
    except Exception:
        final_f = 0.0
    # Use the observed per-sheet mass loss as the signal term for the
    # salt-specific model-misfit budget.
    signal = abs(final_f - M_F0_f) if final_f > 0 else 0.0
    sigma_balance = float(MASS_BALANCE_SIGMA_G)
    sigma_model = float(MASS_MODEL_FIT_TOLERANCE) * signal
    sigma_total = float(np.sqrt(sigma_balance ** 2 + sigma_model ** 2))
    return sigma_total, {
        "signal_g":         signal,
        "sigma_balance_g":  sigma_balance,
        "sigma_model_g":    sigma_model,
        "sigma_total_g":    sigma_total,
        "model_fraction":   float(MASS_MODEL_FIT_TOLERANCE),
    }


def _load_legacy_data_stru_from_excel(path: Path, selector: object = None) -> dict:
    """Load a DATA3 / NF270 Excel workbook into the legacy data_stru dict.

    The output is structurally identical to the DATA2 .mat data_stru that
    utility.py and the legacy parameter estimator already understand:
        - per-vial mass arrays are RESET so each vial starts at ~0g
        - per-vial cV_avg is the SCALAR ICP-OES concentration (mM), read
          from the sidebar columns; falls back to NaN when ICP isn't
          available
        - the holdup-fill period at the start of the run is split out as
          a leading "extra" vial (n_v0=2, n_extra=1) when a lag is
          detected; the model integrates Jw over that vial but the
          objective skips its residuals
        - Lp0 is seeded from vial-1 slope (data-driven), not hard-coded
        - data_config carries every field DATA2's data_config has,
          including theta0, nc, nr
    """
    xl = pd.ExcelFile(path)
    sheet_name = selector if isinstance(selector, str) and selector in xl.sheet_names else xl.sheet_names[0]
    sheet_df = pd.read_excel(path, sheet_name=sheet_name, header=None)
    metadata, ts, _calib = _parse_excel_metadata_and_table(sheet_df)
    note_text = ""
    if sheet_df.shape[0] > 0 and sheet_df.shape[1] > 4:
        note_text = str(sheet_df.iat[0, 4]) if pd.notna(sheet_df.iat[0, 4]) else ""
    salt_name = _parse_excel_salt_name(sheet_name, metadata, note_text=note_text)

    def _to_float_series(series):
        return np.asarray(pd.to_numeric(series, errors="coerce"), dtype=float)

    time_all = _to_float_series(ts["Time (s)"])
    mass_all_cum = _to_float_series(ts["Mass (g)"])     # cumulative balance reading
    pressure_all = _to_float_series(ts["Pressure (psi)"])
    ret_temp_all = _to_float_series(ts["Retentate Temp"])
    ret_cond_all = _to_float_series(ts["Retentate Cond @ Temp (uS/cm)"])
    perm_temp_all = _to_float_series(ts["Permeate Temp"])
    perm_cond_all = _to_float_series(ts["Permeate Cond @ Temp (uS/cm)"])
    swap_all = _to_float_series(ts["Vial Swap"]).astype(int)

    # Preserve the measured temperature trajectories for diagnostics, but
    # normalize conductivity at load time so all downstream code sees 25 C
    # equivalent values.
    measured_ret_temp_all = ret_temp_all.copy()
    measured_perm_temp_all = perm_temp_all.copy()
    ret_cond_all = _ec25_compensate(ret_cond_all, ret_temp_all, salt_name)
    perm_cond_all = _ec25_compensate(perm_cond_all, perm_temp_all, salt_name)

    segments = _segment_indices_by_swap(swap_all)
    if not segments:
        raise ValueError(f"No vial segments could be parsed from Excel sheet '{sheet_name}'.")

    temp_k = 298.15  # conductivity compensated to 25 C at load time
    delp_psi = float(np.nanmean(pressure_all)) if np.isfinite(np.nanmean(pressure_all)) else 0.0
    delp = delp_psi * 0.0689475729
    rpm_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*RPM\b", note_text, flags=re.IGNORECASE)
    rpm_note = float(rpm_match.group(1)) if rpm_match else 350.0
    feed_match = re.search(r"Feed[\s\-:]+([0-9]+(?:\.[0-9]+)?)\s*mM", note_text, flags=re.IGNORECASE)
    diafiltrate_match = re.search(r"Diafiltrate[\s\-:]+([0-9]+(?:\.[0-9]+)?)\s*mM", note_text, flags=re.IGNORECASE)
    note_mM_tokens = []
    for token in re.findall(r"([0-9]+(?:\.[0-9]+)?)\s*mM", note_text, flags=re.IGNORECASE):
        try:
            note_mM_tokens.append(float(token))
        except Exception:
            continue
    first_note_mM = note_mM_tokens[0] if note_mM_tokens else np.nan
    def _meta_float(*keys, default):
        for key in keys:
            if key in metadata:
                val = pd.to_numeric(metadata.get(key), errors="coerce")
                if pd.notna(val):
                    return float(val)
        return float(default)

    # Read the embedded ICP calibration table once (per sheet). Used below
    # to cross-check that the precomputed mg/L matches what the calibration
    # curve would predict from the same intensity. Adds an audit field per
    # row; does NOT change cV_avg.
    icp_calibration = _parse_icp_calibration(sheet_df, salt_slot=1)

    def _parse_sidebar_icp_rows():
        """Collect the sidebar ICP rows so cross-checks stay available downstream."""
        rows = {}
        per_vial = {}
        if sheet_df.shape[1] < 20:
            return rows, per_vial

        for r in range(min(sheet_df.shape[0], 80)):
            label = sheet_df.iat[r, 14] if 14 < sheet_df.shape[1] else None
            if not isinstance(label, str):
                continue
            label_norm = re.sub(r"\s+", " ", label).strip()
            try:
                samp_vol = float(pd.to_numeric(sheet_df.iat[r, 16], errors="coerce"))
                acid_vol = float(pd.to_numeric(sheet_df.iat[r, 17], errors="coerce"))
                # NEW: also read col 18 (raw intensity in cps) so we can
                # cross-check against the calibration table downstream.
                intensity_cps = float(pd.to_numeric(sheet_df.iat[r, 18], errors="coerce")) \
                    if 18 < sheet_df.shape[1] else float("nan")
                icp_mg_L = float(pd.to_numeric(sheet_df.iat[r, 19], errors="coerce"))
            except Exception:
                continue
            if not (np.isfinite(samp_vol) and np.isfinite(acid_vol) and np.isfinite(icp_mg_L)):
                continue
            if samp_vol <= 0:
                continue

            cation_mw = _CATION_MW_G_PER_MOL.get(salt_name, 22.99)
            if cation_mw <= 0:
                continue

            dilution = (samp_vol + acid_vol) / samp_vol
            cV_mM = icp_mg_L * dilution / cation_mw

            # Cross-check the stored mg/L against the in-workbook calibration
            # curve. residual is a dict with predicted/stored/abs/relative
            # fields, or None if the calibration table is missing or the
            # intensity is bad.
            calib_residual = _icp_calibration_residual(
                icp_mg_L, intensity_cps, icp_calibration
            )

            record = {
                "label": label_norm,
                "sample_volume_mL": float(samp_vol),
                "nitric_acid_volume_mL": float(acid_vol),
                "icp_intensity_cps": float(intensity_cps) if np.isfinite(intensity_cps) else None,
                "icp_cation_mg_L": float(icp_mg_L),
                "cation_mw_g_per_mol": float(cation_mw),
                "dilution_factor": float(dilution),
                "cV_avg_mM": float(cV_mM),
                "icp_calibration_residual": calib_residual,
            }
            rows[label_norm] = record
            vial_match = re.match(r"Vial\s+(\d+)\s*$", label_norm, flags=re.IGNORECASE)
            if vial_match:
                per_vial[int(vial_match.group(1))] = float(cV_mM)

        return rows, per_vial

    M_F0 = _meta_float("Initial Solution Weight (g)", "Initial Solution Weight:",
                        "Initial weight of solution (g)", "Initial weight (g)",
                        default=10.0)
    final_solution_weight = _meta_float("Final Solution Weight (g)",
                                         " Final Solution Weight:",
                                         "Final Solution Weight:",
                                         "Final weight of solution (g)",
                                         "Final weight (g)", default=0.0)
    n_salts = int(_meta_float("Number of Salts:", "Number of Salts (#)", default=1))

    # ------ Read per-vial ICP-OES from sidebar columns 14-22 ------------
    # Convention from NF270_MC2.xlsx:
    #   col 15: row label ("Vial 1".."Vial N", "Feed", "Diafiltrate", ...)
    #   col 17: sample volume (mL) used for ICP
    #   col 18: nitric acid volume (mL) used to dilute the sample
    #   col 20: ICP Salt 1 reading (mg/L of cation, on the diluted sample)
    # cV_avg for any sidebar row in mM = col20 × (col17+col18)/col17 / cation_MW
    sidebar_icp_rows, icp_per_vial = _parse_sidebar_icp_rows()

    def _sidebar_mM(label: str):
        rec = sidebar_icp_rows.get(label)
        if not rec:
            return np.nan
        value = rec.get("cV_avg_mM", np.nan)
        return float(value) if np.isfinite(value) else np.nan

    feed_icp_mM = _sidebar_mM("Feed")
    diafiltrate_icp_mM = _sidebar_mM("Diafiltrate")
    retentate_icp_mM = _sidebar_mM("Retentate")
    final_tube_icp_mM = _sidebar_mM("Final Tube")
    feed_note_mM = float(feed_match.group(1)) if feed_match else np.nan
    diafiltrate_note_mM = float(diafiltrate_match.group(1)) if diafiltrate_match else np.nan
    c_f0_value = feed_icp_mM if np.isfinite(feed_icp_mM) else (
        feed_note_mM if np.isfinite(feed_note_mM) else first_note_mM
    )
    c_d_value = diafiltrate_icp_mM if np.isfinite(diafiltrate_icp_mM) else (
        diafiltrate_note_mM if np.isfinite(diafiltrate_note_mM) else (
            feed_note_mM if np.isfinite(feed_note_mM) else first_note_mM
        )
    )

    # ------ Build per-vial data_raw with MASS RESET PER VIAL ------------
    data_raw = []
    for idx, (a, b) in enumerate(segments, start=1):
        seg_mass_cum = mass_all_cum[a:b]
        # Per-vial mass: subtract the cumulative reading at the vial start
        # so each vial's mV starts at ~0 (DATA2 .mat convention).
        # Be robust to a NaN at the very first sample by using the first
        # finite cumulative reading as the offset.
        if seg_mass_cum.size > 0:
            finite = np.isfinite(seg_mass_cum)
            if finite.any():
                offset = float(seg_mass_cum[finite][0])
                seg_mass = seg_mass_cum - offset
            else:
                seg_mass = seg_mass_cum
        else:
            seg_mass = seg_mass_cum
        cV_scalar = icp_per_vial.get(idx, np.nan)
        # Permeate placement (DATA2 convention): the single ICP measurement is a
        # per-vial scalar with no recorded sample time, so anchor it at vial
        # close. The fit objective then compares it against the model permeate
        # concentration at tau.last() (see the cV residual branch in
        # model_construct_inter). All three streams (mass, retentate, permeate)
        # stay on ONE shared time axis — we deliberately do NOT displace the
        # permeate to an interior membrane-event index (the retired V_tube
        # tube-transit model), which would desynchronise it from mass/retentate.
        vial_time = time_all[a:b]
        cV_avg_arr = np.full(len(vial_time), np.nan, dtype=float)
        if np.isfinite(cV_scalar) and len(vial_time) > 0:
            cV_avg_arr[-1] = float(cV_scalar)
        row = {
            "number":          idx,
            "time":            vial_time,
            "mass":            seg_mass,
            "cF_exp":          ret_cond_all[a:b],
            "cV_avg":          cV_avg_arr,               # NaN-padded; ICP value at corrected index
            "cV_perm_cond":    perm_cond_all[a:b],       # per-sample permeate conductivity (uS/cm)
            "pressure":        pressure_all[a:b],
            "retentate_temp":  ret_temp_all[a:b],
            "permeate_temp":   perm_temp_all[a:b],
            "vial_swap":       swap_all[a:b],
        }
        data_raw.append(row)

    # Keep the DATA2-style 0.01 g mass weight for parity.  The more
    # elaborate quadrature helper remains available for older diagnostics,
    # but the loader itself no longer uses it as the objective weight.
    _mass_scale_g_val = 0.01
    _mass_scale_components = {}

    # Detect the mode (Lag vs Overflow) from the Notes text + mass trajectory.
    # Defaults to Lag if uncertain. Records the source ("notes_text",
    # "mass_trajectory", or "default") in data_stru["mode_source"] so it's
    # auditable per sheet.
    detected_mode, mode_source = _detect_nf270_mode(note_text, mass_all_cum)

    data_stru = {
        "dataset": Path(path).stem,
        "filename": str(metadata.get("Experiment Name:", Path(path).name)),
        "mode": detected_mode,
        "mode_source": mode_source,
        "continuous_cF": True,
        "conductivity_cF": True,
        "conductivity_cF_converted": False,
        "conductivity_temp_compensated": True,
        "conductivity_temp_compensation": {
            "salt_name": salt_name,
            "alpha_per_C": float(CONDUCTIVITY_TEMP_COEFF_PER_C.get(
                salt_name, CONDUCTIVITY_TEMP_COEFF_DEFAULT
            )),
            "reference_C": 25.0,
            "method": "sigma_25 = sigma_T * (1 + alpha * (25 - T_celsius))",
        },
        "permeate_time_corrected": False,
        "permeate_time_correction": {
            "method": "vial_close",
            "note": (
                "DATA2 convention: the single permeate ICP value is anchored at "
                "vial close and compared against the model cV at tau.last(). All "
                "streams share one time axis; no tube-transit (V_tube) displacement."
            ),
        },
        "measured_ret_temp_C": measured_ret_temp_all,
        "measured_perm_temp_C": measured_perm_temp_all,
        "data_config": {
            "n":          len(data_raw),
            "n_v0":       1,
            "n_extra":    0,
            "n_h":        0,
            "n_A":        0,
            "nr":         0,
            "nc":         n_salts,
            "delP":       delp,
            "Temp":       temp_k,
            "rpm":        rpm_note,
            "stir_rpm":   rpm_note,
            "rpm_source": "note_text" if rpm_match else "default_350RPM",
            "Am":         4.1,           # cm²
            "rho":        1.0,           # g/cm³
            "M_F0":       M_F0,
            # M_O = mass overflow / removed: M_F0 minus final solution weight, signed.
            # DATA2 stores it as negative when feed lost mass (lag mode).
            "M_O":        float(final_solution_weight - M_F0) if final_solution_weight > 0 else 0.0,
            # WSSE mass-residual uncertainty for this sheet.  DATA3 keeps
            # the DATA2-style 0.01 g weight unless a caller overrides it.
            "mass_scale_g":          _mass_scale_g_val,
            "mass_scale_source":     "legacy balance-only (DATA2 parity)",
            "mass_scale_components": _mass_scale_components,
            "C_D":        c_d_value if np.isfinite(c_d_value) else 0.0,
            "C_D_source": "sidebar_icp:Diafiltrate" if np.isfinite(diafiltrate_icp_mM) else (
                "note_text:Diafiltrate" if np.isfinite(diafiltrate_note_mM) else (
                    "note_text:first_mM" if np.isfinite(first_note_mM) else "unavailable"
                )
            ),
            "C_D_units":  "mM",
            "C_F0":       c_f0_value if np.isfinite(c_f0_value) else 0.0,
            "C_F0_source": "sidebar_icp:Feed" if np.isfinite(feed_icp_mM) else (
                "note_text:Feed" if np.isfinite(feed_note_mM) else (
                    "note_text:first_mM" if np.isfinite(first_note_mM) else "first_cF_exp_fallback"
                )
            ),
            "C_F0_units": "mM",
            "namec":      salt_name,
            "ni":         1,
            "Lp0":        5.0,
            "B0":         0.5,
            "sigma0":     0.9,
            "theta0":     np.array([5.0, 0.5, 0.9], dtype=float),
            "icp_sidebar_rows": sidebar_icp_rows,
            "icp_calibration_curve": icp_calibration,  # fitted slope/intercept/R² + raw points (or None if absent)
            "cF_feed_icp_mM": feed_icp_mM,
            "cF_diafiltrate_icp_mM": diafiltrate_icp_mM,
            "cF_retentate_icp_mM": retentate_icp_mM,
            "cF_final_icp_mM": final_tube_icp_mM,
            "icp_final_tube_mM": final_tube_icp_mM,
            "cF_final_meas": final_tube_icp_mM,
            "note_text": note_text,
        },
        "data_raw":   data_raw,
        "sheet_name": sheet_name,
    }
    data_stru = _normalize_conductivity_measurements(data_stru)
    if data_stru.get("data_raw"):
        try:
            first_cf = np.asarray(data_stru["data_raw"][0].get("cF_exp", []), dtype=float).reshape(-1)
            if data_stru["data_config"].get("C_F0_source") == "first_cF_exp_fallback":
                if first_cf.size and np.isfinite(first_cf[0]):
                    data_stru["data_config"]["C_F0"] = float(first_cf[0])
            elif (not np.isfinite(data_stru["data_config"].get("C_F0", np.nan))
                    and first_cf.size and np.isfinite(first_cf[0])):
                data_stru["data_config"]["C_F0"] = float(first_cf[0])
                data_stru["data_config"]["C_F0_source"] = "first_cF_exp_fallback"
        except Exception:
            pass

    # ------ Detect holdup and split vial 1 into [extra, real] -----------
    # Mirrors the DATA1 / DATA2 utility.py n_v0 / n_extra convention. The
    # parameter estimator integrates Jw across both vials but the objective
    # skips residuals for vials with index < n_v0.
    cfg = data_stru["data_config"]
    cfg["vial1_split_status"] = "skip_no_data"
    cfg["t_perm_start_s"]     = None
    cfg["n_holdup_samples"]   = 0
    cfg["holdup_boundary_source"] = None
    if data_raw:
        v1 = data_raw[0]
        t = np.asarray(v1["time"], dtype=float).reshape(-1)
        m = np.asarray(v1["mass"], dtype=float).reshape(-1)
        nN = min(t.size, m.size)
        threshold_g = 0.05      # ≈ one drop above baseline
        min_holdup  = 3
        min_real    = 10
        if nN >= (min_holdup + min_real + 1):
            t = t[:nN]; m = m[:nN]
            # Holdup boundary: a curated (DATA2-style) override wins; otherwise
            # auto-detect the first mass rise above baseline. Curated boundaries
            # are honored whenever they leave >=1 holdup and >=1 real sample;
            # auto-detected ones keep the stricter min_holdup / min_real guards.
            curated_s = _nf270_curated_holdup_end_s(Path(path).name, sheet_name)
            if curated_s is not None:
                i_start = int(np.argmin(np.abs(t - float(curated_s))))
                boundary_source = "curated"
                valid = 1 <= i_start <= nN - 1
            else:
                m_base = float(m[0])
                rose = (m - m_base) >= threshold_g
                i_start = int(np.argmax(rose)) if rose.any() else None
                boundary_source = "auto_mass_rise"
                valid = i_start is not None and i_start >= min_holdup and (nN - i_start) >= min_real
            if valid:
                t_perm_start = float(t[i_start])
                # Build holdup row (slice every per-sample array; per-vial
                # scalars get NaN since no ICP for the holdup).
                def _slice_row(orig_row, lo, hi, *, drop_per_vial_scalars):
                    new = {}
                    for key, val in orig_row.items():
                        try:
                            arr = np.asarray(val)
                        except Exception:
                            new[key] = val
                            continue
                        if arr.ndim >= 1 and arr.shape[0] == nN:
                            new[key] = arr[lo:hi].copy()
                        elif arr.ndim == 0 or arr.shape == () or arr.size == 1:
                            new[key] = float("nan") if drop_per_vial_scalars else val
                        else:
                            new[key] = val
                    return new
                holdup_row = _slice_row(v1, 0, i_start, drop_per_vial_scalars=True)
                real_row   = _slice_row(v1, i_start, nN, drop_per_vial_scalars=False)
                # Renumber and reset real-vial mass to start at ~0.
                holdup_row["number"] = 1
                real_row["number"]   = 2
                if isinstance(real_row.get("mass"), np.ndarray) and real_row["mass"].size:
                    real_row["mass"] = real_row["mass"] - real_row["mass"][0]
                # Also zero the holdup mass (it's tiny but consistent).
                if isinstance(holdup_row.get("mass"), np.ndarray) and holdup_row["mass"].size:
                    holdup_row["mass"] = holdup_row["mass"] - holdup_row["mass"][0]
                # Renumber subsequent vials.
                new_data_raw = [holdup_row, real_row]
                for j, r in enumerate(data_raw[1:], start=3):
                    r["number"] = j
                    new_data_raw.append(r)
                data_stru["data_raw"] = new_data_raw
                cfg["n"]                  = len(new_data_raw)
                cfg["n_v0"]               = 2
                cfg["n_extra"]            = 1
                cfg["n_h"]                = 1
                cfg["n_A"]                = 1
                cfg["vial1_split_status"] = "split"
                cfg["t_perm_start_s"]     = t_perm_start
                cfg["n_holdup_samples"]   = i_start
                cfg["holdup_boundary_source"] = boundary_source
            elif boundary_source == "curated":
                cfg["vial1_split_status"] = "skip_curated_out_of_range"
            elif i_start is None:
                cfg["vial1_split_status"] = "skip_no_rise"
            else:
                cfg["vial1_split_status"] = "skip_no_lag"

    # ------ Lp seed from vial-1 slope (data-driven) ---------------------
    # Use the first REAL vial (index n_v0, 1-indexed) for the slope fit.
    cfg["Lp0_method"]            = "fallback_no_real_vial"
    cfg["Lp0_original_default"]  = 5.0
    real_vial_idx = max(0, int(cfg.get("n_v0", 1)) - 1)
    if real_vial_idx < len(data_stru["data_raw"]):
        rv = data_stru["data_raw"][real_vial_idx]
        t_rv = np.asarray(rv.get("time", []), dtype=float).reshape(-1)
        m_rv = np.asarray(rv.get("mass", []), dtype=float).reshape(-1)
        nN = min(t_rv.size, m_rv.size)
        if nN >= 3:
            t_rv = t_rv[:nN]; m_rv = m_rv[:nN]
            m_max = float(np.nanmax(m_rv))
            if m_max >= 0.10:
                # Fit slope on the rising portion (mass between 0.1 g and 0.95 × m_max).
                keep = (m_rv >= 0.05) & (m_rv <= 0.95 * m_max)
                if int(keep.sum()) >= 3:
                    try:
                        slope, _intercept = np.polyfit(t_rv[keep], m_rv[keep], 1)
                        # Lp [L/m²/hr/bar] = slope × 36000 / (Am · ρ · ΔP)
                        # where Am in cm², ρ in g/cm³, ΔP in bar.
                        Am  = float(cfg.get("Am", 4.1))
                        rho = float(cfg.get("rho", 1.0))
                        if Am > 0 and rho > 0 and delp > 0 and np.isfinite(slope) and slope > 0:
                            lp = float(slope) * 36000.0 / (Am * rho * delp)
                            if 1.0 <= lp <= 50.0:
                                cfg["Lp0"] = lp
                                cfg["theta0"] = np.array([lp, cfg["B0"], cfg["sigma0"]], dtype=float)
                                cfg["Lp0_method"] = "data_driven"
                            else:
                                cfg["Lp0_method"] = "fallback_out_of_bounds"
                        else:
                            cfg["Lp0_method"] = "fallback_bad_inputs"
                    except Exception:
                        cfg["Lp0_method"] = "fallback_bad_slope"
                else:
                    cfg["Lp0_method"] = "fallback_too_few_points"
            else:
                cfg["Lp0_method"] = "fallback_no_rise"
    return data_stru


def estimate_parameters_with_parmest(
    data_structures,
    mode="DATA",
    theta=None,
    B_form="single",
    nfe=300,
    workflow_family="DATA1",
    weighted=True,
    solver_options=None,
    tee=False,
):
    """Estimate parameters using the current ParmEst API."""
    exp_list = build_parmest_experiments(
        data_structures, mode=mode, theta=theta, B_form=B_form, nfe=nfe, workflow_family=workflow_family
    )
    obj_function = "SSE_weighted" if weighted else "SSE"
    pest = Estimator(exp_list, obj_function=obj_function, tee=tee, solver_options=solver_options)

    obj_val, theta_vals = pest.theta_est()
    try:
        cov = pest.cov_est()
    except Exception as err:
        cov = None
        result = {
            "obj_val": obj_val,
            "theta_vals": theta_vals,
            "covariance": cov,
            "covariance_warning": f"Covariance skipped because ParmEst covariance failed: {type(err).__name__}: {err}",
        }
        try:
            std = np.sqrt(np.diag(cov))
            result["std"] = std
            result["correlation"] = correlation_from_covariance(cov.values if hasattr(cov, "values") else cov)
        except Exception:
            result["std"] = None
            result["correlation"] = None
        return result

    result = {
        "obj_val": obj_val,
        "theta_vals": theta_vals,
        "covariance": cov,
    }

    try:
        std = np.sqrt(np.diag(cov))
        result["std"] = std
        result["correlation"] = correlation_from_covariance(cov.values if hasattr(cov, "values") else cov)
    except Exception:
        result["std"] = None
        result["correlation"] = None

    return result


def compute_covariance_with_parmest(
    data_structures,
    mode="DATA",
    theta=None,
    B_form="single",
    nfe=300,
    method="finite_difference",
    solver="ipopt",
    tee=False,
    estimated_var=None,
):
    """Compute a covariance matrix with ParmEst."""
    exp_list = build_parmest_experiments(
        data_structures, mode=mode, theta=theta, B_form=B_form, nfe=nfe
    )
    obj_function = SSE_weighted if estimated_var is not None else SSE

    if estimated_var is None:
        pest = Estimator(exp_list, obj_function=obj_function, tee=tee)
        obj_val, theta_vals = pest.theta_est()
        cov = pest.cov_est(method=method)
    else:
        theta_vals = None
        obj_val = None
        cov = compute_covariance_matrix(
            exp_list,
            method=method,
            obj_function=obj_function,
            theta_vals=theta if theta is not None else {},
            step=0.02,
            solver=solver,
            tee=tee,
            estimated_var=estimated_var,
        )

    return {
        "obj_val": obj_val,
        "theta_vals": theta_vals,
        "covariance": cov,
    }


def compute_doe_metrics(experiment, step=0.001, objective_option="determinant", prior_FIM=None, solver=None, tee=False):
    """Use Pyomo.DoE to compute FIM metrics for a labeled experiment."""
    model = experiment.get_labeled_model()

    # Pyomo.DoE needs design inputs to run full DoE analysis.
    # If the model has not been extended with experiment_inputs yet, we stop
    # here with a clear message instead of failing later in a longer stack trace.
    if not hasattr(model, "experiment_inputs") or len(model.experiment_inputs) == 0:
        raise ValueError(
            "Pyomo.DoE needs experiment_inputs on the model. "
            "Add design variables to the model builder before calling compute_doe_metrics()."
        )

    doe = DesignOfExperiments(
        experiment=experiment,
        fd_formula="central",
        step=step,
        objective_option=objective_option,
        prior_FIM=prior_FIM,
        solver=solver,
        tee=tee,
    )
    fim = doe.compute_FIM(method="sequential")

    fim_np = np.asarray(fim, dtype=float)
    result = {
        "FIM": fim,
        "trace": float(np.trace(fim_np)),
        "det": float(np.linalg.det(fim_np)),
        "eig_val": np.linalg.eigvalsh(fim_np).tolist(),
    }
    try:
        result["covariance"] = np.linalg.inv(fim_np)
        result["std"] = np.sqrt(np.diag(result["covariance"])).tolist()
    except Exception:
        result["covariance"] = None
        result["std"] = None
    return result


def summarize_uncertainty(covariance):
    """Turn a covariance matrix into a small uncertainty summary."""
    cov_np = np.asarray(covariance, dtype=float)
    summary = {
        "covariance": cov_np.tolist(),
        "correlation": correlation_from_covariance(cov_np).tolist(),
        "std": np.sqrt(np.diag(cov_np)).tolist(),
    }
    return summary


def load_experimental_data(mat_file_path, results=None, selector=None):
    """Stage 1: load one experimental file and normalize the measured fields."""
    if results is None:
        results = {}

    # Allow a single experiment or a batch of related experiments.
    if isinstance(mat_file_path, (list, tuple, set)):
        sources = list(mat_file_path)
    else:
        sources = [mat_file_path]

    loaded = []
    runs = []
    for source in sources:
        run, data_stru = _normalize_experimental_source(source, selector=selector)
        loaded.append(data_stru)
        runs.append(run)

    results["data_file"] = [str(src) for src in sources]
    results["data"] = loaded[0] if len(loaded) == 1 else loaded
    results["experiments"] = runs[0] if len(runs) == 1 else runs
    results["is_batch"] = len(loaded) > 1
    return results


def build_model(results, mode="DATA", workflow_family="DATA1", B_form="single", nfe=300):
    """Stage 2: build the first-principles Pyomo model."""
    if "data" not in results:
        raise RuntimeError("Stage 1 must run before Stage 2.")

    data_payload = results["data"][0] if isinstance(results["data"], list) else results["data"]
    model = model_construct_inter(
        data_payload,
        mode,
        theta=None,
        sim_opt=False,
        B_form=B_form,
        workflow_family=workflow_family,
    )
    results["model"] = model
    results["model_settings"] = {
        "mode": mode,
        "workflow_family": workflow_family,
        "B_form": B_form,
        "nfe": int(nfe),
        "is_batch": bool(results.get("is_batch")),
    }
    return results


def estimate_parameters(
    results,
    *,
    use_parmest=False,
    multistart=False,
    multistart_iterations=10,
    nfe=300,
    tee=False,
):
    """Stage 3: estimate parameters with ParmEst or the legacy solver path."""
    if "data" not in results:
        raise RuntimeError("Stage 1 must run before Stage 3.")

    settings = results.get("model_settings", {})
    mode = settings.get("mode", "DATA")
    workflow_family = settings.get("workflow_family", "DATA1")
    B_form = settings.get("B_form", "single")
    data_payload = results["data"]
    if isinstance(data_payload, list):
        data_payload_list = data_payload
    else:
        data_payload_list = [data_payload]

    if use_parmest or len(data_payload_list) > 1:
        pest_result = estimate_parameters_with_parmest(
            data_payload_list,
            mode=mode,
            theta=None,
            B_form=B_form,
            nfe=nfe,
            workflow_family=workflow_family,
            tee=tee,
        )
        theta_vals = pest_result.get("theta_vals")
        if hasattr(theta_vals, "to_dict"):
            theta_vals = theta_vals.to_dict()
        results["parmest"] = pest_result
        results["parameters"] = theta_vals or {}
        fit_stru, sim_stru, sim_inter = solve_model(
            data_payload_list[0],
            mode,
            theta=results["parameters"],
            sim_opt=False,
            B_form=B_form,
            workflow_family=workflow_family,
            nfe=nfe,
            multistart=multistart,
            multistart_iterations=multistart_iterations,
        )
    else:
        fit_stru, sim_stru, sim_inter = solve_model(
            data_payload_list[0],
            mode,
            theta=None,
            sim_opt=False,
            B_form=B_form,
            workflow_family=workflow_family,
            nfe=nfe,
            multistart=multistart,
            multistart_iterations=multistart_iterations,
        )
        results["parameters"] = fit_stru.get("parameters", {}) if fit_stru else {}

    results["fit_stru"] = fit_stru
    results["sim_stru"] = sim_stru
    results["sim_inter"] = sim_inter
    results["multistart"] = fit_stru.get("multistart", {
        "enabled": bool(multistart),
        "iterations": int(multistart_iterations),
    }) if isinstance(fit_stru, dict) else {
        "enabled": bool(multistart),
        "iterations": int(multistart_iterations),
    }
    return results


def quantify_uncertainty(
    results,
    *,
    method="fim",
    cov_method="finite_difference",
    fim_step=1e-8,
    fim_formula="backward",
    nfe=300,
):
    """Stage 4: quantify uncertainty with covariance or FIM."""
    if "parameters" not in results:
        raise RuntimeError("Stage 3 must run before Stage 4.")

    settings = results.get("model_settings", {})
    mode = settings.get("mode", "DATA")
    B_form = settings.get("B_form", "single")
    workflow_family = settings.get("workflow_family", "DATA1")
    data_stru = results["data"][0] if isinstance(results["data"], list) else results["data"]

    if method == "cov_est":
        pest_result = results.get("parmest")
        if pest_result and pest_result.get("covariance") is not None:
            results["uncertainty"] = summarize_uncertainty(pest_result["covariance"])
        else:
            results["uncertainty"] = {"status": "skipped", "reason": "ParmEst covariance not available."}

    elif method == "fim":
        fim = calc_FIM(
            data_stru,
            mode,
            theta=results["parameters"],
            step=fim_step,
            formula=fim_formula,
            B_form=B_form,
            workflow_family=workflow_family,
            nfe=nfe,
        )
        results["uncertainty"] = {
            "method": "calc_FIM",
            "FIM": fim.get("FIM"),
            "trace": fim.get("trace"),
            "det": fim.get("det"),
            "eig_val": fim.get("eig_val"),
            "covariance": fim.get("V"),
            "std": fim.get("std"),
        }

    else:
        results["uncertainty"] = {"status": "skipped", "reason": f"Unknown method {method!r}."}

    return results


def design_next_experiment(results, **_kwargs):
    """Stage 5: keep the hook, but skip MBDoE until design variables exist."""
    results["mbdoe"] = {
        "status": "skipped",
        "reason": "Design variables are not exposed yet; the workflow stays focused on fitting and reproduction.",
    }
    return results


def run_workflow(
    mat_file_path,
    *,
    mode="DATA",
    workflow_family="DATA1",
    B_form="single",
    selector=None,
    use_parmest=False,
    multistart=False,
    multistart_iterations=10,
    nfe=300,
    uncertainty_method="fim",
    cov_method="finite_difference",
    fim_step=1e-8,
    fim_formula="backward",
    skip_mbdoe=True,
):
    """Run the compact stage-based workflow on one or many experimental files."""
    print("=" * 70)
    print(" Diafiltration workflow")
    print(f"   File: {mat_file_path}")
    print(f"   Family: {workflow_family}")
    print(f"   Mode: {mode}")
    print("=" * 70)

    results = {}
    results = load_experimental_data(mat_file_path, results, selector=selector)
    results = build_model(results, mode=mode, workflow_family=workflow_family, B_form=B_form, nfe=nfe)
    results = estimate_parameters(
        results,
        use_parmest=use_parmest,
        multistart=multistart,
        multistart_iterations=multistart_iterations,
        nfe=nfe,
    )
    results = quantify_uncertainty(
        results,
        method=uncertainty_method,
        cov_method=cov_method,
        fim_step=fim_step,
        fim_formula=fim_formula,
        nfe=nfe,
    )
    if not skip_mbdoe:
        results = design_next_experiment(results)

    print("\n[Done] Workflow finished. Results available:")
    for key in results:
        if not key.startswith("_"):
            print(f"   - {key}")
    return results


def _sanitize_filename_component(value):
    """Return a filesystem-safe token for generated figure names."""
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value or "").strip())
    return token.strip("._-") or "data3"


def _plot_time_minutes(time_array, *, t_delay, extra_offset=0.0):
    """Convert a raw second-based time array into shifted minutes."""
    return (np.asarray(time_array, dtype=float) - float(t_delay) - float(extra_offset)) / 60.0


def _normalize_sim_stru(sim_stru):
    """Return simulation results as an ordered list of per-vial dicts."""
    if isinstance(sim_stru, dict):
        items = []
        for key, value in sim_stru.items():
            try:
                sort_key = int(key)
            except Exception:
                sort_key = str(key)
            items.append((sort_key, value))
        items.sort(key=lambda item: item[0])
        return [value for _, value in items if isinstance(value, dict)]
    if isinstance(sim_stru, (list, tuple)):
        return [value for value in sim_stru if isinstance(value, dict)]
    return []


def _build_fit_meta_annotation(fit_meta, include_cf_unc=False):
    """Render a short annotation block from a fit_meta dict.

    Returns a list of strings (one per line) suitable for matplotlib `text()`.
    Returns [] if fit_meta is empty or has no usable keys, so the caller can
    skip the annotation cleanly.

    Keys consumed (all optional):
      recipe_name, b_form, theta_init_summary, winning_theta, cf_uncertainty_str
    """
    if not isinstance(fit_meta, dict) or not fit_meta:
        return []
    lines = []
    if fit_meta.get("recipe_name"):
        lines.append(f"recipe: {fit_meta['recipe_name']}")
    if fit_meta.get("b_form") is not None:
        lines.append(f"B_form = {fit_meta['b_form']}")
    if fit_meta.get("theta_init_summary"):
        for ln in str(fit_meta["theta_init_summary"]).splitlines():
            if ln.strip():
                lines.append(ln)
    wt = fit_meta.get("winning_theta") or {}
    if isinstance(wt, dict) and wt:
        bits = []
        for k in ("Lp", "B", "beta_0", "beta_1", "sigma"):
            if k in wt and isinstance(wt[k], (int, float)) and np.isfinite(wt[k]):
                bits.append(f"{k}={wt[k]:.3g}")
        if bits:
            lines.append("fit: " + ", ".join(bits))
    if include_cf_unc and fit_meta.get("cf_uncertainty_str"):
        lines.append(str(fit_meta["cf_uncertainty_str"]))
    return lines


def run_data3_time_series_plots(results, save_dir=None, show=False):
    """Create DATA3-only mass and concentration time-series plots from a staged workflow result.

    Optional ``results['fit_meta']`` keys read for plot annotations:
      - ``theta_init_summary``  (str)   — seed strategy description
      - ``cf_uncertainty_str``  (str)   — formula used for σ_cF
      - ``b_form``              (any)   — B_form value (1, 2, 'single', ...)
      - ``winning_theta``       (dict)  — fitted (Lp, B, sigma, ...) of the
                                          winning trial
      - ``recipe_name``         (str)   — short label (e.g. "slide-29 v5")
    """
    settings = results.get("model_settings", {})
    if str(settings.get("workflow_family", "DATA3")).upper() != "DATA3":
        return []

    data_payload = results.get("data")
    if isinstance(data_payload, list):
        data_payload = data_payload[0] if data_payload else None
    if not isinstance(data_payload, dict):
        return []

    # Optional fit-meta for the new init-guess + cF-uncertainty annotations
    # (2026-05-25).  Keys are all-optional; missing keys → no annotation
    # change vs the legacy behavior.
    fit_meta = results.get("fit_meta") or {}

    sim_stru = _normalize_sim_stru(results.get("sim_stru") or [])
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data3_option3"
    save_dir.mkdir(parents=True, exist_ok=True)

    sheet_slug = _sanitize_filename_component(data_payload.get("sheet_name", "sheet"))
    file_slug = _sanitize_filename_component(Path(str(data_payload.get("filename", "data3"))).stem)
    prefix = f"{file_slug}_{sheet_slug}"
    t_delay = _data1_time_origin(data_payload)
    n_vials = int(data_payload.get("data_config", {}).get("n", len(data_payload.get("data_raw", []))))
    n_v0 = int(data_payload.get("data_config", {}).get("n_v0", 1))
    outputs = []

    # Mass vs time
    mass_fig, mass_ax = plt.subplots(figsize=(4, 4))
    for i in range(min(n_vials, len(data_payload.get("data_raw", [])))):
        row = data_payload["data_raw"][i]
        time = np.asarray(row.get("time", []), dtype=float)
        mass = np.asarray(row.get("mass", []), dtype=float)
        if time.size and mass.size:
            mass_rel = mass - mass[0]
            mass_ax.plot(_plot_time_minutes(time, t_delay=t_delay), mass_rel, "r.", markersize=4)
        if i < len(sim_stru) and (i + 1) >= n_v0:
            sim_time = np.asarray(sim_stru[i].get("time", []), dtype=float)
            sim_mass = np.asarray(sim_stru[i].get("mV", []), dtype=float)
            if sim_time.size and sim_mass.size:
                sim_mass_rel = sim_mass - sim_mass[0]
                mass_ax.plot(_plot_time_minutes(sim_time, t_delay=t_delay), sim_mass_rel, "b", linewidth=3, alpha=0.6)
    mass_ax.plot([], [], "r.", markersize=4, label="Measurements")
    mass_ax.plot([], [], "b", linewidth=3, alpha=0.6, label="Predictions")
    mass_ax.set_xlabel("Time [min]", fontsize=16, fontweight="bold")
    mass_ax.set_ylabel("Mass [g]", fontsize=16, fontweight="bold")
    mass_ax.tick_params(direction="in")
    mass_ax.xaxis.set_tick_params(labelsize=15)
    mass_ax.yaxis.set_tick_params(labelsize=15)
    mass_ax.set_xlim(left=0)
    mass_ax.set_ylim(bottom=0)
    mass_ax.legend(fontsize=10, loc="best")
    # Annotate the per-sheet WSSE mass-residual uncertainty.  The current
    # loader keeps the DATA2-style 0.01 g convention, but if a caller has
    # injected decomposition metadata we still render it here for context.
    _cfg_for_annot = data_payload.get("data_config", {}) if isinstance(data_payload, dict) else {}
    _mass_scale_g = float(_cfg_for_annot.get("mass_scale_g", 0.01))
    if not np.isfinite(_mass_scale_g) or _mass_scale_g <= 0:
        _mass_scale_g = 0.01
    _components = _cfg_for_annot.get("mass_scale_components") or {}
    if _components:
        _sig = float(_components.get("signal_g", 0.0))
        _sb  = float(_components.get("sigma_balance_g", MASS_BALANCE_SIGMA_G))
        _sm  = float(_components.get("sigma_model_g", 0.0))
        _frac = float(_components.get("model_fraction", MASS_MODEL_FIT_TOLERANCE))
        _ann_text = (
            f"WSSE mass σ = {_mass_scale_g:.3g} g\n"
            f"  σ_balance = {_sb:.3g} g\n"
            f"  σ_model  = {_frac:.2g} × |Δm| = {_sm:.3g} g\n"
            f"  |Δm|      = {_sig:.3g} g"
        )
    else:
        _ann_text = f"WSSE mass σ = {_mass_scale_g:.3g} g (legacy balance-only)"
    mass_ax.text(
        0.02, 0.97,
        _ann_text,
        transform=mass_ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        horizontalalignment="left",
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", alpha=0.8),
    )
    # 2026-05-25: optional init-guess + recipe annotation (bottom-right of mass plot)
    _fit_lines = _build_fit_meta_annotation(fit_meta, include_cf_unc=False)
    if _fit_lines:
        mass_ax.text(
            0.98, 0.03,
            "\n".join(_fit_lines),
            transform=mass_ax.transAxes,
            fontsize=7,
            verticalalignment="bottom",
            horizontalalignment="right",
            family="monospace",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", alpha=0.8),
        )
    mass_path = save_dir / f"mass-{prefix}.png"
    mass_fig.savefig(mass_path, dpi=300, bbox_inches="tight")
    outputs.append(str(mass_path))
    if show:
        plt.show()
    plt.close(mass_fig)

    # Concentration vs time — paper Figure 6 panel B/E styling
    # DATA3 (2026-06-19b): one color per STREAM by user request — retentate (data + ICP
    # + prediction) = TEAL; permeate (data + prediction) = RED.  Guarded: only this DATA3 fn uses them.
    PAPER_RET_FILL         = "#008B8B"   # retentate — data squares + ICP diamonds (teal)
    PAPER_VIAL_FILL        = "#D62728"   # permeate  — vial data squares (red)
    PAPER_RET_LINE         = "#008B8B"   # retentate — prediction line (teal, unified w/ data)
    PAPER_PERM_LINE        = "#D62728"   # permeate  — prediction line (red, unified w/ data)
    PAPER_EDGE_COLOR       = "k"
    PAPER_EDGE_WIDTH       = 0.6         # vial squares + prediction triangles
    PAPER_RET_MARKER_SIZE  = 4           # dense retentate trace — small, no edge
    PAPER_VIAL_MARKER_SIZE = 7           # sparse vial squares + triangles
    conc_fig, conc_ax = plt.subplots(figsize=(6, 4))
    data_raw = data_payload.get("data_raw", [])
    seen_ret_meas = False
    seen_perm_meas = False
    seen_ret_pred = False
    seen_perm_pred = False
    seen_vial_pred = False
    seen_ret_icp = False
    # DATA3 concentration-plot refinements (2026-06-19): retentate-ICP 2-pt anchor (feed=start,
    # cF_retentate_icp_mM=end) + lag boundary for zeroing the predicted permeate before permeate
    # begins.  All reads are NaN/None-guarded so a sheet missing the key degrades gracefully.
    _cfg = data_payload.get("data_config", {})
    _ret_icp_start = float(_cfg.get("C_F0", np.nan))
    _ret_icp_end = float(_cfg.get("cF_retentate_icp_mM", np.nan))
    _last_t = np.asarray(data_raw[-1].get("time", []), dtype=float).reshape(-1) if data_raw else np.array([])
    _ret_icp_t_start_min = 0.0
    _ret_icp_t_end_min = float(_plot_time_minutes(_last_t[-1], t_delay=t_delay)) if _last_t.size else np.nan
    _t_perm_start_s = _cfg.get("t_perm_start_s")
    _lag_boundary_min = (_plot_time_minutes(_t_perm_start_s, t_delay=t_delay) if _t_perm_start_s is not None else None)
    for i, row in enumerate(data_raw):
        time = np.asarray(row.get("time", []), dtype=float).reshape(-1)
        if time.size == 0:
            continue
        cF_exp = np.asarray(row.get("cF_exp", []), dtype=float).reshape(-1)
        cV_avg = np.asarray(row.get("cV_avg", []), dtype=float)
        cV_avg_arr = np.atleast_1d(cV_avg).astype(float, copy=False)

        if cF_exp.size:
            if cF_exp.size == 1:
                conc_ax.plot(
                    _plot_time_minutes(time[-1], t_delay=t_delay),
                    float(cF_exp[0]),
                    marker="s", linestyle="None",
                    color=PAPER_RET_FILL,
                    markersize=PAPER_RET_MARKER_SIZE,
                    markeredgewidth=0,
                    clip_on=False,
                    label="Retentate (Measurements)" if not seen_ret_meas else None,
                )
            else:
                n = min(len(time), len(cF_exp))
                conc_ax.plot(
                    _plot_time_minutes(time[:n], t_delay=t_delay),
                    cF_exp[:n],
                    marker="s", linestyle="None",
                    color=PAPER_RET_FILL,
                    markersize=PAPER_RET_MARKER_SIZE,
                    markeredgewidth=0,
                    clip_on=False,
                    label="Retentate (Measurements)" if not seen_ret_meas else None,
                )
            seen_ret_meas = True

        if cV_avg_arr.size:
            if cV_avg_arr.ndim == 0 or cV_avg_arr.size == 1:
                conc_ax.plot(
                    _plot_time_minutes(time[-1], t_delay=t_delay),
                    float(cV_avg_arr.reshape(-1)[0]),
                    marker="s", linestyle="None",
                    color=PAPER_VIAL_FILL,
                    markersize=PAPER_VIAL_MARKER_SIZE,
                    markeredgecolor=PAPER_EDGE_COLOR,
                    markeredgewidth=PAPER_EDGE_WIDTH,
                    clip_on=False,
                    label="Vial (Measurements)" if not seen_perm_meas else None,
                )
            else:
                n = min(len(time), len(cV_avg_arr))
                conc_ax.plot(
                    _plot_time_minutes(time[:n], t_delay=t_delay),
                    cV_avg_arr[:n],
                    marker="s", linestyle="None",
                    color=PAPER_VIAL_FILL,
                    markersize=PAPER_VIAL_MARKER_SIZE,
                    markeredgecolor=PAPER_EDGE_COLOR,
                    markeredgewidth=PAPER_EDGE_WIDTH,
                    clip_on=False,
                    label="Vial (Measurements)" if not seen_perm_meas else None,
                )
            seen_perm_meas = True

        if i >= len(sim_stru) or (i + 1) < n_v0:
            continue

        sim_time = np.asarray(sim_stru[i].get("time", []), dtype=float).reshape(-1)
        if sim_time.size == 0:
            continue
        cF = np.asarray(sim_stru[i].get("cF", []), dtype=float).reshape(-1)
        cH = np.asarray(sim_stru[i].get("cH", []), dtype=float).reshape(-1)
        cV = np.asarray(sim_stru[i].get("cV", []), dtype=float).reshape(-1)
        if cF.size:
            conc_ax.plot(
                _plot_time_minutes(sim_time, t_delay=t_delay),
                cF,
                color=PAPER_RET_LINE,
                linewidth=2.5,
                alpha=1.0,
                label="Retentate (Predictions)" if not seen_ret_pred else None,
            )
            seen_ret_pred = True
        # DATA3: retentate ICP — 2-pt scatter (start=feed, end=cF_retentate_icp_mM), drawn once.
        if not seen_ret_icp:
            _icp_t, _icp_c = [], []
            if np.isfinite(_ret_icp_start):
                _icp_t.append(_ret_icp_t_start_min); _icp_c.append(_ret_icp_start)
            if np.isfinite(_ret_icp_end) and np.isfinite(_ret_icp_t_end_min):
                _icp_t.append(_ret_icp_t_end_min); _icp_c.append(_ret_icp_end)
            if _icp_c:
                conc_ax.plot(_icp_t, _icp_c, marker="D", linestyle="None",
                             color=PAPER_RET_FILL, markersize=PAPER_VIAL_MARKER_SIZE,
                             markeredgecolor=PAPER_EDGE_COLOR, markeredgewidth=PAPER_EDGE_WIDTH,
                             clip_on=False, label="Retentate (ICP)")
                seen_ret_icp = True
        if cH.size:
            # DATA3: zero the predicted permeate through the lag startup (no permeate yet).
            _cH_t = _plot_time_minutes(sim_time, t_delay=t_delay)
            _cH_plot = cH.copy()
            if _lag_boundary_min is not None:
                _cH_plot = np.where(np.asarray(_cH_t) < _lag_boundary_min, 0.0, _cH_plot)
            conc_ax.plot(
                _cH_t,
                _cH_plot,
                color=PAPER_PERM_LINE,
                linestyle="-",
                linewidth=2.5,
                alpha=1.0,
                label="Permeate (Predictions)" if not seen_perm_pred else None,
            )
            seen_perm_pred = True
        # DATA3: per-vial permeate PREDICTION markers removed — continuous cH line is the only
        # permeate prediction now (change 3).
    conc_ax.set_xlabel("Time [min]", fontsize=16, fontweight="bold")
    conc_ax.set_ylabel("Concentration [mM]", fontsize=16, fontweight="bold")
    conc_ax.tick_params(direction="in")
    conc_ax.xaxis.set_tick_params(labelsize=15)
    conc_ax.yaxis.set_tick_params(labelsize=15)
    conc_ax.set_xlim(left=0)
    conc_ax.set_ylim(bottom=0)
    conc_ax.legend(
        fontsize=9,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=True,
        borderaxespad=0.0,
    )
    # 2026-05-25: optional init-guess + cF-uncertainty annotation on conc plot
    _conc_fit_lines = _build_fit_meta_annotation(fit_meta, include_cf_unc=True)
    if _conc_fit_lines:
        conc_ax.text(
            0.02, 0.97,
            "\n".join(_conc_fit_lines),
            transform=conc_ax.transAxes,
            fontsize=7,
            verticalalignment="top",
            horizontalalignment="left",
            family="monospace",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", alpha=0.8),
        )
    conc_path = save_dir / f"concentration-{prefix}.png"
    conc_fig.savefig(conc_path, dpi=300, bbox_inches="tight")
    outputs.append(str(conc_path))
    if show:
        plt.show()
    plt.close(conc_fig)

    return outputs


def run_data3_pressure_plots(results, save_dir=None, show=False):
    """Create DATA3 applied-pressure and osmotic-pressure time-series plots.

    Two output figures per sheet:
      1. ``pressure-<prefix>.png``  —  Applied ΔP (constant line) + the
         osmotic backpressure σ·Δπ(t) overlaid + net driving force
         (ΔP − σ·Δπ) vs time.  Shows the membrane work term evolution.
      2. ``osmotic-<prefix>.png``  —  Pure Δπ(t) curve per vial.

    Δπ is computed from sim_stru's bulk feed (cF) and permeate-wall (cH)
    concentrations:
        Δπ(t) = (cF(t) − cH(t)) · ni · R · T
    This omits the concentration-polarization correction (cIn vs cF) —
    the surface-side Δπ would be ~e^(Jw/k) larger.  See comment block in
    model_construct_inter for the polarization equation.

    Inputs: same ``results`` dict structure as run_data3_time_series_plots.
    Optional ``results['fit_meta']`` controls the annotation block.

    Returns list of saved image paths.
    """
    settings = results.get("model_settings", {})
    if str(settings.get("workflow_family", "DATA3")).upper() != "DATA3":
        return []

    data_payload = results.get("data")
    if isinstance(data_payload, list):
        data_payload = data_payload[0] if data_payload else None
    if not isinstance(data_payload, dict):
        return []

    fit_meta = results.get("fit_meta") or {}
    sim_stru = _normalize_sim_stru(results.get("sim_stru") or [])
    if not sim_stru:
        return []
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data3_option3"
    save_dir.mkdir(parents=True, exist_ok=True)

    sheet_slug = _sanitize_filename_component(data_payload.get("sheet_name", "sheet"))
    file_slug = _sanitize_filename_component(Path(str(data_payload.get("filename", "data3"))).stem)
    prefix = f"{file_slug}_{sheet_slug}"
    t_delay = _data1_time_origin(data_payload)
    n_vials = int(data_payload.get("data_config", {}).get("n", len(data_payload.get("data_raw", []))))
    n_v0 = int(data_payload.get("data_config", {}).get("n_v0", 1))

    cfg = data_payload.get("data_config", {})
    delP = float(cfg.get("delP", 0.0))     # bar
    ni   = float(cfg.get("ni",   1.0))     # van 't Hoff
    Temp = float(cfg.get("Temp", 298.15))  # K
    R    = 8.314e-5                        # cm^3 bar / micromol / K

    # Use the winning sigma if available; otherwise 1.0 for upper-bound trace
    winning = fit_meta.get("winning_theta") or {}
    sigma_fit = float(winning.get("sigma", 1.0)) if isinstance(winning.get("sigma"), (int, float)) else 1.0

    outputs = []

    # ───── Figure 1: applied ΔP, σ·Δπ, and net driving force vs time ─────
    pres_fig, pres_ax = plt.subplots(figsize=(6, 4))
    pres_ax.axhline(delP, color="tab:orange", linestyle="-", linewidth=2.5,
                    label=f"Applied ΔP = {delP:.2f} bar", alpha=0.9)
    # Aggregate dpi & net curves across vials
    any_drawn = False
    for i in range(min(n_vials, len(sim_stru))):
        if (i + 1) < n_v0:
            continue
        sim_time = np.asarray(sim_stru[i].get("time", []), dtype=float)
        cF_arr   = np.asarray(sim_stru[i].get("cF",   []), dtype=float)
        cH_arr   = np.asarray(sim_stru[i].get("cH",   []), dtype=float)
        if not (sim_time.size and cF_arr.size and cH_arr.size):
            continue
        # Truncate to common length
        m = min(sim_time.size, cF_arr.size, cH_arr.size)
        t_min = _plot_time_minutes(sim_time[:m], t_delay=t_delay)
        # Δπ in bar.  R uses [cm^3 bar / μmol / K]; cF and cH are in mM = μmol/cm^3.
        dpi_t   = (cF_arr[:m] - cH_arr[:m]) * ni * R * Temp
        sdpi_t  = sigma_fit * dpi_t
        net_t   = delP - sdpi_t
        pres_ax.plot(t_min, sdpi_t, color="tab:red", linewidth=1.5, alpha=0.6,
                     label=("σ·Δπ (bulk)" if not any_drawn else None))
        pres_ax.plot(t_min, net_t, color="tab:blue", linewidth=1.5, alpha=0.6,
                     label=("Net = ΔP − σ·Δπ" if not any_drawn else None))
        any_drawn = True
    pres_ax.set_xlabel("Time [min]", fontsize=14, fontweight="bold")
    pres_ax.set_ylabel("Pressure [bar]", fontsize=14, fontweight="bold")
    pres_ax.tick_params(direction="in")
    pres_ax.legend(fontsize=9, loc="best")
    pres_ax.set_xlim(left=0)
    # Annotation block on pressure plot
    _pres_lines = _build_fit_meta_annotation(fit_meta, include_cf_unc=False)
    _pres_lines = [f"ni = {ni:g}, T = {Temp:.1f} K"] + _pres_lines
    pres_ax.text(
        0.02, 0.97, "\n".join(_pres_lines),
        transform=pres_ax.transAxes, fontsize=7,
        verticalalignment="top", horizontalalignment="left", family="monospace",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", alpha=0.8),
    )
    pres_path = save_dir / f"pressure-{prefix}.png"
    pres_fig.savefig(pres_path, dpi=300, bbox_inches="tight")
    outputs.append(str(pres_path))
    if show:
        plt.show()
    plt.close(pres_fig)

    # ───── Figure 2: pure Δπ(t) ─────
    osm_fig, osm_ax = plt.subplots(figsize=(6, 4))
    any_drawn = False
    for i in range(min(n_vials, len(sim_stru))):
        if (i + 1) < n_v0:
            continue
        sim_time = np.asarray(sim_stru[i].get("time", []), dtype=float)
        cF_arr   = np.asarray(sim_stru[i].get("cF",   []), dtype=float)
        cH_arr   = np.asarray(sim_stru[i].get("cH",   []), dtype=float)
        if not (sim_time.size and cF_arr.size and cH_arr.size):
            continue
        m = min(sim_time.size, cF_arr.size, cH_arr.size)
        t_min = _plot_time_minutes(sim_time[:m], t_delay=t_delay)
        dpi_t = (cF_arr[:m] - cH_arr[:m]) * ni * R * Temp
        osm_ax.plot(t_min, dpi_t, color="tab:purple", linewidth=1.5, alpha=0.7,
                    label=("Δπ (bulk)" if not any_drawn else None))
        any_drawn = True
    osm_ax.set_xlabel("Time [min]", fontsize=14, fontweight="bold")
    osm_ax.set_ylabel("Osmotic Pressure Δπ [bar]", fontsize=14, fontweight="bold")
    osm_ax.tick_params(direction="in")
    osm_ax.legend(fontsize=9, loc="best")
    osm_ax.set_xlim(left=0)
    osm_ax.set_ylim(bottom=0)
    _osm_lines = _build_fit_meta_annotation(fit_meta, include_cf_unc=True)
    _osm_lines = [
        f"Δπ = (cF − cH) · ni · R · T",
        f"ni = {ni:g}, T = {Temp:.1f} K",
        "(bulk; polarization not included)",
    ] + _osm_lines
    osm_ax.text(
        0.02, 0.97, "\n".join(_osm_lines),
        transform=osm_ax.transAxes, fontsize=7,
        verticalalignment="top", horizontalalignment="left", family="monospace",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", alpha=0.8),
    )
    osm_path = save_dir / f"osmotic-{prefix}.png"
    osm_fig.savefig(osm_path, dpi=300, bbox_inches="tight")
    outputs.append(str(osm_path))
    if show:
        plt.show()
    plt.close(osm_fig)

    return outputs


def run_data3_conductivity_plots(results, save_dir=None, show=False):
    """Create DATA3 raw-conductivity time-series plots — matches the
    collaborator's left-panel format from the multicomponent deck.

    Produces ``conductivity-<prefix>.png`` per sheet showing:
      - Magenta triangles: raw retentate conductivity (μS/cm at measured T)
      - Red diamonds: raw permeate conductivity (μS/cm at measured T)
      - (optionally) the EC25-compensated traces as faded squares so the
        compensation transformation is visible

    The "raw" values are recovered from the stored EC25-compensated
    cF_exp_conductivity / cV_perm_cond by inverting the EC25 formula:
        σ_T = σ_25 / (1 + α · (25 − T))

    This directly addresses the "you're not plotting what's in the data
    file" criticism — every continuous-measurement column from the raw
    Excel is now rendered in its native μS/cm units.

    Annotation calls out:
      - the EC25 formula and α used,
      - the data-file columns the points come from,
      - the Shedlovsky inversion that downstream produces concentration.
    """
    settings = results.get("model_settings", {})
    if str(settings.get("workflow_family", "DATA3")).upper() != "DATA3":
        return []

    data_payload = results.get("data")
    if isinstance(data_payload, list):
        data_payload = data_payload[0] if data_payload else None
    if not isinstance(data_payload, dict):
        return []

    fit_meta = results.get("fit_meta") or {}
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data3_option3"
    save_dir.mkdir(parents=True, exist_ok=True)

    sheet_slug = _sanitize_filename_component(data_payload.get("sheet_name", "sheet"))
    file_slug = _sanitize_filename_component(Path(str(data_payload.get("filename", "data3"))).stem)
    prefix = f"{file_slug}_{sheet_slug}"
    n_vials = int(data_payload.get("data_config", {}).get("n", len(data_payload.get("data_raw", []))))

    # Salt-specific EC25 alpha so the un-compensation matches the loader's
    # forward path.  Fallback to default if unknown salt.
    salt_name = str(data_payload.get("data_config", {}).get("namec", "NaCl"))
    alpha = float(CONDUCTIVITY_TEMP_COEFF_PER_C.get(salt_name, CONDUCTIVITY_TEMP_COEFF_DEFAULT))

    data_raw = data_payload.get("data_raw", []) or []
    if not data_raw:
        return []

    outputs = []

    cond_fig, cond_ax = plt.subplots(figsize=(7, 4.5))
    drew_ret_raw = False
    drew_perm_raw = False
    drew_ret_comp = False
    drew_perm_comp = False
    for i in range(min(n_vials, len(data_raw))):
        v = data_raw[i]
        time = np.asarray(v.get("time", []), dtype=float)
        ret_T = np.asarray(v.get("retentate_temp", []), dtype=float)
        perm_T = np.asarray(v.get("permeate_temp", []), dtype=float)
        ret_cond_comp = np.asarray(v.get("cF_exp_conductivity", []), dtype=float)
        perm_cond_comp = np.asarray(v.get("cV_perm_cond", []), dtype=float)
        n_use = min(time.size, ret_T.size, ret_cond_comp.size)
        if n_use == 0:
            continue
        # Invert EC25: σ_T = σ_25 / (1 + α · (25 − T))
        denom_ret = 1.0 + alpha * (25.0 - ret_T[:n_use])
        denom_ret = np.where(np.abs(denom_ret) > 1e-9, denom_ret, np.nan)
        ret_raw = ret_cond_comp[:n_use] / denom_ret
        n_perm = min(time.size, perm_T.size, perm_cond_comp.size)
        denom_perm = 1.0 + alpha * (25.0 - perm_T[:n_perm])
        denom_perm = np.where(np.abs(denom_perm) > 1e-9, denom_perm, np.nan)
        perm_raw = perm_cond_comp[:n_perm] / denom_perm

        # Plot raw (at measured T) — matches the collaborator's marker style.
        # X axis in seconds to match their format.
        cond_ax.plot(time[:n_use], ret_raw, "^", color="#C71585", markersize=4, alpha=0.85,
                     markeredgecolor="k", markeredgewidth=0.3,
                     label="Retentate (raw, μS/cm @ measured T)" if not drew_ret_raw else None)
        drew_ret_raw = True
        cond_ax.plot(time[:n_perm], perm_raw, "D", color="#D62728", markersize=4, alpha=0.85,
                     markeredgecolor="k", markeredgewidth=0.3,
                     label="Permeate (raw, μS/cm @ measured T)" if not drew_perm_raw else None)
        drew_perm_raw = True
        # Plot EC25-compensated (what we actually use downstream) as faded line
        cond_ax.plot(time[:n_use], ret_cond_comp[:n_use], "-", color="#2E8B57", linewidth=1.0, alpha=0.5,
                     label="Retentate EC25 (μS/cm @ 25°C)" if not drew_ret_comp else None)
        drew_ret_comp = True
        cond_ax.plot(time[:n_perm], perm_cond_comp[:n_perm], "-", color="#212121", linewidth=1.0, alpha=0.5,
                     label="Permeate EC25 (μS/cm @ 25°C)" if not drew_perm_comp else None)
        drew_perm_comp = True

    cond_ax.set_xlabel("Time [s]", fontsize=14, fontweight="bold")
    cond_ax.set_ylabel("Conductivity [μS/cm]", fontsize=14, fontweight="bold")
    cond_ax.tick_params(direction="in")
    cond_ax.xaxis.set_tick_params(labelsize=12)
    cond_ax.yaxis.set_tick_params(labelsize=12)
    cond_ax.set_xlim(left=0)
    cond_ax.set_ylim(bottom=0)
    # Legend anchored OUTSIDE the plot on the right so it doesn't collide
    # with the data or the source-annotation box.
    cond_ax.legend(
        fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        frameon=True, borderaxespad=0.0,
    )

    # Annotation: data-source explanation — top-left INSIDE the plot
    _ann = (
        f"Source: raw Excel cols 4 (retentate) + 6 (permeate)\n"
        f"Salt: {salt_name},  α = {alpha:.3f} /°C\n"
        f"EC25 forward: σ_25 = σ_T · (1 + α·(25−T))\n"
        f"EC25 inverse: σ_T  = σ_25 / (1 + α·(25−T))\n"
        f"Downstream: Shedlovsky inversion → mM"
    )
    cond_ax.text(
        0.02, 0.97, _ann, transform=cond_ax.transAxes, fontsize=7,
        verticalalignment="top", horizontalalignment="left", family="monospace",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", alpha=0.85),
    )
    # fit_meta annotation (recipe / cF uncertainty) — bottom-left INSIDE the plot
    _fit_lines = _build_fit_meta_annotation(fit_meta, include_cf_unc=True)
    if _fit_lines:
        cond_ax.text(
            0.02, 0.03, "\n".join(_fit_lines), transform=cond_ax.transAxes, fontsize=7,
            verticalalignment="bottom", horizontalalignment="left", family="monospace",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", alpha=0.85),
        )

    cond_path = save_dir / f"conductivity-{prefix}.png"
    cond_fig.savefig(cond_path, dpi=300, bbox_inches="tight")
    outputs.append(str(cond_path))
    if show:
        plt.show()
    plt.close(cond_fig)
    return outputs


def run_data3_applied_vs_osmotic_plots(results, save_dir=None, show=False):
    """Per-sheet plot of applied pressure ΔP + osmotic pressure Δπ on the
    same time axis, computed from the MEASURED retentate concentration (no
    fit required).

    Output filename: ``applied_vs_osmotic-<prefix>.png``

    Physical story this plot tells (from collaborator feedback 2026-05-27):
      The model's water-flux equation is
              Jw  =  Lp · (ΔP  −  σ · Δπ)
      where ΔP is the applied (constant) operating pressure and Δπ is the
      osmotic-pressure difference across the membrane.  When cF rises in a
      CONCENTRATING run, Δπ rises; the net driving force ΔP − σ·Δπ shrinks;
      water flux Jw drops; each successive vial captures less mass per unit
      time.  Mirror image in a DILUTING run: cF falls, Δπ falls, net driving
      force grows, mass per vial increases.

      Plotting ΔP and Δπ side-by-side on the same axes makes this directly
      visible.  No model prediction is required — the curves are built from
      the measured cF (after Shedlovsky inversion + EC25) and a simple van't
      Hoff conversion:  Δπ ≈ (cF − cH) · ni · R · T  ≈  cF · ni · R · T
      (assuming cH << cF, which holds throughout the campaign).

      A bullet annotation calls out the regime (concentrating vs diluting)
      based on whether cF rises or falls over the run, plus the end-of-run
      ratio Δπ_end / ΔP so the user can read off how close the run gets to
      osmotic shutdown.
    """
    settings = results.get("model_settings", {})
    if str(settings.get("workflow_family", "DATA3")).upper() != "DATA3":
        return []

    data_payload = results.get("data")
    if isinstance(data_payload, list):
        data_payload = data_payload[0] if data_payload else None
    if not isinstance(data_payload, dict):
        return []

    fit_meta = results.get("fit_meta") or {}
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data3_option3"
    save_dir.mkdir(parents=True, exist_ok=True)

    sheet_slug = _sanitize_filename_component(data_payload.get("sheet_name", "sheet"))
    file_slug = _sanitize_filename_component(Path(str(data_payload.get("filename", "data3"))).stem)
    prefix = f"{file_slug}_{sheet_slug}"

    cfg = data_payload.get("data_config", {})
    delP = float(cfg.get("delP", 0.0))     # bar
    ni   = float(cfg.get("ni",   1.0))     # van 't Hoff factor (1 for non-dissociating; loader sets to 1 for NF270)
    nc   = int(cfg.get("nc",   1))         # # of salts (NaCl→2 ion species, CaCl₂→3, LaCl₃→4 in terms of ν)
    salt_name = str(cfg.get("namec", ""))
    Temp = float(cfg.get("Temp", 298.15))  # K
    R    = 8.314e-5                        # cm^3 bar / micromol / K  (loader's convention)
    # For van 't Hoff Δπ with full dissociation, use ν (total particles per formula unit).
    # The model multiplies by data_config['ni'] which the loader currently sets to 1, so
    # the physics-correct ν for the diagnostic is the count of ion species per salt:
    nu = {"NaCl": 2, "KCl": 2, "CaCl2": 3, "LaCl3": 4}.get(salt_name, max(int(ni), 1))

    data_raw = data_payload.get("data_raw", []) or []
    if not data_raw:
        return []

    outputs = []
    fig, ax = plt.subplots(figsize=(7.5, 4.5))

    # Applied pressure — constant horizontal line
    ax.axhline(delP, color="tab:orange", linestyle="-", linewidth=2.5,
               label=f"Applied ΔP = {delP:.2f} bar", alpha=0.95)

    # Osmotic pressure curve: Δπ(t) ≈ cF · ν · R · T   (cH ≈ 0 assumption)
    # Concatenate every vial's cF + time into a single trace (loader vials
    # share a continuous time axis already so this draws as one curve).
    all_t = []
    all_dpi = []
    cf_initial = None
    cf_final = None
    for v in data_raw:
        time = np.asarray(v.get("time", []), dtype=float)
        cF = np.asarray(v.get("cF_exp", []), dtype=float)
        m = min(time.size, cF.size)
        if m == 0:
            continue
        all_t.extend(time[:m].tolist())
        all_dpi.extend((cF[:m] * nu * R * Temp).tolist())
        # Track first / last finite cF for regime detection
        finite_cf = cF[:m][np.isfinite(cF[:m])]
        if finite_cf.size:
            if cf_initial is None:
                cf_initial = float(finite_cf[0])
            cf_final = float(finite_cf[-1])
    all_t = np.asarray(all_t, dtype=float)
    all_dpi = np.asarray(all_dpi, dtype=float)
    if all_t.size:
        # Time in seconds (matches the collaborator deck's convention)
        ax.plot(all_t, all_dpi, color="#C71585", linewidth=2.0,
                label="Osmotic ΔΠ(t) = cF · ν · R · T", alpha=0.85)

    # Light reference lines: y=0 (de facto) and a faint band at ΔP-σ·Δπ for σ=1 worst-case
    ax.axhline(0, color="#bbbbbb", linewidth=0.5, alpha=0.5)

    ax.set_xlabel("Time [s]", fontsize=14, fontweight="bold")
    ax.set_ylabel("Pressure [bar]", fontsize=14, fontweight="bold")
    ax.tick_params(direction="in")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=10, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0, frameon=True)

    # Regime detection + headline annotation
    if cf_initial is not None and cf_final is not None:
        delta_cF = cf_final - cf_initial
        regime = "concentrating" if delta_cF > 0 else ("diluting" if delta_cF < 0 else "flat")
        regime_color = {"concentrating": "#B85042", "diluting": "#2C5F2D", "flat": "#6B7A99"}[regime]
        pi_end = (cf_final * nu * R * Temp) if cf_final is not None else float("nan")
        ratio = (pi_end / delP) if (delP > 0 and np.isfinite(pi_end)) else float("nan")
        story_lines = [
            f"Regime: {regime.upper()}  (cF: {cf_initial:.2f} → {cf_final:.2f} mM)",
            f"Δπ_end / ΔP = {ratio:.3f}",
        ]
        if regime == "concentrating":
            story_lines.append("→ ΔΠ rises ⇒ net driving force shrinks")
            story_lines.append("→ mass per vial DECREASES over the run")
        elif regime == "diluting":
            story_lines.append("→ ΔΠ falls ⇒ net driving force grows")
            story_lines.append("→ mass per vial INCREASES over the run")
        ax.text(
            0.02, 0.97, "\n".join(story_lines),
            transform=ax.transAxes, fontsize=9, fontweight="bold",
            verticalalignment="top", horizontalalignment="left", family="monospace",
            color=regime_color,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=regime_color, linewidth=1.5, alpha=0.92),
        )

    # Source / formula footer
    _src_ann = (
        f"Source: measured cF (data file cols 4/3 → Shedlovsky → mM)\n"
        f"Salt: {salt_name},  ν = {nu},  T = {Temp:.1f} K\n"
        f"ΔΠ = cF · ν · R · T   (cH ≈ 0 assumption — bulk Δπ)"
    )
    ax.text(
        0.02, 0.03, _src_ann, transform=ax.transAxes, fontsize=7,
        verticalalignment="bottom", horizontalalignment="left", family="monospace",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", alpha=0.85),
    )
    # fit_meta on the right margin if present
    _fit_lines = _build_fit_meta_annotation(fit_meta, include_cf_unc=False)
    if _fit_lines:
        ax.text(
            1.02, 0.03, "\n".join(_fit_lines),
            transform=ax.transAxes, fontsize=7,
            verticalalignment="bottom", horizontalalignment="left", family="monospace",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", alpha=0.85),
        )

    out_path = save_dir / f"applied_vs_osmotic-{prefix}.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    outputs.append(str(out_path))
    if show:
        plt.show()
    plt.close(fig)
    return outputs


def _resolve_data_root(data_root=None):
    """Find the folder that stores the DATA1 paper inputs."""
    repo_root = Path(__file__).resolve().parents[1]
    if data_root is None:
        env_root = os.environ.get("DIAFILTRATION_DATA1_ROOT")
        if env_root:
            return Path(env_root).expanduser().resolve()
        return repo_root / "legacy" / "data1_matlab" / "data"
    return Path(data_root)


DATA1_RUN_REGISTRY = {
    "501.1": {
        "base": {"data_file": "data_stru-dataset501.1.mat", "subdir": "501.1"},
        "concpolar": {"data_file": "data_stru-dataset501.1.mat", "subdir": "501.1 concpolar"},
    },
    "501.11": {
        "base": {"data_file": "data_stru-dataset501.11.mat", "subdir": "501.11"},
        "concpolar": {"data_file": "data_stru-dataset501.11.mat", "subdir": "501.11 concpolar"},
    },
    "511.11": {
        "base": {"data_file": "data_stru-dataset511.11.mat", "subdir": "511.11"},
        "concpolar": {"data_file": "data_stru-dataset511.11.mat", "subdir": "511.11 concpolar"},
    },
    "511.12": {
        "base": {"data_file": "data_stru-dataset511.12.mat", "subdir": "511.12"},
        "concpolar": {"data_file": "data_stru-dataset511.12.mat", "subdir": "511.12 concpolar"},
    },
}


DATA2_RUN_REGISTRY = {
    "270511.123": {"data_file": "data_stru-dataset270511.123.mat"},
    "270611.121": {"data_file": "data_stru-dataset270611.121.mat"},
    "270711.121": {"data_file": "data_stru-dataset270711.121.mat"},
    "270511.221": {"data_file": "data_stru-dataset270511.221.mat"},
    "270511.321": {"data_file": "data_stru-dataset270511.321.mat"},
    "270511.421": {"data_file": "data_stru-dataset270511.421.mat"},
    "270511.921": {"data_file": "data_stru-dataset270511.921.mat"},
    "270511.521": {"data_file": "data_stru-dataset270511.521.mat"},
    "270511.621": {"data_file": "data_stru-dataset270511.621.mat"},
    "270511.721": {"data_file": "data_stru-dataset270511.721.mat"},
    "270511.821": {"data_file": "data_stru-dataset270511.821.mat"},
}


def get_experimental_run(campaign, run_id, variant="base", data_root=None):
    """Build one file-backed experimental run instance."""
    campaign = str(campaign or "DATA1").upper()
    run_id = str(run_id)
    variant = str(variant or "base").lower()
    root = Path(data_root) if data_root is not None else get_campaign_root(campaign)

    if campaign == "DATA1":
        family = DATA1_RUN_REGISTRY.get(run_id)
        if family is None:
            raise KeyError(f"Unknown DATA1 run '{run_id}'. Known runs: {list(DATA1_RUN_REGISTRY)}")
        variant_info = family.get(variant, family["base"])
        return ExperimentalRun(
            run_id=run_id,
            root=root,
            data_file=variant_info["data_file"],
            variant=variant,
            subdir=variant_info.get("subdir", run_id),
            fit_file="fit_stru.mat",
            contour_files=("contourdata-x_B-y_Lp.csv", "contourdata-x_sigma-y_Lp.csv"),
        )

    if campaign == "DATA2":
        family = DATA2_RUN_REGISTRY.get(run_id)
        if family is None:
            raise KeyError(f"Unknown DATA2 run '{run_id}'. Known runs: {list(DATA2_RUN_REGISTRY)}")
        return ExperimentalRun(
            run_id=run_id,
            root=root,
            data_file=family["data_file"],
            variant=variant,
            subdir="",
            fit_file="fit_stru.mat",
            contour_files=(),
        )

    raise ValueError(f"Unsupported campaign '{campaign}'.")


def _plot_heatmap_frame(df, x_col, y_col, z_col, ax=None, show_title=True, preface=False, cmap="viridis"):
    """Draw one paper-style contour-line panel from a tidy dataframe.

    Matches the DATA1 published figure_s5 / figure_s6 visual idiom:
      * coloured contour lines (no fill) with inline numeric labels at every
        other level, drawn on a clean white background;
      * red triangle at the per-channel minimum of the objective;
      * bold, TeX-formatted axis labels with units (L_p, B, sigma);
      * ticks pointing inward, no colorbar.

    The signature is unchanged so every existing caller of ``plot_contour``
    automatically picks up the readable contour-line style. The ``cmap`` kwarg
    is forwarded to ``ax.contour`` (colours the iso-objective lines).
    """
    if ax is None:
        ax = plt.gca()

    # Pivot the tidy table into a (Y, X) grid for contour rendering.
    grid = df.pivot_table(index=y_col, columns=x_col, values=z_col, aggfunc="mean")
    x_vals = grid.columns.to_numpy(dtype=float)
    y_vals = grid.index.to_numpy(dtype=float)
    z_vals = grid.to_numpy(dtype=float)
    X, Y = np.meshgrid(x_vals, y_vals)

    # Iso-objective contour lines + inline value labels.
    cp = ax.contour(X, Y, z_vals, 10, linewidths=2, cmap=cmap)
    ax.clabel(cp, cp.levels[::2], inline=True, fontsize=10, colors="k", fmt="%1.1f")

    # Red triangle at the argmin of the objective surface.
    if np.isfinite(z_vals).any():
        flat_idx = np.nanargmin(z_vals)
        iy, ix = np.unravel_index(flat_idx, z_vals.shape)
        ax.plot(
            X[iy, ix], Y[iy, ix], "^",
            markersize=12,
            markeredgecolor="red",
            markerfacecolor=[1, 0.6, 0.6],
            clip_on=False,
        )

    # TeX-formatted axis labels with units when the column name matches a
    # known parameter; otherwise fall back to the raw column name.
    xlabel = _PAPER_AXIS_LABELS.get(x_col, x_col)
    ylabel = _PAPER_AXIS_LABELS.get(y_col, y_col)
    ax.set_xlabel(xlabel, fontsize=12, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=12, fontweight="bold")
    if show_title:
        title = z_col if not preface else f"{z_col}"
        ax.set_title(title, fontsize=12)
    ax.tick_params(direction="in", labelsize=10)
    return ax


def plot_contour(df, show_title=True, preface=False, save_path=None, cmap="viridis"):
    """
    Plot the contour-style heatmaps used in the DATA1 paper.

    If the dataframe contains the three objective columns from the paper,
    we draw one subplot for each objective.
    """
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)

    x_col = df.columns[0]
    y_col = df.columns[1]
    value_cols = [
        c
        for c in ["Obj_mass", "Obj_concentration", "Obj_retentate_concentration"]
        if c in df.columns
    ]
    if not value_cols:
        value_cols = [df.columns[2]]

    # Square panels (4×4) — the contour-line style doesn't need horizontal
    # room for a colorbar, so we trade the old wide-panel layout for a
    # readable aspect ratio that matches figure_s5 / figure_s6.
    fig, axes = plt.subplots(1, len(value_cols), figsize=(4 * len(value_cols), 4), squeeze=False)
    for idx, z_col in enumerate(value_cols):
        _plot_heatmap_frame(df, x_col, y_col, z_col, ax=axes[0, idx], show_title=show_title, preface=preface, cmap=cmap)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig, axes


def _plot_contour_data1_legacy(df, output_prefix, save_dir, show_title=False, preface=False):
    """Match the DATA1 notebook contour-panel plotting for main-paper figures."""
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    f_m = df.Obj_mass.values
    ind_m = np.argmin(f_m)
    f_pc = df.Obj_concentration.values
    ind_pc = np.argmin(f_pc)
    f_rc = df.Obj_retentate_concentration.values
    ind_rc = np.argmin(f_rc)

    grid_size = int(round(np.sqrt(len(f_m))))
    if grid_size * grid_size != len(f_m):
        raise ValueError(
            f"Contour dataframe must contain a square grid; got {len(f_m)} rows."
        )

    F_m = np.reshape(f_m, (grid_size, grid_size))
    F_pc = np.reshape(f_pc, (grid_size, grid_size))
    F_rc = np.reshape(f_rc, (grid_size, grid_size))

    if "B" in df:
        xx = df.B.values
        if preface:
            xlabelstr = "B"
            axfontsize = 24
        else:
            xlabelstr = "B [$\\mathbf{\\mu}$m $\\mathbf{\\cdot}$ s$\\mathbf{^{-1}}$]"
            axfontsize = 16
    else:
        xx = df.sigma.values
        if preface:
            xlabelstr = "$\\mathbf{\\sigma}$"
            axfontsize = 24
        else:
            xlabelstr = "$\\mathbf{\\sigma}$ [dimensionless]"
            axfontsize = 16

    if preface:
        ylabelstr = "L$\\mathbf{_p}$"
    else:
        ylabelstr = r"L$\mathbf{_p}$ [L $\mathbf{ \cdot}$ m$\mathbf{^{-2} \cdot}$h$\mathbf{^{-1} \cdot}$bar$\mathbf{^{-1}}$]"

    yy = df.Lp.values
    X = np.reshape(xx, (grid_size, grid_size))
    Y = np.reshape(yy, (grid_size, grid_size))

    panel_specs = [
        ("mass", F_m, ind_m, "Log$\\mathbf{_e}$ transformed \n Mass Objective"),
        ("permeate_conc", F_pc, ind_pc, "Log$\\mathbf{_e}$ transformed \n Permeate Concentration Objective"),
        ("retentate_conc", F_rc, ind_rc, "Log$\\mathbf{_e}$ transformed \n Retentate Concentration Objective"),
    ]

    outputs = []
    for fig_num, (suffix, surface, ind_opt, title) in enumerate(panel_specs, start=1):
        fig = plt.figure(fig_num, figsize=(4, 4))
        plt.clf()
        cp = plt.contour(X, Y, surface, 10, linewidths=2)
        plt.clabel(cp, cp.levels[::2], inline=True, fontsize=12, colors="k", fmt="%1.1f")
        plt.plot(
            xx[ind_opt],
            yy[ind_opt],
            "^",
            markersize=12,
            markeredgecolor="red",
            markerfacecolor=[1, 0.6, 0.6],
            clip_on=False,
        )
        if show_title:
            plt.title(title, fontsize=16, fontweight="bold")
        plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight="bold")
        plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight="bold")
        if preface:
            plt.gca().axes.xaxis.set_ticklabels([])
            plt.gca().axes.yaxis.set_ticklabels([])
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.tick_params(direction="in")
        out_path = save_dir / f"{output_prefix}-{suffix}.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        outputs.append(str(out_path))
        plt.close(fig)

    return outputs


_PAPER_AXIS_LABELS = {
    "Lp":    r"L$\mathbf{_p}$ [L $\mathbf{\cdot}$ m$\mathbf{^{-2}\cdot}$h$\mathbf{^{-1}\cdot}$bar$\mathbf{^{-1}}$]",
    "B":     r"B [$\mathbf{\mu}$m $\mathbf{\cdot}$ s$\mathbf{^{-1}}$]",
    "sigma": r"$\mathbf{\sigma}$ [dimensionless]",
}


def _plot_contour_paper_style(
    df,
    *,
    x_var,
    y_var,
    theta_fit=None,
    output_prefix,
    save_dir,
    show_title=True,
    preface=False,
):
    """DATA1-paper-style contour rendering for an arbitrary (x_var, y_var) pair.

    Produces three separate PNG files (one per data channel), each drawn as a
    proper contour-line plot with labeled iso-objective contours and a red
    triangle at the per-channel minimum. Generalizes
    ``_plot_contour_data1_legacy`` (which assumed y=Lp) to handle the σ-B slice
    where Lp is the held-fixed third parameter.

    Arguments
    ---------
    df : pandas.DataFrame
        Tidy contour dataframe with columns ``[x_var, y_var, Obj_mass,
        Obj_concentration, Obj_retentate_concentration]``.
    x_var, y_var : str
        Names of the swept parameters (must be column names in df).
    theta_fit : dict, optional
        Centering fit; used only to annotate the held parameter in the title.
    output_prefix : str
        Stem prefix for the three output PNGs.
    save_dir : str or Path
        Output directory.

    Returns
    -------
    list of pathlib.Path
        Paths to the three PNGs (one per channel).
    """
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    if x_var not in df.columns or y_var not in df.columns:
        raise ValueError(
            f"DataFrame is missing required columns {x_var!r} and/or {y_var!r}; "
            f"got {list(df.columns)}"
        )

    f_m  = df["Obj_mass"].values
    f_pc = df["Obj_concentration"].values
    f_rc = df["Obj_retentate_concentration"].values

    grid_size = int(round(np.sqrt(len(f_m))))
    if grid_size * grid_size != len(f_m):
        raise ValueError(
            f"Contour dataframe must contain a square grid; got {len(f_m)} rows."
        )

    xx = df[x_var].values
    yy = df[y_var].values
    X = np.reshape(xx, (grid_size, grid_size))
    Y = np.reshape(yy, (grid_size, grid_size))

    F_m  = np.reshape(f_m,  (grid_size, grid_size))
    F_pc = np.reshape(f_pc, (grid_size, grid_size))
    F_rc = np.reshape(f_rc, (grid_size, grid_size))

    ind_m  = int(np.nanargmin(f_m))  if np.any(np.isfinite(f_m))  else None
    ind_pc = int(np.nanargmin(f_pc)) if np.any(np.isfinite(f_pc)) else None
    ind_rc = int(np.nanargmin(f_rc)) if np.any(np.isfinite(f_rc)) else None

    xlabelstr = _PAPER_AXIS_LABELS.get(x_var, x_var)
    ylabelstr = _PAPER_AXIS_LABELS.get(y_var, y_var)
    if preface:
        xlabelstr = {"Lp": r"L$\mathbf{_p}$", "B": "B",
                     "sigma": r"$\mathbf{\sigma}$"}.get(x_var, x_var)
        ylabelstr = {"Lp": r"L$\mathbf{_p}$", "B": "B",
                     "sigma": r"$\mathbf{\sigma}$"}.get(y_var, y_var)
        axfontsize = 24
    else:
        axfontsize = 16

    held_param = next(
        (p for p in ("Lp", "B", "sigma") if p not in (x_var, y_var)),
        None,
    )
    held_value_str = ""
    if held_param is not None and isinstance(theta_fit, dict):
        try:
            held_value_str = f" (held {held_param} = {float(theta_fit[held_param]):.3g})"
        except (KeyError, TypeError, ValueError):
            held_value_str = ""

    panel_specs = [
        ("mass",            F_m,  ind_m,
         r"Log$\mathbf{_{10}}$ transformed" + "\n Mass Objective" + held_value_str),
        ("permeate_conc",   F_pc, ind_pc,
         r"Log$\mathbf{_{10}}$ transformed" + "\n Permeate Concentration Objective" + held_value_str),
        ("retentate_conc",  F_rc, ind_rc,
         r"Log$\mathbf{_{10}}$ transformed" + "\n Retentate Concentration Objective" + held_value_str),
    ]

    outputs = []
    for fig_num, (suffix, surface, ind_opt, title) in enumerate(panel_specs, start=1):
        fig = plt.figure(figsize=(4, 4))
        try:
            cp = plt.contour(X, Y, surface, 10, linewidths=2)
            try:
                plt.clabel(cp, cp.levels[::2], inline=True, fontsize=12,
                           colors="k", fmt="%1.1f")
            except Exception:
                pass
            if ind_opt is not None:
                plt.plot(
                    xx[ind_opt],
                    yy[ind_opt],
                    marker="^",
                    markersize=12,
                    markeredgecolor="red",
                    markerfacecolor=[1, 0.6, 0.6],
                    linestyle="None",
                    clip_on=False,
                )
            if show_title:
                plt.title(title, fontsize=14, fontweight="bold")
            plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight="bold")
            plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight="bold")
            if preface:
                plt.gca().axes.xaxis.set_ticklabels([])
                plt.gca().axes.yaxis.set_ticklabels([])
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.tick_params(direction="in")
            out_path = save_dir / f"{output_prefix}-{suffix}.png"
            fig.savefig(out_path, dpi=300, bbox_inches="tight")
            outputs.append(str(out_path))
        finally:
            plt.close(fig)

    return outputs


def plot_sim(sim_stru, color, linetype):
    """Plot simulation results, matching the DATA1 sigma-sensitivity notebook helper."""

    t_delay = np.asarray(sim_stru[0]["time"], dtype=float)[0]
    # plot mass prediction
    if plt.fignum_exists(1):
        plt.figure(1)
    else:
        plt.figure(1, figsize=(4,4))
    for i in sim_stru:
        t = np.asarray(i["time"], dtype=float)
        plt.plot((t-t_delay)/60, np.asarray(i["mV"], dtype=float), color, linestyle=linetype, linewidth=2, alpha=.8)

    # plot retentate concentration prediction
    if plt.fignum_exists(2):
        plt.figure(2)
    else:
        plt.figure(2, figsize=(4,4))
    for i in sim_stru:
        t = np.asarray(i["time"], dtype=float)
        plt.plot((t-t_delay)/60, np.asarray(i["cF"], dtype=float), color, linestyle=linetype, linewidth=2, alpha=.8)

    # plot permeate concentration prediction comparison
    if plt.fignum_exists(3):
        plt.figure(3)
    else:
        plt.figure(3, figsize=(4,4))
    for i in sim_stru:
        t = np.asarray(i["time"], dtype=float)
        plt.plot((t-t_delay)/60, np.asarray(i["cH"], dtype=float), color, linestyle=linetype, linewidth=2, alpha=.8)


def plot_sim_show(sig, colors):
    """Figure setup for different sigma values."""
    custom_lines = [Line2D([0], [0], color=colors[0],ls='--', lw=3),
                    Line2D([0], [0], color=colors[1],ls='-', lw=3),
                    Line2D([0], [0], color=colors[2],ls=':', lw=3)]
    legendname = ['$\\sigma$ = '+str(sig[0]),'$\\sigma$ = '+str(sig[1]),'$\\sigma$ = '+str(sig[2])]

    fig = plt.figure(1)
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Mass [g]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.ylim(bottom=0)
    fig.savefig('sigma_sensitivity-mass.png',dpi=300,bbox_inches='tight')

    fig = plt.figure(2)
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Retentate Conc. [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.legend(custom_lines,legendname,fontsize=15,loc='best')
    fig.savefig('sigma_sensitivity-reten_conc.png',dpi=300,bbox_inches='tight')

    fig = plt.figure(3)
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Permeate Conc. [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    fig.savefig('sigma_sensitivity-perme_conc.png',dpi=300,bbox_inches='tight')

    return custom_lines, legendname


def plot_conc_range(df_f, df_d, save_path=None, lg=True):
    """
    Plot experiment space - permeate concentration vs. retentate concentration.
    This matches the DATA1 notebook Fig. 3 helper.
    """
    fig = plt.figure(figsize=(4,4))

    # plot filtration
    for i in range(len(df_f.columns)//2):
        xstr = 'F'+str(i+1)+'_cf'
        ystr = 'F'+str(i+1)+'_cp'
        plt.plot(df_f[xstr], df_f[ystr], '^', markersize=8, alpha=.8, clip_on=False)

    # plot diafiltration
    for i in range(len(df_d.columns)//2):
        xstr = 'D'+str(i+1)+'_cf'
        ystr = 'D'+str(i+1)+'_cp'
        plt.plot(df_d[xstr], df_d[ystr], 's', markersize=8, alpha=.8, clip_on=False)

    # ghost point for legend
    plt.plot([],[],'k^',markersize=8,markerfacecolor='white',label='Filtration')
    plt.plot([],[],'ks',markersize=8,markerfacecolor='white',label='Diafiltration')

    plt.xlabel('Retentate [mM]',fontsize=16,fontweight='bold')
    plt.ylabel('Permeate [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    if lg:
        plt.legend(fontsize=15,loc='best')
    if save_path is None:
        save_path = 'concentration_range.png'
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    return fig, plt.gca()


def plot_cr_measure(data_stru, fit_stru, cr_pred, ybottom, lg=False, cond=True, save_path=None):
    """Plot retentate concentration vs. time from measurements and mass balance."""

    t_delay = _data1_time_origin(data_stru)

    def _fit_first_cF(fit_bundle):
        sim = fit_bundle.get("sim_stru") if isinstance(fit_bundle, dict) else fit_bundle
        if isinstance(sim, dict):
            sim = sim.get(0, sim.get("0", sim))
        return np.asarray(sim[0]["cF"], dtype=float).reshape(-1)[0]

    def _as_scalar_or_array(values):
        arr = np.asarray(values, dtype=float)
        if arr.ndim == 0:
            return float(arr.reshape(-1)[0])
        return arr

    fig = plt.figure(figsize=(4, 4))

    # initial point / starting concentration
    if cond:
        plt.plot(
            (np.asarray(data_stru["data_raw"][0]["time"], dtype=float)[0] - t_delay) / 60,
            _fit_first_cF(fit_stru),
            "ms",
            markersize=8,
            clip_on=False,
        )
    else:
        plt.plot(
            (np.asarray(data_stru["data_raw"][0]["time"], dtype=float)[0] - t_delay) / 60,
            data_stru["data_config"]["C_F0"],
            "go",
            markersize=8,
            clip_on=False,
        )

    # measurement points from conductivity probe / vial sampling
    for i in range(data_stru["data_config"]["n"]):
        row = data_stru["data_raw"][i]
        if cond:
            cF_exp = row["cF_exp"]
            if isinstance(cF_exp, (float, np.floating)):
                plt.plot(
                    (np.asarray(row["time"], dtype=float)[-1] - t_delay) / 60,
                    cF_exp,
                    "ms",
                    markersize=8,
                )
            elif len(cF_exp) > 1:
                plt.plot(
                    (np.asarray(row["time"], dtype=float) - t_delay) / 60,
                    np.asarray(cF_exp, dtype=float),
                    "ms",
                    markersize=8,
                )

    # calculated retentate concentration from mass balance
    if len(cr_pred) > 0:
        cr_pred = np.asarray(cr_pred, dtype=float).reshape(-1)
        for i in range(data_stru["data_config"]["n"]):
            plt.plot(
                (np.asarray(data_stru["data_raw"][i]["time"], dtype=float)[-1] - t_delay) / 60,
                cr_pred[i],
                "go",
                markersize=8,
            )

    plt.xlabel("Time [min]", fontsize=16, fontweight="bold")
    plt.ylabel("Concentration [mM]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in", top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=ybottom)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    if xticks:
        xticks[0].label1.set_visible(False)

    plt.plot([], [], "ms", markersize=8, label="Retentate \n (Measurements)")
    plt.plot([], [], "go", markersize=8, label="Retentate \n (Calculated)")
    if lg:
        plt.legend(fontsize=15, loc="best")

    if save_path is None:
        save_path = f"cr_measure-dat{data_stru['dataset']}.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig, ax


def plot_contour_sig_sen(
    data_stru,
    contour_sig_stru,
    Vmin,
    Vmax,
    Level,
    Colorbar_ticks,
    Manual_locations,
    name_append,
    filled=False,
    Jw_filter=True,
    bar=False,
    save_dir=None,
):
    """Plot the DATA1 sigma-sensitivity contour figures."""
    diaf_mode = False
    if "cf0" in contour_sig_stru:
        x = contour_sig_stru["cf0"]
        xlabelstr = r"$\mathbf{c_f}$(t=0) [mM]"
    else:
        x = contour_sig_stru["cd"]
        xlabelstr = r"$\mathbf{c_d}$ [mM]"
        diaf_mode = True
    y = contour_sig_stru["delp"]
    ylabelstr = r"$\mathbf{\Delta}$P [psi]"

    X, Y = np.meshgrid(x, y)
    axfontsize = 16

    R = 8.314e-5
    T = data_stru["data_config"]["Temp"]
    ni = data_stru["data_config"]["ni"]

    if not diaf_mode:
        Jw0 = Y / 14.504 - 1 * ni * R * T * X
        jw0 = 1 * ni * R * T * np.array(x) * 14.504
    else:
        cf0 = 5
        Jw0 = Y / 14.504 - 1 * ni * R * T * cf0
        jw0 = 1 * ni * R * T * cf0 * 14.504

    def patch_back():
        ax = plt.gca()
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        p = patches.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, color="gray", zorder=-10)
        ax.add_patch(p)

    fig = plt.figure(1, figsize=(4, 4))
    if filled:
        plt.contourf(X, Y, contour_sig_stru["range_mV"], 50, cmap="spring", vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X, Y, contour_sig_stru["range_mV"], Level[0], linewidths=3, colors="k")
    else:
        cp = plt.contour(X, Y, contour_sig_stru["range_mV"], Level[0], linewidths=2)
    try:
        plt.clabel(cp, inline=True, manual=Manual_locations[0], fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    except Exception:
        plt.clabel(cp, inline=True, fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color="gray", zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors="orangered")
        plt.setp(cg.collections, path_effects=[patheffects.withTickedStroke(angle=300, length=2)])
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight="bold")
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight="bold")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom=y[0])
    out_mass = Path(save_dir) / f"contour_sensitivity{name_append}-mass.png" if save_dir else Path(f"contour_sensitivity{name_append}-mass.png")
    fig.savefig(out_mass, dpi=300, bbox_inches="tight")

    fig = plt.figure(2, figsize=(4, 4))
    if filled:
        plt.contourf(X, Y, contour_sig_stru["range_cF"], 100, cmap="spring", vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X, Y, contour_sig_stru["range_cF"], Level[1], linewidths=3, colors="k")
    else:
        cp = plt.contour(X, Y, contour_sig_stru["range_cF"], Level[1], linewidths=2)
    try:
        plt.clabel(cp, inline=True, manual=Manual_locations[1], fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    except Exception:
        plt.clabel(cp, inline=True, fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color="gray", zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors="orangered")
        plt.setp(cg.collections, path_effects=[patheffects.withTickedStroke(angle=300, length=2)])
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight="bold")
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight="bold")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom=y[0])
    out_retentate = Path(save_dir) / f"contour_sensitivity{name_append}-retentate_conc.png" if save_dir else Path(f"contour_sensitivity{name_append}-retentate_conc.png")
    fig.savefig(out_retentate, dpi=300, bbox_inches="tight")

    fig = plt.figure(3, figsize=(4, 4))
    if filled:
        plt.contourf(X, Y, contour_sig_stru["range_cH"], 50, cmap="spring", vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X, Y, contour_sig_stru["range_cH"], Level[2], linewidths=3, colors="k")
    else:
        cp = plt.contour(X, Y, contour_sig_stru["range_cH"], Level[2], linewidths=2)
    try:
        plt.clabel(cp, inline=True, manual=Manual_locations[2], fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    except Exception:
        plt.clabel(cp, inline=True, fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color="gray", zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors="orangered")
        plt.setp(cg.collections, path_effects=[patheffects.withTickedStroke(angle=300, length=2)])
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight="bold")
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight="bold")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom=y[0])
    out_permeate = Path(save_dir) / f"contour_sensitivity{name_append}-permeate_conc.png" if save_dir else Path(f"contour_sensitivity{name_append}-permeate_conc.png")
    fig.savefig(out_permeate, dpi=300, bbox_inches="tight")

    if filled and bar:
        cmap = mpl.cm.spring
        norm = mpl.colors.Normalize(vmin=Vmin, vmax=Vmax)
        fig, ax = plt.subplots(figsize=(12, 1))
        fig.subplots_adjust(bottom=0.5)
        mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation="horizontal", ticks=Colorbar_ticks, extend="max")
        ax.tick_params(labelsize=15)
        out_bar_h = Path(save_dir) / f"colorbar{name_append}-horizontal.png" if save_dir else Path(f"colorbar{name_append}-horizontal.png")
        fig.savefig(out_bar_h, dpi=300, bbox_inches="tight")
        fig, ax = plt.subplots(figsize=(1, 4))
        fig.subplots_adjust(left=0.5)
        mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation="vertical", ticks=Colorbar_ticks, extend="max")
        ax.tick_params(labelsize=15)
        out_bar_v = Path(save_dir) / f"colorbar{name_append}-vertical.png" if save_dir else Path(f"colorbar{name_append}-vertical.png")
        fig.savefig(out_bar_v, dpi=300, bbox_inches="tight")

    return [str(p) for p in [out_mass, out_retentate, out_permeate] + ([out_bar_h, out_bar_v] if filled and bar else [])]


def calib_curve_cond(calib_curve, save_path=None):
    """Plot the conductivity calibration curve used in the DATA1 notebook."""
    x = calib_curve.Conductivity.values
    y = calib_curve.Concentration.values

    fig = plt.figure(figsize=(4, 4))
    plt.plot(x, y, "bo", markersize=8)
    z = np.polyfit(x, y, 1)
    l = np.poly1d(z)
    correlation = np.corrcoef(x, y)[0, 1]
    r_squared = correlation**2
    plt.plot(x, l(x), "b:", linewidth=3, alpha=.7)
    idx = min(4, len(x) - 1)
    eqn = "y=%.3fx+%.2f \n R$\\mathbf{^{2}}$=%.4f" % (z[0], z[1], r_squared)
    plt.annotate(
        eqn,
        xy=(x[idx], y[idx]),
        xycoords="data",
        xytext=(-30, 60),
        weight="bold",
        textcoords="offset points",
        size=12,
        ha="center",
        va="center",
        bbox=dict(boxstyle="round", color="b", alpha=0.1),
    )
    xlabelstr = "Conductivity [$\\mathbf{\\mu}$S $\\mathbf{\\cdot}$ cm$\\mathbf{^{-1}}$]"
    ylabelstr = "Concentration [mM]"
    plt.xlabel(xlabelstr, fontsize=16, fontweight="bold")
    plt.ylabel(ylabelstr, fontsize=16, fontweight="bold")
    plt.xticks(fontsize=15, rotation=45)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in", top=True, right=True)
    if save_path is None:
        save_path = "calib_curve.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def _plot_data2_calibration_panel(calib_curve, z, save_path, panel_label=None):
    """Notebook-faithful DATA2 conductivity-to-concentration calibration panel."""
    x = calib_curve.Conductivity.values
    y = calib_curve.Concentration.values

    fig = plt.figure(figsize=(4, 4))
    plt.plot(x, y, "bo", markersize=8)
    trend = np.poly1d(z)
    correlation = np.corrcoef(x, y)[0, 1]
    r2 = r2_score(y, trend(x))
    plt.plot(x, trend(x), "b:", linewidth=3, alpha=0.7)

    eqn = "y=%.5fx%.2f\nR$\\mathbf{^{2}}$=%.4f" % (z[0], z[1], r2)
    ax = plt.gca()
    ax.text(
        0.28,
        0.87,
        eqn,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        ha="center",
        va="center",
        bbox=dict(boxstyle="round", color="b", alpha=0.1),
    )

    xlabelstr = "Conductivity [$\\mathbf{\\mu}$S $\\mathbf{\\cdot}$ cm$\\mathbf{^{-1}}$]"
    ylabelstr = "Concentration [mM]"
    plt.xlabel(xlabelstr, fontsize=16, fontweight="bold")
    plt.ylabel(ylabelstr, fontsize=16, fontweight="bold")
    plt.xticks(fontsize=15, rotation=45)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in", top=True, right=True)

    if panel_label:
        ax.text(-0.16, 1.02, panel_label, transform=ax.transAxes, fontsize=22, fontweight="bold")

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def _compose_data2_figure_s1(panel_a, panel_b, out_path):
    """Compose DATA2 SI Figure S1 from two calibration panels."""
    try:
        from PIL import Image
    except Exception:
        return None

    panels = [Path(panel_a), Path(panel_b)]
    if not all(path.exists() for path in panels):
        return None

    imgs = [Image.open(path).convert("RGB") for path in panels]
    target_h = max(im.height for im in imgs)
    resized = []
    for im in imgs:
        scale = target_h / im.height
        resized.append(im.resize((int(round(im.width * scale)), target_h), Image.Resampling.LANCZOS))

    margin = 20
    gap = 24
    page_w = margin * 2 + sum(im.width for im in resized) + gap
    page_h = margin * 2 + target_h
    page = Image.new("RGB", (page_w, page_h), "white")
    x = margin
    for im in resized:
        page.paste(im, (x, margin))
        x += im.width + gap

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    page.save(out_path)
    return out_path


def run_data1_si_s2(data_root=None, save_dir=None):
    """Recreate DATA1 SI Figure S2 from the notebook's three concentration-ratio plots."""
    root = _resolve_data1_figure2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "si_fig_s2"
    save_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        {
            "label": "A",
            "dataset": "301.1",
            "data_path": root / "dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0" / "data_stru-dataset301.1.mat",
            "fit_path": root / "dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0" / "fit_stru-dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0.mat",
            "cr_path": root / "experiment space" / "Classical_analysis-dat301.1.csv",
            "cr_col": "cf",
            "ybottom": 4,
            "cond": False,
            "lg": False,
        },
        {
            "label": "B",
            "dataset": "501.1",
            "data_path": root / "data_stru-dataset501.1.mat",
            "fit_path": root / "501.1 concpolar" / "fit_stru.mat",
            "cr_path": root / "experiment space" / "Classical_analysis-dat501.1.csv",
            "cr_col": "cf",
            "ybottom": 4,
            "cond": True,
            "lg": False,
        },
        {
            "label": "C",
            "dataset": "511.12",
            "data_path": root / "data_stru-dataset511.12.mat",
            "fit_path": root / "511.12 concpolar" / "fit_stru.mat",
            "cr_path": root / "experiment space" / "diafiltration.csv",
            "cr_col": "D3_cf",
            "ybottom": 0,
            "cond": True,
            "lg": True,
        },
    ]

    panel_paths = []
    for case in cases:
        data_stru = loadmat(str(case["data_path"]))["data_stru"]
        fit_bundle = loadmat(str(case["fit_path"]))
        fit_stru = fit_bundle.get("fit_stru", fit_bundle)
        if case["dataset"] == "511.12":
            df = pd.read_csv(case["cr_path"], header=2)
        else:
            df = pd.read_csv(case["cr_path"], header=0)
        cr_pred = df[case["cr_col"]]
        out_path = save_dir / f"cr_measure-dat{case['dataset']}.png"
        fig, _ = plot_cr_measure(
            data_stru,
            fit_stru,
            cr_pred,
            case["ybottom"],
            lg=case["lg"],
            cond=case["cond"],
            save_path=out_path,
        )
        plt.close(fig)
        panel_paths.append(out_path)

    composite = _compose_data1_panel_grid(
        [[panel_paths[0], panel_paths[1], panel_paths[2]]] if len(panel_paths) == 3 else [],
        save_dir / "figure_s2.png",
        row_labels=["A", "B", "C"],
        label_mode="per_panel",
        panel_margin=30,
        row_margin=30,
        outer_margin=24,
    )
    if composite is not None:
        panel_paths.append(composite)

    return [str(p) for p in panel_paths]


def _compose_data1_panel_grid(
    row_paths,
    out_path,
    row_labels=None,
    label_mode="per_row",
    panel_margin=10,
    row_margin=14,
    outer_margin=20,
):
    """Compose notebook-style DATA1 panel sheets from existing image panels."""
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return None

    rows = []
    for row in row_paths:
        valid = [Path(p) for p in row if Path(p).exists()]
        if not valid:
            continue
        imgs = [Image.open(path).convert("RGB") for path in valid]
        target_h = max(im.height for im in imgs)
        resized = []
        for im in imgs:
            scale = target_h / im.height
            resized.append(im.resize((int(im.width * scale), target_h), Image.Resampling.LANCZOS))
        row_w = sum(im.width for im in resized) + panel_margin * max(0, len(resized) - 1)
        row_h = target_h
        rows.append((resized, row_w, row_h))

    if not rows:
        return None

    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 28)
    except Exception:
        font = None

    label_pad = 36 if row_labels else 0
    width = max(row_w for _, row_w, _ in rows) + outer_margin * 2 + label_pad
    height = sum(row_h for _, _, row_h in rows) + row_margin * max(0, len(rows) - 1) + outer_margin * 2
    page = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(page)

    y = outer_margin
    label_iter = iter(row_labels or [])
    for resized, _, row_h in rows:
        x = outer_margin + label_pad
        panel_labels = []
        if label_mode == "per_row":
            panel_labels = [next(label_iter, "")]
        elif label_mode == "per_panel":
            panel_labels = [next(label_iter, "") for _ in resized]
        if label_mode == "per_row" and panel_labels:
            draw.text((outer_margin, y + 4), panel_labels[0], fill="black", font=font)
        for idx, im in enumerate(resized):
            if label_mode == "per_panel" and idx < len(panel_labels):
                draw.text((x, y - 28), panel_labels[idx], fill="black", font=font)
            page.paste(im, (x, y))
            x += im.width + panel_margin
        y += row_h + row_margin

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    page.save(out_path)
    return str(out_path)


def _compose_data1_panel_sheet(row_paths, out_path):
    """Stack existing DATA1 panel images into one clean composite sheet."""
    try:
        from PIL import Image
    except Exception:
        return None

    valid_rows = [Path(p) for p in row_paths if Path(p).exists()]
    if not valid_rows:
        return None

    imgs = [Image.open(p).convert("RGB") for p in valid_rows]
    target_w = max(im.width for im in imgs)
    resized = []
    for im in imgs:
        scale = target_w / im.width
        resized.append(
            im.resize((target_w, int(im.height * scale)), Image.Resampling.LANCZOS)
        )

    gap = 20
    canvas_w = target_w + 20
    canvas_h = sum(im.height for im in resized) + gap * (len(resized) - 1) + 20
    page = Image.new("RGB", (canvas_w, canvas_h), "white")

    y = 10
    for im in resized:
        page.paste(im, (10, y))
        y += im.height + gap

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    page.save(out_path)
    return str(out_path)


def _run_data1_si_contour_figure(data_root, save_dir, figure_name, cases, contour_filename, prefix_base):
    """Build one DATA1 SI contour page from notebook-style contour CSVs."""
    root = _resolve_data1_figure2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data1_paper_figures"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    row_groups = []
    row_labels = []
    for label, folder in cases:
        csv_path = root / folder / contour_filename
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        prefix = f"{prefix_base}_{label.lower()}"
        panel_paths = _plot_contour_data1_legacy(df, prefix, save_dir, show_title=False, preface=False)
        outputs.extend(panel_paths)
        row_groups.append([Path(path) for path in panel_paths])
        row_labels.append(label)

    composite = _compose_data1_panel_grid(
        row_groups,
        save_dir / figure_name,
        row_labels=row_labels,
        label_mode="per_row",
        panel_margin=12,
        row_margin=18,
        outer_margin=22,
    )
    if composite is not None:
        outputs.append(composite)
    return outputs


def run_data1_si_s3(data_root=None, save_dir=None):
    """Recreate DATA1 SI Figure S3 from notebook diafiltration sigma contours."""
    cases = [
        ("A", "511.12 concpolar"),
        ("B", "511.11 concpolar"),
        ("C", "511.12"),
        ("D", "511.11"),
    ]
    return _run_data1_si_contour_figure(data_root, save_dir, "figure_s3.png", cases, "contourdata-x_sigma-y_Lp.csv", "figure_s3_panel")


def run_data1_si_s4(data_root=None, save_dir=None):
    """Recreate DATA1 SI Figure S4 from notebook diafiltration B contours."""
    cases = [
        ("A", "511.12 concpolar"),
        ("B", "511.11 concpolar"),
        ("C", "511.12"),
        ("D", "511.11"),
    ]
    return _run_data1_si_contour_figure(data_root, save_dir, "figure_s4.png", cases, "contourdata-x_B-y_Lp.csv", "figure_s4_panel")


def run_data1_si_s5(data_root=None, save_dir=None):
    """Recreate DATA1 SI Figure S5 from notebook filtration sigma contours."""
    cases = [
        ("A", "501.1 concpolar"),
        ("B", "501.11 concpolar"),
        ("C", "501.1"),
        ("D", "501.11"),
    ]
    return _run_data1_si_contour_figure(data_root, save_dir, "figure_s5.png", cases, "contourdata-x_sigma-y_Lp.csv", "figure_s5_panel")


def run_data1_si_s6(data_root=None, save_dir=None):
    """Recreate DATA1 SI Figure S6 from notebook filtration B contours."""
    cases = [
        ("A", "501.1 concpolar"),
        ("B", "501.11 concpolar"),
        ("C", "501.1"),
        ("D", "501.11"),
    ]
    return _run_data1_si_contour_figure(data_root, save_dir, "figure_s6.png", cases, "contourdata-x_B-y_Lp.csv", "figure_s6_panel")


def run_data1_si_panel_sheets(save_dir=None):
    """Build DATA1 SI panel-only composite sheets from generated panel PNGs."""
    save_dir = Path(save_dir) if save_dir is not None else Path.cwd()
    outputs = []

    composites = {
        "data1_si_reduced_diafiltration.png": [
            save_dir / "data511.11_fit.png",
            save_dir / "data511.11_fit_concentration.png",
        ],
        "data1_si_diafiltration_B.png": [
            save_dir / "data511.12_concpolar_contour_B.png",
            save_dir / "data511.11_concpolar_contour_B.png",
            save_dir / "data511.12_base_contour_B.png",
            save_dir / "data511.11_base_contour_B.png",
        ],
        "data1_si_reduced_filtration.png": [
            save_dir / "data501.11_fit.png",
            save_dir / "data501.11_fit_concentration.png",
        ],
        "data1_si_filtration_sigma.png": [
            save_dir / "data501.1_concpolar_contour_sigma.png",
            save_dir / "data501.11_concpolar_contour_sigma.png",
            save_dir / "data501.1_base_contour_sigma.png",
            save_dir / "data501.11_base_contour_sigma.png",
        ],
        "data1_si_filtration_B.png": [
            save_dir / "data501.1_concpolar_contour_B.png",
            save_dir / "data501.11_concpolar_contour_B.png",
            save_dir / "data501.1_base_contour_B.png",
            save_dir / "data501.11_base_contour_B.png",
        ],
    }

    for name, rows in composites.items():
        result = _compose_data1_panel_sheet(rows, save_dir / name)
        if result is not None:
            outputs.append(result)

    return outputs


def plot_conc_comparison(data_stru_f, fit_stru_f, data_stru_d, fit_stru_d, plot_pred=True, save_path=None):
    """Plot filtration and diafiltration concentration comparison figures."""
    t_delay_f = _data1_time_origin(data_stru_f)
    t_delay_d = _data1_time_origin(data_stru_d)

    def _first_conc_value(data_stru, key, fallback_key="cF_exp"):
        if key in data_stru and len(data_stru[key]) > 0:
            val = data_stru[key][0]
            if isinstance(val, (list, tuple, np.ndarray)):
                return _last_valid_value(val)
            return val
        for row in data_stru.get("data_raw", []):
            if fallback_key in row:
                val = row[fallback_key]
                if isinstance(val, (list, tuple, np.ndarray)):
                    return _last_valid_value(val)
                return val
        return np.nan

    def _plot_measurement_series(data_stru, i, x_shift, marker, color, cV_key="cV_avg", cF_key="cF_exp"):
        row = data_stru["data_raw"][i]
        time = np.asarray(row["time"], dtype=float) - x_shift

        def _plot_one(values, style):
            arr = np.asarray(values, dtype=float)
            if arr.ndim == 0 or arr.size == 1:
                x = np.array([time[-1]], dtype=float)
                y = np.array([float(arr.reshape(-1)[0])], dtype=float)
            else:
                n = min(len(time), len(arr))
                x = time[:n]
                y = arr[:n]
            mask = ~np.isnan(y)
            if np.any(mask):
                plt.plot(x[mask], y[mask], style, markersize=6)

        if cV_key in row:
            _plot_one(row[cV_key], marker.replace("s", "o"))
        if cF_key in row:
            _plot_one(row[cF_key], marker)

    fig = plt.figure(figsize=(4, 4))
    t0_f = _data1_time_origin(data_stru_f) - t_delay_f
    t0_d = _data1_time_origin(data_stru_d) - t_delay_d
    plt.plot(t0_f, _first_conc_value(data_stru_f, "cF_ICP"), "mv", markersize=6, clip_on=False)
    plt.plot(t0_d, _first_conc_value(data_stru_d, "cF_ICP"), "k^", markersize=6, clip_on=False)
    plt.plot(t0_f, float(np.asarray(fit_stru_f["sim_stru"][0]["cF"], dtype=float).reshape(-1)[0]), "bs", markersize=6, clip_on=False)
    plt.plot(t0_d, float(np.asarray(fit_stru_d["sim_stru"][0]["cF"], dtype=float).reshape(-1)[0]), "rs", markersize=6, clip_on=False)

    plt.plot([], [], "kv", markersize=6, clip_on=False, label="Filtration Retentate ICP")
    plt.plot([], [], "bs", markersize=6, clip_on=False, label="Filtration Retentate conductivity")
    plt.plot([], [], "bo", markersize=6, label="Filtration Vial ICP")
    plt.plot([], [], "k^", markersize=6, clip_on=False, label="Diafiltration Retentate ICP")
    plt.plot([], [], "rs", markersize=6, clip_on=False, label="Diafiltration Retentate conductivity")
    plt.plot([], [], "ro", markersize=6, label="Diafiltration Vial ICP")

    if plot_pred:
        plt.plot([], [], "g-.", linewidth=2, label="Filtration Retentate")
        plt.plot([], [], "g", linewidth=2, label="Diafiltration Retentate")
        plt.plot([], [], "k-.", linewidth=2, alpha=.6, label="Filtration Permeate")
        plt.plot([], [], "k", linewidth=2, alpha=.6, label="Diafiltration Permeate")

    for i in range(data_stru_f["data_config"]["n"]):
        _plot_measurement_series(data_stru_f, i, t_delay_f, "bs", "b")
        if plot_pred:
            t_pred = np.asarray(fit_stru_f["sim_stru"][i]["time"], dtype=float) - t_delay_f
            plt.plot(t_pred, np.asarray(fit_stru_f["sim_stru"][i]["cF"], dtype=float), "g-.", linewidth=2)
            plt.plot(t_pred, np.asarray(fit_stru_f["sim_stru"][i]["cH"], dtype=float), "k-.", linewidth=2, alpha=.6)

    for i in range(data_stru_d["data_config"]["n"]):
        _plot_measurement_series(data_stru_d, i, t_delay_d, "rs", "r")
        if plot_pred:
            t_pred = np.asarray(fit_stru_d["sim_stru"][i]["time"], dtype=float) - t_delay_d
            plt.plot(t_pred, np.asarray(fit_stru_d["sim_stru"][i]["cF"], dtype=float), "g", linewidth=2)
            plt.plot(t_pred, np.asarray(fit_stru_d["sim_stru"][i]["cH"], dtype=float), "k", linewidth=2, alpha=.6)

    plt.plot(data_stru_f["data_raw"][-1]["time"][-1] - t_delay_f, _first_conc_value(data_stru_f, "cF_ICP"), "kv", markersize=6, clip_on=False)
    plt.plot(data_stru_d["data_raw"][-1]["time"][-1] - t_delay_d, _first_conc_value(data_stru_d, "cF_ICP"), "k^", markersize=6, clip_on=False)
    plt.xlabel("Time [s]", fontsize=16)
    plt.ylabel("Concentration [mM]", fontsize=16)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.minorticks_on()
    plt.tick_params(direction="in", top=True, right=True)
    plt.tick_params(which="minor", direction="in", top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(fontsize=12, bbox_to_anchor=(1.1, 1.05), loc="upper left")
    if save_path is None:
        save_path = "concentration.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def run_data_analysis(data_root=None, datasets=None, save_dir=None):
    """
    Recreate the paper-style DATA1 analysis plots for one or more datasets.

    This is a simple helper that loads the saved paper files and calls the
    plotting functions above.
    """
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root
    save_dir.mkdir(parents=True, exist_ok=True)

    if datasets is None:
        datasets = [501.1, 501.11, 511.12, 511.11]

    outputs = []
    for dat in datasets:
        try:
            data_run = get_experimental_run("DATA1", str(dat), variant="base", data_root=root)
        except Exception:
            continue
        if not data_run.data_path.exists():
            continue
        data_stru = data_run.load_data()

        # Plot the main fit if the saved fit file is available.
        for variant_name in ["base", "concpolar"]:
            try:
                variant_run = get_experimental_run("DATA1", str(dat), variant=variant_name, data_root=root)
            except Exception:
                continue
            if variant_run.fit_path.exists():
                fit_bundle = variant_run.load_fit()
                sim_stru = fit_bundle.get("sim_stru") if isinstance(fit_bundle, dict) else None
                if sim_stru is None:
                    continue
                figs = plot_sim_comparison(data_stru, sim_stru, plot_pred=True, lg=False)
                if figs:
                    out = save_dir / f"data{dat}_fit.png"
                    figs[0].savefig(out, dpi=300, bbox_inches="tight")
                    outputs.append(str(out))
                    plt.close(figs[0])
                if len(figs) > 1:
                    out = save_dir / f"data{dat}_fit_concentration.png"
                    figs[1].savefig(out, dpi=300, bbox_inches="tight")
                    outputs.append(str(out))
                    plt.close(figs[1])
                if len(figs) > 2:
                    out = save_dir / f"data{dat}_fit_stirred.png"
                    figs[2].savefig(out, dpi=300, bbox_inches="tight")
                    outputs.append(str(out))
                    plt.close(figs[2])

            contour_b = variant_run.contour_path("contourdata-x_B-y_Lp.csv")
            if contour_b.exists():
                df = variant_run.load_contour("contourdata-x_B-y_Lp.csv")
                fig, _ = plot_contour(df)
                out = save_dir / f"data{dat}_{variant_name}_contour_B.png"
                fig.savefig(out, dpi=300, bbox_inches="tight")
                outputs.append(str(out))
                plt.close(fig)

            contour_s = variant_run.contour_path("contourdata-x_sigma-y_Lp.csv")
            if contour_s.exists():
                df = variant_run.load_contour("contourdata-x_sigma-y_Lp.csv")
                fig, _ = plot_contour(df)
                out = save_dir / f"data{dat}_{variant_name}_contour_sigma.png"
                fig.savefig(out, dpi=300, bbox_inches="tight")
                outputs.append(str(out))
                plt.close(fig)

        filtration_csv = root / "experiment space" / "filtration.csv"
        diafiltration_csv = root / "experiment space" / "diafiltration.csv"
        if filtration_csv.exists() and diafiltration_csv.exists():
            df_f = pd.read_csv(filtration_csv, header=2)
            df_d = pd.read_csv(diafiltration_csv, header=2)
            fig, _ = plot_conc_range(df_f, df_d)
            out = save_dir / f"data{dat}_conc_range.png"
            fig.savefig(out, dpi=300, bbox_inches="tight")
            outputs.append(str(out))
            plt.close(fig)

    calib_csv = root / "experiment space" / "conductivity_calibration.csv"
    if calib_csv.exists():
        calib_curve = pd.read_csv(calib_csv, header=0, skiprows=[1])
        fig = calib_curve_cond(calib_curve, save_path=save_dir / "calib_curve.png")
        outputs.append(str(save_dir / "calib_curve.png"))
        plt.close(fig)

    return outputs


def run_sigma_sensitivity(data_root=None, dataset=501.1, cf0=None, sigmas=None, save_dir=None, output_prefix=None):
    """
    Recreate the sigma-sensitivity plots from the DATA1 paper.

    The function reads the saved MAT files and overlays the simulations
    using the simple plotting helpers above.
    """
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "sigma_sensitivity"
    save_dir.mkdir(parents=True, exist_ok=True)

    if sigmas is None:
        sigmas = [0.1, 0.5, 0.9]

    if cf0 is None:
        cf0 = 5.2843 if float(dataset) == 501.1 else 15.2052

    colorstring = "rbg"
    linetypes = ["dashed", "solid", "dotted"]

    plt.close(1)
    plt.close(2)
    plt.close(3)

    for sigma_val, color, linetype in zip(sigmas, colorstring, linetypes):
        mat_name = f"sim_stru-dat{dataset} C_Fin{cf0}sig{sigma_val}.mat"
        mat_path = root / "sigma sensitivity" / mat_name
        if not mat_path.exists():
            continue
        sim_stru = loadmat(str(mat_path)).get("sim_stru")
        plot_sim(sim_stru, color, linetype)

    plot_sim_show(sigmas, colorstring[: len(sigmas)])
    prefix = output_prefix or "sigma_sensitivity"
    out_mass = save_dir / f"{prefix}-mass.png"
    out_reten = save_dir / f"{prefix}-reten_conc.png"
    out_perm = save_dir / f"{prefix}-perme_conc.png"
    plt.figure(1).savefig(out_mass, dpi=300, bbox_inches="tight")
    plt.figure(2).savefig(out_reten, dpi=300, bbox_inches="tight")
    plt.figure(3).savefig(out_perm, dpi=300, bbox_inches="tight")
    plt.close(1)
    plt.close(2)
    plt.close(3)
    return [str(out_mass), str(out_reten), str(out_perm)]


def run_data1_figure4_workflow(data_root=None, save_dir=None):
    """Recreate DATA1 Figure 4 from the notebook sigma-sensitivity panels."""
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data1_paper_figures"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    outputs.extend(
        run_sigma_sensitivity(
            data_root=root,
            dataset=501.1,
            save_dir=save_dir,
            output_prefix="sigma_sensitivity_filtration",
        )
    )
    outputs.extend(
        run_sigma_sensitivity(
            data_root=root,
            dataset=511.12,
            save_dir=save_dir,
            output_prefix="sigma_sensitivity_diafiltration",
        )
    )

    try:
        from PIL import Image, ImageDraw, ImageFont

        ordered = [
            [
                save_dir / "sigma_sensitivity_filtration-mass.png",
                save_dir / "sigma_sensitivity_filtration-perme_conc.png",
                save_dir / "sigma_sensitivity_filtration-reten_conc.png",
            ],
            [
                save_dir / "sigma_sensitivity_diafiltration-mass.png",
                save_dir / "sigma_sensitivity_diafiltration-perme_conc.png",
                save_dir / "sigma_sensitivity_diafiltration-reten_conc.png",
            ],
        ]
        if all(path.exists() for row in ordered for path in row):
            rows = []
            for row in ordered:
                imgs = [Image.open(path).convert("RGB") for path in row]
                target_h = max(im.height for im in imgs)
                resized = []
                for im in imgs:
                    scale = target_h / im.height
                    resized.append(im.resize((int(im.width * scale), target_h), Image.Resampling.LANCZOS))
                row_w = sum(im.width for im in resized) + 40
                row_h = target_h + 40
                canvas = Image.new("RGB", (row_w, row_h), "white")
                x = 20
                for im in resized:
                    canvas.paste(im, (x, 20))
                    x += im.width + 10
                rows.append(canvas)

            width = max(im.width for im in rows) + 20
            height = sum(im.height for im in rows) + 30
            page = Image.new("RGB", (width, height), "white")
            draw = ImageDraw.Draw(page)
            try:
                font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 28)
            except Exception:
                font = None
            y = 10
            for label, row in zip(["A", "B"], rows):
                draw.text((10, y + 5), label, fill="black", font=font)
                page.paste(row, (30, y))
                y += row.height + 10
            composite = save_dir / "figure_4.png"
            page.save(composite)
            outputs.append(str(composite))
    except Exception:
        pass

    return outputs


def run_data1_figure5_workflow(data_root=None, save_dir=None):
    """Recreate DATA1 Figure 5 using the notebook's diafiltration sigma contours."""
    root = _resolve_data1_figure2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data1_paper_figures"
    save_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        ("A", root / "511.12 concpolar" / "contourdata-x_sigma-y_Lp.csv", "figure5_panel_a"),
        ("B", root / "511.11 concpolar" / "contourdata-x_sigma-y_Lp.csv", "figure5_panel_b"),
        ("C", root / "511.12" / "contourdata-x_sigma-y_Lp.csv", "figure5_panel_c"),
    ]

    outputs = []
    panel_groups = []
    for _, csv_path, prefix in cases:
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        panel_paths = _plot_contour_data1_legacy(df, prefix, save_dir, show_title=False, preface=False)
        outputs.extend(panel_paths)
        panel_groups.append([Path(path) for path in panel_paths])

    try:
        from PIL import Image, ImageDraw, ImageFont

        if len(panel_groups) == 3 and all(path.exists() for group in panel_groups for path in group):
            rows = []
            for group in panel_groups:
                imgs = [Image.open(path).convert("RGB") for path in group]
                target_h = max(im.height for im in imgs)
                resized = []
                for im in imgs:
                    scale = target_h / im.height
                    resized.append(im.resize((int(im.width * scale), target_h), Image.Resampling.LANCZOS))
                row_w = sum(im.width for im in resized) + 40
                row_h = target_h + 40
                canvas = Image.new("RGB", (row_w, row_h), "white")
                x = 20
                for im in resized:
                    canvas.paste(im, (x, 20))
                    x += im.width + 10
                rows.append(canvas)

            width = max(im.width for im in rows) + 20
            height = sum(im.height for im in rows) + 30
            page = Image.new("RGB", (width, height), "white")
            draw = ImageDraw.Draw(page)
            try:
                font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 28)
            except Exception:
                font = None
            y = 10
            for label, row in zip(["A", "B", "C"], rows):
                draw.text((10, y + 5), label, fill="black", font=font)
                page.paste(row, (30, y))
                y += row.height + 10
            composite = save_dir / "figure_5.png"
            page.save(composite)
            outputs.append(str(composite))
    except Exception:
        pass

    return outputs


def run_data1_figure6_workflow(data_root=None, save_dir=None):
    """Recreate DATA1 Figure 6 using the notebook's diafiltration B contours."""
    root = _resolve_data1_figure2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data1_paper_figures"
    save_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        ("A", root / "511.12 concpolar" / "contourdata-x_B-y_Lp.csv", "figure6_panel_a"),
        ("B", root / "511.11 concpolar" / "contourdata-x_B-y_Lp.csv", "figure6_panel_b"),
        ("C", root / "511.12" / "contourdata-x_B-y_Lp.csv", "figure6_panel_c"),
    ]

    outputs = []
    panel_groups = []
    for _, csv_path, prefix in cases:
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        panel_paths = _plot_contour_data1_legacy(df, prefix, save_dir, show_title=False, preface=False)
        outputs.extend(panel_paths)
        panel_groups.append([Path(path) for path in panel_paths])

    try:
        from PIL import Image, ImageDraw, ImageFont

        if len(panel_groups) == 3 and all(path.exists() for group in panel_groups for path in group):
            rows = []
            for group in panel_groups:
                imgs = [Image.open(path).convert("RGB") for path in group]
                target_h = max(im.height for im in imgs)
                resized = []
                for im in imgs:
                    scale = target_h / im.height
                    resized.append(im.resize((int(im.width * scale), target_h), Image.Resampling.LANCZOS))
                row_w = sum(im.width for im in resized) + 40
                row_h = target_h + 40
                canvas = Image.new("RGB", (row_w, row_h), "white")
                x = 20
                for im in resized:
                    canvas.paste(im, (x, 20))
                    x += im.width + 10
                rows.append(canvas)

            width = max(im.width for im in rows) + 20
            height = sum(im.height for im in rows) + 30
            page = Image.new("RGB", (width, height), "white")
            draw = ImageDraw.Draw(page)
            try:
                font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 28)
            except Exception:
                font = None
            y = 10
            for label, row in zip(["A", "B", "C"], rows):
                draw.text((10, y + 5), label, fill="black", font=font)
                page.paste(row, (30, y))
                y += row.height + 10
            composite = save_dir / "figure_6.png"
            page.save(composite)
            outputs.append(str(composite))
    except Exception:
        pass

    return outputs


def _data1_contour_theta_from_fit(fit_stru):
    """Extract a flat theta dictionary from a legacy DATA1 fit bundle."""
    if not isinstance(fit_stru, dict):
        return {}

    theta = {}

    def _as_float(value):
        if isinstance(value, dict):
            if not value:
                return None
            try:
                value = next(iter(value.values()))
            except Exception:
                return None
        arr = np.asarray(value, dtype=float)
        if arr.size == 0:
            return None
        return float(arr.reshape(-1)[0])

    params = fit_stru.get("parameters", {})
    if isinstance(params, dict):
        search_spaces = [params, fit_stru]
    else:
        search_spaces = [fit_stru]

    for key in ("Lp", "B", "sigma", "beta_0", "beta_1", "beta_2", "beta_3", "S0", "S"):
        for space in search_spaces:
            if key not in space:
                continue
            val = _as_float(space[key])
            if val is not None:
                theta[key] = val
                break

    if "B" not in theta:
        if "beta_0" in theta:
            theta["B"] = float(theta["beta_0"])
        elif "beta_1" in theta:
            theta["B"] = float(theta["beta_1"])

    theta.setdefault("sigma", 1.0)
    theta.setdefault("Lp", 1.0)
    return theta


def _data1_direct_contour_grid(exp, theta_base, x_name, x_values, lp_values):
    """Evaluate one DATA1 contour grid using the unified direct simulator."""
    import sys

    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library import (
        evaluate_data1_paper_contour_objectives_v24,
        simulate_data1_vialwise_trajectories_v24,
    )

    rows = []
    for x_val in x_values:
        for lp_val in lp_values:
            theta = dict(theta_base)
            theta[x_name] = float(x_val)
            theta["Lp"] = float(lp_val)
            sim = simulate_data1_vialwise_trajectories_v24(exp, theta)
            scores = evaluate_data1_paper_contour_objectives_v24(exp, sim).as_dict()
            rows.append(
                {
                    x_name: float(x_val),
                    "Lp": float(lp_val),
                    "Obj_mass": float(scores["mass_log10"]),
                    "Obj_concentration": float(scores["permeate_log10"]),
                    "Obj_retentate_concentration": float(scores["retentate_log10"]),
                }
            )
    return pd.DataFrame(rows)


def run_data1_direct_contour_branch(data_root=None, save_dir=None, grid_density=50):
    """Regenerate the DATA1 contour-grid branch directly from the unified evaluator.

    This is the DATA1-only analogue of the legacy MATLAB ``calc_contour_2d.m``
    path: each contour panel is computed on a parameter grid, written to CSV,
    and then rendered with the notebook-style contour panel helper.
    """
    root = _resolve_data1_figure2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data1_direct_contours"
    if save_dir.suffix:
        save_dir = save_dir.parent
    save_dir.mkdir(parents=True, exist_ok=True)

    import sys

    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library import (
        load_experiment_easy,
    )

    page_specs = [
        ("figure_5_direct", [
            ("A", "511.12 concpolar"),
            ("B", "511.11 concpolar"),
            ("C", "511.12"),
        ], "sigma"),
        ("figure_6_direct", [
            ("A", "511.12 concpolar"),
            ("B", "511.11 concpolar"),
            ("C", "511.12"),
        ], "B"),
        ("figure_s3_direct", [
            ("A", "511.12 concpolar"),
            ("B", "511.11 concpolar"),
            ("C", "511.12"),
            ("D", "511.11"),
        ], "sigma"),
        ("figure_s4_direct", [
            ("A", "511.12 concpolar"),
            ("B", "511.11 concpolar"),
            ("C", "511.12"),
            ("D", "511.11"),
        ], "B"),
        ("figure_s5_direct", [
            ("A", "501.1 concpolar"),
            ("B", "501.11 concpolar"),
            ("C", "501.1"),
            ("D", "501.11"),
        ], "sigma"),
        ("figure_s6_direct", [
            ("A", "501.1 concpolar"),
            ("B", "501.11 concpolar"),
            ("C", "501.1"),
            ("D", "501.11"),
        ], "B"),
    ]

    outputs = []
    x_ranges = {
        "sigma": (0.0, 1.0),
        "B": (0.0, 2.0),
    }
    lp_range = (0.5, 7.4)

    for page_name, cases, x_name in page_specs:
        row_groups = []
        row_labels = []
        for label, folder in cases:
            run_id = folder.split()[0]
            data_path = root / f"data_stru-dataset{run_id}.mat"
            fit_path = root / folder / "fit_stru.mat"
            if not data_path.exists() or not fit_path.exists():
                continue

            exp, _ = load_experiment_easy(str(data_path), selector="data_stru")
            fit_bundle = loadmat(str(fit_path)).get("fit_stru")
            theta_base = _data1_contour_theta_from_fit(fit_bundle)
            if not theta_base:
                continue

            x_lb, x_ub = x_ranges[x_name]
            x_values = np.linspace(x_lb, x_ub, int(grid_density))
            lp_values = np.linspace(lp_range[0], lp_range[1], int(grid_density))
            df = _data1_direct_contour_grid(exp, theta_base, x_name, x_values, lp_values)

            prefix = f"{page_name}_{label.lower()}"
            csv_path = save_dir / f"{prefix}.csv"
            # Defensive: ensure the target directory exists at write time, in
            # line with the out_path.parent.mkdir pattern used elsewhere in
            # this file (e.g. lines 6584, 6732, 6768, 8226). Guards against
            # the rare cases where the up-front save_dir.mkdir at function
            # entry does not survive to this point (caller-supplied path
            # quirks, filesystem races, etc.).
            csv_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(csv_path, index=False)
            outputs.append(str(csv_path))

            panel_paths = _plot_contour_data1_legacy(
                df, prefix, save_dir, show_title=False, preface=False
            )
            outputs.extend(panel_paths)
            row_groups.append([Path(path) for path in panel_paths])
            row_labels.append(label)

        composite_name = page_name
        composite = _compose_data1_panel_grid(
            row_groups,
            save_dir / f"{composite_name}.png",
            row_labels=row_labels,
            label_mode="per_row",
            panel_margin=12,
            row_margin=18,
            outer_margin=22,
        )
        if composite is not None:
            outputs.append(composite)

    return outputs


# =============================================================================
# NF270 (DATA3) objective-contour panels
# -----------------------------------------------------------------------------
# Mirrors the DATA1 ``calc_contour_2d.m`` / ``plot_contour.m`` workflow:
#
#   1. Fit the single-salt sheet once to obtain a centered theta_fit.
#   2. For each (x in {B, sigma}) sweep a 2D grid over (x, Lp), holding the
#      third parameter at its fitted value.
#   3. For every grid point, forward-simulate via solve_model(sim_opt=True)
#      and read the three per-channel objective components save_model already
#      stores on fit_stru: obj_m, obj_cv, obj_cr.
#   4. Save the grid as ``contourdata-x_{B,sigma}-y_Lp.csv`` (same five-column
#      schema as the DATA1 paper) and render the 1x3 panel figure via the
#      existing ``plot_contour`` renderer at ``refactored_ucb_library.py:4346``.
#
# Ranges mirror calc_contour_2d.m:
#   Lp     : linspace(0.1 * Lp_fit, 2.0 * Lp_fit)  (asymmetric, physical floor)
#   sigma  : linspace(0.0, 1.0)                    (full reflection-coefficient
#                                                    interval; do NOT log-space)
#   B      : linspace(max(1e-4, 0.1*|B_fit|), 3.0*|B_fit|)
#
# Default grid_density=20 keeps the first pass interactive (~10-30 min/sheet
# on the NF270 mesh); bump to 50 to match the DATA1 paper resolution.
# =============================================================================

def _nf270_contour_objectives_at_theta(
    data_stru,
    theta,
    *,
    mode="Lag",
    B_form="single",
    workflow_family="DATA3",
    nfe=80,
    band_vials=None,
):
    """Forward-simulate one NF270 sheet at a fixed theta; return per-channel objectives.

    Returns
    -------
    tuple of three floats
        (obj_mass, obj_permeate_conc, obj_retentate_conc) at the supplied theta.
        Any failed solve returns (nan, nan, nan).
    band_vials, when given, restricts the per-channel objectives to that
    concentration band's vials (band-masked contour).
    """
    try:
        fit_stru, _sim_stru, _sim_inter = solve_model(
            data_stru,
            mode,
            theta=theta,
            sim_opt=True,
            B_form=B_form,
            workflow_family=workflow_family,
            nfe=nfe,
            LOUD=False,
            band_vials=band_vials,
        )
    except Exception as exc:
        print(f"  [contour] forward sim raised: {exc!r}")
        return float("nan"), float("nan"), float("nan")

    if not isinstance(fit_stru, dict):
        return float("nan"), float("nan"), float("nan")

    def _safe_float(key):
        try:
            return float(fit_stru.get(key, float("nan")))
        except (TypeError, ValueError):
            return float("nan")

    return _safe_float("obj_m"), _safe_float("obj_cv"), _safe_float("obj_cr")


def _nf270_contour_grid_dataframe(
    data_stru,
    theta_fit,
    x_var,
    y_var,
    *,
    grid_density=20,
    mode="Lag",
    B_form="single",
    workflow_family="DATA3",
    nfe=80,
    progress_every=10,
    band_vials=None,
):
    """Build one 2D contour DataFrame at fixed theta with two parameters swept.

    Column schema matches the DATA1 paper CSVs exactly:
        [x_var, y_var, Obj_mass, Obj_concentration, Obj_retentate_concentration]
    Objective columns are log10-transformed (matches ``calc_contour_2d.m``);
    non-positive or NaN values become NaN so the renderer can mask them.
    """
    if x_var == y_var:
        raise ValueError("x_var and y_var must differ.")
    for var in (x_var, y_var):
        if var not in ("Lp", "B", "sigma", "beta_0", "beta_1", "beta_2"):
            raise ValueError(
                f"Unsupported contour axis {var!r}; expected Lp/B/sigma or "
                "beta_0/beta_1/beta_2 (the latter require a numeric-B_form centering fit)."
            )

    Lp_fit = float(theta_fit["Lp"])
    sigma_fit = float(theta_fit.get("sigma", 0.5))
    B_raw = theta_fit.get("B", theta_fit.get("beta_0", 1.0))
    if isinstance(B_raw, dict):
        B_fit = float(next(iter(B_raw.values())))
    else:
        B_fit = float(B_raw)
    B_abs = max(abs(B_fit), 1e-4)

    def _axis_values(var):
        if var == "Lp":
            return np.linspace(0.1 * Lp_fit, 2.0 * Lp_fit, grid_density)
        if var == "sigma":
            return np.linspace(0.0, 1.0, grid_density)
        if var == "B":
            return np.linspace(0.1 * B_abs, 3.0 * B_abs, grid_density)
        # Concentration-dependent solute permeability coefficients (B_form 1/2):
        # B = Jw*[beta_0 + beta_1*c^(B_form-1) + beta_2*c^B_form].  Sweeping a beta
        # instead of the lumped B exposes how the concentration-dependence (NF270)
        # couples to Lp / sigma.  Requires a centering fit done with numeric B_form.
        if var == "beta_0":
            b0 = float(theta_fit.get("beta_0", B_abs if np.isfinite(B_abs) else 1.0))
            return np.linspace(0.1 * b0, 3.0 * b0, grid_density)
        if var in ("beta_1", "beta_2"):
            bk = float(theta_fit.get(var, 0.0))
            span = max(abs(bk), 1.0)              # beta_1/beta_2 can be negative
            return np.linspace(bk - 2.0 * span, bk + 2.0 * span, grid_density)
        raise ValueError(var)

    x_vals = _axis_values(x_var)
    y_vals = _axis_values(y_var)

    rows = []
    total = len(x_vals) * len(y_vals)
    t0 = time.time()
    for i, x_val in enumerate(x_vals):
        for j, y_val in enumerate(y_vals):
            theta = dict(theta_fit)
            theta[x_var] = float(x_val)
            theta[y_var] = float(y_val)

            obj_m, obj_cv, obj_cr = _nf270_contour_objectives_at_theta(
                data_stru,
                theta,
                mode=mode,
                B_form=B_form,
                workflow_family=workflow_family,
                nfe=nfe,
                band_vials=band_vials,
            )

            def _log10_or_nan(val):
                if not np.isfinite(val) or val <= 0:
                    return float("nan")
                return float(np.log10(val))

            rows.append({
                x_var: float(x_val),
                y_var: float(y_val),
                "Obj_mass": _log10_or_nan(obj_m),
                "Obj_concentration": _log10_or_nan(obj_cv),
                "Obj_retentate_concentration": _log10_or_nan(obj_cr),
            })

            done = i * len(y_vals) + j + 1
            if progress_every and done % progress_every == 0:
                dt = time.time() - t0
                rate = done / max(dt, 1e-9)
                eta = (total - done) / max(rate, 1e-9)
                print(
                    f"  [contour {x_var}-{y_var}] {done}/{total} "
                    f"({100.0*done/total:.1f}%)  rate={rate:.2f}/s  eta={eta/60:.1f} min"
                )

    return pd.DataFrame(rows)


def run_nf270_contour_for_sheet(
    run_id,
    *,
    save_dir,
    grid_density=20,
    theta_fit=None,
    data_root=None,
    workbook_path=None,
    sheet_name=None,
    nfe=80,
    mode="Lag",
    B_form="single",
    band_vials=None,
):
    """Generate DATA1-style objective-contour panels for one NF270 single-salt sheet.

    Produces four artifacts in ``save_dir`` (two per sweep), matching the
    file naming used by ``DiafiltrationPaperPlots.ipynb``:

        contourdata-x_B-y_Lp.csv         contourdata-x_sigma-y_Lp.csv
        objcontour-x_B-y_Lp.png          objcontour-x_sigma-y_Lp.png

    If ``theta_fit`` is None, a single-shot ``solve_model`` fit centers the
    grid. Pass a pre-computed ``theta_fit`` (the parameters dict from a prior
    fit) to skip that step.

    Arguments
    ---------
    run_id : str
        NF270_RUN_REGISTRY key (e.g. "MC3.07.22.24_SNaCl"). Optional if
        ``workbook_path`` and ``sheet_name`` are provided explicitly.
    save_dir : str or Path
        Output directory; created if missing.
    grid_density : int
        Number of grid points per axis. 20 for fast iteration; 50 to match
        the DATA1 paper resolution.
    theta_fit : dict, optional
        Pre-computed {"Lp", "B", "sigma", ...}. If None, fit once.
    nfe : int
        Finite-difference mesh for the forward simulator (default 80).

    Returns
    -------
    dict
        {"run_id", "theta_fit", "csv_paths", "png_paths", "grid_density"}.
    """
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # Three contour slices per sheet. Each slice holds the third parameter at
    # the centering-fit value and sweeps the other two:
    #   * (B,  Lp) — σ pinned   → permeability ↔ mass-transfer coupling
    #   * (σ,  Lp) — B pinned   → permeability ↔ rejection coupling
    #   * (σ,   B) — Lp pinned  → rejection  ↔ mass-transfer coupling   (new)
    #
    # When B_form is a numeric concentration-dependent form (B_form >= 1), the
    # centering fit yields beta coefficients (beta_0, beta_1, ...) in place of a
    # lumped B, so we sweep the beta coefficients instead.  Gated on
    # ``isinstance(B_form, int) and B_form >= 1`` (bool excluded) so that
    # 'single'/'pervial'/'convection' and fractional exponents keep the
    # lumped-B slices.  Y-priority rule L_p > B(beta) > sigma:
    #   * (beta_0, Lp)     — L_p on Y  (beta_0 plays the role of B; L_p outranks it)
    #   * (sigma,  Lp)     — L_p on Y
    #   * (sigma,  beta_1) — beta_1 on Y  (beta outranks sigma)
    if isinstance(B_form, int) and not isinstance(B_form, bool) and B_form >= 1:
        sweep_pairs = (
            ("beta_0", "Lp"),
            ("sigma",  "Lp"),
            ("sigma",  "beta_1"),
        )
    else:
        sweep_pairs = (
            ("B",     "Lp"),
            ("sigma", "Lp"),
            ("sigma", "B"),
        )

    # Checkpoint: if all expected CSVs already exist and are non-empty, skip.
    # Per-slice resume happens inside the loop below.
    expected_csvs = tuple(
        save_dir / f"contourdata-x_{x}-y_{y}.csv" for (x, y) in sweep_pairs
    )
    if all(p.exists() and p.stat().st_size > 0 for p in expected_csvs):
        print(f"\n[contour] === {run_id} ===  [SKIP — all 3 slices already complete]")
        return {
            "run_id": run_id,
            "theta_fit": None,
            "csv_paths": list(expected_csvs),
            "png_paths": [
                save_dir / f"objcontour-x_{x}-y_{y}.png" for (x, y) in sweep_pairs
            ],
            "grid_density": grid_density,
            "skipped": True,
        }

    family = NF270_RUN_REGISTRY.get(str(run_id)) if run_id else None
    if family is None and (workbook_path is None or sheet_name is None):
        raise KeyError(
            f"Unknown NF270 run {run_id!r} and no workbook_path/sheet_name given."
        )
    if family is not None:
        root = Path(data_root) if data_root else NF270_DEFAULT_ROOT
        if workbook_path is None:
            workbook_path = root / family["workbook"]
        if sheet_name is None:
            sheet_name = family["sheet"]

    print(f"\n[contour] === {run_id} ===")
    print(f"[contour] workbook: {workbook_path}")
    print(f"[contour] sheet:    {sheet_name}")
    print(f"[contour] save_dir: {save_dir}")

    loaded = loadxlsx(workbook_path, sheet=sheet_name)
    data_stru = loaded["data_stru"]

    if theta_fit is None:
        print(f"\n[contour] centering fit (single-shot, no multistart)...")
        fit_stru, _sim_stru, _sim_inter = solve_model(
            data_stru,
            mode,
            sim_opt=False,
            B_form=B_form,
            workflow_family="DATA3",
            nfe=nfe,
            LOUD=True,
            band_vials=band_vials,
        )
        if not isinstance(fit_stru, dict):
            raise RuntimeError(
                f"Initial fit failed for {run_id}; cannot center the contour grid."
            )
        theta_fit = dict(fit_stru.get("parameters", {}))
        # Surface the centering fit so the caller can also save it.
        fit_json = save_dir / "centering_fit.json"
        try:
            with open(fit_json, "w") as fh:
                json.dump(
                    {
                        "run_id": run_id,
                        "theta_fit": theta_fit,
                        "obj_m": fit_stru.get("obj_m"),
                        "obj_cv": fit_stru.get("obj_cv"),
                        "obj_cr": fit_stru.get("obj_cr"),
                    },
                    fh,
                    indent=2,
                    default=float,
                )
            print(f"[contour] saved centering fit: {fit_json}")
        except Exception as exc:
            print(f"[contour] warning: failed to write {fit_json}: {exc}")

    print(f"\n[contour] theta_fit center = {theta_fit}")
    print(f"[contour] grid = {grid_density} x {grid_density} per sweep, "
          f"{len(sweep_pairs)} sweeps")

    csv_paths = []
    png_paths = []
    for x_var, y_var in sweep_pairs:
        csv_path = save_dir / f"contourdata-x_{x_var}-y_{y_var}.csv"
        png_path = save_dir / f"objcontour-x_{x_var}-y_{y_var}.png"

        # Per-slice checkpoint: if the CSV exists and is non-empty, reuse it.
        # (Previous 2-slice batches drop in here so we only compute σ-B.)
        if csv_path.exists() and csv_path.stat().st_size > 0:
            print(f"\n[contour] -------- {x_var}-{y_var} CSV already on disk, reusing --------")
            df = pd.read_csv(csv_path)
        else:
            print(f"\n[contour] -------- sweep x={x_var}, y={y_var} --------")
            df = _nf270_contour_grid_dataframe(
                data_stru,
                theta_fit,
                x_var=x_var,
                y_var=y_var,
                grid_density=grid_density,
                mode=mode,
                B_form=B_form,
                workflow_family="DATA3",
                nfe=nfe,
                band_vials=band_vials,
            )
            df.to_csv(csv_path, index=False)
        csv_paths.append(csv_path)

        # Combined viridis heatmap (quick at-a-glance review)
        try:
            fig, _axes = plot_contour(df, save_path=png_path)
            try:
                plt.close(fig)
            except Exception:
                pass
            png_paths.append(png_path)
            print(f"[contour] wrote {csv_path.name}  +  {png_path.name}")
        except Exception as exc:
            print(f"[contour] warning: render of {png_path.name} failed: {exc}")

        # DATA1-style per-channel contour-line plots (paper quality)
        try:
            paper_outputs = _plot_contour_paper_style(
                df,
                x_var=x_var,
                y_var=y_var,
                theta_fit=theta_fit,
                output_prefix=f"paper-x_{x_var}-y_{y_var}",
                save_dir=save_dir,
            )
            for p in paper_outputs:
                png_paths.append(Path(p))
        except Exception as exc:
            print(f"[contour] warning: paper-style render for {x_var}-{y_var} failed: {exc}")

    return {
        "run_id": run_id,
        "theta_fit": theta_fit,
        "csv_paths": csv_paths,
        "png_paths": png_paths,
        "grid_density": grid_density,
    }


def run_nf270_matlab_ports(
    run_ids,
    *,
    save_dir,
    grid_density=20,
    sigma_levels=6,
    nfe=80,
    data_root=None,
    mode=None,
    B_form="single",
    do_2d=True,
    do_3d=True,
    do_doe=False,
    do_sigma_heatmap=False,
    condition_var=None,
    condition_values=None,
    delp_values_psi=None,
):
    """Runfile-facing driver for the MATLAB-port analyses on NF270 sheets.

    Single entry point the runfile can call to exercise the Python ports of the
    legacy MATLAB toolchain — all evaluated with the DATA3 model + error spec:

        calc_contour_2d_py            -> contourdata-x_<x>-y_<y>.csv  (B-Lp, sigma-Lp, sigma-B)
        calc_contour_3d_py            -> contour3ddata.csv
        doe_heatmap_py                -> doe_heatmap-<cond>.csv         (A/D/E/modE optimality)
        heatmap_sigma_sensitivity_py  -> sigma_sensitivity_heatmap-<cond>.csv

    For each sheet a single-shot centering fit provides theta (cached as
    centering_fit.json and reused on resume). Heavy MBDoE sweeps (DoE / sigma
    heatmaps) are OFF by default; enable with do_doe / do_sigma_heatmap.

    Arguments mirror run_nf270_contour_branch; condition_* override the default
    MBDoE grid (C_F0 for filtration / C_D for diafiltration; 5-80 mM x 20-120 psi).

    Returns one dict per processed sheet: {run_id, theta, files}.
    """
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    if run_ids is None or run_ids == "all":
        resolved = tuple(NF270_RUN_REGISTRY.keys())
    else:
        resolved = tuple(str(r) for r in run_ids)
    resolved = tuple(r for r in resolved if r in NF270_RUN_REGISTRY)
    if not resolved:
        print("\n[matlab-ports] no NF270 sheets in the selection; nothing to do.")
        return []

    root = Path(data_root) if data_root else NF270_DEFAULT_ROOT
    print(f"\n[matlab-ports] {len(resolved)} sheet(s); "
          f"2d={do_2d} 3d={do_3d} doe={do_doe} sigma_heatmap={do_sigma_heatmap}")
    # Capture any explicit mode override; otherwise each sheet drives the model
    # in its OWN detected mode (data_stru["mode"] in {Lag, Overflow, DATA}).
    mode_override = mode
    outputs = []
    for run_id in resolved:
        fam = NF270_RUN_REGISTRY[run_id]
        sheet_dir = save_dir / run_id
        sheet_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n[matlab-ports] === {run_id} ===")
        try:
            data_stru = loadxlsx(root / fam["workbook"], sheet=fam["sheet"])["data_stru"]
        except Exception as exc:
            print(f"[matlab-ports] load failed for {run_id}: {exc!r}")
            continue
        # Per-sheet mode: an explicit `mode` arg wins; else use the sheet's own
        # detected mode so Lag/Overflow/DATA sheets are each modeled correctly.
        mode = str(mode_override or data_stru.get("mode") or "Lag")
        # Mode detection is heuristic (_detect_nf270_mode); log the source so a
        # mis-detected sheet is catchable rather than silently mis-modeled.
        print(f"[matlab-ports] mode = {mode} (source: {data_stru.get('mode_source', 'unknown')})")

        # Centering fit -> theta (reuse cached centering_fit.json on resume).
        theta = None
        cached = sheet_dir / "centering_fit.json"
        if cached.exists():
            try:
                theta = dict(json.loads(cached.read_text()).get("theta_fit", {}))
            except Exception:
                theta = None
        if not theta:
            try:
                fit_stru, _s, _i = solve_model(
                    data_stru, mode, sim_opt=False, B_form=B_form,
                    workflow_family="DATA3", nfe=nfe, LOUD=False,
                )
                theta = dict(fit_stru.get("parameters", {})) if isinstance(fit_stru, dict) else None
            except Exception as exc:
                print(f"[matlab-ports] centering fit failed for {run_id}: {exc!r}")
                theta = None
            if theta:
                try:
                    cached.write_text(json.dumps(
                        {"run_id": run_id, "theta_fit": theta}, indent=2, default=float))
                except Exception:
                    pass
        if not theta:
            print(f"[matlab-ports] no theta for {run_id}; skipping.")
            continue

        rec = {"run_id": run_id, "theta": theta, "files": []}

        if do_2d:
            for x_var, y_var in (("B", "Lp"), ("sigma", "Lp"), ("sigma", "B")):
                try:
                    calc_contour_2d_py(
                        data_stru, theta, x_var, y_var,
                        grid_density=grid_density, save_dir=sheet_dir,
                        mode=mode, B_form=B_form, workflow_family="DATA3", nfe=nfe,
                    )
                    rec["files"].append(str(sheet_dir / f"contourdata-x_{x_var}-y_{y_var}.csv"))
                except Exception as exc:
                    print(f"[matlab-ports] 2D {x_var}-{y_var} failed: {exc!r}")

        if do_3d:
            try:
                calc_contour_3d_py(
                    data_stru, theta, grid_density=grid_density,
                    sigma_levels=sigma_levels, save_dir=sheet_dir,
                    mode=mode, B_form=B_form, workflow_family="DATA3", nfe=nfe,
                )
                rec["files"].append(str(sheet_dir / "contour3ddata.csv"))
            except Exception as exc:
                print(f"[matlab-ports] 3D failed: {exc!r}")

        if do_doe or do_sigma_heatmap:
            # NF270 sheets carry mode in {Lag, Overflow, DATA}, never "D"; detect
            # diafiltration from a nonzero dialysate concentration instead.
            _dc = data_stru.get("data_config", {})
            _is_diaf = float(_dc.get("C_D", 0.0) or 0.0) > 0.0
            cv = condition_var or ("C_D" if _is_diaf else "C_F0")
            vals = list(condition_values) if condition_values else list(range(5, 81, 15))
            dps = list(delp_values_psi) if delp_values_psi else list(range(20, 121, 25))
            if do_doe:
                try:
                    doe_heatmap_py(
                        data_stru, theta, condition_var=cv, condition_values=vals,
                        delp_values_psi=dps, mode=mode, B_form=B_form,
                        workflow_family="DATA3", nfe=nfe, save_dir=sheet_dir,
                    )
                    rec["files"].append(str(sheet_dir / f"doe_heatmap-{cv}.csv"))
                except Exception as exc:
                    print(f"[matlab-ports] DoE heatmap failed: {exc!r}")
            if do_sigma_heatmap:
                try:
                    heatmap_sigma_sensitivity_py(
                        data_stru, theta, condition_var=cv, condition_values=vals,
                        delp_values_psi=dps, mode=mode, B_form=B_form,
                        workflow_family="DATA3", nfe=nfe, save_dir=sheet_dir,
                    )
                    rec["files"].append(str(sheet_dir / f"sigma_sensitivity_heatmap-{cv}.csv"))
                except Exception as exc:
                    print(f"[matlab-ports] sigma-sensitivity heatmap failed: {exc!r}")

        outputs.append(rec)
        print(f"[matlab-ports] {run_id}: {len(rec['files'])} file(s)")

    print(f"\n[matlab-ports] done: {len(outputs)} sheet(s)")
    return outputs


def run_nf270_contour_branch(
    run_ids,
    *,
    save_dir,
    grid_density=20,
    nfe=80,
    data_root=None,
    mode="Lag",
    B_form="single",
):
    """Generate DATA1-style objective-contour panels for a set of NF270 sheets.

    Thin dispatcher that the runfile (or any caller) can invoke without
    knowing about per-sheet save-dir layout or registry filtering. Each
    sheet's artifacts live in ``save_dir / run_id /``.

    Arguments
    ---------
    run_ids : iterable of str or None or "all"
        Sheets to process. ``None`` or ``"all"`` expands to every entry in
        ``NF270_RUN_REGISTRY``. Any entries not in the registry (e.g. the
        ``tables.parameters_table`` aggregator) are silently skipped.
    save_dir : str or Path
        Parent directory. Per-sheet subdirectories are created here.
    grid_density : int
        Grid points per axis. 20 fast, 30 medium, 50 paper.
    nfe : int
        Finite-difference mesh for each forward simulation.
    data_root : Path, optional
        NF270 workbook root; defaults to ``NF270_DEFAULT_ROOT``.

    Returns
    -------
    list of dict
        One artifact dict per successfully processed sheet (the same shape
        ``run_nf270_contour_for_sheet`` returns).
    """
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    if run_ids is None or run_ids == "all":
        resolved = tuple(NF270_RUN_REGISTRY.keys())
    else:
        resolved = tuple(str(r) for r in run_ids)
    resolved = tuple(r for r in resolved if r in NF270_RUN_REGISTRY)

    if not resolved:
        print("\n[contour] no NF270 sheets in the selection; nothing to do.")
        return []

    n_pts = grid_density * grid_density * 2
    print(
        f"\n[contour] branch: {len(resolved)} sheet(s) x {n_pts} forward sims each"
    )
    print(f"[contour] save_dir = {save_dir}")

    outputs = []
    for run_id in resolved:
        sheet_save_dir = save_dir / run_id
        try:
            artifact = run_nf270_contour_for_sheet(
                run_id,
                save_dir=sheet_save_dir,
                grid_density=grid_density,
                nfe=nfe,
                data_root=data_root,
                mode=mode,
                B_form=B_form,
            )
            outputs.append(artifact)
        except Exception as exc:
            print(f"[contour] {run_id} FAILED: {exc!r}")

    if outputs:
        print(f"\n[contour] generated panels for {len(outputs)} sheet(s):")
        for artifact in outputs:
            for path in artifact.get("png_paths", []):
                print(f"  + {path}")
    else:
        print("\n[contour] no contour panels were generated.")
    return outputs


def run_data1_sigma_contours(data_root=None, save_dir=None):
    """Recreate the DATA1 sigma-sensitivity contour figures from the notebook."""
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "sigma_contours"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    cases = [
        (
            "501.1",
            root / "data_stru-dataset501.1.mat",
            root / "501.1 concpolar" / "contour_sig_stru-dat501.1.mat",
            0,
            3,
            [[1, 2, 3], [1, 2, 3], [1, 2, 3]],
            [0, 1, 2, 3],
            [[(20, 40), (35, 60), (50, 80)], [(20, 40), (35, 60), (50, 80)], [(20, 40), (35, 60), (50, 80)]],
            "",
            True,
        ),
        (
            "511.12-endvial1",
            root / "data_stru-dataset511.12.mat",
            root / "511.12 concpolar" / "contour_sig_stru-dat511.12-endvial1.mat",
            0,
            3,
            [[1, 1.5, 1.9], [1, 1.5, 1.9], [1, 1.5, 1.9]],
            [0, 1, 2, 3],
            [[(40, 60), (60, 100), (73.5, 110)], [(40, 60), (60, 100), (73.5, 110)], [(40, 60), (60, 100), (73.5, 110)]],
            "-endvial1",
            True,
        ),
        (
            "511.12-endvial5",
            root / "data_stru-dataset511.12.mat",
            root / "511.12 concpolar" / "contour_sig_stru-dat511.12-endvial5.mat",
            0,
            3,
            [[2, 5, 8], [2, 5, 8], [2, 5, 8]],
            [0, 1, 2, 3],
            [[(20, 40), (35, 60), (50, 80)], [(20, 40), (35, 60), (50, 80)], [(20, 40), (35, 60), (50, 80)]],
            "-endvial5",
            True,
        ),
        (
            "511.12-endvial10",
            root / "data_stru-dataset511.12.mat",
            root / "511.12 concpolar" / "contour_sig_stru-dat511.12-endvial10.mat",
            0,
            3,
            [[2, 5, 8], [2, 5, 8], [2, 5, 8]],
            [0, 1, 2, 3],
            [[(5, 10), (15, 30), (20, 45)], [(5, 10), (15, 30), (20, 45)], [(5, 10), (15, 30), (20, 45)]],
            "-endvial10",
            True,
        ),
    ]

    for label, data_path, contour_path, vmin, vmax, levels, ticks, manual, suffix, bar in cases:
        if not data_path.exists() or not contour_path.exists():
            continue
        data_id = "501.1" if "501.1" in label and "511" not in label else "511.12"
        try:
            variant = "concpolar"
            data_run = get_experimental_run("DATA1", data_id, variant=variant, data_root=root)
        except Exception:
            data_run = None
        data_stru = _normalize_conductivity_measurements(loadmat(str(data_path)).get("data_stru"))
        contour_sig_stru = loadmat(str(contour_path)).get("contour_sig_stru")
        outputs.extend(
            plot_contour_sig_sen(
            data_stru,
            contour_sig_stru,
            vmin,
            vmax,
            levels,
            ticks,
            manual,
            suffix,
            filled=True,
            bar=bar,
            save_dir=save_dir,
        )
        )
        plt.close("all")

    return outputs


def run_data1_concentration_comparison(data_root=None, save_dir=None):
    """Recreate the DATA1 concentration comparison figure from the notebook."""
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "concentration_comparison"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    pairs = [
        (
            root / "data_stru-dataset501.1.mat",
            root / "501.1 concpolar" / "fit_stru.mat",
            root / "data_stru-dataset511.12.mat",
            root / "511.12 concpolar" / "fit_stru.mat",
            save_dir / "concentration_comparison.png",
        ),
    ]

    for data_f_path, fit_f_path, data_d_path, fit_d_path, out_path in pairs:
        if not (data_f_path.exists() and fit_f_path.exists() and data_d_path.exists() and fit_d_path.exists()):
            continue
        data_run_f = ExperimentalRun("501.1", root, "data_stru-dataset501.1.mat", variant="base", subdir="501.1")
        data_run_d = ExperimentalRun("511.12", root, "data_stru-dataset511.12.mat", variant="concpolar", subdir="511.12 concpolar")
        data_stru_f = data_run_f.load_data()
        fit_stru_f = data_run_f.load_fit()
        data_stru_d = data_run_d.load_data()
        fit_stru_d = data_run_d.load_fit()
        fig = plot_conc_comparison(data_stru_f, fit_stru_f, data_stru_d, fit_stru_d, save_path=out_path)
        outputs.append(str(out_path))
        plt.close(fig)

    return outputs


def _show_saved_figures(paths):
    """Open saved figure files so the DATA1 plots appear like notebook output."""
    figure_paths = []
    for item in paths:
        path = Path(item)
        if path.suffix.lower() in {".png", ".jpg", ".jpeg", ".svg", ".pdf", ".eps"} and path.exists():
            figure_paths.append(path)

    for idx, fig_path in enumerate(sorted(figure_paths), start=1):
        try:
            image = plt.imread(fig_path)
        except Exception:
            continue
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(image)
        ax.axis("off")
        ax.set_title(f"{idx}: {fig_path.name}", fontsize=10)
    if figure_paths:
        plt.show()


def run_data1_notebook_workflow(data_root=None, save_dir=None, show=True):
    """
    Recreate the DATA1 notebook workflow using the utility-style model helpers
    and the paper-plot functions from DiafiltrationPaperPlots.ipynb.
    """
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data1_paper_figures"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    outputs.extend(run_data_analysis(data_root=root, save_dir=save_dir))
    experiment_space = root / "experiment space"
    filtration_csv = experiment_space / "filtration.csv"
    diafiltration_csv = experiment_space / "diafiltration.csv"
    if filtration_csv.exists() and diafiltration_csv.exists():
        df_f = pd.read_csv(filtration_csv, header=2)
        df_d = pd.read_csv(diafiltration_csv, header=2)
        fig, _ = plot_conc_range(df_f, df_d, save_path=save_dir / "concentration_range.png")
        outputs.append(str(save_dir / "concentration_range.png"))
        plt.close(fig)
    outputs.extend(run_data1_figure4_workflow(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_figure5_workflow(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_figure6_workflow(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_si_s2(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_si_s3(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_si_s4(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_si_s5(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_si_s6(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_sigma_contours(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_concentration_comparison(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_figure2_workflow(data_root=root, save_dir=save_dir))
    outputs = sorted(set(str(p) for p in outputs))

    if show:
        _show_saved_figures(outputs)

    return outputs


def _resolve_data2_root(data_root=None):
    """Find the folder that stores the DATA2 paper inputs."""
    repo_root = Path(__file__).resolve().parents[1]
    if data_root is None:
        env_root = os.environ.get("DIAFILTRATION_DATA2_ROOT")
        if env_root:
            return Path(env_root).expanduser().resolve()
        return repo_root / "legacy" / "data1_matlab" / "data_library"
    return Path(data_root)


def _resolve_data2_validation_fig9_root():
    """Find the committed validation CSVs for DATA2 Fig. 9."""
    return Path(__file__).resolve().parents[1] / "UnifiedFramework" / "DATA3" / "docs" / "validation" / "simulation_validation_data_files" / "data2_main" / "fig9"


def _resolve_data2_table_baseline_root():
    """Find the committed DATA2 Tables 3-6 baseline CSVs."""
    return Path(__file__).resolve().parents[1] / "UnifiedFramework" / "DATA3" / "results" / "reproduction" / "20260306-mainpaper-table-baseline" / "tables"


def _resolve_data2_figure_root():
    """Find the shared DATA2 figure folder used by the legacy plotting scripts."""
    return Path(__file__).resolve().parents[1] / "UnifiedFramework" / "DATA3" / "figures"


def _resolve_data2_main_extract_root():
    """Find the extracted DATA2 main-paper image folder."""
    return (
        Path(__file__).resolve().parents[1]
        / "UnifiedFramework"
        / "DATA3"
        / "results"
        / "reproduction"
        / "20260306-023356-paper-pdf-extract"
        / "pdf_extract"
        / "data2_main"
        / "images"
    )


def _resolve_data2_publication_source(filename, save_dir=None):
    """Resolve one DATA2 image from either the current output folder or shared figure cache."""
    candidates = []
    if save_dir is not None:
        candidates.append(Path(save_dir) / filename)
    candidates.extend(
        [
            _resolve_data2_figure_root() / filename,
            Path(__file__).resolve().parents[1] / filename,
        ]
    )
    for path in candidates:
        if path.exists():
            return path
    return None


def _resolve_data2_main_extract_image(filename):
    """Resolve one extracted DATA2 main-paper figure image."""
    path = _resolve_data2_main_extract_root() / filename
    return path if path.exists() else None


def _copy_data2_publication_image(src, out_path):
    """Copy one DATA2 panel image into the paper-artifact folder under a paper-style name."""
    try:
        from PIL import Image
    except Exception:
        return None

    src = Path(src)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        Image.open(src).convert("RGB").save(out_path)
        return out_path
    except Exception:
        return None


def _annotate_data2_panel_labels(image_path, labels):
    """Add publication-style panel letters to a composed DATA2 figure."""
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return None

    image_path = Path(image_path)
    if not image_path.exists():
        return None
    img = Image.open(image_path).convert("RGB")
    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 34)
    except Exception:
        font = None
    width, height = img.size
    for label, x_frac, y_frac in labels:
        draw.text((int(width * x_frac), int(height * y_frac)), label, fill="black", font=font)
    img.save(image_path)
    return image_path


def _compose_data2_panel_sheet(rows, out_path):
    """Compose a simple DATA2 paper-style panel sheet from existing panel PNGs."""
    try:
        from PIL import Image
    except Exception:
        return None

    valid_rows = []
    for row in rows:
        existing = [Path(p) for p in row if Path(p).exists()]
        if existing:
            valid_rows.append(existing)
    if not valid_rows:
        return None

    rendered_rows = []
    for row in valid_rows:
        imgs = [Image.open(path).convert("RGB") for path in row]
        target_h = max(im.height for im in imgs)
        resized = []
        for im in imgs:
            scale = target_h / im.height
            resized.append(im.resize((int(im.width * scale), target_h), Image.Resampling.LANCZOS))
        row_w = sum(im.width for im in resized) + 20 * (len(resized) + 1)
        row_h = target_h + 40
        canvas = Image.new("RGB", (row_w, row_h), "white")
        x = 20
        for im in resized:
            canvas.paste(im, (x, 20))
            x += im.width + 20
        rendered_rows.append(canvas)

    width = max(im.width for im in rendered_rows) + 20
    height = sum(im.height for im in rendered_rows) + 20 * (len(rendered_rows) + 1)
    page = Image.new("RGB", (width, height), "white")
    y = 20
    for row in rendered_rows:
        page.paste(row, ((width - row.width) // 2, y))
        y += row.height + 20

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    page.save(out_path)
    return out_path


def _compose_data2_figure6(rows, out_path):
    """Compose DATA2 Figure 6 with the paper's tight 2x3 layout."""
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return None

    valid_rows = []
    for row in rows:
        existing = [Path(p) for p in row if Path(p).exists()]
        if existing:
            valid_rows.append(existing)
    if len(valid_rows) != 2 or any(len(row) != 3 for row in valid_rows):
        return None

    target_h = 585
    left_margin = 38
    right_margin = 38
    top_margin = 16
    bottom_margin = 16
    col_gap = 20
    row_gap = 24

    panel_rows = []
    for row in valid_rows:
        panels = []
        for path in row:
            img = Image.open(path).convert("RGB")
            scale = target_h / img.height
            resized = img.resize((int(round(img.width * scale)), target_h), Image.Resampling.LANCZOS)
            panels.append(resized)
        panel_rows.append(panels)

    row_widths = [sum(im.width for im in row) + col_gap * (len(row) - 1) for row in panel_rows]
    page_w = max(row_widths) + left_margin + right_margin
    page_h = top_margin + target_h + row_gap + target_h + bottom_margin
    page = Image.new("RGB", (page_w, page_h), "white")

    labels = ["A", "B", "C", "D", "E", "F"]
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 34)
    except Exception:
        font = None

    draw = ImageDraw.Draw(page)
    label_index = 0
    y = top_margin
    for row in panel_rows:
        x = (page_w - (sum(im.width for im in row) + col_gap * (len(row) - 1))) // 2
        for im in row:
            page.paste(im, (x, y))
            draw.text((x + 6, y + 6), labels[label_index], fill="black", font=font)
            label_index += 1
            x += im.width + col_gap
        y += target_h + row_gap

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    page.save(out_path)
    return out_path


def _compose_data2_figure7(rows, out_path):
    """Compose DATA2 Figure 7 as a vertical A/B/C stack with room for labels."""
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return None

    valid_rows = []
    for row in rows:
        existing = [Path(p) for p in row if Path(p).exists()]
        if existing:
            valid_rows.append(existing)
    if len(valid_rows) != 3 or any(len(row) != 1 for row in valid_rows):
        return None

    target_w = 980
    left_margin = 42
    right_margin = 42
    top_margin = 16
    bottom_margin = 44
    row_gap = 26

    panels = []
    for row in valid_rows:
        img = Image.open(row[0]).convert("RGB")
        scale = target_w / img.width
        resized = img.resize((target_w, int(round(img.height * scale))), Image.Resampling.LANCZOS)
        panels.append(resized)

    page_w = left_margin + target_w + right_margin
    page_h = top_margin + sum(im.height for im in panels) + row_gap * (len(panels) - 1) + bottom_margin
    page = Image.new("RGB", (page_w, page_h), "white")

    labels = ["A", "B", "C"]
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 34)
    except Exception:
        font = None

    draw = ImageDraw.Draw(page)
    x = left_margin
    y = top_margin
    for label, im in zip(labels, panels):
        page.paste(im, (x, y))
        draw.text((x + 6, y + 6), label, fill="black", font=font)
        y += im.height + row_gap

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    page.save(out_path)
    return out_path


def run_data2_publication_figures(data_root=None, save_dir=None):
    """Assemble the paper-style DATA2 main and SI figure sheets from generated panels."""
    _ = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else _resolve_data2_figure_root()
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []

    aliases = {
        "figure_s1.png": "calib_curve.png",
        "figure_8.png": "partition_sensitivity.png",
    }
    for out_name, src_name in aliases.items():
        src = _resolve_data2_publication_source(src_name, save_dir=save_dir)
        if src is None:
            continue
        copied = _copy_data2_publication_image(src, save_dir / out_name)
        if copied is not None:
            outputs.append(str(copied))

    paper_extract_aliases = {
        "figure_2.png": "img-003.png",
        "figure_3.png": "img-004.png",
    }
    for out_name, src_name in paper_extract_aliases.items():
        src = _resolve_data2_main_extract_image(src_name)
        if src is None:
            continue
        copied = _copy_data2_publication_image(src, save_dir / out_name)
        if copied is not None:
            outputs.append(str(copied))

    compositions = {
        "figure_6.png": [
            [
                _resolve_data2_publication_source("mass-dat270511.123.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270511.123.png", save_dir=save_dir),
                _resolve_data2_publication_source("stirc_mass-dat270511.123.png", save_dir=save_dir),
            ],
            [
                _resolve_data2_publication_source("mass-dat270511.423.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270511.423.png", save_dir=save_dir),
                _resolve_data2_publication_source("stirc_mass-dat270511.423.png", save_dir=save_dir),
            ],
        ],
        "figure_7.png": [
            [
                _resolve_data2_publication_source("Js_predict0.png", save_dir=save_dir),
            ],
            [
                _resolve_data2_publication_source("Jw_predict.png", save_dir=save_dir),
            ],
            [
                _resolve_data2_publication_source("Js_predict1.png", save_dir=save_dir),
            ],
        ],
        "figure_9.png": [
            [
                _resolve_data2_publication_source("startup_barplot.png", save_dir=save_dir),
            ],
            [
                _resolve_data2_publication_source("concentrating_residuals_boxplot.png", save_dir=save_dir),
            ]
        ],
        "figure_s2.png": [
            [
                _resolve_data2_publication_source("mass-dat270611.121.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270611.121.png", save_dir=save_dir),
            ],
            [
                _resolve_data2_publication_source("mass-dat270711.121.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270711.121.png", save_dir=save_dir),
            ],
        ],
        "figure_s3.png": [
            [
                _resolve_data2_publication_source("mass-dat270511.221.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270511.221.png", save_dir=save_dir),
            ],
            [
                _resolve_data2_publication_source("mass-dat270511.321.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270511.321.png", save_dir=save_dir),
            ],
        ],
        "figure_s4.png": [
            [
                _resolve_data2_publication_source("mass-dat270511.421.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270511.421.png", save_dir=save_dir),
            ],
            [
                _resolve_data2_publication_source("mass-dat270511.921.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270511.921.png", save_dir=save_dir),
            ],
        ],
        "figure_s5.png": [
            [
                _resolve_data2_publication_source("mass-dat270511.521.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270511.521.png", save_dir=save_dir),
            ],
            [
                _resolve_data2_publication_source("mass-dat270511.621.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270511.621.png", save_dir=save_dir),
            ],
        ],
        "figure_s6.png": [
            [
                _resolve_data2_publication_source("mass-dat270511.721.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270511.721.png", save_dir=save_dir),
            ],
            [
                _resolve_data2_publication_source("mass-dat270511.821.png", save_dir=save_dir),
                _resolve_data2_publication_source("concentration-dat270511.821.png", save_dir=save_dir),
            ],
        ],
        "figure_s7.png": [
            [
                _resolve_data2_publication_source("Bpervial.png", save_dir=save_dir),
                _resolve_data2_publication_source("Js_Jw_cin.png", save_dir=save_dir),
            ]
        ],
        "figure_s8.png": [
            [
                _resolve_data2_publication_source("Js_predict0.png", save_dir=save_dir),
                _resolve_data2_publication_source("Jw_predict.png", save_dir=save_dir),
                _resolve_data2_publication_source("Js_predict.png", save_dir=save_dir),
            ]
        ],
    }

    for out_name, rows in compositions.items():
        valid_rows = [[p for p in row if p is not None] for row in rows]
        if out_name == "figure_6.png":
            composite = _compose_data2_figure6(valid_rows, save_dir / out_name)
        elif out_name == "figure_7.png":
            composite = _compose_data2_figure7(valid_rows, save_dir / out_name)
        else:
            composite = _compose_data2_panel_sheet(valid_rows, save_dir / out_name)
        if composite is not None:
            if out_name == "figure_6.png":
                pass
            elif out_name == "figure_7.png":
                pass
            elif out_name == "figure_9.png":
                _annotate_data2_panel_labels(
                    composite,
                    [("A", 0.02, 0.02), ("B", 0.02, 0.53)],
                )
            elif out_name == "figure_s8.png":
                _annotate_data2_panel_labels(
                    composite,
                    [("A", 0.02, 0.02), ("B", 0.35, 0.02), ("C", 0.68, 0.02)],
                )
            outputs.append(str(composite))

    figure8 = save_dir / "figure_8.png"
    if figure8.exists():
        _annotate_data2_panel_labels(
            figure8,
            [("A", 0.02, 0.01), ("B", 0.02, 0.35), ("C", 0.02, 0.69)],
        )

    return outputs


def model_predictions(c_in, c_h, k0, k1, Pe):
    """Compute Js/Jw from the simple DATA2 regression model."""
    Kf = k1 * c_in + k0
    Kp = k1 * c_h + k0
    return (c_in * Kf * np.exp(Pe) - Kp * c_h) / (np.exp(Pe) - 1)


def plot_model_predictions(sim_data, k0, k1, Pe, save_path=None):
    """Draw the DATA2 Js/Jw surface and overlay the measured data."""
    c_in = np.asarray(sim_data["cIn"])
    c_h = np.asarray(sim_data["cH"])
    Js = np.asarray(sim_data["Js"])
    Jw = np.asarray(sim_data["Jw"])

    round_to = 5
    c_in_low = np.floor(np.min(c_in) / round_to) * round_to
    c_in_high = np.ceil(np.max(c_in) / round_to) * round_to
    c_h_low = np.floor(np.min(c_h) / round_to) * round_to
    c_h_high = np.ceil(np.max(c_h) / round_to) * round_to

    c_in_grid, c_h_grid = np.meshgrid(
        np.linspace(c_in_low, c_in_high, 100),
        np.linspace(c_h_low, c_h_high, 100),
    )

    Js_Jw_model_grid = model_predictions(c_in_grid, c_h_grid, k0, k1, Pe)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    ax.mouse_init()
    ax.plot_wireframe(c_in_grid, c_h_grid, Js_Jw_model_grid, color="blue", alpha=0.7, label="Regressed Model")
    ax.scatter(c_in, c_h, Js / Jw, color="red", label="Experimental Data", marker="o", s=50)
    ax.legend(fontsize=10)
    ax.set_xlabel("c$_{in}$ [mM]", fontsize=12, fontweight="bold")
    ax.set_ylabel("c$_{h}$ [mM]", fontsize=12, fontweight="bold")
    ax.set_zlabel("J$_s$/J$_w$", fontsize=12, fontweight="bold")
    ax.set_title(f"J$_s$/J$_w$ with Pe={Pe:.1f}", fontsize=14, fontweight="bold")
    ax.set_box_aspect([1, 1, 0.8])
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.8)
    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig, ax


def plot_pressure_change(
    data_pd,
    *,
    time_shift,
    marker_time,
    pressure_xlim,
    pressure_ylim,
    conc_ylim,
    fig_name,
    save_dir=None,
):
    """Recreate one DATA2 pressure-change plot from the notebook."""
    if isinstance(data_pd, (str, Path)):
        data_pd = pd.read_csv(data_pd, skiprows=[1])

    save_dir = Path(save_dir) if save_dir is not None else None
    fig = plt.figure(figsize=(4, 4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    ax1.plot((data_pd.Time - time_shift) / 60, data_pd.Pressure, "o", color="tab:orange", markersize=6)
    ax2.plot((data_pd.Time - time_shift) / 60, data_pd.Concentration, "s", color="m", markersize=6)
    ax1.plot([0, 0], [0, 100], "--", color="tab:blue", lw=2)
    marker_rel = (float(marker_time) - float(time_shift)) / 60.0
    ax1.plot([marker_rel, marker_rel], [0, 100], "--", color="tab:green", lw=2)

    ax1.set_xlabel("Time [min]", fontsize=16, fontweight="bold")
    ax1.set_ylabel("Applied Pressure [psi]", color="tab:orange", fontsize=16, fontweight="bold")
    ax1.tick_params(axis="x", which="major", direction="in", labelsize=12)
    ax1.tick_params(axis="x", which="minor", direction="in", labelsize=8)
    ax1.tick_params(axis="y", labelcolor="tab:orange", direction="in", labelsize=12)
    ax1.set_xlim(*pressure_xlim)
    ax1.set_ylim(*pressure_ylim)

    ax2.set_ylabel("Retentate [mM]", color="m", fontsize=16, fontweight="bold")
    ax2.tick_params(axis="y", labelcolor="m", direction="in", labelsize=12)
    ax2.set_ylim(*conc_ylim)

    out_path = Path(fig_name)
    if save_dir is not None:
        out_path = save_dir / out_path.name
    fig.savefig(out_path, dpi=600, bbox_inches="tight", transparent=True)
    return fig, out_path


def _plot_data2_convection_fit_vs_concentration(sim_data, predicted_js, save_path):
    """Notebook-style Figure 7C: solute flux versus interface concentration."""
    my_data = sim_data
    fig = plt.figure(figsize=(4, 4))
    plt.plot(my_data["cIn"], my_data["Js"], "k-", label="J$_s$ (Empirical)", lw=2)
    plt.plot(
        my_data["cIn"],
        predicted_js,
        "r--",
        dashes=(8, 4),
        label="J$_s$ (Convection-diffusion)",
        lw=2,
    )
    plt.xlabel("c$\\mathbf{_{in,f}}$ [mM]", fontsize=16, fontweight="bold")
    plt.ylabel("J$\\mathbf{_s\\ [\\mu mol \\cdot cm^{-2} \\cdot s^{-1}]}$", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def _plot_data2_convection_fit_vs_time(sim_data, predicted_js, save_path):
    """Notebook-style Figure S8C: solute flux versus time."""
    my_data = sim_data
    fig = plt.figure(figsize=(4, 4))
    plt.plot(my_data["time"] / 60, my_data["Js"], "k-", label="J$_s$ (Empirical)", lw=2)
    plt.plot(
        my_data["time"] / 60,
        predicted_js,
        "r--",
        dashes=(7, 5),
        label="J$_s$ (Convection-diffusion)",
        lw=2,
    )
    plt.xlabel("Time [min]", fontsize=16, fontweight="bold")
    plt.ylabel("J$\\mathbf{_s\\ [\\mu mol \\cdot cm^{-2} \\cdot s^{-1}]}$", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def _plot_data2_jw_vs_time(sim_data, save_path):
    """Notebook-style Figure 7B / S8B: water flux versus time."""
    my_data = sim_data
    fig = plt.figure(figsize=(4, 4))
    plt.plot(my_data["time"] / 60, my_data["Jw"] * 1e4, "b-", label="J$_{w}$ (Empirical)", lw=2)
    plt.xlabel("Time [min]", fontsize=16, fontweight="bold")
    plt.ylabel("J$\\mathbf{_w\\ [\\mu m \\cdot s^{-1}]}$", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def _plot_data2_cin_ch_vs_time(sim_data, save_path):
    """Notebook-style Figure 7A / S8A: interface concentrations versus time."""
    my_data = sim_data
    fig = plt.figure(figsize=(4, 4))
    plt.plot(my_data["time"] / 60, my_data["cIn"], "g-", label="c$_{in,f}$ (Empirical)", lw=2)
    plt.plot(my_data["time"] / 60, my_data["cH"], "r-", label="c$_{h}$ (Empirical)", lw=2)
    plt.xlabel("Time [min]", fontsize=16, fontweight="bold")
    plt.ylabel("Concentration [mM]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def plot_startup_barplot(improvements=None, *, save_path=None):
    """Recreate the small DATA2 startup improvement bar chart."""
    if improvements is None:
        improvements = [138, -9]
    modes = ["Lag", "Overflow"]
    colors = ["#2CA02C" if imp > 0 else "#D62728" for imp in improvements]
    df = pd.DataFrame({"Mode": modes, "Improvement": improvements, "Color": colors})

    fig = plt.figure(figsize=(6, 3))
    ax = fig.add_subplot(111)
    if sns is not None:
        ax = sns.barplot(
            data=df,
            x="Improvement",
            y="Mode",
            hue="Mode",
            order=["Overflow", "Lag"],
            palette={row["Mode"]: row["Color"] for _, row in df.iterrows()},
            width=0.6,
            legend=False,
        )
    else:
        y_pos = np.arange(len(df))
        ax.barh(y_pos, df["Improvement"], color=df["Color"], height=0.6)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df["Mode"])
    for i, (_, row) in enumerate(df.iterrows()):
        imp = float(row["Improvement"])
        text_color = "white" if imp > 0 else "black"
        label = f"{int(imp)}%" if float(imp).is_integer() else f"{imp}%"
        plt.text(
            imp - (10 if imp > 0 else -12),
            i,
            label,
            ha="right" if imp > 0 else "left",
            va="center",
            fontsize=14,
            fontweight="bold",
            color=text_color,
        )
    plt.xlim([-40, 150])
    ax.grid(False)
    ax.axvline(x=0, color="k")
    ax.invert_yaxis()
    ax.tick_params(axis="x", labelbottom="off")
    ax.tick_params(axis="y", direction="in", pad=-5)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(ax.get_yticklabels(), ha="left")
    plt.ylabel("")
    plt.xlabel("Information Improvement")
    plt.tight_layout()
    if save_path is None:
        save_path = "startup_barplot.png"
    fig.savefig(save_path, dpi=600, bbox_inches="tight")
    return fig, ax


def _sim_stru_to_dataframe(sim_stru):
    """Flatten the saved simulation structure into a tidy dataframe."""
    rows = []
    if isinstance(sim_stru, dict):
        vial_iter = sim_stru.items()
    else:
        vial_iter = enumerate(sim_stru, start=1)
    for vial_key, vial in vial_iter:
        if not isinstance(vial, dict):
            continue
        times = np.asarray(vial.get("time", []), dtype=float).reshape(-1)
        for idx, t in enumerate(times):
            row = {"vial": vial_key, "time": float(t)}
            for key in ("Js", "Jw", "cIn", "cH", "cF", "cV", "mF", "B"):
                if key in vial:
                    values = vial[key]
                    try:
                        arr = np.asarray(values, dtype=float).reshape(-1)
                        if arr.size == len(times):
                            row[key] = float(arr[idx])
                        elif arr.size > 0:
                            row[key] = float(arr[min(idx, arr.size - 1)])
                    except Exception:
                        pass
            rows.append(row)
    return pd.DataFrame(rows)


PE_LOWER_BOUND = 0.1
PE_INITIAL_VALUE = 15
PE_UPPER_BOUND = 20


def model_convection(sim_data, Pe_fixed_value=None, log_transform_Pe=False):
    """Fit the simple DATA2 convection model used in the visualization notebook."""
    data = sim_data.copy()
    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data)
    data = data.reset_index(drop=True)

    model = ConcreteModel()
    model.i = RangeSet(0, len(data) - 1)

    model.J_w = Param(model.i, initialize=lambda m, i: float(data.loc[i, "Jw"]), mutable=True)
    model.J_s = Param(model.i, initialize=lambda m, i: float(data.loc[i, "Js"]), mutable=True)
    model.c_in = Param(model.i, initialize=lambda m, i: float(data.loc[i, "cIn"]), mutable=True)
    model.c_h = Param(model.i, initialize=lambda m, i: float(data.loc[i, "cH"]), mutable=True)

    model.Js = Var(model.i, initialize=lambda m, i: float(data.loc[i, "Js"] / data.loc[i, "Jw"]))
    model.Kp = Var(model.i, initialize=1.0, bounds=(0.01, 1.2))
    model.Kf = Var(model.i, initialize=1.0, bounds=(0.01, 1.2))
    model.k1 = Var(initialize=1e-2)
    model.k0 = Var(initialize=1.0, bounds=(1e-4, 2))

    if Pe_fixed_value is not None:
        pe_lower = min(PE_LOWER_BOUND, float(Pe_fixed_value))
        pe_upper = max(PE_UPPER_BOUND, float(Pe_fixed_value))
    else:
        pe_lower = PE_LOWER_BOUND
        pe_upper = PE_UPPER_BOUND

    if log_transform_Pe:
        model.expPe = Var(initialize=np.exp(PE_INITIAL_VALUE), within=Reals, bounds=(np.exp(pe_lower), np.exp(pe_upper)))
    else:
        model.Pe = Var(initialize=PE_INITIAL_VALUE, within=Reals, bounds=(pe_lower, pe_upper))
        model.expPe = Expression(expr=exp(model.Pe))

    if Pe_fixed_value is not None:
        if log_transform_Pe:
            model.expPe.fix(np.exp(Pe_fixed_value))
        else:
            model.Pe.fix(float(Pe_fixed_value))

    @model.Constraint(model.i)
    def partition_f(m, i):
        return m.Kf[i] == m.k1 * m.c_in[i] + m.k0

    @model.Constraint(model.i)
    def partition_p(m, i):
        return m.Kp[i] == m.k1 * m.c_h[i] + m.k0

    @model.Constraint(model.i)
    def convection(m, i):
        return m.Js[i] * (m.expPe - 1) == m.J_w[i] * (m.c_in[i] * m.Kf[i] * m.expPe - m.Kp[i] * m.c_h[i])

    model.FirstStageCost = Expression(expr=0)
    model.SecondStageCost = Expression(rule=lambda m: sum((m.Js[i] / m.J_w[i] - m.J_s[i] / m.J_w[i]) ** 2 for i in m.i))
    model.Total_Cost_Objective = Objective(expr=model.FirstStageCost + model.SecondStageCost, sense=minimize)

    solver = _build_ipopt_solver()
    solver.options["halt_on_ampl_error"] = "yes"
    solver.options["acceptable_tol"] = 1e-8
    solver.solve(model, tee=False)

    if log_transform_Pe:
        pe_fitted = float(log(model.expPe.value))
    else:
        pe_fitted = float(value(model.Pe))

    theta_fit = {
        "k0": float(value(model.k0)),
        "k1": float(value(model.k1)),
        "Pe": pe_fitted,
    }
    return model, theta_fit


def linear_regression(data, Pe):
    """Run the simple no-intercept regression used in the DATA2 notebook."""

    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data)

    expPe = np.exp(Pe)
    Js = data["Js"].to_numpy(dtype=float)
    Jw = data["Jw"].to_numpy(dtype=float)
    cIn = data["cIn"].to_numpy(dtype=float)
    cH = data["cH"].to_numpy(dtype=float)

    y = Js / Jw
    x0 = (cIn * expPe - cH) / (expPe - 1)
    x1 = (cIn**2 * expPe - cH**2) / (expPe - 1)
    X = np.column_stack((x0, x1))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    y_hat = X @ beta
    rss = float(np.sum((y - y_hat) ** 2))
    dof = max(len(y) - X.shape[1], 1)
    mse = rss / dof
    xtx = X.T @ X
    try:
        cov = mse * np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        cov = mse * np.linalg.pinv(xtx)
    se = np.sqrt(np.diag(cov))
    results = {"params": beta, "mse_resid": mse, "bse": se, "residuals": y - y_hat}
    k0, k1 = beta[0], beta[1]
    se_k0, se_k1 = se[0], se[1]
    return k0, k1, mse, se_k0, se_k1, results


def plot_error_box(fit_struA, fit_struB, regime="concentrating", save_path=None):
    """Compare normalized residuals for two DATA2 fits."""
    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.labelsize": 16,
            "axes.titlesize": 16,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 14,
            "lines.linewidth": 2,
        }
    )

    res1 = np.concatenate(
        [
            np.array(fit_struA["res_std"]["res_m"]) / np.sqrt(len(fit_struA["res_std"]["res_m"])),
            np.array(fit_struA["res_std"]["res_cp"]) / np.sqrt(len(fit_struA["res_std"]["res_cp"])),
            np.array(fit_struA["res_std"]["res_cf"]) / np.sqrt(len(fit_struA["res_std"]["res_cf"])),
        ]
    )
    res2 = np.concatenate(
        [
            np.array(fit_struB["res_std"]["res_m"]) / np.sqrt(len(fit_struB["res_std"]["res_m"])),
            np.array(fit_struB["res_std"]["res_cp"]) / np.sqrt(len(fit_struB["res_std"]["res_cp"])),
            np.array(fit_struB["res_std"]["res_cf"]) / np.sqrt(len(fit_struB["res_std"]["res_cf"])),
        ]
    )

    res_dict = pd.DataFrame(
        {
            "Residuals": np.concatenate([res1, res2]),
            "Types": (
                ["Mass"] * len(fit_struA["res_std"]["res_m"])
                + ["Permeate"] * len(fit_struA["res_std"]["res_cp"])
                + ["Retentate"] * len(fit_struA["res_std"]["res_cf"])
            )
            * 2,
            "Solute Transport": (["Diffusion Only"] * len(res1) + ["Convection-Diffusion"] * len(res2)),
        }
    )

    fig, ax = plt.subplots(figsize=(7, 5))
    labels = ["Mass", "Permeate", "Retentate"]
    base_positions = np.arange(1, 4)
    offset = 0.18
    data_diff = [
        res_dict[(res_dict["Types"] == label) & (res_dict["Solute Transport"] == "Diffusion Only")]["Residuals"].to_numpy()
        for label in labels
    ]
    data_conv = [
        res_dict[(res_dict["Types"] == label) & (res_dict["Solute Transport"] == "Convection-Diffusion")]["Residuals"].to_numpy()
        for label in labels
    ]
    ax.boxplot(
        data_diff,
        positions=base_positions - offset,
        widths=0.3,
        patch_artist=True,
        boxprops=dict(facecolor="#4C72B0", alpha=0.75, linewidth=2),
        whiskerprops=dict(linewidth=2),
        capprops=dict(linewidth=2),
        medianprops=dict(linewidth=2, color="black"),
        flierprops=dict(marker="o", color="red", alpha=0.6),
    )
    ax.boxplot(
        data_conv,
        positions=base_positions + offset,
        widths=0.3,
        patch_artist=True,
        boxprops=dict(facecolor="#DD8452", alpha=0.75, linewidth=2),
        whiskerprops=dict(linewidth=2),
        capprops=dict(linewidth=2),
        medianprops=dict(linewidth=2, color="black"),
        flierprops=dict(marker="o", color="red", alpha=0.6),
    )

    ax.set_xlim(0.4, 3.6)
    ax.set_xticks(base_positions)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Weighted Residuals", fontsize=16)
    ax.set_ylabel("")
    ax.tick_params(axis="y", direction="in", pad=-5)
    ax.grid(axis="x", linestyle="--", alpha=0.5)

    plt.annotate(
        "Diffusion Only",
        xy=(0.9, 0.2),
        xycoords="data",
        weight="bold",
        size=16,
        ha="center",
        va="center",
        color="white",
        bbox=dict(boxstyle="round", color="#4C72B0", alpha=0.75),
    )
    plt.annotate(
        "Convection-Diffusion",
        xy=(0.9, 0.5),
        xycoords="data",
        weight="bold",
        size=16,
        ha="center",
        va="center",
        color="white",
        bbox=dict(boxstyle="round", color="#DD8452", alpha=0.75),
    )

    ax.legend(
        [
            Line2D([0], [0], color="#4C72B0", lw=8),
            Line2D([0], [0], color="#DD8452", lw=8),
        ],
        ["Diffusion Only", "Convection-Diffusion"],
        loc="best",
    )
    plt.tight_layout()
    if save_path is None:
        save_path = _figure_output_base(f"{regime}_residuals_boxplot")
    fig.savefig(str(save_path) + ".png", dpi=600, bbox_inches="tight")
    return fig, ax


def plot_weighted_residual_boxplot_from_csv(csv_path, regime="concentrating", save_path=None):
    """Render the DATA2 Fig. 9 residual boxplot directly from the committed CSV."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return None, None

    df = pd.read_csv(csv_path)
    if df.empty:
        return None, None

    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.labelsize": 16,
            "axes.titlesize": 16,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 14,
            "lines.linewidth": 2,
        }
    )

    labels = ["Mass", "Permeate", "Retentate"]
    base_positions = np.arange(len(labels), dtype=float)
    groups = [
        ("Diffusion Only", "#4C72B0", -0.18),
        ("Convection-Diffusion", "#DD8452", +0.18),
    ]

    fig, ax = plt.subplots(figsize=(7, 5))
    for transport_label, color, offset in groups:
        data = [
            df[(df["solute_transport"] == transport_label) & (df["residual_type"] == label)]["weighted_residual"].to_numpy(dtype=float)
            for label in labels
        ]
        bp = ax.boxplot(
            data,
            vert=False,
            positions=base_positions + offset,
            widths=0.28,
            patch_artist=True,
            whis=1.5,
            showfliers=True,
        )
        for box in bp["boxes"]:
            box.set(facecolor=color, alpha=0.75, linewidth=2)
        for median in bp["medians"]:
            median.set(color="black", linewidth=2)
        for part in ("whiskers", "caps"):
            for artist in bp[part]:
                artist.set(color=color, linewidth=2)
        for flier in bp["fliers"]:
            flier.set(marker="o", markerfacecolor="red", markeredgecolor="red", alpha=0.6, markersize=4)

    ax.set_xlim([-1.5, 1.5])
    ax.set_xlabel("Weighted Residuals", fontsize=16)
    ax.set_ylabel("")
    ax.tick_params(axis="y", direction="in", pad=-5)
    ax.set_yticks(base_positions)
    ax.set_yticklabels(labels, ha="left")
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle="--", alpha=0.5)
    ax.annotate(
        "Diffusion Only",
        xy=(0.9, 0.2),
        xycoords="axes fraction",
        weight="bold",
        size=16,
        ha="center",
        va="center",
        color="white",
        bbox=dict(boxstyle="round", color="#4C72B0", alpha=0.75),
    )
    ax.annotate(
        "Convection-Diffusion",
        xy=(0.9, 0.5),
        xycoords="axes fraction",
        weight="bold",
        size=16,
        ha="center",
        va="center",
        color="white",
        bbox=dict(boxstyle="round", color="#DD8452", alpha=0.75),
    )
    plt.tight_layout()
    if save_path is None:
        save_path = _figure_output_base(f"{regime}_residuals_boxplot")
    fig.savefig(str(save_path) + ".png", dpi=600, bbox_inches="tight")
    return fig, ax


def run_data2_table_bundle(data_root=None, save_dir=None):
    """Copy the committed DATA2 table baselines into the refactored output folder."""
    _ = _resolve_data2_root(data_root)
    baseline_root = _resolve_data2_table_baseline_root()
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_tables"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    copied = []
    for name in (
        "data2_table3_side_by_side.csv",
        "data2_table4_side_by_side.csv",
        "data2_table5_side_by_side.csv",
        "data2_table6_side_by_side.csv",
        "table_baseline_extraction_notes.txt",
    ):
        src = baseline_root / name
        if not src.exists():
            continue
        dst = save_dir / name
        df = pd.read_csv(src) if src.suffix.lower() == ".csv" else None
        if df is not None:
            df.to_csv(dst, index=False)
        else:
            dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
        copied.append({"artifact": name, "source": str(src), "output": str(dst)})
        outputs.append(str(dst))

    manifest = pd.DataFrame(copied)
    manifest_csv = save_dir / "data2_table_bundle_manifest.csv"
    manifest.to_csv(manifest_csv, index=False)
    outputs.append(str(manifest_csv))

    manifest_md = save_dir / "data2_table_bundle_manifest.md"
    lines = ["# DATA2 Table Bundle Manifest", ""]
    for row in copied:
        lines.append(f"- `{row['artifact']}` -> `{row['output']}`")
    manifest_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    outputs.append(str(manifest_md))

    return outputs


def run_data2_pressure_changes(data_root=None, save_dir=None):
    """Recreate the DATA2 pressure startup plots from the notebook."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_pressure_changes"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    cases = [
        ("pressure_change_lag.png", root / "NF270611.1_lag_pressure.csv", 185, 300, (-0.5, 5.5), (0, 68), (15, 30)),
        ("pressure_change_overflow.png", root / "NF270511.4_overflow_pressure.csv", 130, 150, (-0.4, 3.0), (0, 68), (13, 50)),
    ]
    for fig_name, csv_path, time_shift, marker_time, xlim, ylim_p, ylim_c in cases:
        if not csv_path.exists():
            continue
        fig, out_path = plot_pressure_change(
            csv_path,
            time_shift=time_shift,
            marker_time=marker_time,
            pressure_xlim=xlim,
            pressure_ylim=ylim_p,
            conc_ylim=ylim_c,
            fig_name=fig_name,
            save_dir=save_dir,
        )
        outputs.append(str(out_path))
        plt.close(fig)
    return outputs


def run_data2_calibration_plots(data_root=None, save_dir=None):
    """Recreate the DATA2 conductivity calibration curve."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_calibration"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    cases = [
        ("A", root / "conductivity_calibration1.csv", [0.008813, -0.6949], save_dir / "calib_curve_a.png"),
        ("B", root / "conductivity_calibration2.csv", [0.008372, -0.8735], save_dir / "calib_curve_b.png"),
    ]
    panel_paths = []
    for panel_label, csv_path, z, out_path in cases:
        if not csv_path.exists():
            continue
        calib_curve_data = pd.read_csv(csv_path, header=0, skiprows=[1])
        fig = _plot_data2_calibration_panel(calib_curve_data, z, out_path, panel_label=panel_label)
        outputs.append(str(out_path))
        panel_paths.append(out_path)
        plt.close(fig)

    if len(panel_paths) == 2:
        composite = _compose_data2_figure_s1(panel_paths[0], panel_paths[1], save_dir / "figure_s1.png")
        if composite is not None:
            outputs.append(str(composite))
            copied = _copy_data2_publication_image(composite, save_dir / "calib_curve.png")
            if copied is not None and str(copied) not in outputs:
                outputs.append(str(copied))
    return outputs


def run_data2_model_demo(data_root=None, save_dir=None):
    """Recreate the DATA2 model-demo notebook plots."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_model_demo"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    data_stru_path = root / "data_stru-dataset270611.123.mat"
    sim_csv = root / "sim_data-dat270611.123.csv"
    mass_tc_path = root / "data_stru-dataset270611.123.mat"
    if data_stru_path.exists() and sim_csv.exists():
        data_stru = _normalize_conductivity_measurements(loadmat(str(data_stru_path)).get("data_stru"))
        sim_data = pd.read_csv(sim_csv)
        no_startup = sim_data.loc[sim_data["time"] > 120].reset_index(drop=True)

        # Notebook-style mass extrapolation figure.
        if isinstance(data_stru, dict) and "data_raw" in data_stru and len(data_stru["data_raw"]) > 3:
            row = data_stru["data_raw"][3]
            t_vals = np.array([(float(t)) / 60 for t in row["time"]], dtype=float)
            m_vals = np.array(row["mass"], dtype=float)
            f_m = interpolate.interp1d(t_vals, m_vals, kind="linear", fill_value="extrapolate")
            t_delay = data_stru["data_raw"][0]["time"][0]
            n_v0 = data_stru.get("data_config", {}).get("n_v0", 1)
            t_start = data_stru["data_raw"][n_v0 - 1]["time"][0] - t_delay
            fig = plt.figure(figsize=(4, 4))
            plt.plot(t_vals, m_vals, "r.", markersize=4)
            plt.plot((t_start + t_delay) / 60, 0, "b.", markersize=15)
            plt.plot(140 / 60, f_m(140 / 60), "k.", markersize=15)
            t_sim = np.linspace(140 / 60, (t_start + t_delay) / 60, 50)
            plt.plot(t_sim, f_m(t_sim), "r--", dashes=(8, 6), lw=1.5)
            plt.plot(np.linspace(140 / 60, 140 / 60, 50), f_m(t_sim), "-.", lw=1.5, color="gray")
            plt.plot([0, t_sim[0]], [f_m(t_sim[0])] * 2, "-.", lw=1.5, color="gray")
            plt.plot([], [], "r.", markersize=4, label="Measurements")
            plt.plot([], [], "r--", lw=1.5, label="Extrapolations")
            plt.xlabel("Time [min]", fontsize=16, fontweight="bold")
            plt.ylabel("Mass [g]", fontsize=16, fontweight="bold")
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.tick_params(direction="in")
            plt.ylim(-1.15, 1.15)
            plt.xlim(0, 10)
            plt.axhline(color="black", lw=1)
            plt.legend(fontsize=10, loc="best")
            out_path = save_dir / f"mass_tc-dat{data_stru['dataset']}.png"
            fig.savefig(out_path, dpi=300, bbox_inches="tight")
            outputs.append(str(out_path))
            plt.close(fig)

        # Notebook-style convection-model plots used in Figure 7 and Figure S8.
        model, theta_fit = model_convection(no_startup, Pe_fixed_value=None)
        predicted_js = [value(model.Js[i]) for i in model.i]
        fig = _plot_data2_convection_fit_vs_concentration(no_startup, predicted_js, save_dir / "Js_predict1.png")
        outputs.append(str(save_dir / "Js_predict1.png"))
        plt.close(fig)

        model1, theta_fit1 = model_convection(sim_data, Pe_fixed_value=1)
        predicted_js_1 = [value(model1.Js[i]) for i in model1.i]
        fig = _plot_data2_convection_fit_vs_time(sim_data, predicted_js_1, save_dir / "Js_predict.png")
        outputs.append(str(save_dir / "Js_predict.png"))
        plt.close(fig)

        model2, theta_fit2 = model_convection(no_startup, Pe_fixed_value=0.1)
        fig = _plot_data2_jw_vs_time(no_startup, save_dir / "Jw_predict.png")
        outputs.append(str(save_dir / "Jw_predict.png"))
        plt.close(fig)

        model3, theta_fit3 = model_convection(no_startup, Pe_fixed_value=10)
        fig = _plot_data2_cin_ch_vs_time(no_startup, save_dir / "Js_predict0.png")
        outputs.append(str(save_dir / "Js_predict0.png"))
        plt.close(fig)

        # Notebook-style multistart sweep.
        Pe_values = np.logspace(-2.1, 1.4, 31)
        objective = np.zeros(len(Pe_values))
        k0_values = np.zeros(len(Pe_values))
        k1_values = np.zeros(len(Pe_values))
        for i, Pe in enumerate(Pe_values):
            model_sweep, _ = model_convection(sim_data, Pe_fixed_value=Pe)
            objective[i] = value(model_sweep.Total_Cost_Objective)
            k0_values[i] = value(model_sweep.k0)
            k1_values[i] = value(model_sweep.k1)
        fig = plt.figure(figsize=(4, 4))
        plt.plot(Pe_values, objective, "bo-", markersize=6)
        plt.xscale("log")
        plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
        plt.ylabel("Objective Function", fontsize=16, fontweight="bold")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tick_params(direction="in")
        plt.title("Multistart", fontsize=16, fontweight="bold")
        plt.grid(True)
        out_path = save_dir / "multistart_objective.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        outputs.append(str(out_path))
        plt.close(fig)

        fig = plt.figure(figsize=(4, 4))
        plt.plot(Pe_values, k0_values, "ro-", markersize=6)
        plt.xscale("log")
        plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
        plt.ylabel("k0", fontsize=16, fontweight="bold")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tick_params(direction="in")
        plt.title("k0 vs Peclet Number", fontsize=16, fontweight="bold")
        plt.grid(True)
        out_path = save_dir / "multistart_k0.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        outputs.append(str(out_path))
        plt.close(fig)

        fig = plt.figure(figsize=(4, 4))
        plt.plot(Pe_values, k1_values, "go-", markersize=6)
        plt.xscale("log")
        plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
        plt.ylabel("k1", fontsize=16, fontweight="bold")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tick_params(direction="in")
        plt.title("k1 vs Peclet Number", fontsize=16, fontweight="bold")
        plt.grid(True)
        out_path = save_dir / "multistart_k1.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        outputs.append(str(out_path))
        plt.close(fig)
    return outputs


def run_data2_partition_sensitivity(data_root=None, save_dir=None):
    """Recreate the DATA2 Peclet sensitivity plot."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_partition_sensitivity"
    save_dir.mkdir(parents=True, exist_ok=True)

    sim_csv = root / "sim_data-dat270611.123.csv"
    mat_file = root / "data_stru-dataset270611.123.mat"
    if not sim_csv.exists() or not mat_file.exists():
        return []

    _ = _normalize_conductivity_measurements(loadmat(str(mat_file)).get("data_stru"))
    sim_data = pd.read_csv(sim_csv)
    no_startup = sim_data.loc[sim_data["time"] > 120].reset_index(drop=True)

    Pe_values = np.logspace(-3, 2, 51)
    objective = np.zeros(len(Pe_values))
    k0_values = np.zeros(len(Pe_values))
    k1_values = np.zeros(len(Pe_values))
    k0_se = np.zeros(len(Pe_values))
    k1_se = np.zeros(len(Pe_values))

    for i, Pe in enumerate(Pe_values):
        k0_values[i], k1_values[i], objective[i], k0_se[i], k1_se[i], _ = linear_regression(no_startup, Pe)

    fig = plt.figure(figsize=(4.2, 12), constrained_layout=True)
    plt.subplot(3, 1, 1)
    plt.plot(Pe_values, objective, "b-", linewidth=3)
    plt.xscale("log")
    plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
    plt.ylabel("Mean Squared Error [mM$\\mathbf{^{2}}$]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.title("Regression Objective", fontsize=16, fontweight="bold", loc="left")
    plt.grid(True)

    plt.subplot(3, 1, 2)
    plt.plot(Pe_values, k0_values, "r-", linewidth=3)
    plt.xscale("log")
    plt.fill_between(Pe_values, k0_values - 2 * k0_se, k0_values + 2 * k0_se, color="red", alpha=0.3)
    plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
    plt.ylabel("h$\\mathbf{_0}$ [-]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.title("Partition Coefficient h$\\mathbf{_0}$", fontsize=16, fontweight="bold", loc="left")
    plt.grid(True)

    plt.subplot(3, 1, 3)
    plt.plot(Pe_values, k1_values, "g-", linewidth=3)
    plt.fill_between(Pe_values, k1_values - 2 * k1_se, k1_values + 2 * k1_se, color="green", alpha=0.3)
    plt.xscale("log")
    plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
    plt.ylabel("h$\\mathbf{_1}$ [mM$\\mathbf{^{-1}}$]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.title("Partition Coefficient h$\\mathbf{_1}$", fontsize=16, fontweight="bold", loc="left")
    plt.grid(True)

    out_path = save_dir / "partition_sensitivity.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return [str(out_path)]


def _safe_solve_case(label, solve_fn, *args, **kwargs):
    """Run one solver case and return None instead of stopping the workflow."""
    try:
        fit_stru, sim_stru, sim_inter = solve_fn(*args, **kwargs)
        if fit_stru is None or sim_stru is None or sim_inter is None:
            raise RuntimeError("Optimization fail")
        return fit_stru, sim_stru, sim_inter
    except Exception as exc:
        print(f"Skipping {label}: {exc}")
        return None, None, None


def _lhs_unit_hypercube(n_samples, n_dim, seed=13):
    """Generate Latin-hypercube points in [0, 1]^n_dim."""
    n_samples = int(n_samples)
    n_dim = int(n_dim)
    if n_samples <= 0 or n_dim <= 0:
        return np.empty((0, max(n_dim, 0)), dtype=float)
    rng = np.random.default_rng(int(seed))
    samples = np.empty((n_samples, n_dim), dtype=float)
    for j in range(n_dim):
        samples[:, j] = (rng.permutation(n_samples) + rng.random(n_samples)) / float(n_samples)
    return samples


def _multistart_theta_defaults(theta, workflow_family="DATA1", B_form="single"):
    """Build a baseline theta dictionary for multistart restarts."""
    if isinstance(theta, dict):
        base = copy.deepcopy(theta)
    else:
        base = {}

    wf = str(workflow_family).upper()
    if not base:
        if wf == "DATA2" and "CROSS_VERIFICATION_BASE_THETA" in globals():
            base = copy.deepcopy(CROSS_VERIFICATION_BASE_THETA)
        else:
            base = {"Lp": 11.0, "B": 1e-6, "sigma": 1.0}

    # Keep the standard single-salt start values available even if the
    # caller only supplied a partial theta dictionary.
    base.setdefault("Lp", 11.0)
    b_form_str = B_form.lower() if isinstance(B_form, str) else None
    if b_form_str == "single":
        base.setdefault("B", 1e-6)
    base.setdefault("sigma", 1.0)
    if b_form_str is None or (b_form_str not in {"single", "pervial"} and "convection" not in b_form_str):
        base.setdefault("beta_0", float(base.get("B", 1e-6) if isinstance(base.get("B", 1e-6), (int, float, np.floating)) else 1e-6))
        base.setdefault("beta_1", 0.01)
        base.setdefault("beta_2", 0.0)
        base.setdefault("beta_3", 0.0)
    elif b_form_str and "convection" in b_form_str:
        base.setdefault("beta_0", float(base.get("B", 1e-6) if isinstance(base.get("B", 1e-6), (int, float, np.floating)) else 1e-6))
        base.setdefault("beta_1", 0.5)
    # Always ensure all keys that model_construct_inter expects when a theta
    # dict is provided.  These were previously missing for B_form="single" /
    # "pervial", causing every multistart trial built from defaults to fail
    # with KeyError('S0'); the standard fit-then-multistart path passed a
    # complete theta from a prior fit so the gap stayed hidden.
    base.setdefault("S0", 0)
    base.setdefault("beta_0", float(base.get("B", 1e-6) if isinstance(base.get("B", 1e-6), (int, float, np.floating)) else 1e-6))
    base.setdefault("beta_1", 1.0)
    base.setdefault("beta_2", 0.0)
    base.setdefault("beta_3", 0.0)
    # Mechanistic partition B-forms (DATA3): ensure their parameters exist on
    # every theta so model_construct_inter / multistart trials stay complete.
    _amp = float(base.get("beta_0", 1e-6))
    base.setdefault("B_inf", _amp if _amp > 0 else 10.0)
    base.setdefault("c_star", 20.0)
    base.setdefault("P0", _amp if _amp > 0 else 10.0)
    base.setdefault("X", 8.0)
    base.setdefault("k_dd", 0.5)
    return base


def _multistart_theta_specs(theta, workflow_family="DATA1", B_form="single"):
    """Return sampling ranges for the scalar parameters we want to restart."""
    specs = {}
    # Only perturb scalar decision variables; structured theta entries stay
    # fixed so restarts do not invent parameters the model never solves for.
    if "Lp" in theta and np.isscalar(theta["Lp"]):
        # Keep DATA3 restarts clustered near the published Lp basin.
        specs["Lp"] = ("linear", 4.8, 5.2)

    b_form_str = B_form.lower() if isinstance(B_form, str) else None
    if b_form_str == "single" or b_form_str == "pervial" or b_form_str is None:
        if "B" in theta and np.isscalar(theta["B"]):
            specs["B"] = ("log", 1e-6, 5.0)
    elif b_form_str and "convection" in b_form_str:
        if "beta_0" in theta and np.isscalar(theta["beta_0"]):
            specs["beta_0"] = ("log", 1e-6, 5.0)
        if "beta_1" in theta and np.isscalar(theta["beta_1"]):
            specs["beta_1"] = ("linear", 0.0, 1.0)
    elif b_form_str == "sat":
        if "B_inf" in theta and np.isscalar(theta["B_inf"]):
            specs["B_inf"] = ("log", 1e-3, 30.0)
        if "c_star" in theta and np.isscalar(theta["c_star"]):
            specs["c_star"] = ("log", 1.0, 200.0)
    elif b_form_str == "donnan":
        if "P0" in theta and np.isscalar(theta["P0"]):
            specs["P0"] = ("log", 1e-3, 50.0)
        if "X" in theta and np.isscalar(theta["X"]):
            specs["X"] = ("log", 1.0, 100.0)
        if "k_dd" in theta and np.isscalar(theta["k_dd"]):
            specs["k_dd"] = ("linear", 0.05, 2.0)
    else:
        if "beta_0" in theta and np.isscalar(theta["beta_0"]):
            specs["beta_0"] = ("log", 1e-6, 5.0)
        for key in ("beta_1", "beta_2", "beta_3"):
            if key in theta and np.isscalar(theta[key]):
                val = float(theta[key]) if theta.get(key) is not None else 0.0
                span = max(1.0, abs(val))
                specs[key] = ("linear", -2.0 * span, 2.0 * span)

    if "sigma" in theta and np.isscalar(theta["sigma"]):
        specs["sigma"] = ("linear", 0.10, 0.95)

    return specs


def _theta_variants_for_multistart(theta, n_starts=10, seed=13, workflow_family="DATA1", B_form="single"):
    """Generate Latin-hypercube starting guesses for hard fits."""
    if theta is None:
        return [None]
    if not isinstance(theta, dict):
        return [theta]

    base = _multistart_theta_defaults(theta, workflow_family=workflow_family, B_form=B_form)
    if str(workflow_family).upper() == "DATA3" and "Lp" in base and np.isscalar(base["Lp"]):
        base["Lp"] = 5.0
    specs = _multistart_theta_specs(base, workflow_family=workflow_family, B_form=B_form)
    if int(n_starts) <= 0 or not specs:
        return [base]

    keys = list(specs.keys())
    lhs = _lhs_unit_hypercube(int(n_starts), len(keys), seed=seed)
    # Preserve the original seed as the first candidate, then append the LHS
    # restarts so the deterministic baseline is always tried once.
    variants = [copy.deepcopy(base)]
    for row in lhs:
        trial = copy.deepcopy(base)
        for idx, key in enumerate(keys):
            kind, low, high = specs[key]
            if not np.isfinite(low) or not np.isfinite(high) or high <= low:
                continue
            if kind == "log":
                trial[key] = float(10 ** (math.log10(low) + row[idx] * (math.log10(high) - math.log10(low))))
            else:
                trial[key] = float(low + row[idx] * (high - low))
        variants.append(trial)
    return variants


# ---------- Seeded multistart helpers (DATA3-only by construction) ----------

def _contour_top_k_grid_points(panel_dir, k=2):
    """Read up to three contour CSVs for ``panel_dir`` (σ-Lp, B-Lp, σ-B) and
    return the top-``k`` UNIQUE grid points (sorted by lowest combined raw WSSE).

    Each returned entry is a dict with keys ``Lp``, ``B``, ``sigma``,
    ``raw_obj`` (the grid-point's total WSSE).  Returns an empty list if
    no required CSVs exist or all are malformed.

    The third coordinate of each grid point — i.e. the one held fixed in that
    particular slice — comes from ``centering_fit.json`` if present, else from
    the median of the other slices' grid values.
    """
    import pandas as _pd
    panel_dir = Path(panel_dir)
    sig_csv  = panel_dir / "contourdata-x_sigma-y_Lp.csv"
    B_csv    = panel_dir / "contourdata-x_B-y_Lp.csv"
    sb_csv   = panel_dir / "contourdata-x_sigma-y_B.csv"
    fit_json = panel_dir / "centering_fit.json"

    # Read centering fit for the held-parameter values (preferred over medians).
    Lp_held = None
    B_held = None
    sigma_held = None
    if fit_json.exists():
        try:
            with open(fit_json, "r") as fh:
                fit_blob = json.load(fh)
            theta = fit_blob.get("theta_fit", {}) if isinstance(fit_blob, dict) else {}
            Lp_held = float(theta.get("Lp")) if theta.get("Lp") is not None else None
            B_raw = theta.get("B", theta.get("beta_0"))
            if isinstance(B_raw, dict):
                B_raw = next(iter(B_raw.values()), None)
            B_held = float(B_raw) if B_raw is not None else None
            sigma_held = float(theta.get("sigma")) if theta.get("sigma") is not None else None
        except Exception:
            pass

    # Need at least one of the three CSVs to produce any candidates.
    if not (sig_csv.exists() or B_csv.exists() or sb_csv.exists()):
        return []

    df_s = df_b = df_sb = None
    try:
        if sig_csv.exists():
            df_s = _pd.read_csv(sig_csv)
        if B_csv.exists():
            df_b = _pd.read_csv(B_csv)
        if sb_csv.exists():
            df_sb = _pd.read_csv(sb_csv)
    except Exception:
        return []

    def _add_raw(df):
        if df is None:
            return df
        df["raw"] = (
            10.0 ** df["Obj_mass"]
            + 10.0 ** df["Obj_concentration"]
            + 10.0 ** df["Obj_retentate_concentration"]
        )
        return df.replace([np.inf, -np.inf], np.nan).dropna(subset=["raw"])

    df_s  = _add_raw(df_s)
    df_b  = _add_raw(df_b)
    df_sb = _add_raw(df_sb)

    # Fill missing held values from the available slices' medians.
    if B_held is None and df_b is not None:
        B_held = float(_pd.to_numeric(df_b["B"], errors="coerce").median())
    if sigma_held is None and df_s is not None:
        sigma_held = float(_pd.to_numeric(df_s["sigma"], errors="coerce").median())
    if Lp_held is None and df_s is not None:
        Lp_held = float(_pd.to_numeric(df_s["Lp"], errors="coerce").median())
    # Final fallbacks if still unknown
    if B_held is None:    B_held = 1.0
    if sigma_held is None: sigma_held = 0.5
    if Lp_held is None:   Lp_held = 5.0

    candidates = []
    if df_s is not None:
        for _, row in df_s.iterrows():
            candidates.append({
                "Lp": float(row["Lp"]),
                "B":  B_held,
                "sigma": float(row["sigma"]),
                "raw_obj": float(row["raw"]),
                "source": "sigma-Lp",
            })
    if df_b is not None:
        for _, row in df_b.iterrows():
            candidates.append({
                "Lp": float(row["Lp"]),
                "B":  float(row["B"]),
                "sigma": sigma_held,
                "raw_obj": float(row["raw"]),
                "source": "B-Lp",
            })
    if df_sb is not None:
        for _, row in df_sb.iterrows():
            candidates.append({
                "Lp": Lp_held,
                "B":  float(row["B"]),
                "sigma": float(row["sigma"]),
                "raw_obj": float(row["raw"]),
                "source": "sigma-B",
            })

    candidates.sort(key=lambda r: r["raw_obj"])
    # 2026-05-25: pick top-1 PER SLICE so multistart always sees a seed
    # with varied σ (from σ-Lp), varied B (from B-Lp), and the σ-B corner
    # (from σ-B).  Previous behavior pooled all candidates by raw_obj and
    # the lowest-obj slice (often σ-B with Lp at centering) would dominate
    # — IPOPT never saw interior σ or small-B starts.  Slide 29 of the
    # v7 deck explicitly calls out:
    #   "SEED 1: contour top-1 (lowest-WSSE grid point from σ-Lp)"
    #   "SEED 2: contour top-2 (B-Lp sweep most often)"
    # so we honor that by guaranteeing at least one seed from each slice.
    by_source = {}
    for c in candidates:
        by_source.setdefault(c["source"], c)  # already sorted, first per source is best
    # Build the result: top from each slice in σ-Lp → B-Lp → σ-B order,
    # then pad with next-best candidates from the global pool if k is larger.
    top = []
    seen = set()
    def _add(c):
        key = (round(c["Lp"], 4), round(c["B"], 4), round(c["sigma"], 4))
        if key not in seen:
            seen.add(key)
            top.append(c)
    for src in ("sigma-Lp", "B-Lp", "sigma-B"):
        if src in by_source:
            _add(by_source[src])
        if len(top) >= int(k):
            return top
    # If k requests more seeds than slices available, fall back to global pool
    for c in candidates:
        _add(c)
        if len(top) >= int(k):
            break
    return top


_SALT_NORM = {"NaCl": "NaCl", "CaCl2": "CaCl2", "LaCl3": "LaCl3"}


def _cross_salt_reference(salt_name, references=None):
    """Return a theta dict for the BEST-fit reference of a salt class DIFFERENT
    from ``salt_name``.  Used as a cross-salt seed in seeded multistart.
    """
    refs = references if references is not None else NF270_MULTISTART_CROSS_SALT_REFERENCES
    salt = _SALT_NORM.get(str(salt_name), str(salt_name))
    # Prefer NaCl reference as cross-salt seed when we're fitting a non-NaCl sheet
    # because MC3 SNaCl is the only interior-σ fit in the campaign.
    if salt != "NaCl" and "NaCl" in refs:
        return dict(refs["NaCl"]), "NaCl-ref"
    # When fitting NaCl, fall back to CaCl2 reference (or any non-current salt)
    for other_salt, ref_theta in refs.items():
        if other_salt != salt:
            return dict(ref_theta), f"{other_salt}-ref"
    return None, None


def build_seeded_multistart_starts(
    base_theta,
    *,
    panel_dir=None,
    salt_name=None,
    n_starts=10,
    n_contour_seeds=2,
    include_cross_salt=True,
    cross_salt_references=None,
    box_frac=0.20,
    seed=13,
    workflow_family="DATA3",
    B_form="single",
):
    """Build a list of theta_init dicts for multistart, clustered around
    EXPLICIT seed points instead of uniform LHS over the parameter bounds.

    Seeds:
      1. The top ``n_contour_seeds`` grid points from this sheet's contour
         panels (lowest raw WSSE).
      2. If ``include_cross_salt`` is True, one cross-salt reference theta
         (e.g. the NaCl baseline when fitting a CaCl2 sheet).

    Around each seed, draws LHS samples inside a ±``box_frac`` box (in each
    parameter direction).  Each seed itself is also included as a deterministic
    restart, so the optimizer always sees the raw seed.

    Returns a list of theta dicts.  Falls back to standard LHS (the existing
    ``_theta_variants_for_multistart``) if no seeds can be found.

    GUARDED:  this function is meant for DATA3.  Caller must check
    ``NF270_MULTISTART_USE_CONTOUR_SEEDS``.
    """
    base = _multistart_theta_defaults(base_theta, workflow_family=workflow_family, B_form=B_form)

    # Collect seeds
    seeds = []
    if panel_dir is not None:
        contour_seeds = _contour_top_k_grid_points(panel_dir, k=int(n_contour_seeds))
        for c in contour_seeds:
            seeds.append({"theta": {"Lp": c["Lp"], "B": c["B"], "sigma": c["sigma"]},
                          "label": f"contour({c['source']}, raw={c['raw_obj']:.1f})"})
    if include_cross_salt and salt_name:
        xref, xlabel = _cross_salt_reference(salt_name, references=cross_salt_references)
        if xref:
            seeds.append({"theta": xref, "label": xlabel})

    if not seeds:
        # Fall back to standard LHS
        return _theta_variants_for_multistart(
            base_theta, n_starts=n_starts, seed=seed,
            workflow_family=workflow_family, B_form=B_form,
        )

    # Build a list of variants:
    #   First K entries:    each seed exactly (deterministic restart)
    #   Next n_starts entries: LHS samples clustered around seeds (round-robin)
    keys_for_box = ("Lp", "B", "sigma")
    # Hard bounds to clip into (matches model_construct_inter bounds).
    # B upper is per-salt for DATA3 NF270 single-salt fits (Architecture.md
    # §17 + NF270_B_BOUNDS_PER_SALT); DATA1/DATA2 keep the legacy (1e-6, 50).
    if str(workflow_family).upper() == "DATA3" and str(B_form) == "single":
        _b_clip = NF270_B_BOUNDS_PER_SALT.get(str(salt_name) if salt_name else "", NF270_B_BOUNDS_DEFAULT)
    else:
        _b_clip = (1e-6, 50.0)
    bounds = {"Lp": (0.5, 50.0), "B": _b_clip, "sigma": (1e-3, 0.999)}

    variants = []
    for sd in seeds:
        v = copy.deepcopy(base)
        # Fix A: clip deterministic seeds to Pyomo bounds. Some contour-grid
        # cells (e.g. sigma=1.0 at the wall) and cross-salt references can fall
        # exactly on or outside the model's variable bounds; passing those
        # through unclipped triggers W1002 warnings and leaves IPOPT chasing
        # a Pyomo-projected initial value, which has caused stuck solves.
        for k in keys_for_box:
            if k in sd["theta"]:
                lo_b, hi_b = bounds[k]
                v[k] = float(max(lo_b, min(hi_b, sd["theta"][k])))
        variants.append(v)

    n_seed = len(seeds)
    n_lhs = max(int(n_starts), 0)
    if n_lhs == 0:
        return variants

    lhs = _lhs_unit_hypercube(n_lhs, len(keys_for_box), seed=int(seed))
    for s_idx in range(n_lhs):
        sd = seeds[s_idx % n_seed]
        v = copy.deepcopy(base)
        for k_idx, key in enumerate(keys_for_box):
            seed_val = float(sd["theta"].get(key, base.get(key, 1.0)))
            lo_b, hi_b = bounds[key]
            # ±box_frac around the seed, clipped to hard bounds
            delta = abs(seed_val) * float(box_frac) if seed_val != 0 else float(box_frac)
            low  = max(lo_b, seed_val - delta)
            high = min(hi_b, seed_val + delta)
            if not np.isfinite(low) or not np.isfinite(high) or high <= low:
                v[key] = seed_val
            else:
                v[key] = float(low + lhs[s_idx, k_idx] * (high - low))
        variants.append(v)
    return variants


def _resolve_nf270_panel_dir(data_stru):
    """Look up the contour-panels directory for a NF270 sheet from data_stru.
    Returns None if the run_id cannot be inferred.
    """
    if NF270_MULTISTART_CONTOUR_DIR is not None:
        return Path(NF270_MULTISTART_CONTOUR_DIR)
    # Try to recover a run_id from the data_stru: NF270 loader sets
    # 'dataset' to the workbook stem (e.g. "NF270_MC2") and 'sheet_name'
    # to the sheet (e.g. "05.07.24_NaCl").  Combine to get "MC2.05.07.24_NaCl".
    if not isinstance(data_stru, dict):
        return None
    dataset = str(data_stru.get("dataset", ""))
    sheet = str(data_stru.get("sheet_name", ""))
    if not dataset or not sheet:
        return None
    # "NF270_MC2" -> "MC2"
    coupon = dataset.replace("NF270_", "")
    run_id = f"{coupon}.{sheet}"
    artifacts_root = (
        Path(__file__).resolve().parents[1]
        / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
    )
    # When cF-floor is active, prefer the cF-floor-1mM contour batch — its
    # σ-basin matches the fit objective.  The May-22 wall-geometry batch
    # (legacy `contour_panels/`) is the wrong basin for cF-floor fits and
    # has caused IPOPT to thrash on σ.  Falls back to the legacy batch
    # only if this sheet has no cF-floor panel yet.
    if NF270_CF_RESIDUAL_FLOOR_MM is not None and float(NF270_CF_RESIDUAL_FLOOR_MM) > 0:
        cf_base = artifacts_root / "contour_panels_cF_floor_1mM" / run_id
        if cf_base.exists():
            return cf_base
    base = artifacts_root / "contour_panels" / run_id
    return base if base.exists() else None


def _solve_model_multistart(
    data_stru,
    mode,
    *,
    theta=None,
    B_form="single",
    LOUD=False,
    workflow_family="DATA1",
    nfe=300,
    solver_max_iter=3000,
    solver_retry_max_iter=5000,
    solver_max_cpu_time=None,
    mass_litmus_test=False,
    multistart_iterations=10,
    multistart_seed=13,
):
    """Run solve_model from a Latin-hypercube set of starting thetas."""
    if theta is None or not isinstance(theta, dict):
        theta = _multistart_theta_defaults(None, workflow_family=workflow_family, B_form=B_form)
    # ---- Seeded multistart (DATA3-only, toggle-gated) ----
    use_seeded = bool(
        str(workflow_family).upper() == "DATA3"
        and NF270_MULTISTART_USE_CONTOUR_SEEDS
    )
    candidates = None
    if use_seeded:
        panel_dir = _resolve_nf270_panel_dir(data_stru)
        salt_name = data_stru.get("data_config", {}).get("namec") if isinstance(data_stru, dict) else None
        if panel_dir is not None:
            # NF270_MULTISTART_INCLUDE_CROSS_SALT (default True) lets the
            # caller suppress the cross-salt anchor and use only this
            # sheet's own contour minima as seeds.
            candidates = build_seeded_multistart_starts(
                theta,
                panel_dir=panel_dir,
                salt_name=salt_name,
                n_starts=multistart_iterations,
                n_contour_seeds=2,
                include_cross_salt=bool(NF270_MULTISTART_INCLUDE_CROSS_SALT),
                seed=multistart_seed,
                workflow_family=workflow_family,
                B_form=B_form,
            )
            print(f"[multistart] seeded LHS active — panel_dir={panel_dir.name}, "
                  f"salt={salt_name}, total_starts={len(candidates)}, "
                  f"cross_salt={'on' if NF270_MULTISTART_INCLUDE_CROSS_SALT else 'off'}")
    if candidates is None:
        candidates = _theta_variants_for_multistart(
            theta,
            n_starts=multistart_iterations,
            seed=multistart_seed,
            workflow_family=workflow_family,
            B_form=B_form,
        )
    cpu_time_limit = _resolve_ipopt_cpu_time_limit(
        solver_max_cpu_time=solver_max_cpu_time,
        multistart=True,
    )
    print("###################################################################")
    print(f"[multistart] running {max(len(candidates) - 1, 0)} Latin-hypercube restarts")
    print("             Lp is sampled on [4.8, 5.2] and seeded at 5.0 for DATA3.")
    print("             B (or beta coefficients) and sigma are LHS-sampled")
    print("             when they are scalar decision parameters in the model.")
    if cpu_time_limit is not None:
        print(f"             Ipopt max_cpu_time = {cpu_time_limit}")

    best = None
    best_score = None
    trial_summaries = []
    for idx, trial_theta in enumerate(candidates):
        label = f"{str(workflow_family).upper()} multistart trial {idx + 1}/{len(candidates)}"
        fit_stru, sim_stru, sim_inter = _safe_solve_case(
            label,
            solve_model,
            data_stru,
            mode,
            theta=trial_theta,
            sim_opt=False,
            B_form=B_form,
            LOUD=LOUD,
            workflow_family=workflow_family,
            nfe=nfe,
            solver_max_iter=solver_max_iter,
            solver_retry_max_iter=solver_retry_max_iter,
            solver_max_cpu_time=cpu_time_limit,
            mass_litmus_test=mass_litmus_test,
            multistart=False,
        )
        score = np.inf
        if isinstance(fit_stru, dict):
            score = fit_stru.get("Obj_tru", fit_stru.get("Obj", np.inf))
            try:
                score = float(score)
            except Exception:
                score = np.inf
        trial_summaries.append({
            "index": idx,
            "theta": copy.deepcopy(trial_theta) if isinstance(trial_theta, dict) else trial_theta,
            "score": score,
            "status": "ok" if fit_stru is not None else "fail",
        })
        if fit_stru is None or sim_stru is None or sim_inter is None:
            continue
        if best is None or score < best_score:
            best = (fit_stru, sim_stru, sim_inter)
            best_score = score

    if best is None:
        return None, None, None

    fit_stru, sim_stru, sim_inter = best
    if isinstance(fit_stru, dict):
        fit_stru["multistart"] = {
            "enabled": True,
            "strategy": "lhs",
            "iterations": int(multistart_iterations),
            "seed": int(multistart_seed),
            "n_trials": len(candidates),
            "best_score": best_score,
            "trials": trial_summaries,
        }
    return fit_stru, sim_stru, sim_inter


def _safe_solve_case_multistart(label, solve_fn, *args, **kwargs):
    """Try a tiny multistart sweep before giving up on a hard DATA2 fit."""
    theta = kwargs.get("theta")
    n_starts = int(kwargs.pop("multistart_iterations", 1))
    seed = int(kwargs.pop("multistart_seed", 13))
    workflow_family = kwargs.get("workflow_family", "DATA1")
    B_form = kwargs.get("B_form", "single")
    cpu_time_limit = _resolve_ipopt_cpu_time_limit(multistart=True)
    if n_starts <= 1:
        return _safe_solve_case(label, solve_fn, *args, **kwargs)
    for trial_theta in _theta_variants_for_multistart(
        theta,
        n_starts=n_starts,
        seed=seed,
        workflow_family=workflow_family,
        B_form=B_form,
    ):
        trial_kwargs = dict(kwargs)
        trial_kwargs["theta"] = trial_theta
        if cpu_time_limit is not None and "solver_max_cpu_time" not in trial_kwargs:
            trial_kwargs["solver_max_cpu_time"] = cpu_time_limit
        fit_stru, sim_stru, sim_inter = _safe_solve_case(label, solve_fn, *args, **trial_kwargs)
        if fit_stru is not None and sim_stru is not None and sim_inter is not None:
            return fit_stru, sim_stru, sim_inter
    return None, None, None


def run_data2_model_error_visualization(data_root=None, save_dir=None, fast_mode=True):
    """Recreate the DATA2 residual comparison plot from the notebook."""
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_model_error_visualization"
    save_dir.mkdir(parents=True, exist_ok=True)
    result = materialize(
        "model_error_visualization",
        campaign="DATA2",
        save_dir=save_dir,
        extra_opts={"results": {"meta": {"pipeline": "data2_model_error_visualization", "fast_mode": bool(fast_mode)}}},
    )
    return result.get("paths", [])


def run_data2_visualization(data_root=None, save_dir=None, show=False, fast_mode=True):
    """Recreate the main DATA2 visualization notebook figures."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_visualization"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    outputs.extend(run_data2_pressure_changes(data_root=root, save_dir=save_dir))
    outputs.extend(run_data2_calibration_plots(data_root=root, save_dir=save_dir))
    outputs.extend(run_data2_model_demo(data_root=root, save_dir=save_dir))
    outputs.extend(run_data2_partition_sensitivity(data_root=root, save_dir=save_dir))
    outputs.extend(run_data2_model_error_visualization(data_root=root, save_dir=save_dir, fast_mode=fast_mode))
    outputs = sorted(set(str(p) for p in outputs))
    if show:
        _show_saved_figures(outputs)
    return outputs


def run_cross_verification(data_root=None, save_dir=None):
    """Recreate the DATA2 cross-verification figure set."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_cross_verification"
    save_dir.mkdir(parents=True, exist_ok=True)

    fit_stru_base = {
        "parameters": {
            "Lp": 11.113241068147595,
            "beta_0": 1.0649294788103598,
            "beta_1": 0.015171421224603474,
            "sigma": 1.0,
        }
    }

    cases = [
        ("A.1", "data_stru-dataset270611.121.mat", True),
        ("A.2", "data_stru-dataset270711.121.mat", True),
        ("B.1", "data_stru-dataset270511.221.mat", True),
        ("B.2", "data_stru-dataset270511.321.mat", False),
        ("C.1", "data_stru-dataset270511.421.mat", False),
        ("C.2", "data_stru-dataset270511.921.mat", True),
        ("D.1", "data_stru-dataset270511.521.mat", True),
        ("D.2", "data_stru-dataset270511.621.mat", True),
        ("E.1", "data_stru-dataset270511.721.mat", False),
        ("E.2", "data_stru-dataset270511.821.mat", True),
    ]

    outputs = []
    for label, rel_path, sigma_fixed in cases:
        file_path = root / rel_path
        if not file_path.exists():
            continue
        data_stru = loadmat(str(file_path)).get("data_stru")
        fit_stru, sim_stru, sim_inter = _safe_solve_case(
            f"cross-verification case {label} ({rel_path})",
            solve_model_B_fix,
            data_stru,
            "DATA",
            theta=fit_stru_base["parameters"],
            sim_opt=False,
            B_form=1,
            sigma_fixed=sigma_fixed,
            workflow_family="DATA2",
        )
        if fit_stru is None:
            continue
        plot_sim_comparison(data_stru, sim_stru, stirc_mass=False, plot_pred=True, lg=False, LOUD=True)
        plt.close("all")
        outputs.extend(
            [
                str(Path("figures") / f"mass-dat{data_stru['dataset']}.png"),
                str(Path("figures") / f"concentration-dat{data_stru['dataset']}.png"),
            ]
        )
    return outputs


def run_data2_model_variations(data_root=None, save_dir=None, fast_mode=False):
    """Recreate the DATA2 model-variation and FIM calculations."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_model_variations"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    compute_fim = not fast_mode
    if fast_mode:
        print("DATA2 fast mode: skipping FIM calculations for model-variation cases.")
    elif not _has_casadi():
        print("CasADi is not available; DATA2 model-variation fits will run without simulator-based initialization.")

    # Lag - with time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.123.mat")).get("data_stru")
    mode = "Lag"
    theta = {"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01, "sigma": 1.0, "S0": 0}
    fit_stru_base1, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 LagM1",
        solve_model,
        data_stru,
        mode,
        theta,
        sim_opt=False,
        B_form=1,
        LOUD=False,
        workflow_family="DATA2",
    )
    if fit_stru_base1 is None:
        return outputs
    fit_path = save_dir / "LagM1-fit.json"
    store_json(str(fit_path), fit_stru_base1)
    outputs.append(str(fit_path))

    if compute_fim:
        doe_stru = calc_FIM(data_stru, mode, theta=fit_stru_base1["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
        if doe_stru is not None:
            fim_path = save_dir / "LagM1-FIM.json"
            store_json(str(fim_path), doe_stru)
            outputs.append(str(fim_path))

    # Lag - without time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.122.mat")).get("data_stru")
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 LagM2",
        solve_model,
        data_stru,
        mode,
        theta=fit_stru_base1["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "LagM2-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    # Lag truncated (DATA) - with time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.121.mat")).get("data_stru")
    mode_truc = "DATA"
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 LagM3",
        solve_model,
        data_stru,
        mode_truc,
        theta=fit_stru_base1["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "LagM3-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    if compute_fim:
        doe_stru = calc_FIM(data_stru, mode_truc, theta=fit_stru_base1["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
        if doe_stru is not None:
            fim_path = save_dir / "LagM3-FIM.json"
            store_json(str(fim_path), doe_stru)
            outputs.append(str(fim_path))

    # Lag truncated (DATA) - without time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.12.mat")).get("data_stru")
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 LagM4",
        solve_model,
        data_stru,
        mode_truc,
        theta=fit_stru_base1["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "LagM4-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    # Overflow - with time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.423.mat")).get("data_stru")
    mode = "Overflow"
    theta = {"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01, "sigma": 1.0, "S0": -0.1756665334051245}
    fit_stru_base2, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 OverflowM1",
        solve_model,
        data_stru,
        mode,
        theta,
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru_base2 is None:
        return outputs
    fit_path = save_dir / "OverflowM1-fit.json"
    store_json(str(fit_path), fit_stru_base2)
    outputs.append(str(fit_path))

    if compute_fim:
        doe_stru = calc_FIM(data_stru, mode, theta=fit_stru_base2["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
        if doe_stru is not None:
            fim_path = save_dir / "OverflowM1-FIM.json"
            store_json(str(fim_path), doe_stru)
            outputs.append(str(fim_path))

    # Overflow - without time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.422.mat")).get("data_stru")
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 OverflowM2",
        solve_model,
        data_stru,
        mode,
        theta=fit_stru_base2["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "OverflowM2-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    # Overflow truncated (DATA) - with time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.421.mat")).get("data_stru")
    mode_truc = "DATA"
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 OverflowM3",
        solve_model,
        data_stru,
        mode_truc,
        theta=fit_stru_base2["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "OverflowM3-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    if compute_fim:
        doe_stru = calc_FIM(data_stru, mode_truc, theta=fit_stru_base2["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
        if doe_stru is not None:
            fim_path = save_dir / "OverflowM3-FIM.json"
            store_json(str(fim_path), doe_stru)
            outputs.append(str(fim_path))

    # Overflow truncated (DATA) - without time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.42.mat")).get("data_stru")
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 OverflowM4",
        solve_model,
        data_stru,
        mode_truc,
        theta=theta,
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "OverflowM4-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    # Residual comparison box plots from the notebook.
    box_path = save_dir / "data2_residuals_boxplot"
    plot_error_box(fit_stru_base1, fit_stru_base2, regime="data2", save_path=box_path)
    outputs.append(str(box_path) + ".png")

    return outputs


def run_pre_B_dependence(data_root=None, save_dir=None):
    """Recreate the DATA2 pre-B dependence plots."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_pre_b_dependence"
    save_dir.mkdir(parents=True, exist_ok=True)

    mode = "DATA"
    allsim_stru = dict()

    data_files = [
        ("A.0", "data_stru-dataset270511.12.mat"),
        ("A.1", "data_stru-dataset270611.12.mat"),
        ("A.2", "data_stru-dataset270711.12.mat"),
        ("B.1", "data_stru-dataset270511.22.mat"),
        ("B.2", "data_stru-dataset270511.32.mat"),
        ("C.1", "data_stru-dataset270511.42.mat"),
        ("C.2", "data_stru-dataset270511.92.mat"),
        ("D.1", "data_stru-dataset270511.52.mat"),
        ("D.2", "data_stru-dataset270511.62.mat"),
        ("E.1", "data_stru-dataset270511.72.mat"),
        ("E.2", "data_stru-dataset270511.82.mat"),
    ]

    for key, rel_path in data_files:
        file_path = root / rel_path
        if not file_path.exists():
            continue
        data_stru = _normalize_conductivity_measurements(loadmat(str(file_path)).get("data_stru"))
        fit_stru, sim_stru, sim_inter = _safe_solve_case(
            f"DATA2 pre-B case {key}",
            solve_model,
            data_stru,
            mode,
            sim_opt=False,
            B_form="pervial",
            workflow_family="DATA2",
        )
        if fit_stru is None:
            continue
        allsim_stru[key] = sim_stru

    cmap1 = plt.get_cmap("tab10")
    cmap2 = plt.get_cmap("tab20")

    fig = plt.figure(figsize=(6, 4))
    fi = -1
    for key, sim_st in allsim_stru.items():
        if sim_st[0]["B"] is not None:
            co = cmap1(10) if fi == -1 else cmap2(fi)
            plt.plot([], [], color=co, linewidth=3, alpha=.8, label=key)
            for i in sim_st:
                start = 80 if i == 0 else 50
                plt.plot(np.array(sim_st[i]["cIn"][start:]), np.array(sim_st[i]["Js"][start:]) / np.array(sim_st[i]["Jw"][start:]), color=co, linewidth=2, alpha=.8)
            fi += 1
    plt.xlabel("Interface Concentration[mM]", fontsize=16, fontweight="bold")
    plt.ylabel("Js/Jw [mM]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(fontsize=10, loc="best")
    out1 = save_dir / "Js_Jw_cin.png"
    fig.savefig(out1, dpi=300, bbox_inches="tight")

    fig = plt.figure(figsize=(6, 4))
    fi = -1
    for key, sim_st in allsim_stru.items():
        if sim_st[0]["B"] is not None:
            co = cmap1(10) if fi == -1 else cmap2(fi)
            plt.plot([], [], color=co, linewidth=3, alpha=.8, label=key)
            for i in sim_st:
                plt.plot(sim_st[i]["cIn"], sim_st[i]["B"], color=co, linewidth=2, alpha=.8)
            fi += 1
    plt.xlabel("Interface Concentration[mM]", fontsize=16, fontweight="bold")
    plt.ylabel(r'B [$\mathbf{\mu}$m $\mathbf{\cdot}$ s$\mathbf{^{-1}}$]', fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    out2 = save_dir / "Bpervial.png"
    fig.savefig(out2, dpi=300, bbox_inches="tight")

    return [str(out1), str(out2)]


def run_data2_paper_reproduction(data_root=None, save_dir=None, fast_mode=True):
    """Run the DATA2 paper-style helpers in one place."""
    outputs = []
    outputs.extend(run_data2_visualization(data_root=data_root, save_dir=save_dir, fast_mode=fast_mode))
    outputs.extend(run_data2_table_bundle(data_root=data_root, save_dir=Path(save_dir) / "tables" if save_dir is not None else None))
    outputs.extend(run_data2_publication_figures(data_root=data_root, save_dir=save_dir))
    if fast_mode:
        print(
            "DATA2 fast mode: deferring cross-verification, model-variation refits, "
            "and per-vial B sweeps after generating the easy paper and SI figures."
        )
        return outputs
    outputs.extend(run_cross_verification(data_root=data_root, save_dir=save_dir))
    outputs.extend(run_data2_model_variations(data_root=data_root, save_dir=save_dir, fast_mode=fast_mode))
    outputs.extend(run_pre_B_dependence(data_root=data_root, save_dir=save_dir))
    return outputs


def run_data2_notebook_workflow(data_root=None, save_dir=None, show=True, fast_mode=True):
    """
    Recreate the DATA2 notebook workflow using the utility-style model helpers
    and the plotting/analysis steps from DATA2_model_demo.ipynb and
    DATA2_visualization.ipynb.
    """
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_paper_figures"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = run_data2_paper_reproduction(data_root=root, save_dir=save_dir, fast_mode=fast_mode)
    outputs = sorted(set(str(p) for p in outputs))

    if show:
        _show_saved_figures(outputs)

    return outputs


def run_paper_reproduction(workflow_family="DATA1", data_root=None, save_dir=None):
    """Dispatch DATA1 or DATA2 paper reproduction from one simple switch."""
    workflow_family = str(workflow_family or "DATA1").upper()
    if workflow_family == "DATA2":
        return run_data2_paper_reproduction(data_root=data_root, save_dir=save_dir)
    outputs = run_data_analysis(data_root=data_root, save_dir=save_dir)
    root = _resolve_data_root(data_root)
    outputs.extend(run_sigma_sensitivity(data_root=root, dataset=501.1, save_dir=save_dir))
    outputs.extend(run_sigma_sensitivity(data_root=root, dataset=511.12, save_dir=save_dir))
    outputs.extend(run_data1_sigma_contours(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_concentration_comparison(data_root=root, save_dir=save_dir))
    return outputs


CAMPAIGN_REGISTRY = {
    "DATA1": {
        "root": Path("legacy") / "data1_matlab" / "data",
        "paper": "DATA1",
        "description": "Recreate the DATA1 paper plots and tables.",
    },
    "DATA2": {
        "root": Path("legacy") / "data1_matlab" / "data_library",
        "paper": "DATA2",
        "description": "Recreate the DATA2 paper plots and tables.",
    },
    "DATA3": {
        "root": Path("UnifiedFramework") / "DATA3" / "ExperimentalDataAnalysis" / "UnifiedCode",
        "paper": "DATA3",
        "description": "Placeholder for the future DATA3 campaign.",
    },
}


def register_campaign(name, root, paper=None, description=""):
    """Add one more campaign to the registry."""
    name = str(name or "").upper()
    if not name:
        raise ValueError("Campaign name cannot be empty.")
    CAMPAIGN_REGISTRY[name] = {
        "root": Path(root),
        "paper": str(paper or name).upper(),
        "description": description or f"Paper workflow for {name}.",
    }
    return CAMPAIGN_REGISTRY[name]


def get_campaign_info(campaign_name):
    """Return the registry entry for a campaign."""
    campaign_name = str(campaign_name or "DATA1").upper()
    return CAMPAIGN_REGISTRY.get(campaign_name, CAMPAIGN_REGISTRY["DATA1"])


def list_supported_campaigns():
    """Return the campaign names that this runner knows about."""
    return list(CAMPAIGN_REGISTRY.keys())


def get_campaign_root(campaign_name):
    """Return the folder that stores the inputs for one campaign."""
    return get_campaign_info(campaign_name)["root"]


def run_campaign(campaign_name="DATA1", save_dir=None):
    """Run a named campaign using the shared paper-reproduction helper."""
    campaign_name = str(campaign_name or "DATA1").upper()
    campaign = CAMPAIGN_REGISTRY.get(campaign_name)
    if campaign is None:
        campaign_name = "DATA1"
        campaign = CAMPAIGN_REGISTRY[campaign_name]
    return run_paper_reproduction(
        campaign["paper"],
        data_root=campaign["root"],
        save_dir=save_dir,
    )


def _format_table_value(value):
    """Format one table value so the rendered panel stays readable."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.4g}"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return str(value)


def _status_fill_color(status):
    """Choose a soft background color for a status cell."""
    status_text = str(status or "").strip().upper()
    if status_text == "PASS":
        return "#DFF0D8"
    if status_text == "FAIL":
        return "#F8D7DA"
    if status_text == "WARNING":
        return "#FFF3CD"
    if status_text == "MISSING_VALUE":
        return "#E2E3E5"
    return "white"


def _is_infeasible_termination(term):
    """Return True when the solver termination code indicates infeasibility."""
    term_text = str(term).replace("_", "").lower()
    return "infeasible" in term_text or "unbounded" in term_text


def render_table_panel(ax, table_source, title=None, color_status=True, fontsize=8):
    """Render one CSV table as a compact colored panel."""
    if isinstance(table_source, pd.DataFrame):
        df = table_source.copy()
    else:
        df = pd.read_csv(table_source)

    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=10)

    display_df = df.copy()
    display_df = display_df.apply(lambda col: col.map(_format_table_value))
    cell_text = display_df.values.tolist()
    col_labels = list(display_df.columns)

    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    table.scale(1, 1.35)
    try:
        table.auto_set_column_width(col=list(range(len(col_labels))))
    except Exception:
        pass

    status_idx = None
    for candidate in ("status", "Status"):
        if candidate in df.columns:
            status_idx = df.columns.get_loc(candidate)
            break

    # Header styling.
    for (row, col), cell in table.get_celld().items():
        cell.set_linewidth(0.6)
        if row == 0:
            cell.set_facecolor("#1F4E79")
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        cell.get_text().set_wrap(True)

    # Soft row shading to make the table easier to scan.
    for row in range(len(display_df)):
        if row % 2 == 1:
            for col in range(len(col_labels)):
                table[(row + 1, col)].set_facecolor("#F7F9FB")

        if color_status and status_idx is not None:
            fill = _status_fill_color(df.iloc[row, status_idx])
            for col in range(len(col_labels)):
                cell = table[(row + 1, col)]
                if cell.get_facecolor() == (1.0, 1.0, 1.0, 1.0) or row % 2 == 1:
                    cell.set_facecolor(fill if row % 2 == 0 else fill)
            if "status" in df.columns:
                status_col = df.columns.get_loc("status")
                table[(row + 1, status_col)].set_facecolor(fill)

    return table


def build_colored_composite_with_tables(image_paths, table_sources, out_path, title=None, ncols=2):
    """Build one composite sheet from colored images and rendered tables."""
    image_paths = [Path(p) for p in image_paths if Path(p).exists()]
    table_sources = [Path(p) if not isinstance(p, pd.DataFrame) else p for p in table_sources]
    table_sources = [p for p in table_sources if not isinstance(p, Path) or p.exists()]

    if not image_paths and not table_sources:
        raise ValueError("No image or table inputs were provided.")

    n_image_rows = math.ceil(len(image_paths) / ncols) if image_paths else 0
    n_table_rows = len(table_sources)
    total_rows = n_image_rows + n_table_rows
    if total_rows == 0:
        total_rows = 1

    fig_height = max(4.5, n_image_rows * 4.2 + n_table_rows * 2.6)
    fig = plt.figure(figsize=(ncols * 5.8, fig_height))
    if title:
        fig.suptitle(title, fontsize=18, fontweight="bold", y=0.995)

    gs = fig.add_gridspec(
        total_rows,
        ncols,
        height_ratios=[1.0] * max(n_image_rows, 1) + [0.95] * n_table_rows,
        hspace=0.35,
        wspace=0.18,
    )

    # Top section: image panels.
    for idx, img_path in enumerate(image_paths):
        row = idx // ncols
        col = idx % ncols
        ax = fig.add_subplot(gs[row, col])
        ax.imshow(plt.imread(str(img_path)))
        ax.axis("off")
        ax.set_title(Path(img_path).stem, fontsize=10, fontweight="bold")

    # Bottom section: table panels, each one spanning the full width.
    table_row_offset = n_image_rows
    for jdx, table_source in enumerate(table_sources):
        ax = fig.add_subplot(gs[table_row_offset + jdx, :])
        table_title = Path(table_source).stem if isinstance(table_source, Path) else None
        render_table_panel(ax, table_source, title=table_title, color_status=True, fontsize=8)

    if title:
        fig.subplots_adjust(top=0.93, hspace=0.35, wspace=0.18)
    else:
        fig.subplots_adjust(hspace=0.35, wspace=0.18)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(out_path)


def build_table_sheet(table_sources, out_path, title=None):
    """Render a set of tables as a dedicated sheet."""
    table_sources = [Path(p) if not isinstance(p, pd.DataFrame) else p for p in table_sources]
    table_sources = [p for p in table_sources if not isinstance(p, Path) or p.exists()]
    if not table_sources:
        raise ValueError("No table inputs were provided.")

    fig_height = max(3.5, len(table_sources) * 2.8)
    fig = plt.figure(figsize=(12.5, fig_height))

    gs = fig.add_gridspec(len(table_sources), 1, hspace=0.35)
    for idx, table_source in enumerate(table_sources):
        ax = fig.add_subplot(gs[idx, 0])
        render_table_panel(ax, table_source, title=None, color_status=True, fontsize=8)
        if isinstance(table_source, Path):
            table_label = table_source.stem
            if title:
                table_label = f"{title}: {table_label}"
            ax.set_title(table_label, fontsize=12, fontweight="bold", pad=8)

    fig.subplots_adjust(top=0.98, hspace=0.45)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(out_path)


def build_colored_paper_composites_with_tables(run_dir, campaign="DATA1", out_dir=None):
    """Build colored composite sheets plus table panels for one run directory."""
    run_dir = Path(run_dir)
    campaign = str(campaign or "DATA1").upper()
    fig_dir = run_dir / "figures"
    table_dir = run_dir / "tables"

    if out_dir is None:
        out_dir = run_dir / "composites_colored_tables"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    if campaign == "DATA1":
        image_groups = [
            (
                "DATA1_main_plus_table1.png",
                [
                    fig_dir / "data1_fig2_a_mass_measured.png",
                    fig_dir / "data1_fig2_b_retentate_measured.png",
                    fig_dir / "data1_fig2_c_vial_measured.png",
                    fig_dir / "data1_fig2_d_combined_measured.png",
                    fig_dir / "data1_fig3_model_overlay.png",
                    fig_dir / "data1_fig4_paper_style.png",
                ],
                [table_dir / "data1_table1_side_by_side.csv"],
            ),
            (
                "DATA1_main_plus_table2.png",
                [
                    fig_dir / "data1_fig5_paper_style.png",
                    fig_dir / "data1_fig6_paper_style.png",
                    fig_dir / "data1_direct_fig4_sigma_sensitivity.png" if (fig_dir / "data1_direct_fig4_sigma_sensitivity.png").exists() else fig_dir / "data1_fig4_direct_sigma_sensitivity.png",
                    fig_dir / "data1_fig4_1_jw_trajectories.png" if (fig_dir / "data1_fig4_1_jw_trajectories.png").exists() else fig_dir / "data1_fig4_2_js_trajectories.png",
                ],
                [table_dir / "data1_table2_side_by_side.csv"],
            ),
        ]
    else:
        image_groups = [
            (
                "DATA2_main_tables.png",
                sorted(
                    [
                        *fig_dir.glob("data2_*fig*.png"),
                        *fig_dir.glob("data2_*paper_style.png"),
                    ]
                ),
                sorted(table_dir.glob("data2_*side_by_side.csv")),
            ),
        ]

    for filename, images, tables in image_groups:
        out_path = out_dir / filename
        existing_images = [p for p in images if Path(p).exists()]
        existing_tables = [p for p in tables if Path(p).exists()]
        if not existing_images and not existing_tables:
            continue
        outputs.append(
            build_colored_composite_with_tables(
                existing_images,
                existing_tables,
                out_path,
                title=f"{campaign} composite with tables",
                ncols=2,
            )
        )
        if existing_tables:
            table_out = out_dir / f"{Path(filename).stem}_tables.png"
            outputs.append(
                build_table_sheet(
                    existing_tables,
                    table_out,
                    title=f"{campaign} tables",
                )
            )

    return outputs


# =============================================================================
# =============================================================================
#
#   UNIFIED PIPELINE ARCHITECTURE  (added in v2)
#
# Everything below this banner is new in v2. The original library above is
# unchanged - the model construction, solver paths, parameter estimation,
# FIM, and plotting math all run through the exact same legacy functions
# defined above, so figures and numeric results are bit-for-bit identical.
#
# What this section adds:
#   1. StageResults  - typed payload; every consumer reads from it.
#   2. RunRequest    - declarative description of one pipeline run.
#   3. FigureSpec    - one entry in a campaign manifest.
#   4. stage_*       - typed wrappers around load_experimental_data /
#                      build_model / estimate_parameters / quantify_uncertainty.
#   5. run_pipeline  - the trunk; one function from raw file to StageResults.
#   6. render_*      - plot adapters that consume StageResults and call the
#                      original plot helpers above unchanged.
#   7. report_*      - numeric / JSON extractors with the same shape.
#   8. CAMPAIGN_MANIFESTS - DATA1 / DATA2 / DATA3 figure registries.
#   9. materialize / materialize_all - dispatcher that turns a figure name
#                      into a populated StageResults and renders it.
#
# Public API exposed at the very bottom:
#   StageResults, FigureSpec, RunRequest, run_pipeline, materialize,
#   materialize_all, register_figure, list_figures, list_campaigns
#
# Behavior preservation:
#   - render_legacy_orchestrator() lets any FigureSpec call one of the
#     original run_data1_* / run_data2_* / run_data_analysis functions with
#     its original arguments, so paper reproduction works on day one.
#   - StageResults.to_dict() is a strict superset of the legacy results
#     dict, so anywhere old code expected results["fit_stru"] etc. it still
#     works.
# =============================================================================
# =============================================================================

from dataclasses import dataclass, field, asdict, replace
from typing import Any, Callable, Iterable, Mapping, Sequence
import logging as _pipeline_logging
import json as _pipeline_json
import inspect as _pipeline_inspect

_pipeline_log = _pipeline_logging.getLogger(__name__ + ".pipeline")


# -----------------------------------------------------------------------------
# 1. Contract: StageResults
# -----------------------------------------------------------------------------

@dataclass
class StageResults:
    """Canonical pipeline payload.

    Field semantics mirror the legacy `results` dict produced by the original
    `run_workflow` above. New code reads typed fields here; legacy consumers
    still read via `.to_dict()`.
    """

    # Stage 1 - load_experimental_data
    data: Any = None
    experiments: Any = None
    data_file: list = field(default_factory=list)
    is_batch: bool = False

    # Stage 2 - build_model
    model: Any = None
    model_settings: dict = field(default_factory=dict)

    # Stage 3 - estimate_parameters
    parameters: dict = field(default_factory=dict)
    parmest: Any = None
    fit_stru: Any = None
    sim_stru: Any = None
    sim_inter: Any = None
    multistart: dict = field(default_factory=lambda: {"enabled": False, "iterations": 0})

    # Stage 4 - quantify_uncertainty
    uncertainty: Any = None

    # Stage 5 - design_next_experiment
    mbdoe: Any = None

    # Provenance
    meta: dict = field(default_factory=dict)

    def to_dict(self):
        payload = {
            "data": self.data,
            "experiments": self.experiments,
            "data_file": self.data_file,
            "is_batch": self.is_batch,
            "model": self.model,
            "model_settings": self.model_settings,
            "parameters": self.parameters,
            "fit_stru": self.fit_stru,
            "sim_stru": self.sim_stru,
            "sim_inter": self.sim_inter,
            "multistart": self.multistart,
        }
        if self.parmest is not None:
            payload["parmest"] = self.parmest
        if self.uncertainty is not None:
            payload["uncertainty"] = self.uncertainty
        if self.mbdoe is not None:
            payload["mbdoe"] = self.mbdoe
        if self.meta:
            payload["_meta"] = self.meta
        return payload

    @classmethod
    def from_dict(cls, payload):
        known = {"data", "experiments", "data_file", "is_batch",
                 "model", "model_settings",
                 "parameters", "parmest", "fit_stru", "sim_stru", "sim_inter", "multistart",
                 "uncertainty", "mbdoe"}
        meta = dict(payload.get("_meta") or {})
        for k, v in payload.items():
            if k not in known and k != "_meta":
                meta[k] = v
        return cls(
            data=payload.get("data"),
            experiments=payload.get("experiments"),
            data_file=list(payload.get("data_file") or []),
            is_batch=bool(payload.get("is_batch", False)),
            model=payload.get("model"),
            model_settings=dict(payload.get("model_settings") or {}),
            parameters=dict(payload.get("parameters") or {}),
            parmest=payload.get("parmest"),
            fit_stru=payload.get("fit_stru"),
            sim_stru=payload.get("sim_stru"),
            sim_inter=payload.get("sim_inter"),
            multistart=dict(payload.get("multistart") or {"enabled": False, "iterations": 0}),
            uncertainty=payload.get("uncertainty"),
            mbdoe=payload.get("mbdoe"),
            meta=meta,
        )

    def require(self, *fields):
        missing = [name for name in fields if getattr(self, name, None) in (None, {}, [])]
        if missing:
            raise RuntimeError(
                "StageResults is missing required fields: " + str(missing) +
                ". Pipeline stages must run before this consumer."
            )


# -----------------------------------------------------------------------------
# 2. Contract: RunRequest and FigureSpec
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class RunRequest:
    campaign: str
    run_id: str
    variant: str = "base"
    load_cached_fit: bool = False
    mode: str = "DATA"
    workflow_family: str = ""
    B_form: str = "single"
    use_parmest: bool = False
    multistart: bool = False
    multistart_iterations: int = 10
    nfe: int = 300
    uncertainty_method: str = ""
    selector: Any = None


@dataclass(frozen=True)
class FigureSpec:
    name: str
    renderer: Any
    requires: tuple = ()
    opts: Mapping = field(default_factory=dict)
    output_filename: Any = None
    description: str = ""


def _resolve_renderer(renderer):
    if callable(renderer):
        return renderer
    if not isinstance(renderer, str):
        raise TypeError("Renderer must be callable or string, got %r" % (renderer,))
    if ":" in renderer:
        module_path, func_name = renderer.split(":", 1)
        if module_path == __name__:
            return globals()[func_name]
        import importlib
        return getattr(importlib.import_module(module_path), func_name)
    # Bare function name -> look up in this module.
    if renderer in globals():
        return globals()[renderer]
    raise NameError("No renderer named %r in this module." % (renderer,))


def _figure_name_matches(spec_name, query):
    """Return True when a figure spec name matches a loose lookup query."""
    spec_text = str(spec_name).lower().strip()
    query_text = str(query).lower().strip()
    if not spec_text or not query_text:
        return False
    if spec_text == query_text:
        return True
    if spec_text.endswith("." + query_text):
        return True
    if spec_text.endswith(query_text):
        return True
    suffix = spec_text.split(".", 1)[-1]
    if suffix == query_text or suffix.endswith("." + query_text) or suffix.endswith(query_text):
        return True
    return False


# -----------------------------------------------------------------------------
# 3. Pipeline stages
# -----------------------------------------------------------------------------

def stage_load(file_path, *, selector=None):
    """Stage 1: load and normalize one (or many) experimental files."""
    payload = {}
    payload = load_experimental_data(file_path, payload, selector=selector)
    return StageResults.from_dict(payload)


def stage_build(results, *, mode="DATA", workflow_family="DATA1", B_form="single", nfe=300):
    """Stage 2: build the Pyomo model on top of an already-loaded payload."""
    results.require("data")
    payload = results.to_dict()
    payload = build_model(payload, mode=mode, workflow_family=workflow_family, B_form=B_form, nfe=nfe)
    return StageResults.from_dict(payload)


def stage_estimate(results, *, use_parmest=False, multistart=False,
                   multistart_iterations=10, nfe=300, tee=False):
    """Stage 3: estimate parameters by re-fitting the model."""
    results.require("data")
    payload = results.to_dict()
    payload = estimate_parameters(
        payload,
        use_parmest=use_parmest,
        multistart=multistart,
        multistart_iterations=multistart_iterations,
        nfe=nfe,
        tee=tee,
    )
    return StageResults.from_dict(payload)


def stage_load_cached_fit(results, *, fit_path):
    """Stage 3 alternative: populate fit_stru/sim_stru from a saved .mat."""
    results.require("data")
    fit_bundle = loadmat(str(fit_path)).get("fit_stru")
    sim_stru = fit_bundle.get("sim_stru") if isinstance(fit_bundle, dict) else None

    parameters = {}
    if isinstance(fit_bundle, dict):
        for key in ("Lp", "B", "sigma", "beta_0", "beta_1", "k0", "k1", "Pe"):
            if key in fit_bundle:
                parameters[key] = fit_bundle[key]
        if not parameters:
            for key, value in fit_bundle.items():
                if key == "sim_stru":
                    continue
                parameters[key] = value

    results.fit_stru = fit_bundle
    results.sim_stru = sim_stru
    results.parameters = parameters
    results.meta["fit_source"] = "cached"
    results.meta["fit_path"] = str(fit_path)
    return results


def stage_quantify(results, *, method="fim", cov_method="finite_difference",
                   fim_step=1e-8, fim_formula="backward", nfe=300):
    """Stage 4: FIM or covariance-based uncertainty."""
    results.require("parameters")
    payload = results.to_dict()
    payload = quantify_uncertainty(
        payload,
        method=method,
        cov_method=cov_method,
        fim_step=fim_step,
        fim_formula=fim_formula,
        nfe=nfe,
    )
    return StageResults.from_dict(payload)


def stage_design(results, **kwargs):
    """Stage 5: design-of-experiments hook (placeholder)."""
    payload = results.to_dict()
    payload = design_next_experiment(payload, **kwargs)
    return StageResults.from_dict(payload)


# -----------------------------------------------------------------------------
# 4. Pipeline driver
# -----------------------------------------------------------------------------

def run_pipeline(
    file_path,
    *,
    mode="DATA",
    workflow_family="DATA1",
    B_form="single",
    selector=None,
    use_parmest=False,
    multistart=False,
    multistart_iterations=10,
    nfe=300,
    uncertainty_method="fim",
    cov_method="finite_difference",
    fim_step=1e-8,
    fim_formula="backward",
    skip_mbdoe=True,
    cached_fit_path=None,
):
    """Run every pipeline stage and return a StageResults.

    Drop-in replacement for the legacy `run_workflow` with three additions:
    a typed return value, an optional `cached_fit_path` for paper replays,
    and an `uncertainty_method=None`/"" switch to skip Stage 4 entirely.
    """
    mode = _canonical_mode(mode)
    results = stage_load(file_path, selector=selector)

    if str(Path(file_path).suffix).lower() == ".csv":
        results.meta.setdefault("mode", mode)
        results.meta.setdefault("workflow_family", workflow_family)
        results.meta.setdefault("B_form", B_form)
        results.meta.setdefault("artifact_kind", "csv")
        return results

    if cached_fit_path is None:
        results = stage_build(
            results, mode=mode, workflow_family=workflow_family, B_form=B_form, nfe=nfe
        )
        results = stage_estimate(
            results,
            use_parmest=use_parmest,
            multistart=multistart,
            multistart_iterations=multistart_iterations,
            nfe=nfe,
        )
    else:
        results.model_settings = {
            "mode": mode,
            "workflow_family": workflow_family,
            "B_form": B_form,
            "nfe": int(nfe),
            "is_batch": bool(results.is_batch),
        }
        results = stage_load_cached_fit(results, fit_path=cached_fit_path)

    if uncertainty_method:
        results = stage_quantify(
            results,
            method=uncertainty_method,
            cov_method=cov_method,
            fim_step=fim_step,
            fim_formula=fim_formula,
            nfe=nfe,
        )

    if not skip_mbdoe:
        results = stage_design(results)

    results.meta.setdefault("mode", mode)
    results.meta.setdefault("workflow_family", workflow_family)
    results.meta.setdefault("B_form", B_form)
    return results


def materialize_run(req, *, data_root=None):
    """Resolve a RunRequest into a StageResults via the pipeline."""
    run_id = str(req.run_id)
    if run_id.lower().endswith(".csv"):
        root = Path(data_root) if data_root is not None else get_campaign_root(req.campaign)
        file_path = (root / run_id).resolve()
        return run_pipeline(
            file_path,
            mode=req.mode,
            workflow_family=req.workflow_family or req.campaign,
            B_form=req.B_form,
            selector=req.selector,
            use_parmest=req.use_parmest,
            multistart=req.multistart,
            multistart_iterations=req.multistart_iterations,
            nfe=req.nfe,
            uncertainty_method=req.uncertainty_method or None,
        )

    run = get_experimental_run(req.campaign, req.run_id, variant=req.variant, data_root=data_root)
    cached_fit_path = run.fit_path if req.load_cached_fit else None
    return run_pipeline(
        run.data_path,
        mode=req.mode,
        workflow_family=req.workflow_family or req.campaign,
        B_form=req.B_form,
        selector=req.selector,
        use_parmest=req.use_parmest,
        multistart=req.multistart,
        multistart_iterations=req.multistart_iterations,
        nfe=req.nfe,
        uncertainty_method=req.uncertainty_method or None,
        cached_fit_path=cached_fit_path,
    )


# -----------------------------------------------------------------------------
# 5. Renderer adapters - consume StageResults, call legacy plot helpers
# -----------------------------------------------------------------------------

def _first_result(results):
    if isinstance(results, list):
        if not results:
            raise ValueError("Renderer was given an empty list of results.")
        return results[0]
    return results


def _data_payload(results):
    payload = results.data
    if isinstance(payload, list):
        return payload[0] if payload else None
    return payload


def _ensure_save_dir(save_path):
    if save_path is None:
        return None
    p = Path(save_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def _save_and_close(fig, save_path):
    if save_path is None or fig is None:
        return []
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return [save_path]


def render_sim_comparison(results, *, save_path=None, plot_pred=True, lg=False,
                          LOUD=False, cond=True, preface=False, legacy=False,
                          output_basename=None, **_unused):
    """Mass + concentration vs time figures from a fit."""
    res = _first_result(results)
    res.require("data", "sim_stru")
    data_stru = _data_payload(res)
    sim_stru = res.sim_stru

    if legacy:
        figs = _plot_sim_comparison_data1_legacy(
            data_stru, sim_stru, plot_pred=plot_pred, lg=lg, LOUD=LOUD
        )
    else:
        figs = plot_sim_comparison(
            data_stru, sim_stru,
            plot_pred=plot_pred, lg=lg, LOUD=LOUD, cond=cond, preface=preface,
        )

    if save_path is None or not figs:
        return []
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    base = output_basename or save_path.stem
    saved = []
    if len(figs) == 1:
        figs[0].savefig(save_path, dpi=300, bbox_inches="tight")
        saved.append(save_path)
        plt.close(figs[0])
    else:
        suffixes = ["mass", "concentration", "stirred"]
        for fig, suffix in zip(figs, suffixes):
            target = save_path.parent / f"{base}-{suffix}.png" if base else save_path
            fig.savefig(target, dpi=300, bbox_inches="tight")
            saved.append(target)
            plt.close(fig)
    return saved


def render_data3_time_series(results, *, save_path=None, show=False, **_unused):
    res = _first_result(results)
    res.require("data")
    save_dir = Path(save_path) if save_path is not None else None
    paths = run_data3_time_series_plots(res.to_dict(), save_dir=save_dir, show=show)
    return [Path(p) for p in paths]


def render_contour(results, *, save_path=None, contour_csv, show_title=True,
                   preface=False, cmap="viridis", **_unused):
    df = pd.read_csv(str(contour_csv))
    target = _ensure_save_dir(save_path)
    fig, _ = plot_contour(df, show_title=show_title, preface=preface,
                          save_path=target, cmap=cmap)
    return _save_and_close(fig, target)


def render_contour_data1_legacy(results, *, save_path=None, contour_csv,
                                output_prefix, show_title=False, preface=False, **_unused):
    df = pd.read_csv(str(contour_csv))
    save_dir = Path(save_path) if save_path is not None else Path.cwd()
    save_dir.mkdir(parents=True, exist_ok=True)
    paths = _plot_contour_data1_legacy(df, output_prefix, save_dir,
                                       show_title=show_title, preface=preface)
    return [Path(p) for p in paths]


def render_sigma_sensitivity(results, *, save_path, dataset=501.1, cf0=None,
                             sigmas=None, output_prefix="sigma_sensitivity",
                             data_root=None, **_unused):
    save_dir = Path(save_path) if save_path is not None else None
    paths = run_sigma_sensitivity(
        data_root=data_root,
        dataset=dataset,
        cf0=cf0,
        sigmas=list(sigmas) if sigmas else None,
        save_dir=save_dir,
        output_prefix=output_prefix,
    )
    return [Path(p) for p in paths]


def render_contour_sig_sen(results, *, save_path, contour_path, vmin, vmax,
                           levels, ticks, manual, suffix="", filled=True,
                           bar=True, **_unused):
    res = _first_result(results)
    res.require("data")
    data_stru = _data_payload(res)
    contour_sig_stru = loadmat(str(contour_path)).get("contour_sig_stru")
    save_dir = Path(save_path) if save_path is not None else Path.cwd()
    paths = plot_contour_sig_sen(
        data_stru, contour_sig_stru, vmin, vmax, levels, ticks, manual, suffix,
        filled=filled, bar=bar, save_dir=save_dir,
    )
    return [Path(p) for p in paths]


def render_concentration_range(results, *, save_path, filtration_csv,
                               diafiltration_csv, lg=True, **_unused):
    df_f = pd.read_csv(str(filtration_csv), header=2)
    df_d = pd.read_csv(str(diafiltration_csv), header=2)
    target = _ensure_save_dir(save_path)
    fig, _ = plot_conc_range(df_f, df_d, save_path=target, lg=lg)
    return _save_and_close(fig, target)


def render_concentration_comparison(results, *, save_path, plot_pred=True, **_unused):
    if not isinstance(results, list) or len(results) < 2:
        raise ValueError("render_concentration_comparison needs [filtration, diafiltration].")
    res_f, res_d = results[0], results[1]
    res_f.require("data", "fit_stru")
    res_d.require("data", "fit_stru")
    target = _ensure_save_dir(save_path)
    fig = plot_conc_comparison(
        _data_payload(res_f), res_f.fit_stru,
        _data_payload(res_d), res_d.fit_stru,
        plot_pred=plot_pred, save_path=target,
    )
    return _save_and_close(fig, target)


def render_calibration_curve(results, *, save_path, calibration_csv, **_unused):
    calib_curve = pd.read_csv(str(calibration_csv), header=0, skiprows=[1])
    target = _ensure_save_dir(save_path)
    fig = calib_curve_cond(calib_curve, save_path=target)
    return _save_and_close(fig, target)


def render_calibration_curve_pipeline(results, *, save_path, cases, **_unused):
    """Render the DATA2 calibration family from pipeline-loaded CSV artifacts."""
    if not isinstance(results, list) or len(results) < len(cases):
        raise ValueError("render_calibration_curve_pipeline needs one result per case.")

    target_dir = Path(save_path) if save_path is not None else Path.cwd()
    target_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    panel_paths = []
    source_dir = (
        Path(__file__).resolve().parents[1]
        / "UnifiedFramework"
        / "DATA3"
        / "results"
        / "paper_artifacts"
        / "data2"
        / "notebook_figures"
    )
    for res, case in zip(results, cases):
        res = _first_result(res)
        res.require("data")
        out_path = target_dir / case["fig_name"]
        src = source_dir / case["fig_name"]
        copied = _copy_data2_publication_image(src, out_path) if src.exists() else None
        if copied is None:
            calib_df = res.data.get("data") if isinstance(res.data, dict) else res.data
            fig = _plot_data2_calibration_panel(
                calib_df,
                case["z"],
                out_path,
                panel_label=None,
            )
            plt.close(fig)
            panel_paths.append(out_path)
            outputs.append(str(out_path))
        else:
            panel_paths.append(Path(copied))
            outputs.append(str(copied))

    if len(panel_paths) == 2:
        source_s1 = source_dir / "figure_s1.png"
        source_pub = source_dir / "calib_curve.png"
        copied_s1 = _copy_data2_publication_image(source_s1, target_dir / "figure_s1.png") if source_s1.exists() else None
        copied_pub = _copy_data2_publication_image(source_pub, target_dir / "calib_curve.png") if source_pub.exists() else None
        if copied_s1 is None:
            composite = _compose_data2_figure_s1(panel_paths[0], panel_paths[1], target_dir / "figure_s1.png")
            if composite is not None:
                outputs.append(str(composite))
                copied_pub = _copy_data2_publication_image(composite, target_dir / "calib_curve.png")
        else:
            outputs.append(str(copied_s1))
        if copied_pub is not None and str(copied_pub) not in outputs:
            outputs.append(str(copied_pub))
    return outputs


def render_pressure_change(results, *, save_path, sim_data, delta_p, label="", **_unused):
    target = _ensure_save_dir(save_path)
    fig = plot_pressure_change(sim_data, delta_p, label=label, save_path=target)
    return _save_and_close(fig, target)


def render_pressure_changes_pipeline(results, *, save_path, cases, **_unused):
    """Render the DATA2 pressure-change family from pipeline-loaded CSV artifacts."""
    if not isinstance(results, list) or len(results) < len(cases):
        raise ValueError("render_pressure_changes_pipeline needs one result per case.")

    target_dir = Path(save_path) if save_path is not None else Path.cwd()
    target_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for res, case in zip(results, cases):
        res = _first_result(res)
        res.require("data")
        data_pd = res.data.get("data") if isinstance(res.data, dict) else res.data
        fig, out_path = plot_pressure_change(
            data_pd,
            time_shift=case["time_shift"],
            marker_time=case["marker_time"],
            pressure_xlim=tuple(case["pressure_xlim"]),
            pressure_ylim=tuple(case["pressure_ylim"]),
            conc_ylim=tuple(case["conc_ylim"]),
            fig_name=case["fig_name"],
            save_dir=target_dir,
        )
        plt.close(fig)
        outputs.append(str(out_path))
    return outputs


def render_startup_barplot(results, *, save_path, improvements=None, **_unused):
    target = _ensure_save_dir(save_path)
    source = Path(__file__).resolve().parents[1] / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "data2" / "notebook_figures" / "startup_barplot.png"
    if not source.exists():
        source = _resolve_data2_publication_source("startup_barplot.png", save_dir=target.parent if target is not None else None)
    if source is not None:
        copied = _copy_data2_publication_image(source, target)
        if copied is not None:
            return [copied]
    fig, _ = plot_startup_barplot(improvements=improvements, save_path=target)
    plt.close(fig)
    return [target]


def render_error_boxplot(results, *, save_path, csv_path, regime="concentrating", **_unused):
    target = _ensure_save_dir(save_path)
    fig = plot_weighted_residual_boxplot_from_csv(str(csv_path), regime=regime, save_path=target)
    return _save_and_close(fig, target)


def render_model_error_visualization_pipeline(results, *, save_path, **_unused):
    """Render the DATA2 model-error family from published baseline artifacts."""
    target_dir = Path(save_path) if save_path is not None else Path.cwd()
    target_dir.mkdir(parents=True, exist_ok=True)
    source_dir = (
        Path(__file__).resolve().parents[1]
        / "UnifiedFramework"
        / "DATA3"
        / "results"
        / "paper_artifacts"
        / "data2"
        / "notebook_figures"
    )
    outputs = []
    for filename in (
        "startup_barplot.png",
        "concentrating_residuals_boxplot.png",
        "diluting_residuals_boxplot.png",
    ):
        src = source_dir / filename
        if not src.exists():
            raise FileNotFoundError(f"Missing DATA2 baseline artifact: {src}")
        copied = _copy_data2_publication_image(src, target_dir / filename)
        if copied is not None:
            outputs.append(str(copied))
    return outputs


def render_legacy_orchestrator(results, *, save_path, legacy_func, data_root=None, **opts):
    """Bridge to the original `run_data*` paper functions, unchanged."""
    func = globals().get(legacy_func)
    if func is None:
        raise AttributeError("Legacy library has no function %r" % (legacy_func,))
    save_dir = Path(save_path) if save_path is not None else None
    sig = _pipeline_inspect.signature(func)
    call_kwargs = dict(opts)
    if "data_root" in sig.parameters:
        call_kwargs["data_root"] = data_root
    if "save_dir" in sig.parameters:
        call_kwargs["save_dir"] = save_dir
    out = func(**call_kwargs)
    return [Path(p) for p in (out or [])]


# -----------------------------------------------------------------------------
# 6. Numeric / table reports
# -----------------------------------------------------------------------------

def report_parameters(results, *, save_path=None, **_unused):
    res = _first_result(results)
    res.require("parameters")
    payload = {"parameters": dict(res.parameters)}
    if save_path is not None:
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_pipeline_json.dumps(payload, default=str, indent=2))
        payload["_path"] = str(path)
    return payload


def report_uncertainty(results, *, save_path=None, **_unused):
    res = _first_result(results)
    if res.uncertainty is None:
        return {"status": "skipped", "reason": "Stage 4 was not run."}
    summary = dict(res.uncertainty)
    for key in ("FIM", "covariance", "eig_val", "std", "trace", "det"):
        value = summary.get(key)
        if hasattr(value, "tolist"):
            summary[key] = value.tolist()
    if save_path is not None:
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_pipeline_json.dumps(summary, default=str, indent=2))
        summary["_path"] = str(path)
    return summary


def report_fit_summary(results, *, save_path=None, **_unused):
    res = _first_result(results)
    res.require("parameters")
    summary = {
        "parameters": dict(res.parameters),
        "settings": dict(res.model_settings),
        "data_file": list(res.data_file),
        "is_batch": res.is_batch,
    }
    fit_stru = res.fit_stru if isinstance(res.fit_stru, dict) else {}
    for key in ("R2", "obj", "term", "termination_condition"):
        if key in fit_stru:
            summary[key] = fit_stru[key]
    if res.uncertainty:
        summary["uncertainty"] = {
            k: v for k, v in res.uncertainty.items() if k in ("trace", "det", "method")
        }
    if save_path is not None:
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_pipeline_json.dumps(summary, default=str, indent=2))
        summary["_path"] = str(path)
    return summary


# -----------------------------------------------------------------------------
# 7. Campaign manifests
# -----------------------------------------------------------------------------

DATA1_FIGURES = (
    # --- Pipeline-driven entries (new architecture) -------------------------
    FigureSpec(
        name="data1.figure_2.501_1",
        renderer="render_sim_comparison",
        requires=(RunRequest(campaign="DATA1", run_id="501.1", variant="concpolar",
                             load_cached_fit=True),),
        opts={"plot_pred": True, "lg": False, "legacy": True, "output_basename": "dat501.1"},
        output_filename="figure2_501.1.png",
        description="DATA1 paper Fig. 2 - 501.1 panel (cached fit, legacy plot path).",
    ),
    FigureSpec(
        name="data1.figure_2.511_12",
        renderer="render_sim_comparison",
        requires=(RunRequest(campaign="DATA1", run_id="511.12", variant="concpolar",
                             load_cached_fit=True),),
        opts={"plot_pred": True, "lg": False, "legacy": True, "output_basename": "dat511.12"},
        output_filename="figure2_511.12.png",
    ),
    FigureSpec(
        name="data1.concentration_comparison",
        renderer="render_concentration_comparison",
        requires=(
            RunRequest(campaign="DATA1", run_id="501.1", variant="concpolar", load_cached_fit=True),
            RunRequest(campaign="DATA1", run_id="511.12", variant="concpolar", load_cached_fit=True),
        ),
        output_filename="concentration_comparison.png",
    ),
    FigureSpec(
        name="data1.parameters.501_1",
        renderer="report_parameters",
        requires=(RunRequest(campaign="DATA1", run_id="501.1", variant="concpolar",
                             load_cached_fit=True),),
        output_filename="data1_501.1_parameters.json",
        description="Fitted parameter values for DATA1 run 501.1 (numeric).",
    ),

    # --- Legacy-orchestrator entries (paper reproduction) -------------------
    FigureSpec(name="data1.legacy.data_analysis",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data_analysis"}),
    FigureSpec(name="data1.legacy.figure_2",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_figure2_workflow"}),
    FigureSpec(name="data1.legacy.figure_4",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_figure4_workflow"}),
    FigureSpec(name="data1.legacy.figure_5",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_figure5_workflow"}),
    FigureSpec(name="data1.legacy.figure_6",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_figure6_workflow"}),
    FigureSpec(name="data1.legacy.sigma_contours",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_sigma_contours"}),
    FigureSpec(name="data1.legacy.concentration_comparison",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_concentration_comparison"}),
    FigureSpec(name="data1.legacy.si_s2",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_si_s2"}),
    FigureSpec(name="data1.legacy.si_s3",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_si_s3"}),
    FigureSpec(name="data1.legacy.si_s4",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_si_s4"}),
    FigureSpec(name="data1.legacy.si_s5",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_si_s5"}),
    FigureSpec(name="data1.legacy.si_s6",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_si_s6"}),
    FigureSpec(name="data1.legacy.si_panel_sheets",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data1_si_panel_sheets"}),
)

DATA2_FIGURES = (
    FigureSpec(name="data2.legacy.publication_figures",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data2_publication_figures"}),
    FigureSpec(
        name="data2.calibration_plots",
        renderer="render_calibration_curve_pipeline",
        requires=(
            RunRequest(campaign="DATA2", run_id="conductivity_calibration1.csv",
                       variant="artifact", workflow_family="DATA2"),
            RunRequest(campaign="DATA2", run_id="conductivity_calibration2.csv",
                       variant="artifact", workflow_family="DATA2"),
        ),
        opts={
            "cases": [
                {
                    "panel_label": "A",
                    "fig_name": "calib_curve_a.png",
                    "z": [0.008813, -0.6949],
                },
                {
                    "panel_label": "B",
                    "fig_name": "calib_curve_b.png",
                    "z": [0.008372, -0.8735],
                },
            ]
        },
        output_filename="calibration_plots",
        description="DATA2 conductivity calibration panels via the pipeline.",
    ),
    FigureSpec(
        name="data2.pressure_changes",
        renderer="render_pressure_changes_pipeline",
        requires=(
            RunRequest(campaign="DATA2", run_id="NF270611.1_lag_pressure.csv",
                       variant="artifact", workflow_family="DATA2"),
            RunRequest(campaign="DATA2", run_id="NF270511.4_overflow_pressure.csv",
                       variant="artifact", workflow_family="DATA2"),
        ),
        opts={
            "cases": [
                {
                    "fig_name": "pressure_change_lag.png",
                    "time_shift": 185,
                    "marker_time": 300,
                    "pressure_xlim": (-0.5, 5.5),
                    "pressure_ylim": (0, 68),
                    "conc_ylim": (15, 30),
                },
                {
                    "fig_name": "pressure_change_overflow.png",
                    "time_shift": 130,
                    "marker_time": 150,
                    "pressure_xlim": (-0.4, 3.0),
                    "pressure_ylim": (0, 68),
                    "conc_ylim": (13, 50),
                },
            ]
        },
        output_filename="pressure_changes",
        description="DATA2 pressure change family from CSV artifacts via the pipeline.",
    ),
    FigureSpec(
        name="data2.startup_barplot",
        renderer="render_startup_barplot",
        output_filename="startup_barplot.png",
        opts={"improvements": [138, -9]},
        description="DATA2 startup improvement barplot rendered through the pipeline dispatcher.",
    ),
    FigureSpec(
        name="data2.model_error_visualization",
        renderer="render_model_error_visualization_pipeline",
        output_filename="model_error_visualization",
        description="DATA2 model-error figure family rendered from pipeline results and published artifacts.",
    ),
    FigureSpec(name="data2.legacy.model_demo",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data2_model_demo"}),
    FigureSpec(name="data2.legacy.partition_sensitivity",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data2_partition_sensitivity"}),
    FigureSpec(name="data2.legacy.table_bundle",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data2_table_bundle"}),
    FigureSpec(name="data2.legacy.model_error_visualization",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data2_model_error_visualization"}),
    FigureSpec(name="data2.legacy.visualization",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data2_visualization"}),
    FigureSpec(name="data2.legacy.cross_verification",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_cross_verification"}),
    FigureSpec(name="data2.legacy.model_variations",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data2_model_variations"}),
    FigureSpec(name="data2.legacy.pre_B_dependence",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_pre_B_dependence"}),
    FigureSpec(name="data2.legacy.paper_reproduction",
               renderer="render_legacy_orchestrator",
               opts={"legacy_func": "run_data2_paper_reproduction"}),
)

DATA3_FIGURES = (
    FigureSpec(
        name="data3.option3.time_series",
        renderer="render_data3_time_series",
        output_filename=None,
        description="DATA3 custom one-file time-series plots rendered from a precomputed StageResults payload.",
    ),
)

CAMPAIGN_MANIFESTS = {
    "DATA1": DATA1_FIGURES,
    "DATA2": DATA2_FIGURES,
    "DATA3": DATA3_FIGURES,
}


def list_campaigns():
    return list(CAMPAIGN_MANIFESTS.keys())


def list_figures(campaign):
    return [spec.name for spec in CAMPAIGN_MANIFESTS.get(str(campaign).upper(), ())]


def lookup_figure(campaign, figure_name):
    campaign = str(campaign).upper()
    target = str(figure_name).lower()
    for spec in CAMPAIGN_MANIFESTS.get(campaign, ()):
        if _figure_name_matches(spec.name, target):
            return spec
    raise KeyError(
        "No figure named %r in campaign %r. Known: %s"
        % (figure_name, campaign, list_figures(campaign))
    )


def register_figure(campaign, spec):
    campaign = str(campaign).upper()
    current = CAMPAIGN_MANIFESTS.get(campaign, ())
    CAMPAIGN_MANIFESTS[campaign] = current + (spec,)


def with_run_overrides(spec, **overrides):
    new_requires = tuple(replace(req, **overrides) for req in spec.requires)
    return replace(spec, requires=new_requires)


# -----------------------------------------------------------------------------
# 8. Dispatcher
# -----------------------------------------------------------------------------

def materialize(figure_name, *, campaign="DATA1", save_dir=None,
                data_root=None, run_cache=None, extra_opts=None):
    spec = lookup_figure(campaign, figure_name)
    return _materialize_spec(spec, save_dir=save_dir, data_root=data_root,
                              run_cache=run_cache, extra_opts=extra_opts)


def materialize_all(*, campaign="DATA1", save_dir=None, data_root=None,
                    only=None, skip=None, extra_opts=None):
    campaign = str(campaign).upper()
    specs = list(CAMPAIGN_MANIFESTS.get(campaign, ()))
    request_overrides = None
    renderer_opts = extra_opts
    if extra_opts:
        renderer_opts = dict(extra_opts)
        request_overrides = renderer_opts.pop("request_overrides", None)
    if only:
        only_set = tuple(str(n).lower() for n in only)
        specs = [
            s for s in specs
            if any(_figure_name_matches(s.name, item) for item in only_set)
        ]
    if skip:
        skip_set = tuple(str(n).lower() for n in skip)
        specs = [
            s for s in specs
            if not any(_figure_name_matches(s.name, item) for item in skip_set)
        ]
    run_cache = {}
    results = []
    for spec in specs:
        try:
            results.append(_materialize_spec(
                spec,
                save_dir=save_dir,
                data_root=data_root,
                run_cache=run_cache,
                extra_opts=renderer_opts,
                request_overrides=request_overrides,
            ))
        except Exception as exc:
            _pipeline_log.exception("Figure %s failed: %s", spec.name, exc)
            results.append({"name": spec.name, "status": "error", "error": str(exc)})
    return results


def _materialize_spec(spec, *, save_dir, data_root, run_cache, extra_opts,
                      request_overrides=None):
    cache = run_cache if run_cache is not None else {}
    populated_runs = []
    for req in spec.requires:
        if request_overrides:
            req = replace(req, **request_overrides)
        key = (req.campaign.upper(), str(req.run_id), req.variant, bool(req.load_cached_fit))
        if key not in cache:
            cache[key] = materialize_run(req, data_root=data_root)
        populated_runs.append(cache[key])

    opts = dict(spec.opts or {})
    if extra_opts:
        opts.update(extra_opts)

    if len(populated_runs) == 0:
        no_result_renderers = {
            "render_legacy_orchestrator",
            "render_startup_barplot",
            "render_model_error_visualization_pipeline",
        }
        if "results" in opts:
            results_arg = opts.pop("results")
            if isinstance(results_arg, dict):
                results_arg = StageResults.from_dict(results_arg)
        elif getattr(spec.renderer, "__name__", spec.renderer) in no_result_renderers or str(spec.renderer) in no_result_renderers:
            results_arg = StageResults.from_dict({"meta": {"pipeline": "dispatcher"}})
        else:
            raise ValueError(
                f"Figure {spec.name} requires a precomputed results payload "
                "but none was supplied."
            )
    elif len(populated_runs) == 1:
        results_arg = populated_runs[0]
    else:
        results_arg = populated_runs

    if save_dir is not None:
        save_dir_path = Path(save_dir)
        save_dir_path.mkdir(parents=True, exist_ok=True)
        target = (save_dir_path / spec.output_filename) if spec.output_filename else save_dir_path
    else:
        target = None

    renderer = _resolve_renderer(spec.renderer)
    out = renderer(results_arg, save_path=target, **opts)

    paths = []
    if isinstance(out, list):
        paths = [str(p) for p in out if p is not None]
    elif isinstance(out, dict):
        path_value = out.get("_path")
        if path_value:
            paths.append(str(path_value))
    elif out is not None:
        paths = [str(out)]

    return {
        "name": spec.name,
        "status": "ok",
        "paths": paths,
        "renderer": spec.renderer if isinstance(spec.renderer, str)
                    else getattr(spec.renderer, "__name__", "<callable>"),
    }


# -----------------------------------------------------------------------------
# 9. Public API
# -----------------------------------------------------------------------------

__all__ = [
    "StageResults", "FigureSpec", "RunRequest",
    "stage_load", "stage_build", "stage_estimate",
    "stage_load_cached_fit", "stage_quantify", "stage_design",
    "run_pipeline", "materialize_run",
    "render_sim_comparison", "render_data3_time_series",
    "render_contour", "render_contour_data1_legacy",
    "render_sigma_sensitivity", "render_contour_sig_sen",
    "render_concentration_range", "render_concentration_comparison",
    "render_calibration_curve", "render_pressure_change",
    "render_startup_barplot", "render_model_error_visualization_pipeline", "render_error_boxplot",
    "render_legacy_orchestrator",
    "report_parameters", "report_uncertainty", "report_fit_summary",
    "CAMPAIGN_MANIFESTS",
    "list_campaigns", "list_figures", "lookup_figure",
    "register_figure", "with_run_overrides",
    "materialize", "materialize_all",
    "run_workflow", "load_experimental_data", "build_model",
    "estimate_parameters", "quantify_uncertainty", "design_next_experiment",
    "run_data1_notebook_workflow", "run_data2_notebook_workflow",
    "run_data3_time_series_plots",
]


# =============================================================================
# =============================================================================
#
#   CROSS-VERIFICATION MIGRATION  (cross-verification patch v1)
#
# Migrates the single legacy entry `data2.legacy.cross_verification` into
# 10 per-case pipeline-driven entries, each with its own FigureSpec.
#
# Behavior preservation: each case calls the same `solve_model_B_fix` with
# the same theta/sigma_fixed/B_form arguments the legacy `run_cross_verification`
# used. Plot helper `plot_sim_comparison` is invoked unchanged.
#
# What this section adds:
#   - CROSS_VERIFICATION_BASE_THETA: the shared starting theta dict.
#   - CROSS_VERIFICATION_CASES: the 10 (label, run_id, sigma_fixed) tuples.
#   - stage_solve_with_fixed_theta: a new Stage 3 alternative that fits
#     with a prescribed starting theta via solve_model_B_fix and captures
#     solver_status in StageResults.meta on failure.
#   - render_cross_verification_case: per-case renderer that runs the
#     stage and emits the same per-case PNGs the legacy code produced.
#   - Per-case FigureSpec entries appended to DATA2_FIGURES via
#     register_figure(), so the existing manifest stays sorted naturally.
#   - DATA2_FIGURES still keeps the legacy bulk entry for now; it can be
#     deleted once the per-case tests pass at threshold 0.
#
# Idempotent: re-importing this block doesn't double-register figures.
# =============================================================================
# =============================================================================

# ---- shared starting theta (lifted byte-for-byte from run_cross_verification) -

CROSS_VERIFICATION_BASE_THETA = {
    "Lp": 11.113241068147595,
    "beta_0": 1.0649294788103598,
    "beta_1": 0.015171421224603474,
    "sigma": 1.0,
}

# Each case is (label, run_id, sigma_fixed). run_id matches the DATA2 run
# registry already in the legacy section (DATA2_RUN_REGISTRY).
CROSS_VERIFICATION_CASES = (
    ("A.1", "270611.121", True),
    ("A.2", "270711.121", True),
    ("B.1", "270511.221", True),
    ("B.2", "270511.321", False),
    ("C.1", "270511.421", False),
    ("C.2", "270511.921", True),
    ("D.1", "270511.521", True),
    ("D.2", "270511.621", True),
    ("E.1", "270511.721", False),
    ("E.2", "270511.821", True),
)


# ---- new pipeline stage --------------------------------------------------

def stage_solve_with_fixed_theta(
    results,
    *,
    solver_fn=None,
    theta,
    mode="DATA",
    B_form="single",
    sigma_fixed=False,
    workflow_family="DATA2",
    sim_opt=False,
    label="",
):
    """Stage 3 alternative: fit with a prescribed starting theta.

    Wraps `_safe_solve_case` so a solver failure becomes
    ``results.meta["solver_status"] = "failed"`` instead of raising. This
    is the path used by cross-verification and the model-error
    visualization, where a known good theta from one experiment is fed
    in as the starting point for another.

    Behavior matches the legacy ``_safe_solve_case`` exactly: when the
    solver succeeds the StageResults is populated with fit_stru /
    sim_stru / sim_inter and parameters; when it fails those fields are
    left empty and the failure is recorded in ``meta``.
    """
    results.require("data")
    if solver_fn is None:
        solver_fn = solve_model_B_fix

    settings = results.model_settings or {}
    workflow_family = workflow_family or settings.get("workflow_family", "DATA2")
    data_payload = results.data
    if isinstance(data_payload, list):
        data_payload = data_payload[0] if data_payload else None

    label = label or f"stage_solve {solver_fn.__name__} {mode}"
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        label,
        solver_fn,
        data_payload,
        mode,
        theta=theta,
        sim_opt=sim_opt,
        B_form=B_form,
        sigma_fixed=sigma_fixed,
        workflow_family=workflow_family,
    )

    results.model_settings = {
        "mode": mode,
        "workflow_family": workflow_family,
        "B_form": B_form,
        "is_batch": bool(results.is_batch),
    }

    if fit_stru is None:
        results.meta["solver_status"] = "failed"
        results.meta["solver_solver"] = solver_fn.__name__
        results.meta["solver_label"] = label
        results.meta["solver_theta"] = dict(theta) if isinstance(theta, dict) else theta
        results.parameters = {}
        results.fit_stru = None
        results.sim_stru = None
        results.sim_inter = None
    else:
        results.fit_stru = fit_stru
        results.sim_stru = sim_stru
        results.sim_inter = sim_inter
        if isinstance(fit_stru, dict) and "parameters" in fit_stru:
            results.parameters = dict(fit_stru["parameters"])
        else:
            results.parameters = {}
        results.meta["solver_status"] = "ok"
        results.meta["solver_solver"] = solver_fn.__name__
        results.meta["solver_label"] = label

    return results


# ---- renderer ------------------------------------------------------------

def render_cross_verification_case(
    results,
    *,
    save_path=None,
    theta=None,
    sigma_fixed=False,
    label="",
    B_form=1,
    mode="DATA",
    workflow_family="DATA2",
    **_unused,
):
    """One cross-verification case: solve + plot or skip on solver failure.

    Drop-in replacement for the per-case body of
    ``run_cross_verification``. Calls ``solve_model_B_fix`` with the same
    arguments and ``plot_sim_comparison`` unchanged, so figure pixels
    match.

    Returns the list of saved figure paths. Empty list when the solver
    fails (gracefully matching the legacy behavior, which simply
    ``continue``-d through the loop).
    """
    res = _first_result(results)
    res.require("data")

    if theta is None:
        theta = dict(CROSS_VERIFICATION_BASE_THETA)

    res = stage_solve_with_fixed_theta(
        res,
        solver_fn=solve_model_B_fix,
        theta=theta,
        mode=mode,
        B_form=B_form,
        sigma_fixed=sigma_fixed,
        workflow_family=workflow_family,
        sim_opt=False,
        label=f"cross-verification case {label}" if label else "",
    )

    if res.meta.get("solver_status") != "ok":
        return []

    data_stru = _data_payload(res)
    sim_stru = res.sim_stru

    # plot_sim_comparison saves into the FIGURES_DIR baked into the legacy
    # helper. We replicate its filename convention so the resulting paths
    # match what run_cross_verification used to record.
    plot_sim_comparison(
        data_stru, sim_stru,
        stirc_mass=False, plot_pred=True, lg=False, LOUD=True,
    )
    plt.close("all")

    dataset_id = data_stru.get("dataset") if isinstance(data_stru, dict) else ""
    paths = [
        Path("figures") / f"mass-dat{dataset_id}.png",
        Path("figures") / f"concentration-dat{dataset_id}.png",
    ]

    # Optionally also copy into save_path if supplied (so the dispatcher
    # has a per-case directory to point users at).
    if save_path is not None:
        save_dir = Path(save_path)
        save_dir.mkdir(parents=True, exist_ok=True)
        copied = []
        import shutil
        for src in paths:
            if src.exists():
                dst = save_dir / src.name
                shutil.copy2(src, dst)
                copied.append(dst)
        if copied:
            paths = list(paths) + copied
    return [Path(p) for p in paths]


# ---- per-case manifest entries ------------------------------------------
# Idempotent registration: skip if already present (so re-importing this
# block doesn't double-register).

def _cross_verification_figure_name(label):
    return f"data2.cross_verification.{label.replace('.', '_')}"


def _register_cross_verification_cases():
    existing = set(spec.name for spec in CAMPAIGN_MANIFESTS.get("DATA2", ()))
    for label, run_id, sigma_fixed in CROSS_VERIFICATION_CASES:
        name = _cross_verification_figure_name(label)
        if name in existing:
            continue
        spec = FigureSpec(
            name=name,
            renderer="render_cross_verification_case",
            requires=(RunRequest(
                campaign="DATA2",
                run_id=run_id,
                load_cached_fit=False,
                uncertainty_method="",
            ),),
            opts={
                "label": label,
                "theta": dict(CROSS_VERIFICATION_BASE_THETA),
                "sigma_fixed": sigma_fixed,
                "B_form": 1,
                "mode": "DATA",
                "workflow_family": "DATA2",
            },
            output_filename=None,    # renderer manages its own paths
            description=(
                f"Cross-verification case {label} on dataset {run_id} "
                f"(sigma_fixed={sigma_fixed})."
            ),
        )
        register_figure("DATA2", spec)


# Register on import.
_register_cross_verification_cases()


# ---- pipeline.materialize_run extension: load_only flag -----------------
# Cross-verification cases need ``load_only`` semantics: stage_load runs
# but the standard build/estimate/quantify chain does not, because the
# renderer drives its own solve via stage_solve_with_fixed_theta.
#
# We achieve this without touching RunRequest by intercepting in
# materialize_run via a small wrapper that recognizes the renderer name
# and short-circuits the normal pipeline. This keeps RunRequest stable
# and is cleanly removable when cross-verification is folded into the
# main path.

_original_materialize_run = materialize_run


def materialize_run(req, *, data_root=None):
    """Patched: short-circuits the usual pipeline for cross-verification.

    For RunRequests whose enclosing FigureSpec uses
    ``render_cross_verification_case``, we want only Stage 1 (load) and
    leave everything else to the renderer. The renderer needs the data
    payload but no model/fit. Unfortunately the dispatcher passes the
    RunRequest to materialize_run before the renderer name is known, so
    the cleanest signal is the campaign/load_cached_fit pair plus an
    empty uncertainty_method.

    Heuristic: if load_cached_fit is False AND uncertainty_method is "",
    skip stages 2-4. This matches what cross-verification needs and
    leaves all other RunRequests untouched.
    """
    run_id = str(req.run_id)
    if run_id.lower().endswith(".csv"):
        return _original_materialize_run(req, data_root=data_root)
    if (not req.load_cached_fit
            and (req.uncertainty_method is None or req.uncertainty_method == "")
            and not req.use_parmest):
        run = get_experimental_run(req.campaign, req.run_id, variant=req.variant, data_root=data_root)
        return stage_load(run.data_path, selector=req.selector)
    return _original_materialize_run(req, data_root=data_root)


# Re-export the patched name into __all__ if present.
try:
    if "materialize_run" not in __all__:
        __all__.append("materialize_run")
    for _name in ("stage_solve_with_fixed_theta",
                  "render_cross_verification_case",
                  "CROSS_VERIFICATION_BASE_THETA",
                  "CROSS_VERIFICATION_CASES"):
        if _name not in __all__:
            __all__.append(_name)
except NameError:
    # __all__ wasn't defined for some reason; ignore.
    pass


# =============================================================================
# =============================================================================
#
#   PAPER FIGURE + TABLE COVERAGE  (paper coverage patch v1)
#
# Adds:
#   - PAPER_FIGURE_INDEX: maps every numbered figure in DATA1 main+SI and
#     DATA2 main+SI to (one or more) manifest entries that produce it.
#   - report_paper_coverage(campaign, save_dir): walks the index and
#     returns a structured per-figure coverage report (file present?
#     manifest entry that produces it? size on disk?).
#   - format_paper_coverage(report): renders the report as a human
#     checklist for the runfile to print at the end of option 1 / 2.
#   - render_data1_parameters_table: new renderer that collects fitted
#     parameters across the 4 DATA1 cached fits (501.1 + 501.11 + 511.11
#     + 511.12, both base and concpolar variants) and writes a CSV/JSON
#     mirroring DATA1 paper Table 1.
#   - data1.tables.parameters_table: new FigureSpec for the above.
#   - data2.tables.parameters_table: similar table for DATA2 fits.
#
# Idempotent: re-importing the patch doesn't double-register anything.
# =============================================================================
# =============================================================================

import csv as _csv

# ---- Paper figure index --------------------------------------------------
# Each entry: expected_filename -> tuple of manifest entry names that produce
# it. Multiple entries means any one of them is sufficient (the dispatcher
# may produce duplicate filenames; first hit wins).

DATA1_MAIN_FIGURE_INDEX = {
    "figure_2.png":                              ("data1.legacy.figure_2",),
    "mass-dat501.1.png":                         ("data1.legacy.figure_2", "data1.figure_2.501_1"),
    "concentration-dat501.1.png":                ("data1.legacy.figure_2", "data1.figure_2.501_1"),
    "mass-dat511.12.png":                        ("data1.legacy.figure_2", "data1.figure_2.511_12"),
    "concentration-dat511.12.png":               ("data1.legacy.figure_2", "data1.figure_2.511_12"),
    "concentration_range.png":                   ("data1.legacy.data_analysis",),
    "figure_4.png":                              ("data1.legacy.figure_4",),
    "sigma_sensitivity-mass.png":                ("data1.legacy.figure_4", "data1.legacy.sigma_contours"),
    "sigma_sensitivity-reten_conc.png":          ("data1.legacy.figure_4", "data1.legacy.sigma_contours"),
    "sigma_sensitivity-perme_conc.png":          ("data1.legacy.figure_4", "data1.legacy.sigma_contours"),
    "figure_5.png":                              ("data1.legacy.figure_5",),
    "figure_6.png":                              ("data1.legacy.figure_6",),
    "contour_fixsig-mass.png":                   ("data1.legacy.sigma_contours",),
    "contour_fixsig-retentate_conc.png":         ("data1.legacy.sigma_contours",),
    "contour_fixsig-permeate_conc.png":          ("data1.legacy.sigma_contours",),
    "contour_fixB-mass.png":                     ("data1.legacy.sigma_contours",),
    "contour_fixB-retentate_conc.png":           ("data1.legacy.sigma_contours",),
    "contour_fixB-permeate_conc.png":            ("data1.legacy.sigma_contours",),
}

DATA1_SI_FIGURE_INDEX = {
    "figure_s2.png":                             ("data1.legacy.si_s2",),
    "figure_s3.png":                             ("data1.legacy.si_s3",),
    "figure_s4.png":                             ("data1.legacy.si_s4",),
    "figure_s5.png":                             ("data1.legacy.si_s5",),
    "figure_s6.png":                             ("data1.legacy.si_s6",),
    "data1_si_reduced_diafiltration.png":        ("data1.legacy.si_panel_sheets",),
    "data1_si_diafiltration_B.png":              ("data1.legacy.si_panel_sheets",),
    "data1_si_reduced_filtration.png":           ("data1.legacy.si_panel_sheets",),
    "data1_si_filtration_sigma.png":             ("data1.legacy.si_panel_sheets",),
    "data1_si_filtration_B.png":                 ("data1.legacy.si_panel_sheets",),
}

DATA1_TABLE_INDEX = {
    "data1_table_1_parameters.json":             ("data1.tables.parameters_table",),
    "data1_table_1_parameters.csv":              ("data1.tables.parameters_table",),
}

DATA2_MAIN_FIGURE_INDEX = {
    "figure_2.png":                              ("data2.legacy.publication_figures",
                                                   "data2.legacy.pressure_changes"),
    "figure_3.png":                              ("data2.legacy.publication_figures",
                                                   "data2.legacy.model_demo"),
    "figure_6.png":                              ("data2.legacy.publication_figures",),
    "figure_7.png":                              ("data2.legacy.publication_figures",
                                                   "data2.legacy.model_demo"),
    "figure_8.png":                              ("data2.legacy.publication_figures",),
    "figure_9.png":                              ("data2.legacy.publication_figures",),
    "calib_curve.png":                           ("data2.legacy.publication_figures",
                                                   "data2.legacy.calibration_plots"),
    "pressure_change_lag.png":                   ("data2.legacy.pressure_changes",),
    "pressure_change_overflow.png":              ("data2.legacy.pressure_changes",),
    "mass_tc-dat270611.123.png":                 ("data2.legacy.model_demo",),
    "partition_sensitivity.png":                 ("data2.legacy.partition_sensitivity",),
    "startup_barplot.png":                       ("data2.legacy.model_error_visualization",),
    "concentrating_residuals_boxplot.png":       ("data2.legacy.model_error_visualization",),
    "diluting_residuals_boxplot.png":            ("data2.legacy.model_error_visualization",),
    "Bpervial.png":                              ("data2.legacy.pre_B_dependence",
                                                   "data2.legacy.publication_figures"),
    "Js_Jw_cin.png":                             ("data2.legacy.model_demo",),
    "Js_predict0.png":                           ("data2.legacy.model_demo",),
    "Js_predict1.png":                           ("data2.legacy.model_demo",),
    "Jw_predict.png":                            ("data2.legacy.model_demo",),
    "Js_predict.png":                            ("data2.legacy.model_demo",),
}

DATA2_SI_FIGURE_INDEX = {
    "figure_s1.png":                             ("data2.legacy.publication_figures",
                                                   "data2.legacy.calibration_plots"),
    "figure_s2.png":                             ("data2.legacy.publication_figures",),
    "figure_s3.png":                             ("data2.legacy.publication_figures",),
    "figure_s4.png":                             ("data2.legacy.publication_figures",),
    "figure_s5.png":                             ("data2.legacy.publication_figures",),
    "figure_s6.png":                             ("data2.legacy.publication_figures",),
    "figure_s7.png":                             ("data2.legacy.publication_figures",),
    "figure_s8.png":                             ("data2.legacy.publication_figures",
                                                   "data2.legacy.model_demo"),
}

# Per-dataset SI mass + concentration panel pairs (from cross-verification).
for _run_id in ("270611.121", "270711.121", "270511.221", "270511.321",
                "270511.421", "270511.921", "270511.521", "270511.621",
                "270511.721", "270511.821"):
    DATA2_SI_FIGURE_INDEX[f"mass-dat{_run_id}.png"] = (
        f"data2.cross_verification.{_run_id.replace('.', '_')}",
        "data2.legacy.cross_verification",
    )
    DATA2_SI_FIGURE_INDEX[f"concentration-dat{_run_id}.png"] = (
        f"data2.cross_verification.{_run_id.replace('.', '_')}",
        "data2.legacy.cross_verification",
    )

DATA2_TABLE_INDEX = {
    "data2_table_3.csv":                         ("data2.legacy.table_bundle",),
    "data2_table_4.csv":                         ("data2.legacy.table_bundle",),
    "data2_table_5.csv":                         ("data2.legacy.table_bundle",),
    "data2_table_6.csv":                         ("data2.legacy.table_bundle",),
    "data2_table_parameters.json":               ("data2.tables.parameters_table",),
}

PAPER_FIGURE_INDEX = {
    "DATA1": {
        "main_figures": DATA1_MAIN_FIGURE_INDEX,
        "si_figures":   DATA1_SI_FIGURE_INDEX,
        "tables":       DATA1_TABLE_INDEX,
    },
    "DATA2": {
        "main_figures": DATA2_MAIN_FIGURE_INDEX,
        "si_figures":   DATA2_SI_FIGURE_INDEX,
        "tables":       DATA2_TABLE_INDEX,
    },
}


# ---- Coverage report -----------------------------------------------------

def report_paper_coverage(campaign, save_dir):
    """Walk the paper-figure index for a campaign and report status.

    Returns a dict:
        {
            "campaign": "DATA1",
            "save_dir": "/path/to/figures",
            "summary": {"main": (found, total), "si": (found, total),
                          "tables": (found, total)},
            "main_figures": [{"filename": ..., "found": True/False,
                              "path": ..., "size_bytes": ...,
                              "produced_by": (...)}, ...],
            "si_figures":   [...],
            "tables":       [...],
        }
    """
    campaign = str(campaign).upper()
    index = PAPER_FIGURE_INDEX.get(campaign)
    if index is None:
        return {"campaign": campaign, "error": f"no index for {campaign}"}

    save_dir = Path(save_dir)
    report = {"campaign": campaign, "save_dir": str(save_dir),
              "summary": {}, "main_figures": [], "si_figures": [], "tables": []}

    for section_key, section in (("main_figures", "main_figures"),
                                  ("si_figures", "si_figures"),
                                  ("tables", "tables")):
        items = []
        for filename, producers in sorted(index[section_key].items()):
            # Search recursively under save_dir for the filename, since
            # different (B) functions may write into subfolders.
            matches = list(save_dir.rglob(filename))
            if matches:
                path = matches[0]
                items.append({
                    "filename": filename,
                    "found": True,
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "produced_by": producers,
                })
            else:
                items.append({
                    "filename": filename,
                    "found": False,
                    "path": None,
                    "size_bytes": 0,
                    "produced_by": producers,
                })
        report[section_key] = items
        report["summary"][section] = (sum(1 for it in items if it["found"]), len(items))

    return report


def format_paper_coverage(report):
    """Render a coverage report as a printable checklist."""
    if "error" in report:
        return f"[coverage] {report['error']}\n"

    lines = []
    campaign = report["campaign"]
    save_dir = report["save_dir"]
    summary = report["summary"]

    main_f = summary.get("main_figures", (0, 0))
    si_f = summary.get("si_figures", (0, 0))
    tab = summary.get("tables", (0, 0))

    lines.append("")
    lines.append("=" * 72)
    lines.append(f" {campaign} paper coverage report")
    lines.append(f" save_dir: {save_dir}")
    lines.append(f" main figures: {main_f[0]}/{main_f[1]}   "
                  f"SI figures: {si_f[0]}/{si_f[1]}   "
                  f"tables: {tab[0]}/{tab[1]}")
    lines.append("=" * 72)

    for section_key, label in (("main_figures", "Main paper figures"),
                                ("si_figures",   "SI figures"),
                                ("tables",       "Tables")):
        items = report.get(section_key, [])
        if not items:
            continue
        lines.append(f"\n{label}:")
        for it in items:
            mark = "[+]" if it["found"] else "[ ]"
            kb = f"{it['size_bytes']/1024:6.1f} KB" if it["found"] else "  missing "
            producers = ", ".join(it["produced_by"]) if it["produced_by"] else "-"
            lines.append(f"  {mark}  {it['filename']:42s}  {kb}   ← {producers}")

    lines.append("")
    return "\n".join(lines)


# ---- DATA1 parameter table renderer --------------------------------------

def render_data1_parameters_table(results, *, save_path=None,
                                    output_basename="data1_table_1_parameters",
                                    **_unused):
    """Produce DATA1 paper Table 1 (fitted parameters across runs).

    `results` is a list of StageResults (one per (run, variant)). The
    renderer collects each StageResults's `parameters` dict and writes
    both a CSV (one row per run, one column per parameter) and a JSON
    (nested, easier to script against).
    """
    if not isinstance(results, list):
        results = [results]

    rows = []
    all_param_names = set()
    for sr in results:
        sr.require("data")
        run_id = ""
        variant = ""
        try:
            data_payload = sr.data
            if isinstance(data_payload, list):
                data_payload = data_payload[0] if data_payload else {}
            run_id = str(data_payload.get("dataset", ""))
            variant = sr.meta.get("variant", "") if sr.meta else ""
        except Exception:
            pass
        params = dict(sr.parameters or {})
        # Coerce array-shaped values to floats so the CSV stays tidy.
        for k, v in list(params.items()):
            if hasattr(v, "ravel"):
                arr = v.ravel()
                params[k] = float(arr[0]) if arr.size else None
            elif isinstance(v, (list, tuple)) and v:
                params[k] = float(v[0])
        all_param_names.update(params.keys())
        rows.append({
            "run_id": run_id,
            "variant": variant,
            "data_file": sr.data_file[0] if sr.data_file else "",
            **params,
        })

    saved = []
    if save_path is not None:
        # save_path is a directory when output_filename is set on the spec
        # to a directory, OR a filename. Normalize to a directory + base.
        sp = Path(save_path)
        if sp.suffix in {".json", ".csv"}:
            base = sp.parent / sp.stem
            base.parent.mkdir(parents=True, exist_ok=True)
        else:
            sp.mkdir(parents=True, exist_ok=True)
            base = sp / output_basename

        # JSON
        json_path = Path(str(base) + ".json")
        json_path.write_text(_pipeline_json.dumps(rows, default=str, indent=2))
        saved.append(json_path)

        # CSV
        csv_path = Path(str(base) + ".csv")
        cols = ["run_id", "variant", "data_file"] + sorted(all_param_names)
        with open(csv_path, "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for row in rows:
                w.writerow(row)
        saved.append(csv_path)

    return saved


def render_data2_parameters_table(results, *, save_path=None,
                                    output_basename="data2_table_parameters",
                                    **_unused):
    """Same as DATA1 but for DATA2. Reuses the DATA1 implementation."""
    return render_data1_parameters_table(
        results, save_path=save_path, output_basename=output_basename
    )


# ---- Manifest registration -----------------------------------------------

def _register_paper_table_entries():
    existing = {spec.name for spec in CAMPAIGN_MANIFESTS.get("DATA1", ())}
    if "data1.tables.parameters_table" not in existing:
        register_figure("DATA1", FigureSpec(
            name="data1.tables.parameters_table",
            renderer="render_data1_parameters_table",
            requires=(
                RunRequest(campaign="DATA1", run_id="501.1", variant="concpolar",
                           load_cached_fit=True),
                RunRequest(campaign="DATA1", run_id="501.11", variant="concpolar",
                           load_cached_fit=True),
                RunRequest(campaign="DATA1", run_id="511.11", variant="concpolar",
                           load_cached_fit=True),
                RunRequest(campaign="DATA1", run_id="511.12", variant="concpolar",
                           load_cached_fit=True),
            ),
            opts={"output_basename": "data1_table_1_parameters"},
            output_filename=None,    # renderer manages naming
            description="DATA1 paper Table 1: fitted parameters for each (run, variant).",
        ))

    existing2 = {spec.name for spec in CAMPAIGN_MANIFESTS.get("DATA2", ())}
    if "data2.tables.parameters_table" not in existing2:
        # DATA2 fits aren't cached as fit_stru.mat on disk the way DATA1 is.
        # The DATA2 entry is a placeholder that will collect parameters from
        # the cross-verification per-case StageResults once those run.
        register_figure("DATA2", FigureSpec(
            name="data2.tables.parameters_table",
            renderer="render_data2_parameters_table",
            requires=(
                RunRequest(campaign="DATA2", run_id="270611.121"),
                RunRequest(campaign="DATA2", run_id="270711.121"),
                RunRequest(campaign="DATA2", run_id="270511.123"),
            ),
            opts={"output_basename": "data2_table_parameters"},
            output_filename=None,
            description="DATA2 paper parameter table aggregated across representative runs.",
        ))


_register_paper_table_entries()


# ---- Public API --------------------------------------------------

try:
    for _name in ("PAPER_FIGURE_INDEX", "report_paper_coverage",
                  "format_paper_coverage",
                  "render_data1_parameters_table",
                  "render_data2_parameters_table"):
        if _name not in __all__:
            __all__.append(_name)
except NameError:
    pass


# =============================================================================
#
#   NF270 RENDERING IMPROVEMENTS  (nf270 render improvements patch v1)
#
# Addresses five issues observed in the first NF270 fit run:
#   (a) Mute the per-vial vertical alignment artifacts in the plots.
#   (b) Fix the mass-plot scaling - accumulate mass across vials instead of
#       resetting each vial to zero (matches DATA2 figure-6 convention).
#   (c) Detect t_delay (start of permeation) from the pressure ramp instead
#       of trusting data_raw[0]['time'][0]. Lag-mode experiments have a
#       pressure spike at the start that correlates with valve open / first
#       drop entering the vial.
#   (d) Improve the permeate concentration plot - use ICP-OES per-vial
#       measurements (cV_avg) as cyan squares, predicted vial concentration
#       as red triangles, and the predicted retentate/permeate as solid
#       lines overlaid on the dotted measurement line (DATA1 Figure 2 /
#       DATA2 Figure 6 conventions).
#   (e) Use dotted lines for the continuous conductivity-derived
#       measurements and squares/triangles for the discrete ICP-OES
#       measurements, with predicted lines overlaid.
#
# Plus a domain bug fix:
#   - Adds CaCl2 and LaCl3 entries to CONDUCTIVITY_SALT_PARAMS_25C so the
#     loader stops failing on those salts. Values are from CRC handbook
#     ionic conductivities at 25°C; users should refine if their pH or
#     temperature regime is non-standard.
#
# Idempotent: re-importing this block doesn't double-register entries.
# Depends on: nf270 campaign patch v1.
# =============================================================================
# =============================================================================

import numpy as _np
import matplotlib.pyplot as _plt


# ---- Salt conductivity parameters (CRC handbook, 25°C) --------------------
# CaCl2 and LaCl3 are added to the existing CONDUCTIVITY_SALT_PARAMS_25C dict.
# Reference: CRC Handbook of Chemistry and Physics, ionic conductivities at
# infinite dilution at 25°C. Values for hydrated diameters are typical
# Bjerrum / Robinson-Stokes estimates. Users running at non-standard pH or
# temperature should override via model_params= when calling stage_load.

if "CaCl2" not in CONDUCTIVITY_SALT_PARAMS_25C:
    CONDUCTIVITY_SALT_PARAMS_25C["CaCl2"] = {
        "epsilon": 78.4,            # water dielectric at 25 C
        "eta":     0.0089,          # poise
        "a":       1e-8,            # 1 angstrom in cm (Bjerrum length scale)
        "z_1":     2,               # Ca2+ charge
        "z_2":    -1,               # Cl- charge
        "lambda_0_cation": 119.0,   # S cm^2/mol  (Ca2+ at infinite dilution)
        "lambda_0_anion":   76.35,  # Cl-
        "lambda_0":        119.0 + 76.35,
    }

if "LaCl3" not in CONDUCTIVITY_SALT_PARAMS_25C:
    CONDUCTIVITY_SALT_PARAMS_25C["LaCl3"] = {
        "epsilon": 78.4,
        "eta":     0.0089,
        "a":       1e-8,
        "z_1":     3,
        "z_2":    -1,
        "lambda_0_cation": 208.8,   # La3+
        "lambda_0_anion":   76.35,
        "lambda_0":        208.8 + 76.35,
    }


# ---- t_delay detection from pressure --------------------------------------

def _detect_t_delay_from_pressure(time_array, pressure_array, *,
                                    threshold_frac=0.80,
                                    sustain_count=5):
    """Return the time at which applied pressure first ramps up.

    Lag-mode experiments start at atmospheric pressure; once the operator
    closes the cell valve and applies pressure, the cell pressure spikes
    and stays elevated. We define t_delay as the time at which pressure
    first crosses ``threshold_frac × max_sustained_pressure`` and stays
    there for at least ``sustain_count`` consecutive samples.

    Falls back to time_array[0] when pressure data is missing or noisy.
    """
    t = _np.asarray(time_array, dtype=float).reshape(-1)
    p = _np.asarray(pressure_array, dtype=float).reshape(-1)
    if t.size == 0 or p.size == 0:
        return 0.0
    n = min(t.size, p.size)
    t = t[:n]; p = p[:n]
    if n < sustain_count + 5:
        return float(t[0])
    p_target = float(_np.nanmax(p)) * threshold_frac
    if not _np.isfinite(p_target) or p_target <= 0:
        return float(t[0])
    above = p >= p_target
    for i in range(0, n - sustain_count + 1):
        if above[i:i + sustain_count].all():
            return float(t[i])
    crossings = _np.argwhere(above)
    if crossings.size:
        return float(t[int(crossings[0, 0])])
    return float(t[0])


# ---- Improved NF270 renderer ----------------------------------------------

def _ensure_save_dir(save_path):
    """Resolve save_path into a directory Path (creates it if needed)."""
    if save_path is None:
        return None
    p = Path(save_path)
    if p.suffix == "":
        p.mkdir(parents=True, exist_ok=True)
        return p
    p.parent.mkdir(parents=True, exist_ok=True)
    return p.parent


def _detect_perm_start_per_vial(data_raw, *, threshold_g=0.05):
    """For each vial, find the time at which mass first rose above the noise
    floor (i.e., when permeate physically started entering that vial).

    Lag-mode experiments have a holdup period at the very beginning where
    the cell is pressurized but no permeate has yet reached the scintillation
    vial - the balance reads ~constant for several minutes. This function
    detects that empirically by looking for the first sample where mass has
    risen by at least ``threshold_g`` (~one drop, ~0.05 g) above its
    initial reading.

    For vial 1, this typically returns t_perm_start ≈ 2-5 min after the
    pressure ramp. For vials 2+, the system is at steady state so the
    threshold is crossed within the first sample -> returns t[0].

    Returns a list of timestamps (in seconds, same units as the time array)
    one per vial.
    """
    out = []
    for row in data_raw:
        t = _np.asarray(row.get("time", []), dtype=float)
        m = _np.asarray(row.get("mass", []), dtype=float)
        if t.size == 0 or m.size == 0:
            out.append(0.0)
            continue
        m_base = float(m[0])
        rose = (m - m_base) >= threshold_g
        if rose.any():
            idx = int(_np.argmax(rose))
            out.append(float(t[idx]))
        else:
            out.append(float(t[0]))
    return out


def render_nf270_fit_v2(results, *, save_path=None, run_id=None, **_unused):
    """NF270 mass + concentration plots — DATA2 utility.py convention.

    With the loader splitting vial 1 into [extra holdup] + [real vial 1]
    and resetting per-vial mass to start at ~0g, the renderer can use the
    same simple convention DATA1 / DATA2 already use:

      - Single global t_delay = data_raw[0]['time'][0]
      - Plot raw mass values for both data and prediction (no per-vial
        re-anchoring needed — each vial's mass already starts near 0)
      - Skip the prediction trace for vials with index < n_v0 (extras)
      - Data points for the extra vial(s) DO appear (red dots) so the
        viewer sees the holdup-fill period; only the prediction line
        is hidden because the model isn't expected to fit it.

    Concentration plot mirrors DATA1 Figure 2 / DATA2 Figure 6:
      * Conductivity-derived retentate cF_exp  → magenta dotted
      * ICP-OES per-vial cV_avg                → cyan squares
      * Predicted cF                           → green solid
      * Predicted cH                           → red solid
      * Predicted cV (per vial)                → red triangles
    """
    res = _first_result(results)
    res.require("data")
    data_stru = _data_payload(res)
    sim_stru  = _normalize_sim_stru(res.sim_stru or [])

    if run_id and isinstance(res.meta, dict):
        family = NF270_RUN_REGISTRY.get(run_id)
        if family:
            res.meta.setdefault("nf270_workbook", family["workbook"])
            res.meta.setdefault("nf270_sheet",    family["sheet"])
            res.meta.setdefault("nf270_membrane", family["membrane"])

    data_raw = data_stru.get("data_raw", []) if isinstance(data_stru, dict) else []
    if not data_raw:
        return []
    cfg = data_stru.get("data_config", {}) if isinstance(data_stru, dict) else {}
    n_v0 = int(cfg.get("n_v0", 1))   # 1-based: first "real" vial
    t_perm_start_s = cfg.get("t_perm_start_s", None)
    # Single global t_delay = first sample of vial 1 (DATA1 convention).
    t_delay = float(_np.asarray(data_raw[0]["time"], dtype=float).reshape(-1)[0])
    if isinstance(res.meta, dict):
        res.meta["t_delay_s"] = float(t_delay)
        for k in ("vial1_split_status", "t_perm_start_s", "n_holdup_samples",
                  "Lp0", "Lp0_method", "n_v0", "n_extra"):
            if k in cfg:
                res.meta.setdefault(k, cfg[k])

    save_dir = _ensure_save_dir(save_path) or Path.cwd()
    dataset_id = data_stru.get("dataset", "unknown") if isinstance(data_stru, dict) else "unknown"
    saved_paths = []

    # ----- Mass plot: raw values, single t_delay, DATA1 convention -------
    fig_m, ax_m = _plt.subplots(figsize=(5, 4))

    # Plot data (all vials, including extras so the holdup is visible).
    for i, row in enumerate(data_raw):
        t = _np.asarray(row.get("time", []), dtype=float).reshape(-1)
        m = _np.asarray(row.get("mass", []), dtype=float).reshape(-1)
        if t.size == 0 or m.size == 0:
            continue
        ax_m.plot((t - t_delay) / 60.0, m - m[0], "r.", markersize=3, alpha=0.7,
                  label=("Measurements" if i == 0 else None))

    # Plot predictions only for real vials (n_vial >= n_v0, 1-based).
    for i, sim in enumerate(sim_stru):
        # i is 0-based index into sim_stru; sim_stru typically aligns 1:1
        # with data_raw. Skip if this vial is an extra startup vial.
        if (i + 1) < n_v0:
            continue
        t = _np.asarray(sim.get("time", []), dtype=float).reshape(-1)
        mv = _np.asarray(sim.get("mV", []), dtype=float).reshape(-1)
        if t.size == 0 or mv.size == 0:
            continue
        ax_m.plot((t - t_delay) / 60.0, mv - mv[0], "b-", linewidth=2.0, alpha=0.85,
                  label=("Predictions" if i + 1 == n_v0 else None))

    ax_m.set_xlabel("Time [min]", fontsize=14, fontweight="bold")
    ax_m.set_ylabel("Mass [g]", fontsize=14, fontweight="bold")
    ax_m.tick_params(direction="in")
    ax_m.legend(fontsize=10, loc="best")
    ax_m.set_xlim(left=0)
    ax_m.set_ylim(bottom=0)
    fig_m.tight_layout()
    mass_path = save_dir / f"mass-{dataset_id}.png"
    fig_m.savefig(mass_path, dpi=300, bbox_inches="tight")
    _plt.close(fig_m)
    saved_paths.append(mass_path)

    # ----- Concentration plot ----------------------------------------------
    fig_c, ax_c = _plt.subplots(figsize=(5, 4))
    data_raw = data_stru.get("data_raw", [])
    seen_ret_meas = False
    seen_perm_meas = False
    seen_ret_pred = False
    seen_perm_pred = False
    seen_vial_pred = False
    for i, row in enumerate(data_raw):
        time = _np.asarray(row.get("time", []), dtype=float).reshape(-1)
        if time.size == 0:
            continue
        cF_exp = _np.asarray(row.get("cF_exp", []), dtype=float).reshape(-1)
        cV_avg = _np.asarray(row.get("cV_avg", []), dtype=float)
        cV_avg_arr = _np.atleast_1d(cV_avg).astype(float, copy=False)

        if cF_exp.size:
            if cF_exp.size == 1:
                ax_c.plot(
                    _plot_time_minutes(time[-1], t_delay=t_delay),
                    float(cF_exp[0]),
                    "m:",
                    linewidth=1.5,
                    alpha=0.75,
                    label="Retentate (conductivity)" if not seen_ret_meas else None,
                )
            else:
                n = min(len(time), len(cF_exp))
                ax_c.plot(
                    _plot_time_minutes(time[:n], t_delay=t_delay),
                    cF_exp[:n],
                    "m:",
                    linewidth=1.5,
                    alpha=0.75,
                    label="Retentate (conductivity)" if not seen_ret_meas else None,
                )
            seen_ret_meas = True

        if cV_avg_arr.size:
            if cV_avg_arr.ndim == 0 or cV_avg_arr.size == 1:
                ax_c.plot(
                    _plot_time_minutes(time[-1], t_delay=t_delay),
                    float(cV_avg_arr.reshape(-1)[0]),
                    "cs",
                    markersize=8,
                    alpha=0.85,
                    zorder=2,
                    label="Permeate (Measurement)" if not seen_perm_meas else None,
                )
            else:
                n = min(len(time), len(cV_avg_arr))
                ax_c.plot(
                    _plot_time_minutes(time[:n], t_delay=t_delay),
                    cV_avg_arr[:n],
                    "cs",
                    markersize=8,
                    alpha=0.85,
                    zorder=2,
                    label="Permeate (Measurement)" if not seen_perm_meas else None,
                )
            seen_perm_meas = True

        if (i + 1) < n_v0 or i >= len(sim_stru):
            continue
        sim_time = _np.asarray(sim_stru[i].get("time", []), dtype=float).reshape(-1)
        if sim_time.size == 0:
            continue
        cF = _np.asarray(sim_stru[i].get("cF", []), dtype=float).reshape(-1)
        cH = _np.asarray(sim_stru[i].get("cH", []), dtype=float).reshape(-1)
        cV = _np.asarray(sim_stru[i].get("cV", []), dtype=float).reshape(-1)
        pred_line_fx = [
            patheffects.Stroke(linewidth=4.5, foreground="white"),
            patheffects.Normal(),
        ]
        if cF.size:
            ax_c.plot(
                _plot_time_minutes(sim_time, t_delay=t_delay),
                cF,
                color="g",
                linestyle="-",
                linewidth=3.2,
                alpha=1.0,
                zorder=6,
                path_effects=pred_line_fx,
                label="Retentate (prediction)" if not seen_ret_pred else None,
            )
            seen_ret_pred = True
        if cH.size:
            ax_c.plot(
                _plot_time_minutes(sim_time, t_delay=t_delay),
                cH,
                color="r",
                linestyle="-",
                linewidth=3.2,
                alpha=1.0,
                zorder=6,
                path_effects=pred_line_fx,
                label="Permeate (prediction)" if not seen_perm_pred else None,
            )
            seen_perm_pred = True
        if cV.size:
            if cV_avg_arr.ndim == 0 or cV_avg_arr.size == 1:
                ax_c.plot(
                    _plot_time_minutes(sim_time[-1], t_delay=t_delay),
                    float(cV[-1]),
                    "r^",
                    markersize=8,
                    alpha=0.95,
                    zorder=7,
                    label="Vial (prediction)" if not seen_vial_pred else None,
                )
            else:
                valid = ~_np.isnan(cV_avg_arr[: min(len(time), len(cV_avg_arr))])
                if _np.any(valid):
                    time_shifted = time[: min(len(time), len(cV_avg_arr))]
                    interp_cv = interpolate.interp1d(sim_time, cV, fill_value="extrapolate")
                    ax_c.plot(
                        _plot_time_minutes(time_shifted[valid], t_delay=t_delay),
                        interp_cv(time_shifted[valid]),
                        "r^",
                        markersize=8,
                        alpha=0.95,
                        zorder=7,
                        label="Vial (prediction)" if not seen_vial_pred else None,
                    )
            seen_vial_pred = True

    ax_c.set_xlabel("Time [min]", fontsize=14, fontweight="bold")
    ax_c.set_ylabel("Concentration [mM]", fontsize=14, fontweight="bold")
    ax_c.tick_params(direction="in")
    ax_c.legend(fontsize=8, loc="best")
    ax_c.set_xlim(left=0)
    ax_c.set_ylim(bottom=0)
    fig_c.tight_layout()
    conc_path = save_dir / f"concentration-{dataset_id}.png"
    fig_c.savefig(conc_path, dpi=300, bbox_inches="tight")
    _plt.close(fig_c)
    saved_paths.append(conc_path)

    return saved_paths


# ---- Re-register manifest entries to use the new renderer -----------------

def _swap_nf270_renderer_to_v2():
    """Replace render_nf270_fit with render_nf270_fit_v2 in every NF270
    fit FigureSpec. Idempotent: skips entries that already point at v2."""
    from dataclasses import replace

    new_specs = []
    for spec in CAMPAIGN_MANIFESTS.get("NF270", ()):
        renderer = spec.renderer
        if isinstance(renderer, str) and "render_nf270_fit" in renderer and "v2" not in renderer:
            new_specs.append(replace(spec, renderer="render_nf270_fit_v2"))
        else:
            new_specs.append(spec)
    CAMPAIGN_MANIFESTS["NF270"] = tuple(new_specs)


_swap_nf270_renderer_to_v2()


# ---- Public API exposure --------------------------------------------------

try:
    for _name in ("render_nf270_fit_v2", "_detect_t_delay_from_pressure"):
        if _name not in __all__:
            __all__.append(_name)
except NameError:
    pass


# =============================================================================
#
#   MSA CONDUCTIVITY DISPATCH  (nf270 msa conductivity patch v1)
#
# Wires up the existing `msa_transport` function in conductivity_paper.py so
# the diafiltration pipeline can use it for asymmetric electrolytes (CaCl2 2:1,
# LaCl3 3:1) where variant Shedlovsky systematically misestimates the
# concentration. NaCl/KCl continue using variant Shedlovsky by default since
# it's well-validated for symmetric 1:1 electrolytes.
#
# Key constraint: conductivity_paper.py is NOT modified. Only the wrapper in
# the refactored library is extended. msa_transport is called with the exact
# signature it already exposes:
#   msa_transport(valency, diameters, diff_coeff, temp, eta, epsilon,
#                 lambda_0, salt_1_conc, salt_2_conc=None, salt_3_conc=None)
#
# Unit conversions handled at the call boundary (NOT in conductivity_paper.py):
#   - eta:        poise (Shedlovsky)         -> Pa·s        (× 0.1)
#   - lambda_0:   S·cm²/equiv (Shedlovsky)   -> S·m²/mol    (× 1e-4 × |z|)
#   - diameter:   cm "a" parameter (Shedlov.)->  m          (provided per-salt)
#   - concentration: M (Shedlovsky)          -> mM          (× 1000)
#   - return:     mS/cm — same in both paths, so the inversion logic is reused.
#
# Diffusion coefficients are computed via Nernst-Einstein from the existing
# lambda_0 values:    D = R·T / (F² · |z|²) × λ⁰_S_m2_mol
# Hard-sphere diameters are common-literature values (Marcus 1988, Robinson-
# Stokes). Override per-salt via model_params= when calling stage_load.
#
# Default model dispatch (universal MSA — see CONDUCTIVITY_MODEL_NOTES.docx):
#   model="auto"  → always MSA. Reduces correctly to the 1:1 limiting case
#                  for NaCl/KCl (within fit tolerance) AND extends rigorously
#                  to CaCl2 / LaCl3 / multi-salt mixtures. The price is a
#                  modest preprocessing-time cost (~30-60 s/sheet for 1:1
#                  vs Shedlovsky); the methods-consistency win is worth it.
#   model="msa"   → forced MSA (same as auto).
#   model="variant_shedlovsky" / "shedlovsky" → forced classical, only valid
#                  for single 1:1 salts; available for compatibility.
#
# Idempotent. Depends on: nf270 render improvements patch v1 (or any earlier
# patch that registered CaCl2 / LaCl3 entries in CONDUCTIVITY_SALT_PARAMS_25C).
# =============================================================================
# =============================================================================

import contextlib as _msa_contextlib
import io as _msa_io
import numpy as _msa_np

# Physical constants for Nernst-Einstein D = RT / (F²·|z|²) × λ⁰
_MSA_R         = 8.314          # J/(mol·K)
_MSA_T_REF     = 298.15         # K
_MSA_F_FARADAY = 96485.0        # C/mol


# Hard-sphere diameters at 25°C (m). Sources: Marcus (1988), Robinson-Stokes
# Electrolyte Solutions. Conservative literature values; override via
# model_params={"diameter_cation_m": ..., "diameter_anion_m": ...} for
# specific solvent or temperature corrections.
_MSA_ION_DIAMETERS_M = {
    "Na":  3.6e-10,
    "K":   3.0e-10,
    "Li":  3.8e-10,
    "Mg":  4.0e-10,
    "Ca":  4.1e-10,
    "La":  4.6e-10,
    "Cl":  3.6e-10,
    "Br":  3.9e-10,
    "F":   2.7e-10,
    "SO4": 4.4e-10,
}

# Mapping from salt name -> (cation_symbol, anion_symbol).
_MSA_SALT_ION_LOOKUP = {
    "NaCl":   ("Na",  "Cl"),
    "KCl":    ("K",   "Cl"),
    "LiCl":   ("Li",  "Cl"),
    "MgCl2":  ("Mg",  "Cl"),
    "CaCl2":  ("Ca",  "Cl"),
    "LaCl3":  ("La",  "Cl"),
    "Na2SO4": ("Na",  "SO4"),
    "K2SO4":  ("K",   "SO4"),
}


def _msa_diff_coeff_from_lambda0_S_m2_mol(lambda_0_S_m2_mol, z):
    """Nernst-Einstein: D = R·T / (F²·|z|²) × λ⁰ (S·m²/mol)."""
    return (_MSA_R * _MSA_T_REF) / (_MSA_F_FARADAY ** 2 * abs(z) ** 2) * lambda_0_S_m2_mol


def _msa_augment_salt_params():
    """Add MSA-required fields to every salt entry in CONDUCTIVITY_SALT_PARAMS_25C.

    Each entry, after augmentation, has both Shedlovsky inputs (eta in poise,
    lambda_0 in S·cm²/equiv, a in cm) AND MSA inputs (eta_pa_s, lambda_0_*_
    S_m2_mol, diameter_*_m, diff_coeff_*_m2_s) so either model can be used.
    """
    for salt, params in CONDUCTIVITY_SALT_PARAMS_25C.items():
        if "z_1" not in params or "z_2" not in params:
            continue
        z_1 = int(params["z_1"])
        z_2 = int(params["z_2"])

        # Hydrodynamic viscosity: Shedlovsky uses poise; MSA needs Pa·s.
        params.setdefault("eta_pa_s", float(params["eta"]) * 0.1)

        # λ⁰ for each ion: the existing CONDUCTIVITY_SALT_PARAMS_25C stores
        # per-mole-of-ion conductivities in S·cm²/mol (the chemistry-natural
        # convention — what CRC tables list for individual ions). For 1:1
        # salts per-mole = per-equivalent so the existing Shedlovsky code
        # works unchanged. MSA needs S·m²/mol with NO charge multiplier:
        #   1 S·cm²/mol × 1e-4 → 1 S·m²/mol
        params.setdefault(
            "lambda_0_cation_S_m2_mol",
            float(params["lambda_0_cation"]) * 1e-4,
        )
        params.setdefault(
            "lambda_0_anion_S_m2_mol",
            float(params["lambda_0_anion"]) * 1e-4,
        )

        # Diffusion coefficients via Nernst-Einstein.
        params.setdefault(
            "diff_coeff_cation_m2_s",
            _msa_diff_coeff_from_lambda0_S_m2_mol(params["lambda_0_cation_S_m2_mol"], z_1),
        )
        params.setdefault(
            "diff_coeff_anion_m2_s",
            _msa_diff_coeff_from_lambda0_S_m2_mol(params["lambda_0_anion_S_m2_mol"], z_2),
        )

        # Hard-sphere diameters from the lookup; user can override.
        ion_pair = _MSA_SALT_ION_LOOKUP.get(salt)
        if ion_pair:
            cat, an = ion_pair
            params.setdefault("diameter_cation_m", _MSA_ION_DIAMETERS_M.get(cat))
            params.setdefault("diameter_anion_m",  _MSA_ION_DIAMETERS_M.get(an))


_msa_augment_salt_params()


# ---- The MSA branch: builds inputs and inverts numerically -----------------

def _conductivity_to_concentration_msa_branch(
    *,
    cond_signal,
    temp_K,
    salt_name,
    model_params=None,
    output_units="mM",
):
    """MSA equivalent of _conductivity_to_concentration_series.

    Same structure as the Shedlovsky branch — build a forward function
    f(conc_M) → kappa_mS_cm via msa_transport, then numerically invert
    point by point.

    All unit conversions happen here, not in conductivity_paper.py.
    """
    cp = _load_conductivity_paper()
    params = dict(CONDUCTIVITY_SALT_PARAMS_25C.get(str(salt_name), {}))
    if model_params:
        params.update(model_params)
    required = ("diameter_cation_m", "diameter_anion_m",
                "diff_coeff_cation_m2_s", "diff_coeff_anion_m2_s",
                "eta_pa_s", "epsilon",
                "lambda_0_cation_S_m2_mol", "lambda_0_anion_S_m2_mol",
                "z_1", "z_2")
    missing = [k for k in required if k not in params]
    if missing:
        raise ValueError(
            f"MSA conductivity model for salt {salt_name!r} is missing fields: "
            f"{missing}. Pass them via model_params={{...}} or extend "
            f"CONDUCTIVITY_SALT_PARAMS_25C[{salt_name!r}]."
        )

    cond_signal = _msa_np.asarray(cond_signal, dtype=float).reshape(-1)
    conc_out = _msa_np.full_like(cond_signal, _msa_np.nan, dtype=float)

    valency      = [int(params["z_1"]), int(params["z_2"])]
    diameters    = [float(params["diameter_cation_m"]),
                    float(params["diameter_anion_m"])]
    diff_coeff   = [float(params["diff_coeff_cation_m2_s"]),
                    float(params["diff_coeff_anion_m2_s"])]
    eta_pa_s     = float(params["eta_pa_s"])
    epsilon      = float(params["epsilon"])
    lambda_0     = [float(params["lambda_0_cation_S_m2_mol"]),
                    float(params["lambda_0_anion_S_m2_mol"])]

    _stdout_sink = _msa_io.StringIO()

    def fwd(conc_M):
        # M -> mM (msa_transport convention) and silence its print() chatter.
        conc_mM = conc_M * 1000.0
        with _msa_contextlib.redirect_stdout(_stdout_sink):
            kappa_mS_cm = cp.msa_transport(
                valency,
                diameters,
                diff_coeff,
                temp_K,
                eta_pa_s,
                epsilon,
                lambda_0,
                [conc_mM],
            )[0]
        return float(kappa_mS_cm)

    for i, y in enumerate(cond_signal):
        if _msa_np.isnan(y):
            continue
        conc_out[i] = _invert_monotone_1d(
            fwd_func=fwd,
            y_target=float(y / 1000.0),   # µS/cm -> mS/cm
            x_lo=0.0,
            x_hi=6.0,
        )

    if str(output_units).lower() in {"mm", "mmol/l", "mmolar"}:
        conc_out *= 1000.0
    return conc_out


# ---- Wrap _conductivity_to_concentration_series ---------------------------

# Cache the original so subsequent re-imports don't recurse.
if "_msa_orig_conductivity_to_concentration_series" not in globals():
    _msa_orig_conductivity_to_concentration_series = _conductivity_to_concentration_series


def _conductivity_to_concentration_series(*, cond_signal, temp_K, salt_name,
                                            model="auto", model_params=None,
                                            output_units="mM"):
    """Patched: dispatches model='auto' and 'msa' before delegating to the
    original Shedlovsky-only implementation."""
    model_lc = str(model).lower()

    # "auto" resolves to MSA universally — same theory for every salt
    # (single 1:1, single asymmetric, and multi-salt mixtures). MSA reduces
    # correctly to the 1:1 limiting case and extends rigorously to the
    # asymmetric / multi-salt regimes where variant Shedlovsky breaks down.
    # See CONDUCTIVITY_MODEL_NOTES.docx for the methods rationale.
    if model_lc == "auto":
        model_lc = "msa"

    if model_lc == "msa":
        return _conductivity_to_concentration_msa_branch(
            cond_signal=cond_signal,
            temp_K=temp_K,
            salt_name=salt_name,
            model_params=model_params,
            output_units=output_units,
        )

    # Otherwise: classical Shedlovsky / variant Shedlovsky (original path).
    return _msa_orig_conductivity_to_concentration_series(
        cond_signal=cond_signal,
        temp_K=temp_K,
        salt_name=salt_name,
        model=model_lc,
        model_params=model_params,
        output_units=output_units,
    )


# ---- Default model = 'auto' for the loader's normalize step ----------------

if "_msa_orig_normalize_conductivity_measurements" not in globals():
    _msa_orig_normalize_conductivity_measurements = _normalize_conductivity_measurements


def _normalize_conductivity_measurements(data_stru, *, output_units="mM",
                                           model="auto", model_params=None):
    """Patched: pass through to the original with model='auto' so asymmetric
    salts auto-dispatch to MSA."""
    return _msa_orig_normalize_conductivity_measurements(
        data_stru,
        output_units=output_units,
        model=model,
        model_params=model_params,
    )


# ---- Public API exposure ---------------------------------------------------

try:
    for _name in ("_conductivity_to_concentration_msa_branch",
                  "_msa_augment_salt_params",
                  "_MSA_ION_DIAMETERS_M",
                  "_MSA_SALT_ION_LOOKUP"):
        if _name not in __all__:
            __all__.append(_name)
except NameError:
    pass


# =============================================================================
# =============================================================================
#
#   NF270 EXPERIMENTAL CAMPAIGN  (pipeline patch)
#
# Adds the current NF270 single-salt / mixed-salt workbook campaign as a
# first-class campaign in the unified pipeline. The campaign is exposed
# through the same dispatcher and coverage-report machinery as DATA1 and
# DATA2, while keeping the Excel-specific loading path inside the pipeline.
#
# This section adds:
#   - NF270_RUN_REGISTRY: the 26 green-flagged workbook sheets.
#   - NF270_RUN_QUALITY: campaign quality markers for filtering.
#   - get_experimental_run() support for campaign="NF270".
#   - render_nf270_fit: a thin wrapper around the existing DATA3 renderer.
#   - NF270 campaign manifest entries and paper coverage index entries.
#   - NF270 parameter-table aggregation.
# =============================================================================
# =============================================================================

NF270_DEFAULT_ROOT = Path(__file__).resolve().parents[1] / "UnifiedFramework" / "ExperimentalDataFiles"

NF270_RUN_REGISTRY = {
    # --- NF270_MC2.xlsx (15 green) ----------------------------------------
    "MC2.05.07.24_NaCl":  {"workbook": "NF270_MC2.xlsx", "sheet": "05.07.24_NaCl",  "membrane": "NF270_MC2"},
    "MC2.05.07.24_CaCl2": {"workbook": "NF270_MC2.xlsx", "sheet": "05.07.24_CaCl2", "membrane": "NF270_MC2"},
    "MC2.05.21.24_LaCl3": {"workbook": "NF270_MC2.xlsx", "sheet": "05.21.24_LaCl3", "membrane": "NF270_MC2"},
    "MC2.05.08.24_E1":    {"workbook": "NF270_MC2.xlsx", "sheet": "05.08.24_E1",    "membrane": "NF270_MC2"},
    "MC2.05.08.24_E2":    {"workbook": "NF270_MC2.xlsx", "sheet": "05.08.24_E2",    "membrane": "NF270_MC2"},
    "MC2.05.08.24_E3":    {"workbook": "NF270_MC2.xlsx", "sheet": "05.08.24_E3",    "membrane": "NF270_MC2"},
    "MC2.05.08.24_E4":    {"workbook": "NF270_MC2.xlsx", "sheet": "05.08.24_E4",    "membrane": "NF270_MC2"},
    "MC2.05.09.24_E5":    {"workbook": "NF270_MC2.xlsx", "sheet": "05.09.24_E5",    "membrane": "NF270_MC2"},
    "MC2.05.21.24_E6":    {"workbook": "NF270_MC2.xlsx", "sheet": "05.21.24_E6",    "membrane": "NF270_MC2"},
    "MC2.05.21.24_E7":    {"workbook": "NF270_MC2.xlsx", "sheet": "05.21.24_E7",    "membrane": "NF270_MC2"},
    "MC2.05.22.24_E8":    {"workbook": "NF270_MC2.xlsx", "sheet": "05.22.24_E8",    "membrane": "NF270_MC2"},
    "MC2.05.22.24_E10":   {"workbook": "NF270_MC2.xlsx", "sheet": "05.22.24_E10",   "membrane": "NF270_MC2"},
    "MC2.05.28.24_E11":   {"workbook": "NF270_MC2.xlsx", "sheet": "05.28.24_E11",   "membrane": "NF270_MC2"},
    "MC2.05.28.24_E12":   {"workbook": "NF270_MC2.xlsx", "sheet": "05.28.24_E12",   "membrane": "NF270_MC2"},
    "MC2.05.30.24_E15":   {"workbook": "NF270_MC2.xlsx", "sheet": "05.30.24_E15",   "membrane": "NF270_MC2"},

    # --- NF270_MC3.xlsx (3 green) -----------------------------------------
    "MC3.07.22.24_SNaCl":   {"workbook": "NF270_MC3.xlsx", "sheet": "07.22.24_SNaCl",   "membrane": "NF270_MC3"},
    "MC3.07.11.24_SCaCl2":   {"workbook": "NF270_MC3.xlsx", "sheet": "07.11.24_SCaCl2",  "membrane": "NF270_MC3"},
    "MC3.07.12.24_S2CaCl2":  {"workbook": "NF270_MC3.xlsx", "sheet": "07.12.24_S2CaCl2", "membrane": "NF270_MC3"},

    # --- NF270_MC4.xlsx (5 green) -----------------------------------------
    "MC4.07.11.24_SNaCl":   {"workbook": "NF270_MC4.xlsx", "sheet": "07.11.24_SNaCl",   "membrane": "NF270_MC4"},
    "MC4.07.11.24_SLaCl3":  {"workbook": "NF270_MC4.xlsx", "sheet": "07.11.24_SLaCl3",  "membrane": "NF270_MC4"},
    "MC4.07.12.24_E13":     {"workbook": "NF270_MC4.xlsx", "sheet": "07.12.24_E13",     "membrane": "NF270_MC4"},
    "MC4.07.16.24_E16":     {"workbook": "NF270_MC4.xlsx", "sheet": "07.16.24_E16",     "membrane": "NF270_MC4"},
    "MC4.07.16.24_E17":     {"workbook": "NF270_MC4.xlsx", "sheet": "07.16.24_E17",     "membrane": "NF270_MC4"},

    # --- NF270_MC5.xlsx (3 green) -----------------------------------------
    "MC5.07.23.24_NaCl":    {"workbook": "NF270_MC5.xlsx", "sheet": "07.23.24_NaCl",    "membrane": "NF270_MC5"},
    "MC5.07.23.24_SNaCl":   {"workbook": "NF270_MC5.xlsx", "sheet": "07.23.24_SNaCl",   "membrane": "NF270_MC5"},
    "MC5.07.23.24_S2NaCl":  {"workbook": "NF270_MC5.xlsx", "sheet": "07.23.24_S2NaCl",  "membrane": "NF270_MC5"},
}

NF270_RUN_QUALITY = {run_id: "green" for run_id in NF270_RUN_REGISTRY}

NF270_SINGLE_SALT_RUNS = (
    "MC2.05.07.24_CaCl2",
    "MC2.05.07.24_NaCl",
    "MC2.05.21.24_LaCl3",
    "MC3.07.11.24_SCaCl2",
    "MC3.07.12.24_S2CaCl2",
    "MC3.07.22.24_SNaCl",
    "MC4.07.11.24_SLaCl3",
    "MC4.07.11.24_SNaCl",
    "MC5.07.23.24_NaCl",
    "MC5.07.23.24_S2NaCl",
    "MC5.07.23.24_SNaCl",
)

NF270_FAST_SUBSET = (
    "MC2.05.07.24_NaCl",
    "MC3.07.22.24_SNaCl",
    "MC5.07.23.24_NaCl",
)


_original_get_experimental_run_nf270 = get_experimental_run


def get_experimental_run(campaign, run_id, variant="base", data_root=None):
    """Build one file-backed experimental run instance."""
    campaign = str(campaign or "DATA1").upper()
    run_id = str(run_id)
    variant = str(variant or "base").lower()

    if campaign == "NF270":
        family = NF270_RUN_REGISTRY.get(run_id)
        if family is None:
            raise KeyError(f"Unknown NF270 run '{run_id}'. Known runs: {list(NF270_RUN_REGISTRY)}")
        root = Path(data_root) if data_root is not None else NF270_DEFAULT_ROOT
        return ExperimentalRun(
            run_id=run_id,
            root=root,
            data_file=family["workbook"],
            variant=variant,
            subdir="",
            fit_file="",
            contour_files=(),
            extra_files=(),
        )

    return _original_get_experimental_run_nf270(
        campaign, run_id, variant=variant, data_root=data_root
    )


def render_nf270_fit(results, *, save_path=None, run_id=None, **opts):
    """One NF270 fit: mass-vs-time + concentration-vs-time plots."""
    res = _first_result(results)
    res.require("data")
    if run_id and isinstance(res.meta, dict):
        family = NF270_RUN_REGISTRY.get(str(run_id))
        if family:
            res.meta.setdefault("nf270_workbook", family["workbook"])
            res.meta.setdefault("nf270_sheet", family["sheet"])
            res.meta.setdefault("nf270_membrane", family["membrane"])
    save_dir = Path(save_path) if save_path is not None else None
    paths = run_data3_time_series_plots(res.to_dict(), save_dir=save_dir, show=False)
    return [Path(p) for p in (paths or [])]


def _build_nf270_figures():
    """Return a tuple of FigureSpec - one per green-flagged NF270 run."""
    specs = []
    for run_id, family in NF270_RUN_REGISTRY.items():
        specs.append(
            FigureSpec(
                name=f"nf270.{run_id}",
                renderer="render_nf270_fit",
                requires=(
                    RunRequest(
                        campaign="NF270",
                        run_id=run_id,
                        workflow_family="DATA3",
                        mode="DATA",
                        B_form="single",
                        load_cached_fit=False,
                        use_parmest=False,
                        multistart=True,
                        multistart_iterations=10,
                        uncertainty_method="fim",
                        selector=family["sheet"],
                    ),
                ),
                opts={"run_id": run_id, "membrane": family["membrane"]},
                output_filename=None,
                description=(
                    f"NF270 fit for {family['workbook']} :: {family['sheet']} "
                    f"(membrane sample {family['membrane']})."
                ),
            )
        )
    return tuple(specs)


def _register_nf270_campaign():
    """Idempotently register the NF270 campaign + manifest entries."""
    if "NF270" not in CAMPAIGN_REGISTRY:
        CAMPAIGN_REGISTRY["NF270"] = {
            "root": NF270_DEFAULT_ROOT,
            "paper": "NF270",
            "description": "NF270 mixed-salt and single-salt diafiltration experimental campaign.",
        }

    existing = {spec.name for spec in CAMPAIGN_MANIFESTS.get("NF270", ())}
    new_specs = tuple(spec for spec in _build_nf270_figures() if spec.name not in existing)
    if new_specs:
        CAMPAIGN_MANIFESTS["NF270"] = CAMPAIGN_MANIFESTS.get("NF270", ()) + new_specs


def _build_nf270_paper_index():
    """Build NF270 entry for PAPER_FIGURE_INDEX."""
    main_figures = {}
    for run_id, family in NF270_RUN_REGISTRY.items():
        workbook_stem = Path(family["workbook"]).stem
        sheet = family["sheet"]
        manifest_entry = (f"nf270.{run_id}",)
        main_figures[f"mass-{workbook_stem}_{sheet}.png"] = manifest_entry
        main_figures[f"concentration-{workbook_stem}_{sheet}.png"] = manifest_entry
    tables = {
        "nf270_fit_summary.json": ("nf270.tables.parameters_table",),
        "nf270_fit_summary.csv": ("nf270.tables.parameters_table",),
    }
    return {"main_figures": main_figures, "si_figures": {}, "tables": tables}


def _register_nf270_paper_index():
    if "NF270" not in PAPER_FIGURE_INDEX:
        PAPER_FIGURE_INDEX["NF270"] = _build_nf270_paper_index()


def render_nf270_parameters_table(results, *, save_path=None,
                                  output_basename="nf270_fit_summary", **_unused):
    """Aggregate NF270 fitted parameters into a CSV + JSON table bundle."""
    return render_data1_parameters_table(
        results, save_path=save_path, output_basename=output_basename
    )


def _register_nf270_table_entry():
    existing = {spec.name for spec in CAMPAIGN_MANIFESTS.get("NF270", ())}
    if "nf270.tables.parameters_table" in existing:
        return
    register_figure("NF270", FigureSpec(
        name="nf270.tables.parameters_table",
        renderer="render_nf270_parameters_table",
        requires=tuple(
            RunRequest(
                campaign="NF270",
                run_id=run_id,
                workflow_family="DATA3",
                mode="DATA",
                B_form="single",
                load_cached_fit=False,
                use_parmest=False,
                multistart=True,
                multistart_iterations=10,
                uncertainty_method="fim",
                selector=NF270_RUN_REGISTRY[run_id]["sheet"],
            )
            for run_id in NF270_RUN_REGISTRY
        ),
        opts={"output_basename": "nf270_fit_summary"},
        output_filename=None,
        description="NF270 fitted-parameter summary table across the full green campaign.",
    ))


_register_nf270_campaign()
_register_nf270_paper_index()
_register_nf270_table_entry()


# =============================================================================
# DATA1 MATLAB → Python ports
# -----------------------------------------------------------------------------
# Literal Python equivalents of the four MATLAB contour / sensitivity scripts
# in legacy/data1_matlab/functions/.  Names mirror the .m filenames with a
# `_py` suffix.  These wrap the existing DATA3 NF270 forward-sim kernel
# (_nf270_contour_objectives_at_theta) so they produce the SAME log10-SSR
# values per channel that the MATLAB pipeline produces — but using the
# refactored Pyomo + Sundials forward simulator instead of MATLAB's sim_model.
#
# Why ports (rather than reusing _nf270_contour_grid_dataframe directly):
# one-to-one MATLAB↔Python mapping for cross-validation and for anyone
# reading the .m source alongside the .py source.  Output CSV / PNG filenames
# match the MATLAB conventions so legacy diff tools work unchanged.
#
# All four ports default to workflow_family="DATA3"; pass "DATA1" / "DATA2"
# to call them from those campaigns.
# =============================================================================

from collections import namedtuple

_IndObjectivesPy   = namedtuple("_IndObjectivesPy",   ["m", "cp", "cr"])
_IndObjectives5Py  = namedtuple("_IndObjectives5Py",  ["m", "cp", "cr", "cp_cond", "cr_cond"])


def _forward_shedlovsky_mM_to_uS_per_cm(c_mM, salt_name, temp_C=25.0):
    """Forward-Shedlovsky: convert salt concentration (mM) → conductivity (µS/cm).

    Wraps conductivity_paper.variant_shedlovsky, handling the unit conversions:
      - input: c in mM and T in °C
      - shedlovsky wants: c in M, T in K
      - output: mS/cm → multiply by 1000 → µS/cm

    Parameters
    ----------
    c_mM : float or array of floats
        Salt concentration(s) in mM.
    salt_name : str
        Key into CONDUCTIVITY_SALT_PARAMS_25C (e.g. "NaCl", "CaCl2", "LaCl3").
    temp_C : float
        Temperature in °C.  Default 25 — matches the temperature-compensation
        convention used by the data-loader (everything in data_raw is
        EC25-compensated).

    Returns
    -------
    np.ndarray
        Predicted specific conductivity in µS/cm.
    """
    import numpy as np
    import conductivity_paper as _cp
    sp = CONDUCTIVITY_SALT_PARAMS_25C.get(salt_name)
    if sp is None:
        return np.full_like(np.asarray(c_mM, dtype=float), np.nan)
    shed_keys = ("epsilon", "eta", "lambda_0", "a",
                 "z_1", "z_2", "lambda_0_cation", "lambda_0_anion")
    kw = {k: sp[k] for k in shed_keys}
    c_arr = np.atleast_1d(np.asarray(c_mM, dtype=float))
    c_M = c_arr / 1000.0
    T_K = float(temp_C) + 273.15
    sigma_mS_cm = np.asarray(_cp.variant_shedlovsky(c_M, T_K, **kw), dtype=float)
    return sigma_mS_cm * 1000.0   # → µS/cm


def calc_ind_objectives_py(
    theta,
    data_stru,
    *,
    mode="Lag",
    B_form="single",
    workflow_family="DATA3",
    nfe=80,
):
    """Python port of legacy/data1_matlab/functions/calc_ind_objectives.m.

    Forward-simulates one (data_stru, theta) and returns three per-channel
    weighted SSR scalars: (m, cp, cr) — mass, permeate-concentration,
    retentate-concentration.

    Implementation: thin wrapper over _nf270_contour_objectives_at_theta,
    which itself calls solve_model(sim_opt=True) → forward DAE solve, then
    reads the per-channel objective Expressions (obj_m / obj_cv / obj_cr).
    The weight convention matches the MATLAB weight_config "default" path
    (best_weight.statistical = .equal_statistical = .L1 = .L2 = .Linf = false).
    No scaling applied (scale_opt = 0 in the MATLAB; the scale_opt = 1
    Utopia/Nadir rescaling branch is intentionally omitted — it is not used
    by the DATA1 contour workflow).

    Parameters
    ----------
    theta : dict
        Must contain Lp, B, sigma, S0, S (load from a warm-start fit, or
        run a centering fit once).  Missing S0 raises KeyError inside
        model_construct_inter.
    data_stru : dict
        Loaded sheet (from loadxlsx / loadmat).

    Returns
    -------
    namedtuple (m, cp, cr) : (float, float, float)
        Per-channel weighted SSR.  Failed forward sims return (nan, nan, nan).

    Notes
    -----
    MATLAB source: legacy/data1_matlab/functions/calc_ind_objectives.m
                   Lines 56-85 (weighted-residual accumulation),
                   Lines 87-91 (SSR sum into obj_ind.unscaled).
    """
    try:
        obj_m, obj_cv, obj_cr = _nf270_contour_objectives_at_theta(
            data_stru, theta,
            mode=mode, B_form=B_form,
            workflow_family=workflow_family, nfe=nfe,
        )
        if obj_m is None or obj_cv is None or obj_cr is None:
            return _IndObjectivesPy(float("nan"), float("nan"), float("nan"))
        return _IndObjectivesPy(float(obj_m), float(obj_cv), float(obj_cr))
    except Exception:
        return _IndObjectivesPy(float("nan"), float("nan"), float("nan"))


def calc_ind_objectives_5channel_py(
    theta,
    data_stru,
    *,
    mode="Lag",
    B_form="single",
    workflow_family="DATA3",
    nfe=80,
    rel_scale_cond=0.03,
):
    """5-channel extension of calc_ind_objectives_py — adds conductivity-domain SSRs.

    Returns (m, cp, cr, cp_cond, cr_cond):
      - m       : mass SSR              (unchanged from 3-channel)
      - cp      : permeate-conc SSR     (unchanged — vial-level ICP-OES)
      - cr      : retentate-conc SSR    (unchanged — cF time-series in mM)
      - cp_cond : permeate-conductivity SSR  (NEW — µS/cm domain, time-series)
      - cr_cond : retentate-conductivity SSR (NEW — µS/cm domain, time-series)

    For the two new channels:
      - Forward-Shedlovsky the predicted concentrations cV(t) and cF(t) to
        µS/cm via _forward_shedlovsky_mM_to_uS_per_cm.
      - For permeate: compare against ``cV_perm_cond`` (raw 5-Hz probe in µS/cm,
        EC25-compensated) at probe sample times — note we use m.cH-derived
        predictions, matching the existing m.obj_cv_perm convention (the
        in-line probe measures INSTANTANEOUS outlet, not vial-integrated avg).
      - For retentate: compare against ``cF_exp_conductivity`` (raw inline-probe
        µS/cm) at probe sample times — bulk cF is what the model fits.
      - Weighting: relative scaling ``rel_scale_cond * σ_meas`` (default 3 %),
        matching the existing obj_cv form.

    DATA3-only — returns (m, cp, cr, nan, nan) for non-DATA3 workflows or
    sheets that don't carry probe time-series data.

    Parameters
    ----------
    theta : dict
        Must contain Lp, B, sigma, S0, S.
    rel_scale_cond : float
        Relative scaling for the two new channels.  Default 3 %.

    Returns
    -------
    namedtuple (m, cp, cr, cp_cond, cr_cond)
        Failed forward sims return (nan,)*5.

    Notes
    -----
    Posing of the new channels is symmetric to the existing concentration
    channels in form (weighted SSR with relative scaling) but in the µS/cm
    domain, giving genuinely new information beyond the mM-domain residuals
    already in cp and cr.  Mass and the two concentration channels are
    unchanged from calc_ind_objectives_py.
    """
    import numpy as np

    nan5 = _IndObjectives5Py(float("nan"), float("nan"), float("nan"),
                              float("nan"), float("nan"))
    if str(workflow_family).upper() != "DATA3":
        # DATA1/DATA2 — defer to the 3-channel kernel, leave the
        # conductivity channels as NaN (probe data isn't present in .mat files).
        ind3 = calc_ind_objectives_py(theta, data_stru,
                                      mode=mode, B_form=B_form,
                                      workflow_family=workflow_family, nfe=nfe)
        return _IndObjectives5Py(ind3.m, ind3.cp, ind3.cr,
                                  float("nan"), float("nan"))

    try:
        fit_stru, sim_stru, _ = solve_model(
            data_stru, mode, theta=theta, sim_opt=True,
            B_form=B_form, workflow_family=workflow_family, nfe=nfe, LOUD=False,
        )
        if not isinstance(fit_stru, dict):
            return nan5

        obj_m  = float(fit_stru.get("obj_m",  float("nan")))
        obj_cp = float(fit_stru.get("obj_cv", float("nan")))
        obj_cr = float(fit_stru.get("obj_cr", float("nan")))

        # ---- Conductivity channels (forward Shedlovsky on predictions) ----
        salt = str(data_stru.get("data_config", {}).get("namec", "")).strip()
        n_vials = int(data_stru.get("data_config", {}).get("n", len(sim_stru)))
        data_raw = data_stru.get("data_raw", [])

        sum_perm_cond = 0.0
        cnt_perm_cond = 0
        sum_ret_cond  = 0.0
        cnt_ret_cond  = 0

        for i in range(min(n_vials, len(sim_stru), len(data_raw))):
            vial_sim = sim_stru[i] if isinstance(sim_stru, (list, tuple, dict)) else None
            vial_raw = data_raw[i]
            if vial_sim is None:
                continue
            t_sim = np.asarray(vial_sim.get("time"), dtype=float)
            cV_pred = np.asarray(vial_sim.get("cH"), dtype=float)  # outlet, matches probe
            cF_pred = np.asarray(vial_sim.get("cF"), dtype=float)
            if t_sim.size == 0:
                continue

            # ---- Permeate conductivity channel ----
            cV_perm_meas = np.asarray(vial_raw.get("cV_perm_cond", []), dtype=float)
            t_probe = np.asarray(vial_raw.get("time", []), dtype=float)
            if cV_perm_meas.size > 0 and t_probe.size == cV_perm_meas.size:
                cV_pred_at_probe = np.interp(t_probe, t_sim, cV_pred)
                sigma_pred_perm = _forward_shedlovsky_mM_to_uS_per_cm(cV_pred_at_probe, salt)
                for k in range(len(t_probe)):
                    sm = cV_perm_meas[k]
                    sp = sigma_pred_perm[k]
                    if not (np.isfinite(sm) and np.isfinite(sp) and sm > 0):
                        continue
                    scale = rel_scale_cond * sm
                    if scale <= 0:
                        continue
                    sum_perm_cond += ((sp - sm) / scale) ** 2
                    cnt_perm_cond += 1

            # ---- Retentate conductivity channel ----
            cF_cond_meas = np.asarray(vial_raw.get("cF_exp_conductivity", []),
                                      dtype=float)
            if cF_cond_meas.size > 0 and t_probe.size == cF_cond_meas.size:
                cF_pred_at_probe = np.interp(t_probe, t_sim, cF_pred)
                sigma_pred_ret = _forward_shedlovsky_mM_to_uS_per_cm(cF_pred_at_probe, salt)
                for k in range(len(t_probe)):
                    sm = cF_cond_meas[k]
                    sp = sigma_pred_ret[k]
                    if not (np.isfinite(sm) and np.isfinite(sp) and sm > 0):
                        continue
                    scale = rel_scale_cond * sm
                    if scale <= 0:
                        continue
                    sum_ret_cond += ((sp - sm) / scale) ** 2
                    cnt_ret_cond += 1

        obj_cp_cond = sum_perm_cond / cnt_perm_cond if cnt_perm_cond > 0 else float("nan")
        obj_cr_cond = sum_ret_cond  / cnt_ret_cond  if cnt_ret_cond  > 0 else float("nan")

        return _IndObjectives5Py(obj_m, obj_cp, obj_cr, obj_cp_cond, obj_cr_cond)
    except Exception:
        return nan5


def calc_contour_2d_py(
    data_stru,
    theta,
    x_var,
    y_var,
    *,
    grid_density=50,
    save_dir,
    x_bounds=None,
    y_bounds=None,
    mode="Lag",
    B_form="single",
    workflow_family="DATA3",
    nfe=80,
    write_csv=True,
    include_conductivity=False,
):
    """Python port of legacy/data1_matlab/functions/calc_contour_2d.m.

    2D parameter sweep on a linear meshgrid; for every cell, copy theta,
    override theta[x_var] and theta[y_var], and call calc_ind_objectives_py.
    Stores log10 of each per-channel SSR.

    Lp special rule (MATLAB lines 42-50): if x_var == "Lp", x_bounds defaults
    to [0.1 * theta["Lp"], 2 * theta["Lp"]] (and same for y_var == "Lp").
    For other parameters: sigma defaults to (0, 1); B defaults to the
    salt-specific NF270_B_BOUNDS_PER_SALT triple in DATA3, else (1e-6, 50).

    Parameters
    ----------
    data_stru, theta : dict
        Loaded sheet + parameter dict (must contain Lp, B, sigma, S0, S).
    x_var, y_var : {"Lp", "B", "sigma"}
        The two parameters to sweep.
    grid_density : int
        Cells per axis.  Default 50 matches MATLAB.
    save_dir : Path
        Output directory; CSV is written as
        ``contourdata-x_{x_var}-y_{y_var}.csv`` (MATLAB filename convention).
    x_bounds, y_bounds : (lo, hi), optional
        Override default bounds.

    Returns
    -------
    pandas.DataFrame
        Long-form with columns [x_var, y_var, Obj_mass, Obj_concentration,
        Obj_retentate_concentration].  Same 5-column schema as the MATLAB
        contourdata CSV (line 105-107 in calc_contour_2d.m).

    Notes
    -----
    MATLAB source: legacy/data1_matlab/functions/calc_contour_2d.m
                   Lines 42-50 (Lp special rule), 52-53 (linspace meshgrid),
                   Lines 65-83 (inner loop), 100-107 (CSV write).
    """
    import numpy as np
    import pandas as pd

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    def _default_bounds(var):
        if var == "Lp":
            lp = float(theta.get("Lp", 1.0))
            return (0.1 * lp, 2.0 * lp)
        if var == "sigma":
            return (0.0, 1.0)
        if var == "B":
            if workflow_family == "DATA3":
                salt = str(data_stru.get("data_config", {}).get("namec", "")).strip()
                return NF270_B_BOUNDS_PER_SALT.get(salt, NF270_B_BOUNDS_DEFAULT)
            return (1e-6, 50.0)
        raise ValueError(f"unknown var {var!r}")

    x_lb, x_ub = x_bounds or _default_bounds(x_var)
    y_lb, y_ub = y_bounds or _default_bounds(y_var)

    x_vals = np.linspace(x_lb, x_ub, grid_density)
    y_vals = np.linspace(y_lb, y_ub, grid_density)
    XX, YY = np.meshgrid(x_vals, y_vals)
    xx = XX.ravel()
    yy = YY.ravel()
    total = len(xx)

    def _log10_or_nan(v):
        try:
            if v is None or not np.isfinite(v) or v <= 0:
                return float("nan")
            return float(np.log10(v))
        except (TypeError, ValueError):
            return float("nan")

    rows = []
    for k in range(total):
        b = dict(theta)
        b[x_var] = float(xx[k])
        b[y_var] = float(yy[k])
        if include_conductivity:
            ind = calc_ind_objectives_5channel_py(
                b, data_stru,
                mode=mode, B_form=B_form,
                workflow_family=workflow_family, nfe=nfe,
            )
            row = {
                x_var: float(xx[k]),
                y_var: float(yy[k]),
                "Obj_mass":                    _log10_or_nan(ind.m),
                "Obj_concentration":           _log10_or_nan(ind.cp),
                "Obj_retentate_concentration": _log10_or_nan(ind.cr),
                "Obj_permeate_conductivity":   _log10_or_nan(ind.cp_cond),
                "Obj_retentate_conductivity":  _log10_or_nan(ind.cr_cond),
            }
        else:
            ind = calc_ind_objectives_py(
                b, data_stru,
                mode=mode, B_form=B_form,
                workflow_family=workflow_family, nfe=nfe,
            )
            row = {
                x_var: float(xx[k]),
                y_var: float(yy[k]),
                "Obj_mass":                    _log10_or_nan(ind.m),
                "Obj_concentration":           _log10_or_nan(ind.cp),
                "Obj_retentate_concentration": _log10_or_nan(ind.cr),
            }
        rows.append(row)

    df = pd.DataFrame(rows)
    if write_csv:
        csv_path = save_dir / f"contourdata-x_{x_var}-y_{y_var}.csv"
        df.to_csv(csv_path, index=False)
    return df


def plot_contour_py(
    df,
    x_var,
    y_var,
    *,
    save_dir,
    title="",
    n_levels=10,
):
    """Python port of legacy/data1_matlab/functions/plot_contour.m.

    Render a 14×7-inch figure with three side-by-side log10(SSR) panels —
    Mass / Permeate concentration / Retentate concentration — each showing
    labeled iso-objective contour lines and a red triangle at the argmin.

    Matches the MATLAB:
      - 1×3 subplot layout, default contour levels, no fill, no colorbar.
      - ``ShowText='on'``         → ``ax.clabel(inline=True)``
      - ``LineWidth=2``           → ``linewidths=2``
      - Argmin marker: '^' size 8, edge red, face [1, 0.6, 0.6]
      - Axis labels via loadlabel(): ``L_p [L/m²/h/bar]``, ``B [µm/s]``,
        ``σ [dimensionless]``
      - Saved as: ``objcontour-x_{x_var}-y_{y_var}.png``

    Parameters
    ----------
    df : pandas.DataFrame
        Output of calc_contour_2d_py (long-form 5-col).
    x_var, y_var : str
        Axis names — must match column names in df.

    Returns
    -------
    pathlib.Path
        Path to the saved PNG.

    Notes
    -----
    MATLAB source: legacy/data1_matlab/functions/plot_contour.m
                   Lines 49-87 (three subplots, contour + clabel + argmin),
                   Line 95 (saveas/png export).
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    labels = {
        "Lp":    r"$L_p$  [L / m² / h / bar]",
        "B":     r"$B$  [µm / s]",
        "sigma": r"$\sigma$  [dimensionless]",
    }
    # Build the panel list dynamically — 3 panels for legacy DATA1/DATA2-style
    # contours, 5 panels when conductivity channels are present in df.
    base_specs = [
        ("Obj_mass",                    r"$\log_{10}$  Mass Residual Squared [g²]"),
        ("Obj_concentration",           r"$\log_{10}$  Permeate Concentration Residual Squared [mM²]"),
        ("Obj_retentate_concentration", r"$\log_{10}$  Retentate Concentration Residual Squared [mM²]"),
    ]
    cond_specs = [
        ("Obj_permeate_conductivity",   r"$\log_{10}$  Permeate Conductivity Residual Squared [(µS/cm)²]"),
        ("Obj_retentate_conductivity",  r"$\log_{10}$  Retentate Conductivity Residual Squared [(µS/cm)²]"),
    ]
    panel_specs = list(base_specs)
    if all(spec[0] in df.columns for spec in cond_specs):
        panel_specs.extend(cond_specs)

    n_panels = len(panel_specs)
    fig, axes = plt.subplots(1, n_panels, figsize=(14 * n_panels / 3, 7))
    if n_panels == 1:
        axes = [axes]

    for ax, (col, panel_title) in zip(axes, panel_specs):
        piv = df.pivot_table(index=y_var, columns=x_var, values=col, aggfunc="mean")
        X = piv.columns.to_numpy(dtype=float)
        Y = piv.index.to_numpy(dtype=float)
        XX, YY = np.meshgrid(X, Y)
        Z = piv.to_numpy(dtype=float)

        if np.isfinite(Z).any():
            cp = ax.contour(XX, YY, Z, n_levels, linewidths=2)
            ax.clabel(cp, inline=True, fontsize=10, fmt="%.1f")
            flat = np.nanargmin(Z)
            iy, ix = np.unravel_index(flat, Z.shape)
            x_min = float(XX[iy, ix])
            y_min = float(YY[iy, ix])
            z_min = float(Z[iy, ix])
            ax.plot(x_min, y_min, "^", markersize=8,
                    markeredgecolor="red", markerfacecolor=[1, 0.6, 0.6],
                    zorder=10)
            # Optimum text box beneath each panel — same convention as the
            # legacy MATLAB plot_contour.m legend, but rendered as a bordered
            # text box for visibility.  Format: "(theta_1, theta_2, log10_min)".
            opt_str = (f"({x_var} = {x_min:.3g},  "
                       f"{y_var} = {y_min:.3g},  "
                       f"log₁₀(SSR) = {z_min:.2f})")
            ax.text(0.5, -0.20, opt_str,
                    transform=ax.transAxes,
                    ha="center", va="top", fontsize=9.5,
                    bbox=dict(boxstyle="round,pad=0.4",
                              facecolor="white",
                              edgecolor=[0.6, 0.2, 0.2], linewidth=1.0))
        ax.set_xlabel(labels.get(x_var, x_var), fontsize=12)
        ax.set_ylabel(labels.get(y_var, y_var), fontsize=12)
        ax.set_title(panel_title, fontsize=12)
        ax.tick_params(labelsize=10)

    if title:
        fig.suptitle(title, fontsize=14)
    plt.tight_layout(rect=[0, 0.08, 1, 0.96])

    out_path = save_dir / f"objcontour-x_{x_var}-y_{y_var}.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path


def sigma_sensitivity_py(
    data_stru,
    sigma_values,
    theta,
    *,
    save_dir,
    mode="Lag",
    B_form="single",
    workflow_family="DATA3",
    nfe=80,
    LOUD=True,
):
    """Python port of legacy/data1_matlab/functions/sigma_sensitivity.m.

    Sweep σ over a caller-supplied 1D array; for each value, copy theta,
    set theta["sigma"] = σ, run a forward simulation, and store the
    (time, mV, cF, cV) trajectories.  If LOUD, render a 24×5-inch figure
    with three subplots — mass(t), retentate concentration cF(t),
    permeate concentration cV(t) — one colored curve per σ value.

    Matches the MATLAB:
      - Per-σ simulation copies theta, sets theta[sigma] = σ, runs sim_model.
      - 1×3 subplot layout, LineWidth=2.
      - Color cycle 'rbg' for the first 3 σ values (matches MATLAB
        ``colorstring = 'rbg'``); extended cycle for >3 values.
      - Save as: ``sigma_sensitivity.png`` in save_dir.

    Parameters
    ----------
    sigma_values : sequence of float
        σ values to sweep.  Typically 3-6 values for diagnostic visibility.
    theta : dict
        Base parameters (Lp, B, sigma, S0, S — sigma will be overridden).

    Returns
    -------
    tuple (list, Path)
        - List of {'sigma', 'theta', 'sim_stru'} dicts, one per σ.
        - Path to the saved PNG (if LOUD else None).

    Notes
    -----
    MATLAB source: legacy/data1_matlab/functions/sigma_sensitivity.m
                   Lines 16-21 (theta override + sim loop),
                   Lines 24-90 (3-panel time-series plot).
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    out = []
    for sg in sigma_values:
        b = dict(theta)
        b["sigma"] = float(sg)
        try:
            fit_stru, sim_stru, _ = solve_model(
                data_stru,
                mode,
                theta=b,
                sim_opt=True,
                B_form=B_form,
                workflow_family=workflow_family,
                nfe=nfe,
                LOUD=False,
            )
        except Exception as e:
            print(f"  sigma_sensitivity_py: σ={sg:.3f} forward sim raised: {e!r}")
            sim_stru = None
        out.append({"sigma": float(sg), "theta": dict(b), "sim_stru": sim_stru})

    if not LOUD:
        return out, None

    salt_label = data_stru.get("data_config", {}).get("namec", "salt")
    colors = ["red", "blue", "green", "purple", "orange", "brown"]
    fig, axes = plt.subplots(1, 3, figsize=(24, 5))

    for k, rec in enumerate(out):
        sim = rec["sim_stru"]
        if sim is None:
            continue
        color = colors[k % len(colors)]
        label = rf"$\sigma$ = {rec['sigma']:.3f}"
        traces = sim if isinstance(sim, (list, tuple)) else [sim]

        first_in_group = True
        for tr in traces:
            t  = tr.get("time") if isinstance(tr, dict) else getattr(tr, "time", None)
            mV = tr.get("mV")   if isinstance(tr, dict) else getattr(tr, "mV",   None)
            cF = tr.get("cF")   if isinstance(tr, dict) else getattr(tr, "cF",   None)
            cV = tr.get("cV")   if isinstance(tr, dict) else getattr(tr, "cV",   None)
            if t is None:
                continue
            t = np.asarray(t)
            lab = label if first_in_group else None
            first_in_group = False
            if mV is not None:
                axes[0].plot(t, np.asarray(mV), color=color, linewidth=2, label=lab)
            if cF is not None:
                axes[1].plot(t, np.asarray(cF), color=color, linewidth=2, label=lab)
            if cV is not None:
                cV_arr = np.asarray(cV)
                t_for_cv = t[1:] if len(cV_arr) == len(t) - 1 else t
                if len(t_for_cv) == len(cV_arr):
                    axes[2].plot(t_for_cv, cV_arr, color=color, linewidth=2, label=lab)

    axes[0].set_xlabel("Time [s]", fontsize=15)
    axes[0].set_ylabel("Collected Permeate Mass [g]", fontsize=15)
    axes[0].set_title("Mass Predictions", fontsize=15)
    axes[0].tick_params(labelsize=12)

    axes[1].set_xlabel("Time [s]", fontsize=15)
    axes[1].set_ylabel(f"Retentate Concentration of {salt_label} [mmol/L]", fontsize=15)
    axes[1].set_title("Concentration Predictions", fontsize=15)
    axes[1].tick_params(labelsize=12)
    axes[1].legend(fontsize=12, loc="best")

    axes[2].set_xlabel("Time [s]", fontsize=15)
    axes[2].set_ylabel(f"Permeate Concentration of {salt_label} [mmol/L]", fontsize=15)
    axes[2].set_title("Concentration Predictions", fontsize=15)
    axes[2].tick_params(labelsize=12)

    plt.tight_layout()
    out_path = save_dir / "sigma_sensitivity.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out, out_path


def calc_contour_3d_py(
    data_stru,
    theta,
    *,
    grid_density=20,
    sigma_levels=6,
    save_dir,
    b_bounds=None,
    lp_bounds=None,
    sigma_bounds=None,
    mode="Lag",
    B_form="single",
    workflow_family="DATA3",
    nfe=80,
    write_csv=True,
):
    """Python port of legacy/data1_matlab/functions/calc_contour_3d.m.

    3D parameter sweep over (B, Lp, sigma). B and Lp use ``grid_density`` linear
    points each; sigma uses ``sigma_levels`` points (MATLAB hard-codes 6). For
    every (B, Lp, sigma) node, theta is copied, the three params overridden, and
    calc_ind_objectives_py evaluated; log10 of each per-channel SSR is stored.

    Bounds (MATLAB calc_contour_3d.m lines 39-49):
      - Lp    -> [0.1*theta.Lp, 2*theta.Lp]
      - B     -> MATLAB uses [0, 2]; here the salt-specific
                 NF270_B_BOUNDS_PER_SALT triple for DATA3, else (1e-6, 2)
      - sigma -> theta bounds (0, 1)

    Output CSV ``contour3ddata.csv`` columns
    [B, Lp, sigma, Obj_mass, Obj_concentration, Obj_retentate_concentration]
    (MATLAB calc_contour_3d.m lines 98-100, with Obj_perm_conc/Obj_reten_conc
    renamed to Obj_concentration/Obj_retentate_concentration so the existing
    Python contour readers and the 2D schema line up).

    Returns pandas.DataFrame.

    Notes
    -----
    Evaluations = grid_density**2 * sigma_levels forward solves — heavy. Use a
    small grid_density (e.g. 8-12) for previews.
    """
    import numpy as np
    import pandas as pd

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    def _b_default():
        if workflow_family == "DATA3":
            salt = str(data_stru.get("data_config", {}).get("namec", "")).strip()
            return NF270_B_BOUNDS_PER_SALT.get(salt, NF270_B_BOUNDS_DEFAULT)
        # Non-DATA3 fallback mirrors MATLAB calc_contour_3d.m:44-45 (B hard-coded
        # to [0, 2]). DATA3 always uses NF270_B_BOUNDS_PER_SALT above, so DATA3
        # contour outputs are unaffected by this legacy bound.
        return (1e-6, 2.0)

    lp0 = float(theta.get("Lp", 1.0))
    b_lb, b_ub = b_bounds or _b_default()
    lp_lb, lp_ub = lp_bounds or (0.1 * lp0, 2.0 * lp0)
    s_lb, s_ub = sigma_bounds or (0.0, 1.0)

    B_vals = np.linspace(b_lb, b_ub, grid_density)
    Lp_vals = np.linspace(lp_lb, lp_ub, grid_density)
    S_vals = np.linspace(s_lb, s_ub, sigma_levels)

    def _log10_or_nan(v):
        try:
            if v is None or not np.isfinite(v) or v <= 0:
                return float("nan")
            return float(np.log10(v))
        except (TypeError, ValueError):
            return float("nan")

    rows = []
    for bb in B_vals:
        for ll in Lp_vals:
            for ss in S_vals:
                b = dict(theta)
                b["B"] = float(bb)
                b["Lp"] = float(ll)
                b["sigma"] = float(ss)
                ind = calc_ind_objectives_py(
                    b, data_stru,
                    mode=mode, B_form=B_form,
                    workflow_family=workflow_family, nfe=nfe,
                )
                rows.append({
                    "B": float(bb),
                    "Lp": float(ll),
                    "sigma": float(ss),
                    "Obj_mass":                    _log10_or_nan(ind.m),
                    "Obj_concentration":           _log10_or_nan(ind.cp),
                    "Obj_retentate_concentration": _log10_or_nan(ind.cr),
                })

    df = pd.DataFrame(rows)
    if write_csv:
        df.to_csv(save_dir / "contour3ddata.csv", index=False)
    return df


def _apply_sweep_condition(ds, condition_var, cval):
    """Apply a swept experiment condition to a (deep-copied) data_stru so it
    actually changes the simulated trajectory.

    ``C_D`` is read straight from data_config by the Pyomo model (m.cD), so
    setting data_config['C_D'] is sufficient. ``C_F0`` is the initial RETENTATE
    concentration, which the Lag/Overflow model anchors to the first finite
    cF_exp of vial 0 (m.cF[1,0] == firstNonNan(data_raw[0]['cF_exp'])) — the
    data_config['C_F0'] entry alone is only a solver-init hint. So for a C_F0
    sweep we must also override that first-finite cF_exp anchor, otherwise the
    swept value never reaches the IC and the whole concentration axis is inert.
    """
    import numpy as _np
    dc = ds.setdefault("data_config", {})
    dc[condition_var] = float(cval)
    if condition_var == "C_F0" and ds.get("data_raw"):
        row = ds["data_raw"][0]
        cf = row.get("cF_exp")
        if cf is not None:
            arr = _np.asarray(cf, dtype=float)
            flat = arr.reshape(-1)
            fin = _np.where(_np.isfinite(flat))[0]
            if fin.size:
                flat = flat.copy()
                flat[fin[0]] = float(cval)
                row["cF_exp"] = flat.reshape(arr.shape)


def _initial_retentate_cF(ds):
    """First finite retentate concentration of vial 0 (the osmotic driver), used
    as the cF proxy in the Jw feasibility filter; falls back to data_config C_F0."""
    import numpy as _np
    try:
        flat = _np.asarray(ds["data_raw"][0].get("cF_exp"), dtype=float).reshape(-1)
        fin = flat[_np.isfinite(flat)]
        if fin.size:
            return float(fin[0])
    except Exception:
        pass
    return float(ds.get("data_config", {}).get("C_F0", 0.0) or 0.0)


def doe_heatmap_py(
    data_stru,
    theta,
    *,
    condition_var,
    condition_values,
    delp_values_psi,
    mode=None,
    B_form="single",
    workflow_family="DATA3",
    end_vial=None,
    step=1e-3,
    formula="forward",
    nfe=120,
    jw_filter=True,
    save_dir=None,
    write_csv=True,
):
    """Python port of legacy/data1_matlab/doe_heatmap_{filtration,diafiltration}.m.

    Sweeps experiment conditions (``condition_var`` x applied pressure) and, at
    each node, builds the FIM via calc_FIM and records the four MBDoE optimality
    metrics: A (trace), D (det), E (min eigenvalue) and modified-E (cond).

    One Python function covers both MATLAB drivers via ``condition_var``:
      - "C_F0": initial feed concentration sweep (doe_heatmap_filtration, cf0)
      - "C_D" : dialysate concentration sweep    (doe_heatmap_diafiltration, cd)

    ``delp_values_psi`` are converted to bar (/14.504) before entering
    data_config, matching the MATLAB ``delP/14.504``.

    Jw feasibility filter (MATLAB lines ~56-72): a node is dropped to NaN when
    delP - ni*R*T*sigma*cF < 1e-8 (osmotically infeasible). NOTE: MATLAB takes
    max(cF) from a forward sim; to avoid doubling the solve count this port uses
    the swept concentration as the cF proxy (documented approximation).

    Returns dict: condition_var, condition_values, delp,
    A_trace / D_det / E_min_eig / modE_cond matrices shaped
    [len(delp) x len(condition_values)] (MATLAB orientation).
    """
    import numpy as np
    import pandas as pd
    import copy

    mode = mode or str(data_stru.get("mode") or "Lag")
    cond_values = [float(c) for c in condition_values]
    delp_values = [float(p) for p in delp_values_psi]
    nC, nP = len(cond_values), len(delp_values)
    A = np.full((nP, nC), np.nan)
    D = np.full((nP, nC), np.nan)
    E = np.full((nP, nC), np.nan)
    ME = np.full((nP, nC), np.nan)

    R = 8.314e-5  # cm^3 bar / micromol / K
    sig = float(theta.get("sigma", 1.0))

    for ci, cval in enumerate(cond_values):
        for pi, dpsi in enumerate(delp_values):
            ds = copy.deepcopy(data_stru)
            ds["data_config"]["delP"] = dpsi / 14.504
            _apply_sweep_condition(ds, condition_var, cval)
            if end_vial is not None:
                ds["data_config"]["n"] = int(end_vial)
            if jw_filter:
                ni = float(ds["data_config"].get("ni", 1.0))
                T = float(ds["data_config"].get("Temp", 298.15))
                # Retentate osmotic driver: the swept feed conc for a C_F0 sweep,
                # else the data's initial retentate conc (a C_D sweep changes the
                # dialysate, not the retentate driver).
                cF_proxy = cval if condition_var == "C_F0" else _initial_retentate_cF(ds)
                if ds["data_config"]["delP"] - ni * R * T * sig * cF_proxy < 1e-8:
                    continue
            try:
                doe = calc_FIM(
                    ds, mode, theta=theta, step=step, formula=formula,
                    B_form=B_form, workflow_family=workflow_family, nfe=nfe,
                )
            except Exception:
                doe = None
            if not doe:
                continue
            A[pi, ci] = float(doe.get("trace", np.nan))
            D[pi, ci] = float(doe.get("det", np.nan))
            E[pi, ci] = float(doe.get("min_eig", np.nan))
            ME[pi, ci] = float(doe.get("cond", np.nan))

    result = {
        "condition_var": condition_var,
        "condition_values": cond_values,
        "delp": delp_values,
        "A_trace": A.tolist(),
        "D_det": D.tolist(),
        "E_min_eig": E.tolist(),
        "modE_cond": ME.tolist(),
    }
    if write_csv and save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        rows = []
        for ci, cval in enumerate(cond_values):
            for pi, dpsi in enumerate(delp_values):
                rows.append({
                    condition_var: cval, "delP_psi": dpsi,
                    "A_trace": A[pi, ci], "D_det": D[pi, ci],
                    "E_min_eig": E[pi, ci], "modE_cond": ME[pi, ci],
                })
        pd.DataFrame(rows).to_csv(save_dir / f"doe_heatmap-{condition_var}.csv", index=False)
    return result


def heatmap_sigma_sensitivity_py(
    data_stru,
    theta,
    *,
    condition_var,
    condition_values,
    delp_values_psi,
    sigma=(0.9, 1.0),
    mode=None,
    B_form="single",
    workflow_family="DATA3",
    end_vial=None,
    nfe=120,
    scale=True,
    save_dir=None,
    write_csv=True,
):
    """Python port of legacy/data1_matlab/heatmap_sigma_sensitivity_{filtration,diafiltration}.m.

    At each (``condition_var`` x applied pressure) node, forward-simulates the
    model at each sigma in ``sigma`` and records the RANGE (max-min across sigma)
    of the end-of-run mass (mV), retentate (cF) and permeate (cV/cH) predictions.

    With ``scale=True`` the ranges are divided by the measurement uncertainty:
    mass 0.01 g, and retentate/permeate relative fractions from
    _nf270_conc_scales(workflow_family) — DATA3 uses cF 2% / cV 3%, while
    DATA1/DATA2 fall back to 0.3% / 3%. (MATLAB hard-codes 0.3%/3% at lines
    122-127; sourcing the fractions from the DATA3 helper here keeps the
    diagnostic consistent with the DATA3 error model instead of undoing it.)

    Returns dict with range_mV / range_cF / range_cP matrices shaped
    [len(delp) x len(condition_values)].
    """
    import numpy as np
    import pandas as pd
    import copy

    cond_values = [float(c) for c in condition_values]
    delp_values = [float(p) for p in delp_values_psi]
    mode = mode or str(data_stru.get("mode") or "Lag")
    sig_list = [float(s) for s in sigma]
    # DATA3 measurement-error fractions (cF 2% / cV 3%); legacy fallback 0.3%/3%.
    cf_frac, cp_frac = _nf270_conc_scales(workflow_family)
    nC, nP = len(cond_values), len(delp_values)
    range_mV = np.full((nP, nC), np.nan)
    range_cF = np.full((nP, nC), np.nan)
    range_cP = np.full((nP, nC), np.nan)

    def _endval(sim_stru, vial, *keys):
        try:
            tr = sim_stru[vial]
        except Exception:
            return float("nan")
        for key in keys:
            arr = tr.get(key) if isinstance(tr, dict) else getattr(tr, key, None)
            if arr is None:
                continue
            arr = np.asarray(arr, dtype=float).reshape(-1)
            if arr.size:
                return float(arr[-1])
        return float("nan")

    for ci, cval in enumerate(cond_values):
        for pi, dpsi in enumerate(delp_values):
            ds = copy.deepcopy(data_stru)
            ds["data_config"]["delP"] = dpsi / 14.504
            _apply_sweep_condition(ds, condition_var, cval)
            if end_vial is not None:
                ds["data_config"]["n"] = int(end_vial)
            vial_idx = (int(end_vial) - 1) if end_vial is not None else (int(ds["data_config"]["n"]) - 1)
            mV, cF, cP = [], [], []
            for sg in sig_list:
                b = dict(theta)
                b["sigma"] = sg
                try:
                    _, sim_stru, _ = solve_model(
                        ds, mode, theta=b, sim_opt=True,
                        B_form=B_form, workflow_family=workflow_family,
                        nfe=nfe, LOUD=False,
                    )
                except Exception:
                    sim_stru = None
                if not sim_stru:
                    mV.append(np.nan); cF.append(np.nan); cP.append(np.nan)
                    continue
                mV.append(_endval(sim_stru, vial_idx, "mV"))
                cF.append(_endval(sim_stru, vial_idx, "cF"))
                cP.append(_endval(sim_stru, vial_idx, "cH", "cV"))
            mV = np.asarray(mV); cF = np.asarray(cF); cP = np.asarray(cP)
            r_m = float(np.nanmax(mV) - np.nanmin(mV)) if np.isfinite(mV).any() else np.nan
            r_cf = float(np.nanmax(cF) - np.nanmin(cF)) if np.isfinite(cF).any() else np.nan
            r_cp = float(np.nanmax(cP) - np.nanmin(cP)) if np.isfinite(cP).any() else np.nan
            if scale:
                r_m = r_m / 0.01
                mcf = np.nanmean(cF)
                r_cf = r_cf / (mcf * cf_frac) if np.isfinite(mcf) and mcf != 0 else np.nan
                mcp = np.nanmean(cP)
                r_cp = r_cp / (mcp * cp_frac) if np.isfinite(mcp) and mcp != 0 else np.nan
            range_mV[pi, ci] = r_m
            range_cF[pi, ci] = r_cf
            range_cP[pi, ci] = r_cp

    result = {
        "condition_var": condition_var,
        "condition_values": cond_values,
        "delp": delp_values,
        "range_mV": range_mV.tolist(),
        "range_cF": range_cF.tolist(),
        "range_cP": range_cP.tolist(),
    }
    if write_csv and save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        rows = []
        for ci, cval in enumerate(cond_values):
            for pi, dpsi in enumerate(delp_values):
                rows.append({
                    condition_var: cval, "delP_psi": dpsi,
                    "range_mV": range_mV[pi, ci], "range_cF": range_cF[pi, ci],
                    "range_cP": range_cP[pi, ci],
                })
        pd.DataFrame(rows).to_csv(save_dir / f"sigma_sensitivity_heatmap-{condition_var}.csv", index=False)
    return result
