"""
Copolymer membranes (Estrada's COOH-functionalized "spacer arm" family):
NaCl (1:1) and Na2SO4 (1:2, divalent SO4^2- co-ion) -- a non-adsorbing-salt
positive control (docs/prompts/copolymer_nacl_na2so4.md).

Repeats the NF270 single-salt workflow (raw-signal preprocessing with a
per-experiment fitted conductivity calibration, thin-film concentration
polarization, bounded Donnan regression, multi-start reliability, AR(1)
autocorrelation correction, adsorption-hypothesis diagnostics) on four new
datasets, generalizing rather than forking the existing machinery:
  - analysis/preprocessing/preprocess_nacl.py: reused unmodified except a
    (backward-compatible) generalization of the `Salt` dataclass to support
    a divalent ANION (Na2SO4's SO4^2- co-ion) instead of only a divalent
    CATION (CaCl2/LaCl3) -- see that module's `SALT_NA2SO4` and the
    `z_anion`/`cation_stoich` fields. Verified byte-identical on the
    existing NaCl/CaCl2/LaCl3 regression numbers (z_anion defaults to 1).
  - analysis/step2_nacl/regress_nacl.py: NaCl reuses `f_corrected`
    (the quadratic 1:1 co-ion model) UNCHANGED -- only the transformed
    response `deltaC_m` differs (copolymer ell, Mackie-Meares membrane
    diffusivities from the SI phi_w, not NF270's).
  - analysis/cacl2/regress_cacl2.py: template for the bounded (|chi|,K)
    fit, multi-start reliability check, profile-likelihood CI, sign-runs
    test, and sliding-window effective-|chi| diagnostic -- re-implemented
    here parameterized by an arbitrary model_fn (cacl2's versions are
    hardwired to f_cacl2; per the project guardrails, existing NF270 code
    is not modified, so these are new, generic siblings, not forks of the
    physics).
  - analysis/step8_autocorrelation/ar1_correction.py: `lag1_autocorr`/
    `whiten` (already dataset-agnostic utilities) reused directly; the
    Cochrane-Orcutt loop is re-implemented generically here (ar1_correction's
    own loop is hardwired to the NaCl deltaC_m model), matching the pattern
    already used in analysis/step9_raw_measurement/raw_regression.py and
    analysis/chi_spline/chispline_ar1_whitened.py.

The Na2SO4 (1:2) co-ion cubic (Estrada's reworked equations, NEW to this
project): membrane negatively charged (COOH), co-ion = SO4^2- (divalent),
counter-ion = Na+:
    (c_m^co)^3 + |chi|(c_m^co)^2 + (1/4)|chi|^2 c_m^co - K (c_s^co)^3 = 0,
    K := delta* * delta_counter_circ (lumped, as CaCl2/LaCl3 lump K).
Unlike CaCl2's cubic (no linear term), this one has all three coefficients
present because sulfate enters the Donnan/Boltzmann balance with an extra
|chi| cross-term. Solved by a robust vectorized bisection (see
`co_root_na2so4` for the monotonicity argument that makes bisection --
rather than CaCl2's closed-form Cardano branch -- the right tool here).

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/copolymer/regress_copolymer.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import least_squares
from scipy.stats import f as f_dist

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "preprocessing"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "cacl2"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step8_autocorrelation"))
import plot_style  # noqa: E402
import preprocess_nacl as pp  # noqa: E402
import regress_nacl as step2  # noqa: E402
import regress_cacl2 as cacl2  # noqa: E402
import ar1_correction as s8  # noqa: E402

plot_style.apply_style()

DATA_DIR = REPO_ROOT / "data" / "copolymer"
FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05
RNG_SEED = 20260707
N_MULTISTART = 100
APPLY_CP = True  # standing convention throughout this project; column B of
# the copolymer processed sheets is confirmed BULK (non-CP) by the same
# c_int/kappa_r coefficient-of-variation diagnostic used for NF270 (CV
# 0.020-0.033 here, well inside the "bulk" cluster, vs. ~0.25 for an
# already-CP-corrected NF270 sheet) -- see report-back / README.

# ---------------------------------------------------------------------------
# Membrane constants (from Estrada's spacer-arm manuscript SI,
# Spacer_arm_SI.docx: Table S1 zeta-potential chi, Section S2; Table S3
# water volume fraction phi_w, Section S5). The SI's membrane names have no
# decimal suffix (PEG0, PEG6, ALK1); the DATA3 dataset names (PEG0.1,
# PEG6.5, ALK1.2) are read as batch/replicate identifiers of those base
# membranes (not confirmed in the SI itself -- flagged).
# ---------------------------------------------------------------------------
ELL_ASSUMED_M = 80e-9  # membrane thickness -- NOT reported anywhere in the SI
# for any membrane (checked: Section S5 tabulates only the wet/dry
# thickness RATIO phi_w, not the thicknesses themselves). Working
# assumption: reuse NF270's assumed 80 nm effective transport-length value
# (docs/reports/main.tex sec:reproMATLAB) for all four copolymer datasets,
# flagged explicitly -- biases the absolute Dm/ell-derived K/deltaC_m scale,
# not the (|chi|,K) fit shape or the adsorption-control diagnostics.

PHI_W = {
    "PEG0": 0.41,   # Table S3, measured directly
    "PEG6": 0.51,   # NOT in Table S3 (only Azide/PEG0/PEG4/ALK1 measured);
                    # PEG4's value used as the nearest-family proxy --
                    # flagged, not a measured PEG6 number.
    "ALK1": 0.06,   # Table S3, measured directly
}
PHI_W_IS_PROXY = {"PEG0": False, "PEG6": True, "ALK1": False}

CHI_SI_MM = {
    # Table S1, zeta-potential-derived |chi| [mM] at pH 6.08 (closest to
    # this project's "unadjusted pH ~6" datasets). Sign is negative in the
    # SI (COOH deprotonation); we report the magnitude to match this
    # project's |chi| convention throughout.
    "PEG0": 12.60,
    # PEG6 and ALK1 have NO entry in Table S1 -- only Azide/PEG0/PEG10 were
    # characterized by zeta potential. Absence is real, not a lookup miss.
}

DATASETS = [
    {"key": "PEG0.1-NaCl", "membrane": "PEG0.1", "membrane_si": "PEG0",
     "salt": pp.SALT_NACL, "raw_file": "PEG0.1.xlsx", "raw_sheet": "6.18.24NaCl"},
    {"key": "PEG6.5-NaCl", "membrane": "PEG6.5", "membrane_si": "PEG6",
     "salt": pp.SALT_NACL, "raw_file": "PEG6.5.xlsx", "raw_sheet": "6.27.24NaCl (2)"},
    {"key": "PEG0.1-Na2SO4", "membrane": "PEG0.1", "membrane_si": "PEG0",
     "salt": pp.SALT_NA2SO4, "raw_file": "PEG0.1.xlsx", "raw_sheet": "8.5.24Na2SO4"},
    {"key": "ALK1.2-Na2SO4", "membrane": "ALK1.2", "membrane_si": "ALK1",
     "salt": pp.SALT_NA2SO4, "raw_file": "ALK1.2.xlsx", "raw_sheet": "5.7.25Na2SO4",
     # This raw sheet has NO per-vial conductivity at all (every vial's cond
     # cell is blank) -- fall back to a calibration derived from the
     # processed sheet's own c_int-vs-kappa_r relationship (R^2=0.998).
     "processed_sheet_for_calib": "ALK1.2 5.7.25Na2SO4"},
]
DATASET_COLORS = plot_style.get_dataset_colors([d["key"] for d in DATASETS])

BETA0 = np.array([10.0, 0.3])          # (|chi| [mM], K or delta* [-])
FIT_BOUNDS = ([0.0, 0.0], [np.inf, np.inf])
CHI_LO, CHI_HI = 1e-2, 200.0            # multi-start global search box
K_LO, K_HI = 1e-4, 50.0

N_WINDOW_POINTS = 80
WINDOW_STRIDE = 40
WINDOW_CHI_BOUND = 300.0
WINDOW_RELSE_UNRELIABLE = 2.0


# ---------------------------------------------------------------------------
# Mackie-Meares membrane-phase diffusivities + ambipolar D from a Salt
# ---------------------------------------------------------------------------
def mackie_meares(d_solution, phi_w):
    """D_membrane = D_solution * (phi_w / (2 - phi_w))^2 (Mackie & Meares,
    1955 -- standard tortuosity/obstruction correction for a swollen gel
    membrane phase, parametrized by the water volume fraction)."""
    return d_solution * (phi_w / (2.0 - phi_w)) ** 2


def membrane_ambipolar_D(salt: pp.Salt, phi_w: float) -> float:
    """Ambipolar (Nernst-Hartley) MEMBRANE-phase diffusivity for `salt`,
    built by Mackie-Meares-correcting each ion's solution diffusivity then
    reusing Salt's own ambipolar-D formula (identical algebraic form in
    either phase; only the D's fed in change) -- reuses
    `pp.Salt.solution_D` on a throwaway Salt carrying the membrane-phase
    D's, rather than duplicating the ambipolar formula."""
    d_cation_m = mackie_meares(salt.d_cation_m2_s, phi_w)
    d_anion_m = mackie_meares(salt.d_anion_m2_s, phi_w)
    membrane_salt = pp.Salt(name=salt.name + " (membrane)", z_cation=salt.z_cation,
                             icp_molar_mass=salt.icp_molar_mass, d_cation_m2_s=d_cation_m,
                             z_anion=salt.z_anion, d_anion_m2_s=d_anion_m)
    return membrane_salt.solution_D


PROCESSED_WORKBOOK = DATA_DIR / "Analysis of Select Data.xlsx"


def derive_calibration_from_processed(processed_sheet: str) -> pp.Calibration:
    """Fallback for a raw sheet with NO per-vial conductivity at all (ALK1.2
    5.7.25Na2SO4: every vial's `cond` cell is blank -- confirmed via
    load_raw_sheet, not a parsing bug). Recover an equivalent linear
    calibration kappa = intercept + slope*c_bulk by regressing the
    PROCESSED sheet's own (already-bulk, non-CP -- see the module-level CP
    status note) c_int [col B] against the retentate conductivity [col N]
    across the whole run, the same "recover the embedded/implied
    calibration" move already used for NF270 sheets lacking usable vial
    data (analysis/preprocessing/README.md)."""
    df = pd.read_excel(PROCESSED_WORKBOOK, sheet_name=processed_sheet, header=1)
    df = df.iloc[:, :15]
    c_int = df.iloc[:, 1].to_numpy(dtype=float)
    kappa_r = df.iloc[:, 13].to_numpy(dtype=float)
    mask = np.isfinite(c_int) & np.isfinite(kappa_r)
    c_int, kappa_r = c_int[mask], kappa_r[mask]
    slope, intercept = np.polyfit(c_int, kappa_r, 1)
    pred = intercept + slope * c_int
    ss_res = np.sum((kappa_r - pred) ** 2)
    ss_tot = np.sum((kappa_r - kappa_r.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    return pp.Calibration(slope=slope, intercept=intercept, n_points=len(c_int),
                           r2=r2, source="derived from processed sheet (no vial cond.)")


# ---------------------------------------------------------------------------
# Data loading: raw -> processed, our own per-experiment calibration + CP
# ---------------------------------------------------------------------------
def load_and_process(ds: dict) -> pd.DataFrame:
    salt = ds["salt"]
    raw = pp.load_raw_sheet(DATA_DIR / ds["raw_file"], ds["raw_sheet"])
    calibration = None
    if ds.get("processed_sheet_for_calib"):
        calibration = derive_calibration_from_processed(ds["processed_sheet_for_calib"])
    result = pp.process_experiment(raw, apply_cp=APPLY_CP, salt=salt, calibration=calibration)
    if result.flags:
        raise RuntimeError(f"{ds['key']}: {result.flags}")
    df = result.df.rename(columns={"c_int_mM": "c_int", "c_p_mM": "c_p",
                                    "Jw_m3_m2_s": "Jw"})
    valid = df[["c_int", "c_p", "Jw"]].notna().all(axis=1) & (df["c_int"] > df["c_p"]) \
        & (df["c_int"] > 0) & (df["Jw"] > 0)
    df_valid = df[valid].copy()
    df_valid.attrs["calibration"] = result.calibration
    return df_valid


def build_deltaC_m(df: pd.DataFrame, Dm: float, ell: float = ELL_ASSUMED_M) -> pd.DataFrame:
    """deltaC_m [mM] = J_s*ell/Dm, J_s = c_p*Jw (co-ion flux; no
    stoichiometric multiplier -- both NaCl's Cl- and Na2SO4's SO4^2- are
    exactly 1-per-formula-unit co-ions here, per the module docstring)."""
    df = df.copy()
    J_s = df["c_p"] * df["Jw"]
    df["J_s"] = J_s
    df["deltaC_m"] = J_s * ell / Dm
    return df


# ---------------------------------------------------------------------------
# Models: NaCl reuses step2.f_corrected unmodified; Na2SO4 is new (cubic
# with a linear term, unlike CaCl2's cubic -- see module docstring)
# ---------------------------------------------------------------------------
def co_root_na2so4(beta, c_s, n_iter=60):
    """Unique nonnegative real root of
        t^3 + chi*t^2 + 0.25*chi^2*t - K*c_s^3 = 0
    via vectorized bisection on [0, hi]. This cubic's derivative
    g'(t) = 3t^2 + 2*chi*t + 0.25*chi^2 has discriminant
    (2chi)^2 - 4*3*(0.25chi^2) = chi^2 >= 0, with roots
    t = (-2chi +/- chi)/6 = -chi/6 or -chi/2 -- BOTH <= 0 for chi >= 0 --
    so g is strictly increasing (hence convex-monotonic, well-conditioned
    for bisection) on t >= 0 for any chi >= 0: no closed-form Cardano
    branch is needed here (contrast CaCl2's cubic, which needed one purely
    for speed, not correctness -- its g is also monotonic on t>=0)."""
    chi, K = beta
    c_s = np.atleast_1d(np.asarray(c_s, dtype=float))
    chi = max(float(chi), 0.0)
    K = max(float(K), 0.0)

    def g(t):
        return t ** 3 + chi * t ** 2 + 0.25 * chi ** 2 * t - K * c_s ** 3

    lo = np.zeros_like(c_s)
    hi = np.where(c_s > 0, c_s, 1.0).astype(float)
    for _ in range(50):
        bad = g(hi) < 0
        if not np.any(bad):
            break
        hi = np.where(bad, hi * 2.0, hi)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        pos = g(mid) > 0
        hi = np.where(pos, mid, hi)
        lo = np.where(pos, lo, mid)
    return 0.5 * (lo + hi)


def f_na2so4(beta, c_int, c_p):
    """Prediction = c_co_m(c_int) - c_co_m(c_p); co-ion solution
    concentration = the salt concentration directly (1 SO4^2- per Na2SO4,
    stoichiometry `cs_co = c`, per docs/prompts/copolymer_nacl_na2so4.md --
    unlike CaCl2's `cs_co = 2*c_CaCl2`)."""
    c_co_int = co_root_na2so4(beta, np.asarray(c_int))
    c_co_p = co_root_na2so4(beta, np.asarray(c_p))
    return c_co_int - c_co_p


def model_for(salt: pp.Salt):
    return step2.f_corrected if salt.name == "NaCl" else f_na2so4


def validate_na2so4_cubic(rng, n=200):
    """Sanity check: for random (chi,K,c_s), verify co_root_na2so4 actually
    solves the cubic to high precision (residual ~0) and matches an
    independent numpy.roots solve."""
    chis = rng.uniform(0.0, 100.0, n)
    Ks = rng.uniform(1e-6, 5.0, n)
    css = rng.uniform(1e-2, 200.0, n)
    max_resid = 0.0
    max_rel_err_vs_roots = 0.0
    for chi, K, cs in zip(chis[:50], Ks[:50], css[:50]):
        t = co_root_na2so4((chi, K), np.array([cs]))[0]
        resid = t ** 3 + chi * t ** 2 + 0.25 * chi ** 2 * t - K * cs ** 3
        max_resid = max(max_resid, abs(resid) / max(K * cs ** 3, 1e-12))
        roots = np.roots([1.0, chi, 0.25 * chi ** 2, -K * cs ** 3])
        real_pos = [r.real for r in roots if abs(r.imag) < 1e-6 * max(abs(r.real), 1) and r.real >= -1e-9]
        if real_pos:
            t_ref = max(real_pos)
            if t_ref > 1e-12:
                max_rel_err_vs_roots = max(max_rel_err_vs_roots, abs(t - t_ref) / t_ref)
    return max_resid, max_rel_err_vs_roots


# ---------------------------------------------------------------------------
# Fitting / multi-start / profile-likelihood (generic siblings of cacl2's,
# parameterized by model_fn -- see module docstring for why these are new
# rather than edits to regress_cacl2.py)
# ---------------------------------------------------------------------------
def fit_model_bounded(model_fn, c_int, c_p, y, beta0=BETA0, label="fit"):
    def resid(beta):
        return model_fn(beta, c_int, c_p) - y

    result = least_squares(resid, beta0, method="trf", bounds=FIT_BOUNDS)
    beta_hat = result.x
    residuals = result.fun
    sse = float(np.sum(residuals ** 2))
    n = len(y)
    p = len(beta_hat)
    dof = max(n - p, 1)
    J = result.jac
    JTJ = J.T @ J
    try:
        JTJ_inv = np.linalg.inv(JTJ)
    except np.linalg.LinAlgError:
        JTJ_inv = np.linalg.pinv(JTJ)
    sigma2 = sse / dof
    cov = sigma2 * JTJ_inv
    se = np.sqrt(np.diag(cov))
    return step2.FitResult(label=label, beta=beta_hat, se=se, sse=sse, n=n, p=p, cov=cov)


def multistart_reliability(model_fn, c_int, c_p, y, sse_global, rng, n_starts=N_MULTISTART):
    log_chi = rng.uniform(np.log(CHI_LO), np.log(CHI_HI), n_starts)
    log_K = rng.uniform(np.log(K_LO), np.log(K_HI), n_starts)
    chi0s, K0s = np.exp(log_chi), np.exp(log_K)
    n_global = 0
    for chi0, K0 in zip(chi0s, K0s):
        def resid(beta):
            return model_fn(beta, c_int, c_p) - y
        try:
            res = least_squares(resid, np.array([chi0, K0]), method="trf",
                                 bounds=([1e-8, 1e-8], [1e5, 1e5]))
            sse_hat = float(np.sum(res.fun ** 2))
            if res.success and sse_hat <= sse_global * (1 + 1e-6) + 1e-8:
                n_global += 1
        except Exception:
            pass
    return 100.0 * n_global / n_starts


def f_test_threshold(sse_min, n, p, alpha):
    Fcrit = f_dist.ppf(1 - alpha, p, n - p)
    return sse_min * (1 + (p / (n - p)) * Fcrit)


def profile_likelihood(model_fn, beta_hat, se, sse_min, n, p, c_int, c_p, y, param_idx,
                        alpha=0.05, n_grid=41, span_se=8.0):
    lo_bound = max(beta_hat[param_idx] - span_se * se[param_idx], 1e-6)
    grid = np.linspace(lo_bound, beta_hat[param_idx] + span_se * se[param_idx], n_grid)
    other_idx = 1 - param_idx
    profile_sse = np.zeros(len(grid))
    for i, val in enumerate(grid):
        def resid_fixed(other_val):
            beta = np.zeros(2)
            beta[param_idx] = val
            beta[other_idx] = other_val[0]
            return model_fn(beta, c_int, c_p) - y
        x0 = np.array([max(beta_hat[other_idx], 1e-7)])
        res = least_squares(resid_fixed, x0, method="trf", bounds=([1e-8], [1e5]))
        profile_sse[i] = float(np.sum(res.fun ** 2))
    thresh = f_test_threshold(sse_min, n, p, alpha)
    return grid, profile_sse, thresh


def profile_ci_from_curve(grid, profile_sse, thresh):
    imin = int(np.argmin(profile_sse))
    lo = grid[0]
    for i in range(imin, 0, -1):
        if profile_sse[i - 1] > thresh >= profile_sse[i]:
            lo = np.interp(thresh, [profile_sse[i], profile_sse[i - 1]],
                            [grid[i], grid[i - 1]])
            break
    else:
        lo = grid[0]
    hi = grid[-1]
    for i in range(imin, len(grid) - 1):
        if profile_sse[i] <= thresh < profile_sse[i + 1]:
            hi = np.interp(thresh, [profile_sse[i], profile_sse[i + 1]],
                            [grid[i], grid[i + 1]])
            break
    else:
        hi = grid[-1]
    return lo, hi


# ---------------------------------------------------------------------------
# Adsorption-hypothesis diagnostics (the control test) -- generic siblings
# of cacl2's sign_runs_test/effective_chi_vs_conc
# ---------------------------------------------------------------------------
def sign_runs_test(resid):
    signs = np.sign(resid)
    signs = signs[signs != 0]
    n = len(signs)
    n_pos = int(np.sum(signs > 0))
    n_neg = n - n_pos
    runs = 1 + int(np.sum(signs[1:] != signs[:-1]))
    if n_pos == 0 or n_neg == 0:
        return runs, np.nan, np.nan
    mean_runs = 2 * n_pos * n_neg / n + 1
    var_runs = (2 * n_pos * n_neg * (2 * n_pos * n_neg - n)) / (n ** 2 * (n - 1))
    z = (runs - mean_runs) / np.sqrt(var_runs) if var_runs > 0 else np.nan
    return runs, mean_runs, z


def effective_chi_vs_conc(model_fn, c_int, c_p, y, K_fixed, n_window=N_WINDOW_POINTS,
                           stride=WINDOW_STRIDE, chi0=10.0):
    order = np.argsort(c_int)
    c_int_s, c_p_s, y_s = c_int[order], c_p[order], y[order]
    n = len(c_int_s)
    centers, chi_locals, chi_ses, reliable = [], [], [], []
    start = 0
    while start + n_window <= n:
        sl = slice(start, start + n_window)
        ci, cp, yy = c_int_s[sl], c_p_s[sl], y_s[sl]

        def resid(chi):
            return model_fn((chi[0], K_fixed), ci, cp) - yy

        x0_win = min(max(chi0, 1e-3), WINDOW_CHI_BOUND / 2)
        res = least_squares(resid, np.array([x0_win]), method="trf",
                             bounds=([1e-6], [WINDOW_CHI_BOUND]))
        sse = float(np.sum(res.fun ** 2))
        dof = max(len(yy) - 1, 1)
        sigma2 = sse / dof
        J = res.jac
        JTJ = J.T @ J
        try:
            se = np.sqrt(sigma2 * np.linalg.inv(JTJ)[0, 0])
        except np.linalg.LinAlgError:
            se = np.nan
        chi_hat = float(res.x[0])
        at_bound = chi_hat > 0.995 * WINDOW_CHI_BOUND
        se_ok = np.isfinite(se) and se < WINDOW_RELSE_UNRELIABLE * max(chi_hat, 1e-3)
        centers.append(float(np.mean(ci)))
        chi_locals.append(chi_hat)
        chi_ses.append(float(se))
        reliable.append(bool(se_ok and not at_bound))
        start += stride
    return (np.array(centers), np.array(chi_locals), np.array(chi_ses), np.array(reliable))


# ---------------------------------------------------------------------------
# AR(1) whitening (generic Cochrane-Orcutt, reusing s8.lag1_autocorr/whiten
# -- see module docstring)
# ---------------------------------------------------------------------------
def cochrane_orcutt_generic(model_fn, theta0, c_int, c_p, y, max_iter=50, tol=1e-10):
    def resid_fn(theta):
        return model_fn(theta, c_int, c_p) - y

    theta = np.asarray(theta0, dtype=float)
    r = resid_fn(theta)
    phi = s8.lag1_autocorr(r)
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        def resid_whitened(th):
            return s8.whiten(resid_fn(th), phi)
        res = least_squares(resid_whitened, theta, method="trf", bounds=FIT_BOUNDS)
        theta_new = res.x
        r_new = resid_fn(theta_new)
        phi_new = s8.lag1_autocorr(r_new)
        d_theta = np.max(np.abs(theta_new - theta) / np.maximum(np.abs(theta), 1e-12))
        d_phi = abs(phi_new - phi)
        theta, phi = theta_new, phi_new
        if d_theta < tol and d_phi < tol:
            break
    return theta, phi, n_iter


def whitened_fit_result_generic(model_fn, theta_hat, phi_hat, c_int, c_p, y, label="GLS/AR(1)"):
    def resid(th):
        return s8.whiten(model_fn(th, c_int, c_p) - y, phi_hat)

    res = least_squares(resid, theta_hat, method="trf", bounds=FIT_BOUNDS)
    beta_hat = res.x
    sse = float(np.sum(res.fun ** 2))
    n = len(res.fun)
    p = len(beta_hat)
    dof = max(n - p, 1)
    J = res.jac
    JTJ = J.T @ J
    try:
        JTJ_inv = np.linalg.inv(JTJ)
    except np.linalg.LinAlgError:
        JTJ_inv = np.linalg.pinv(JTJ)
    sigma_a2 = sse / dof
    cov = sigma_a2 * JTJ_inv
    se = np.sqrt(np.diag(cov))
    return step2.FitResult(label=label, beta=beta_hat, se=se, sse=sse, n=n, p=p, cov=cov)


# ---------------------------------------------------------------------------
# Per-dataset analysis
# ---------------------------------------------------------------------------
def compute_sse_surface(model_fn, c_int, c_p, y, chi_grid, k_grid):
    """Raw (non-profiled) 2D SSE(chi,K) grid -- both parameters fixed at
    each grid point, unlike the profile-likelihood curves elsewhere (which
    optimize over the other parameter). This is the surface that directly
    shows a flat ridge/valley as evidence of weak identifiability, per
    docs/prompts/copolymer_nacl_na2so4.md's request for an SSE heatmap."""
    SSE = np.empty((len(chi_grid), len(k_grid)))
    for i, chi in enumerate(chi_grid):
        for j, k in enumerate(k_grid):
            r = model_fn((chi, k), c_int, c_p) - y
            SSE[i, j] = np.nansum(r ** 2)
    return SSE


def analyze_dataset(ds: dict, rng):
    salt = ds["salt"]
    model_fn = model_for(salt)
    phi_w = PHI_W[ds["membrane_si"]]
    Dm = membrane_ambipolar_D(salt, phi_w)

    df_raw = load_and_process(ds)
    df = build_deltaC_m(df_raw, Dm=Dm)
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()
    y = df["deltaC_m"].to_numpy()
    n = len(y)

    beta0 = BETA0 if salt.name != "NaCl" else np.array([10.0, 0.1])
    fr = fit_model_bounded(model_fn, c_int, c_p, y, beta0=beta0,
                            label=f"{ds['key']}: {salt.name}")

    resid = model_fn(fr.beta, c_int, c_p) - y
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    pseudo_r2 = 1.0 - fr.sse / ss_tot if ss_tot > 0 else float("nan")
    rmse = float(np.sqrt(np.mean(resid ** 2)))
    rmse_over_range = rmse / (y.max() - y.min())

    pct_global = multistart_reliability(model_fn, c_int, c_p, y, fr.sse, rng)

    profile_cis = []
    for idx in range(2):
        grid, prof_sse, thresh = profile_likelihood(
            model_fn, fr.beta, fr.se, fr.sse, n, fr.p, c_int, c_p, y, idx, alpha=ALPHA,
        )
        profile_cis.append(profile_ci_from_curve(grid, prof_sse, thresh))

    order = np.argsort(c_int)
    runs, mean_runs, z_runs = sign_runs_test(resid[order])

    centers, chi_locals, chi_ses, reliable = effective_chi_vs_conc(
        model_fn, c_int, c_p, y, K_fixed=fr.beta[1], chi0=fr.beta[0],
    )

    # AR(1)
    phi_ols = s8.lag1_autocorr(resid)
    theta_w, phi_hat, n_iter = cochrane_orcutt_generic(model_fn, fr.beta, c_int, c_p, y)
    fr_w = whitened_fit_result_generic(model_fn, theta_w, phi_hat, c_int, c_p, y)
    n_eff = n * (1 - phi_hat) / (1 + phi_hat)
    z95 = 1.959963984540054
    wald_iid = [(fr.beta[k] - z95 * fr.se[k], fr.beta[k] + z95 * fr.se[k]) for k in range(2)]
    wald_w = [(fr_w.beta[k] - z95 * fr_w.se[k], fr_w.beta[k] + z95 * fr_w.se[k]) for k in range(2)]
    fold_wald = [(wald_w[k][1] - wald_w[k][0]) / (wald_iid[k][1] - wald_iid[k][0]) for k in range(2)]

    # SSE(chi,K) surface: a wide, adaptive log-log grid centered on the
    # fitted point (+/- 3 decades) -- wide enough to show whether the
    # optimum sits in a sharp bowl (well-identified) or a flat ridge
    # (weakly/non-identified), per docs/prompts/copolymer_nacl_na2so4.md.
    n_grid = 45
    chi_grid = np.logspace(np.log10(max(fr.beta[0], 1e-3)) - 3,
                            np.log10(max(fr.beta[0], 1e-3)) + 3, n_grid)
    k_grid = np.logspace(np.log10(max(fr.beta[1], 1e-8)) - 3,
                          np.log10(max(fr.beta[1], 1e-8)) + 3, n_grid)
    sse_surface = compute_sse_surface(model_fn, c_int, c_p, y, chi_grid, k_grid)

    return dict(
        ds=ds, df=df, fr=fr, resid=resid, pseudo_r2=pseudo_r2, rmse_over_range=rmse_over_range,
        pct_global=pct_global, profile_cis=profile_cis, runs=runs, mean_runs=mean_runs,
        z_runs=z_runs, chi_window_centers=centers, chi_window_locals=chi_locals,
        chi_window_ses=chi_ses, chi_window_reliable=reliable, phi_w=phi_w, Dm=Dm,
        phi_ols=phi_ols, phi_hat=phi_hat, n_eff=n_eff, fr_w=fr_w, wald_iid=wald_iid,
        wald_w=wald_w, fold_wald=fold_wald, n_iter_co=n_iter,
        chi_si_mM=CHI_SI_MM.get(ds["membrane_si"]),
        phi_w_is_proxy=PHI_W_IS_PROXY[ds["membrane_si"]],
        chi_grid=chi_grid, k_grid=k_grid, sse_surface=sse_surface,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def print_summary_table(results):
    print("=" * 130)
    print("SUMMARY: Copolymer NaCl/Na2SO4 regression, non-adsorbing-salt control")
    print("=" * 130)
    header = (f"{'Dataset':<18}{'Salt':<8}{'n':>6}{'|chi| (mM)':>13}{'95% CI':>20}"
              f"{'K/delta*':>11}{'SSE':>12}{'R2':>8}{'MS%':>7}{'runs(z)':>12}")
    print(header)
    for r in results:
        fr = r["fr"]
        lo1, hi1 = r["profile_cis"][0]
        print(f"{r['ds']['key']:<18}{r['ds']['salt'].name:<8}{fr.n:>6}{fr.beta[0]:>13.4f}"
              f"  ({lo1:6.3f},{hi1:6.3f})"
              f"{fr.beta[1]:>11.5f}{fr.sse:>12.4f}{r['pseudo_r2']:>8.4f}{r['pct_global']:>6.1f}%"
              f"  {r['runs']:>3d}({r['z_runs']:+.2f})")
    print()

    print("-" * 130)
    print("AR(1) autocorrelation correction (Method A, GLS/whitening, n-p df)")
    print("-" * 130)
    for r in results:
        print(f"{r['ds']['key']:<18} phi_OLS={r['phi_ols']:.4f} -> CO phi_hat={r['phi_hat']:.4f} "
              f"({r['n_iter_co']} iters), n_eff={r['n_eff']:.1f}/{r['fr'].n} "
              f"({100*r['n_eff']/r['fr'].n:.1f}%), Wald fold |chi|={r['fold_wald'][0]:.2f}x "
              f"K/delta*={r['fold_wald'][1]:.2f}x")
    print()

    print("-" * 130)
    print("Validation vs. independent zeta-potential |chi| (SI Table S1)")
    print("-" * 130)
    for r in results:
        chi_hat = r["fr"].beta[0]
        chi_si = r["chi_si_mM"]
        if chi_si is None:
            print(f"{r['ds']['key']:<18} regressed |chi|={chi_hat:.3f} mM; "
                  f"NO zeta-potential chi in the SI for {r['ds']['membrane_si']} "
                  f"(only Azide/PEG0/PEG10 characterized).")
        else:
            ratio = chi_hat / chi_si
            print(f"{r['ds']['key']:<18} regressed |chi|={chi_hat:.3f} mM vs. "
                  f"zeta-potential |chi|={chi_si:.2f} mM ({r['ds']['membrane_si']}, pH 6.08) "
                  f"-> ratio={ratio:.2f}x")
    print()

    print("-" * 130)
    print("Effective |chi| drift (sliding window, K/delta* fixed at global estimate)")
    print("-" * 130)
    for r in results:
        centers = r["chi_window_centers"]
        locals_ = r["chi_window_locals"]
        rel = r["chi_window_reliable"]
        n_rel = int(rel.sum())
        if n_rel < 2:
            print(f"{r['ds']['key']:<18} insufficient reliable windows ({n_rel}/{len(locals_)})")
            continue
        slope = np.polyfit(centers[rel], locals_[rel], 1)[0]
        print(f"{r['ds']['key']:<18} n_windows={len(locals_):>3} ({n_rel} reliable)  "
              f"chi range [{locals_[rel].min():.2f}, {locals_[rel].max():.2f}] mM  "
              f"slope={slope:+.4f} mM/mM  "
              f"({'DRIFTS with concentration' if abs(slope) > 0.05 else 'roughly flat (as expected, non-adsorbing)'})")
    print()


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_all_fits(results, outfile):
    w = plot_style.fig_width("full")
    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(w, w / n * 1.05))
    for ax, r in zip(axes, results):
        fr = r["fr"]
        model_fn = model_for(r["ds"]["salt"])
        df = r["df"]
        c_int = df["c_int"].to_numpy()
        c_p = df["c_p"].to_numpy()
        y = df["deltaC_m"].to_numpy()
        order = np.argsort(c_int)
        y_fit = model_fn(fr.beta, c_int[order], c_p[order])
        color = DATASET_COLORS[r["ds"]["key"]]
        ax.plot(c_int, y, "o", color=color, markersize=3, alpha=0.4, label="data")
        ax.plot(c_int[order], y_fit, "--", color="k", linewidth=2, label="fit")
        ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
        ax.set_title(f"{r['ds']['key']}", fontsize=9)
        ax.set_box_aspect(1)
    axes[0].set_ylabel(r"$\Delta \boldsymbol{c_m}$ [mM]")
    shared_handles = [
        Line2D([], [], marker="o", linestyle="none", color="gray", alpha=0.6, label="data"),
        Line2D([], [], color="k", linestyle="--", linewidth=2, label="fit"),
    ]
    fig.legend(handles=shared_handles, loc="lower center", bbox_to_anchor=(0.5, -0.12),
               ncol=2, fontsize=9, borderaxespad=0)
    fig.suptitle("Copolymer NaCl/Na$_2$SO$_4$ Donnan model fit vs. data")
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_estimates_vs_zeta(results, outfile):
    """Log-scale y-axis: the fitted |chi| and its 95% CI span nearly SIX
    orders of magnitude across these four datasets (PEG0.1-Na2SO4's ~1 mM
    to PEG6.5-NaCl's ~5x10^5 mM CI upper bound) -- a linear axis makes
    every point but PEG6.5-NaCl collapse to the baseline, which would
    misleadingly hide rather than show the identifiability story."""
    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 1.1))
    labels = [r["ds"]["key"] for r in results]
    x = np.arange(len(labels))
    colors = [DATASET_COLORS[label] for label in labels]
    for xi, r, c in zip(x, results, colors):
        fr = r["fr"]
        lo_ci = max(r["profile_cis"][0][0], 1e-3)  # floor: log axis needs > 0
        lo = max(fr.beta[0] - lo_ci, 1e-6)
        hi = max(r["profile_cis"][0][1] - fr.beta[0], 1e-6)
        ax.errorbar(xi, fr.beta[0], yerr=[[lo], [hi]], fmt="o", color=c,
                    capsize=5, markersize=8, elinewidth=2, label="regressed" if xi == 0 else None)
        if r["chi_si_mM"] is not None:
            ax.plot(xi, r["chi_si_mM"], "*", color="k", markersize=16,
                    label=r"zeta-potential $\chi$" if xi == 0 else None)
    ax.set_yscale("log")
    ax.set_ylabel(r"$\boldsymbol{|\chi|}$ [mM] (log scale)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=1, fontsize=8)
    fig.suptitle("Regressed vs. zeta-potential $\\boldsymbol{|\\chi|}$\n"
                 "(log scale: PEG6.5-NaCl's CI spans nearly 2 decades)")
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_residuals_and_drift(results, outfile):
    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(1, 2, figsize=(w, w * 0.5))
    ax = axes[0]
    for r in results:
        label = r["ds"]["key"]
        color = DATASET_COLORS[label]
        c_int = r["df"]["c_int"].to_numpy()
        ax.plot(c_int, r["resid"], "o", color=color, markersize=3, alpha=0.4, label=label)
    ax.axhline(0, color="k", lw=1, ls=":")
    ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"residual [mM]")
    ax.set_title("(a) residual vs. concentration", fontsize=10)
    ax.set_box_aspect(1)

    ax = axes[1]
    for r in results:
        label = r["ds"]["key"]
        color = DATASET_COLORS[label]
        centers = r["chi_window_centers"]
        locals_ = r["chi_window_locals"]
        ses = r["chi_window_ses"]
        rel = r["chi_window_reliable"]
        if len(centers) == 0:
            continue
        if np.any(rel):
            ax.errorbar(centers[rel], locals_[rel], yerr=ses[rel], fmt="o-", color=color,
                        markersize=4, capsize=2, elinewidth=1, label=label)
        if np.any(~rel):
            ax.plot(centers[~rel], locals_[~rel], "x", color=color, markersize=5, alpha=0.35)
        ax.axhline(r["fr"].beta[0], color=color, ls=":", lw=1, alpha=0.6)
    ax.set_xlabel(r"window-center $\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"local $\boldsymbol{|\chi|}$ [mM]")
    ax.set_title("(b) sliding-window effective $|\\chi|$", fontsize=10)
    ax.set_box_aspect(1)

    # panel (a) always has a labeled artist per dataset; panel (b) drops
    # PEG6.5-NaCl (0/59 reliable windows), so pull the shared legend from (a).
    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="lower center", bbox_to_anchor=(0.5, -0.18),
               ncol=4, fontsize=8, borderaxespad=0)
    fig.suptitle("Adsorption-hypothesis diagnostics (control test: expect flat, "
                 "structure-free residuals)")
    fig.tight_layout(rect=(0, 0.12, 1, 0.93))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_sse_heatmaps(results, outfile):
    """SSE(chi,K) heatmaps, log-log axes, one panel per dataset: log10(SSE)
    color-mapped, global optimum marked with a star. A sharp, localized
    bowl means (chi,K) is well-identified; a broad, flat valley/ridge --
    visually elongated with little color change along it -- means the data
    cannot pin the two parameters down individually (only some combination
    of them), which is the direct visual counterpart of the SSE-profile
    check reported in the text (e.g. PEG6.5-NaCl's SSE changes <1% over a
    500x range in chi)."""
    w = plot_style.fig_width("full")
    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(w, w / n * 1.15))
    for ax, r in zip(axes, results):
        chi_grid, k_grid, sse = r["chi_grid"], r["k_grid"], r["sse_surface"]
        sse_min = np.nanmin(sse)
        log_ratio = np.log10(np.maximum(sse / sse_min, 1.0))
        pcm = ax.pcolormesh(chi_grid, k_grid, log_ratio.T, shading="auto",
                             cmap="viridis_r", vmin=0, vmax=max(np.nanmax(log_ratio), 0.5))
        ax.contour(chi_grid, k_grid, log_ratio.T, levels=[np.log10(1.05)],
                   colors="white", linewidths=1.2)
        ax.plot(r["fr"].beta[0], r["fr"].beta[1], "*", color="red", markersize=14,
                markeredgecolor="white", markeredgewidth=0.5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\boldsymbol{|\chi|}$ [mM]")
        ax.set_title(r["ds"]["key"], fontsize=9)
        ax.set_box_aspect(1)
        cbar = fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.06)
        cbar.set_label(r"$\log_{10}(\mathrm{SSE}/\mathrm{SSE}_{\min})$", fontsize=7)
    axes[0].set_ylabel(r"$\boldsymbol{K}$ (or $\boldsymbol{\delta^*}$) [$-$]")
    fig.suptitle("SSE$(|\\chi|,K)$ surfaces: red star = fitted optimum, white contour = "
                 "5% above minimum\n(a broad pale band along the contour signals weak "
                 "identifiability, not a fitting failure)")
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_ar1_regions(results, outfile):
    w = plot_style.fig_width("full")
    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(w, w / n * 1.1))
    for ax, r in zip(axes, results):
        fr, fr_w = r["fr"], r["fr_w"]
        color = DATASET_COLORS[r["ds"]["key"]]
        ax.errorbar([0], [fr.beta[0]],
                    yerr=[[fr.beta[0] - r["wald_iid"][0][0]], [r["wald_iid"][0][1] - fr.beta[0]]],
                    fmt="o", color=color, capsize=5, markersize=7, label="i.i.d. Wald")
        ax.errorbar([1], [fr_w.beta[0]],
                    yerr=[[max(fr_w.beta[0] - r["wald_w"][0][0], 0)],
                          [max(r["wald_w"][0][1] - fr_w.beta[0], 0)]],
                    fmt="s", color="k", capsize=5, markersize=7, label="AR(1)-whitened Wald")
        ax.set_xlim(-0.5, 1.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["i.i.d.", "AR(1)"])
        ax.set_title(f"{r['ds']['key']}\n$\\hat\\phi$={r['phi_hat']:.2f}, "
                     f"fold={r['fold_wald'][0]:.2f}x", fontsize=8)
        ax.set_box_aspect(1)
    axes[0].set_ylabel(r"$\boldsymbol{|\chi|}$ [mM] $\pm$ 95% Wald")
    fig.suptitle("AR(1)-whitened vs. i.i.d. 95% Wald intervals")
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 130)
    print("Copolymer membranes: NaCl (1:1) and Na2SO4 (1:2) non-adsorbing-salt control")
    print("=" * 130)
    print(f"APPLY_CP = {APPLY_CP}")

    rng = np.random.default_rng(RNG_SEED)

    max_resid, max_rel_err = validate_na2so4_cubic(rng)
    print(f"Na2SO4 cubic root-solver sanity check: max relative cubic residual = "
          f"{max_resid:.2e}, max relative error vs. numpy.roots = {max_rel_err:.2e}")
    print()

    results = []
    for ds in DATASETS:
        print(f"--- {ds['key']} ({ds['raw_file']} / {ds['raw_sheet']}) ---")
        t0 = time.time()
        r = analyze_dataset(ds, rng)
        print(f"  n={r['fr'].n}, phi_w={r['phi_w']:.2f}"
              f"{' [PROXY]' if r['phi_w_is_proxy'] else ''}, Dm={r['Dm']:.3e} m^2/s, "
              f"done in {time.time()-t0:.1f}s")
        results.append(r)
    print()

    print_summary_table(results)

    plot_all_fits(results, FIG_DIR / "copolymer_fits.png")
    plot_estimates_vs_zeta(results, FIG_DIR / "copolymer_zeta_validation.png")
    plot_residuals_and_drift(results, FIG_DIR / "copolymer_adsorption_diagnostics.png")
    plot_ar1_regions(results, FIG_DIR / "copolymer_ar1.png")
    plot_sse_heatmaps(results, FIG_DIR / "copolymer_sse_heatmaps.png")

    return results


if __name__ == "__main__":
    main()
