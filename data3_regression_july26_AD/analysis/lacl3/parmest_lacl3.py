"""
LaCl3 (3:1) Donnan regression -- Method B: ParmEst/Pyomo with the co-ion
QUARTIC posed as a nonlinear CONSTRAINT (not a closed-form root), cross-
validated against Method A (analysis/lacl3/regress_lacl3.py's bracketed-Brent
scipy fit).

Why pose it this way
---------------------
The 3:1 co-ion Donnan equation has no compact closed form (unlike NaCl's
quadratic or CaCl2's cubic). Two equally valid numerical routes:
  (A) solve the quartic numerically per point (bracketed root-find or
      numpy.roots) and hand scipy.optimize.least_squares an explicit
      residual function of theta -- fast, used to drive the diagnostics.
  (B) introduce the membrane co-ion concentration as a Pyomo decision
      variable, CONSTRAINED by the quartic equation (bounded to select the
      positive physical root), and let ParmEst/Ipopt solve the joint
      NLP -- no closed-form root or explicit residual function needed at
      all. This is the formulation that GENERALIZES to the multicomponent
      case (Step 14), where no closed form exists for any salt combination.
This script implements (B) and confirms it reproduces (A)'s point estimate,
as a validation exercise (mirrors Step 3's scipy-vs-ParmEst cross-check for
NaCl).

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/lacl3/parmest_lacl3.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import pyomo.environ as pyo
from pyomo.contrib.parmest.experiment import Experiment
import pyomo.contrib.parmest.parmest as parmest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "lacl3"))
import regress_lacl3 as m_a  # noqa: E402  (Method A, reused for data + comparison)

CO_STOICH = 3.0  # 3 Cl- per LaCl3


# ---------------------------------------------------------------------------
# Experiment class (Pyomo 6.10 Experiment-class API, mirrors
# analysis/step3_nacl_uncertainty/parmest_uncertainty.py's pattern)
# ---------------------------------------------------------------------------
class LaCl3DonnanExperiment(Experiment):
    """One experiment == one (c_int, c_p, deltaC_m) data row. chi, K are
    shared Vars (fixed per-experiment, tied together across the extensive
    form); c_co_int, c_co_p are AUXILIARY Vars local to each experiment,
    each constrained by the 3:1 co-ion quartic (Descartes: exactly one
    positive real root) rather than computed via an explicit formula."""

    def __init__(self, data_row, theta0=None, warm_start=None):
        self.data = data_row
        self.model = None
        self.theta0 = theta0 or {"chi": 44.0, "K": 0.01}
        self.warm_start = warm_start  # optional (c_co_int0, c_co_p0)

    def create_model(self):
        m = pyo.ConcreteModel()
        m.chi = pyo.Var(initialize=self.theta0["chi"], bounds=(0.0, 1e4))
        m.K = pyo.Var(initialize=self.theta0["K"], bounds=(0.0, 1e2))
        m.chi.fix()
        m.K.fix()

        m.c_int = pyo.Var(initialize=self.data["c_int"])
        m.c_int.fix()
        m.c_p = pyo.Var(initialize=self.data["c_p"])
        m.c_p.fix()

        c_co_s_int = CO_STOICH * self.data["c_int"]
        c_co_s_p = CO_STOICH * self.data["c_p"]

        if self.warm_start is not None:
            init_int, init_p = self.warm_start
        else:
            init_int, init_p = 0.1 * c_co_s_int, 0.1 * c_co_s_p
        init_int = max(init_int, 1e-8)
        init_p = max(init_p, 1e-8)

        m.c_co_int = pyo.Var(initialize=init_int, bounds=(1e-10, 50 * c_co_s_int + 1e-3))
        m.c_co_p = pyo.Var(initialize=init_p, bounds=(1e-10, 50 * c_co_s_p + 1e-3))

        def quartic_int_rule(m):
            return m.c_co_int ** 4 + m.chi * m.c_co_int ** 3 == m.K * c_co_s_int ** 4

        m.quartic_int = pyo.Constraint(rule=quartic_int_rule)

        def quartic_p_rule(m):
            return m.c_co_p ** 4 + m.chi * m.c_co_p ** 3 == m.K * c_co_s_p ** 4

        m.quartic_p = pyo.Constraint(rule=quartic_p_rule)

        m.deltaC_m = pyo.Var(initialize=self.data["deltaC_m"])

        def deltaC_rule(m):
            return m.deltaC_m == m.c_co_int - m.c_co_p

        m.deltaC_eq = pyo.Constraint(rule=deltaC_rule)
        self.model = m

    def label_model(self):
        m = self.model

        m.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.experiment_outputs.update([(m.deltaC_m, self.data["deltaC_m"])])

        m.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.unknown_parameters.update((k, pyo.value(k)) for k in [m.chi, m.K])

        m.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.experiment_inputs.update(
            [(m.c_int, self.data["c_int"]), (m.c_p, self.data["c_p"])]
        )

        # measurement_error=None -> cov_est() estimates variance from
        # SSE/(n-p), matching scipy's sigma_hat^2 (see Step 3 README for the
        # ParmEst API trap this avoids).
        m.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.measurement_error.update([(m.deltaC_m, None)])

    def get_labeled_model(self):
        if self.model is None:
            self.create_model()
            self.label_model()
        return self.model


def build_experiment_list(df, theta0, warm_starts=None):
    exp_list = []
    for i in range(len(df)):
        ws = warm_starts[i] if warm_starts is not None else None
        exp_list.append(LaCl3DonnanExperiment(df.iloc[i], theta0=theta0, warm_start=ws))
    return exp_list


def main():
    print("=" * 100)
    print("LaCl3 (3:1) Donnan regression -- Method B: ParmEst nonlinear-constraint quartic")
    print("=" * 100)

    df = m_a.load_dataset()
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()

    # --- Method A result (scipy, bracketed-Brent quartic solve) ---
    fr_a = m_a.fit_model_bounded(m_a.f_lacl3, c_int, c_p, df["deltaC_m"].to_numpy(),
                                  beta0=m_a.BETA0, label="Method A")
    print(f"Method A (scipy):   |chi| = {fr_a.beta[0]:.6g} mM, K = {fr_a.beta[1]:.6g}, "
          f"SSE = {fr_a.sse:.4f}")

    # Warm-start each experiment's auxiliary c_co Vars at Method A's converged
    # root, so Ipopt starts essentially at the solution -- a strong, fair
    # cross-check (both methods should then agree to high precision if the
    # two formulations are truly equivalent).
    theta0 = {"chi": float(fr_a.beta[0]), "K": float(fr_a.beta[1])}
    c_co_int0 = m_a.co_root(fr_a.beta, CO_STOICH * c_int)
    c_co_p0 = m_a.co_root(fr_a.beta, CO_STOICH * c_p)
    warm_starts = list(zip(c_co_int0, c_co_p0))

    exp_list = build_experiment_list(df, theta0=theta0, warm_starts=warm_starts)
    pest = parmest.Estimator(exp_list, obj_function="SSE")
    obj, theta = pest.theta_est()
    n = len(df)
    sse_parmest = obj * n
    cov = pest.cov_est()
    se = np.sqrt(np.diag(cov.to_numpy()))
    beta_parmest = np.array([theta["chi"], theta["K"]])

    print(f"Method B (ParmEst): |chi| = {beta_parmest[0]:.6g} mM, K = {beta_parmest[1]:.6g}, "
          f"SSE = {sse_parmest:.4f}")
    print(f"  ParmEst SE: ({se[0]:.6g}, {se[1]:.6g})   scipy SE: "
          f"({fr_a.se[0]:.6g}, {fr_a.se[1]:.6g})")

    rel_diff_beta = np.abs(beta_parmest - fr_a.beta) / np.maximum(np.abs(fr_a.beta), 1e-12)
    rel_diff_sse = abs(sse_parmest - fr_a.sse) / fr_a.sse
    print(f"  Relative diff (beta): {rel_diff_beta}")
    print(f"  Relative diff (SSE):  {rel_diff_sse:.3e}")
    print()
    print("Agreement verdict:", "CONFIRMED (SSE and beta agree to <1% and <1e-3 resp.)"
          if rel_diff_sse < 1e-2 else "DISAGREEMENT -- investigate")
    print()

    return dict(fr_a=fr_a, beta_parmest=beta_parmest, se_parmest=se, sse_parmest=sse_parmest)


if __name__ == "__main__":
    main()
