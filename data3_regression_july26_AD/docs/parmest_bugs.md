# ParmEst issues found during DATA3 regression work

This file collects issues, silent footguns, and enhancement ideas we hit while
using `pyomo.contrib.parmest` (Pyomo **6.10.1**) for the NaCl Donnan regression
(`analysis/step3_nacl_uncertainty/parmest_uncertainty.py`). Each entry is written
as a **self-contained prompt**: hand one to a fresh Claude Code (Sonnet) session
pointed at a clone of `pyomo`, and it should be able to reproduce the issue, write
a test, and draft a patch/PR.

## Versions (all issues below were observed with these)

| Package | Version |
|---|---|
| **Pyomo** (includes `pyomo.contrib.parmest`) | **6.10.1** |
| Python | 3.11.15 |
| Ipopt (conda-forge) | 3.14.19 |
| NumPy | 2.4.6 |
| SciPy | 1.17.1 |
| OS | macOS (Darwin arm64) |

ParmEst ships inside Pyomo, so its version is the Pyomo version above. All entries
are reproducible in the `data3-regression` conda env (`environment.yml`).

Severity legend: **[bug]** wrong/misleading result · **[usability]** correct but a
silent footgun · **[enhancement]** performance/feature · **[verify]** needs
confirmation against current API.

---

## Issue 1 — `cov_est()` silently uses a placeholder `measurement_error` as the residual variance  **[usability, CONFIRMED]**

### Problem
ParmEst's `cov_est()` computes the parameter covariance as
`cov = sigma^2 * (J^T J)^{-1}`. Which `sigma^2` it uses depends entirely on the
`measurement_error` Suffix in `label_model()`, and the two modes are **not obvious
and produce silently very different standard errors with no warning**:

- `measurement_error[output] = None` → ParmEst **estimates** the residual variance
  from the fit, `sigma_hat^2 = SSE/(n-p)` (this matches scipy's
  `curve_fit`/`least_squares` convention).
- `measurement_error[output] = <a number>` (e.g. a placeholder `1.0`) → ParmEst
  **uses that number as the variance**, so every reported SE is scaled by
  `sqrt(supplied_variance / sigma_hat^2)`.

A user who sets `measurement_error = 1.0` as a "required placeholder" (a natural
guess from some examples) gets standard errors wrong by a factor of
`1/sqrt(sigma_hat^2)` with no error, warning, or log message. In our NaCl fit this
was a **~3.5x inflation**; in the minimal example below it is **4.26x**.

### Minimum working example (confirmed to reproduce)
Model `y = a*x`, one parameter, 6 noisy points. The point estimate and SSE are
identical for both settings; only `cov_est()` differs.

```python
import numpy as np, pandas as pd
import pyomo.environ as pyo
from pyomo.contrib.parmest.experiment import Experiment
import pyomo.contrib.parmest.parmest as parmest

DATA = pd.DataFrame({"x": [1.,2.,3.,4.,5.,6.],
                     "y": [2.1,3.9,6.2,7.8,10.3,11.7]})  # ~ y = 2x + noise

class LinExp(Experiment):
    def __init__(self, row, meas_err):
        self.row, self.meas_err, self.model = row, meas_err, None
    def create_model(self):
        m = pyo.ConcreteModel()
        m.a = pyo.Var(initialize=2.0, bounds=(-1e3, 1e3)); m.a.fix()
        m.x = pyo.Var(initialize=self.row["x"]); m.x.fix()
        m.y = pyo.Var(initialize=self.row["y"])
        m.con = pyo.Constraint(rule=lambda m: m.y == m.a*m.x)
        self.model = m
    def label_model(self):
        m = self.model
        m.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.experiment_outputs[m.y] = self.row["y"]
        m.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.unknown_parameters[m.a] = pyo.value(m.a)
        m.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.experiment_inputs[m.x] = self.row["x"]
        m.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.measurement_error[m.y] = self.meas_err          # None  OR  1.0
        return m
    def get_labeled_model(self):
        if self.model is None:
            self.create_model(); self.label_model()
        return self.model

def run(meas_err):
    exps = [LinExp(DATA.iloc[i], meas_err) for i in range(len(DATA))]
    pest = parmest.Estimator(exps, obj_function="SSE")
    obj, theta = pest.theta_est()
    se = float(np.sqrt(np.diag(pest.cov_est().to_numpy()))[0])
    return theta["a"], obj*len(DATA), se

for me in (None, 1.0):
    a, sse, se = run(me)
    print(f"measurement_error={me!r}:  a={a:.6f}  SSE={sse:.6f}  SE(a)={se:.6f}")
```

**Observed output:**
```
measurement_error=None:  a=1.993407  SSE=0.276044  SE(a)=0.024631
measurement_error=1.0 :  a=1.993407  SSE=0.276044  SE(a)=0.104828
```
`SE ratio = 4.26 = 1/sqrt(SSE/(n-p)) = 1/sqrt(0.0552)`.

### Expected vs. actual
Both behaviors are individually defensible (estimate the variance vs. use a
supplied one). The problem is that the switch is **silent** and keys off a value
that looks like a placeholder. A user cannot tell from the output which mode ran.

### Fix plan
1. In `cov_est()` (parmest `parmest.py`), detect when `measurement_error` values
   are all equal and were not clearly user-specified, or at minimum **log/emit a
   message** stating which variance mode is in effect
   (`"estimating sigma^2 = SSE/(n-p) = ..."` vs `"using supplied measurement
   variance = ..."`).
2. Document the two modes prominently in the `cov_est()` docstring and the ParmEst
   covariance user guide, with this exact scaling relationship.
3. Consider a keyword argument (e.g. `variance='estimate'|'supplied'`) so the
   intent is explicit rather than inferred from a Suffix value.
4. Add a regression test asserting the `1/sqrt(sigma_hat^2)` scaling relationship
   between the two modes (the MWE above is the test skeleton).

---

## Issue 2 — `cov_est()` rejects a custom `obj_function` callable that `theta_est()` accepts  **[bug, verify]**

### Problem (reported from Step 3; needs a clean MWE)
`parmest.Estimator(exp_list, obj_function=<callable>)` works for `theta_est()`
(as in the older `rooney_biegler` examples), but calling `cov_est()` on the same
estimator raises `ValueError` — `cov_est()` accepts only the two built-in
objective names `'SSE'` and `'SSE_weighted'`. This is an inconsistency: an
estimator that can fit cannot report covariance without being rebuilt with a
string objective.

### Task for the Sonnet session
1. Write a minimal reproducer: build one `Estimator` with a custom
   `obj_function` callable, confirm `theta_est()` succeeds and `cov_est()` raises
   `ValueError`; contrast with `obj_function='SSE'`.
2. Decide the intended contract: either (a) make `cov_est()` support custom
   objectives (documenting the variance assumption), or (b) fail *early* in the
   `Estimator` constructor / `theta_est()` with a clear message that covariance
   requires a built-in objective — not a late `ValueError` from `cov_est()`.
3. Implement and add a test.

---

## Issue 3 — `theta_est(calc_cov=..., cov_n=...)` deprecated in favor of `cov_est()`  **[note]**

Not a bug — recording so we use the current API. In 6.10.1, passing `calc_cov`/
`cov_n` to `theta_est()` raises a `DeprecationWarning`; covariance now lives in
the separate `cov_est()` method. If any of our tutorials/notebooks still call the
old form, update them. No patch needed unless docs elsewhere still show the old
signature (worth a docs sweep).

---

## Issue 4 — `objective_at_theta()` re-solves the full NLP per grid point  **[enhancement]**

### Problem
Building a confidence-region grid with `objective_at_theta(theta_vals)` re-solves
the entire extensive-form NLP (all `n` experiment blocks) once **per grid point**
(~450 ms/point for our 692-experiment NaCl model → several minutes for a 21×21
grid). When the unknown parameters are fixed, evaluating the objective is a
square-system / function evaluation, not an optimization; the per-point solve is
avoidable.

### Task for the Sonnet session
1. Profile `objective_at_theta()` and confirm it invokes the solver per point.
2. Add a fast path: with all `unknown_parameters` fixed, evaluate the objective by
   (a) a single square-system solve of the model constraints, or (b) direct
   expression evaluation where the outputs are explicit — avoiding a full NLP
   optimization per point.
3. Benchmark against the current path on the MWE and add a test asserting numeric
   agreement between fast and slow paths.

---

## Issue 5 — Built-in multi-start?  **[verify]**

Bill/Alex indicated ParmEst **recently added** a multi-start capability. In
6.10.1 this session found none and emulated multi-start by rebuilding the
`Estimator` with different per-experiment initial values and re-running
`theta_est()`. 

### Task for the Sonnet session
1. Check the current ParmEst API (main/dev branch, not just 6.10.1) for a
   built-in multi-start / initialization-sweep feature (search `parmest.py` and
   recent PRs/changelog).
2. If it exists: document the correct call and update our
   `parmest_uncertainty.py` to use it. If it does not: this becomes an
   enhancement request (a `theta_est(multistart=k, bounds=...)` option that
   returns the best of `k` solves and the outcome distribution).
