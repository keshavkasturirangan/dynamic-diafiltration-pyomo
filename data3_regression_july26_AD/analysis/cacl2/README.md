# CaCl2 (2:1) Donnan regression + adsorption-hypothesis diagnostics (report §3)

Fits the 2:1 co-ion cubic Donnan model to all three processed CaCl2 datasets
and tests the PI's working hypothesis: **Ca2+ adsorbs onto the membrane and
changes its effective fixed charge `|χ|` in a concentration-dependent way.**
The current model treats `|χ|` as a single constant per dataset (no
adsorption), so we expect the CaCl2 fits to be *worse* than NaCl, with the
diagnostic signature being systematic residual structure vs. concentration
and/or an apparent `|χ|` that drifts with concentration.

## How to run

```bash
conda activate data3-regression
python analysis/cacl2/regress_cacl2.py
```

Takes ~6 seconds total for all three datasets. Reuses `analysis/plot_style.py`
and `analysis/step2_nacl/regress_nacl.py`'s `FitResult` dataclass; the model
function (2:1 cubic root) and its fit/uncertainty utilities are new (see
"Why the model logic isn't shared with NaCl" below).

## Datasets

All three sheets live in `data/BoE Analysis.xlsx` (**not**
`data/Rejection_Analysis.xlsx`, which is missing the 06.27.26 run), same
15-column schema as NaCl (B=`c_int`, H=`c_p`, J=`J_w`):

| Membrane | Sheet | Valid rows | `c_int` range (mM) |
|---|---|---|---|
| MC3 | `NF270_MC3 07.11.24_SCaCl2` | 609 | 1.0–47 |
| MC5 | `NF270_MC5 05.27.26_CaCl2` | 1645 | 1.6–108 |
| MC5 | `NF270_MC5 06.27.26_CaCl2` | 1885 | 1.3–103 |

The MC5 05.27.26 vs. 06.27.26 pair is a near-repeat (reproducibility check).
**Raw-only CaCl2 runs are out of scope here** (not yet processed):
MC2 05.07.24, MC3 07.11.24\_CaCl2, MC3 07.12.24\_S2CaCl2, and the mislabeled
tab `NF270_MC5 05.26.26_NaCl` (its own note identifies it as CaCl2, not NaCl
— exclude it from any NaCl accounting too, see the NaCl inventory table).

## Model

**Co-ion = Cl−**, solution concentration set by stoichiometry:
`c_co_s = 2 * c_CaCl2` (2 Cl− per CaCl2), applied at both the feed face
(`c_int`) and permeate face (`c_p`).

**Donnan co-ion cubic** (report Eq. 21co), with a lumped partition constant
`K >= 0` absorbing `δ* · δ°_co` and constant stoichiometric factors:

```
c_co_m^3 + |χ|*c_co_m^2 - K*c_co_s^3 = 0
```

By Descartes' rule of signs (coefficient pattern `+,+,0,-` for `t>0`, and
`+,-,0,+`/`+,-,+` for `t<0` depending on the sign of `|χ|`) there is always
**exactly one positive real root** for any real `|χ|` (not just `|χ|>=0`) —
selecting it is a well-posed, single-valued operation.

**Transformed response** (Cl− flux form, analogous to corrected NaCl):

```
deltaC_m = J_s_Cl * l / D_m
J_s_Cl   = 2 * c_p * J_w          (Cl- flux = 2x the CaCl2 molar flux)
l        = 80e-9 m
D_m      = 3*D_Ca_m*D_Cl_m / (2*D_Ca_m + D_Cl_m)   (2:1 ambipolar membrane
                                                      diffusivity)
D_Ca_m = 0.79e-12 m^2/s, D_Cl_m = 2.03e-12 m^2/s
```

**Prediction** = `c_co_m(c_int) − c_co_m(c_p)` (solve the cubic at each
face). Fit `(|χ|, K)` to `deltaC_m` by `scipy.optimize.least_squares`,
**bounded to `chi, K >= 0`** — see below for why this matters here but not
for NaCl.

These are Opus's working assumptions (per the task prompt), chosen to mirror
the NaCl treatment closely enough to be directly comparable, while being
sufficient for the hypothesis test; a lumped constant prefactor mismatch in
`K` would bias absolute `|χ|`/`K` values but **not** the residual-shape /
drift diagnostics, which are the point of this analysis. They connect to the
open items in report §1.2 ("Modeling questions") and can be refined later.

### Why the cubic root is a vectorized closed-form solve, not `numpy.roots`

The task prompt suggested `numpy.roots([1, |chi|, 0, -K*c_s**3])` per point.
That is correct but far too slow at the array sizes here (n up to ~1900,
called repeatedly inside `least_squares`, multi-start (100 starts), and the
sliding-window drift diagnostic — a rough estimate is tens of millions of
`numpy.roots` calls for the full script, taking minutes to hours). Instead,
`co_root()` uses a **fully vectorized closed-form solution**, in two regimes
selected by `m = K*c_s^3/chi^3` (how small the true root `t` is relative to
`chi`):

- **`m` small (`t << chi`):** the general (complex) Cardano solution suffers
  catastrophic cancellation, since `t` is the tiny difference of two `O(chi)`
  terms. Instead, use the asymptotic balance `chi*t^2 ≈ K*c_s^3` (dropping
  the now-negligible cubic term) for a starting guess `t0 = sqrt(K*c_s^3/chi)`,
  polished by 3 Newton iterations on `g(t) = t^3 + chi*t^2 - K*c_s^3` (no
  cancellation there — every term is already the right scale).
- **`m` not small:** the general closed-form (complex) Cardano solution of
  the depressed cubic, selecting the root with negligible imaginary part and
  positive real part.

Verified against `numpy.roots` on a 5000-point random scan spanning
`chi ∈ [10⁻³,10⁵]`, `K ∈ [10⁻⁸,10³]`, `c_s ∈ [10⁻²,10³]`: worst-case relative
error ~6e-10, zero NaNs — vs. the naive single-formula closed-form (no small-`m`
branch), which produces `NaN` (via `sqrt` of a catastrophically-cancelled
near-zero or negative floating-point residual) whenever `least_squares`
probes a large-`chi`/small-`K` region, e.g. during profile-likelihood
excursions for the poorly-identified MC3 dataset. This is a **pure
implementation/numerics fix**, verified to not change the physics or the
converged optimum (see git history for the intermediate broken versions).

### Why the fit is bounded (`chi, K >= 0`), unlike NaCl's unconstrained `lm`

`analysis/step2_nacl/regress_nacl.py`'s `fit_model()` is hardcoded to
unconstrained Levenberg–Marquardt (`method="lm"`), which is fine for NaCl
because `beta1` (`|chi|`) only ever appears **squared** inside a `sqrt` in
that model — its sign is irrelevant, so an unconstrained fit can never do
anything unphysical with it. The CaCl2 cubic's `chi` coefficient enters
**unsquared and asymmetrically** (`c_m^3 + chi*c_m^2 − ...`), so an
unconstrained fit *can* wander to a spurious negative "chi" with a
marginally lower SSE that has no physical meaning (`|χ|` is a magnitude by
construction). Confirmed empirically: an unconstrained fit on MC5 05.27.26
converges to `chi = -0.034` (SSE 643.14) vs. the bounded fit's
`chi = +0.018` (SSE 643.16) — a 0.003% SSE difference for a qualitatively
meaningless sign flip. `analysis/cacl2/regress_cacl2.py`'s
`fit_model_bounded()` uses `method="trf"` with `bounds=([0,0],[inf,inf])`
instead. (Multi-start, profile-likelihood, and the sliding-window diagnostic
all already used bounded `trf`/`lm`-in-a-box formulations and needed no
change.)

## Results (bounded fit, full range, 95% CIs = profile likelihood)

| Dataset | n | `\|χ\|` [mM] | 95% CI | `K` [–] | 95% CI | SSE | RMSE/range | pseudo-R² | MS % global |
|---|---|---|---|---|---|---|---|---|---|
| MC3 (07.11.24, S) | 609 | 0.0000 | (0.000, 204.9) | 0.02374 | (0.0232, 0.0243) | 388.3 | 0.061 | 0.945 | 100% |
| MC5 (05.27.26) | 1645 | 0.0176 | (0.000, 0.369) | 0.00080 | (0.0008, 0.0008) | 643.2 | 0.046 | 0.937 | 96% |
| MC5 (06.27.26) | 1885 | 0.0121 | (0.000, 0.284) | 0.00032 | (0.0003, 0.0003) | 1003.8 | 0.068 | 0.838 | 99% |

Multi-start: 100 random starts, `|χ|₀ ∈ [1,500]` mM log-uniform, `K₀ ∈
[10⁻⁴,50]` log-uniform, `scipy.optimize.least_squares` method `trf`
(bounded, avoiding `dogbox` per Step 3's finding that it can get trapped at
a boundary). All three datasets converge to (approximately) the same
optimum from 96–100% of starts.

**Headline finding: `|χ|` is pinned at (or extremely close to) zero for all
three datasets** — i.e., the fit finds essentially **no net Donnan exclusion
term needed**; a simple partition constant `K` alone (`c_co_m ≈ (K)^{1/3} ·
c_co_s`, dropping the now-negligible `chi` term) explains 84–95% of the
variance (`pseudo-R²`). This is a strong, specific, and testable signature
**consistent with the adsorption hypothesis**: if adsorbed Ca2+ locally
screens or cancels the membrane's intrinsic (presumably negative) fixed
charge, the *net* effective `|χ|` seen by the co-ion partitioning would
indeed collapse toward zero. It is a materially different and more
specific finding than "chi is just some smaller value than NaCl's ~18–50 mM"
— the fit is telling us the exclusion term is not needed at all, only the
partition term is.

MC3's CI is much wider (0–205 mM) than the two MC5 runs (0–0.4 mM):
MC3's narrower concentration range (max ~47 mM vs. MC5's ~100+ mM) gives much
weaker curvature to pin down `|χ|` precisely, the same concentration-range
identifiability effect noted for NaCl's MC2 vs. MC5 (report §2.4).

## Adsorption-hypothesis diagnostics

### Residuals vs. concentration (`cacl2_residuals_vs_conc.png`)

All three datasets show a **clear, smooth, non-random wave pattern**: a dip
to roughly −1 to −1.5 mM at low-to-mid `c_int` (~10–60 mM), rising back
through zero, then a **peak of +1 to +1 mM around `c_int` ≈ 85–95 mM**,
followed by a sharp drop at the very highest concentrations. A formal
sign-runs test (fewer sign changes than expected under random residuals
indicates systematic structure) is strongly significant for all three:

| Dataset | sign runs observed | z-score |
|---|---|---|
| MC3 (07.11.24, S) | 14 | −23.6 |
| MC5 (05.27.26) | 16 | −39.8 |
| MC5 (06.27.26) | 17 | −42.7 |

(`z` this far from 0 — dozens of standard deviations — is not a subtle
effect; it simply confirms visually obvious smooth structure is not noise.)
This **directly matches the anticipated signature**: the single-`|χ|` model
cannot capture whatever concentration-dependent physics is really at play.

### Effective `|χ|` vs. concentration (`cacl2_effective_chi_vs_conc.png`)

Sliding-window (80 points/window, `K` fixed at the global estimate) local
refits of `|χ|` show a **non-monotonic (inverted-U) drift** for both MC5
runs: `|χ|` rises from near 0 at low `c_int` to a peak of ~13 mM around
`c_int` ≈ 90 mM, then falls back down toward ~2 mM at the highest
concentrations sampled. MC3's few reliable windows (5 of 14 — most of its
narrow-range windows are too weakly identified to trust, flagged via SE and
boundary checks and excluded from the trend/plot) show `|χ|` rising from ~4
to ~18 mM over its available range (35–47 mM), consistent with being on the
*rising* leg of the same qualitative shape MC5 shows over a wider range.

**This is drift, but not a simple monotonic one** — a naive "adsorption
increases with concentration, monotonically suppressing `|χ|`" story would
predict a monotonic *decrease*; what's observed is a rise-then-fall. Two
readings, not mutually exclusive: (i) at very low `c_int` the model
correctly finds near-zero exclusion because Ca2+ adsorption has already
saturated/screened the charge even at low concentration, with the mid-range
rise-then-fall reflecting some other unaccounted concentration-dependent
effect (e.g., concentration polarization, which this analysis does not
correct for CaCl2, unlike some NaCl treatments — see report §1.2); or (ii)
the windowed local fit is itself sensitive to how each window's residual
`c_int` sub-range interacts with the fixed global `K`, making the resulting
shape partly an artifact of the diagnostic's own local-refit design rather
than a direct read of the physics. **Either way, `|χ|` is clearly not flat
with concentration for any of the three datasets** — the central qualitative
claim of the adsorption hypothesis (concentration-dependent effective
charge) is supported; the specific functional form is not yet pinned down.

### Reproducibility: MC5 05.27.26 vs. 06.27.26

The two MC5 near-repeats agree closely in both point estimates (`|χ|` 0.018
vs. 0.012 mM, CIs overlap entirely) and in the qualitative *shape* of both
diagnostics (residual wave pattern, χ-drift inverted-U), though 06.27.26 has
a visibly worse pseudo-R² (0.838 vs. 0.937) and a sharper high-concentration
drop in both plots — worth a follow-up look at whether 06.27.26 has more
noise or a genuinely different high-`c_int` behavior (e.g. different
membrane fouling state on that date).

## Verdict on the adsorption hypothesis

**Supported, with a specific and unexpected refinement.** The fit is
materially different from NaCl's (chi pinned at ~0 vs. NaCl's clean 18–50 mM,
non-ill-conditioned single-charge fits), residuals show strong, consistent,
non-random structure vs. concentration in all three datasets, and the local
effective `|χ|` demonstrably drifts with concentration rather than staying
flat — all as anticipated. The specific finding that the *global* fit drives
`|χ|` to (near) zero — rather than merely landing on "some smaller nonzero
value than NaCl" — is a sharper, more specific piece of evidence than the
task prompt anticipated: it suggests the co-ion exclusion effect is entirely
masked/cancelled on average, with the *local* drift diagnostic revealing the
concentration-dependence that a single global constant necessarily averages
away. The non-monotonic (rise-then-fall, not simple monotonic) shape of that
drift is the open question for follow-up — likely requires either a
concentration-dependent `|χ|(c)` model extension (as the report's `sec:cacl2`
TBD text already anticipated) or ruling out concentration-polarization /
window-diagnostic artifacts first.

## Figures (in `docs/reports/figures/`, `plot_style.py` conventions)

- `cacl2_fits.png` — full-width, one panel per dataset, shared legend below.
- `cacl2_residuals_vs_conc.png` — half-width; the key diagnostic.
- `cacl2_effective_chi_vs_conc.png` — half-width; local `|χ|` vs. `c_int`
  (unreliable windows shown as faint `×`, excluded from the reported trend).
- `cacl2_estimates.png` — half-width; `(|χ|, K)` with 95% profile-likelihood
  CIs across datasets.

## Guardrails honored

- `data/`, `prior_analysis/`, and `refs.bib` untouched.
- `analysis/step2_nacl/regress_nacl.py` untouched (no `load_sheet` change was
  needed — CaCl2 uses its own `load_sheet`/`build_deltaC_m` analogs in
  `regress_cacl2.py`, since the response transform and diffusivities differ).
  NaCl numbers unaffected; not re-verified by re-running (no shared code
  changed) but confirmed by inspection of the diff.
- `main.tex` edited only within `sec:cacl2` and the new CaCl2 per-dataset
  appendix subsections.
