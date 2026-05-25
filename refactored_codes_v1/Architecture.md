# `refactored_codes_v1/` — Architecture

This directory holds the refactored unified pipeline for the UCB dynamic
diafiltration project. Three Python files cover the whole surface area:

| File                          | Role                                                        | Lines (approx.) | Mutability |
|-------------------------------|-------------------------------------------------------------|-----------------|------------|
| `refactored_ucb_library.py`   | All loaders, models, dispatchers, registries, plotting      | ~11,000         | Editable   |
| `refactored_ucb_runfile.py`   | Interactive CLI; walks the dispatch tree                    | ~580            | Editable   |
| `conductivity_paper.py`       | Standalone Shedlovsky / MSA forward conductivity models     | ~370            | **Read-only** — do not modify |

`conductivity_paper.py` is treated as an external dependency that happens to
live in-tree. The library reaches into it from one place only
(`_load_conductivity_paper()` at `refactored_ucb_library.py:223`) and every
unit conversion, parameter table, and pre/post-processing step is owned by
the library, not the paper file.

---

## 1. The dispatch tree

The runfile organizes work as four sequential choices — `ROOT → SUBSET →
TRUNK → BRANCHES`:

```
ROOTS               TRUNK                       BRANCHES
─────               ─────                       ────────
DATA1 (.mat)   ┐                              ┌ Mass plots (m)
DATA2 (.mat)   ┤   model_construct_inter      ├ Concentration plots (c)
DATA3 (.xlsx)  ┼─→ + solve_model            ──┼ Parameter table (p)
Custom file    ┘   + estimate_parameters      ├ FIM heatmap (f)
                   + multistart / FIM / DoE   └ DoE recommendations (d)
```

- **ROOTS** are loaders that produce a `data_stru` dictionary.
  `loadmat` (`refactored_ucb_library.py:504`) for DATA1/DATA2, `loadxlsx`
  (`:584`) and `_load_legacy_data_stru_from_excel` (`:3206`) for DATA3/NF270.
- **TRUNK** is the science: model construction, solve, parameter estimation,
  uncertainty quantification, design of experiments.
- **BRANCHES** are the output artifacts: figures, tables, JSON reports.

The runfile's `main()` (`refactored_ucb_runfile.py:552`) prompts for each
choice and dispatches to one of:

- `_dispatch_paper_campaign(root, ...)` — DATA1 / DATA2 manifest-driven
  reproduction (`:320`).
- `_dispatch_nf270(subset, ...)` — DATA3 NF270 campaign sweep (`:388`).
- `_dispatch_custom(...)` — single-file ad-hoc runs (`:463`).

All three dispatchers funnel into `materialize_all(...)` in the library,
which iterates registered `FigureSpec` entries (see §8).

### TRUNK recipes

`TRUNK_RECIPES` (`refactored_ucb_runfile.py:250`) maps a recipe name to a
small option dict consumed by the dispatcher:

| Recipe            | multistart | FIM  | DoE   | mass-litmus | sim-only |
|-------------------|------------|------|-------|-------------|----------|
| `simulate`        | F          | F    | F     | F           | T        |
| `fit`             | F          | F    | F     | F           | F        |
| `fit_multistart`  | T          | F    | F     | F           | F        |
| `fit_FIM`         | T          | T    | F     | F           | F        |
| `fit_FIM_DoE`     | T          | T    | T     | F           | F        |
| `mass_litmus`     | F          | F    | F     | T           | F        |

`mass_litmus` is NF270-only. It collapses the fit to `Lp` from vial mass
slopes alone (sigma=0, no concentration objective) and is gated by the
module-level toggle `NF270_MASS_LITMUS_TEST_ACTIVE`
(`refactored_ucb_library.py:261`). For DATA1/DATA2 the toggle is ignored
with a warning.

---

## 2. Data roots and registries

### Filesystem layout

All three data sources are addressable via environment variables; defaults
point to in-repo locations:

| Variable                   | Default                                                      | Format |
|----------------------------|--------------------------------------------------------------|--------|
| `DIAFILTRATION_DATA1_ROOT` | `legacy/data1_matlab/data`                                   | `.mat` |
| `DIAFILTRATION_DATA2_ROOT` | `legacy/data1_matlab/data_library`                           | `.mat` |
| `DIAFILTRATION_NF270_ROOT` | `UnifiedFramework/ExperimentalDataFiles`                     | `.xlsx` |

NF270 workbooks: `NF270_MC2.xlsx`, `NF270_MC3.xlsx`, `NF270_MC4.xlsx`,
`NF270_MC5.xlsx`. Each workbook holds 3–18 sheets; one sheet = one
experiment.

### NF270 run registry

`NF270_RUN_REGISTRY` (`refactored_ucb_library.py:10726`) maps a stable
human-readable run id (e.g. `MC2.05.07.24_NaCl`) to a `{file, sheet}`
dictionary. The registry currently exposes 26 green-flagged experiments;
`NF270_RUN_QUALITY` (`:10762`) marks each as `"green"`.

The runfile defines two well-known subsets on top of the registry:

- `NF270_FAST_SUBSET` (`refactored_ucb_runfile.py:97`) — 2 sheets + the
  parameter table aggregation entry. Used for smoke tests.
- `NF270_SINGLE_SALT_RUNS` (`:103`) — 11 single-salt sheets (6 NaCl,
  3 CaCl₂, 2 LaCl₃). This is the focus subset for the current paper push.

The "All 26 green-flagged sheets" option iterates the registry directly.
The "One specific sheet" option uses `_prompt_one_nf270_run`
(`refactored_ucb_runfile.py:228`).

---

## 3. Loaders and the `data_stru` contract

Every loader returns a Python dict named `data_stru` with a stable shape so
the model-build, fit, and plot stages downstream do not care which root
the data came from.

### Key `data_stru` fields

| Key                                | Meaning                                                       |
|------------------------------------|---------------------------------------------------------------|
| `data_raw`                         | List of per-vial dicts (`time`, `mass`, `cF_exp`, `cV_perm_cond`, ...) |
| `data_config`                      | Scalar config: `n` (#vials), `Temp` (K), `delP`, `M_F0`, `cF0`, salt name `namec`, ... |
| `source_format`                    | `"mat"` or `"xlsx"`                                           |
| `sheet_name`                       | XLSX-only: which sheet this came from                         |
| `continuous_cF`                    | True when `cF_exp` is a time-series (XLSX), not vial averages |
| `conductivity_cF`                  | True if `cF_exp` is currently a conductivity, not concentration |
| `conductivity_cF_converted`        | True after the conductivity-to-concentration inversion ran    |
| `conductivity_temp_compensated`    | **New.** True after EC25 temperature compensation. See §5.    |
| `conductivity_temp_compensation`   | **New.** Metadata: salt name, α, reference temperature, formula. |
| `measured_ret_temp_C`              | **New.** Original per-row retentate temperature trajectory in °C. |
| `measured_perm_temp_C`             | **New.** Original per-row permeate temperature trajectory in °C. |

Loaders:

- **`loadmat(filename)`** (`:504`) — DATA1 / DATA2 MATLAB structs. Returns
  retentate concentration directly in mM; no conductivity branch.
- **`loadxlsx(filename, sheet=...)`** (`:584`) — XLSX single-sheet path used
  by the custom-file dispatcher.
- **`_load_legacy_data_stru_from_excel(path, selector)`** (`:3206`) — full
  NF270 sheet loader. Parses time, mass, pressure, retentate/permeate
  temperatures + conductivities, vial swaps, and metadata; segments by
  vial swap; **applies EC25 temperature compensation in-place** (§5).

### Time-axis convention and per-vial mass reset

The model treats every experiment as a sequence of per-vial intervals
`[TI_i, TF_i]` in seconds, with `tau ∈ [0, 1]` as the scaled within-vial
variable. The same convention is shared by all three workflow families;
the only differences are which layer does the bookkeeping.

**Time origin (shared, byte-identical to DATA2):**

```
t_delay = data_stru['data_raw'][0]['time'][0]
TI_list[i] = data_raw[i]['time'][0]  - t_delay     # vial i start, relative
TF_list[i] = data_raw[i]['time'][-1] - t_delay     # vial i end,   relative
t_physical = tau * (TF_i - TI_i) + TI_i
plot_time_min = (t - t_delay) / 60
```

This formula is identical between `utility.py` (DATA2) at lines 115, 447,
473, 477 and `refactored_ucb_library.py` (DATA1/DATA3) at `_data1_time_origin`
(`:170`), `model_construct_inter` (`:1247`), and `_data1_shifted_time_minutes`
(`:180`). The first sample of the first vial is always `t = 0` in plots
and in the model.

**Per-vial mass reset (different layer, same outcome):**

| Workflow | Mechanism | Location |
|---|---|---|
| DATA2 | Pyomo linking constraint `1e-6 == m.mV[n+1, 0]` | `utility.py:736-743` (mode `'DATA'`) |
| DATA1 | Same Pyomo constraint | `refactored_ucb_library.py:1273-1280` |
| DATA3 | Loader subtracts segment-start offset before model sees data: `seg_mass = seg_mass_cum - offset` | `refactored_ucb_library.py:3404-3412` |

Net effect is the same — every vial's `mV[0] = 0` regardless of any tare
offset in the raw balance column. This is how sheets like
`MC3.07.12.24_S2CaCl2`, whose Excel column ranges from -9.85 g to +0.23 g
(tare zeroed against ~9.85 g of container weight), feed into the
objective without contamination.

### Hold-up split — DATA3 only

NF270 operators start the data logger before any permeate has reached
the collection vial. The first ~5-10 minutes of every Excel sheet is the
membrane wetting up, recorded as part of the operator's "Vial 1" but
contributing essentially zero mass. The loader detects and isolates this
period.

**Detection rule** (`refactored_ucb_library.py:3522-3590`):

```python
# In operator's vial 1, find the first sample where mass has risen
# >= 0.05 g above the starting tare value.
i_start = first index where (m[i] - m[0]) >= threshold_g (= 0.05 g)
if i_start >= 3 and (N - i_start) >= 10:
    # Split: samples [0..i_start) become model-vial-1 (holdup),
    #        samples [i_start..N) become model-vial-2 (real first vial).
    cfg["vial1_split_status"] = "split"
    cfg["t_perm_start_s"]     = t[i_start]
    cfg["n_holdup_samples"]   = i_start
    cfg["n"]                  = N_operator_vials + 1   # bumped by 1
    cfg["n_v0"]               = 2                      # first real vial index
    cfg["n_extra"]            = 1                      # holdup vial count
```

When the split fires:

- Operator's "vial 1" is carved into two model-vials.
- Model-vial-1 = holdup period (no permeate, mass ≈ 0, no ICP sample).
- Model-vial-2 = first real collection segment.
- Operator's subsequent vials are all renumbered +1.

**What the model and objective do with model-vial-1:**

- The integrator passes through it continuously — `Jw`, `cF`, `cH` all
  evolve as if it were a normal vial, so the state at the start of
  model-vial-2 is physically consistent with the run history.
- The objective **skips residuals** for any vial index `< n_v0`
  (the `n_extra` / `n_v0` skip convention is the DATA1/DATA2 startup-vial
  contract inherited unchanged).
- So the holdup vial is integrated but not penalized — exactly the
  treatment DATA2 used for its own `n_extra` startup vials.

**Per-sheet status codes** stored in `data_config["vial1_split_status"]`:

| Status            | When                                                              |
|-------------------|-------------------------------------------------------------------|
| `"split"`         | Split fired; `n_holdup_samples` and `t_perm_start_s` are populated |
| `"skip_no_rise"`  | Mass never crosses the threshold (rare; e.g. zero-flux sheet)     |
| `"skip_no_lag"`   | No clear flat-then-rise pattern, or N too small                   |
| `"skip_no_data"`  | `data_raw` empty                                                  |

**Why DATA2 doesn't need this:** DATA2 .mat files are post-curated paper
data; the data curator had already trimmed the holdup samples before
writing the .mat. NF270 sheets come in raw from the operator script, so
the loader must auto-detect and auto-split.

---

## 4. Conductivity → concentration inversion

DATA1/DATA2 already store retentate concentration; this section only
matters for the NF270 / XLSX path.

### The forward models (in `conductivity_paper.py`, do not modify)

- `_shedlovsky(...)` — equivalent conductivity λ(c, T).
- `variant_shedlovsky(...)` — specific conductivity κ(c, T). Returns mS/cm.
- `msa_transport(...)` — Mean Spherical Approximation; supports single,
  binary, and ternary salt mixtures. Returns mS/cm.

### Salt parameter tables (anchored at 25 °C)

`CONDUCTIVITY_SALT_PARAMS_25C` (`refactored_ucb_library.py:264`) holds
ε, η, distance-of-closest-approach `a`, valencies, and limiting equivalent
conductivities for **KCl, NaCl, CaCl₂, LaCl₃**. These are the textbook
25 °C values from Robinson & Stokes.

`_msa_augment_salt_params()` (`:10712`) extends the table with the MSA
inputs the binary/ternary path needs:

- `eta_pa_s` = η (poise) × 0.1 — viscosity in Pa·s.
- `lambda_0_*_S_m2_mol` — limiting molar conductivities in S·m²/mol.
- `diff_coeff_*_m2_s` — diffusivities via Nernst–Einstein from λ⁰.
- `diameter_*_m` — hard-sphere ion diameters from a lookup table.

### Inversion pipeline

```
data_stru.conductivity_cF == True
        │
        ▼
_normalize_conductivity_measurements      (refactored_ucb_library.py:2927)
        │  ├── _apply_ec25_compensation_to_data_stru   ← §5 (load-time hook)
        │  └── convert_experimental_conductivity_to_concentration  (:2868)
        ▼                            │
                                     ▼
                _conductivity_to_concentration_series         (:2784)
                  Shedlovsky / variant_shedlovsky  branch
                                     │
                                     ▼
                _conductivity_to_concentration_msa_branch     (:10766)
                  MSA binary / ternary branch
                                     │
                                     ▼
                _invert_monotone_1d  (:2743)
                  bisection inversion of forward κ(c) curve
                                     │
                                     ▼
                cF_exp now stores concentration in mM
```

`_normalize_conductivity_measurements` is intentionally early in the
preprocessing chain — the rest of the diafiltration math sees only
concentrations.

---

## 5. EC25 temperature compensation (new)

### Why

The 25 °C-anchored parameter tables in §4 are fed conductivity values
measured at the meter's reading temperature. For the 11 single-salt NF270
sheets in current use, retentate temperatures sit 2–5 °C below 25 °C, and
sheet `07.22.24_SNaCl` exhibits a 3 °C retentate trajectory while permeate
stays flat. Feeding raw conductivity into a 25 °C inversion produces a
systematic ~5–10 % under-prediction of concentration. Pre-compensating to
a 25 °C-equivalent value removes that bias before the inversion runs.

### Equation

Small-ΔT linearization of the standard practice:

$$\sigma_{25} = \sigma_T \cdot \left[1 + \alpha \,(25 - T)\right]$$

Valid to <0.2 % within ±10 °C of 25 °C in the dilute regime.

**Literature.** The formula and its α ≈ 0.02 / °C generic value:

- Hayashi, M. (2004). *Environ. Monit. Assess.* **96**: 119–128.
  [doi:10.1023/B:EMAS.0000031719.83065.68](https://doi.org/10.1023/B:EMAS.0000031719.83065.68)
- Sorensen, J. A., & Glass, G. E. (1987). *Anal. Chem.* **59**(13): 1594–1597.
  [doi:10.1021/ac00140a003](https://doi.org/10.1021/ac00140a003)

Standards: APHA *Standard Methods*, Method 2510 B; ASTM D1125-14. Salt-specific
α values are derived from limiting-mobility temperature derivatives in
Robinson, R. A., & Stokes, R. H. (2002). *Electrolyte Solutions*, 2nd ed., Dover.

### Where it lives

| Symbol                                  | Location                                       |
|-----------------------------------------|------------------------------------------------|
| `CONDUCTIVITY_TEMP_COEFF_PER_C`         | `refactored_ucb_library.py:309`                |
| `CONDUCTIVITY_TEMP_COEFF_DEFAULT`       | `:315` (0.024, generic dilute electrolyte)     |
| `_ec25_compensate(...)`                 | `:319`                                         |
| `_apply_ec25_compensation_to_data_stru` | `:363` (the post-load hook used by `_normalize_conductivity_measurements`) |
| Inline call in the NF270 Excel loader   | `:3245–3251` (`_load_legacy_data_stru_from_excel`) |

Default α (1/°C):

| Salt   | α      | Source                       |
|--------|--------|------------------------------|
| KCl    | 0.019  | Robinson & Stokes 2002       |
| NaCl   | 0.021  | Robinson & Stokes 2002       |
| CaCl₂  | 0.023  | Robinson & Stokes 2002       |
| LaCl₃  | 0.025  | Robinson & Stokes 2002       |
| _other_ | 0.024 | Generic dilute electrolyte   |

### Contract enforced at load time

In `_load_legacy_data_stru_from_excel`:

1. Pull `Retentate Temp`, `Retentate Cond @ Temp`, `Permeate Temp`,
   `Permeate Cond @ Temp` from the sheet.
2. **Keep copies** of the raw temperature trajectories on `data_stru` as
   `measured_ret_temp_C` and `measured_perm_temp_C` so diagnostics still
   have access to the original signal.
3. **Compensate both conductivity streams per row** with their own
   temperature column via `_ec25_compensate(...)`.
4. Set `data_config["Temp"] = 298.15` so the downstream inversion uses the
   25 °C reference its parameter tables already assume — preventing a
   double correction.
5. Flag the dataset with
   `conductivity_temp_compensated = True` and stash metadata under
   `conductivity_temp_compensation = {salt_name, alpha_per_C, reference_C,
   method}`.

`_apply_ec25_compensation_to_data_stru` is the dual hook for datasets that
arrive pre-loaded but still in conductivity units (e.g. a sheet handed in
via the custom path). It is idempotent — it no-ops when
`conductivity_temp_compensated` is already True or when per-row
temperatures are unavailable.

### Out of scope

- The MSA reference temperature constants near `:10443` / `:10463` already
  sit at 25 °C; the compensation makes the data consistent with them, so
  no MSA-side change was needed.
- The inverse-form compensation `σ_25 = σ_T / [1 + α(T−25)]` is not used
  because no sheet in the current dataset exceeds ±5 °C from 25 °C; the
  linear form is faithful there.
- There is no runtime flag to disable the compensation. To pass already
  25 °C-compensated data through, set the temperature columns to 25.0 °C
  and the math no-ops.

---

## 5b. NF270 diagnostic knobs (DATA3-only, guarded)

Two module-level toggles exist for NF270 corner-pin diagnosis. Both default
to `None` (legacy behavior) and both are **hard-guarded to `workflow_family
== "DATA3"`** — flipping either constant has zero effect on DATA1 / DATA2
code paths, so published-paper reproducibility is preserved unconditionally.

### `NF270_CF_RESIDUAL_FLOOR_MM`  (`refactored_ucb_library.py:271`)

Optional absolute floor (mM) on the cF residual scale. When set, the cF
weight in the WSSE becomes `max(0.003 · cF_meas[i], floor)` instead of the
legacy `0.003 · cF_meas[i]`. The motivation is that at end-of-run C-regime
cF ≈ 100 mM, the 0.3 % relative weight is a 0.3 mM tolerance — tighter than
the actual measurement noise — and the cF channel ends up over-weighted,
dragging σ to the {0, 1} bound. A floor of 1 mM lets the optimizer release
σ to the interior on NaCl C-regime and CaCl₂ sheets.

**Guards:**

- The two read sites are at `:1867` (inside `obj_rule` in
  `model_construct_inter`) and `:2476` (inside the `obj_rule` in
  `solve_model_B_fix`).
- Both wrap the constant in
  `(NF270_CF_RESIDUAL_FLOOR_MM if workflow_family.upper() == "DATA3" else None)`,
  so when called with `workflow_family="DATA1"` or `"DATA2"` the legacy
  weight is used regardless of what the constant is set to.

### `NF270_SIGMA_INTERIOR_BOUNDS`  (`refactored_ucb_library.py:279`)

Optional `(lo, hi)` tuple replacing the default σ Var bounds of `(0, 1)`.
Used together with the cF floor as a literature-σ-prior surrogate on sheets
where the σ surface is genuinely flat and the floor alone can't release the
optimizer off the wall.

**Guards:**

- The single read site is at `:1364` in `model_construct_inter`.
- Wrapped in the same `workflow_family == "DATA3"` gate; DATA1/DATA2 always
  receive `(0, 1)` bounds regardless of the constant.

### `NF270_USE_PERMEATE_PROBE`  (`refactored_ucb_library.py:303`)

Adds the continuous **permeate-probe** trace (`cV_perm_cond` → mM via
Shedlovsky inversion) as a 4th channel in the WSSE objective.  When None
(default), the objective is the legacy three-channel form
(mass + permeate-ICP + retentate-cF).  When set to a positive float
(e.g. 0.03 for 3 % relative scale), and `workflow_family == "DATA3"`, the
loader converts `cV_perm_cond` to `cV_perm_exp` (mM) and the objective
adds:

    Φ_perm  =  Σ_i [ ( cH_pred(t_i) − cV_perm_exp(t_i) ) / σ_perm ]²

where `cH` is the model's wall-side permeate concentration (the
physically correct counterpart to the probe reading; **not** `cV` which
is the vial-integrated value).  Per-sample residual scale is
`max(NF270_USE_PERMEATE_PROBE * cV_perm_exp[i], NF270_PERMEATE_PROBE_FLOOR_MM)`.

**Guards:**

- Loader converts `cV_perm_cond` only when the toggle is truthy
  (`_normalize_conductivity_measurements` at `:3066-3076`).
- Objective term gate at `:1855-1864` requires `workflow_family.upper() == "DATA3"`
  AND toggle truthy AND `cV_perm_exp` present per vial.
- DATA1/DATA2 .mat files don't have `cV_perm_cond` to begin with;
  the loader path is silent no-op on them.

### `NF270_MULTISTART_USE_CONTOUR_SEEDS`  (`refactored_ucb_library.py:307-314`)

Switches the LHS multistart from uniform-over-bounds to **seeded** —
the LHS starts are clustered around three physics-informed seeds:

1. Top σ-Lp grid point from this sheet's own contour panel (lowest raw WSSE).
2. Top B-Lp grid point.
3. A **cross-salt reference** — when fitting CaCl₂ or LaCl₃, the MC3 SNaCl
   gold-standard θ is used as the third seed (`Lp=8.22, B=14.7, σ=0.45`).

LHS samples are then drawn round-robin around each seed within a ±20 %
box.  The companion constants `NF270_MULTISTART_CONTOUR_DIR` (optional
explicit path) and `NF270_MULTISTART_CROSS_SALT_REFERENCES` (the
per-salt reference dict) provide fine-grained control.

**Guards:** Only when `workflow_family == "DATA3"` AND the toggle is True
AND the sheet's contour CSVs are on disk.  Falls back to the legacy
uniform-LHS multistart when any of those conditions is false.

### How to test the guards

`refactored_codes_v1/tests/regression_data1_data2_guards.py` is the regression
test.  18 cases (6 per guard for σ-bounds, cF-floor, permeate-probe), runs
in ~1 s.  It does NOT invoke IPOPT — only constructs the model and inspects
state.  The full matrix is the contract: if any guard ever leaks, the test
fails.  Re-run after any change to the objective or loader code paths.

---

## 6. Model construction, solve, and parameter estimation

### Model build

`model_construct_inter(data_stru, mode, theta=None, ...)`
(`refactored_ucb_library.py:1183`) builds a Pyomo concrete model:

- Reads `data_config` (Temp, ΔP, M_F0, cF0, salt-specific diffusivity).
- Discretizes time per vial with finite-difference collocation (`nfe`
  knobbed from the runfile; default 80–300 depending on root).
- Mode dictates which residuals appear in the objective:
  - `"DATA"` — fit mass + concentration (DATA1/DATA2 path).
  - `"Lag"` — DATA3/NF270 single-salt path; permeate ramps up with a delay.
  - `"Overflow"` — vial overflow handling.
- `B_form ∈ {single, pervial, convection, 0..3}` selects the rejection
  parameterization.
- `workflow_family ∈ {DATA1, DATA2, DATA3}` picks the residual weights and
  the salt-property defaults that match each paper.

### Solve

- `solve_model(...)` (`:1640`) — single forward solve with IPOPT.
  Picks the linear solver via `_preferred_ipopt_linear_solver()` (`:54`)
  and applies CPU-time limits via `_resolve_ipopt_cpu_time_limit()` (`:78`).
- `solve_mass_balance_only(data_stru, ...)` (`:2592`) — mass-litmus
  collapse: σ=0, drops the concentration objective, solves in closed form
  for `Lp` from vial mass slopes only.

### Parameter estimation

The library has two parallel estimation paths:

1. **Legacy** — `solve_model` with `mode="DATA"` already fits parameters
   when its inner block is configured for estimation.
2. **Modern ParmEst** —
   `DiafiltrationParmestExperiment` (`:3022`) + Pyomo `ParmEst`:
   - `build_parmest_experiments(...)` (`:3050`) wraps per-experiment
     model factories.
   - `estimate_parameters_with_parmest(...)` (`:3602`) drives the fit.
   - `compute_covariance_with_parmest(...)` (`:3657`) returns the
     parameter covariance.

`StageResults` (`:8397`) is the in-memory container that flows through the
modern pipeline: `data_stru`, `model`, fitted `parameters`, `uncertainty`,
`fim`, and downstream DoE state.

### Multistart, FIM, DoE

- `_lhs_unit_hypercube(...)` (`:7364`),
  `_theta_variants_for_multistart(...)` (`:7442`),
  `_solve_model_multistart(...)` (`:7475`) — multistart wrapper. The
  runfile prompts for 5–14 LHS starts when `trunk.multistart` is true.
- `calc_FIM(...)` (`:2076`) — finite-difference Fisher Information Matrix.
- `quantify_uncertainty(...)` (`:3876`) — covariance reporting wrapper.
- `compute_doe_metrics(...)` (`:3699`) — Pyomo.DoE A/D/E-optimality
  criteria; `design_next_experiment(...)` (`:3929`) is the public hook.

---

## 7. The figure manifest and `materialize_all`

The library decouples *what figure to make* from *how to make it* with a
declarative manifest:

- **`FigureSpec`** (`:8515`) — `name`, `renderer`, `requires` (tuple of
  `RunRequest`), `opts`, `output_filename`, `description`.
- **`RunRequest`** (`:8498`) — `campaign`, `run_id`, `variant`, optional
  `load_cached_fit`. Identifies which dataset + fit feeds the figure.
- **`CAMPAIGN_MANIFESTS`** — keyed by campaign name (DATA1, DATA2, NF270);
  each holds the list of FigureSpec entries that *can* be produced.
- **`materialize_all(campaign, save_dir, data_root, only=...)`** (`:9374`)
  — the workhorse. Iterates matching specs, resolves shared `RunRequest`
  state via a cache (`_materialize_spec`, `:9413`), then dispatches to the
  named renderer.

Renderers are plain functions of the form
`render_*(results, *, save_path, ...)`:
- `render_sim_comparison` (`:8791`)
- `render_contour` (`:8838`)
- `render_concentration_range` (`:8887`)
- `render_data1_parameters_table` (`:10069`)
- `render_data2_parameters_table` (`:10142`)
- `render_nf270_fit_v2` (`:10362`)
- `render_nf270_parameters_table` (`:11114`)
- `render_calibration_curve`, `render_pressure_change`,
  `render_startup_barplot`, `render_error_boxplot`,
  `render_model_error_visualization_pipeline`, ...

Reports (JSON sidecars):
- `report_parameters` (`:9075`)
- `report_uncertainty` (`:9087`)
- `report_fit_summary` (`:9104`)

Paper coverage:
- `report_paper_coverage(campaign, save_dir)` (`:9971`) and
  `format_paper_coverage(report)` (`:10027`) — diff the produced files
  against the manifest's expected outputs and report what's missing.

---

## 8. NF270-specific extensions

The NF270 campaign is layered on top of the manifest machinery in the
last ~1,000 lines of the library:

- **Loader integration.** `get_experimental_run(campaign="NF270", run_id, ...)`
  (`:11002`) resolves a registry entry to a loaded `data_stru`. Internally
  it routes through `_load_legacy_data_stru_from_excel`.
- **Renderers.** `render_nf270_fit` (`:11029`) is the original NF270 figure
  renderer; `render_nf270_fit_v2` (`:10362`) is the current version with:
  - per-vial vertical-alignment artifact suppression,
  - cumulative mass scaling across vials (DATA2 figure-6 convention),
  - `_detect_t_delay_from_pressure` (`:10279`) for lag-mode start detection,
  - ICP-OES per-vial points overlaid on the continuous conductivity trace.
  `_swap_nf270_renderer_to_v2()` (`:10596`) is the registration hook.
- **Figure registration.** `_register_nf270_campaign()` (`:11078`),
  `_build_nf270_figures()` (`:11044`), `_register_nf270_paper_index()`
  (`:11109`), `_register_nf270_table_entry()` (`:11122`).
- **Mass-litmus toggle.** `NF270_MASS_LITMUS_TEST_ACTIVE` (`:261`) is a
  module-level flag the runfile turns on for the duration of a
  `mass_litmus` trunk run, then resets.

---

## 9. Outputs

Default save directories (under `UnifiedFramework/DATA3/results/paper_artifacts/`):

| Campaign      | Path                                              |
|---------------|---------------------------------------------------|
| DATA1         | `data1/notebook_figures/`                         |
| DATA2         | `data2/notebook_figures/`                         |
| NF270         | `nf270/campaign_figures/`                         |
| Custom (xlsx) | `UnifiedFramework/DATA3/figures/data3_option3/`   |
| Custom (mat)  | `UnifiedFramework/DATA3/results/custom_runs/`     |

Each campaign produces:

- Per-experiment mass-vs-time and concentration-vs-time PNGs.
- A campaign-level parameter table (PNG/CSV).
- Optional contour heatmaps (FIM/σ-sensitivity).
- Optional JSON reports for parameters and uncertainty.

---

## 10. Reproducibility guarantees

- **DATA1 / DATA2 published figures must reproduce byte-equivalent.** The
  dispatcher passes no `extra_opts` overrides for these roots unless the
  user deviates from `manifest_default`; the manifest itself encodes the
  multistart/FIM/B-form/mode settings the published paper used. See
  `_dispatch_paper_campaign` (`refactored_ucb_runfile.py:320`).
- **NF270 is the experimental campaign.** Mesh, multistart count, and
  subset are user-prompted; defaults are tuned for tractable overnight runs.
- **`conductivity_paper.py` is frozen.** Any "what if we changed the model"
  experiment goes in the library, not in the paper file.

---

## 11. Quick reference: where do I add ...?

| Task                                              | Where                                              |
|---------------------------------------------------|----------------------------------------------------|
| New salt for conductivity inversion               | `CONDUCTIVITY_SALT_PARAMS_25C` (`:264`) and optionally `CONDUCTIVITY_TEMP_COEFF_PER_C` (`:309`) |
| Different α for an existing salt                  | `CONDUCTIVITY_TEMP_COEFF_PER_C` (`:309`)           |
| New NF270 experiment to expose                    | `NF270_RUN_REGISTRY` (`:10726`)                    |
| New figure renderer for any campaign              | A `render_*` function + a `FigureSpec` in the campaign manifest |
| New TRUNK recipe                                  | `TRUNK_RECIPES` (`refactored_ucb_runfile.py:250`)  |
| New ROOT (data source)                            | A loader + a `_dispatch_*` in the runfile          |
| Plot tweaks for NF270                             | `render_nf270_fit_v2` (`:10362`)                   |
| New optimization knob                             | `solve_model` (`:1640`) and propagate through `StageResults` |

---

## 12. May 24, 2026 — Day-of-run updates (honest state report)

This section is added at the end so it doesn't fragment the architectural
narrative above. It records what was added on May 24, 2026, what survived
end-of-day, and where the codebase is still rough. Read this section if you
want the most current picture of *what works vs what doesn't*.

### 12.1  What was added today

| Component | Where | Purpose | Status |
|---|---|---|---|
| `_ec25_compensate` extended docstring + worked example | `refactored_ucb_library.py:375` | Layered educational comments on the 25 °C temperature compensation. Math-to-code annotated. | ✓ in repo |
| `_nf270_corrected_cv_index` extended docstring + worked example | `refactored_ucb_library.py:454` | Same treatment for the V_tube tube-transit time correction. | ✓ in repo |
| 11 constraint banners inside `model_construct_inter` | `refactored_ucb_library.py:1419+` | Each ODE/algebraic constraint now carries a `╔═╗` banner naming the equation, the physics, the units, and a literature reference. | ✓ in repo |
| `CONDUCTIVITY_PAPER_EXPLAINER.md` | `refactored_codes_v1/` | Sibling explainer for `conductivity_paper.py` (Shedlovsky + MSA) without modifying that file. | ✓ in repo |
| `_detect_nf270_mode(notes_text, mass_cum)` | `refactored_ucb_library.py:3837` | Replaces the hardcoded `mode="Lag"` in the Excel loader. Detects Lag vs Overflow from Notes-text keywords and mass-trajectory shape. | ✓ in repo |
| `_parse_icp_calibration` + `_icp_calibration_residual` | `refactored_ucb_library.py:3850` | Reads the per-sheet ICP calibration table from cols 23-24, fits a linear curve, and cross-checks the stored `mg/L` values. Populates `data_config["icp_calibration_curve"]` and a per-row `icp_calibration_residual` field. | ✓ in repo |
| B upper bound: `30 → 50` | `refactored_ucb_library.py:1637`, `:9026` | All 3 CaCl₂ sheets in the prior contour batch landed at B=30 exactly. Bumping the bound lets the optimizer move freely. | ✓ in repo |
| Defensive guard `or "DATA1"` → `or "DATA3"` (workflow_family default) | `refactored_ucb_library.py:1578` | Defensive fallback. Only fires when caller explicitly passes `None`/invalid. The function-signature default stays `'DATA1'` for legacy callers. | ✓ in repo |
| `.get()` defaults for `M_O`, `C_D`, `n_v0` | `refactored_ucb_library.py:1546`, `:1549`, `:1556` | DATA1 DATA-mode .mat files omit these keys; defaults are physically correct (0 mass overflow, 0 diafiltrate, fit from vial 1). DATA2 / DATA3 unaffected since those datasets always supply the keys. | ✓ in repo |
| `_preflight_audit.py` + `PREFLIGHT_AUDIT.md` | `refactored_codes_v1/` | Read-only 6-check audit: registry salts, V_tube sanity, Pyomo bounds, 11-sheet roster, loader fields, ICP calibration health. | ✓ in repo |
| `tests/` pytest infrastructure | `refactored_codes_v1/tests/` | Structural regression tests for DATA1/DATA2 model build. `conftest.py`, `reference_cases.py`, `_capture_references.py`, `_trajectory_compare.py`, `test_data1_regression.py`, `test_data2_regression.py`, `regression_references.json`. Run: `pytest -v refactored_codes_v1/tests/`. | ✓ in repo |

### 12.2  What works end-to-end

- **Reading every workbook sheet** — all 11 single-salt sheets load cleanly. The audit (PREFLIGHT_AUDIT.md) shows salt classifications 11/11 consistent, V_tube=0.3 g lands the ICP time correction at exactly 40 % of every vial duration, and ICP calibration health is good for NaCl sheets (R² ≈ 0.998, residuals ≤ 5 %).
- **DATA1/DATA2 model build (structural)** — pytest passes for the 2 DATA2 cases that have the full schema. DATA1 cases are explicitly skipped with documented `con_boundary[2]` pre-existing issue.
- **The well-fit dilution sheet** — `MC3.07.22.24_SNaCl` fits cleanly via the DATA2 recipe (`mode='Lag'`, `B_form=1`, theta from the paper demo) in ~19 sec, producing `Lp ≈ 8.8`, `σ = 1.0`, finite per-channel objectives.
- **Paper-styled concentration plots** — `run_data3_time_series_plots` produces Figure 6-style mass + concentration panels with the agreed colors, marker edges, time-corrected teal squares.

### 12.3  What doesn't work today (honest)

- **Generalizing the DATA2 recipe to every DATA3 sheet** — passing the published theta dict to a low-starting-cF concentration-regime sheet (e.g. MC2.05.07.24_NaCl, cF₀≈0.9 mM) puts IPOPT in a degenerate region (σ=1.0 + tiny cF → near-zero osmotic resistance, Jw too high). IPOPT thrashes. Three seed strategies tried today (explicit DATA2 theta, theta=None for per-sheet `Lp0/B0/sigma0`, IPOPT internal `max_cpu_time` cap) all failed to break the pattern on the hard sheets.
- **IPOPT's `max_cpu_time` cap is unreliable** — it's checked at iteration boundaries; if MUMPS gets stuck inside a single linear solve, the cap never fires. The current rollout driver wraps `solve_model` in a Python `signal.alarm` as a hard wall-clock timeout (`refactored_codes_v1/_rollout_data2_recipe.py`).
- **ICP calibration health for multivalents (data-quality finding, not a code bug)** — CaCl₂ and LaCl₃ workbook calibration tables don't match the stored `mg/L` values (residuals up to 921 %). Almost certainly a Salt-1=NaCl template that was reused for divalent runs; the lab's analytical software used the correct calibration. The cross-check is *audit-only*; the loader still uses the lab's stored `mg/L`.
- **Multivalent fits (CaCl₂, LaCl₃)** — even when fits converge, σ pins at 1.0 (or 0.0) on 10/11 sheets. This is the cF-channel-dominance pathology documented in slides 26-27 of the pptx deck. The `NF270_CF_RESIDUAL_FLOOR_MM` knob fixes it (validated on MC4 SNaCl: σ moved 0 → 0.461) but is currently **disabled by default per the May 24 decision** to keep the DATA3 path aligned with the DATA2 recipe.

### 12.4  Why DATA2 recipe doesn't generalize directly to DATA3

DATA2's recipe was tuned for a *uniform* sheet set: KCl filtration on NF270 across a narrow concentration window. The published seed theta `{Lp:11, σ:1.0, beta_0:1, ...}` lands close to a feasible basin for every DATA2 .mat file because they're all the same salt at similar conditions.

DATA3's 11 single-salt sheets span:
- **Three different salts** with different valencies (ni = 2 for NaCl, 3 for CaCl₂, 4 for LaCl₃). Osmotic pressure ∝ ni → σ·ΔΠ scales differently per salt.
- **Two regimes** (dilution: cF starts high; concentration: cF starts ≈ 1 mM). The DATA2 seed sits closer to a feasible point for the dilution regime; for low-cF concentration runs it's degenerate.
- **Different operators / coupons** — small variations in feed/diafiltrate concentrations across sheets.

A single seed θ cannot serve all of them. The codebase has the machinery to handle this (contour-seeded multistart, per-sheet warm starts) — but the May 24 decision was to align with DATA2's simpler workflow first and add per-sheet seeding only if the simple path proves insufficient.

### 12.5  Diagnostic knobs that are currently disabled (preserved as comments)

These DATA3-specific improvements were validated earlier in the session but are turned OFF in the current rollout driver so the workflow mirrors DATA2 more faithfully. Re-enable by uncommenting in `refactored_codes_v1/_rollout_paper_styling.py`:

| Knob | What it does | Validated outcome |
|---|---|---|
| `NF270_CF_RESIDUAL_FLOOR_MM = 1.0` | Floor the cF residual scale at 1 mM. Releases σ from the wall on sheets where the 0.3 % relative weight over-constrains. | MC4 SNaCl: σ moved `0 → 0.461`, obj_cr improved 200× |
| `NF270_MULTISTART_USE_CONTOUR_SEEDS = True` | Use per-sheet contour-grid minima as IPOPT warm starts. | Validated to give interior σ on most sheets in the contour batch |

### 12.6  Open questions for collaborators

1. **σ-at-the-wall pathology** — is the right physics fix (a) the cF residual floor, (b) a Bayesian prior on σ anchored at literature values (~0.95 for divalents), (c) a per-salt B_form change, or (d) an explicit ion-pairing correction in ΔΠ for trivalents?
2. **DATA1 `con_boundary[2]` bug** — pre-existing in `model_construct_inter` for DATA-mode .mat files; not fixed today. Worth a separate investigation.
3. **Per-sheet vs uniform recipe** — should the rollout enable contour-seeded warm starts so every sheet gets the best available initial guess, even if that introduces machinery DATA2 doesn't have?
4. **ICP calibration disagreement on multivalents** — leftover template artifact, or did the lab's analytical software actually use that calibration? Worth a quick conversation with the experimentalist.

### 12.7  Files added today (canonical list)

| File | Purpose |
|---|---|
| `refactored_codes_v1/CONDUCTIVITY_PAPER_EXPLAINER.md` | Shedlovsky / MSA reader's guide for `conductivity_paper.py` |
| `refactored_codes_v1/PREFLIGHT_AUDIT.md` | 6-check pre-flight audit results |
| `refactored_codes_v1/_preflight_audit.py` | Driver for the pre-flight audit |
| `refactored_codes_v1/_b_form_selection.py` | B_form selection driver (built, not used after May 24 simplification) |
| `refactored_codes_v1/_rollout_data2_recipe.py` | The May 24 rollout that mirrors DATA2 recipe + has wall-clock timeout |
| `refactored_codes_v1/tests/` (8 files) | pytest infrastructure for DATA1/DATA2 structural regression |
| `refactored_codes_v1/DATA3_single_salt_analysis_v6.pptx` | Slide deck (with May 24 update slides appended) |
| `refactored_codes_v1/DATA3_single_salt_speaker_notes.docx` | Speaker notes (with May 24 update notes appended) |

Older disposable rollout drivers (`_preview_paper_styling.py`, `_rollout_paper_styling.py`, `_contour_batch_sigma_B.py`) are also in `refactored_codes_v1/`; they survived for compatibility with the contour batch and the earlier paper-styling rollout. None are blockers.
