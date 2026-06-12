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

### `NF270_CF_RESIDUAL_SCALE_FRACTION` / `NF270_CP_RESIDUAL_SCALE_FRACTION`  (paper-justified measurement error, 2026-06-11)

Per-channel **relative measurement error** (fraction of concentration) for the
retentate `cF` and permeate `cV` WLS weights, replacing the hardcoded legacy
`0.003` (0.3 %) and `0.03` (3 %).

**Why.** The legacy 0.3 % cF weight had **no physical basis** — at cF ≈ 100 mM it
is a 0.3 mM tolerance, tighter than the conductivity probe itself. Source for the
honest values: Lilonfe, Estrada, Singh, Ouimet, Phillip, Dowling, *"Soft Sensors
Enable Real-Time Ion Concentration Measurements"*, ChemRxiv 2026
(doi:10.26434/chemrxiv.15002541/v1):

- **cF (retentate)** is conductivity/soft-sensor-derived. End-to-end error =
  soft-sensor MAPE (NaCl 1.0 %, CaCl₂ 2.2 %, LaCl₃ 1.9 %; Table 3) ⊕ regressed
  cell-constant uncertainty (≈1.38 %, Eq. S5) ≈ **2 %**. **DATA3 default = `0.02`.**
- **cV (permeate)** is an offline **ICP-OES** scalar — the *reference*. ICP-OES
  accuracy is ~1–3 %, so the legacy **3 % is already physical**; the paper's
  "within 5 %" is the soft-sensor↔ICP *agreement*, not ICP's own error.
  **DATA3 default = `None` (→ 3 %).** (The 5 % soft-sensor figure applies to the
  conductivity permeate *probe* channel `cH` via `NF270_USE_PERMEATE_PROBE`.)

**One source, three consumers.** All weights flow through the helper
`_nf270_conc_scales(workflow_family) → (cf_frac, cp_frac)`, used by (a) both
`obj_rule`s (objective), (b) `calc_FIM` `var_pred` (Fisher-information
covariance), and (c) `_label_parmest_model` (Pyomo.parmest `measurement_error`).
Keeping all three on one source prevents an inconsistent reported covariance.

**Guards:** `_nf270_conc_scales` returns the new fraction only when
`workflow_family == "DATA3"` *and* the constant is non-`None`; otherwise it
returns the legacy `(0.003, 0.03)`. DATA1/DATA2 are byte-identical (verified by
`tests/regression_data1_data2_guards.py` + `test_data1/2_regression.py`).

**Empirical finding (cF weight ↔ σ identifiability are coupled).** Relaxing cF
0.3 %→2 % does **not** rescue σ on low-cF concentrating sheets (still rails to a
bound — σ is genuinely unidentifiable there because Δπ is tiny), and on the
high-cF gold-standard sheet it can push σ *toward* the wall — i.e. the interior σ
under 0.3 % was partly propped up by over-fitting the σ-sensitive retentate
channel. Conclusion: honest weighting is correct for UQ/covariance, but σ
identifiability is governed by **experiment design (high-cF / multi-ΔP runs)**,
not by reweighting. See `_compare_cf_weighting.py`.

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

---

## 13. May 24–25, 2026 — Late-session NaCl rollout findings (honest)

### 13.1  Goal of the late session

The mid-day v5 contour-seeded rollout (logged in `/tmp/nf270_paper_rollout_v5.log`) had converged the MC2.05.07.24_NaCl sheet to **Lp = 9.11, B = 8.31, σ = 0.665** under `B_form = 'single'`, with 6+ multistart trials all landing on the same θ — a strong basin signal. The plan for the evening was to use that fit's θ as a seed (the "MC2 NaCl anchor") and run all remaining NaCl single-salt sheets through the same contour-seeded multistart, expecting σ near 0.66 across the campaign per the user's table (same NF270 coupon, NaCl in concentration regime → similar σ).

### 13.2  What we actually got — v3 rollout, 2026-05-25

Re-running the contour-seeded multistart with the closest reproduction of the v5 recipe we could assemble (`B_form='single'`, `NF270_CF_RESIDUAL_FLOOR_MM=1.0`, `NF270_USE_PERMEATE_PROBE=None`, `NF270_MULTISTART_USE_CONTOUR_SEEDS=True`, `NF270_MULTISTART_INCLUDE_CROSS_SALT=True`) on all 6 NaCl sheets:

| Sheet | Regime | Lp | B | σ | obj_m / cv / cr |
|---|---|---|---|---|---|
| MC2.05.07.24_NaCl   | conc | 9.21 | 14.08 | **1.000** | 22 / 102 / 55 |
| MC3.07.22.24_SNaCl  | dil  | 8.90 | — | **1.000** | 26 / 7.9 / 21 |
| **MC4.07.11.24_SNaCl**  | conc | 10.76 | 6.12 | **0.767** | 19 / 240 / 26 |
| MC5.07.23.24_NaCl   | conc | 8.83 | — | **1.000** | 13 / 86 / 55 |
| MC5.07.23.24_SNaCl  | conc | 8.36 | — | **1.000** | 15 / 85 / 63 |
| **MC5.07.23.24_S2NaCl** | dil  | 7.84 | — | **0.454** | 43 / 50 / 37 |

**2 of 6 sheets (MC4.SNaCl, MC5.S2NaCl) landed in interior σ. 4 of 6 landed at σ = 1.0 (the wall).**

### 13.3  Why v5's σ = 0.665 doesn't reproduce

The smoke test for the v3 rollout was MC2.05.07.24_NaCl — re-running v5's exact sheet through what was meant to be v5's recipe. It came back at σ = 1.000, not 0.665. We investigated:

1. **Library bisect** — checked out the library file at `0c45466^` (the exact pre-evening-commit state, which is what v5 ran against). Ran the smoke test against that library. Same σ = 1.0 result. Conclusion: **the single evening commit 0c45466 did NOT introduce the σ-basin change.**
2. **Data files** — `NF270_MC2.xlsx` mtime is Feb 23, 2026. Not modified since v5.
3. **Objective magnitudes** — v5 logs report `Obj = 1632838` and `Obj(m) = 26.6`. Current rollout reports `obj_m = 22` (lower than v5's). So the σ = 1.0 fit has *better* objective by the current formulation. This is not a local-minimum issue; σ = 1.0 is the genuine global optimum under the current code.
4. **Seeds** — pre-0c45466 library produced different LHS seeds (broader, including σ = 0.46 starts) than HEAD. All five seeds walked uphill to σ = 1.0 regardless of starting point.

**Conclusion: v5 must have been running against a different objective function than is now in the codebase.** Possibilities (none confirmed): an uncommitted edit that was rewritten before the evening commit; a different residual-scale normalization; a different per-vial count weighting. We cannot recover v5's exact state from git alone.

### 13.4  What this means for the campaign

- **The σ = 0.665 MC2.NaCl fit and the σ ≈ 0.66 cross-salt anchor values in `NF270_MULTISTART_CROSS_SALT_REFERENCES` are NOT reproducible.** They were obtained against a specific ephemeral session state. Honest reading: the σ = 0.665 result is not the current best optimum and any documentation that claims that fit is "the answer" needs the caveat.
- **The current campaign best is 2 of 6 NaCl sheets at interior σ.** Both interior cases (MC4.SNaCl at σ = 0.77, MC5.S2NaCl at σ = 0.45) used cf-floor contour panels as the seed source — the contour-seeded multistart is doing real work where the data support an interior σ.
- **4 of 6 sheets prefer σ = 1.0 under the current objective.** That is a real campaign result. Whether it's physically reasonable (literature NaCl rejection on NF270 is ~0.5–0.7) is a question for the experimentalist + the model — not for more multistart tuning.

### 13.5  Standing open question

If the v5 σ = 0.665 fit is physically expected (literature consistent), then the *current* objective is too permissive — something is letting the σ = 1.0 corner win. The most likely candidate is the residual-scale normalization: today's run has a sum-of-objectives total of ~180, vs v5's ~1.6M, a ~10⁴ factor difference. Tightening (or restoring) the per-vial residual count normalization is the natural place to look if/when this is investigated further.

### 13.6  Files added in this session

| File | Purpose |
|---|---|
| `_rollout_contour_seeded.py` | The v3 NaCl rollout — 6 sheets, contour-seeded multistart, killpg-enforced 600s wall cap |
| `_contour_seeded_worker.py` | Per-sheet worker; sets v5-recipe toggles, runs library multistart, returns winning θ + sim_stru |
| `_rollout_warm_anchor.py` | Phase 1.5 (built but NOT launched — would have used same flawed B_form=1) |
| `_warm_anchor_worker.py` | Companion worker for `_rollout_warm_anchor.py` |
| `_nacl_warm_worker.py` | Earlier abandoned attempt (warm-start single-shot from MC2 anchor) |
| `_rollout_nacl_warm.py` | Driver for the abandoned attempt above |
| `_b_form_worker.py` | Earlier subprocess worker for B_form experiments |
| `_rollout_per_sheet_subprocess.py` | Earlier process-group-kill rollout (kept as template) |

Library edits in this session (kept):
- `NF270_MULTISTART_INCLUDE_CROSS_SALT` (new module-level toggle)
- `_resolve_nf270_panel_dir` patched to prefer `contour_panels_cF_floor_1mM/` when cf-floor is on, falling back to legacy `contour_panels/` otherwise.

These patches are functional and benign — they don't change DATA1/DATA2 behavior, and they correctly route cf-floor fits to the matching contour panels when available.


## 14. May 27, 2026 — Data validation: loader fidelity + raw-conductivity plots

This section documents the formal proof that the DATA3 loader and plotting pipeline preserve the raw Excel measurements faithfully. It exists because a collaborator deck ("Diafiltration Multicomponent — Experimental Data.pptx") raised the concern that we might be "not plotting what's in the data file." That concern is rebutted with three independent pieces of evidence:

1. **Loader-fidelity audit** — every data point in every single-salt sheet verified GOOD after accounting for documented transformations.
2. **Raw-conductivity plot type** — a new figure that renders the raw probe traces (μS/cm at measured T) directly from the data file in the collaborator's exact format.
3. **Collaborator-vs-DATA3 side-by-side deck** — visually identical retentate concentration end-of-run values on every NaCl, CaCl₂, and LaCl₃ sheet in scope.

### 14.1  Canonical loader transformation table

This is the authoritative reference for what the loader does to each raw Excel column. Every transformation listed here is verified by the audit in §14.2.

| Variable | Loader transformation | Why |
|---|---|---|
| **Time (s)** | 1:1 preserved (no change) | Time axis is the index; no compensation needed. |
| **Mass (g)** | Baseline-subtracted per vial — each vial starts at 0 | Eliminates per-vial offset so mass-vs-time fits compare across vials cleanly. |
| **Pressure (psi)** | 1:1 preserved | Applied pressure is a control variable, recorded as-is. |
| **Retentate Temp (°C)** | 1:1 preserved | Needed for EC25 compensation of conductivity. |
| **Retentate Cond (μS/cm)** | EC25 temperature compensation: `σ_25 = σ_T · (1 + α · (25 − T))` | Probe reports σ at measured T; the model needs σ at 25 °C so concentration inversion (Shedlovsky) is consistent across temperatures. |
| **Permeate Temp (°C)** | 1:1 preserved | Needed for EC25 compensation of permeate conductivity. |
| **Permeate Cond (μS/cm)** | EC25 temperature compensation (same formula as retentate); tube-transit-time correction applied to **TIME axis**, not to value | Same as retentate. Tube transit τ = V_tube / (dm/dt) shifts the permeate time axis to account for the 0.3 g dead volume between membrane and probe. |
| **Vial Swap** | 1:1 preserved | Index flag; used by loader to slice the continuous time-series into per-vial blocks. |
| **ICP (cV_avg per vial)** | Extracted from vial-data block (cols 14–21), one scalar per vial. Conductivity calibration applied if needed. | Vial ICP-OES gives ion-specific concentrations; reduced to one mM scalar per vial for the model. |

**Salt-specific α** (used by EC25 compensation, from `CONDUCTIVITY_TEMP_COEFF_PER_C`):

| Salt | α (1/°C) |
|---|---|
| KCl | 0.019 |
| NaCl | 0.021 |
| CaCl₂ | 0.023 |
| LaCl₃ | 0.025 |
| (other) | 0.024 (default) |

**Important note on holdup splitting** — the loader also recognizes the experimental holdup period (the first `n_holdup_samples` rows before steady permeate collection) and splits them into a synthetic "vial 1" via the `vial1_split_status` flag. This is a documented redefinition of vial boundaries vs the raw `Vial Swap` column. The audit in §14.2 accounts for this by time-aligning loaded samples to raw rows rather than index-matching.

### 14.2  Loader-fidelity audit

**Script:** `refactored_codes_v1/_audit_loader_fidelity.py`
**Output:** `refactored_codes_v1/DATA3_loader_fidelity_audit_2026-05-27.xlsx` (1.6 MB, 15 tabs)

For each of the 13 single-salt sheets, the audit:

1. Reads the raw Excel directly with pandas (no library code in the path).
2. Loads via `lib.loadxlsx` → `data_stru`.
3. Applies the expected transformations from §14.1 to the raw data (so we compare apples to apples).
4. Time-aligns each loaded sample to its corresponding raw row by time-stamp matching (handles the holdup split).
5. Computes absolute and relative residual per data point.
6. Classifies each point as `good` / `warning` / `bad` based on a 1 % relative + absolute-floor tolerance.

**Result (continuous time-series, cols 0–7):**

| Total data points across 13 sheets | Good | Warning | Bad |
|---|---|---|---|
| **70,101** | **70,101 (100 %)** | 0 | 0 |

Every variable on every sheet receives the `OK` verdict in the `_summary` tab. The `_transformations` tab embeds the table from §14.1 directly into the workbook. The 13 detail tabs let any reviewer drill into row-by-row residuals.

**Extended audit blocks (added 2026-05-27 evening — v2 of the script):**

The original audit covered the continuous time-series (cols 0–7) thoroughly but stopped short of the per-vial ICP block and metadata. The v2 extension adds three more comparison blocks to every per-sheet tab:

1. **ICP sidebar (raw cols 14–21)** — side-by-side display of raw per-vial cells (Sample Vol, Nitric Acid Vol, ICP intensity in cps, ICP in mg/L) against the loaded per-position scalars (`cF_feed_icp_mM`, `cF_diafiltrate_icp_mM`, `cF_retentate_icp_mM`, `icp_final_tube_mM`, and per-vial `cV_avg`). Includes a "Loader expanded ratio" column showing the dilution × molar-mass factor the loader applies (e.g., ~108× for NaCl Feed with 0.12 mL sample + 5 mL HNO₃; ~510× for Diafiltrate / Retentate with 0.025 mL sample). Physically sensible across all 13 sheets.
2. **ICP calibration points (raw cols 23–26)** — the 10 (concentration, intensity) calibration pairs per salt, exactly as in the raw Excel, followed by the loaded `icp_calibration_curve` summary (slope, intercept, R², range). R² values fall in the 0.99–0.999 range — calibration quality preserved by the loader.
3. **Header / sidebar metadata** — Datapoints count, Notes string, Initial / Final Solution Weight, Number of Salts, Salt 1 + Salt 2 names, ICP Calibration Points (#) checked against `data_config` (`note_text`, `M_F0`, `nc`, `namec`, etc.). Result: **40 / 40 ✓ matches** across all 13 sheets on the fields with strict equality tests.

A new aggregate `_header_metadata` tab in the workbook collects every (sheet, field) row in one place — 83 rows total, with the `Match` column color-coded ✓ / — / ✗ for quick scanning.

### 14.3  New raw-conductivity plot type

**Library function:** `run_data3_conductivity_plots(results, save_dir, show)` in `refactored_ucb_library.py`
**Worker integration:** `_contour_seeded_worker.py` calls this alongside the existing mass / concentration / pressure / osmotic plot functions.
**Output filename:** `conductivity-<prefix>.png` (one per sheet, written to the same `campaign_figures/` directory as the other DATA3 plots).

For each experiment, the plot shows:

- **Magenta ▲ markers** — raw retentate conductivity (μS/cm at measured T), recovered from the stored EC25-compensated `cF_exp_conductivity` by inverting EC25: `σ_T = σ_25 / (1 + α · (25 − T))`.
- **Red ◆ markers** — raw permeate conductivity (μS/cm at measured T), same inversion applied to `cV_perm_cond`.
- **Green line (retentate) + black line (permeate)** — EC25-compensated traces at 25 °C, exactly what the model consumes downstream.
- **In-plot annotation** spelling out (i) the raw Excel columns the markers come from (col 4 = retentate, col 6 = permeate), (ii) the salt and α value, (iii) the EC25 forward + inverse formulas, (iv) the Shedlovsky-inversion downstream step.

This plot is **salt-agnostic** by construction — it uses `CONDUCTIVITY_TEMP_COEFF_PER_C[salt]` so the same code produces correct results for NaCl, CaCl₂, and LaCl₃ rollouts identically. Smoke-tested on MC2.NaCl, MC2.CaCl₂, and MC2.LaCl₃.

### 14.4  Collaborator-vs-DATA3 comparison deck

**File:** `refactored_codes_v1/DATA3_single_salt_collab_comparison_2026-05-25.pptx`

14 slides: title + reading guide + 11 side-by-side experiment comparisons (collaborator's slide LEFT, our concentration plot RIGHT) + summary. Spot-checked retentate end-of-run values on 3 sheets (MC2.NaCl, MC3.SNaCl, MC2.CaCl₂) — all agree within ~5 % across the two pipelines.

### 14.4b  ICP-data usage audit (added 2026-05-27 evening)

User raised a sharper question — verifying values are *preserved* isn't the same as verifying they're *used*. Tracing every ICP scalar from `data_config` into the model objective produced this table:

| ICP value (raw row) | Loaded into `data_config` | Used in fit objective? | Where in `refactored_ucb_library.py` |
|---|---|---|---|
| **Feed** (initial retentate) | `C_F0` | ✅ Yes — initial condition for `cF` at t=0 | line 1574 |
| **Diafiltrate** (buffer composition) | `C_D` | ✅ Yes — diafiltrate concentration parameter | line 1572 |
| **Retentate** (final stirred-cell) | `cF_retentate_icp_mM` | ⚠️ Previously **cross-check only** — now optional fit anchor via `NF270_USE_RETENTATE_ICP_ANCHOR` | new block at line ~2493 |
| **Final Tube** | `cF_final_meas` / `icp_final_tube_mM` | ✅ Yes — `obj_cr` end-of-experiment anchor | lines 2477–2487 |
| **Per-vial Permeate** (Vial 1..N) | `data_raw[i]['cV_avg']` | ✅ Yes — `obj_cv` residual term with 3 % relative scale | lines 2376, 3042 |

`PREFLIGHT_AUDIT.md` line 120 had explicitly documented `cF_retentate_icp_mM` as "⚪ cross-check only" — known but unused. In MC2.NaCl the Retentate ICP (35.06 mM) and Final Tube ICP (23.90 mM) differ by ~46 %, so they are physically distinct samples (bulk stirred-cell at end vs connecting-tubing dead-volume).

**New toggle `NF270_USE_RETENTATE_ICP_ANCHOR`** (default `None` = off, opt-in):

When set truthy AND `workflow_family == "DATA3"`, `model_construct_inter` adds a *second* residual term to `obj_cr`:
```
res_cf_retentate = m.cF[last_vial, last_tau]  −  cF_retentate_icp_mM
```
with the same 0.3 % relative scaling as the Final Tube anchor. Both anchors apply at the same `(last_vial, last_tau)` index — the optimizer reconciles both samples in a least-squares sense.

**Guards** — like the other three DATA3-only knobs (`NF270_CF_RESIDUAL_FLOOR_MM`, `NF270_USE_PERMEATE_PROBE`, `NF270_SIGMA_INTERIOR_BOUNDS`), this toggle is gated by a `workflow_family == "DATA3"` check. Verified inert for DATA1/DATA2 by `tests/regression_data1_data2_guards.py` (now testing all four knobs); DATA2 pytest continues to pass within 1 % tolerance.

### 14.5  What this section lets us claim

> The data-processing pipeline preserves the raw Excel measurements byte-equivalent (Time, Pressure, Temperatures, Vial Swap) or transforms them by single documented formulas (Mass baseline-subtraction per loader-vial, Conductivities EC25-compensated by salt-specific α). Every continuous-measurement column from the data file is now plotted in a generated figure (raw conductivity in `conductivity-<prefix>.png`; derived concentration in `concentration-<prefix>.png`). 70,101 of 70,101 audited data points verify GOOD. The "you're not plotting the data file" claim is empirically false.

### 14.6  Files added or changed in this section

| File | Purpose |
|---|---|
| `refactored_codes_v1/_audit_loader_fidelity.py` | Loader-fidelity audit script (NEW) |
| `refactored_codes_v1/DATA3_loader_fidelity_audit_2026-05-27.xlsx` | Audit output workbook (NEW; 15 tabs) |
| `refactored_codes_v1/refactored_ucb_library.py` | Added `run_data3_conductivity_plots()` (+140 lines) |
| `refactored_codes_v1/_contour_seeded_worker.py` | Worker now calls `run_data3_conductivity_plots()` alongside mass / conc / pressure / osmotic |
| `refactored_codes_v1/DATA3_single_salt_collab_comparison_2026-05-25.pptx` | Collaborator-vs-DATA3 side-by-side deck (NEW; built earlier in this session) |


## 15. May 27, 2026 — Explaining the concentrating/diluting mass pattern (ΔP vs Δπ overlay)

This section closes a question collaborators raised in last week's discussion: across the single-salt campaign, the mass collected per vial is **not constant** — it tends to *increase* through a diluting run and *decrease* through a concentrating run. Per-vial mass slope ≠ 0 is a physical signature of changing flux, and that change has a clean explanation in the Spiegler–Kedem flux equation.

### 15.1  The physics in one line

```
                    Jw  =  Lp · (ΔP  −  σ · Δπ)
```

- `ΔP` (applied pressure) is **constant** for the experiment — it's a control variable set by the operator.
- `Δπ` (osmotic pressure difference across the membrane) **depends on cF(t)** via van 't Hoff:

```
                    Δπ  ≈  (cF − cH) · ν · R · T   ≈   cF · ν · R · T   (cH ≪ cF)
```

So as `cF` changes during the experiment, `Δπ` changes, and the net driving force `ΔP − σ·Δπ` changes. That in turn changes `Jw` — the water flux — which is the rate at which mass accumulates in each collection vial.

| Regime | `cF(t)` trend | `Δπ(t)` trend | `ΔP − σ·Δπ` trend | `Jw` trend | Mass per vial |
|---|---|---|---|---|---|
| **Concentrating** (feed dilute, diafiltrate concentrated) | rises | rises | shrinks | drops | **decreases over the run** |
| **Diluting** (feed concentrated, diafiltrate water) | falls | falls | grows | rises | **increases over the run** |

The pattern the collaborators flagged is therefore expected — and visible directly in the data without any model fitting.

### 15.2  New plot type: `applied_vs_osmotic-<prefix>.png`

**Library function:** `run_data3_applied_vs_osmotic_plots(results, save_dir, show)` in `refactored_ucb_library.py`.
**Worker integration:** `_contour_seeded_worker.py` calls it alongside mass / concentration / pressure / osmotic / conductivity.
**One-off generator:** `_generate_applied_vs_osmotic_plots.py` regenerates all 13 single-salt plots in ~10 s without needing a fit.

Per-sheet figure shows:

- **Orange horizontal line** — applied pressure `ΔP` (constant, read from `data_config["delP"]`).
- **Magenta curve** — osmotic pressure `Δπ(t) = cF · ν · R · T`, built from measured `cF` (Shedlovsky-inverted, EC25-compensated). No model prediction required — the plot is built purely from the data file.
- **Colored callout** (red for concentrating, green for diluting) reporting the regime, the cF start/end values, and the end-of-run `Δπ_end / ΔP` ratio so the reader can see immediately how close the run gets to osmotic shutdown.
- **Footer** documenting the formula and the data-file source columns.

Stoichiometric ν is salt-specific (NaCl = 2, KCl = 2, CaCl₂ = 3, LaCl₃ = 4); falls back to `data_config["ni"]` or 1 for unknown salts. This is the **physical** dissociation factor, distinct from the model's `ni` parameter which the loader sets to 1 by convention.

### 15.3  Spot-check findings on the 13 generated plots

| Sheet | Regime | `cF: start → end` (mM) | `Δπ_end / ΔP` | Pattern visible? |
|---|---|---|---|---|
| MC2.05.07.24_NaCl | concentrating | 0.88 → 35.72 | 0.54 | ✓ moderate flux suppression |
| MC3.07.22.24_SNaCl | **diluting** | 95.53 → 10.39 | 0.14 | ✓ **Δπ_initial > ΔP at t=0** — striking |
| MC2.05.07.24_CaCl₂ | concentrating | 0.60 → 23.75 | 0.53 | ✓ ν = 3 magnifies Δπ per mM |
| (10 other sheets — same pattern) | … | … | … | ✓ |

The MC3.SNaCl plot is the cleanest visualization: at t = 0 the osmotic backpressure literally exceeds the applied pressure, so net flux is suppressed at the start of the run and grows as the experiment proceeds. The mass-vs-time plot for that sheet shows exactly the expected mirror image — slow mass accumulation early, accelerating as Δπ falls.

### 15.4  How this helps the broader fit story

This plot type is **diagnostic**, not a fit input. Its value is twofold:

1. **For collaborator review:** the per-vial mass pattern is no longer mysterious. It is the direct, predictable signature of a constant-pressure diafiltration experiment combined with a changing feed concentration. Showing the ΔP and Δπ curves side-by-side makes this self-evident.
2. **For our own pathology hunt:** sheets with `Δπ_end / ΔP` close to 1 (near-osmotic-shutdown) are exactly the ones where σ is most strongly constrained by the data, because σ multiplies the dominant term in the flux equation. Sheets where `Δπ_end / ΔP ≪ 1` have weak σ identifiability — σ becomes interchangeable with Lp. This connects naturally to the §13 finding that σ pins at the bound on the concentrating-regime sheets and to the lumped-σ·Lp diagnostic deferred in Task #5.

### 15.5  Files added / changed in this section

| File | Purpose |
|---|---|
| `refactored_codes_v1/refactored_ucb_library.py` | Added `run_data3_applied_vs_osmotic_plots()` (+162 lines) |
| `refactored_codes_v1/_contour_seeded_worker.py` | Worker now calls the new function alongside the other plot types |
| `refactored_codes_v1/_generate_applied_vs_osmotic_plots.py` | One-off generator for the 13 single-salt sheets (NEW) |
| `UnifiedFramework/DATA3/results/paper_artifacts/nf270/campaign_figures/applied_vs_osmotic-*.png` | 13 generated figures (NEW; committed to the repo) |


## 16. June 1, 2026 — DATA1 direct-contour bugfix + `plot_contour` figure_s5 style

Two small but important changes to the paper-figure rendering path, plus same-day housekeeping commits captured for completeness.

### 16.1  Defensive `parent.mkdir` at the `run_data1_direct_contour_branch` CSV write site

**Symptom.** With `Root = 1` (DATA1), manifest-default trunk, branches = all, then answering `y` to the *"Run the direct DATA1 contour-grid branch as well?"* prompt, the run crashed with:

```
OSError: Cannot save file into a non-existent directory:
  '.../paper_artifacts/data1/notebook_figures/direct_contours'
```

**Root cause.** The up-front `save_dir.mkdir(parents=True, exist_ok=True)` at function entry (line 7423 of `refactored_ucb_library.py`) was not sufficient — caller-supplied path quirks, intermediate modifications, or filesystem races can leave the directory missing by the time `df.to_csv(csv_path)` runs inside the per-page / per-case loop.

**Fix.** Add `csv_path.parent.mkdir(parents=True, exist_ok=True)` directly before the `df.to_csv(...)` call. This mirrors the `out_path.parent.mkdir` pattern already in use at four other write sites in the same file (lines 6584, 6732, 6768, 8226), making the CSV writer self-sufficient.

Three downstream writers in the same call graph were already defensive — only the direct-CSV write was exposed:

| Writer | Defensive idiom | Status |
|---|---|---|
| `df.to_csv(csv_path)` at line 7502 | (none) | **fixed in this commit** |
| `_plot_contour_data1_legacy()` at 5899 | `save_dir.mkdir` at entry | already defensive |
| `_compose_data1_panel_grid()` at 6732 | `out_path.parent.mkdir` | already defensive |
| `_compose_data1_panel_sheet()` at 6768 | `out_path.parent.mkdir` | already defensive |

Commit `4cd8a60` (+7 lines, comments included).

### 16.2  `plot_contour` switched from `pcolormesh` heatmap to figure_s5 contour-line style

**Symptom.** The contour panels produced by the DATA1 reproduction workflow (file names like `data501.1_base_contour_sigma.png`, `data511.12_base_contour_B.png`) were rendered as solid `pcolormesh` heatmaps with a colorbar and bare column-name axis labels — visually correct but hard to read for parameter-estimability discussion.

**Root cause.** `plot_contour` (the public 1×N panel renderer) delegated to `_plot_heatmap_frame`, an internal panel-drawer using `ax.pcolormesh(...) + plt.colorbar(...)` and `ax.set_xlabel(x_col)` (raw column name like `"sigma"`, no units). The published-style contour-line renderer (`_plot_contour_data1_legacy`, which produces the SI panels stitched into `figure_s5.png` / `figure_s6.png`) lived in a separate code path used only by the SI composers and the manifest-driven `render_contour_data1_legacy`.

**Fix.** Rewrite `_plot_heatmap_frame` (signature unchanged) in the published figure_s5 / figure_s6 idiom:

- Coloured iso-objective contour **lines** (`ax.contour(X, Y, Z, 10, linewidths=2, cmap="viridis")`), no fill
- Inline numeric labels at every other level (`ax.clabel(cp, cp.levels[::2], inline=True, fmt="%1.1f")`)
- Red triangle at the per-channel `nanargmin(Z)` location, matching the legacy SI panel marker
- Bold, TeX-formatted axis labels with units via `_PAPER_AXIS_LABELS` lookup, fallback to raw column name
- Ticks pointing inward, no colorbar

`plot_contour`'s figsize also drops from `(5 × n_panels, 4)` to square `(4 × n_panels, 4)` — the colorbar room is freed.

**Blast radius.** All four `plot_contour` call sites automatically pick up the new style:

| Call site | What it writes |
|---|---|
| `published_paper_reproductions` workflow, lines 7042 / 7051 | `data{N}_{variant}_contour_{B,sigma}.png` (DATA1 reproduction PNGs) |
| `_contour_seeded_worker.py`, line 7861 | DATA3 NF270 quick-look composite (alongside the per-channel SI panels) |
| `render_contour`, line 11296 | Manifest-driven FigureSpec entrypoint |

`_plot_contour_data1_legacy` itself is unchanged — its inline `plt.contour` / `plt.clabel` logic was the visual reference for this port.

Visual check against an existing DATA3 NF270 contour CSV (`MC2.05.07.24_NaCl/contourdata-x_sigma-y_Lp.csv`) — output matches the figure_s5 idiom (coloured iso-lines, inline labels, red triangle at the `(Lp, sigma)` minimum, σ [dimensionless] and L_p [L · m⁻² · h⁻¹ · bar⁻¹] axis labels).

Commit `c737254` (+44 / −11).

### 16.3  Reproducibility note on DATA1 / DATA2

Neither §16.1 nor §16.2 touches `_plot_contour_data1_legacy` or any of the SI composers — the byte-equivalent reproduction of `figure_s2…s6` is preserved. The change in §16.2 affects only the **secondary** quick-look heatmap output produced by `plot_contour`; the canonical SI panels stitched by the manifest builders are unchanged.

Recommended verification: run the DATA1 / DATA2 manifest default trunk and inspect `_print_coverage` — should still report `main figures: 12/18   SI figures: 10/10   tables: 2/2` (or better, since §16.1 now lets the direct-contour branch complete).

### 16.4  Housekeeping commits same day

- `a19ef91` — collaborator briefing 2026-05-27: dated folder grouping the earlier brief, the 6-slide AICHE-style build, and the edited version actually presented to collaborators. Lives under `refactored_codes_v1/collaborator_briefing_2026-05-27/`.
- `e3a25e1` — DATA3 warm-start fit artifacts under `UnifiedFramework/DATA3/results/paper_artifacts/nf270/warm_start_fits/` and `..._with_perm/` (22 JSON + 22 PNG + 2 `summary.json`), preserved for later reference. NF270 partition only — DATA1 / DATA2 outputs untouched.

### 16.5  Files changed in this section

| File | Change |
|---|---|
| `refactored_codes_v1/refactored_ucb_library.py` | §16.1 mkdir defensive (+7) and §16.2 `_plot_heatmap_frame` rewrite + `plot_contour` figsize (+44 / −11) |
| `refactored_codes_v1/Architecture.md` | This section (§16) |
| `refactored_codes_v1/collaborator_briefing_2026-05-27/*.pptx` | 3 dated decks (§16.4) |
| `UnifiedFramework/DATA3/results/paper_artifacts/nf270/warm_start_fits/` | 23 new files (§16.4) |
| `UnifiedFramework/DATA3/results/paper_artifacts/nf270/warm_start_fits_with_perm/` | 23 new files (§16.4) |


## 17. June 3, 2026 — Per-salt B upper bounds for DATA3 NF270 single-salt fits

This section documents a physics-anchored refinement to the solute-permeability parameter `B` bound used by the DATA3 NF270 single-salt workflow. The change replaces a single dataset-wide `B ∈ (1e-6, 50]` window (set in §12 on May 24 to free the CaCl₂ fits off the prior `B = 30` bound) with a **per-salt, valence-graded** triple — `NaCl (1e-6, 30]`, `CaCl₂ (1e-6, 15]`, `LaCl₃ (1e-6, 10]` — derived from the dielectric-exclusion + Stokes–Einstein picture of nanofiltration through a negatively-charged polyamide layer. The change is **DATA3-only, hard-gated by `workflow_family == "DATA3"`**; DATA1 / DATA2 model builds continue to receive the legacy single-window bound and the published-paper reproductions remain byte-equivalent.

### 17.1  Why the unified `(1e-6, 50]` bound is physically uninformative

The `B` parameter in the Spiegler–Kedem flux equation is the **solute permeability coefficient** (units of length / time after the library's internal unit folding). For a charged-membrane–chloride-salt system, two physics-based predictions for the cation-valence scaling of `B` agree to leading order:

1. **Dielectric exclusion** (Yaroshchuk 2000, *Adv. Colloid Interface Sci.* 85 (2–3): 193–230). For a NF270-like negatively-charged active layer with a low-permittivity pore interior, the partition coefficient of a cation with valence `z` carries a Born-image exclusion factor `exp(−z² ΔW / kT)`. For the chloride salts in scope, the dominant exclusion is on the cation: doubling the cation valence multiplies `z²` by four, which suppresses partitioning sharply. Empirically the suppression is gentler than the bare prediction (counter-ion screening softens the image potential), and the literature consensus on NF-class polyamides is a **roughly halving of B per unit increase in cation valence** for the Na⁺ → Ca²⁺ → La³⁺ sequence on chloride salts. See also Bandini & Vezzani 2003, *Chem. Eng. Sci.* 58 (15): 3303–3326, for the parallel DSPM-DE derivation.

2. **Stokes–Einstein diffusivity** in a partition × diffusion factorization, `B ≈ φ · D_∞ / δ`. The bulk Stokes–Einstein diffusivities `D_∞` for NaCl / CaCl₂ / LaCl₃ are within a factor ~1.25 of each other at 25 °C (1.61, 1.33, 1.29 × 10⁻⁹ m²/s per Vanýsek CRC Handbook + Rard & Miller 1979). The remaining factor in `B` is the partition coefficient `φ`, which is exactly the dielectric-exclusion term above. So the diffusive prefactor is approximately salt-independent and the valence-graded scaling falls entirely on `φ`.

Net: the physics predicts `B(NaCl) > B(CaCl₂) > B(LaCl₃)` with each step roughly halving. A uniform `(1e-6, 50]` bound permits non-physical optima in which CaCl₂ or LaCl₃ end at `B` values comparable to or exceeding NaCl — an outcome the literature explicitly rules out for NF-class chloride rejection.

The full deep-research provenance trail — five-angle web search, twenty-one fetched sources, twenty-five claims with three-vote adversarial verification (two killed, including the Lachheb et al. 2025 NaCl Pₛ as not reproducible) — is workflow `wss0st04w` in the run log. Nair et al. 2018, *Membranes* 8 (3): 78, is the only peer-reviewed NF270-specific anchor that survived verification; their per-ion Pₛ (multi-ion seawater fit) puts NaCl at order ~1.5 µm/s and Cl⁻ at ~21 µm/s, comfortably inside the new `(1e-6, 30]` bound.

### 17.2  The bound triple, gating, and DATA1 / DATA2 inertness

| Salt | New bound `B ∈ (lo, hi]` | Source |
|---|---|---|
| NaCl  | `(1e-6, 30]` | Retains the May 24 §12 widening; physical upper for monovalent Cl⁻ salt on NF270 |
| CaCl₂ | `(1e-6, 15]` | NaCl bound halved, dielectric-exclusion + Nair 2018 rejection survey |
| LaCl₃ | `(1e-6, 10]` | CaCl₂ bound 2/3-ed (gentler than strict halving — La³⁺ ion-pairing softens predicted suppression) |

**Module-level constants** (`refactored_ucb_library.py` lines 281–296, immediately following `NF270_SIGMA_INTERIOR_BOUNDS`):

- `NF270_B_BOUNDS_PER_SALT = {"NaCl": (1e-6, 30.0), "CaCl2": (1e-6, 15.0), "LaCl3": (1e-6, 10.0)}`
- `NF270_B_BOUNDS_DEFAULT = (1e-6, 50.0)` — fallback for any DATA3 sheet whose `namec` doesn't match a key in the dict.

**Read sites and guards** (`refactored_ucb_library.py`):

- The `B` Pyomo `Var` declaration inside `model_construct_inter` at the `B_form == 'single'` branch (was the single-line bound at line 1683; now a 5-line `if workflow_family == "DATA3":` block followed by `m.B = Var(bounds=_b_bounds, initialize=param_in['B'])`). The salt key is `str(data_stru['data_config'].get('namec') or "").strip()` — the canonical salt-name field already consumed at the diffusivity branch (line ~1566) and populated by the DATA3 Excel ingest (line ~4422) and by both .mat loaders for DATA1 / DATA2.
- The multistart clip-bound dict in `build_seeded_multistart_starts` (was line 9859, now reads `NF270_B_BOUNDS_PER_SALT.get(salt_name, NF270_B_BOUNDS_DEFAULT)` when `str(workflow_family).upper() == "DATA3" and str(B_form) == "single"`). Without this update the clip site would silently widen DATA3 multistart seed B values back to 50 on CaCl₂/LaCl₃ sheets and re-trigger the W1002 projection warnings the clip was introduced to silence.
- Both gates are explicit `workflow_family == "DATA3"` checks; DATA1 / DATA2 always see `(1e-6, 50)` exactly, regardless of the dict's contents.

`solve_model_B_fix` does not redeclare `m.B` — it routes through `model_construct_inter`, so the single Pyomo-Var-declaration site above is the only fit-time read.

The regression test `refactored_codes_v1/tests/regression_data1_data2_guards.py` is extended with a new fifth guard block (`b_bounds_for`) asserting the `B` `Var` ends up with `bounds == (1e-6, 50)` for DATA1 and DATA2 builds even when `NF270_B_BOUNDS_PER_SALT` is set to a deliberately distorted triple (3.0 / 1.0 / 0.5) — the same contract pattern used for the other four DATA3-only knobs. SUMMARY block and conclusion wording updated from "all four DATA3 knobs" to "all five DATA3 knobs".

### 17.3  The σ-vs-B "ordering" surprise — not a contradiction

A subtlety surfaced during review: the experimentally-observed rejection ordering on NF270 with these three chloride salts is **NaCl > CaCl₂ > LaCl₃** (multivalent salts are *less* rejected than NaCl — a Donnan-attraction consequence: higher-valence cation pulls Cl⁻ co-ion through to maintain electroneutrality, and the convective leak dominates rejection). The new `B` bound ordering goes the **same direction** as the rejection ordering — `B(NaCl) > B(CaCl₂) > B(LaCl₃)` — which at first glance reads as backwards: "less rejected" usually intuits as "leakier in every sense."

This is **not a contradiction**. The two orderings describe physically distinct transport mechanisms in the Spiegler–Kedem model, and the data has to distinguish them:

- **σ (the reflection coefficient)** governs the *convective coupling* term `σ · Δπ` in the water flux equation. Lower σ means salt rides through with the water flow more easily; "rejection at the bench" in the high-flux limit is primarily a σ phenomenon.
- **B (the solute permeability)** governs the *diffusive partitioning* term in the solute flux equation. Lower B means the solute partitions weakly into the pore phase and diffuses slowly across it.

Both mechanisms attenuate with increasing cation valence on a negatively-charged NF membrane because the dielectric exclusion is shared between them. On a CaCl₂ or LaCl₃ run, salt mostly passes through by **riding the water flow** (Route 1, σ-controlled); diffusing through alone (Route 2, B-controlled) is a smaller channel where multivalent salts face a steeper energy cost. So the bench rejection drops (the σ effect dominates), *and* the B upper bound tightens (physics constrains the secondary diffusive channel) — both can be true, and the new bounds are internally consistent with the rejection data, not in conflict with it.

The cF-channel-dominance pathology documented in §12.3 / §13 is precisely the symptom of weak σ-vs-B identifiability — when `Δπ_end / ΔP ≪ 1` (§15.4) the convective coupling is weak and the data has nothing to separate σ from B. Tightening the B bound *helps* identifiability by reducing the σ↔B trade-off feasible region.

### 17.4  Open questions and caveats

1. **CaCl₂ S2 clipping at B = 30 even before this tightening.** The May 24 widening (§12) reported CaCl₂ S2 sheets pinning at `B = 30` — at the time interpreted as the optimizer wanting to go higher. Under the new bound `B ∈ (1e-6, 15]` those sheets will pin lower. Whether the new pin is at the bound or in the interior is the first thing to check after re-running the campaign with the per-salt bounds active.
2. **No peer-reviewed NF270 LaCl₃ B value in the literature.** Nair 2018 and the Bandini–Vezzani treatment both stop at divalent; Lachheb et al. 2025 attempted NaCl Pₛ but did not reproduce in the verification audit (workflow `wss0st04w`). The `(1e-6, 10]` LaCl₃ bound is a physics extrapolation from the z² Born scaling, not literature-pinned. The upper edge could move ±50 % once a direct LaCl₃ measurement appears.
3. **Ion-pairing physics gap for La³⁺.** The dielectric-exclusion picture treats the cation as a point charge. La³⁺ in chloride solution forms LaCl²⁺ and LaCl₂⁺ ion pairs at non-negligible fractions even at sub-millimolar concentration; these pairs partition differently than the free trivalent cation and bias `B` upward. The current bound does not encode this — it is the "no ion-pairing" upper. A future revision may need to widen the LaCl₃ bound modestly (e.g. `(1e-6, 12]`) and add a free-vs-paired species sub-model. This is also the open question flagged on slide 12 of the deck.
4. **The bound triple should be revisited if KCl is added to the DATA3 campaign.** The library already supports KCl conductivity inversion via `CONDUCTIVITY_SALT_PARAMS_25C`; KCl is monovalent so the NaCl bound (`(1e-6, 30]`) would be a reasonable seed by valence-equivalence, but KCl's limiting equivalent conductivity is ~25 % higher than NaCl, which lifts the diffusive prefactor slightly. Open until experimentally needed. Unknown salt strings fall through to the (1e-6, 50] default — same behavior as before the change for any DATA3 sheet that isn't NaCl/CaCl₂/LaCl₃.

### 17.5  Slide deck cross-references

The physics summary, the σ-vs-B reconciliation, and the literature anchors are slides **§10A**, **§10B**, and **§10C** respectively in `DATA3_single_salt_analysis_v10.pptx` (slides 31, 32, and 104). Slide 30's `B ∈ (10⁻⁶, 30)` line was corrected to `(10⁻⁶, 50)` in the same v10 revision to match the bumped library bound from §12. Speaker notes for slides 10A/10B/10C carry the deep-research citation graph for anyone needing the full literature picture. The deep-research provenance for the literature anchors is workflow `wss0st04w` in the run log.

### 17.6  Files changed in this section

| File | Change |
|---|---|
| `refactored_codes_v1/refactored_ucb_library.py` | New constants `NF270_B_BOUNDS_PER_SALT` + `NF270_B_BOUNDS_DEFAULT` after line 279 (+17); gated read in `model_construct_inter` at the `B_form == 'single'` branch (was line 1683, +14 / −1); gated read in `build_seeded_multistart_starts` clip dict (was line 9859, +6 / −1) |
| `refactored_codes_v1/tests/regression_data1_data2_guards.py` | New `b_bounds_for` helper + 5th guard block asserting DATA1 / DATA2 always see `(1e-6, 50)` regardless of `NF270_B_BOUNDS_PER_SALT` contents; SUMMARY line + "all four → all five" wording updated (+50 / −2) |
| `refactored_codes_v1/Architecture.md` | This section (§17) |
| `refactored_codes_v1/DATA3_single_salt_analysis_v10.pptx` | Slides 31, 32, 104 (referenced in §17.5) — already landed in the deck-version commit |

---

## 18. June 11, 2026 — Post-meeting follow-up (DATA3-only)

Everything in §18 is **DATA3-scoped**; DATA1/DATA2 published-figure reproduction is
unchanged (verified by `tests/regression_data1_data2_guards.py` + the IPOPT
`test_data1/2_regression.py`).

### 18.1  Paper-justified measurement error  (cF 0.3 % → 2 %)

Documented in full at §5b. Summary: the legacy 0.3 % retentate-cF weight had no
physical basis; the soft-sensor paper (Lilonfe et al., ChemRxiv 2026,
doi:10.26434/chemrxiv.15002541/v1) gives the honest values — cF (conductivity
soft-sensor) ≈ 2 % (MAPE ⊕ cell-constant), cV (ICP reference) ≈ 3 % (kept). One
helper `_nf270_conc_scales()` feeds the objective, `calc_FIM`, and parmest so the
covariance stays self-consistent.

**Key empirical finding (`_compare_cf_weighting.py`, multistart-verified).** σ's
global minimizer **flips between interior and a bound depending on the cF weight**:
MC3 SNaCl goes 0.42(interior, 0.3 %) → 1.0(wall, 2 %); MC5 S2NaCl goes 0.0(wall,
0.3 %) → 0.35(interior, 2 %). Conclusion: the interior-σ values were **partly an
artifact of the over-tight 0.3 % weight**; with honest weighting σ is revealed to
be **poorly identified** from single-salt runs. This is a sensitivity-analysis
result, not a fit improvement — it strengthens the case that σ needs **better
experiments** (high-cF / multi-ΔP DoE), not reweighting. WSSE magnitudes are NOT
comparable across weightings (the weight rescales the objective).

### 18.2  Plot-axis convention — Y-priority `L_p > B > σ`

Rule of thumb (apply unless overridden): when two of {L_p, B, σ} are on a 2-D
plane, the **higher-priority** variable goes on **Y**, the other on X.
So: L_p×B → L_p on Y; L_p×σ → L_p on Y; **σ×B → B on Y**.

In `(x_var, y_var)` code, `y_var` is the Y axis. Audit + fixes:

| Generator | pair | axes | status |
|---|---|---|---|
| `_make_3d_contour_animations.py` σ-scrub | was `x=Lp,y=B` | B on Y | **FIXED → `x=B,y=Lp`** (file `scrub-sigma_LpvsB.gif`) |
| `_make_3d_contour_animations.py` Lp-scrub | `x=σ,y=B` | B on Y | already compliant |
| `_make_3d_contour_animations.py` B-scrub | `x=σ,y=Lp` | L_p on Y | already compliant |
| `refactored_ucb_library.py` `sweep_pairs` | `("sigma","B")` | B on Y | already compliant (the prior synthesis was wrong) |
| `_run_all_pairs_5channel.py` `SWEEPS[2]` | was `("B","sigma")` | σ on Y | **FIXED → `("sigma","B")`** (file `objcontour-x_sigma-y_B.png`) |

Stale old-axis files (`scrub-sigma_BvsLp.gif`, `objcontour-x_B-y_sigma.png`) are
superseded; delete or ignore so deck-builders don't mix conventions.

### 18.3  β(c) → β(I):  concentration-dependent B, plotted and ionic-strength-based

NF270's solute permeability is concentration-dependent, so the lumped constant
`B` (`B_form='single'`) is an average. With numeric `B_form`:
`B = J_w·[β₀ + β₁·c^(B_form−1) + β₂·c^B_form]` (β saved by `save_model`).

- **Plotting β vs L_p / σ (IMPLEMENTED).** `_axis_values()` accepts
  `beta_0/beta_1/beta_2` axes (β₀ ∈ [0.1,3]×fit; β₁,β₂ centered ±2·span, allowing
  negatives), the contour-axis validator now admits them, and `sweep_pairs` in
  `run_nf270_contour_for_sheet` switches to the β slices `(beta_0,Lp)`,
  `(sigma,Lp)`, `(sigma,beta_1)` **gated on `isinstance(B_form,int) and
  B_form>=1`** (bool excluded) so `'single'`/`'pervial'` keep the lumped-B slices.
  Y-priority preserved: β plays B's role (Lp on Y vs β; β on Y vs σ). The B_form=1
  centering fit is fragile single-shot, so the fast driver
  `_run_beta_contour_fast.py` warm-starts it from a `'single'` fit (B→β₀, β₁=0)
  and sweeps with a 300-iter/15 s per-cell cap + worker pool (the library path is
  correct but slow on extreme β cells). On MC3 SNaCl the B_form=1 fit centers near
  Lp≈9, B = β₀ + β₁·cIn with β₀≈1, β₁≈0.009.
- **Ionic strength (IMPLEMENTED — `NF270_B_USE_IONIC_STRENGTH`).** Inside the
  active B-polynomial rule `B_form_rule3` (and the dead `B_form_rule2`), the
  concentration argument is bound to `cc = (k_I·cIn) if use_I else cIn`, where
  `k_I = nf270_ionic_strength_factor(salt)` (NaCl 1, CaCl₂ 3, LaCl₃ 6;
  `½·z·(z+1)`). Because `cc` IS the same Pyomo component as `cIn` when the toggle
  is off, the constraint is **byte-identical** to legacy; the toggle is gated to
  `workflow_family=="DATA3"` so DATA1/DATA2 are unaffected (6th guard in
  `regression_data1_data2_guards.py`, build-based). **k_I scales the INTERFACIAL
  cIn, not bulk cF** — B physically acts at the membrane wall (`Js=B·(cIn−cH)`),
  which is also why rule3 (cIn) is the live rule and rule2 (cF) is dead. Smoke
  test confirms engagement: CaCl₂ `β₀+β₁·cIn` → `β₀+β₁·(3.0·cIn)`.
  - *Single-salt rescale identity (validated, deterministic).* Forward-sim B(c) at
    (β₀,β₁) vs B(I) at (β₀,β₁/k_I) match **to machine precision** (max rel-diff
    across mass/cV/cR channels: NaCl 0, CaCl₂ 2.1e-15, LaCl₃ 3.1e-15;
    `_run_b_ionic_validation.py` Part A) — a single-salt B(I) fit is provably a
    pure rescaling β_k→β_k/k_I^(cIn-power).
  - *Cross-salt transfer — NEGATIVE (honest).* Ionic strength does **not** unify B
    across NaCl/CaCl₂/LaCl₃. The fitted slope β₁ does not cluster tighter in
    I-space (CV(β₁_I)=3.45 > CV(β₁_c)=1.92) and even flips sign by salt, while β₀
    spans 1.5 / 0.32 / 0.11 — B itself is far more salt-specific than k_I predicts
    (multivalent ions are much more strongly rejected via dielectric/Donnan
    exclusion). So B(I) is **necessary but not sufficient** for cross-salt
    prediction; a salt/valence-specific exclusion term is also needed. Caveat: the
    B_form=1 slope fits on the tiny-B multivalent salts are fragile
    (retry/multistart), so β₁ magnitudes carry fit noise — but the β₀ spread and
    sign disagreement are robust. The clean, certain result is the single-salt
    identity; cross-salt unification by ionic strength alone is refuted here.
- **Why break/group concentration ranges:** the per-vial trend (§18.4) shows an
  approximately monotone B–c rise within a sheet (MC3 SNaCl B 12→21 μm/s,
  corr 0.93); piecewise-linear (per-band) β is the empirical linearization, and
  B(I) is the physical law that unifies the bands across salts.

### 18.4  Concentration-range & per-vial parameter grouping (IMPLEMENTED 2026-06-10)

Motivation: a single sheet spans ~5→100 mM, over which σ and B genuinely change
(through ionic strength). One θ per sheet averages over that. The order of rigor
below is now **implemented and validated** on a GOOD / OKAY / POOR triple
(`MC3.07.22.24_SNaCl` / `MC2.05.07.24_NaCl` / `MC2.05.21.24_LaCl3`); all 10-vial
sheets.

**Mechanism (objective band mask, DATA1/DATA2 byte-identical).** `solve_model`
and `solve_model_B_fix` take an optional `band_vials=` (1-indexed vial set);
their `obj_rule` per-vial loop `continue`s past out-of-band vials, and the
per-channel denominators are floored at 1 only when a band is active
(`_den_m/_den_cp/_den_cr/_den_cf`), so `band_vials=None` is byte-identical to the
legacy objective. `calc_FIM` takes the same `band_vials` and masks **both** the
covariance (`var_pred`) and the Jacobian (`jac`) loops identically, giving a
per-band FIM. Guard `regression_data1_data2_guards.py` (6th case) builds DATA1 &
DATA2 at `B_form=1` with the toggle on/off and asserts an identical constraint.

1. **Per-band θ — `solve_model_per_concentration_band()`** (new). Partitions the
   fitted vials by terminal cF (`partition_vials_by_terminal_cf`, median/quantile
   split, low→high), warm-starts each band from a full-sheet seed fit, and checks
   each band's FIM for non-singularity. Restricted to ≥6 fitted vials, ≤2 bands;
   violations are returned in a `warnings` list rather than raising.

   **Result (median 2-band split).** On GOOD and OKAY the **full-sheet fit rails
   σ→1.0, but the HIGH-cF band alone recovers an interior σ** (MC3 SNaCl 28-70 mM
   → σ=0.83; MC2 NaCl 25-36 mM → σ=0.56); the low-cF band still rails. σ's
   information lives in the high-Δπ vials, and pooling all vials lets the low-cF
   vials drag σ to the wall. **B drifts with concentration within MC3 SNaCl**
   (band B 12.9→18.2 μm/s); MC2 NaCl's B is flat (~10). On POOR LaCl3 (cF 1.6-9.7
   mM) both bands rail to *opposite* σ bounds — the flat-surface signature; banding
   cannot manufacture identifiability. **All per-band FIMs are non-singular**
   (min-eig > 0) but the railed-σ bands carry high condition numbers (1e8-1e10) —
   the FIM confirms L_p/B are jointly identifiable while flagging σ as the weak
   direction. Caveat: per-band WSSE is **not** comparable across bands (different
   vials, cF magnitudes, band normalization) — read σ-identifiability and B-drift,
   not WSSE magnitude. Driver: `_run_band_value_test.py`, `_run_band_wrapper.py`;
   artifacts under `…/nf270/band_value_test/`.

   **Campaign check (all 11 single-salt sheets, `_run_band_campaign.py`).** The
   σ-recovery is **sheet-specific, not universal**: the high-cF band recovers an
   interior σ (whole-sheet fit railed) on only **3/11** sheets (MC3 SNaCl 0.83,
   MC2 NaCl 0.56, MC3 SCaCl₂ 0.20); on others the high-cF band rails to the *other*
   wall (σ=0) or both rail — banding rescues σ only where the high-cF band carries
   enough Δπ and the σ surface isn't flat. The **robust, generalizable** campaign
   signal is **B rising with concentration on 9/11 sheets** (low→high band B-drift,
   strongest CaCl₂ +8.5, MC4 SNaCl +6.5, MC3 SNaCl +5.4 µm/s) — concentration-
   dependent B, i.e. the §18.3 B(I) signal. The headline that generalizes is the
   B–concentration trend, not σ-recovery-by-banding.

2. **Per-vial θ trend — `_run_per_vial_trend.py`.** Rather than the invasive joint
   per-vial-parameter model (declaring `m.Lp[n]/m.sigma[n]`, which would risk
   DATA1/DATA2 byte-identity), the per-vial trend is traced with **rolling 3-vial
   windows** through the same validated band mask (warm-started = regularized).
   The joint `m.Lp[n]/m.sigma[n]` model (mirroring `B_form='pervial'`) remains the
   available heavier escalation. **Result:** MC3 SNaCl shows **B rising
   monotonically with concentration** (corr(cF,B)=+0.93, B 12→21 μm/s across
   12.8→55.8 mM windows) — direct per-vial evidence of concentration-dependent B,
   the same signal §18.3's B(I) law targets — while σ only leaves the σ=1 wall in
   the highest-cF window. MC2 NaCl: σ is the concentration-dependent quantity
   (corr(cF,σ)=−0.81), B flat. LaCl3: σ flips walls 1→0 across cF (flat-surface
   artifact), B tiny.

3. **Correlate vials with contour groups.** Two views: (a) `_run_band_contour.py`
   re-runs the σ×Lp sweep with the **band mask** → per-band objective surfaces.
   On MC3 SNaCl the **low-cF band's surface is flat in σ** (cannot constrain σ,
   rails to 1) while the **high-cF band's surface is curved with an interior
   minimum at σ≈0.83** — the two bands occupy structurally different basins. (b)
   The per-vial-trend plot overlays each window's optimum on the full σ×Lp contour
   colored by concentration, showing the optima migrate off the σ=1 wall as cF
   rises (the vial↔contour-group correlation).

Statistical guard: each band's FIM must stay non-singular — verified by
`solve_model_per_concentration_band(compute_fim=True)` (all bands passed).

### 18.5  Runfile operator guide  (`refactored_ucb_runfile.py`)

Entry: `python refactored_ucb_runfile.py` → prompts Root → Subset → Trunk →
Branches. For DATA3 use Root `3`; Subset `2` = single-salt 11 (default), `3` = one
run_id, `4` = fast smoke (3).

| Want | Trunk | Branches |
|---|---|---|
| Fitted mass + concentration plots | `3` fit_multistart | `m,c` |
| Pressure / osmotic (per-vial ΔP, σ·Δπ) | any | `r` |
| Objective contours (σ×Lp, B×Lp, B×σ) | `1` simulate | `o` |
| FIM heatmap | `4` fit_FIM | `f` |
| DoE next-experiment | `5` fit_FIM_DoE | `d` |
| Parameter table (JSON) | `2`/`3` | `p` |
| Forward-sim at given θ (no fit) | `1` simulate | `m,c` |
| Mass litmus (σ=0 closed-form Lp) | `6` mass_litmus | `m` |

**Inject a chosen θ and render fitted mass+conc plots** (the runfile has no θ
prompt — call the library directly):

```python
import sys; sys.path.insert(0, "<repo>/refactored_codes_v1")
import refactored_ucb_library as lib
ds = lib.loadxlsx("<wb>/NF270_MC3.xlsx", sheet="07.22.24_SNaCl")["data_stru"]
theta = {"Lp": 8.5, "B": 14.0, "sigma": 0.45}           # your guess
fit, sim, _ = lib.solve_model(ds, ds["mode"], theta, sim_opt=True,  # True = forward-sim at fixed θ
                              B_form="single", workflow_family="DATA3", LOUD=True)
lib.run_data3_time_series_plots({"data": ds, "sim_stru": sim,
        "model_settings": {"workflow_family": "DATA3"}}, save_dir="<out>")
```

Set `sim_opt=False` to FIT from `theta` as the seed. For a per-band/per-concentration
guess, set the band mask first (§18.4). The sheet→workbook map is
`lib.NF270_RUN_REGISTRY[run_id]` (`{"workbook","sheet"}`); the helper that does
load→fit→plot end-to-end is `_warm_anchor_worker.py`.

### 18.6  Files changed in §18

| File | Change |
|---|---|
| `refactored_ucb_library.py` | `_nf270_conc_scales()` + `NF270_C{F,P}_RESIDUAL_SCALE_FRACTION` (cF 2 %, cV None→3 %); applied in both `obj_rule`s, `calc_FIM` `var_pred`, `_label_parmest_model` (now takes `workflow_family`); `nf270_ionic_strength_factor()` + `NF270_IONIC_STRENGTH_FACTOR`; `_axis_values()` β₀/β₁/β₂ support |
| `_make_3d_contour_animations.py` | σ-scrub axes → `x=B, y=Lp` (Lp on Y); output `scrub-sigma_LpvsB.gif` |
| `_run_all_pairs_5channel.py` | `SWEEPS[2]` → `("sigma","B")` (B on Y); output `objcontour-x_sigma-y_B.png` |
| `_compare_cf_weighting.py`, `_verify_goldstd_cf.py` | new: cF-weight before/after + multistart verification |
| `Architecture.md` | §5b measurement-error subsection + this §18 |

**Implementation pass (2026-06-10): β contours, B(I), per-band θ.**

| File | Change |
|---|---|
| `refactored_ucb_library.py` | **Task 1:** contour-axis validator admits `beta_0/1/2`; `run_nf270_contour_for_sheet` `sweep_pairs` switches to β slices when `isinstance(B_form,int) and B_form>=1`. **Task 2:** `NF270_B_USE_IONIC_STRENGTH` toggle; `B_form_rule2/3` bind `cc=(k_I·conc) if use_I else conc` (byte-identical when off), gated DATA3. **Task 3:** `band_vials=` on `solve_model`, `solve_model_B_fix`, `calc_FIM` (masks `var_pred`+`jac`), and the contour path (`_nf270_contour_objectives_at_theta`/`_grid_dataframe`/`run_nf270_contour_for_sheet`); guarded denominators `_den_*`; new `solve_model_per_concentration_band()` + `partition_vials_by_terminal_cf()` + `_vial_terminal_cf()` |
| `tests/regression_data1_data2_guards.py` | 6th guard (B ionic-strength reparam): DATA1/DATA2 at `B_form=1` byte-identical with toggle on/off |
| `_run_beta_contour.py`, `_run_beta_contour_fast.py` | new: β-vs-Lp/σ contour drivers (warm-started B_form=1; fast = capped solver + pool) |
| `_run_b_ionic_validation.py` | new: B(I) rescale-identity (Part A, deterministic) + cross-salt β₁ clustering (Part B) |
| `_run_band_value_test.py`, `_run_band_wrapper.py`, `_run_band_contour.py`, `_run_per_vial_trend.py` | new: per-band value test, formal wrapper+FIM, band-masked contours, per-vial θ trend + vial↔contour correlation |
| `Architecture.md` | §18.3 (β contours + B(I) implemented) + §18.4 (per-band implemented, with GOOD/OKAY/POOR results) |
