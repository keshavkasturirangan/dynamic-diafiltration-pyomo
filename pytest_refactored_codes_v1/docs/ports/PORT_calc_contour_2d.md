# Port #1 — `calc_contour_2d.m` → `calc_contour_2d_py`

> **Status: the Python port is already done. This proposal is for the REWIRING that connects DATA1 figure workflows to the already-existing Python function instead of reading MATLAB-precomputed CSVs.**
>
> Parent audit: [`../REPRODUCIBILITY_AUDIT.md`](../REPRODUCIBILITY_AUDIT.md)
> Operator approval gate: NO code change lands until this proposal is approved.

---

## 1. Discovery (was not in the audit at first pass)

`refactored_codes_v1/refactored_ucb_library.py` already contains:

- **`calc_contour_2d_py`** (line 13943) — explicit Python port of `legacy/data1_matlab/functions/calc_contour_2d.m`. Its docstring cites MATLAB lines 42-50 (the Lp special rule), 52-53 (the linspace meshgrid), 65-83 (the inner loop), and 100-107 (the CSV write).
- **`plot_contour_py`** (line 14077) — Python port of `plot_contour.m` (the renderer that emits the 14×7-inch 3-panel `objcontour-x_{x_var}-y_{y_var}.png` figure).
- **`calc_ind_objectives_py`** (line 13744) and **`calc_ind_objectives_5channel_py`** (line 13801) — Python equivalents of MATLAB `calc_ind_objectives.m`.

So the math + renderer + CSV writer are already in Python. **No new function is needed.** The verdict in the audit (`PRESENT_WITH_MAT_DEPENDENCY`) is caused by the WORKFLOW wrappers — `run_data1_figure5_workflow` and `run_data1_figure6_workflow` — reading pre-computed MATLAB CSVs instead of calling `calc_contour_2d_py` to recompute them.

## 2. What the MATLAB does (math + I/O)

`legacy/data1_matlab/functions/calc_contour_2d.m` (104 lines):

- **Input:** `data_stru`, `model_stru`, `theta`, `x_var`, `x_ind`, `y_var`, `y_ind`, `grid_density`.
- **Bounds:** read `model_stru.theta_lb/theta_ub` for each variable. Lp special rule overrides to `[0.1·θ_Lp, 2·θ_Lp]` (lines 42-50).
- **Grid:** `meshgrid(linspace(x_lb, x_ub, grid_density), linspace(y_lb, y_ub, grid_density))`.
- **Inner loop:** for every (i) of `grid_density²` cells, set `beta(x_ind)=xx(i)`, `beta(y_ind)=yy(i)`, call `calc_ind_objectives(0, beta, data_stru, model_stru)` (scale_opt=0 = no scaling), take `log10` of the unscaled objectives (mass, perm_conc, reten_conc).
- **Output:** writes `[model_stru.filenamestr]/contourdata-x_{x_var}-y_{y_var}.csv` with columns `{x_var, y_var, Obj_mass, Obj_concentration, Obj_retentate_concentration}`.

Default `grid_density` used by `run_data_analysis.m` is **50**, producing **2,500 cells = 2,501 CSV rows including header**. Verified by `wc -l` on existing CSVs.

## 3. What `calc_contour_2d_py` already does

`refactored_codes_v1/refactored_ucb_library.py:13943`. Excerpt of signature:

```python
def calc_contour_2d_py(
    data_stru, theta, x_var, y_var, *,
    grid_density=50,
    save_dir,
    x_bounds=None, y_bounds=None,
    mode="Lag", B_form="single", workflow_family="DATA3",
    nfe=80,
    write_csv=True,
    include_conductivity=False,
) -> pd.DataFrame
```

Behavior:
- Same Lp special rule (`x_lb=0.1*Lp`, `x_ub=2*Lp` when `var=="Lp"`).
- Same `np.linspace + np.meshgrid` grid generation.
- Calls `calc_ind_objectives_py` (or `calc_ind_objectives_5channel_py` if `include_conductivity=True`) — Python analogue of MATLAB `calc_ind_objectives.m`.
- Writes `contourdata-x_{x_var}-y_{y_var}.csv` with the EXACT same 5-column schema and filename convention as MATLAB.
- Default `grid_density=50` matches MATLAB.

Default `workflow_family="DATA3"`. For DATA1 reproductions the caller must pass `workflow_family="DATA1"`.

## 4. What's actually broken — the wrappers don't call the port

`run_data1_figure5_workflow` (line 7267) currently does:

```python
cases = [
    ("A", root / "511.12 concpolar" / "contourdata-x_sigma-y_Lp.csv", "figure5_panel_a"),
    ("B", root / "511.11 concpolar" / "contourdata-x_sigma-y_Lp.csv", "figure5_panel_b"),
    ("C", root / "511.12"            / "contourdata-x_sigma-y_Lp.csv", "figure5_panel_c"),
]
for _, csv_path, prefix in cases:
    if not csv_path.exists():
        continue
    df = pd.read_csv(csv_path)
    panel_paths = _plot_contour_data1_legacy(df, prefix, save_dir, ...)
```

It reads CSVs that were produced by the MATLAB `run_data_analysis.m` script *months ago* and committed to `legacy/data1_matlab/data/*/contourdata-*.csv`. Deleting those CSVs breaks the workflow.

The same pattern applies to:
- `run_data1_figure6_workflow` (Lp × B contours, line 7332)
- `run_data1_si_s3` / `si_s4` / `si_s5` / `si_s6` (4-model SI variants, lines 6871-6912)

All 6 figures of the audit's "PRESENT_WITH_MAT_DEPENDENCY" rows for contour figures share this same blocker.

## 5. Proposed rewiring (CAREFULLY SCOPED)

For EACH of the 6 figure-workflow functions, replace the "read CSV" with a "compute CSV (or read if cached)" branch. Keep the existing MATLAB CSVs as the FALLBACK so byte-equivalence is provable, and add the Python-native generation as the PRIMARY path:

```python
def run_data1_figure5_workflow(data_root=None, save_dir=None, *, regenerate=False):
    root = _resolve_data1_figure2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data1_paper_figures"
    save_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        # (panel_label, dataset_id, variant, prefix)
        ("A", "511.12", "concpolar", "figure5_panel_a"),
        ("B", "511.11", "concpolar", "figure5_panel_b"),
        ("C", "511.12", None,        "figure5_panel_c"),
    ]
    outputs, panel_groups = [], []
    for label, ds, variant, prefix in cases:
        # NEW: prefer native generation, fall back to cached MATLAB CSV
        df = _data1_contour_df(
            data_root=root,
            dataset_id=ds, variant=variant,
            x_var="sigma", y_var="Lp",
            grid_density=50,
            regenerate=regenerate,
        )
        if df is None:
            continue
        panel_paths = _plot_contour_data1_legacy(df, prefix, save_dir,
                                                   show_title=False, preface=False)
        outputs.extend(panel_paths)
        panel_groups.append([Path(p) for p in panel_paths])
    # ... existing composite logic unchanged ...
```

with new helper:

```python
def _data1_contour_df(*, data_root, dataset_id, variant,
                      x_var, y_var, grid_density=50, regenerate=False) -> pd.DataFrame | None:
    """Return the contour DataFrame for a DATA1 (dataset, variant, x, y) case.

    Behavior:
      1. If regenerate=False AND the legacy MATLAB CSV exists, read it. (byte-identical to old behavior)
      2. Else, load the dataset .mat, load the warm-start theta from the fit_stru.mat, and call
         calc_contour_2d_py(..., workflow_family="DATA1", mode="DATA", B_form="single",
                             grid_density=grid_density).
      3. If the Pyomo build raises, return None (let the caller skip the panel) — same as
         the current "if not csv_path.exists(): continue".
    """
    variant_dir = f"{dataset_id} {variant}" if variant else dataset_id
    csv_path = Path(data_root) / variant_dir / f"contourdata-x_{x_var}-y_{y_var}.csv"
    if not regenerate and csv_path.exists():
        return pd.read_csv(csv_path)

    # Native path
    mat_path = Path(data_root) / f"data_stru-dataset{dataset_id}.mat"
    fit_path = Path(data_root) / variant_dir / "fit_stru.mat"
    if not mat_path.exists() or not fit_path.exists():
        return None
    data_stru = loadmat(str(mat_path))["data_stru"]
    fit_stru = loadmat(str(fit_path)).get("fit_stru", {})
    theta = _data1_contour_theta_from_fit(fit_stru)  # already exists, line ~7397
    try:
        return calc_contour_2d_py(
            data_stru, theta, x_var, y_var,
            grid_density=grid_density,
            save_dir=csv_path.parent,
            workflow_family="DATA1", mode="DATA", B_form="single",
            write_csv=True,
        )
    except Exception:
        return None
```

`_data1_contour_theta_from_fit` already exists in `refactored_ucb_library.py:7397` (I read it during the audit). It already extracts Lp, B, σ, β₀, β₁, S₀, S from a legacy `fit_stru` dict.

## 6. Acceptance test (the gate for landing this)

The proposed rewiring **does not change current behavior** when the cached MATLAB CSV is present (`regenerate=False` default = read CSV). Behavior only changes when the CSV is absent OR `regenerate=True`.

Acceptance test (added to `pytest_refactored_codes_v1/tests/`):

```python
@pytest.mark.regression
@pytest.mark.parametrize("dataset_id,variant", [
    ("501.1", "concpolar"), ("501.11", "concpolar"),
    ("511.11", "concpolar"), ("511.12", "concpolar"),
    ("501.1", None), ("501.11", None), ("511.11", None), ("511.12", None),
])
@pytest.mark.parametrize("xy", [("sigma","Lp"), ("B","Lp")])
def test_data1_contour_native_matches_matlab(dataset_id, variant, xy, lib, tmp_path):
    """Python-native contour computation matches the committed MATLAB CSV
    column-for-column within tolerance."""
    x_var, y_var = xy
    variant_dir = f"{dataset_id} {variant}" if variant else dataset_id
    csv_matlab = LEGACY_DATA1_ROOT / variant_dir / f"contourdata-x_{x_var}-y_{y_var}.csv"
    if not csv_matlab.exists():
        pytest.skip(f"MATLAB reference CSV not present: {csv_matlab.name}")
    df_matlab = pd.read_csv(csv_matlab)
    df_native = lib._data1_contour_df(
        data_root=LEGACY_DATA1_ROOT,
        dataset_id=dataset_id, variant=variant,
        x_var=x_var, y_var=y_var,
        grid_density=50,
        regenerate=True,  # force native
    )
    assert df_native is not None, "native generation failed (model build raised?)"
    assert set(df_native.columns) == set(df_matlab.columns)
    # Coordinate columns: exact match (linspace is deterministic)
    np.testing.assert_allclose(df_native[x_var], df_matlab[x_var], rtol=1e-12)
    np.testing.assert_allclose(df_native[y_var], df_matlab[y_var], rtol=1e-12)
    # Objective columns: tighter than visual-similarity threshold but loose enough
    # to absorb IPOPT float-noise. log10-SSR values typically range 0.5-3.0.
    for col in ("Obj_mass", "Obj_concentration", "Obj_retentate_concentration"):
        # Drop rows where either side is NaN (Pyomo build may fail on boundary cells)
        mask = df_matlab[col].notna() & df_native[col].notna()
        np.testing.assert_allclose(
            df_native.loc[mask, col], df_matlab.loc[mask, col],
            rtol=1e-3, atol=1e-3,
            err_msg=f"{col} drift between native and MATLAB",
        )
```

**16 parametrized cases** (8 dataset/variant × 2 axis pairs). Wall-clock estimate: 50×50 = 2,500 cells per case × ~1.5 s per Pyomo build = ~62 min per case sequential, OR ~6 min per case with 10 workers via `multiprocessing.Pool`. Total at 10 workers: ~100 min for all 16 cases. **Tag `nightly`**, not `regression`.

A LIGHTER smoke variant runs `grid_density=8` (64 cells) per case → ~6 sec each → 16 × 6 = ~100 sec total. Tag this `regression`. The full 50×50 stays `nightly`.

## 7. What changes in `refactored_ucb_library.py`

| Function | Change |
|---|---|
| `run_data1_figure5_workflow` (L7267) | Replace `cases` list of CSV paths with `(label, dataset, variant, prefix)` tuples + call new helper. Add optional `regenerate=False` kwarg. |
| `run_data1_figure6_workflow` (L7332) | Same pattern, for B × Lp. |
| `run_data1_si_s3` (L6871) | Same pattern (4-model variant, dataset list expands to include `511.11`). |
| `run_data1_si_s4` (L6882) | Same. |
| `run_data1_si_s5` (L6893) | Same, FILTRATION cases (501.1 / 501.11). |
| `run_data1_si_s6` (L6904) | Same. |
| (new) `_data1_contour_df(...)` | New private helper, ~30 LOC, near the existing `_data1_contour_theta_from_fit` (L7397). |

**Estimated diff:** ~80 lines net add, ~30 lines net delete. All within `refactored_ucb_library.py`. No new public API.

## 8. Risks + mitigations

| Risk | Mitigation |
|---|---|
| Python-native output drifts from MATLAB reference | Default `regenerate=False` reads cached CSV → bit-identical old behavior. Native path only fires when CSV missing or `regenerate=True`. |
| `calc_ind_objectives_py` produces different log10(SSR) than MATLAB `calc_ind_objectives.m` because of IPOPT float-noise on the inner DAE solve | Test tolerance `rtol=1e-3, atol=1e-3` absorbs typical noise. If a real divergence shows, that's a genuine model regression — flag, don't silence. |
| Loader `loadmat` returns slightly different `data_stru` shape from MATLAB | `data_stru` is already loaded the same way in `run_data1_si_s2` (L6697); reuse that path. |
| `_data1_contour_theta_from_fit` returns θ in wrong unit/scale | Existing function is already used by `run_data1_direct_contour_branch` and has a tested θ extraction. |
| Native path takes 50× longer than CSV read | Why we keep CSV read as default. `regenerate=True` is opt-in. |

## 9. What this rewiring removes from the audit

After this lands, the audit verdicts shift:

| Item | Before | After |
|---|---|---|
| DATA1 Fig. 5 | `PRESENT_WITH_MAT_DEPENDENCY` | `PRESENT` (when `regenerate=True`) |
| DATA1 Fig. 6 | `PRESENT_WITH_MAT_DEPENDENCY` | `PRESENT` |
| DATA1 Fig. S3 | `PRESENT_WITH_MAT_DEPENDENCY` | `PRESENT` |
| DATA1 Fig. S4 | `PRESENT_WITH_MAT_DEPENDENCY` | `PRESENT` |
| DATA1 Fig. S5 | `PRESENT_WITH_MAT_DEPENDENCY` | `PRESENT` |
| DATA1 Fig. S6 | `PRESENT_WITH_MAT_DEPENDENCY` | `PRESENT` |

**6 of 21 `PRESENT_WITH_MAT_DEPENDENCY` items resolved** → total `PRESENT` count goes from 0 to 6.

## 10. What this rewiring does NOT do

- Does NOT touch `calc_contour_2d_py` or `plot_contour_py` (already correct).
- Does NOT touch `calc_ind_objectives_py` (the inner objective evaluator).
- Does NOT modify the existing MATLAB CSVs in `legacy/data1_matlab/data/`.
- Does NOT remove the legacy `run_data_analysis.m` MATLAB code path.
- Does NOT change `materialize_all` dispatch logic.
- Does NOT touch DATA2.
- Does NOT touch `conductivity_paper.py`.

## 11. Operator decision needed

Approve any subset:

- [ ] **(a)** Approve the rewiring of `run_data1_figure5_workflow` (Lp × σ, diafiltration main paper Fig. 5)
- [ ] **(b)** Approve `run_data1_figure6_workflow` (Lp × B, diafiltration main paper Fig. 6)
- [ ] **(c)** Approve `run_data1_si_s3` + `s4` (diafiltration SI 4-model variants)
- [ ] **(d)** Approve `run_data1_si_s5` + `s6` (filtration SI 4-model variants)
- [ ] **(e)** Approve the acceptance test additions to `pytest_refactored_codes_v1/tests/`

Recommended order: **(a) + (b) first** (smallest scope, main paper figures, 8 acceptance cases). Validate. Then **(c) + (d)** (SI figures, 8 more cases). Then **(e)** lands the regression tests.

If you'd like a different acceptance-test tolerance or want to keep `regenerate=False` as the default (so the suite gains the option but doesn't enforce native generation unless asked), say so before I write the code.
