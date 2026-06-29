# Port #2.5 — DATA1 model-build fix-ups (firstNonNan + n_extra)

> **Operator-approved 2026-06-09 via option (B) after Port #2 surfaced 2 pre-existing DATA1 bugs.**
>
> Parent audit: [`../REPRODUCIBILITY_AUDIT.md`](../REPRODUCIBILITY_AUDIT.md)
> Trigger: [`PORT_sigma_sensitivity.md`](./PORT_sigma_sensitivity.md) found that `solve_model(workflow_family="DATA1")` crashes on both committed DATA1 datasets due to bugs unrelated to Port #2 itself.

---

## 1. Why this exists

`test_data1_regression.py` has had both DATA1 cases SKIPPED since the suite was built, with reason:
> "pre-existing DATA1 model-build bug (con_boundary[2] resolves to False)"

When Port #2 wired `solve_model(workflow_family="DATA1", sim_opt=True)` into a native simulation path, the same root cause surfaced as 2 distinct exceptions on the two committed DATA1 datasets:

| Dataset | mode | Failure | Library line |
|---|---|---|---|
| `dataset501.1.mat` (filtration) | DATA | `TypeError: 'float' object is not iterable` (firstNonNan on scalar) | 2222 |
| `dataset511.12.mat` (diafiltration) | DATA | `KeyError: 'n_extra'` (DATA2-only data_config field) | 2376, 2418, 3068, 3110 |

Both are pre-existing — caused by the DATA1 vs DATA2 data shape divergence:

- **DATA1 filtration** stores `cF_exp` as a **scalar float** per vial (one end-of-vial ICP-OES sample). DATA1 diafiltration stores it as a list (inline probe).
- **DATA1** has **no startup vials** (every vial is real data). **DATA2** has `n_extra` startup vials that aren't fit (modeled separately by the loader).

The library *already* contains the right idiomatic guards for both cases — they just aren't applied everywhere:

- `np.ndim(...) > 0` guard for scalar/list shape mismatch (used for `cV_avg` at lines 2215-2218; missing for `cF_exp` at line 2222)
- `.get('n_extra', 0)` semantic — line 1613 is already guarded by `if mode != 'DATA':` (so DATA1 skips it); the 4 other callsites are not guarded

## 2. The 11 surgical changes (Bugs A + B + C)

> **Scope extended 2026-06-09** during implementation. The initial proposal
> identified Bugs A (firstNonNan scalar) and B (n_extra KeyError, 4 callsites).
> When those landed, the next DATA2-only field surfaced: **Bug C — n_v0
> KeyError at 6 callsites**. Operator approved the Bug C inclusion under the
> same proposal because it follows the EXACT same pattern as Bug B (bare
> `data_stru['data_config']['n_v0']` → `.get('n_v0', 1)`; default `1` matches
> the pre-existing pattern already used at line 1611).
>
> Final count: **1 (Bug A) + 4 (Bug B) + 6 (Bug C) = 11 surgical edits.**



### Bug A — `cF_exp` scalar guard at line 2222

**Current**:
```python
yield m.cF[1,0]==firstNonNan(data_stru['data_raw'][0]['cF_exp'])
```

**Fix** (mirrors the existing pattern at lines 2215-2218):
```python
cf0_raw = data_stru['data_raw'][0]['cF_exp']
if np.ndim(cf0_raw) > 0:
    cF_init = firstNonNan(cf0_raw)
else:
    cF_init = cf0_raw
yield m.cF[1,0] == cF_init
```

### Bug B — `n_extra` defaulted to 0 for DATA1

DATA1 datasets have no `n_extra` field in `data_config`. The library has already documented `"n_extra": 0` as the DATA1-appropriate default at line 4418 (in the dataclass field definition for `MAT_GLOBAL_DEFAULTS`). The fix replaces 4 unguarded direct lookups with `.get('n_extra', 0)`:

- **Line 2376** (`obj_rule` inside `solve_model`):
  - **Before**: `collect_vial = data_stru['data_config']['n'] - data_stru['data_config']['n_extra']`
  - **After**:  `collect_vial = data_stru['data_config']['n'] - data_stru['data_config'].get('n_extra', 0)`
- **Line 2418** (same `obj_rule`, vial-loop conditional):
  - **Before**: `if n_vial > data_stru['data_config']['n_extra']:`
  - **After**:  `if n_vial > data_stru['data_config'].get('n_extra', 0):`
- **Line 3068** (`obj_rule` inside `solve_model_B_fix`):
  - **Before**: `collect_vial = data_stru['data_config']['n'] - data_stru['data_config']['n_extra']`
  - **After**:  same `.get('n_extra', 0)` substitution
- **Line 3110** (same `obj_rule`, vial-loop conditional):
  - **Before**: `if n_vial > data_stru['data_config']['n_extra']:`
  - **After**:  same `.get('n_extra', 0)` substitution

**Line 1613** is NOT touched — it's already guarded by `if mode != 'DATA':`.

### Bug C — `n_v0` defaulted to 1 for DATA1 (added 2026-06-09 during impl)

DATA1 datasets also have no `n_v0` field (verified empirically: DATA1
501.1 / 511.12 `data_config` lacks `n_v0`). Semantically `n_v0` is the
"first vial whose residuals enter the objective"; default value `1` means
"start from vial 1" — every vial counts. Line **1611** already documents
this with the pattern `N_V0 = data_stru['data_config'].get('n_v0', 1)`.

**6 callsites** need the same `.get` substitution:

- **Line 1444** (`model_construct_inter` time-grid setup):
  `if i >= data_stru['data_config']['n_v0']-1:` → `... .get('n_v0', 1)-1`
- **Line 2422** (`obj_rule` inside `solve_model`, mass residual gate)
- **Line 2488** (`obj_rule` inside `solve_model`, cp residual gate)
- **Line 2509** (`obj_rule` inside `solve_model`, perm-probe gate)
- **Line 3118** (`obj_rule` inside `solve_model_B_fix`, mass residual gate)
- **Line 3180** (`obj_rule` inside `solve_model_B_fix`, cp residual gate)

All 6 share the identical shape `if <something> >= data_stru['data_config']['n_v0']`
→ `.get('n_v0', 1)`. Implementation used `replace_all=true` on the bare
lookup substring; the pre-existing `.get('n_v0', 1)` at line 1611 uses
different syntax (`.get()` not `[]`) so it was unaffected.

For DATA1 (`n_v0` missing → default 1): every vial 1..n is included in the
residual sums. ✓ correct (DATA1 has no skipped startup vials).
For DATA2 (`n_v0=3`): vials 1-2 are skipped, vials 3+ included. ✓ unchanged.

## 3. Semantic correctness

For DATA1 (`n_extra` missing → default 0):
- `collect_vial = n - 0 = n` → all vials count toward the objective. ✓ correct (DATA1 has no startup vials)
- `if n_vial > 0` → always true since `n_vial ∈ {1, 2, ..., n}` → all vials enter the residual sum. ✓ correct

For DATA2 (`n_extra=3`):
- `collect_vial = n - 3 = 10` → only the 10 real vials count. ✓ unchanged behavior
- `if n_vial > 3` → vials 1-3 (startup) are skipped, vials 4-13 count. ✓ unchanged behavior

For DATA3/NF270: the NF270 loader (line 4418 onwards) populates `n_extra` explicitly. ✓ unchanged behavior

**Bottom line: the change is byte-equivalent for DATA2 + DATA3 (which already supplied `n_extra`), and unlocks DATA1 (where the field is absent).**

## 4. Acceptance test

A new file `pytest_refactored_codes_v1/tests/test_data1_model_build_fixups.py`:

```python
@pytest.mark.regression
@pytest.mark.parametrize("dataset", [501.1, 511.12])
def test_data1_model_builds_and_simulates(dataset, lib):
    """DATA1 model_construct_inter + solve_model(sim_opt=True) succeed.

    Catches the two pre-existing bugs Port #2.5 fixed:
      Bug A: firstNonNan on scalar cF_exp (501.1 filtration)
      Bug B: KeyError 'n_extra' on DATA1 (no startup-vial concept)
    """
    data_stru = lib.loadmat(f"legacy/data1_matlab/data/data_stru-dataset{dataset}.mat")["data_stru"]
    fit_stru = lib.loadmat(f"legacy/data1_matlab/data/{dataset} concpolar/fit_stru.mat").get("fit_stru")
    theta = lib._data1_contour_theta_from_fit(fit_stru)
    fit, sim_stru, _ = lib.solve_model(
        data_stru, mode="DATA", theta=theta, sim_opt=True,
        B_form="single", workflow_family="DATA1", LOUD=False,
    )
    assert len(sim_stru) == int(data_stru['data_config']['n'])
    for vial in sim_stru:
        for key in ("time", "mV", "cF"):
            assert key in vial


@pytest.mark.smoke
def test_solve_model_n_extra_uses_get_default(lib):
    """Quick guard: source contains the .get('n_extra', 0) pattern, not direct lookups."""
    src = inspect.getsource(lib.solve_model)
    # All 4 lookups inside solve_model + solve_model_B_fix obj_rule should be .get
    assert "data_stru['data_config']['n_extra']" not in src + inspect.getsource(lib.solve_model_B_fix)


@pytest.mark.smoke
def test_first_non_nan_callsite_guarded(lib):
    """Quick guard: the cF_exp callsite at line 2222 wraps firstNonNan with np.ndim check."""
    src = inspect.getsource(lib.solve_model)
    # Look for the pattern that hold the new guard
    assert "np.ndim(cf0_raw)" in src or "ndim(cf0_raw)" in src
```

## 5. What this rewiring unlocks

| Test/figure | Before | After |
|---|---|---|
| `test_data1_regression.py::test_data1_structural_matches_reference[data1__501.11__single]` | SKIPPED ("pre-existing DATA1 model-build bug") | should now run (and if MATLAB structural fingerprint differs, that's a real signal to investigate) |
| `test_data1_regression.py::test_data1_structural_matches_reference[data1__511.12__single]` | SKIPPED | same |
| Port #2 native-vs-MATLAB acceptance tests (6 cases, currently XFAIL) | XFAIL | should PASS (within rtol=5e-3 trajectory tolerance) |
| DATA1 Fig. 4 audit verdict | `PRESENT_WITH_MAT_DEPENDENCY` | **`PRESENT`** (via `regenerate=True`) |

Audit total `PRESENT` count moves: 6 (after Port #1) → **7** of 41 published items.

## 6. What this fix does NOT touch

- Does NOT change `firstNonNan` itself (it's correct for arrays — we only need to guard the scalar call site).
- Does NOT remove `n_extra` semantics from DATA2 or DATA3 — they still work exactly as before.
- Does NOT touch `solve_mass_balance_only`, `model_construct_inter` initial-conditions code (already handles scalar `cV_avg`).
- Does NOT modify any test file beyond ADDING the new `test_data1_model_build_fixups.py`.
- Does NOT touch `conductivity_paper.py`.
- Does NOT touch loader-stage code (DATA1 data_config shape unchanged on disk).

## 7. Risk + mitigation

| Risk | Mitigation |
|---|---|
| `.get('n_extra', 0)` semantically wrong for DATA1 | Cross-checked: DATA1 has no startup vials, so `n_extra=0` means "all vials are real data" — that's correct. Library's own `MAT_GLOBAL_DEFAULTS` at line 4418 documents `"n_extra": 0` as the DATA1 default. |
| `cf0_raw` scalar guard breaks DATA2/DATA3 | DATA2/DATA3 `cF_exp` is always a list; `np.ndim() > 0` is True → original `firstNonNan` path. Byte-equivalent. |
| Existing test_data1_regression cases now FAIL instead of SKIP | If they fail, that's REAL information — the structural fingerprint may have drifted since capture, OR the model build differs from MATLAB. Either way it's a true signal, not a Port #2.5 regression. |

## 8. After the fix lands — follow-up checks

1. Re-run `test_data1_sigma_sensitivity_native_matches_matlab.py` → expect all 6 cases to PASS instead of XFAIL.
2. Re-run `test_data1_regression.py::test_data1_structural_matches_reference` → expect SKIPPED cases to now RUN; report PASS/FAIL.
3. If the structural-fingerprint cases fail, that's a separate diagnosis — report findings to operator before any other change.
4. Update `REPRODUCIBILITY_AUDIT.md` to flip DATA1 Fig. 4 from `PRESENT_WITH_MAT_DEPENDENCY` to `PRESENT`.
