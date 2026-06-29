# `pytest_refactored_codes_v1/` — Pytest Suite for the Refactored Diafiltration Codebase

> Built overnight on 2026-06-08 → 2026-06-09 as a self-contained test harness for
> `refactored_codes_v1/refactored_ucb_library.py` and `refactored_codes_v1/refactored_ucb_runfile.py`.
> Every test file, fixture input, captured baseline, allowlist, and per-session
> run-output lives inside this folder. No edits are made to the production
> library or to the existing `refactored_codes_v1/tests/` tree (that tree is
> COPY-migrated here; see §11 for the migration history).

---

## 1. What this suite does and why

The refactored library is the active code path for DATA1 (KCl closed-cell paper
reproductions), DATA2 (NF270 / KCl Lag + Overflow paper reproductions), and
DATA3 (NF270 single-salt experimental campaign — the live work that drives the
advisor deck and the 3D contour analyses). Every commit to
`refactored_ucb_library.py` or `refactored_ucb_runfile.py` risks silent drift in
one of three places:

| Drift surface | Symptom |
|---|---|
| Loader (`loadmat`, `loadxlsx`) | `data_stru` shape changes → DAE constraints reference wrong indices → fits silently shift |
| Model build (`model_construct_inter`) | Pyomo Var/Constraint set changes → fitted θ moves by 1–5 % with no apparent cause |
| Dispatcher (`materialize_all`, `run_pipeline`) | Subset / branch / nfe defaults change → published-paper figures stop reproducing byte-for-byte |

This suite catches all three with a four-layer architecture (§3) that runs in
under 30 seconds for the every-commit smoke layer, and under ~10 minutes for the
full regression layer (one cached forward simulation per DATA1/DATA2 reference
case, no IPOPT fits).

The goal is not perfect coverage of every code path — it's **fast, deterministic
detection of structural drift on the four entry-point families** documented by
`refactored_ucb_runfile.py`: DATA1 / DATA2 / DATA3 / Custom (plus the
`DATA2_workflow_for_DATA3` compatibility bridge, which is currently a non-goal
but reachable through the manifest).

---

## 2. Folder layout

```
pytest_refactored_codes_v1/
├── README.md                            ← you are here
├── pytest.ini                           ← markers + addopts + testpaths
├── conftest.py                          ← sys.path injection, session fixtures
├── _capture_references.py               ← refresh baselines (migrated)
├── .gitignore                           ← keeps _runs/ out of git
├── helpers/
│   ├── __init__.py
│   ├── allowlist.py                     ← read_exceptions, evaluate (failures/warnings/resolved)
│   ├── reference_cases.py               ← migrated catalog of DATA1/DATA2 cases
│   └── trajectory_compare.py            ← migrated structural-fingerprint comparator
├── baselines/
│   ├── regression_references.json       ← canonical structural fingerprints (migrated)
│   ├── known_nonpass.csv                ← allowlist (warning-only by default)
│   └── manifest.csv                     ← Layer 4 parametrization driver
├── tests/
│   ├── __init__.py
│   ├── test_data1_smoke.py              ← Layer 1 — DATA1 wiring
│   ├── test_data2_smoke.py              ← Layer 1 — DATA2 wiring
│   ├── test_data3_smoke.py              ← Layer 1 — DATA3 / NF270 wiring
│   ├── test_custom_smoke.py             ← Layer 1 — Custom dispatcher
│   ├── test_data1_artifacts.py          ← Layer 2 (regression marker, skipped if DATA roots absent)
│   ├── test_data2_artifacts.py          ← Layer 2
│   ├── test_data3_artifacts.py          ← Layer 2
│   ├── test_data1_regression.py         ← Layer 4 (migrated)
│   ├── test_data2_regression.py         ← Layer 4 (migrated)
│   ├── test_data3_guards.py             ← Layer 4 — guard inertness (converted from print-script)
│   └── test_validation_gate.py          ← Layer 3
├── _runs/                               ← OUTPUTS — gitignored; per-session tmp dirs land here
│   └── .gitkeep
└── docs/                                ← supplementary docs (architecture diagrams, etc.)
```

---

## 3. Four-layer architecture

Inspired by `pytest_starter_template.md` in your `Downloads/` archive, mapped onto
the active codebase:

### Layer 1 — Smoke (`smoke` marker)

Imports, constants, type checks, registry presence. No solver, no model build,
no file I/O beyond reading text on disk. Target: every commit, < 5 seconds total.

What it catches: a function being deleted or renamed, a constant being moved,
a registry entry being silently dropped, an environment-variable indirection
breaking.

### Layer 2 — Artifact (`regression` marker)

Runs one cached `materialize_all` invocation per entry-point family at the
smallest subset that exercises the dispatch (DATA1_FAST_SUBSET, DATA2_FAST_SUBSET,
one NF270 sheet). Then inspects the produced file tree: are the expected PNGs
present, the parameter JSON written, the side-by-side CSV emitted, etc.

What it catches: a dispatcher change that drops an output, a save-path
regression, a `branches` selection bug that skips a figure family.

These tests are auto-skipped if the underlying data roots
(`DIAFILTRATION_DATA1_ROOT` etc.) are not on disk — so the suite stays green on
a developer machine without the full data library.

### Layer 3 — Validation gate (`regression` marker)

A helper (`helpers/allowlist.py`) parses `baselines/known_nonpass.csv` into a
`{target_id: meta}` map, filters by `active=true` and `review_after >= today`,
and returns a `(failures, warnings, resolved)` triple from any side-by-side
DataFrame produced by Layer 2.

- `failures` — rows with `status == "FAIL"` that are not in the allowlist
- `warnings` — rows with `status == "FAIL"` that ARE allowlisted (expected_status matches)
- `resolved` — rows with `status == "PASS"` that are still in the allowlist (stale)

Default mode is **warning-only** (matches your archived `nightly_ops.md` semantics).
Switch to fail-hard by setting `VALIDATION_GATE_FAIL_HARD=1` in the environment
(see §9).

### Layer 4 — Parametrized regression (`regression` + `nightly` markers)

Two complementary halves:

**(a) Structural fingerprint regression** — migrated from
`refactored_codes_v1/tests/`. For each reference case
(`helpers/reference_cases.py`), build the Pyomo model via `model_construct_inter`
WITHOUT solving it, and assert every loader fingerprint scalar + Pyomo Var name
+ constraint name matches the captured baseline in
`baselines/regression_references.json`. Runs in ~1 second per case.

**(b) Numeric regression** (extensible) — parametrized over `baselines/manifest.csv`.
Per row: load a baseline CSV, run the corresponding library invocation through a
cached fixture, compute the named metric (`nrmse`, `mae`, `rmse`, `param_rmse`),
assert it is within the row's `tolerance`. Empty tolerance → `pytest.mark.xfail`.

What this catches: silent shifts in `n_atomic_vars` (a new Var leaks into the
build), `n_atomic_cons` (a constraint is duplicated), loader scalar changes
(an Excel cell shifts), or numeric drift in committed CSV baselines.

---

## 4. How baselines are captured

The structural-fingerprint baselines in `baselines/regression_references.json` were
captured ONCE with the committed code state, by running:

```bash
cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo
python3 pytest_refactored_codes_v1/_capture_references.py
```

This script (migrated from `refactored_codes_v1/tests/_capture_references.py`,
paths re-pointed at the new home) iterates every case in
`helpers/reference_cases.py`, builds the Pyomo model via `model_construct_inter`
without solving, records:

- Loader fingerprint: numeric `data_config` scalars + per-vial array lengths
- Build fingerprint: Pyomo Var/Constraint counts + names (atomic + component-level)

…and writes `baselines/regression_references.json`.

**When to re-capture:**

| Situation | Action |
|---|---|
| Deliberate physics change (new ODE, new constraint) | Re-capture, review the JSON diff, commit both changes together |
| Loader fix that corrects a real data shape bug | Re-capture, commit with a note explaining what scalar changed |
| Random IPOPT noise causing fit-level drift | NOT applicable — these tests don't depend on the solver. If a fingerprint moves, something structural changed. |

**When NOT to re-capture:**

- Tests are failing after a refactor → first understand why. The capture script
  doesn't know whether your change was intentional or a bug.
- Want to silence a noisy test → use the allowlist (§5) instead. Re-capturing
  to silence a real regression destroys the historical record.

---

## 5. Allowlist semantics (`baselines/known_nonpass.csv`)

Schema:

```csv
target_id,expected_status,active,review_after,reason
DATA1.fig2_A.vial3,FAIL,true,2026-09-01,known mass-trace drift; see commit abc123
DATA2.calib_plots,NOT_APPLICABLE,true,2027-01-01,upstream calibration CSV pending
```

| Column | Meaning |
|---|---|
| `target_id` | Stable identifier matching the side-by-side CSV's `target_id` column |
| `expected_status` | One of `FAIL`, `NOT_APPLICABLE`, `PASS_WITH_EXPLANATION` |
| `active` | `true` to honor the entry, `false` to soft-delete |
| `review_after` | ISO date; entry is auto-ignored if `today > review_after` (forces cleanup) |
| `reason` | One-line human-readable justification; ALWAYS link the commit or issue |

**Workflow:**

1. A test starts failing on a row the team has agreed to accept (e.g. a Fig 2
   panel that's known to be off by 1.5 % because the wet-lab calibration
   shifted). Add a row with `active=true` and a near-term `review_after`.
2. The test passes again (as a warning), but a `WARN` line shows up in pytest output.
3. When the underlying issue is fixed, the test starts emitting `resolved`
   warnings; the team removes the stale row.
4. When `today > review_after`, the entry is ignored and the test fails again —
   forcing a re-evaluation. Defaults: `review_after = capture_date + 90 days`.

**Default mode:** warning-only (`failures` array asserted empty;
`warnings` and `resolved` are logged but don't fail). To switch to fail-hard:

```bash
VALIDATION_GATE_FAIL_HARD=1 pytest -m regression pytest_refactored_codes_v1/
```

---

## 6. Markers

Defined in `pytest.ini`:

| Marker | Use for | Run with |
|---|---|---|
| `smoke` | Imports, constants, registry checks; < 5 s total | `-m smoke` |
| `regression` | One cached run per family + structural fingerprints; ~1 min | `-m "regression and not nightly"` |
| `nightly` | Full numeric regression, multi-sheet NF270 sweeps; ~30 min | `-m nightly` |
| `slow` | Anything > 60 s wall-clock | usually stacked with nightly |

`addopts = -ra --strict-markers` — unknown markers fail collection rather than
silently being treated as expressions.

---

## 7. Exact commands

```bash
# Fast smoke (every commit; runs even without the .mat data on disk):
MPLCONFIGDIR=/tmp/mplconfig pytest -q -m smoke pytest_refactored_codes_v1/

# Regression layer (one cached run per family — needs DATA roots on disk):
MPLCONFIGDIR=/tmp/mplconfig pytest -q -m "regression and not nightly" pytest_refactored_codes_v1/

# Nightly layer (heavier — full NF270 single-salt sweep):
MPLCONFIGDIR=/tmp/mplconfig pytest -q -m nightly pytest_refactored_codes_v1/

# Whole suite:
MPLCONFIGDIR=/tmp/mplconfig pytest -q pytest_refactored_codes_v1/

# Single test while iterating:
MPLCONFIGDIR=/tmp/mplconfig pytest -q pytest_refactored_codes_v1/tests/test_data2_regression.py -k Lag_with_time_correction

# Force fail-hard validation gate:
VALIDATION_GATE_FAIL_HARD=1 MPLCONFIGDIR=/tmp/mplconfig pytest -q -m regression pytest_refactored_codes_v1/

# Refresh baselines after a deliberate physics change:
python3 pytest_refactored_codes_v1/_capture_references.py
```

`MPLCONFIGDIR=/tmp/mplconfig` is the convention from your archived `nightly_ops.md`
— it isolates matplotlib's font cache from `~/.matplotlib/` so parallel test runs
don't fight over the same lockfile.

---

## 8. Automation surface

This suite is designed to slot into existing automation in three ways:

### 8.1 Every-commit pre-push hook (recommended)

```bash
# In .git/hooks/pre-push:
MPLCONFIGDIR=/tmp/mplconfig pytest -q -m smoke pytest_refactored_codes_v1/ || exit 1
```

Smoke layer is fast enough (< 5 s) to run on every push without slowing down the
developer loop. Catches the 90 % case (a public function was deleted, a constant
was moved, a registry entry was dropped).

### 8.2 Local nightly cron

```bash
# crontab:
0 2 * * *  cd ~/GitHub/keshav-dev-dynamic-diafiltration-pyomo && \
           MPLCONFIGDIR=/tmp/mplconfig pytest -q -m "regression and not nightly" \
             pytest_refactored_codes_v1/ \
             --junitxml=pytest_refactored_codes_v1/_runs/latest_junit.xml \
             >> pytest_refactored_codes_v1/_runs/nightly.log 2>&1
```

Writes a JUnit XML for trend analysis and appends to a rolling log. Pair with
`caffeinate -i -m` if the cron runs while the laptop sleeps.

### 8.3 GitHub Actions (when you want CI back)

The archived `.github/workflows/nightly-validation.yml` is retired — when you
re-enable a CI workflow, the entry-point is `pytest -m smoke
pytest_refactored_codes_v1/`. No new workflow file is created by this overnight
build (per §12 non-goals); you author the workflow whenever you choose to
re-enable CI.

---

## 9. Environment variables that affect this suite

| Variable | Default | Effect |
|---|---|---|
| `DIAFILTRATION_DATA1_ROOT` | `legacy/data1_matlab/data` | Where DATA1 `.mat` files live (smoke skips if absent) |
| `DIAFILTRATION_DATA2_ROOT` | `legacy/data1_matlab/data_library` | DATA2 `.mat` location |
| `DIAFILTRATION_NF270_ROOT` | `UnifiedFramework/ExperimentalDataFiles` | DATA3 `.xlsx` location |
| `MPLCONFIGDIR` | (unset) | `MPLCONFIGDIR=/tmp/mplconfig` recommended to avoid font-cache races |
| `VALIDATION_GATE_FAIL_HARD` | `0` | `1` makes Layer 3 fail on unexpected FAILs instead of warning |
| `PYTEST_REFACTORED_RUNS_RETAIN` | `3` | How many session `_runs/` subdirs to keep before rotating |

These match the conventions in `refactored_ucb_runfile.py` so the same env vars
work for both running the production code and running the tests.

---

## 10. Adding a new reference case

1. Edit `helpers/reference_cases.py`, append a dict to `DATA1_CASES` or `DATA2_CASES`:

   ```python
   {
       "case_id": "data2__270511.121__truncated_DATA",
       "data_file": "data_stru-dataset270511.121.mat",
       "data_dir":  "legacy/data1_matlab/data_library",
       "mode":      "DATA",
       "B_form":    1,
       "theta":     {"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01,
                     "sigma": 1.0, "S0": 0, "S": 0.0},
       "workflow_family": "DATA2",
   },
   ```

2. Re-run the capture script:

   ```bash
   python3 pytest_refactored_codes_v1/_capture_references.py
   ```

3. The new pytest case is auto-picked-up via `@pytest.mark.parametrize` —
   no test-file edit needed.

4. Commit `helpers/reference_cases.py` and `baselines/regression_references.json`
   together so the reference and the case definition are in lockstep.

---

## 11. Migration history (what came from where)

This suite was built on 2026-06-08 → 2026-06-09 by COPY-migrating the existing
`refactored_codes_v1/tests/` scaffold into the new sibling folder. **The
original files were not deleted** — they remain in place as a safety net. You
can compare and delete the old tree once you're satisfied this suite covers
everything.

| Old path (still on disk) | New path |
|---|---|
| `refactored_codes_v1/tests/conftest.py` | Split into `pytest_refactored_codes_v1/conftest.py` + `helpers/__init__.py` (paths re-rooted) |
| `refactored_codes_v1/tests/reference_cases.py` | `pytest_refactored_codes_v1/helpers/reference_cases.py` |
| `refactored_codes_v1/tests/_capture_references.py` | `pytest_refactored_codes_v1/_capture_references.py` (paths re-rooted) |
| `refactored_codes_v1/tests/_trajectory_compare.py` | `pytest_refactored_codes_v1/helpers/trajectory_compare.py` |
| `refactored_codes_v1/tests/regression_data1_data2_guards.py` | Converted from print-script to pytest module: `tests/test_data3_guards.py` |
| `refactored_codes_v1/tests/regression_references.json` | `pytest_refactored_codes_v1/baselines/regression_references.json` |
| `refactored_codes_v1/tests/test_data1_regression.py` | `pytest_refactored_codes_v1/tests/test_data1_regression.py` |
| `refactored_codes_v1/tests/test_data2_regression.py` | `pytest_refactored_codes_v1/tests/test_data2_regression.py` |

The new test files (Layer 1 smoke + Layer 2 artifacts + Layer 3 validation gate)
have NO predecessor — they were written fresh for this suite.

---

## 12. Non-goals (what this suite deliberately doesn't do)

- **No edits to `refactored_codes_v1/`.** The library is the system under test.
  Tests that surface real bugs report them; they don't paper over.
- **No edits to `UnifiedFramework/.../UnifiedCode/`** (the V24 stack). Out of
  scope per the architectural decision in this conversation.
- **No new baselines generated by the test suite at runtime.** Baselines are
  refreshed only via the explicit `_capture_references.py` script (§4).
- **No CI workflow re-enabled.** The retired `nightly-validation.yml` stays
  retired. You author CI re-enablement when you choose to.
- **No solver-level fit-result regression** at the every-commit layer. Structural
  fingerprints are what catch the changes you actually care about; full IPOPT
  fits are too flaky to gate on at commit time.
- **No `DATA2_workflow_for_DATA3` coverage** at first; can be added in a
  follow-up.

---

## 13. Defaults chosen overnight (subject to your review tomorrow)

The plan I presented before bed asked you four open questions. I picked these
defaults overnight; flip any of them with a one-line change.

| Question | Default chosen | Where to override |
|---|---|---|
| 1. Migration strategy | COPY (old tree intact) | Delete `refactored_codes_v1/tests/` after audit |
| 2. Validation gate failure mode | Warning-only | `VALIDATION_GATE_FAIL_HARD=1` env var |
| 3. DATA3 regression scope | Single sheet `MC2.05.07.24_NaCl` | Add sheets to `helpers/reference_cases.py:NF270_CASES` |
| 4. `_runs/` retention | Keep last 3 | `PYTEST_REFACTORED_RUNS_RETAIN=N` env var |

---

## 14. What runs when you wake up

I deliberately did NOT trigger any solver-touching pytest layers overnight,
because `bj211xovt` (the 3-sheet 3D contour remainder) is consuming 8 IPOPT
workers and the regression tests would either fight for cores or wait for the
GIL. The smoke layer was the only test execution; it verifies imports and
constants, no model build, no solve.

In the morning, to fully exercise the suite:

```bash
# After bj211xovt completes (check task #11):
MPLCONFIGDIR=/tmp/mplconfig pytest -q pytest_refactored_codes_v1/
```

This runs all four layers. Expected wall-clock: smoke ~5 s, structural
regression ~5 s, Layer 2 artifacts ~3–10 min depending on which DATA roots
are configured, validation gate ~1 s. Total < 15 min on a quiet machine.

---

## 15. One-line index for the impatient

```
pytest -q -m smoke pytest_refactored_codes_v1/   # 5 seconds, runs on every push
```

Everything else is documented above. Questions or pushback: edit this file or
delete a section — the suite has no implicit dependency on the prose here.
