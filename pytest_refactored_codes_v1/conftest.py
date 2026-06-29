"""Session-level conftest for pytest_refactored_codes_v1.

Responsibilities
----------------
1. Put the SUT (refactored_codes_v1/) on sys.path so `import refactored_ucb_library`
   and `import refactored_ucb_runfile` work without modifying the production
   package.
2. Put this folder's `helpers/` on sys.path so tests can do
   `from helpers.reference_cases import ALL_CASES`.
3. Set `MPLCONFIGDIR` if unset, to avoid matplotlib font-cache races.
4. chdir to the repo root, because the library's loaders use relative paths
   like `"data_library/data_stru-...mat"`.
5. Expose session-scoped fixtures: `repo_root`, `lib`, `runfile`, `references`,
   `allowlist`, `session_run_dir`.
6. Rotate `_runs/` so we keep last N session dirs (default 3).

Every cross-test contract lives here. Individual test files import from
`helpers/` and from fixtures defined below.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import time
from pathlib import Path

import pytest


# ----------------------------------------------------------------------
# Path anchors
# ----------------------------------------------------------------------

HERE = Path(__file__).resolve().parent                  # pytest_refactored_codes_v1/
REPO_ROOT = HERE.parent                                 # repo root
LIB_DIR = REPO_ROOT / "refactored_codes_v1"             # the SUT
HELPERS_DIR = HERE / "helpers"
BASELINES_DIR = HERE / "baselines"
RUNS_DIR = HERE / "_runs"

# Put the SUT and helpers on sys.path BEFORE pytest collects any tests.
# Both have to be ahead of any site-packages with the same names.
for p in (str(LIB_DIR), str(HELPERS_DIR), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

# Loaders in the SUT use relative paths — anchor at repo root.
os.chdir(str(REPO_ROOT))

# Isolate the matplotlib font cache so parallel test runs don't race.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")


# ----------------------------------------------------------------------
# Constants exported for tests to import
# ----------------------------------------------------------------------

# Where the .mat data lives — same convention as the runfile.
DATA1_DIR = Path(os.environ.get(
    "DIAFILTRATION_DATA1_ROOT",
    REPO_ROOT / "legacy" / "data1_matlab" / "data"
)).expanduser().resolve()

DATA2_DIR = Path(os.environ.get(
    "DIAFILTRATION_DATA2_ROOT",
    REPO_ROOT / "legacy" / "data1_matlab" / "data_library"
)).expanduser().resolve()

NF270_DIR = Path(os.environ.get(
    "DIAFILTRATION_NF270_ROOT",
    REPO_ROOT / "UnifiedFramework" / "ExperimentalDataFiles"
)).expanduser().resolve()

# Baseline JSONs and CSVs.
REF_JSON = BASELINES_DIR / "regression_references.json"
ALLOWLIST_CSV = BASELINES_DIR / "known_nonpass.csv"
MANIFEST_CSV = BASELINES_DIR / "manifest.csv"

# Same numeric tolerance as the original test suite.
REL_TOLERANCE = 1e-2


# ----------------------------------------------------------------------
# Session housekeeping — rotate _runs/
# ----------------------------------------------------------------------

def _rotate_runs(keep: int = 3) -> None:
    """Keep only the most recent N session dirs under _runs/."""
    if not RUNS_DIR.exists():
        return
    sessions = sorted(
        (d for d in RUNS_DIR.iterdir() if d.is_dir() and d.name.startswith("run-")),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    for d in sessions[keep:]:
        shutil.rmtree(d, ignore_errors=True)


def pytest_sessionstart(session):
    """Rotate _runs/ at the start of every session."""
    keep = int(os.environ.get("PYTEST_REFACTORED_RUNS_RETAIN", "3"))
    _rotate_runs(keep=keep)


# ----------------------------------------------------------------------
# Session-scoped fixtures
# ----------------------------------------------------------------------

@pytest.fixture(scope="session")
def repo_root() -> Path:
    """Repo root (parent of pytest_refactored_codes_v1/)."""
    return REPO_ROOT


@pytest.fixture(scope="session")
def lib():
    """Import refactored_ucb_library once per session.

    Heavy import — Pyomo, CasADi, parmest, matplotlib, NF270 registries — so
    we cache it at session scope. Every test that needs the library imports
    via this fixture instead of doing its own `import refactored_ucb_library`.
    """
    import refactored_ucb_library as _lib
    return _lib


@pytest.fixture(scope="session")
def runfile():
    """Import refactored_ucb_runfile once per session.

    The runfile defines DATA1_FAST_SUBSET, DATA2_FAST_SUBSET, NF270_FAST_SUBSET,
    NF270_SINGLE_SALT_RUNS, TRUNK_RECIPES, DATA1_ROOT, DATA2_ROOT, NF270_ROOT —
    every smoke test needs at least one of these.
    """
    import refactored_ucb_runfile as _rf
    return _rf


@pytest.fixture(scope="session")
def references() -> dict:
    """Load the captured structural fingerprints. Fails loudly if missing."""
    if not REF_JSON.exists():
        pytest.fail(
            f"Reference fingerprints not captured yet at {REF_JSON}.\n"
            f"Run: python3 {HERE / '_capture_references.py'}"
        )
    with open(REF_JSON) as fh:
        return json.load(fh)


@pytest.fixture(scope="session")
def allowlist() -> dict:
    """Parse known_nonpass.csv into {target_id: meta}. Empty dict if absent."""
    from helpers.allowlist import read_exceptions
    if not ALLOWLIST_CSV.exists():
        return {}
    return read_exceptions(ALLOWLIST_CSV)


@pytest.fixture(scope="session")
def session_run_dir() -> Path:
    """Per-session output directory under _runs/.

    Named run-<unixtime>-<pid> so consecutive sessions don't collide. Created
    lazily — fixtures that don't write outputs don't trigger directory creation.
    """
    name = f"run-{int(time.time())}-{os.getpid()}"
    p = RUNS_DIR / name
    p.mkdir(parents=True, exist_ok=True)
    return p


@pytest.fixture(scope="session")
def data_roots(runfile) -> dict:
    """Resolved DATA{1,2,3} roots plus a 'present' flag for each.

    Tests that need a root call `data_roots["DATA1"]["path"]` and
    `data_roots["DATA1"]["present"]` to decide whether to skip.
    """
    return {
        "DATA1": {"path": DATA1_DIR, "present": DATA1_DIR.exists()},
        "DATA2": {"path": DATA2_DIR, "present": DATA2_DIR.exists()},
        "DATA3": {"path": NF270_DIR, "present": NF270_DIR.exists()},
    }


# ----------------------------------------------------------------------
# Public helper — used by both old and new regression tests
# ----------------------------------------------------------------------

def assert_params_match(fitted: dict, reference: dict, tol: float = REL_TOLERANCE):
    """Assert each numeric parameter in `reference` matches `fitted` within
    relative tolerance `tol`. Migrated from refactored_codes_v1/tests/conftest.py.
    """
    mismatches = []
    for key, ref_val in reference.items():
        if key not in fitted:
            mismatches.append(f"  missing key in fitted: {key!r}")
            continue
        if isinstance(ref_val, dict) or isinstance(ref_val, str):
            continue
        if ref_val is None:
            continue
        try:
            ref_f = float(ref_val)
            fit_f = float(fitted[key])
        except (TypeError, ValueError):
            continue
        denom = max(abs(ref_f), 1e-9)
        rel = abs(fit_f - ref_f) / denom
        if rel > tol:
            mismatches.append(
                f"  {key!r}: fitted={fit_f:.6g}, ref={ref_f:.6g}, "
                f"rel_diff={rel:.3%} (tol={tol:.1%})"
            )
    if mismatches:
        pytest.fail("Parameters drifted beyond tolerance:\n" + "\n".join(mismatches))
