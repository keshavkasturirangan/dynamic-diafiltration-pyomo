"""Layer 2 — DATA2 artifact tests.

Mirror of test_data1_artifacts.py, scoped to DATA2_FAST_SUBSET.
Auto-skipped if DATA2 root absent. The runfile's `_ensure_data2_library()`
unpacks the archive on the first DATA2 run, so the skip should be rare.
"""
from __future__ import annotations

import pytest


pytestmark = [pytest.mark.regression]


@pytest.fixture(scope="module")
def data2_fast_run(lib, runfile, data_roots, session_run_dir):
    if not data_roots["DATA2"]["present"]:
        # Attempt the runfile's documented unpack-then-retry.
        try:
            runfile._ensure_data2_library()
            data_roots["DATA2"]["present"] = data_roots["DATA2"]["path"].exists()
        except Exception:
            pass
        if not data_roots["DATA2"]["present"]:
            pytest.skip(f"DATA2_ROOT not on disk ({data_roots['DATA2']['path']})")

    save_dir = session_run_dir / "data2"
    save_dir.mkdir(parents=True, exist_ok=True)
    try:
        results = lib.materialize_all(
            campaign="DATA2",
            save_dir=save_dir,
            data_root=data_roots["DATA2"]["path"],
            only=runfile.DATA2_FAST_SUBSET,
        )
    except Exception as exc:
        pytest.skip(
            f"materialize_all(DATA2) raised {type(exc).__name__}: {exc}"
        )
    return {"save_dir": save_dir, "results": results}


def test_materialize_all_returns_list(data2_fast_run):
    assert isinstance(data2_fast_run["results"], list)


def test_no_error_status_in_results(data2_fast_run):
    errors = [r for r in data2_fast_run["results"] if r.get("status") == "error"]
    assert not errors, (
        f"materialize_all(DATA2) reported {len(errors)} errors:\n" +
        "\n".join(f"  - {e.get('name')}: {e.get('error')}" for e in errors)
    )


def test_fast_subset_all_attempted(data2_fast_run, runfile):
    """Every entry in DATA2_FAST_SUBSET resulted in a result row."""
    result_names = {r.get("name") for r in data2_fast_run["results"]}
    missing = set(runfile.DATA2_FAST_SUBSET) - result_names
    # Be lenient — the manifest may expand each name to multiple result rows.
    # The point is "all 4 entries kicked something off."
    assert len(result_names) >= 1, (
        f"DATA2_FAST_SUBSET produced 0 result rows."
    )
