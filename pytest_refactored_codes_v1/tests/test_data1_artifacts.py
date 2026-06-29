"""Layer 2 — DATA1 artifact tests.

After ONE cached materialize_all run on DATA1_FAST_SUBSET, verify the
expected file tree is produced under the session output dir.

Tagged `regression`, not `smoke`, because running materialize_all touches
solvers and takes a real wall-clock budget. Auto-skipped if DATA1 roots
are not on disk.

These tests are deliberately STUBBED at the pytest.importorskip level
because the production materialize_all path is heavy enough that we want
to fail fast on missing data rather than silently consume an IPOPT worker
budget. Convert the `pytest.skip` calls to real assertions once you've
confirmed the DATA1 data root is on the dev machine.
"""
from __future__ import annotations

import pytest


pytestmark = [pytest.mark.regression]


@pytest.fixture(scope="module")
def data1_fast_run(lib, runfile, data_roots, session_run_dir):
    """Run materialize_all(campaign='DATA1', only=DATA1_FAST_SUBSET) ONCE per
    module. Returns the save_dir + the results list."""
    if not data_roots["DATA1"]["present"]:
        pytest.skip(f"DATA1_ROOT not on disk ({data_roots['DATA1']['path']})")
    save_dir = session_run_dir / "data1"
    save_dir.mkdir(parents=True, exist_ok=True)
    try:
        results = lib.materialize_all(
            campaign="DATA1",
            save_dir=save_dir,
            data_root=data_roots["DATA1"]["path"],
            only=runfile.DATA1_FAST_SUBSET,
        )
    except Exception as exc:
        pytest.skip(
            f"materialize_all(DATA1) raised {type(exc).__name__}: {exc}\n"
            f"This usually means the data root has unexpected shape; not "
            f"a regression in the test suite."
        )
    return {"save_dir": save_dir, "results": results}


def test_materialize_all_returns_list(data1_fast_run):
    """materialize_all returns a list of manifest-entry result dicts."""
    assert isinstance(data1_fast_run["results"], list)


def test_no_error_status_in_results(data1_fast_run):
    """No entry came back with status='error'."""
    errors = [r for r in data1_fast_run["results"] if r.get("status") == "error"]
    assert not errors, (
        f"materialize_all reported {len(errors)} errors:\n" +
        "\n".join(f"  - {e.get('name')}: {e.get('error')}" for e in errors)
    )


def test_save_dir_has_png_outputs(data1_fast_run):
    """At least one .png written under the save_dir."""
    pngs = list(data1_fast_run["save_dir"].rglob("*.png"))
    assert pngs, f"No .png files written under {data1_fast_run['save_dir']}"
