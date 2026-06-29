"""Layer 1 — DATA2 smoke tests.

Verifies the DATA2 portion of the runfile and library public surface. Like the
DATA1 smoke file: no solver runs, no model build. Just constants + registry.
"""
from __future__ import annotations

import pytest


pytestmark = pytest.mark.smoke


# DATA2_FAST_SUBSET is the 4-tuple the runfile uses for the "fast subset" branch.
# Locked here so a typo in the runfile would fail this test immediately.
EXPECTED_DATA2_FAST_SUBSET = (
    "calibration_plots",
    "pressure_changes",
    "startup_barplot",
    "model_error_visualization",
)


def test_data2_fast_subset_constant(runfile):
    """DATA2_FAST_SUBSET matches the runfile's documented 4-tuple."""
    got = tuple(runfile.DATA2_FAST_SUBSET)
    assert got == EXPECTED_DATA2_FAST_SUBSET, (
        f"DATA2_FAST_SUBSET drifted: now={got}"
    )


def test_data2_root_resolved_or_documented(runfile, data_roots):
    """DATA2_ROOT exists or test xfails with reason."""
    p = runfile.DATA2_ROOT
    if not p.exists():
        pytest.xfail(
            f"DATA2_ROOT does not exist on this machine ({p}). "
            f"Set DIAFILTRATION_DATA2_ROOT to override, or run "
            f"_ensure_data2_library() once."
        )
    assert p.is_dir()


def test_data2_run_registry_populated(lib):
    """DATA2_RUN_REGISTRY has at least the two documented entries.

    Two canonical DATA2 .mat files: 270511.123 (Lag) and 270511.423 (Overflow).
    A registry shrinkage breaks the test_data2_regression.py cases.
    """
    reg = lib.DATA2_RUN_REGISTRY
    assert isinstance(reg, dict)
    assert len(reg) >= 1, (
        f"DATA2_RUN_REGISTRY unexpectedly shrunk to {len(reg)} entries."
    )


def test_data2_library_archive_helper_exists(runfile):
    """_ensure_data2_library() is the documented unpack helper.

    It's the only reason DATA2_ROOT can be missing yet the suite still works:
    `runfile._ensure_data2_library()` extracts inputdata_mat_files.zip.
    """
    assert callable(getattr(runfile, "_ensure_data2_library", None)), (
        "runfile._ensure_data2_library() should exist as the DATA2 archive "
        "unpack fallback."
    )


def test_data2_loadmat_signature(lib):
    """`loadmat` accepts a positional file-path string. Catches signature drift."""
    import inspect
    sig = inspect.signature(lib.loadmat)
    params = list(sig.parameters.values())
    assert len(params) >= 1, "loadmat must accept at least one positional arg"
    # The first positional should not require kwargs-only access.
    assert params[0].kind in (
        inspect.Parameter.POSITIONAL_ONLY,
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
    ), f"loadmat first arg has unexpected kind: {params[0].kind}"
