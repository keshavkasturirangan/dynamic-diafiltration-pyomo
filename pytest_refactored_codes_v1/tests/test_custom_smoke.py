"""Layer 1 — Custom dispatcher smoke tests.

The Custom path in `refactored_ucb_runfile._dispatch_custom()` is the
single-file invocation: the user gives a `.mat` or `.xlsx` path and the runfile
routes it through `loadmat` / `loadxlsx` + `run_pipeline` accordingly.

These tests assert the dispatcher's exposed surface — without actually running
solvers — so a refactor that breaks the Custom contract surfaces here.
"""
from __future__ import annotations

import inspect

import pytest


pytestmark = pytest.mark.smoke


def test_dispatch_custom_function_exists(runfile):
    """`_dispatch_custom` is the documented Custom-root handler."""
    assert callable(getattr(runfile, "_dispatch_custom", None)), (
        "runfile._dispatch_custom is the public Custom-file dispatcher."
    )


def test_dispatch_custom_uses_run_pipeline(runfile, lib):
    """_dispatch_custom calls `ucb.run_pipeline` somewhere in its body."""
    src = inspect.getsource(runfile._dispatch_custom)
    assert "ucb.run_pipeline" in src, (
        "_dispatch_custom() no longer routes through ucb.run_pipeline()."
    )


def test_dispatch_custom_handles_xlsx_branch(runfile):
    """The xlsx branch detection + DATA3-style time-series plotting is still wired."""
    src = inspect.getsource(runfile._dispatch_custom)
    assert ".xlsx" in src.lower() or "xls" in src.lower(), (
        "_dispatch_custom() no longer detects .xlsx files."
    )
    assert "DATA3" in src, (
        "_dispatch_custom() no longer routes xlsx through DATA3 workflow_family."
    )


def test_dispatch_custom_litmus_path(runfile):
    """The mass-litmus closed-form branch in Custom dispatch is still present."""
    src = inspect.getsource(runfile._dispatch_custom)
    assert "solve_mass_balance_only" in src, (
        "_dispatch_custom() no longer dispatches to solve_mass_balance_only "
        "for the litmus recipe."
    )


def test_run_pipeline_signature(lib):
    """`run_pipeline` accepts the documented keyword set.

    Catches a refactor that renames `mode`, `workflow_family`, `B_form`, etc.
    Customizing _dispatch_custom would silently break without this assertion.
    """
    sig = inspect.signature(lib.run_pipeline)
    documented = {
        "mode", "workflow_family", "B_form", "selector",
        "use_parmest", "uncertainty_method", "nfe",
        "multistart", "multistart_iterations",
    }
    missing = documented - set(sig.parameters.keys())
    assert not missing, (
        f"run_pipeline() lost documented kwargs: {missing}\n"
        f"current signature: {list(sig.parameters)}"
    )


def test_loadxlsx_signature(lib):
    """`loadxlsx` accepts a `sheet=` selector. Catches signature drift."""
    sig = inspect.signature(lib.loadxlsx)
    assert "sheet" in sig.parameters, (
        f"loadxlsx() lost the documented `sheet` kwarg. "
        f"signature: {list(sig.parameters)}"
    )
