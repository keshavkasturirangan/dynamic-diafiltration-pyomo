"""Layer 1 — DATA1 smoke tests.

No solver, no model build, no file I/O beyond reading constants. Verifies
that the public surface of `refactored_ucb_library` and `refactored_ucb_runfile`
exposes the DATA1 contract the runfile expects.

Catches: a public function being deleted/renamed, a constant being moved out
of the library, an environment-variable indirection breaking, DATA1_RUN_REGISTRY
being silently shrunk.
"""
from __future__ import annotations

import pytest


pytestmark = pytest.mark.smoke


# ----------------------------------------------------------------------
# Library import + public surface
# ----------------------------------------------------------------------

# Names the runfile and analysis scripts actually call into. If any of these
# disappear, the active code path breaks. Catching renames here is the whole
# point of the smoke layer.
REQUIRED_PUBLIC_NAMES = (
    "loadmat",
    "loadxlsx",
    "model_construct_inter",
    "solve_model",
    "solve_mass_balance_only",
    "estimate_parameters",
    "calc_FIM",
    "calc_contour_2d_py",
    "calc_ind_objectives_py",
    "calc_ind_objectives_5channel_py",
    "plot_sim_comparison",
    "plot_contour",
    "report_parameters",
    "report_uncertainty",
    "run_pipeline",
    "materialize_all",
    "DATA1_RUN_REGISTRY",
    "DATA2_RUN_REGISTRY",
    "NF270_RUN_REGISTRY",
    "CAMPAIGN_REGISTRY",
)


def test_library_imports_cleanly(lib):
    """The library imports without raising."""
    assert lib is not None


@pytest.mark.parametrize("name", REQUIRED_PUBLIC_NAMES, ids=lambda n: n)
def test_library_public_surface(lib, name):
    """Every name in REQUIRED_PUBLIC_NAMES is present on the library module."""
    assert hasattr(lib, name), (
        f"refactored_ucb_library is missing public name {name!r} — "
        f"a refactor likely deleted or renamed it. Check git log."
    )


# ----------------------------------------------------------------------
# DATA1 runfile contract
# ----------------------------------------------------------------------

def test_runfile_imports_cleanly(runfile):
    """The runfile imports without raising."""
    assert runfile is not None
    assert hasattr(runfile, "main"), "runfile.main() is the documented CLI entry."


def test_data1_fast_subset_constant(runfile):
    """DATA1_FAST_SUBSET is the documented 3-tuple from refactored_ucb_runfile.py."""
    assert runfile.DATA1_FAST_SUBSET == ("figure_2", "figure_3", "data_analysis"), (
        f"DATA1_FAST_SUBSET drifted: now={runfile.DATA1_FAST_SUBSET}"
    )


def test_data1_root_resolved_or_documented(runfile, data_roots):
    """DATA1_ROOT resolves to an existing dir OR the test xfails with a clear
    'data root not on this machine' reason. Either is fine — both prove the
    indirection didn't break, just that the dev box may not have the data."""
    p = runfile.DATA1_ROOT
    if not p.exists():
        pytest.xfail(
            f"DATA1_ROOT does not exist on this machine ({p}). "
            f"Set DIAFILTRATION_DATA1_ROOT to override."
        )
    assert p.is_dir()


# ----------------------------------------------------------------------
# TRUNK_RECIPES contract — defines all named workflows the runfile dispatches
# ----------------------------------------------------------------------

EXPECTED_RECIPE_NAMES = {
    "simulate", "fit", "fit_multistart", "fit_FIM", "fit_FIM_DoE", "mass_litmus",
}

EXPECTED_RECIPE_KEYS = {"multistart", "fim", "doe", "litmus", "sim_only"}


def test_trunk_recipes_complete(runfile):
    """Every named recipe is present and carries the documented bool flags."""
    recipes = runfile.TRUNK_RECIPES
    missing = EXPECTED_RECIPE_NAMES - set(recipes)
    assert not missing, f"TRUNK_RECIPES missing: {missing}"
    for name in EXPECTED_RECIPE_NAMES:
        rec = recipes[name]
        bad_keys = EXPECTED_RECIPE_KEYS - set(rec.keys())
        assert not bad_keys, f"{name} missing keys: {bad_keys}"
        for k in EXPECTED_RECIPE_KEYS:
            assert isinstance(rec[k], bool), (
                f"{name}.{k} should be bool, got {type(rec[k]).__name__}"
            )


def test_trunk_recipe_simulate_is_sim_only(runfile):
    """The 'simulate' recipe has sim_only=True and everything else False."""
    rec = runfile.TRUNK_RECIPES["simulate"]
    assert rec["sim_only"] is True
    assert rec["multistart"] is False
    assert rec["fim"] is False
    assert rec["doe"] is False
    assert rec["litmus"] is False


def test_trunk_recipe_mass_litmus_is_litmus(runfile):
    """The 'mass_litmus' recipe has litmus=True (NF270-specific)."""
    rec = runfile.TRUNK_RECIPES["mass_litmus"]
    assert rec["litmus"] is True
    assert rec["sim_only"] is False
