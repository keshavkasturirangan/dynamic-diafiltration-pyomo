"""Layer 1 — DATA3 / NF270 smoke tests.

The NF270 campaign is the live experimental work that drives the advisor deck
and the 3D contour analyses. These smoke checks ensure the runfile constants,
the library registry, and the NF270 knob constants are all in sync.

No XLSX loaded (would be too slow for the smoke layer); just registry + constants.
"""
from __future__ import annotations

import pytest


pytestmark = pytest.mark.smoke


# ----------------------------------------------------------------------
# Single-salt subset — the 11 sheets the active campaign uses
# ----------------------------------------------------------------------

EXPECTED_SINGLE_SALT_RUNS = (
    # NaCl (6 sheets)
    "MC2.05.07.24_NaCl",
    "MC3.07.22.24_SNaCl",
    "MC4.07.11.24_SNaCl",
    "MC5.07.23.24_NaCl",
    "MC5.07.23.24_SNaCl",
    "MC5.07.23.24_S2NaCl",
    # CaCl2 (3 sheets)
    "MC2.05.07.24_CaCl2",
    "MC3.07.11.24_SCaCl2",
    "MC3.07.12.24_S2CaCl2",
    # LaCl3 (2 sheets)
    "MC2.05.21.24_LaCl3",
    "MC4.07.11.24_SLaCl3",
)


def test_nf270_single_salt_constant(runfile):
    """NF270_SINGLE_SALT_RUNS matches the documented 11-sheet tuple."""
    got = tuple(runfile.NF270_SINGLE_SALT_RUNS)
    assert got == EXPECTED_SINGLE_SALT_RUNS, (
        f"NF270_SINGLE_SALT_RUNS drifted:\n"
        f"  now: {got}\n"
        f"  expected: {EXPECTED_SINGLE_SALT_RUNS}"
    )


def test_nf270_fast_subset_constant(runfile):
    """NF270_FAST_SUBSET is the 3-tuple smoke-test entry."""
    got = tuple(runfile.NF270_FAST_SUBSET)
    assert got == (
        "MC2.05.08.24_E1",
        "MC4.07.12.24_E13",
        "tables.parameters_table",
    ), f"NF270_FAST_SUBSET drifted: {got}"


# ----------------------------------------------------------------------
# NF270_RUN_REGISTRY — every documented sheet must be a registered run
# ----------------------------------------------------------------------

@pytest.mark.parametrize("sheet", EXPECTED_SINGLE_SALT_RUNS, ids=lambda s: s)
def test_nf270_single_salt_in_registry(lib, sheet):
    """Each single-salt sheet is a registered run in NF270_RUN_REGISTRY."""
    assert sheet in lib.NF270_RUN_REGISTRY, (
        f"{sheet!r} dropped from NF270_RUN_REGISTRY."
    )


def test_nf270_registry_size(lib):
    """NF270_RUN_REGISTRY has at least the 11 single-salt sheets + green-flagged extras."""
    n = len(lib.NF270_RUN_REGISTRY)
    assert n >= 11, f"NF270_RUN_REGISTRY shrunk to {n} entries (expected >= 11)."


def test_nf270_root_resolved_or_documented(runfile):
    """NF270_ROOT exists or test xfails with reason."""
    p = runfile.NF270_ROOT
    if not p.exists():
        pytest.xfail(
            f"NF270_ROOT does not exist on this machine ({p}). "
            f"Set DIAFILTRATION_NF270_ROOT to the directory with NF270_MC*.xlsx."
        )
    assert p.is_dir()


# ----------------------------------------------------------------------
# Per-salt B bound contract — must halve per +1 valence
# ----------------------------------------------------------------------

def test_nf270_b_bounds_per_salt_shape(lib):
    """NF270_B_BOUNDS_PER_SALT has the documented NaCl/CaCl2/LaCl3 entries with
    upper bounds 30, 15, 10 respectively (halve per +1 valence)."""
    bounds = lib.NF270_B_BOUNDS_PER_SALT
    assert set(bounds.keys()) >= {"NaCl", "CaCl2", "LaCl3"}, (
        f"NF270_B_BOUNDS_PER_SALT missing salts: {set(bounds.keys())}"
    )
    # Upper bounds — these are the live values the active campaign uses.
    assert bounds["NaCl"][1] == 30.0, f"NaCl upper bound drifted: {bounds['NaCl']}"
    assert bounds["CaCl2"][1] == 15.0, f"CaCl2 upper bound drifted: {bounds['CaCl2']}"
    assert bounds["LaCl3"][1] == 10.0, f"LaCl3 upper bound drifted: {bounds['LaCl3']}"


# ----------------------------------------------------------------------
# Branches set — pick_branches accepts the documented 8 letters
# ----------------------------------------------------------------------

def test_pick_branches_documented_letters(runfile):
    """pick_branches() emits a set drawn from {m, c, r, p, f, d, o, l}."""
    # We don't actually call pick_branches() (it does input()), but we can
    # inspect the literal 'all' set in the function source via the documented
    # behavior: BRANCHES exposes 8 letters m/c/r/p/f/d/o/l. Catching one being
    # dropped is what this test exists for. Read the source as a sanity check.
    import inspect
    src = inspect.getsource(runfile.pick_branches)
    for letter in ("m", "c", "r", "p", "f", "d", "o", "l"):
        # Either present as a quoted single letter in the 'valid' set or in
        # the 'all' return.
        assert f'"{letter}"' in src, (
            f"pick_branches() no longer references branch letter {letter!r}."
        )
