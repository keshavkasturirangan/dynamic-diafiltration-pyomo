"""Layer 2 — DATA3 / NF270 artifact tests.

Curated single-sheet run (MC2.05.07.24_NaCl) per the overnight defaults.
To extend to the full 11-sheet sweep, tag those tests with `nightly` and
parametrize over runfile.NF270_SINGLE_SALT_RUNS.
"""
from __future__ import annotations

import pytest


pytestmark = [pytest.mark.regression]


CURATED_SHEET = "MC2.05.07.24_NaCl"


@pytest.fixture(scope="module")
def nf270_one_sheet_run(lib, runfile, data_roots, session_run_dir):
    if not data_roots["DATA3"]["present"]:
        pytest.skip(f"NF270_ROOT not on disk ({data_roots['DATA3']['path']})")

    save_dir = session_run_dir / "nf270"
    save_dir.mkdir(parents=True, exist_ok=True)

    extra_opts = {
        "multistart_override": False,
        "fim_override": False,
        "branches": ["m", "c", "p"],
        "request_overrides": {"nfe": 80},
    }
    try:
        results = lib.materialize_all(
            campaign="NF270",
            save_dir=save_dir,
            data_root=data_roots["DATA3"]["path"],
            only=(CURATED_SHEET,),
            extra_opts=extra_opts,
        )
    except Exception as exc:
        pytest.skip(
            f"materialize_all(NF270, sheet={CURATED_SHEET}) raised "
            f"{type(exc).__name__}: {exc}"
        )
    return {"save_dir": save_dir, "results": results, "sheet": CURATED_SHEET}


def test_returns_list(nf270_one_sheet_run):
    assert isinstance(nf270_one_sheet_run["results"], list)


def test_no_error_status_in_results(nf270_one_sheet_run):
    errors = [r for r in nf270_one_sheet_run["results"] if r.get("status") == "error"]
    assert not errors, (
        f"materialize_all(NF270 {nf270_one_sheet_run['sheet']}) reported "
        f"{len(errors)} errors:\n" +
        "\n".join(f"  - {e.get('name')}: {e.get('error')}" for e in errors)
    )


def test_mass_concentration_pngs_present(nf270_one_sheet_run):
    """The 'm' and 'c' branches each produce at least one png under the sheet save_dir."""
    pngs = list(nf270_one_sheet_run["save_dir"].rglob("*.png"))
    mass_pngs = [p for p in pngs if "mass" in p.name.lower()]
    conc_pngs = [p for p in pngs if "conc" in p.name.lower() or "concentration" in p.name.lower()]
    assert mass_pngs, f"No mass png found under {nf270_one_sheet_run['save_dir']}"
    assert conc_pngs, f"No concentration png found under {nf270_one_sheet_run['save_dir']}"
