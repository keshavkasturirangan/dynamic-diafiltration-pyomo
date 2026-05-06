from __future__ import annotations

from pathlib import Path

import pytest

import refactored_codes_v1.refactored_ucb_library as ucb


NF270_ROOT = Path(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/ExperimentalDataFiles"
)


@pytest.mark.slow
def test_nf270_single_salt_smoke_materialize_and_coverage(tmp_path):
    save_dir = tmp_path / "nf270_smoke"
    results = ucb.materialize_all(
        campaign="NF270",
        save_dir=save_dir,
        data_root=NF270_ROOT,
        only=("MC2.05.07.24_NaCl",),
        extra_opts={
            "request_overrides": {
                "multistart": False,
                "multistart_iterations": 1,
                "uncertainty_method": "",
            }
        },
    )

    assert results, "Expected a non-empty NF270 smoke result list"
    assert results[0]["status"] == "ok", results[0]
    assert len(results[0]["paths"]) == 2

    report = ucb.report_paper_coverage("NF270", save_dir)
    assert report["campaign"] == "NF270"
    assert report["summary"]["main_figures"][0] >= 2
    assert report["summary"]["main_figures"][1] >= 26
    assert "parameters_table" in ucb.lookup_figure("NF270", "tables.parameters_table").name


def test_nf270_registration_and_table_contract():
    assert "NF270" in ucb.list_campaigns()
    assert len(ucb.list_figures("NF270")) == 27
    table_spec = ucb.lookup_figure("NF270", "tables.parameters_table")
    assert table_spec.name == "nf270.tables.parameters_table"

