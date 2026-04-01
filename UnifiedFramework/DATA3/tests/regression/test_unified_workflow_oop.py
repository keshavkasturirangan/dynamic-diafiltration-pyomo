"""Regression coverage for the object-oriented unified workflow module."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest


os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

REPO_ROOT = Path(__file__).resolve().parents[4]
UNICODE_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "ExperimentalDataAnalysis" / "UnifiedCode"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(UNICODE_ROOT) not in sys.path:
    sys.path.insert(0, str(UNICODE_ROOT))

from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode import unified_codebase_runfile as runfile  # noqa: E402
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library import (  # noqa: E402
    ExperimentMode,
    ModelOptions,
    ParameterTreatmentMode,
    ProcessModelProfile,
)
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.workflow.unified_workflow import (  # noqa: E402
    DataLoader,
    DatasetRequest,
    UQEngine,
)


@pytest.mark.regression
def test_object_oriented_loader_supports_mat_and_xlsx_requests() -> None:
    """The new DataLoader should treat MAT and XLSX as first-class requests."""
    loader = DataLoader()
    requests = [
        DatasetRequest(file_path=runfile.PRESET_DATA1_PATHS[0], dataset_id="DATA1_501_1"),
        DatasetRequest(
            file_path=runfile.PRESET_DATA3_PATH,
            selector=runfile.PRESET_DATA3_SELECTOR,
            dataset_id="DATA3_MC2",
        ),
    ]

    datasets, validations = loader.load_many(requests)

    assert len(datasets) == 2
    assert validations[0][0] is True
    assert validations[1][0] is True
    assert datasets[0].source.value == "mat"
    assert datasets[1].filename.endswith(".xlsx")
    assert datasets[1].sheet_name == runfile.PRESET_DATA3_SELECTOR
    assert len(datasets[1].vials) >= 1


@pytest.mark.regression
def test_uqengine_compare_models_returns_information_criteria() -> None:
    """The OOP UQEngine should expose AIC-style model comparison results."""
    loader = DataLoader()
    dataset, validation = loader.load(
        DatasetRequest(file_path=runfile.PRESET_DATA1_PATHS[0], dataset_id="DATA1_501_1")
    )
    assert validation[0] is True

    engine = UQEngine(data_loader=loader)
    comparison = engine.compare_models(
        [dataset],
        [
            ModelOptions(
                mode=ExperimentMode.DATA,
                parameter_treatment_mode=ParameterTreatmentMode.ESTIMATION,
                b_form="single",
                nfe=10,
                use_sigma_logit_transform=False,
                process_model_profile=ProcessModelProfile.DATA1,
                paper_profile="DATA1_PAPER",
            )
        ],
    )

    assert not comparison.empty
    assert {"candidate_name", "objective", "aic", "aicc", "bic", "delta_aic"}.issubset(comparison.columns)
    assert np.isfinite(float(comparison.loc[0, "objective"]))
