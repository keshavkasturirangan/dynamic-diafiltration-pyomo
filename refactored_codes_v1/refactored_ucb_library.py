"""
Library of functions for diafiltration experiment modeling
Xinhong Liu
University of Notre Dame
"""

import numpy as np
import pandas as pd
import scipy.io as spio
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.patheffects as patheffects
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import idaes
import time
import copy
import os
import json
import subprocess
import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union
from sklearn.metrics import r2_score

from pyomo.environ import *
from pyomo.dae import *
import idaes.core.util.scaling as iscale


class SourceType(str, Enum):
    """What kind of file the loader is reading."""

    MAT = "mat"
    XLSX = "xlsx"


class WorkflowFamily(str, Enum):
    """Which model family to use."""

    DATA1 = "DATA1"
    DATA2 = "DATA2"
    CUSTOM = "CUSTOM"


@dataclass(frozen=True)
class ExperimentalSource:
    """One input file and its small bit of metadata."""

    path: Path
    selector: Optional[Union[str, int]] = None
    label: Optional[str] = None


@dataclass
class ExperimentalBundle:
    """A group of loaded experiments and where they came from."""

    sources: List[ExperimentalSource] = field(default_factory=list)
    experiments: List[Dict[str, Any]] = field(default_factory=list)
    validation: List[Tuple[bool, List[Tuple[str, str]]]] = field(default_factory=list)
    source_types: List[Optional[SourceType]] = field(default_factory=list)
    profile: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.experiments)


@dataclass(frozen=True)
class ModelOptions:
    """Simple settings that control one run."""

    model_family: WorkflowFamily = WorkflowFamily.DATA1
    backend_name: str = "auto"
    model_builder_name: Optional[str] = None
    mode: str = 'DATA'
    theta: Optional[Dict[str, Any]] = None
    B_form: Union[int, str] = 1
    sim_opt: bool = False
    sigma_fixed: bool = True
    LOUD: bool = False
    log_transform_Pe: bool = False
    plot_pred: bool = True
    cond: bool = True
    lg: bool = False
    preface: bool = False
    output_dir: Optional[Union[str, Path]] = None
    formula: str = 'backward'
    step: float = 1e-8
    profile: Optional[str] = None


@dataclass
class DataLoader:
    """Load files and do small input conversions."""

    repo_root: Path = field(default_factory=lambda: Path.cwd())
    data_root: Optional[Union[str, Path]] = None

    def __post_init__(self):
        self.repo_root = Path(self.repo_root).expanduser().resolve()
        self.data_root = Path(self.data_root).expanduser().resolve() if self.data_root is not None else self.repo_root

    def resolve(self, relative_path: Union[str, Path]) -> Path:
        path = Path(relative_path).expanduser()
        if path.is_absolute():
            return path.resolve()
        return (self.data_root / path).resolve()

    def load_mat(self, relative_path: Union[str, Path], key: str = 'data_stru') -> Dict[str, Any]:
        data = loadmat(str(self.resolve(relative_path)))
        return data[key] if key in data else data

    def load_csv(self, relative_path: Union[str, Path], **kwargs) -> pd.DataFrame:
        return pd.read_csv(self.resolve(relative_path), **kwargs)

    def convert_conductivity_frame(
        self,
        frame: pd.DataFrame,
        *,
        conductivity_column: str = "Conductivity",
        output_column: str = "Concentration",
        temp_K: float = 298.15,
        model: str = "variant_shedlovsky",
        model_params: Optional[Dict[str, object]] = None,
        output_units: str = "mM",
        salt_name: Optional[str] = None,
    ) -> pd.DataFrame:
        """Convert conductivity to concentration using the paper model."""
        converted = frame.copy()
        if conductivity_column not in converted.columns:
            raise KeyError(f"Missing conductivity column: {conductivity_column!r}")
        # Call the shared DATA3 converter directly so we avoid duplicating the
        # conductivity physics in this refactored layer.
        from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library import (
            _conductivity_to_concentration_series,
        )

        conc = _conductivity_to_concentration_series(
            cond_uS_cm=converted[conductivity_column].to_numpy(dtype=float),
            temp_K=temp_K,
            model=model,
            model_params=dict(model_params or {}),
            output_units=output_units,
            salt_name=salt_name,
        )
        converted[output_column] = conc
        converted.attrs["conductivity_model"] = model
        converted.attrs["conductivity_output_units"] = output_units
        return converted

    def load_bundle(self, file_paths: Sequence[Union[str, Path]], selectors: Optional[Sequence[Optional[Union[str, int]]]] = None, profile: Optional[str] = None) -> ExperimentalBundle:
        return load_experiment_bundle([self.resolve(fp) for fp in file_paths], selectors=selectors, profile=profile)

    @staticmethod
    def data1_registry() -> Dict[str, List[str]]:
        return {
            'core_mat': [
                'data/data_stru-dataset501.1.mat',
                'data/data_stru-dataset501.11.mat',
                'data/data_stru-dataset511.11.mat',
                'data/data_stru-dataset511.12.mat',
            ],
            'experimental_space_csv': [
                'data/experiment space/Classical_analysis-dat301.1.csv',
                'data/experiment space/Classical_analysis-dat501.1.csv',
                'data/experiment space/conductivity_calibration.csv',
                'data/experiment space/diafiltration.csv',
                'data/experiment space/filtration.csv',
            ],
            'contour_data': [
                'data/501.1/contourdata-x_B-y_Lp.csv',
                'data/501.1/contourdata-x_sigma-y_Lp.csv',
                'data/501.1/fit_stru.mat',
                'data/501.1 concpolar/contourdata-x_B-y_Lp.csv',
                'data/501.1 concpolar/contourdata-x_sigma-y_Lp.csv',
                'data/501.1 concpolar/fit_stru.mat',
                'data/501.1 concpolar/contour_sig_stru-dat501.1.mat',
                'data/501.11/contourdata-x_B-y_Lp.csv',
                'data/501.11/contourdata-x_sigma-y_Lp.csv',
                'data/501.11/fit_stru.mat',
                'data/501.11 concpolar/contourdata-x_B-y_Lp.csv',
                'data/501.11 concpolar/contourdata-x_sigma-y_Lp.csv',
                'data/501.11 concpolar/fit_stru.mat',
                'data/511.11/contourdata-x_B-y_Lp.csv',
                'data/511.11/contourdata-x_sigma-y_Lp.csv',
                'data/511.11/fit_stru.mat',
                'data/511.11 concpolar/contourdata-x_B-y_Lp.csv',
                'data/511.11 concpolar/contourdata-x_sigma-y_Lp.csv',
                'data/511.11 concpolar/fit_stru.mat',
                'data/511.12/contourdata-x_B-y_Lp.csv',
                'data/511.12/contourdata-x_sigma-y_Lp.csv',
                'data/511.12/fit_stru.mat',
                'data/511.12 concpolar/contourdata-x_B-y_Lp.csv',
                'data/511.12 concpolar/contourdata-x_sigma-y_Lp.csv',
                'data/511.12 concpolar/fit_stru.mat',
                'data/511.12 concpolar/contour_sig_stru-dat511.12-endvial1.mat',
                'data/511.12 concpolar/contour_sig_stru-dat511.12-endvial5.mat',
                'data/511.12 concpolar/contour_sig_stru-dat511.12-endvial10.mat',
            ],
            'special_case': [
                'data/dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0/data_stru-dataset301.1.mat',
                'data/dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0/fit_stru-dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0.mat',
            ],
            'sigma_sensitivity': [
                'data/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.1.mat',
                'data/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.5.mat',
                'data/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.9.mat',
                'data/sigma sensitivity/sim_stru-dat511.12 C_Fin15.2052sig0.1.mat',
                'data/sigma sensitivity/sim_stru-dat511.12 C_Fin15.2052sig0.5.mat',
                'data/sigma sensitivity/sim_stru-dat511.12 C_Fin15.2052sig0.9.mat',
            ],
        }


class DiafiltrationExperiment:
    """One experiment, its settings, and its results."""

    def __init__(
        self,
        data_loader: DataLoader,
        data_source: Union[str, Path, Dict[str, Any]],
        options: Optional[ModelOptions] = None,
        *,
        label: Optional[str] = None,
    ):
        self.data_loader = data_loader
        self.data_source = data_source
        self.options = options or ModelOptions()
        self.label = label
        self.data_stru: Optional[Dict[str, Any]] = None
        self.fit_stru: Optional[Dict[str, Any]] = None
        self.sim_stru: Optional[Any] = None
        self.sim_inter: Optional[Any] = None
        self.fim_stru: Optional[Dict[str, Any]] = None

    @property
    def dataset(self) -> Optional[str]:
        if self.data_stru is None:
            return None
        return str(self.data_stru.get('dataset', self.label or 'unknown'))

    def load(self) -> Dict[str, Any]:
        if isinstance(self.data_source, dict):
            self.data_stru = self.data_source
        else:
            loaded = self.data_loader.load_mat(self.data_source)
            self.data_stru = loaded['data_stru'] if isinstance(loaded, dict) and 'data_stru' in loaded else loaded
        return self.data_stru

    def fit(self) -> Tuple[Dict[str, Any], Any, Any]:
        if self.data_stru is None:
            self.load()
        model_builder = resolve_model_builder(self.options.model_family, self.options.model_builder_name or self.options.backend_name)
        self.fit_stru, self.sim_stru, self.sim_inter = solve_experiment(self.data_stru, self.options, model_builder=model_builder)
        return self.fit_stru, self.sim_stru, self.sim_inter

    def calc_fim(self, **overrides) -> Dict[str, Any]:
        if self.data_stru is None:
            self.load()
        model_builder = resolve_model_builder(self.options.model_family, self.options.model_builder_name or self.options.backend_name)
        params = {
            'theta': self.options.theta,
            'step': self.options.step,
            'formula': self.options.formula,
            'B_form': self.options.B_form,
        }
        params.update(overrides)
        self.fim_stru = calc_FIM(
            self.data_stru,
            self.options.mode,
            theta=params['theta'],
            step=params['step'],
            formula=params['formula'],
            B_form=params['B_form'],
            model_builder=model_builder,
        )
        return self.fim_stru

    def plot(self, output_dir: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
        if self.data_stru is None or self.sim_stru is None:
            self.fit()
        plot_dir = output_dir or self.options.output_dir
        return plot_sim_comparison(
            self.data_stru,
            self.sim_stru,
            plot_pred=self.options.plot_pred,
            cond=self.options.cond,
            lg=self.options.lg,
            preface=self.options.preface,
            output_dir=plot_dir,
        )

    def summarize(self) -> Dict[str, Any]:
        return {
            'dataset': self.dataset,
            'label': self.label,
            'mode': self.options.mode,
            'has_fit': self.fit_stru is not None,
            'has_fim': self.fim_stru is not None,
        }


class UQEngine:
    """Small helper that runs experiments for the user."""

    def __init__(self, loader: DataLoader):
        self.loader = loader

    def build_experiment(
        self,
        data_source: Union[str, Path, Dict[str, Any]],
        options: Optional[ModelOptions] = None,
        *,
        label: Optional[str] = None,
    ) -> DiafiltrationExperiment:
        return DiafiltrationExperiment(self.loader, data_source, options=options, label=label)

    def fit_and_plot(self, experiment: DiafiltrationExperiment) -> Dict[str, Any]:
        experiment.load()
        experiment.fit()
        plot_result = experiment.plot()
        return {'fit': experiment.fit_stru, 'plot': plot_result, 'summary': experiment.summarize()}

    def compute_fim(self, experiment: DiafiltrationExperiment, **overrides) -> Dict[str, Any]:
        return experiment.calc_fim(**overrides)

    def run_data1_bundle(
        self,
        data_sources: Sequence[Union[str, Path, Dict[str, Any]]],
        options: Optional[ModelOptions] = None,
    ) -> List[Dict[str, Any]]:
        results = []
        for idx, src in enumerate(data_sources):
            exp = self.build_experiment(src, options=options, label=f'data1-{idx}')
            exp.load()
            exp.fit()
            results.append({'experiment': exp.summarize(), 'plot': exp.plot()})
        return results


PAPER_REPRODUCTION_TARGETS: Dict[str, Dict[str, Any]] = {
    "DATA1": {
        "script": "UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py",
        # This is the canonical DATA1 paper reproduction target used by the
        # unified validation script.
        "datasets": ["DATA1_511.12"],
    },
    "DATA2": {
        "script": "UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py",
        # These are the canonical DATA2 paper reproduction targets used by the
        # unified validation script.
        "datasets": ["DATA2_270511.123", "DATA2_270611.123"],
    },
}


def run_paper_reproduction(
    paper_family: Union[str, WorkflowFamily],
    *,
    repo_root: Optional[Union[str, Path]] = None,
    output_root: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Run the exact paper reproduction shortcut for DATA1 or DATA2.

    The helper keeps the runtime simple: the user only chooses the paper
    family, and this function knows which datasets the validation script uses.
    """
    family_name = paper_family.value if isinstance(paper_family, WorkflowFamily) else str(paper_family)
    if family_name not in PAPER_REPRODUCTION_TARGETS:
        raise ValueError(f"Unknown paper family: {family_name!r}")

    root = Path(repo_root).expanduser().resolve() if repo_root is not None else Path.cwd().resolve()
    target = PAPER_REPRODUCTION_TARGETS[family_name]
    script = root / target["script"]
    out_root = Path(output_root).expanduser().resolve() if output_root is not None else (root / "UnifiedFramework" / "DATA3" / "results" / "reproduction")

    cmd = [
        sys.executable,
        str(script),
        "--repo-root",
        str(root),
        "--out-root",
        str(out_root),
        "--datasets",
        *target["datasets"],
    ]
    subprocess.run(cmd, check=True, cwd=root)
    return {
        "paper_family": family_name,
        "script": str(script),
        "output_root": str(out_root),
        "datasets": list(target["datasets"]),
    }


_MODEL_BUILDER_REGISTRY: Dict[Tuple[WorkflowFamily, str], Any] = {}


def register_model_builder(family: Union[WorkflowFamily, str], builder_name: str, builder_fn) -> None:
    """Register a model builder for a workflow family."""
    family_enum = family if isinstance(family, WorkflowFamily) else WorkflowFamily(family)
    _MODEL_BUILDER_REGISTRY[(family_enum, builder_name.lower())] = builder_fn


def resolve_model_builder(family: Union[WorkflowFamily, str], builder_name: Optional[str] = None):
    """Resolve the requested model builder for a workflow family."""
    family_enum = family if isinstance(family, WorkflowFamily) else WorkflowFamily(family)
    chosen_name = (builder_name or "auto").lower()
    if (family_enum, chosen_name) in _MODEL_BUILDER_REGISTRY:
        return _MODEL_BUILDER_REGISTRY[(family_enum, chosen_name)]
    if (family_enum, "auto") in _MODEL_BUILDER_REGISTRY:
        return _MODEL_BUILDER_REGISTRY[(family_enum, "auto")]
    raise ValueError(f"No model builder registered for family={family_enum.value!r} builder={chosen_name!r}")


def solve_experiment(data_stru: Dict[str, Any], options: ModelOptions, *, model_builder=None):
    """Run the common solver flow with a family-specific model builder."""
    model_builder = model_builder or resolve_model_builder(options.model_family, options.model_builder_name or options.backend_name)
    if options.sigma_fixed and options.B_form == 1:
        return solve_model_B_fix(
            data_stru,
            options.mode,
            theta=options.theta,
            sim_opt=options.sim_opt,
            B_form=options.B_form,
            sigma_fixed=options.sigma_fixed,
            LOUD=options.LOUD,
            model_builder=model_builder,
        )
    return solve_model(
        data_stru,
        options.mode,
        theta=options.theta,
        sim_opt=options.sim_opt,
        B_form=options.B_form,
        LOUD=options.LOUD,
        model_builder=model_builder,
    )


def _build_model_data1(data_stru: Dict[str, Any], mode, theta=None, sim_opt=False, B_form='single'):
    return model_construct_inter(data_stru, mode, theta=theta, sim_opt=sim_opt, B_form=B_form)


def _build_model_data2(data_stru: Dict[str, Any], mode, theta=None, sim_opt=False, B_form='single'):
    return model_construct_inter(data_stru, mode, theta=theta, sim_opt=sim_opt, B_form=B_form)


register_model_builder(WorkflowFamily.DATA1, "auto", _build_model_data1)
register_model_builder(WorkflowFamily.DATA2, "auto", _build_model_data2)


def detect_source_type(file_path: Union[str, Path]) -> Optional[SourceType]:
    """Detect the source type from the file extension."""
    suffix = Path(file_path).suffix.lower()
    if suffix == ".mat":
        return SourceType.MAT
    if suffix in {".xlsx", ".xls"}:
        return SourceType.XLSX
    return None


def _load_xlsx_minimal(xlsx_path: Path, selector: Optional[Union[str, int]] = None) -> Dict[str, Any]:
    """Minimal XLSX loader for refactored multi-file support.

    The full unified loader has a richer sheet parser. This lightweight adapter
    keeps the refactored library able to accept XLSX inputs while the refactor
    work continues.
    """
    xls = pd.ExcelFile(xlsx_path)
    sheet_name = xls.sheet_names[0] if selector is None else selector
    sheet = xls.parse(sheet_name)
    return {
        "source": SourceType.XLSX.value,
        "filename": xlsx_path.name,
        "path": str(xlsx_path),
        "sheet_name": sheet_name,
        "dataframe": sheet,
        "sheet_names": list(xls.sheet_names),
    }


def load_experiment_file(
    file_path: Union[str, Path],
    selector: Optional[Union[str, int]] = None,
) -> Tuple[Dict[str, Any], Tuple[bool, List[Tuple[str, str]]], Optional[SourceType]]:
    """Load a single experimental file into a normalized dictionary.

    MAT files use the existing `loadmat` helper.
    XLSX files use a minimal adapter that keeps the workbook and selected sheet
    available for downstream processing.
    """
    path = Path(file_path).expanduser().resolve()
    source_type = detect_source_type(path)
    if source_type is None:
        raise ValueError(f"Unsupported file type: {path.suffix}")

    issues: List[Tuple[str, str]] = []
    ok_fatal = True

    if source_type == SourceType.MAT:
        data = loadmat(str(path))
        if selector is not None and selector != "data_stru":
            data = data.get(selector, data)
        if not isinstance(data, dict):
            ok_fatal = False
            issues.append(("FATAL", f"MAT file {path.name} did not load into a dictionary."))
            data = {"raw": data}
        data.setdefault("source", SourceType.MAT.value)
        data.setdefault("path", str(path))
        data.setdefault("filename", path.name)
    else:
        data = _load_xlsx_minimal(path, selector=selector)
        issues.append(("INFO", f"Loaded XLSX workbook {path.name} with sheet selector {selector!r}."))

    return data, (ok_fatal, issues), source_type


def load_experiment_bundle(
    file_paths: Sequence[Union[str, Path]],
    selectors: Optional[Sequence[Optional[Union[str, int]]]] = None,
    *,
    profile: Optional[str] = None,
) -> ExperimentalBundle:
    """Load multiple experimental files and preserve provenance.

    This is the new multi-file entrypoint for the refactored library.
    """
    bundle = ExperimentalBundle(profile=profile)
    selectors = list(selectors) if selectors is not None else [None] * len(file_paths)
    if len(selectors) != len(file_paths):
        raise ValueError("selectors must match file_paths length")

    for file_path, selector in zip(file_paths, selectors):
        path = Path(file_path).expanduser().resolve()
        source = ExperimentalSource(path=path, selector=selector, label=path.stem)
        data, validation, source_type = load_experiment_file(path, selector=selector)
        bundle.sources.append(source)
        bundle.experiments.append(data)
        bundle.validation.append(validation)
        bundle.source_types.append(source_type)

    bundle.metadata["source_count"] = len(bundle.sources)
    bundle.metadata["profile"] = profile
    bundle.metadata["has_xlsx"] = any(st == SourceType.XLSX for st in bundle.source_types)
    bundle.metadata["has_mat"] = any(st == SourceType.MAT for st in bundle.source_types)
    return bundle


def load_experiment_files(
    file_paths: Sequence[Union[str, Path]],
    selectors: Optional[Sequence[Optional[Union[str, int]]]] = None,
    *,
    profile: Optional[str] = None,
) -> Tuple[List[Dict[str, Any]], List[Tuple[bool, List[Tuple[str, str]]]]]:
    """Convenience wrapper returning just experiments and validation results."""
    bundle = load_experiment_bundle(file_paths, selectors, profile=profile)
    return bundle.experiments, bundle.validation


def _ensure_output_dir(output_dir: Optional[Union[str, Path]] = None) -> Path:
    """Resolve an output directory for generated DATA1 figures."""
    path = Path(output_dir).expanduser().resolve() if output_dir is not None else Path("figures").resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path


def _save_current_figure(path: Union[str, Path], *, close: bool = True) -> str:
    """Save the active matplotlib figure and optionally close it."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    if close:
        plt.close()
    return str(out_path)


def plot_sim_comparison_paper(
    data_stru,
    fit_stru,
    plot_pred: bool = True,
    cond: bool = True,
    lg: bool = False,
    preface: bool = False,
    output_dir: Optional[Union[str, Path]] = None,
):
    """Plot DATA1 paper-style simulation comparisons.

    This is the Python translation of the `DiafiltrationPaperPlots.ipynb`
    helper used throughout the DATA1 figure-generation workflow.
    """
    outdir = _ensure_output_dir(output_dir)
    t_delay = data_stru['data_raw'][0]['time'][0]

    fig = plt.figure(figsize=(4, 4))
    for i in range(data_stru['data_config']['n']):
        plt.plot((data_stru['data_raw'][i]['time'] - t_delay) / 60,
                 data_stru['data_raw'][i]['mass'], 'r.', markersize=4)
        if plot_pred:
            plt.plot((fit_stru['sim_stru'][i]['time'] - t_delay) / 60,
                     fit_stru['sim_stru'][i]['mV'], 'b', linewidth=3, alpha=.6)
    plt.plot([], [], 'r.', markersize=4, label='Measurements')
    if plot_pred:
        plt.plot([], [], 'b', linewidth=3, alpha=.6, label='Predictions')
    if preface:
        plt.xlabel('Time', fontsize=24, fontweight='bold')
        plt.ylabel('Mass', fontsize=24, fontweight='bold')
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    else:
        plt.xlabel('Time [min]', fontsize=16, fontweight='bold')
        plt.ylabel('Mass [g]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    if xticks:
        xticks[0].label1.set_visible(False)
    if lg:
        plt.legend(fontsize=10, loc='best')
    fig.savefig(outdir / (('mass_preface-dat' if preface else 'mass-dat') + str(data_stru['dataset']) + '.png'),
                dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(figsize=(4, 4))
    plt.plot((data_stru['data_raw'][0]['time'][0] - t_delay) / 60,
             fit_stru['sim_stru'][0]['cF'][0],
             'ms', markersize=8, clip_on=False)
    plt.plot([], [], 'ms', markersize=8, clip_on=False, label='Retentate (Measurement)')
    if plot_pred:
        plt.plot([], [], 'g', linewidth=3, label='Retentate (Prediction)')
    plt.plot([], [], 'cs', markersize=8, label='Vial (Measurement)')
    if plot_pred:
        plt.plot([], [], 'r^', markersize=8, label='Vial (Prediction)')
        plt.plot([], [], 'r-', linewidth=3, alpha=.6, label='Permeate (Prediction)')
    for i in range(data_stru['data_config']['n']):
        plt.plot((data_stru['data_raw'][i]['time'][-1] - t_delay) / 60,
                 data_stru['data_raw'][i]['cV_avg'],
                 'cs', markersize=8)
        if cond:
            cf_exp = data_stru['data_raw'][i]['cF_exp']
            if isinstance(cf_exp, float):
                plt.plot((data_stru['data_raw'][i]['time'][-1] - t_delay) / 60, cf_exp, 'ms', markersize=8)
            elif len(cf_exp) > 1:
                plt.plot((data_stru['data_raw'][i]['time'] - t_delay) / 60, cf_exp, 'ms', markersize=8)
        if plot_pred:
            plt.plot((fit_stru['sim_stru'][i]['time'] - t_delay) / 60,
                     fit_stru['sim_stru'][i]['cF'], 'g', linewidth=3, alpha=.6)
            plt.plot((fit_stru['sim_stru'][i]['time'] - t_delay) / 60,
                     fit_stru['sim_stru'][i]['cH'], 'r-', linewidth=3, alpha=.6)
            if not preface:
                plt.plot((fit_stru['sim_stru'][i]['time'][-1] - t_delay) / 60,
                         fit_stru['sim_stru'][i]['cV'][-1],
                         'r^', markersize=8, alpha=.6)
    if not cond:
        plt.plot((data_stru['data_raw'][-1]['time'][-1] - t_delay) / 60,
                 data_stru['data_raw'][-1]['cF_exp'], 'ms', markersize=8)
    if preface:
        plt.xlabel('Time', fontsize=24, fontweight='bold')
        plt.ylabel('Concentration', fontsize=24, fontweight='bold')
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
        plt.annotate('Retentate', xy=((fit_stru['sim_stru'][5]['time'][-1] - t_delay) / 60,
                     fit_stru['sim_stru'][5]['cF'][-1]),  xycoords='data',
                     xytext=(36, 50), weight='bold', textcoords='offset points',
                     size=20, ha='right', va="center",
                     bbox=dict(boxstyle="round", color="green", alpha=0.1),
                     arrowprops=dict(arrowstyle="wedge,tail_width=0.3", color="green", alpha=0.1))
        plt.annotate('Permeate', xy=((data_stru['data_raw'][7]['time'][-1] - t_delay) / 60,
                 data_stru['data_raw'][7]['cV_avg']),  xycoords='data',
                     xytext=(36, 50), weight='bold', textcoords='offset points',
                     size=20, ha='right', va="center",
                     bbox=dict(boxstyle="round", color="red", alpha=0.1),
                     arrowprops=dict(arrowstyle="wedge,tail_width=0.3", color="red", alpha=0.1))
    else:
        plt.xlabel('Time [min]', fontsize=16, fontweight='bold')
        plt.ylabel('Concentration [mM]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in", top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    if xticks:
        xticks[0].label1.set_visible(False)
    if lg:
        plt.legend(fontsize=10, loc='best')
    fig.savefig(outdir / (('concentration_preface-dat' if preface else 'concentration-dat') + str(data_stru['dataset']) + '.png'),
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    return {'output_dir': str(outdir)}


def plot_contour(df, show_title: bool = False, preface: bool = False, output_dir: Optional[Union[str, Path]] = None):
    """Plot objective contours used in DATA1 paper figures."""
    outdir = _ensure_output_dir(output_dir)
    f_m = df.Obj_mass.values
    ind_m = np.argmin(f_m)
    f_pc = df.Obj_concentration.values
    ind_pc = np.argmin(f_pc)
    f_rc = df.Obj_retentate_concentration.values
    ind_rc = np.argmin(f_rc)
    F_m = np.reshape(f_m, (50, 50))
    F_pc = np.reshape(f_pc, (50, 50))
    F_rc = np.reshape(f_rc, (50, 50))
    if 'B' in df:
        xx = df.B.values
        if preface:
            xlabelstr = 'B'
            axfontsize = 24
            fixstr = 'preface_fixsig'
        else:
            xlabelstr = 'B [$\\mathbf{\\mu}$m $\\mathbf{\\cdot}$ s$\\mathbf{^{-1}}$]'
            axfontsize = 16
            fixstr = 'fixsig'
    else:
        xx = df.sigma.values
        if preface:
            xlabelstr = '$\\mathbf{\\sigma}$'
            axfontsize = 24
            fixstr = 'preface_fixB'
        else:
            xlabelstr = '$\\mathbf{\\sigma}$ [dimensionless]'
            axfontsize = 16
            fixstr = 'fixB'
    ylabelstr = 'L$\\mathbf{_p}$' if preface else r'L$\mathbf{_p}$ [L $\mathbf{ \cdot}$ m$\mathbf{^{-2} \cdot}$h$\mathbf{^{-1} \cdot}$bar$\mathbf{^{-1}}$]'
    yy = df.Lp.values
    X = np.reshape(xx, (50, 50))
    Y = np.reshape(yy, (50, 50))

    fig = plt.figure(1, figsize=(4, 4))
    cp = plt.contour(X, Y, F_m, 10, linewidths=2)
    plt.clabel(cp, cp.levels[::2], inline=True, fontsize=12, colors='k', fmt='%1.1f')
    plt.plot(xx[ind_m], yy[ind_m], '^', markersize=12, markeredgecolor='red', markerfacecolor=[1, .6, .6], clip_on=False)
    if show_title:
        plt.title('Log$\\mathbf{_e}$ transformed \n Mass Objective', fontsize=16, fontweight='bold')
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight='bold')
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight='bold')
    if preface:
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    fig.savefig(outdir / f'contour_{fixstr}-mass.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(2, figsize=(4, 4))
    cp = plt.contour(X, Y, F_pc, 10, linewidths=2)
    plt.clabel(cp, cp.levels[::2], inline=True, fontsize=12, colors='k', fmt='%1.1f')
    plt.plot(xx[ind_pc], yy[ind_pc], '^', markersize=12, markeredgecolor='red', markerfacecolor=[1, .6, .6], clip_on=False)
    if show_title:
        plt.title('Log$\\mathbf{_e}$ transformed \n Permeate Concentration Objective', fontsize=16, fontweight='bold')
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight='bold')
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight='bold')
    if preface:
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    fig.savefig(outdir / f'contour_{fixstr}-permeate_conc.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(3, figsize=(4, 4))
    cp = plt.contour(X, Y, F_rc, 10, linewidths=2)
    plt.clabel(cp, cp.levels[::2], inline=True, fontsize=12, colors='k', fmt='%1.1f')
    plt.plot(xx[ind_rc], yy[ind_rc], '^', markersize=12, markeredgecolor='red', markerfacecolor=[1, .6, .6], clip_on=False)
    if show_title:
        plt.title('Log$\\mathbf{_e}$ transformed \n Retentate Concentration Objective', fontsize=16, fontweight='bold')
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight='bold')
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight='bold')
    if preface:
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    fig.savefig(outdir / f'contour_{fixstr}-retentate_conc.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return {'output_dir': str(outdir)}


def plot_sim(sim_stru, color, linetype):
    """Plot sigma-sensitivity simulation series."""
    t_delay = sim_stru[0]['time'][0]
    plt.figure(1, figsize=(4, 4))
    for i in sim_stru:
        plt.plot((i['time'] - t_delay) / 60, i['mV'], color, linestyle=linetype, linewidth=2, alpha=.8)
    plt.figure(2, figsize=(4, 4))
    for i in sim_stru:
        plt.plot((i['time'] - t_delay) / 60, i['cF'], color, linestyle=linetype, linewidth=2, alpha=.8)
    plt.figure(3, figsize=(4, 4))
    for i in sim_stru:
        plt.plot((i['time'] - t_delay) / 60, i['cH'], color, linestyle=linetype, linewidth=2, alpha=.8)


def plot_sim_show(sig, colors, output_dir: Optional[Union[str, Path]] = None):
    """Finalize sigma-sensitivity plot figures and legends."""
    outdir = _ensure_output_dir(output_dir)
    custom_lines = [Line2D([0], [0], color=colors[0], ls='--', lw=3),
                    Line2D([0], [0], color=colors[1], ls='-', lw=3),
                    Line2D([0], [0], color=colors[2], ls=':', lw=3)]
    legendname = ['$\\sigma$ = ' + str(sig[0]), '$\\sigma$ = ' + str(sig[1]), '$\\sigma$ = ' + str(sig[2])]
    fig = plt.figure(1)
    plt.xlabel('Time [min]', fontsize=16, fontweight='bold')
    plt.ylabel('Mass [g]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.ylim(bottom=0)
    fig.savefig(outdir / 'sigma_sensitivity-mass.png', dpi=300, bbox_inches='tight')
    fig = plt.figure(2)
    plt.xlabel('Time [min]', fontsize=16, fontweight='bold')
    plt.ylabel('Retentate Conc. [mM]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.legend(custom_lines, legendname, fontsize=15, loc='best')
    fig.savefig(outdir / 'sigma_sensitivity-reten_conc.png', dpi=300, bbox_inches='tight')
    fig = plt.figure(3)
    plt.xlabel('Time [min]', fontsize=16, fontweight='bold')
    plt.ylabel('Permeate Conc. [mM]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    fig.savefig(outdir / 'sigma_sensitivity-perme_conc.png', dpi=300, bbox_inches='tight')
    plt.close('all')
    return {'output_dir': str(outdir)}


def plot_contour_sig_sen(data_stru, contour_sig_stru, Vmin, Vmax, Level, Colorbar_ticks, Manual_locations, name_append, filled=False, Jw_filter=True, bar=False, output_dir: Optional[Union[str, Path]] = None):
    """Plot sigma-sensitivity contour figures."""
    outdir = _ensure_output_dir(output_dir)
    diaf_mode = False
    if 'cf0' in contour_sig_stru:
        x = contour_sig_stru['cf0']
        xlabelstr = '$\\mathbf{c_f}$(t=0) [mM]'
    else:
        x = contour_sig_stru['cd']
        xlabelstr = '$\\mathbf{c_d}$ [mM]'
        diaf_mode = True
    y = contour_sig_stru['delp']
    ylabelstr = '$\\mathbf{\\Delta}$P [psi]'
    X, Y = np.meshgrid(x, y)
    axfontsize = 16
    R = 8.314e-5
    T = data_stru['data_config']['Temp']
    ni = data_stru['data_config']['ni']
    if not diaf_mode:
        Jw0 = Y / 14.504 - 1 * ni * R * T * X
        jw0 = 1 * ni * R * T * np.array(x) * 14.504
    else:
        cf0 = 5
        Jw0 = Y / 14.504 - 1 * ni * R * T * cf0
        jw0 = 1 * ni * R * T * cf0 * 14.504

    def patch_back():
        ax = plt.gca()
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        p = patches.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, color='gray', zorder=-10)
        ax.add_patch(p)

    fig = plt.figure(1, figsize=(4, 4))
    if filled:
        plt.contourf(X, Y, contour_sig_stru['range_mV'], 50, cmap='spring', vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X, Y, contour_sig_stru['range_mV'], Level[0], linewidths=3, colors='k')
    else:
        cp = plt.contour(X, Y, contour_sig_stru['range_mV'], Level[0], linewidths=2)
    plt.clabel(cp, inline=True, manual=Manual_locations[0], fontsize=15, colors='k', fmt='%1.1f', zorder=2)
    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color='gray', zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors='orangered')
        plt.setp(cg.collections, path_effects=[patheffects.withTickedStroke(angle=300, length=2)])
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight='bold')
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom=y[0])
    fig.savefig(outdir / f'contour_sensitivity{name_append}-mass.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(2, figsize=(4, 4))
    if filled:
        plt.contourf(X, Y, contour_sig_stru['range_cF'], 100, cmap='spring', vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X, Y, contour_sig_stru['range_cF'], Level[1], linewidths=3, colors='k')
    else:
        cp = plt.contour(X, Y, contour_sig_stru['range_cF'], Level[1], linewidths=2)
    plt.clabel(cp, inline=True, manual=Manual_locations[1], fontsize=15, colors='k', fmt='%1.1f', zorder=2)
    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color='gray', zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors='orangered')
        plt.setp(cg.collections, path_effects=[patheffects.withTickedStroke(angle=300, length=2)])
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight='bold')
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom=y[0])
    fig.savefig(outdir / f'contour_sensitivity{name_append}-retentate_conc.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(3, figsize=(4, 4))
    if filled:
        plt.contourf(X, Y, contour_sig_stru['range_cH'], 50, cmap='spring', vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X, Y, contour_sig_stru['range_cH'], Level[2], linewidths=3, colors='k')
    else:
        cp = plt.contour(X, Y, contour_sig_stru['range_cH'], Level[2], linewidths=2)
    plt.clabel(cp, inline=True, manual=Manual_locations[2], fontsize=15, colors='k', fmt='%1.1f', zorder=2)
    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color='gray', zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors='orangered')
        plt.setp(cg.collections, path_effects=[patheffects.withTickedStroke(angle=300, length=2)])
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight='bold')
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom=y[0])
    fig.savefig(outdir / f'contour_sensitivity{name_append}-permeate_conc.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    if filled and bar:
        cmap = mpl.cm.spring
        norm = mpl.colors.Normalize(vmin=Vmin, vmax=Vmax)
        fig, ax = plt.subplots(figsize=(12, 1))
        fig.subplots_adjust(bottom=0.5)
        mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal', ticks=Colorbar_ticks, extend='max')
        ax.tick_params(labelsize=15)
        fig.savefig(outdir / f'colorbar{name_append}-horizontal.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(1, 4))
        fig.subplots_adjust(left=0.5)
        mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical', ticks=Colorbar_ticks, extend='max')
        ax.tick_params(labelsize=15)
        fig.savefig(outdir / f'colorbar{name_append}-vertical.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
    return {'output_dir': str(outdir)}


def plot_conc_range(df_f, df_d, lg: bool = True, output_dir: Optional[Union[str, Path]] = None):
    """Plot the DATA1 experiment space concentration range."""
    outdir = _ensure_output_dir(output_dir)
    fig = plt.figure(figsize=(4, 4))
    for i in range(len(df_f.columns) // 2):
        plt.plot(df_f[f'F{i+1}_cf'], df_f[f'F{i+1}_cp'], '^', markersize=8, alpha=.8, clip_on=False)
    for i in range(len(df_d.columns) // 2):
        plt.plot(df_d[f'D{i+1}_cf'], df_d[f'D{i+1}_cp'], 's', markersize=8, alpha=.8, clip_on=False)
    plt.plot([], [], 'k^', markersize=8, markerfacecolor='white', label='Filtration')
    plt.plot([], [], 'ks', markersize=8, markerfacecolor='white', label='Diafiltration')
    plt.xlabel('Retentate [mM]', fontsize=16, fontweight='bold')
    plt.ylabel('Permeate [mM]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    if lg:
        plt.legend(fontsize=15, loc='best')
    fig.savefig(outdir / 'concentration_range.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return {'output_dir': str(outdir)}


def plot_cr_measure(data_stru, fit_stru, cr_pred, ybottom, lg: bool = False, cond: bool = True, output_dir: Optional[Union[str, Path]] = None):
    """Plot retentate concentration comparison against time."""
    outdir = _ensure_output_dir(output_dir)
    t_delay = data_stru['data_raw'][0]['time'][0]
    fig = plt.figure(figsize=(4, 4))
    if cond:
        plt.plot((data_stru['data_raw'][0]['time'][0] - t_delay) / 60, fit_stru['sim_stru'][0]['cF'][0], 'ms', markersize=8, clip_on=False)
    else:
        plt.plot((data_stru['data_raw'][0]['time'][0] - t_delay) / 60, data_stru['data_config']['C_F0'], 'go', markersize=8, clip_on=False)
    for i in range(data_stru['data_config']['n']):
        if cond:
            cf_exp = data_stru['data_raw'][i]['cF_exp']
            if isinstance(cf_exp, float):
                plt.plot((data_stru['data_raw'][i]['time'][-1] - t_delay) / 60, cf_exp, 'ms', markersize=8)
            elif len(cf_exp) > 1:
                plt.plot((data_stru['data_raw'][i]['time'] - t_delay) / 60, cf_exp, 'ms', markersize=8)
    if len(cr_pred) > 0:
        for i in range(data_stru['data_config']['n']):
            plt.plot((data_stru['data_raw'][i]['time'][-1] - t_delay) / 60, cr_pred[i], 'go', markersize=8)
    plt.xlabel('Time [min]', fontsize=16, fontweight='bold')
    plt.ylabel('Concentration [mM]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in", top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=ybottom)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    if xticks:
        xticks[0].label1.set_visible(False)
    plt.plot([], [], 'ms', markersize=8, label='Retentate \n (Measurements)')
    plt.plot([], [], 'go', markersize=8, label='Retentate \n (Calculated)')
    if lg:
        plt.legend(fontsize=15, loc='best')
    fig.savefig(outdir / f'cr_measure-dat{data_stru["dataset"]}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return {'output_dir': str(outdir)}


def calib_curve_cond(calib_curve, output_dir: Optional[Union[str, Path]] = None):
    """Plot the conductivity calibration curve."""
    outdir = _ensure_output_dir(output_dir)
    x = calib_curve.Conductivity.values
    y = calib_curve.Concentration.values
    fig = plt.figure(figsize=(4, 4))
    plt.plot(x, y, 'bo', markersize=8)
    z = np.polyfit(x, y, 1)
    l = np.poly1d(z)
    r_squared = np.corrcoef(x, y)[0, 1] ** 2
    plt.plot(x, l(x), 'b:', linewidth=3, alpha=.7)
    eqn = 'y=%.3fx+%.2f \n R$\\mathbf{^{2}}$=%.4f' % (z[0], z[1], r_squared)
    plt.annotate(eqn, xy=(x[4], y[4]), xycoords='data', xytext=(-30, 60), weight='bold', textcoords='offset points',
                 size=12, ha='center', va="center", bbox=dict(boxstyle="round", color="b", alpha=0.1))
    plt.xlabel('Conductivity [$\\mathbf{\\mu}$S $\\mathbf{\\cdot}$ cm$\\mathbf{^{-1}}$]', fontsize=16, fontweight='bold')
    plt.ylabel('Concentration [mM]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in", top=True, right=True)
    fig.savefig(outdir / 'calib_curve.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return {'output_dir': str(outdir), 'r_squared': float(r_squared), 'coefficients': [float(z[0]), float(z[1])]}


def plot_conc_comparison(data_stru_f, fit_stru_f, data_stru_d, fit_stru_d, plot_pred: bool = True, output_dir: Optional[Union[str, Path]] = None):
    """Compare concentration measurements between filtration and diafiltration datasets."""
    outdir = _ensure_output_dir(output_dir)
    t_delay_f = data_stru_f['data_raw'][0]['time'][0]
    t_delay_d = data_stru_d['data_raw'][0]['time'][0]
    fig = plt.figure(figsize=(4, 4))
    plt.plot(data_stru_f['data_raw'][0]['time'][0] - t_delay_f, data_stru_f['cF_ICP'][0], 'mv', markersize=6, clip_on=False)
    plt.plot(data_stru_d['data_raw'][0]['time'][0] - t_delay_d, data_stru_d['cF_ICP'][0], 'k^', markersize=6, clip_on=False)
    plt.plot(data_stru_f['data_raw'][0]['time'][0] - t_delay_f, fit_stru_f['sim_stru'][0]['cF'][0], 'bs', markersize=6, clip_on=False)
    plt.plot(data_stru_d['data_raw'][0]['time'][0] - t_delay_d, fit_stru_d['sim_stru'][0]['cF'][0], 'rs', markersize=6, clip_on=False)
    plt.plot([], [], 'kv', markersize=6, clip_on=False, label='Filtration Retentate ICP')
    plt.plot([], [], 'bs', markersize=6, clip_on=False, label='Filtration Retentate conductivity')
    plt.plot([], [], 'bo', markersize=6, label='Filtration Vial ICP')
    plt.plot([], [], 'k^', markersize=6, clip_on=False, label='Diafiltration Retentate ICP')
    plt.plot([], [], 'rs', markersize=6, clip_on=False, label='Diafiltration Retentate conductivity')
    plt.plot([], [], 'ro', markersize=6, label='Diafiltration Vial ICP')
    if plot_pred:
        plt.plot([], [], 'g-.', linewidth=2, label='Filtration Retentate')
        plt.plot([], [], 'g', linewidth=2, label='Diafiltration Retentate')
        plt.plot([], [], 'k-.', linewidth=2, alpha=.6, label='Filtration Permeate')
        plt.plot([], [], 'k', linewidth=2, alpha=.6, label='Diafiltration Permeate')
    for i in range(data_stru_f['data_config']['n']):
        plt.plot(data_stru_f['data_raw'][i]['time'][-1] - t_delay_f, data_stru_f['data_raw'][i]['cV_avg'], 'bo', markersize=6)
        plt.plot(data_stru_f['data_raw'][i]['time'][-1] - t_delay_f, data_stru_f['data_raw'][i]['cF_exp'], 'bs', markersize=6)
        if plot_pred:
            plt.plot(fit_stru_f['sim_stru'][i]['time'] - t_delay_f, fit_stru_f['sim_stru'][i]['cF'], 'g-.', linewidth=2)
            plt.plot(fit_stru_f['sim_stru'][i]['time'] - t_delay_f, fit_stru_f['sim_stru'][i]['cH'], 'k-.', linewidth=2, alpha=.6)
    for i in range(data_stru_d['data_config']['n']):
        plt.plot(data_stru_d['data_raw'][i]['time'][-1] - t_delay_d, data_stru_d['data_raw'][i]['cV_avg'], 'ro', markersize=6)
        plt.plot(data_stru_d['data_raw'][i]['time'][-1] - t_delay_d, data_stru_d['data_raw'][i]['cF_exp'], 'rs', markersize=6)
        if plot_pred:
            plt.plot(fit_stru_d['sim_stru'][i]['time'] - t_delay_d, fit_stru_d['sim_stru'][i]['cF'], 'g', linewidth=2)
            plt.plot(fit_stru_d['sim_stru'][i]['time'] - t_delay_d, fit_stru_d['sim_stru'][i]['cH'], 'k', linewidth=2, alpha=.6)
    plt.plot(data_stru_f['data_raw'][-1]['time'][-1] - t_delay_f, data_stru_f['cF_ICP'][-1], 'kv', markersize=6, clip_on=False)
    plt.plot(data_stru_d['data_raw'][-1]['time'][-1] - t_delay_d, data_stru_d['cF_ICP'][-1], 'k^', markersize=6, clip_on=False)
    plt.xlabel('Time [s]', fontsize=16)
    plt.ylabel('Concentration [mM]', fontsize=16)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.minorticks_on()
    plt.tick_params(direction="in", top=True, right=True)
    plt.tick_params(which="minor", direction="in", top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(fontsize=12, bbox_to_anchor=(1.1, 1.05), loc='upper left')
    fig.savefig(outdir / 'concentration.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return {'output_dir': str(outdir)}


def plot_sim_comparison(*args, **kwargs):
    """Compatibility wrapper for both the legacy and paper-style helpers.

    Legacy DATA2 scripts call ``plot_sim_comparison(data_stru, sim_stru, ...)``.
    The newer DATA1 paper workflow calls ``plot_sim_comparison(data_stru, fit_stru, ...)``.
    This wrapper dispatches based on the shape of the second positional argument.
    """
    if len(args) >= 2:
        second = args[1]
        if isinstance(second, list) or 'stirc_mass' in kwargs or 'LOUD' in kwargs:
            return plot_sim_comparison_legacy(*args, **kwargs)
        if isinstance(second, dict) and 'sim_stru' in second:
            return plot_sim_comparison_paper(*args, **kwargs)
    if 'preface' in kwargs or 'cond' in kwargs:
        return plot_sim_comparison_paper(*args, **kwargs)
    return plot_sim_comparison_legacy(*args, **kwargs)


def _data2_load_frame(data: Union[str, Path, pd.DataFrame]) -> pd.DataFrame:
    """Load a DATA2 helper frame from a CSV path or return a copy of the frame."""
    if isinstance(data, pd.DataFrame):
        return data.copy()
    return pd.read_csv(Path(data).expanduser().resolve())


def plot_data2_pressure_change(
    data_pd: Union[str, Path, pd.DataFrame],
    *,
    time_offset: float,
    transition_time: Optional[float] = None,
    output_dir: Optional[Union[str, Path]] = None,
    filename: str = 'pressure_change.png',
    pressure_ylim: Optional[Tuple[float, float]] = None,
    concentration_ylim: Optional[Tuple[float, float]] = None,
    xlim: Optional[Tuple[float, float]] = None,
) -> Dict[str, Any]:
    """Plot the DATA2 pressure and concentration profile with twin y-axes."""
    outdir = _ensure_output_dir(output_dir)
    df = _data2_load_frame(data_pd)
    color_pressure = 'tab:orange'
    color_conc = 'm'
    fig, ax1 = plt.subplots(figsize=(4, 4))
    ax2 = ax1.twinx()
    ax1.plot((df['Time'] - time_offset) / 60, df['Pressure'], 'o', color=color_pressure, markersize=6)
    ax2.plot((df['Time'] - time_offset) / 60, df['Concentration'], 's', color=color_conc, markersize=6)
    ax1.plot([0, 0], [0, 100], '--', color='tab:blue', lw=2)
    if transition_time is not None:
        ax1.plot([(transition_time - time_offset) / 60, (transition_time - time_offset) / 60], [0, 100], '--', color='tab:green', lw=2)
    ax1.set_xlabel('Time [min]', fontsize=16, fontweight='bold')
    ax1.set_ylabel('Applied Pressure [psi]', color=color_pressure, fontsize=16, fontweight='bold')
    ax1.tick_params(axis='x', which='major', direction='in', labelsize=12)
    ax1.tick_params(axis='x', which='minor', direction='in', labelsize=8)
    ax1.tick_params(axis='y', labelcolor=color_pressure, direction='in', labelsize=12)
    if xlim is not None:
        ax1.set_xlim(*xlim)
    ax1.set_ylim(*(pressure_ylim if pressure_ylim is not None else (0, 68)))
    ax2.set_ylabel('Retentate [mM]', color=color_conc, fontsize=16, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor=color_conc, direction='in', labelsize=12)
    ax2.set_ylim(*(concentration_ylim if concentration_ylim is not None else (15, 30)))
    out_path = outdir / filename
    fig.savefig(out_path, dpi=600, bbox_inches='tight', transparent=True)
    plt.close(fig)
    return {'output_dir': str(outdir), 'figure': str(out_path), 'rows': int(len(df))}


def data2_model_predictions(c_in, c_h, k0, k1, Pe):
    """Calculate Js/Jw for the DATA2 NF270 regression model."""
    Kf = k1 * c_in + k0
    Kp = k1 * c_h + k0
    return (c_in * Kf * np.exp(Pe) - Kp * c_h) / (np.exp(Pe) - 1)


def data2_model_convection(sim_data, Pe_fixed_value=None, log_transform_Pe=False):
    """Fit the DATA2 convection-diffusion regression model from the notebook."""
    model = ConcreteModel()
    model.i = RangeSet(sim_data.index[0] + 1, sim_data.index[-1] + 1)
    model.J_w = Param(model.i, initialize=lambda m, i: sim_data['Jw'][i - 1], mutable=True)
    model.J_s = Param(model.i, initialize=lambda m, i: sim_data['Js'][i - 1], mutable=True)
    model.c_in = Param(model.i, initialize=lambda m, i: sim_data['cIn'][i - 1], mutable=True)
    model.c_h = Param(model.i, initialize=lambda m, i: sim_data['cH'][i - 1], mutable=True)
    model.Js = Var(model.i, initialize=lambda m, i: sim_data['Jw'][i - 1])
    model.Kp = Var(model.i, initialize=1.0, bounds=(0.01, 1.2))
    model.Kf = Var(model.i, initialize=1.0, bounds=(0.01, 1.2))
    model.k1 = Var(initialize=1e-2)
    model.k0 = Var(initialize=1.0, bounds=(1e-4, 2))

    pe_lower = min(0.1, Pe_fixed_value) if Pe_fixed_value is not None else 0.1
    pe_upper = max(20, Pe_fixed_value) if Pe_fixed_value is not None else 20
    if log_transform_Pe:
        model.expPe = Var(initialize=np.exp(15), within=Reals, bounds=(np.exp(pe_lower), np.exp(pe_upper)))
    else:
        model.Pe = Var(initialize=15, within=Reals, bounds=(pe_lower, pe_upper))
        model.expPe = Expression(expr=exp(model.Pe))
    if Pe_fixed_value is not None:
        if log_transform_Pe:
            model.expPe.fix(np.exp(Pe_fixed_value))
        else:
            model.Pe.fix(Pe_fixed_value)

    @model.Constraint(model.i)
    def partition_f(m, i):
        return m.Kf[i] == m.k1 * m.c_in[i] + m.k0

    @model.Constraint(model.i)
    def partition_p(m, i):
        return m.Kp[i] == m.k1 * m.c_h[i] + m.k0

    @model.Constraint(model.i)
    def convection(m, i):
        return m.Js[i] * (m.expPe - 1) == m.J_w[i] * (m.c_in[i] * m.Kf[i] * m.expPe - m.Kp[i] * m.c_h[i])

    model.FirstStageCost = Expression(expr=0)
    model.SecondStageCost = Expression(expr=sum((model.Js[i] / model.J_w[i] - model.J_s[i] / model.J_w[i]) ** 2 for i in model.i))
    model.Total_Cost_Objective = Objective(expr=model.FirstStageCost + model.SecondStageCost, sense=minimize)

    solver = SolverFactory('ipopt')
    solver.options['linear_solver'] = 'ma97'
    solver.options['halt_on_ampl_error'] = 'yes'
    solver.options['acceptable_tol'] = 1e-8
    solver.solve(model, tee=True)

    pe_value = np.log(value(model.expPe)) if log_transform_Pe else value(model.Pe)
    theta_fit = {'k0': value(model.k0), 'k1': value(model.k1), 'Pe': float(pe_value)}
    return model, theta_fit


def plot_data2_model_predictions(
    sim_data: pd.DataFrame,
    k0: float,
    k1: float,
    Pe: float,
    *,
    output_dir: Optional[Union[str, Path]] = None,
    filename: Optional[str] = None,
) -> Dict[str, Any]:
    """Render the DATA2 3D wireframe prediction plot."""
    outdir = _ensure_output_dir(output_dir)
    c_in = sim_data['cIn'].values
    c_h = sim_data['cH'].values
    Js = sim_data['Js'].values
    Jw = sim_data['Jw'].values
    round_to = 5
    c_in_low = np.floor(np.min(c_in) / round_to) * round_to
    c_in_high = np.ceil(np.max(c_in) / round_to) * round_to
    c_h_low = np.floor(np.min(c_h) / round_to) * round_to
    c_h_high = np.ceil(np.max(c_h) / round_to) * round_to
    c_in_grid, c_h_grid = np.meshgrid(np.linspace(c_in_low, c_in_high, 100), np.linspace(c_h_low, c_h_high, 100))
    Js_Jw_model_grid = data2_model_predictions(c_in_grid, c_h_grid, k0, k1, Pe)
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d', computed_zorder=False)
    ax.mouse_init()
    ax.plot_wireframe(c_in_grid, c_h_grid, Js_Jw_model_grid, color='blue', alpha=0.7, label='Regressed Model')
    ax.scatter(c_in, c_h, Js / Jw, color='red', label='Experimental Data', marker='o', s=50)
    ax.legend(fontsize=10)
    ax.set_xlabel('c$_{in}$ [mM]', fontsize=12, fontweight='bold')
    ax.set_ylabel('c$_{h}$ [mM]', fontsize=12, fontweight='bold')
    ax.set_zlabel('J$_s$/J$_w$', fontsize=12, fontweight='bold')
    ax.set_title(f'J$_s$/J$_w$ with Pe={Pe:.1f}', fontsize=14, fontweight='bold')
    ax.set_box_aspect([1, 1, 0.8])
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.8)
    out_name = filename or f'data2_model_predictions_pe_{Pe:.2f}.png'
    out_path = outdir / out_name
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return {'output_dir': str(outdir), 'figure': str(out_path)}


def data2_linear_regression(data: pd.DataFrame, Pe: float):
    """Fit k0 and k1 for a fixed Peclet number using ordinary least squares."""
    import statsmodels.api as sm
    expPe = np.exp(Pe)
    y = data['Js'].values / data['Jw'].values
    x0 = (data['cIn'].values * expPe - data['cH'].values) / (expPe - 1)
    x1 = (data['cIn'].values ** 2 * expPe - data['cH'].values ** 2) / (expPe - 1)
    X = np.column_stack((x0, x1))
    results = sm.OLS(y, X).fit()
    k0, k1 = results.params[0], results.params[1]
    return k0, k1, results.mse_resid, results.bse[0], results.bse[1], results


def plot_data2_regression_sensitivity(
    data: pd.DataFrame,
    pe_values: Optional[Sequence[float]] = None,
    *,
    output_dir: Optional[Union[str, Path]] = None,
    prefix: str = 'data2_regression_sensitivity',
) -> Dict[str, Any]:
    """Recreate the notebook's Peclet-number sensitivity plots."""
    outdir = _ensure_output_dir(output_dir)
    pe_values = np.array(pe_values if pe_values is not None else np.logspace(-3, 2, 51))
    objective = np.zeros(len(pe_values))
    k0_values = np.zeros(len(pe_values))
    k1_values = np.zeros(len(pe_values))
    k0_se = np.zeros(len(pe_values))
    k1_se = np.zeros(len(pe_values))
    for i, pe in enumerate(pe_values):
        k0_values[i], k1_values[i], objective[i], k0_se[i], k1_se[i], _ = data2_linear_regression(data, pe)

    fig = plt.figure(figsize=(4, 4))
    plt.plot(pe_values, objective, 'b-', linewidth=3)
    plt.xscale('log')
    plt.xlabel('Peclet Number', fontsize=16, fontweight='bold')
    plt.ylabel('Mean Squared Error [mM$\\mathbf{^{2}}$]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction='in')
    plt.title('Regression Objective', fontsize=16, fontweight='bold', loc='left')
    plt.grid(True)
    fig.savefig(outdir / f'{prefix}_objective.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(figsize=(4, 4))
    plt.plot(pe_values, k0_values, 'r-', linewidth=3)
    plt.xscale('log')
    plt.xlabel('Peclet Number', fontsize=16, fontweight='bold')
    plt.ylabel('h$\\mathbf{_0}$ [-]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction='in')
    plt.title('k0 vs Peclet Number', fontsize=16, fontweight='bold', loc='left')
    plt.grid(True)
    fig.savefig(outdir / f'{prefix}_k0.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(figsize=(4, 4))
    plt.plot(pe_values, k1_values, 'g-', linewidth=3)
    plt.xscale('log')
    plt.xlabel('Peclet Number', fontsize=16, fontweight='bold')
    plt.ylabel('h$\\mathbf{_1}$ [mM$\\mathbf{^{-1}}$]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction='in')
    plt.title('k1 vs Peclet Number', fontsize=16, fontweight='bold', loc='left')
    plt.grid(True)
    fig.savefig(outdir / f'{prefix}_k1.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    return {
        'output_dir': str(outdir),
        'objective': objective.tolist(),
        'k0': k0_values.tolist(),
        'k1': k1_values.tolist(),
        'k0_se': k0_se.tolist(),
        'k1_se': k1_se.tolist(),
        'pe_values': pe_values.tolist(),
    }


def run_data2_model_demo(
    *,
    lag_dataset: Optional[Union[str, Path]] = None,
    overflow_dataset: Optional[Union[str, Path]] = None,
    lag_theta: Optional[Dict[str, float]] = None,
    overflow_theta: Optional[Dict[str, float]] = None,
    output_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Bundle the DATA2 model demo notebook into a callable workflow."""
    outdir = _ensure_output_dir(output_dir)
    summary: Dict[str, Any] = {'output_dir': str(outdir), 'cases': []}
    lag_theta = lag_theta or {'Lp': 11, 'beta_c': 15, 'beta_0': 1, 'beta_1': 0.01, 'sigma': 1.0, 'S0': 0}
    overflow_theta = overflow_theta or {'Lp': 11, 'beta_c': 15, 'beta_0': 1, 'beta_1': 0.01, 'sigma': 1.0, 'S0': -0.1756665334051245}
    if lag_dataset is not None:
        data_stru = loadmat(str(lag_dataset))['data_stru']
        plot_sim_comparison(data_stru, [], plot_pred=False, lg=True, LOUD=False)
        fit_stru, sim_stru, sim_inter = solve_experiment(
            data_stru,
            ModelOptions(model_family=WorkflowFamily.DATA2, mode='Lag', theta=lag_theta, sim_opt=False, B_form=1, LOUD=False, sigma_fixed=False),
        )
        plot_sim_comparison(data_stru, sim_stru, stirc_mass=True, plot_pred=True, lg=False, LOUD=False)
        summary['cases'].append({'case': 'lag', 'dataset': str(lag_dataset), 'fit': fit_stru})
    if overflow_dataset is not None:
        data_stru = loadmat(str(overflow_dataset))['data_stru']
        plot_sim_comparison(data_stru, [], plot_pred=False, lg=True, LOUD=False)
        fit_stru, sim_stru, sim_inter = solve_experiment(
            data_stru,
            ModelOptions(model_family=WorkflowFamily.DATA2, mode='Overflow', theta=overflow_theta, sim_opt=False, B_form=1, LOUD=False, sigma_fixed=False),
        )
        plot_sim_comparison(data_stru, sim_stru, stirc_mass=True, plot_pred=True, lg=False, LOUD=False)
        summary['cases'].append({'case': 'overflow', 'dataset': str(overflow_dataset), 'fit': fit_stru})
    return summary


def run_data2_cross_verification(
    *,
    dataset_paths: Optional[Sequence[Union[str, Path]]] = None,
    base_parameters: Optional[Dict[str, float]] = None,
    mode: str = 'DATA',
    output_dir: Optional[Union[str, Path]] = None,
    b_form: Union[int, str] = 1,
) -> Dict[str, Any]:
    """Run the DATA2 empirical-B cross-verification workflow."""
    outdir = _ensure_output_dir(output_dir)
    dataset_paths = list(dataset_paths) if dataset_paths is not None else [
        'data_library/data_stru-dataset270611.121.mat',
        'data_library/data_stru-dataset270711.121.mat',
        'data_library/data_stru-dataset270511.221.mat',
        'data_library/data_stru-dataset270511.321.mat',
        'data_library/data_stru-dataset270511.421.mat',
        'data_library/data_stru-dataset270511.921.mat',
        'data_library/data_stru-dataset270511.521.mat',
        'data_library/data_stru-dataset270511.621.mat',
        'data_library/data_stru-dataset270511.721.mat',
        'data_library/data_stru-dataset270511.821.mat',
    ]
    base_parameters = base_parameters or {'Lp': 11.113241068147595, 'beta_0': 1.0649294788103598, 'beta_1': 0.015171421224603474, 'sigma': 1.0}
    cases = []
    for ds in dataset_paths:
        data_stru = loadmat(str(ds))['data_stru']
        fit_stru, sim_stru, sim_inter = solve_experiment(
            data_stru,
            ModelOptions(model_family=WorkflowFamily.DATA2, mode=mode, theta=base_parameters, sim_opt=False, B_form=b_form, sigma_fixed=True, LOUD=False),
        )
        plot_sim_comparison(data_stru, sim_stru, stirc_mass=False, plot_pred=True, lg=False, LOUD=True)
        cases.append({'dataset': str(ds), 'fit': fit_stru})
    return {'output_dir': str(outdir), 'cases': cases}


def run_data2_model_variations(
    *,
    output_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Bundle the DATA2 model-variation/FIM script into a function."""
    outdir = _ensure_output_dir(output_dir)
    results: List[Dict[str, Any]] = []

    def _run_case(dataset: Union[str, Path], mode: str, theta: Dict[str, float], label: str, fim_formula: str = 'Backward'):
        data_stru = loadmat(str(dataset))['data_stru']
        fit_stru, sim_stru, sim_inter = solve_experiment(
            data_stru,
            ModelOptions(model_family=WorkflowFamily.DATA2, mode=mode, theta=theta, sim_opt=False, B_form=1, LOUD=False, sigma_fixed=False),
        )
        fit_path = outdir / f'{label}-fit.json'
        store_json(str(fit_path), fit_stru)
        doe_stru = calc_FIM(
            data_stru,
            mode,
            theta=fit_stru['parameters'],
            step=1e-8,
            formula=fim_formula,
            B_form=1,
            model_builder=resolve_model_builder(WorkflowFamily.DATA2),
        )
        fim_path = outdir / f'{label}-FIM.json'
        store_json(str(fim_path), doe_stru)
        results.append({'label': label, 'fit': fit_stru, 'fit_path': str(fit_path), 'fim_path': str(fim_path)})

    base_lag = {'Lp': 11, 'beta_c': 15, 'beta_0': 1, 'beta_1': 0.01, 'sigma': 1.0, 'S0': 0}
    _run_case('data_library/data_stru-dataset270511.123.mat', 'Lag', base_lag, 'LagM1')
    _run_case('data_library/data_stru-dataset270511.122.mat', 'Lag', {'Lp': 11, 'beta_c': 15, 'beta_0': 1, 'beta_1': 0.01, 'sigma': 1.0, 'S0': 0}, 'LagM2')
    _run_case('data_library/data_stru-dataset270511.121.mat', 'DATA', {'Lp': 11, 'beta_c': 15, 'beta_0': 1, 'beta_1': 0.01, 'sigma': 1.0, 'S0': 0}, 'LagM3')
    _run_case('data_library/data_stru-dataset270511.12.mat', 'DATA', {'Lp': 11, 'beta_c': 15, 'beta_0': 1, 'beta_1': 0.01, 'sigma': 1.0, 'S0': 0}, 'LagM4')
    base_overflow = {'Lp': 11, 'beta_c': 15, 'beta_0': 1, 'beta_1': 0.01, 'sigma': 1.0, 'S0': -0.1756665334051245}
    _run_case('data_library/data_stru-dataset270511.423.mat', 'Overflow', base_overflow, 'OverflowM1')
    _run_case('data_library/data_stru-dataset270511.422.mat', 'Overflow', base_overflow, 'OverflowM2')
    _run_case('data_library/data_stru-dataset270511.421.mat', 'DATA', base_overflow, 'OverflowM3')
    _run_case('data_library/data_stru-dataset270511.42.mat', 'DATA', base_overflow, 'OverflowM4')
    return {'output_dir': str(outdir), 'cases': results}


def run_data2_pre_b_dependence(
    *,
    dataset_paths: Optional[Sequence[Union[str, Path]]] = None,
    mode: str = 'DATA',
    output_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Bundle the DATA2 pre-B dependence investigation into a function."""
    outdir = _ensure_output_dir(output_dir)
    dataset_paths = list(dataset_paths) if dataset_paths is not None else [
        ('A.0', 'data_library/data_stru-dataset270511.12.mat'),
        ('A.1', 'data_library/data_stru-dataset270611.12.mat'),
        ('A.2', 'data_library/data_stru-dataset270711.12.mat'),
        ('B.1', 'data_library/data_stru-dataset270511.22.mat'),
        ('B.2', 'data_library/data_stru-dataset270511.32.mat'),
        ('C.1', 'data_library/data_stru-dataset270511.42.mat'),
        ('C.2', 'data_library/data_stru-dataset270511.92.mat'),
        ('D.1', 'data_library/data_stru-dataset270511.52.mat'),
        ('D.2', 'data_library/data_stru-dataset270511.62.mat'),
        ('E.1', 'data_library/data_stru-dataset270511.72.mat'),
        ('E.2', 'data_library/data_stru-dataset270511.82.mat'),
    ]
    all_sim: Dict[str, Any] = {}

    if dataset_paths and isinstance(dataset_paths[0], tuple):
        items = dataset_paths  # type: ignore[assignment]
    else:
        items = [(Path(p).stem, p) for p in dataset_paths]

    for label, data_file in items:
        data_stru = loadmat(str(data_file))['data_stru']
        _, sim_stru, _ = solve_experiment(
            data_stru,
            ModelOptions(model_family=WorkflowFamily.DATA2, mode=mode, sim_opt=False, B_form='pervial', sigma_fixed=False),
        )
        all_sim[label] = sim_stru

    cmap1 = plt.get_cmap('tab10')
    cmap2 = plt.get_cmap('tab20')
    fig = plt.figure(figsize=(6, 4))
    fi = -1
    for key, sim_st in all_sim.items():
        if sim_st[0]['B'] is not None:
            co = cmap1(10) if fi == -1 else cmap2(fi)
            plt.plot([], [], color=co, linewidth=3, alpha=.8, label=key)
            for i in sim_st:
                if i == 0:
                    plt.plot(np.array(sim_st[i]['cIn'][80:]), np.array(sim_st[i]['Js'][80:]) / np.array(sim_st[i]['Jw'][80:]), color=co, linewidth=2, alpha=.8)
                else:
                    plt.plot(np.array(sim_st[i]['cIn'][50:]), np.array(sim_st[i]['Js'][50:]) / np.array(sim_st[i]['Jw'][50:]), color=co, linewidth=2, alpha=.8)
            fi += 1
    plt.xlabel('Interface Concentration[mM]', fontsize=16, fontweight='bold')
    plt.ylabel('Js/Jw [mM]', fontsize=16, fontweight='bold')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction='in')
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(fontsize=10, loc='best')
    fig.savefig(outdir / 'Js_Jw_cin.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return {'output_dir': str(outdir), 'cases': list(all_sim.keys())}


def data2_function_index() -> Dict[str, List[str]]:
    """Return a structured index of the DATA2 notebook and script helpers."""
    return {
        'visualization_helpers': [
            'plot_data2_pressure_change',
            'data2_calibration_curve',
            'data2_model_predictions',
            'data2_model_convection',
            'plot_data2_model_predictions',
            'data2_linear_regression',
            'plot_data2_regression_sensitivity',
        ],
        'workflow_wrappers': [
            'run_data2_model_demo',
            'run_data2_cross_verification',
            'run_data2_model_variations',
            'run_data2_pre_b_dependence',
        ],
        'compatibility_helpers': [
            'plot_sim_comparison',
            'plot_sim_comparison_paper',
            'plot_sim_comparison_legacy',
        ],
        'legacy_sources': [
            'DATA2_model_demo.ipynb',
            'DATA2_visualization.ipynb',
            'run_cross_verification.py',
            'run_DATA2_model_variations.py',
            'run_pre_B_dependence.py',
        ],
    }


def run_data_analysis(
    datasets: Sequence[Union[str, Path, float, int]],
    *,
    repo_root: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Python wrapper for the top-level DATA1 run_data_analysis.m workflow.

    The function currently focuses on loading, plotting, and returning a
    structured execution summary so the refactor can grow around it.
    """
    root = Path(repo_root).expanduser().resolve() if repo_root is not None else Path.cwd()
    outdir = _ensure_output_dir(output_dir or (root / 'figures'))
    summary: Dict[str, Any] = {'datasets': [], 'output_dir': str(outdir)}
    for item in datasets:
        ds = str(item)
        mat_path = root / 'legacy' / 'data1_matlab' / 'data' / f'data_stru-dataset{ds}.mat'
        if not mat_path.exists():
            mat_path = root / 'legacy' / 'data1_matlab' / 'data' / f'dat {ds} oneCPNT holdup concpolar cvmv fixed ch0' / f'data_stru-dataset{ds}.mat'
        if not mat_path.exists():
            raise FileNotFoundError(f"Could not locate dataset MAT file for {ds}")
        data_stru = loadmat(str(mat_path))['data_stru']
        summary['datasets'].append({'dataset': ds, 'mat_path': str(mat_path), 'loaded': True})
        # The plotting helpers can be called once fit/sim structures exist.
        # Here we return the load plan and inputs so downstream code can drive it.
    return summary


def run_sigma_sensitivity(
    dataset: Union[str, float, int],
    sigma: Sequence[float],
    theta: Sequence[float],
    *,
    repo_root: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Python wrapper for the DATA1 sigma sensitivity MATLAB workflow."""
    root = Path(repo_root).expanduser().resolve() if repo_root is not None else Path.cwd()
    outdir = _ensure_output_dir(output_dir or (root / 'sigma_sensitivity'))
    ds = str(dataset)
    mat_path = root / 'legacy' / 'data1_matlab' / 'data' / f'data_stru-dataset{ds}.mat'
    if not mat_path.exists():
        raise FileNotFoundError(f"Could not locate dataset MAT file for {ds}")
    data_stru = loadmat(str(mat_path))['data_stru']
    return {
        'dataset': ds,
        'mat_path': str(mat_path),
        'sigma': list(sigma),
        'theta': list(theta),
        'output_dir': str(outdir),
        'data_config_keys': list(data_stru['data_config'].keys()) if isinstance(data_stru, dict) and 'data_config' in data_stru else [],
    }


def do_data1_notebook_figures(
    *,
    repo_root: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Placeholder orchestration for the DATA1 notebook figure set.

    The function is intentionally data-driven: it returns the file plan for the
    notebook-based figures so future refactor steps can wire actual fit outputs
    into these plotting helpers.
    """
    root = Path(repo_root).expanduser().resolve() if repo_root is not None else Path.cwd()
    outdir = _ensure_output_dir(output_dir or (root / 'figures'))
    return {
        'output_dir': str(outdir),
        'figures': [
            'mass-dat501.1.png',
            'mass-dat501.11.png',
            'mass-dat511.12.png',
            'concentration-dat501.1.png',
            'concentration-dat501.11.png',
            'concentration-dat511.12.png',
            'contour_fixsig-mass.png',
            'contour_fixsig-permeate_conc.png',
            'contour_fixsig-retentate_conc.png',
            'sigma_sensitivity-mass.png',
            'sigma_sensitivity-reten_conc.png',
            'sigma_sensitivity-perme_conc.png',
            'concentration_range.png',
            'calib_curve.png',
            'cr_measure-dat511.12.png',
        ],
    }


def data1_matlab_function_index() -> Dict[str, List[str]]:
    """Return a structured index of the MATLAB scripts/functions that were ported."""
    return {
        'top_level_scripts': [
            'run_data_analysis',
            'run_sigma_sensitivity',
            'doe_heatmap_diafiltration',
            'doe_heatmap_filtration',
            'heatmap_sigma_sensitivity_diafiltration',
            'heatmap_sigma_sensitivity_filtration',
        ],
        'plotting_helpers': [
            'plot_sim_comparison',
            'plot_contour',
            'plot_sim',
            'plot_sim_show',
            'plot_contour_sig_sen',
            'plot_conc_range',
            'plot_cr_measure',
            'calib_curve_cond',
            'plot_conc_comparison',
        ],
        'legacy_support_functions': [
            'loadmat',
            'load_experiment_bundle',
            'load_experiment_files',
            'load_experiment_file',
            'detect_source_type',
        ],
    }


def loadmat(filename):
    '''
    Read in nested structure(mat file) generated from MATLAB and output dictionaries.

    This function is called instead of using spio.loadmat directly as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries which are still mat-objects.
    Adapted from https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
    
    Arguments:
        filename: str, location + filename
    
    Returns:
        dictionary contains structured data
    '''

    print("\nLoading data file =",filename,"\n")

    def _check_keys(d):
        '''
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        '''
        for key in d:
            if isinstance(d[key], spio.matlab.mat_struct):
                d[key] = _todict(d[key])
            elif isinstance(d[key], np.ndarray):
                d[key] = _tolist(d[key])

        return d

    def _todict(matobj):
        '''
        A recursive function which constructs from matobjects nested dictionaries
        '''
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, spio.matlab.mat_struct):
                d[strg] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                d[strg] = _tolist(elem)
            else:
                d[strg] = elem
        return d

    def _tolist(ndarray):
        '''
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        '''
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, spio.matlab.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(_tolist(sub_elem))
            elif isinstance(sub_elem, int):
                elem_list.append(float(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    
    return _check_keys(data)


def plot_sim_comparison_legacy(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=False):
    '''
    Plot simulation results comparing with measurements
    
    Arguments:
        data_stru: dict, experimental data dictionary
        sim_stru: dict, simulation results dictionary
        stirc_mass: boolean, if plot mass change in stirred cell
        plot_pred: boolean, if plot model predictions/simualtions
        lg: boolean, if plot legends
        LOUD: boolean, if store the figures
    
    Actions:
        create plots and store (optional)
    '''
    t_delay = data_stru['data_raw'][0]['time'][0]
    vial1 = data_stru['data_config']['n_v0']-1
    t_start = 0

    # plot mass data/prediction comparison
    fig = plt.figure(figsize=(4,4))
    for i in range(data_stru['data_config']['n']):
        plt.plot([(float(t)-t_delay-t_start)/60 for t in data_stru['data_raw'][i]['time']],
                 data_stru['data_raw'][i]['mass'],'r.',markersize=4)
        if plot_pred and i >= data_stru['data_config']['n_v0']-1:
            plt.plot([(time-t_start)/60 for time in sim_stru[i]['time']],
                     sim_stru[i]['mV'],'b',linewidth=2,
                     alpha=.6)

    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Mass in Vial [g]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.ylim(bottom=0)
    #plt.show()
    
    if LOUD:
        fname = 'figures/mass-dat'+str(data_stru['dataset'])
        fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')

    # plot concentration data/prediction comparison
    fig = plt.figure(figsize=(4,4))

    for i in range(data_stru['data_config']['n']):
        if type(data_stru['data_raw'][i]['cV_avg']) != list:
            plt.plot((float(data_stru['data_raw'][i]['time'][-1])-t_delay-t_start)/60,
                     data_stru['data_raw'][i]['cV_avg'],
                     'cs',markersize=6)
        else:
            plt.plot([(float(t)-t_delay-t_start)/60 for t in data_stru['data_raw'][i]['time']],
                     data_stru['data_raw'][i]['cV_avg'],'cs',markersize=6,clip_on=False)
        plt.plot([(float(t)-t_delay-t_start)/60 for t in data_stru['data_raw'][i]['time']],
                 data_stru['data_raw'][i]['cF_exp'],'ms',markersize=6,clip_on=False)
        if plot_pred:
            plt.plot([(time-t_start)/60 for time in sim_stru[i]['time']],
                     sim_stru[i]['cF'],
                     'g',linewidth=2)
            plt.plot([(time-t_start)/60 for time in sim_stru[i]['time']],
                     sim_stru[i]['cH'],
                     'r-',linewidth=2,alpha=.6)
            if type(data_stru['data_raw'][i]['cV_avg']) != list:
                plt.plot((sim_stru[i]['time'][-1]-t_start)/60,
                         sim_stru[i]['cV'][-1],
                         'r^',linewidth=2,alpha=.6)
            else:
                time = data_stru['data_raw'][i]['time']-t_delay
                index = ~np.isnan(data_stru['data_raw'][i]['cV_avg'])
                t = [time[i] for i, val in enumerate(index) if val]
                f = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['cV'], fill_value='extrapolate')  
                plt.plot([(ti-t_start)/60 for ti in t],
                         f(t),
                         'r^',linewidth=2,alpha=.6)
                
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Concentration [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in",top=True, right=True)
    plt.ylim(bottom=0)
    #plt.show()
    
    if LOUD:
        fname = 'figures/concentration-dat'+str(data_stru['dataset'])
        fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
        
    # plot mass of stirred cell
    if stirc_mass:
        fig = plt.figure(figsize=(4,4))
        for i in range(data_stru['data_config']['n']):
            if plot_pred:
                plt.plot([(time-t_start)/60 for time in sim_stru[i]['time']],
                         sim_stru[i]['mF'],'b',linewidth=2,
                         alpha=.6)

        # ghost point for legend
        plt.plot([],[],'r.',markersize=4,label='Mass (Measurements)')
        if plot_pred:
            plt.plot([],[],'b',linewidth=2,alpha=.6,label='Mass (Predictions)')

        plt.plot([],[],'ms',markersize=6,clip_on=False,label='Retentate (Measurements)')
        if plot_pred:
            plt.plot([],[],'g',linewidth=2,label='Retentate (Predictions)')
        if type(data_stru['data_raw'][i]['cV_avg']) != list:
            plt.plot([],[],'cs',markersize=6,label='Vial (Measurements)')
        else:
            plt.plot([],[],'cs',markersize=6,label='Vial (Measurements)')
        if plot_pred:
            plt.plot([],[],'r^',markersize=6,label='Vial (Predictions)')        
            plt.plot([],[],'r-',linewidth=2,alpha=.6,label='Permeate (Predictions)')

        plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
        plt.ylabel('Mass in Stirred Cell [g]',fontsize=16,fontweight='bold')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tick_params(direction="in")
        ytop = max([max(sim_stru[i]['mF']) for i in range(data_stru['data_config']['n'])])
        plt.ylim(bottom=0,top=ytop+0.5)
        if lg:
            plt.legend(fontsize=12.5,loc='best')#bbox_to_anchor=(1.02, 0.3),borderaxespad=0,ncol=3)
        #plt.show()
        if LOUD:
            fname = 'figures/stirc_mass-dat'+str(data_stru['dataset'])
            fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')


def save_model(m, LO=True, B_form='single', LOUD=False):
    '''
    Save model results
    
    Arguments:
        m: pyomo model instance
        LO: boolean, if mode is lag/overflow
        B_form: float/str, different solute permeability coefficient B formula 
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration              
        LOUD: boolean, if print out parameter results
    
    Returns:
        fit_stru: dict, parameter fit results
        sim_stru: dict, model predictions
    '''
    fit_stru = dict()
    fit_stru['parameters'] = dict()
    fit_stru['sideparameters'] = dict()
    
    fit_stru['parameters']['Lp'] = value(m.Lp) 
    
    if B_form=='single':
        fit_stru['parameters']['B'] = value(m.B)
    elif B_form=='pervial':
        fit_stru['parameters']['B'] = dict()
        for i in m.n_vial:
            fit_stru['parameters']['B'][i] = value(m.B[i])
    elif B_form=='convection':   
        fit_stru['parameters']['beta_0'] = value(m.beta_0)
        fit_stru['parameters']['beta_1'] = value(m.beta_1)
    else:
        fit_stru['parameters']['beta_0'] = value(m.beta_0)
        if not isinstance(B_form, str) and B_form != 0:
            fit_stru['parameters']['beta_1'] = value(m.beta_1)
        if not isinstance(B_form, str) and B_form > 1:
            fit_stru['parameters']['beta_2'] = value(m.beta_2)
        if not isinstance(B_form, str) and B_form > 2:
            fit_stru['parameters']['beta_3'] = value(m.beta_3)
            
    fit_stru['parameters']['sigma'] = value(m.sigma)
    if LO:
        fit_stru['parameters']['S0'] = value(m.S0)
        fit_stru['parameters']['S'] = value(m.S)   
    fit_stru['Obj'] = value(m.Obj)
    fit_stru['obj_m'] = value(m.obj_m)
    fit_stru['obj_cv'] = value(m.obj_cv)
    fit_stru['obj_cr'] = value(m.obj_cr)
    fit_stru['Obj_tru'] = value(m.Obj_tru)
    fit_stru['obj_cr_tru'] = value(m.obj_cr_tru)
    fit_stru['res_std'] = {'res_m': [value(res) for res in m.res_m],
                 'res_cp': [value(res) for res in m.res_cp],
                 'res_cf': [value(res) for res in m.res_cf]}   
        
    if LOUD:
        # print parameters
        print('Lp = ',value(m.Lp),' L / m / m / hr / bar')
        if isinstance(B_form, str):
            if B_form=='single':
                print('B = ',value(m.B),' micrometers / s')
            elif B_form=='pervial':    
                print('B = ',[value(m.B[i]) for i in m.n_vial],' micrometers / s')
            elif B_form=='convection':   
                print('D/l = ',value(m.beta_0))
                print('H = ',value(m.beta_1))
        else:
            if B_form > 2:
                print('B = Jw*[',value(m.beta_0),'+',value(m.beta_1),'* c_in ^',B_form-2,'+',value(m.beta_2),'* c_in ^',B_form-1,'+',value(m.beta_3),'* c_in ^',B_form, ']')
            elif B_form > 1:
                print('B = Jw*[',value(m.beta_0),'+',value(m.beta_1),'* c_in ^',B_form-1,'+',value(m.beta_2),'* c_in ^',B_form,']')
            elif B_form == 0:
                print('B = Jw*[',value(m.beta_0), ']')
            else:
                print('B = Jw*[',value(m.beta_0),'+',value(m.beta_1),'* c_in ^',B_form,']')

        print('sigma = ',value(m.sigma),' dimensionless')
        print('cD = ',value(m.cD),' mM')

        if LO:
            print('S0 = ',value(m.S0),' g / s')
            print('S = ',value(m.S),' g / hr')
        print('Obj = ' ,value(m.Obj))
        print('Obj(m) = ' ,value(m.obj_m))
        print('Obj(cv) = ' ,value(m.obj_cv))
        print('Obj(cr) = ' ,value(m.obj_cr))
        if LO:
            print('Obj(truncate) = ' ,value(m.obj_tru))
            print('Obj(cr_truncate) = ' ,value(m.obj_cr_tru))
            print('count(cr0) = ' ,value(m.count_cr0))
        print('count = ' ,value(m.count))
        print('count(m) = ' ,value(m.count_m))
        print('count(cv) = ' ,value(m.count_cv))
        print('count(cr) = ' ,value(m.count_cr))
        print('likelihood1 = ' ,value(m.llh1))
        print('likelihood2 = ' ,value(m.llh2)) 
    
    sim_stru = dict()
    
    for n_vial in m.n_vial:
        time = [t * (m.tf[n_vial]-m.ti[n_vial]) + m.ti[n_vial] for t in m.tau]
        cF = [value(m.cF[n_vial,t]) for t in m.tau]
        cIn = [value(m.cIn[n_vial,t]) for t in m.tau]
        cH = [value(m.cH[n_vial,t]) for t in m.tau]
        mV = [value(m.mV[n_vial,t]) for t in m.tau]
        cV = [value(m.cV[n_vial,t]) for t in m.tau]
        Jw = [value(m.Jw[n_vial,t]) for t in m.tau]
        Js = [value(m.Js[n_vial,t]) for t in m.tau]
        
        sim_stru[n_vial-1] = {'time': time,
                'cF': cF,
                'cIn': cIn,
                'cH': cH,
                'mV': mV,
                'cV': cV,
                'Jw': Jw,
                'Js': Js}
        if LO:
            mF = [value(m.mF[n_vial,t]) for t in m.tau]
            sim_stru[n_vial-1]['mF'] = mF
        
        if isinstance(B_form, str) and 'convection' in B_form:
            B = None       
        elif B_form=='single':
            B = [value(m.B) for t in m.tau]
        elif B_form=='pervial':
            B = [value(m.B[n_vial]) for t in m.tau]
        else:
            B = [value(m.B[n_vial,t]) for t in m.tau]        
        sim_stru[n_vial-1]['B'] = B

    return fit_stru, sim_stru


def inter_model(data_stru, sim_stru, fit_stru):
    '''
    Interpolate model predictions at same time stamps of experimental measurements
    
    Arguments:
        data_stru: dict, experimental data dictionary
        sim_stru: dict, model predictions
        fit_stru: dict, parameter fit results
    
    Returns:
        sim_inter: dict, model predictions for experimental measurements 
    '''
    sim_inter = dict()
    for i in range(data_stru['data_config']['n']):
        t_delay = data_stru['data_raw'][0]['time'][0]
        t_meas = data_stru['data_raw'][i]['time'] - t_delay
        f_m = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['mV'], kind='linear', fill_value='extrapolate')
        if i >= data_stru['data_config']['n_v0']-1:            
            ind_m = ~np.isnan(data_stru['data_raw'][i]['mass'])
            t_m = [t_meas[i] for i, val in enumerate(ind_m) if val]
        else:
            t_m=[]
        if type(data_stru['data_raw'][i]['cV_avg']) == list:
            f_cv = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['cV'], kind='linear', fill_value='extrapolate')
            ind_cv = ~np.isnan(data_stru['data_raw'][i]['cV_avg'])
            t_cv = [t_meas[i] for i, val in enumerate(ind_cv) if val]
            cV = f_cv(t_cv)
        else:
            cV = sim_stru[i]['cV'][-1]
        f_cf = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['cF'], kind='linear', fill_value='extrapolate')
        ind_cf = ~np.isnan(data_stru['data_raw'][i]['cF_exp'])
        t_cf = [t_meas[i] for i, val in enumerate(ind_cf) if val]
        sim_inter[i] = {'time': t_meas,
                        'mV': f_m(t_m),
                        'cF': f_cf(t_cf),
                        'cV': cV}       

    return sim_inter


def model_construct_inter(data_stru, mode, theta=None, sim_opt=False, B_form='single'):
    '''
    Build diafiltration model in pyomo
    
    Arguments:
        data_stru: dict, experimental data dictionary
        mode: str, experiment mode, {DATA, lag, overflow}
        theta: dict, preset parameter values
        sim_opt: boolean, if run simulation with fixed parameter values
        B_form: float/str, different solute permeability coefficient formula
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration 
    
    Returns:
        m: pyomo model instance
    '''
    # known parameters
    # applied pressure
    delP = data_stru['data_config']['delP'] #[bar]
    # gas constant
    R = 8.314e-5 #[cm^3 bar / micromol / K]
    # temperature
    T = data_stru['data_config']['Temp'] #[K]
    # membrane area
    Am = data_stru['data_config']['Am'] #[cm^2]
    # density
    rho = data_stru['data_config']['rho'] #[g/cm^3]
    # number of dissolved species
    ni = data_stru['data_config']['ni']
    # kinematic viscosity
    nu = 8.927e-3 #[cm^2/s]
    # diameter of stirred cell
    b = 2.2860 #[cm]
    # average velocity within the system
    v = 350/60 * np.pi * b #[cm/s]
    # diffusion coefficient
    if isinstance(data_stru['data_config']['namec'], str) and 'K' in data_stru['data_config']['namec']:
        D = 1.960e-5 #[cm^2/s]  - K+
    else:
        raise NotImplementedError
    # mass transfer coefficien
    k = 0.23 * v**0.57 * D**0.67 / (nu**0.24 * b**0.43)

    # known constants
    t_delay = data_stru['data_raw'][0]['time'][0]
    M_F0 = data_stru['data_config']['M_F0']
    M_O = data_stru['data_config']['M_O']
    cD = data_stru['data_config']['C_D']
    mH = 0.25 #ml
    C_F0 = data_stru['data_config']['C_F0']
    C_V0 = 1e-6
    
    N_VIAL = data_stru['data_config']['n']
    N_V0 = data_stru['data_config']['n_v0']
    if mode !='DATA':
        N_extra = data_stru['data_config']['n_extra']
        N_H = data_stru['data_config']['n_h']
        N_A = data_stru['data_config']['n_A']       
        if mode == 'Overflow':
            N_A0 = 1
            
    Tauf = 1 #scaled ending time

    # create a model object
    m = ConcreteModel()

    # define the independent variable
    m.n_vial = RangeSet(N_VIAL) # number of vials
    m.tau = ContinuousSet(bounds=(0, Tauf))#scaled time
    
    TF_list = [data_stru['data_raw'][i]['time'][-1]-t_delay for i in range(N_VIAL)]
    TF_dict = dict(zip(m.n_vial,TF_list)) # unscaled time elapse for each vial 
    m.tf = Set(initialize=TF_list)
    
    TI_list = [data_stru['data_raw'][i]['time'][0]-t_delay for i in range(N_VIAL)]
    TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time per vial

    # parameter initialization
    param_in = dict()
    if theta is None:
        param_in['Lp'] = data_stru['data_config']['Lp0']
        param_in['B'] = data_stru['data_config']['B0']
        param_in['sigma'] = data_stru['data_config']['sigma0']
        param_in['beta_0'] = param_in['B']*36000/param_in['Lp']/delP #param_in['B']
        param_in['beta_1'] = 1
        param_in['beta_2'] = 0
        param_in['beta_3'] = 0
        if mode == 'Lag':
            param_in['S0'] = 0      
        elif mode == 'Overflow':
            param_in['S0'] = -M_O/10 *3600
    else:
        param_in = theta
    # define model parameters
    if sim_opt:        
        m.cD = cD
        m.C_H0 = 1e-6
        #simulation parameters
        m.Lp = Param(initialize=param_in['Lp'],mutable=True)
        
        if B_form=='single':
            m.B = Param(initialize=param_in['B'],mutable=True)
        elif B_form=='pervial':    
            m.B = Param(m.n_vial,initialize=param_in['B'],mutable=True)
        elif B_form=='convection':    
            m.beta_0 = Param(initialize=param_in['beta_0'])
            m.H = Var(m.n_vial, m.tau,initialize=0.5)    
        else:
            m.beta_0 = Param(initialize=param_in['beta_0'],mutable=True)
            if not isinstance(B_form, str):
                if B_form != 0:
                    m.beta_1 = Param(initialize=param_in['beta_1'],mutable=True)
                if B_form > 1:
                    m.beta_2 = Param(initialize=param_in['beta_2'],mutable=True)
                if B_form > 2:
                    m.beta_3 = Param(initialize=param_in['beta_3'],mutable=True) 
                m.B = Var(m.n_vial, m.tau, bounds=(1e-6,50))
        m.sigma = Param(initialize=param_in['sigma'],mutable=True)

    else: 
        #parameter estimation
        m.C_H0 = 1e-6
        m.cD = cD
        m.Lp = Var(bounds=(0.5,50),initialize=param_in['Lp'])
        if B_form=='single':
            m.B = Var(bounds=(1e-6,30),initialize=param_in['B'])
        elif B_form=='pervial':    
            m.B = Var(m.n_vial, bounds=(1e-6,30),initialize=param_in['B'])
        elif isinstance(B_form, str) and 'convection' in B_form:
            m.beta_0 = Var(bounds=(1+1e-6,50),initialize=param_in['beta_0'])
            m.H = Var(m.n_vial, m.tau,initialize=0.5)  
            m.beta_1 = Var(bounds=(0.0,1.0),initialize=0.5)         
        else:    
            m.beta_0 = Var(initialize=param_in['beta_0'])
            m.B = Var(m.n_vial, m.tau)
            if B_form != 0:
                m.beta_1 = Var(bounds=(-20,20),initialize=param_in['beta_1'])
            if B_form > 1:
                m.beta_2 = Var(bounds=(-20,20),initialize=param_in['beta_2'])
            if B_form > 2:
                m.beta_3 = Var(bounds=(-20,20),initialize=param_in['beta_3'])
                
        m.sigma = Var(bounds=(0.0,1.0), initialize=param_in['sigma'])
    
    if mode == 'Lag':
        if sim_opt:
            m.S0 = Param(initialize=param_in['S0'],mutable=True) 
            m.S = Param(initialize=param_in['S'],mutable=True)
        else:
            # g/s
            m.S0 = Param(initialize=param_in['S0'])
            # g/hr
            m.S = Var(initialize=M_O/sum(m.tf-m.ti)*3600)
    elif mode == 'Overflow':
        if sim_opt:
            m.S0 = Param(initialize=param_in['S0'],mutable=True)
            m.S = Param(initialize=param_in['S'],mutable=True)
        else:
            # g/s
            m.S = Var(initialize=0)
            # g/hr
            m.S0 = Var(initialize=-M_O/10)
            
    # define dependent variables 
    if mode !='DATA':
        m.mF = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=M_F0)
    m.cF = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=C_F0)
    m.cIn = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=C_F0)
    m.cH = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=1e-6)
    m.mV = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=1e-6)
    m.cVmV = Var(m.n_vial, m.tau,initialize=C_V0 * 1e-6)
    m.cV = Var(m.n_vial,  m.tau, domain=NonNegativeReals,initialize=1e-6)

    # define intermediate variables
    m.Jw = Var(m.n_vial, m.tau)
    m.Js = Var(m.n_vial, m.tau)
    if isinstance(B_form, str) and 'convection' in B_form:
        m.Js_exp = Var(m.n_vial, m.tau, bounds=(1+1e-6,1e4))
    elif B_form=='K':   
        m.Js_exp = Var(m.n_vial, m.tau)

    # define derivatives
    if mode !='DATA':
        m.dmF = DerivativeVar(m.mF,wrt=m.tau)
    m.dcF = DerivativeVar(m.cF,wrt=m.tau)
    m.dcH = DerivativeVar(m.cH,wrt=m.tau)
    m.dmV = DerivativeVar(m.mV,wrt=m.tau)
    m.dcVmV = DerivativeVar(m.cVmV,wrt=m.tau)
    
    # define the differential equation as constraints
    def ode_mF_rule(m, n, t):
        if mode == 'Overflow':
            if n < N_A0:
                # (dmF_dt = 0 - Jw * Am * rho) * tf
                return m.dmF[n,t] == (0 - Am * rho * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
            elif N_A0-1 < n <= N_A:
                # (dmF_dt = -S0 - Jw * Am * rho) * tf
                return m.dmF[n,t] == (- m.S0  - Am * rho * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
            else:
                # (dmF_dt = S) * tf
                return m.dmF[n,t] == m.S / 3600  * (TF_dict[n]-TI_dict[n])/Tauf
        elif mode == 'Lag':
            if n <= N_A:
                # (dmF_dt = -S0 - Jw * Am * rho) * tf
                return m.dmF[n,t] == (- m.S0 - Am * rho * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf   
            else:
                # (dmF_dt = S) * tf
                return m.dmF[n,t] == m.S / 3600  * (TF_dict[n]-TI_dict[n])/Tauf
    if mode !='DATA':
        m.ode_mF = Constraint(m.n_vial, m.tau, rule=ode_mF_rule)
     
    def ode_cF_rule(m, n, t):
        if mode == 'Overflow':
            if n < N_A0:
                return m.dcF[n,t] == 1 / m.mF[n,t] * (Am * rho * (m.cF[n,t] * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf        
            elif N_A0-1 < n <= N_A:  
                # (dcF/dt = (cF - cD) * S0 / mF + Am * rho / mF * (cF * Jw - Js)) * tf
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cF[n,t] - m.cD ) * m.S0 + Am * rho * (m.cF[n,t] * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf                     
            else:
                # (dcF/dt = (cD - cF) * S / mF + Am * rho / mF * (cD * Jw - Js)) * tf
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cD - m.cF[n,t]) * m.S / 3600 + Am * rho * (m.cD * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf
        elif mode == 'Lag':
            if n <= N_A:
                # (dcF/dt = (cF - cD) * S0 / mF + Am * rho / mF * (cF * Jw - Js)) * tf
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cF[n,t] - m.cD ) * m.S0 + Am * rho * (m.cF[n,t] * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf                     
            else:
                # (dcF/dt = (cD - cF) * S / mF + Am * rho / mF * (cD * Jw - Js)) * tf
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cD - m.cF[n,t]) * m.S / 3600 + Am * rho * (m.cD * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf
        elif mode == 'DATA':
            return m.dcF[n,t] == Am * rho / M_F0 * (m.cD * m.Jw[n,t] - m.Js[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf

    m.ode_cF = Constraint(m.n_vial, m.tau, rule=ode_cF_rule)

    # (dcH_dt = Am * rho / mH *(Js - Jw * cH)) * tf
    def ode_cH_rule(m, n, t):
        if mode !='DATA' and n <= N_H:
            return m.dcH[n,t] == Am * rho / m.mV[n,t] * (m.Js[n,t] - m.cH[n,t] * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
        else:
            return m.dcH[n,t] == Am * rho / mH * (m.Js[n,t] - m.cH[n,t] * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_cH = Constraint(m.n_vial, m.tau, rule=ode_cH_rule)
    
    # (dmV_dt = Jw * Am * rho) * tf
    def ode_mV_rule(m, n, t):
        return m.dmV[n,t] == m.Jw[n,t] * Am * rho * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_mV = Constraint(m.n_vial, m.tau, rule=ode_mV_rule)
    
    # (dmVcV_dt = Jw * Am * rho * cH) * tf
    def ode_cVmV_rule(m, n, t):
        if mode !='DATA' and n <= N_H:
            return m.dcVmV[n,t] == Am * rho * m.Js[n,t] * (TF_dict[n]-TI_dict[n])/Tauf
        else:                
            return m.dcVmV[n,t] == m.Jw[n,t] * m.cH[n,t] * Am * rho * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_cVmV = Constraint(m.n_vial, m.tau, rule=ode_cVmV_rule)
    
    # intermediate variables constraints
    # mF(end) - mF(0) = m_o
    def eqn_S_rule(m):
        return m.mF[m.n_vial.last(),Tauf] - M_F0 == M_O
    if mode !='DATA' and not sim_opt:
        m.eqn_S = Constraint(rule=eqn_S_rule)    
    # cIn = (cF - cH) .* exp(Jw./ k) + cH
    def eqn_cIn_rule(m, n, t):
        return m.cIn[n,t] == (m.cF[n,t] - m.cH[n,t]) * exp(m.Jw[n,t]/ k) + m.cH[n,t]
    m.eqn_cIn = Constraint(m.n_vial, m.tau, rule=eqn_cIn_rule)
    
    # Jw = Lp*(delP - (cIn - cH).' *(ni.*sigma)*R*T)
    # cm/s
    def eqn_Jw_rule(m, n, t):
        return m.Jw[n,t]*36000 == m.Lp *(delP - (m.cIn[n,t] - m.cH[n,t])*ni*m.sigma*R*T)
    m.eqn_Jw = Constraint(m.n_vial, m.tau, rule=eqn_Jw_rule)
    
    # Js = B*(cIn - cH)
    # micromole/cm^2/s
    def eqn_Js_rule(m, n, t):
        if B_form=='single':
            return m.Js[n,t]*10000 == m.B * (m.cIn[n,t] - m.cH[n,t])
        elif B_form=='pervial':    
            return m.Js[n,t]*10000 == m.B[n] * (m.cIn[n,t] - m.cH[n,t])
        elif B_form=='convection':
            return m.Js[n,t] == m.Jw[n,t] * m.H[n,t] * (m.cIn[n,t] * m.Js_exp[n,t] - m.cH[n,t]) / (m.Js_exp[n,t]-1)
        else:
            #return m.Js[n,t]*10000 ==  m.B[n,t] * (m.cIn[n,t] - m.cH[n,t]) # No Jw in Js
            return m.Js[n,t]*10000 ==  (m.Jw[n,t] * m.B[n,t] * 10000) * (m.cIn[n,t] - m.cH[n,t]) #Jw in Js
            
    m.eqn_Js = Constraint(m.n_vial, m.tau, rule=eqn_Js_rule)
    
    def eqn_Js_exp_rule(m, n, t):
        if B_form=='convection':
            return m.Js_exp[n,t]  == exp(m.Jw[n,t]/m.beta_0*10000)   
        else:
            return Constraint.Skip
    m.eqn_Js_exp = Constraint(m.n_vial, m.tau, rule=eqn_Js_exp_rule)
    
    def eqn_H_rule(m,n,t):
        if B_form=='convection':
            return m.H[n,t]  == m.beta_1
        else:
            return Constraint.Skip
    m.eqn_H = Constraint(m.n_vial, m.tau, rule=eqn_H_rule)
    
    # cV = (mV*cV) /mV
    def eqn_cV_rule(m, n, t):
        return m.mV[n,t] * m.cV[n,t] == m.cVmV[n,t]
    m.eqn_cV = Constraint(m.n_vial, m.tau, rule=eqn_cV_rule)# numercally better without division
    
    # link the variables from different vials 
    def mF_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        else:
            return m.mF[n,Tauf]==m.mF[n+1,0]  
    if mode !='DATA':
        m.mF_linking = Constraint(m.n_vial,rule=mF_linking_rule) 
    
    def cF_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        else:
            return m.cF[n,Tauf]==m.cF[n+1,0]                                                                  
    m.cF_linking = Constraint(m.n_vial,rule=cF_linking_rule)   

    def cH_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        else:
            return m.cH[n,Tauf]==m.cH[n+1,0]                                                                  
    m.cH_linking = Constraint(m.n_vial,rule=cH_linking_rule)      

    def mV_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        elif mode !='DATA'and N_H < N_V0 and n < N_H:
            return m.mV[n,Tauf]==m.mV[n+1,0]
        elif mode !='DATA'and N_H < N_V0 and N_H < n < N_V0-1:
            return m.mV[n,Tauf]==m.mV[n+1,0]
        elif mode !='DATA'and N_V0-1 < N_extra and N_V0 <= n < N_extra+1:
            return m.mV[n,Tauf]==m.mV[n+1,0]
        else:
            return 1e-6==m.mV[n+1,0]                                                                  
    m.mV_linking = Constraint(m.n_vial,rule=mV_linking_rule) 
    
    def cVmV_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        elif mode !='DATA'and N_H < N_V0 and n < N_H:
            return m.cVmV[n,Tauf]==m.cVmV[n+1,0]
        elif mode !='DATA'and N_H < N_V0 and N_H < n < N_V0-1:
            return m.cVmV[n,Tauf]==m.cVmV[n+1,0]
        elif mode !='DATA'and N_V0-1 < N_extra and N_V0 <= n < N_extra+1:
            return m.cVmV[n,Tauf]==m.cVmV[n+1,0]
        else:
            return m.cV[n,Tauf]*1e-6==m.cVmV[n+1,0]                                                                  
    m.cVmV_linking = Constraint(m.n_vial,rule=cVmV_linking_rule)  
    
   # B per vial, skip startup vials
    def B_form_rule1(m,n):
        if n < N_V0-1:
            return m.B[n] == m.B[n+1]
        else:
            return Constraint.Skip

    # B-cf dependence
    def B_form_rule2(m,n,t):
        if not isinstance(B_form, str):
            if B_form > 2:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cF[n,t]**(B_form - 2) + m.beta_2 * m.cF[n,t]**(B_form - 1) + m.beta_3 * m.cF[n,t]**B_form #+ m.beta_4 * m.cH[n,t]        
            elif B_form > 1:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cF[n,t]**(B_form - 1) + m.beta_2 * m.cF[n,t]**B_form #+ m.beta_4 * m.cH[n,t]
            else:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cF[n,t]**B_form
        else:
            return Constraint.Skip

    # B-cin dependence
    def B_form_rule3(m,n,t):
        if not isinstance(B_form, str):
            if B_form > 2:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cIn[n,t]**(B_form - 2) + m.beta_2 * m.cIn[n,t]**(B_form - 1) + m.beta_3 * m.cIn[n,t]**B_form# + m.beta_4 * m.cH[n,t]        
            elif B_form > 1:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cIn[n,t]**(B_form - 1) + m.beta_2 * m.cIn[n,t]**B_form# + m.beta_4 * m.cH[n,t]
            elif B_form == 0:
                return m.B[n,t] == m.beta_0
            elif B_form < 0:
                return m.B[n,t] * m.cIn[n,t]**(-B_form) == m.beta_0 * m.cIn[n,t]**(-B_form) + m.beta_1
            else:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cIn[n,t]**B_form
        else:
            return Constraint.Skip
    
    if B_form=='pervial'and not sim_opt: 
        m.b_form = Constraint(m.n_vial,rule=B_form_rule1)
        pass
    else:
        m.b_form = Constraint(m.n_vial,m.tau,rule=B_form_rule3) 

    #initial conditions
    def _init(m): 
        def firstNonNan(listfloats):
            for item in listfloats:
                if np.isnan(item) == False:
                    return item
        if mode !='DATA':
            yield m.mF[1,0]==data_stru['data_config']['M_F0']
            yield m.cH[1,0]==m.C_H0 #1e-6
            yield m.cF[1,0]==firstNonNan(data_stru['data_raw'][0]['cF_exp'])
        else:
            if type(data_stru['data_raw'][0]['cV_avg']) == list:
                C_H0 = firstNonNan(data_stru['data_raw'][0]['cV_avg'])
            else:
                C_H0 = data_stru['data_raw'][0]['cV_avg']
            yield m.cH[1,0]==C_H0 * 0.8
            yield m.cF[1,0]==firstNonNan(data_stru['data_raw'][0]['cF_exp'])
               
        yield m.cVmV[1,0]==1e-6*1e-6
        yield m.mV[1,0]==1e-6
        yield ConstraintList.End
        
    m.con_boundary = ConstraintList(rule=_init)
        
    return m


def solve_model(data_stru, mode, theta=None, sim_opt=False, B_form='single', LOUD=False, model_builder=None):
    """
    Solve pyomo model
    
    Arguments:
        data_stru: dict, experimental data dictionary
        mode: str, experiment mode, {DATA, lag, overflow}
        theta: dict, preset parameter values
        sim_opt: boolean, if run simulation with fixed parameter values
        B_form: float/str, different solute permeability coefficient formula
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration 
        LOUD: boolean, if print out parameter results and store model predictions
    
    Returns:
        fit_stru: dict, parameter fit results
        sim_stru: dict, model predictions
        sim_inter: dict, model predictions for experimental measurements 
    """

    print("###################################################################")
    print("Creating and solving the Pyomo model with the following settings: ")
    print("mode =", mode)
    print("theta =",theta)
    print("sim_opt =", sim_opt)
    print("B_form =", B_form)
    print(" ")

    # interpolation function
    def interpolation(m,var,n_vial,t):
        """
        Interpolation function for pyomo model
        """
        tp = list(m.tau)
        for j in range(0,len(tp)-1):
            if tp[j+1]>=t and tp[j]<=t:
                var_inter = (var[n_vial,tp[j+1]]-var[n_vial,tp[j]]) * (t - tp[j])/(tp[j+1]-tp[j]) + var[n_vial,tp[j]]
        return var_inter
    
    # objective function
    def obj_rule(m):
        """
        pyomo Objective function (residual calculated using interpolation from model)
        """
        obj_m = 0
        obj_cp = 0
        obj_cf0 = 0
        obj_cf = 0
        ob_m = 0
        ob_cp = 0
        ob_cf0 = 0
        ob_cf = 0
        res_m_assemble = []
        res_cp_assemble = []
        res_cf_assemble = []
        
        collect_vial = data_stru['data_config']['n'] - data_stru['data_config']['n_extra']
        Count_m = 0
        Count_cp = 0
        Count_cf = 0
        count_cf0 = 0
        t_delay = data_stru['data_raw'][0]['time'][0]
        TF_list = [data_stru['data_raw'][i]['time'][-1]-t_delay for i in range(data_stru['data_config']['n'])]
        TF_dict = dict(zip(m.n_vial,TF_list)) # unscaled time elapse for each vial 
        TI_list = [data_stru['data_raw'][i]['time'][0]-t_delay for i in range(data_stru['data_config']['n'])]
        TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time for each vial
    
        for n_vial in m.n_vial:
            t_meas = data_stru['data_raw'][n_vial-1]['time'] - t_delay
            mv_meas = data_stru['data_raw'][n_vial-1]['mass']
            cp_meas = data_stru['data_raw'][n_vial-1]['cV_avg']
            cf_meas = data_stru['data_raw'][n_vial-1]['cF_exp']
    
            t_meas_scaled = [(t-TI_dict[n_vial])/(TF_dict[n_vial]-TI_dict[n_vial]) for t in t_meas]
            
            mv_pred=[]
            res_m_=[]
            obj_mi = 0
            ob_mi = 0
            count_m = 0
            for i in range(0,len(t_meas_scaled)):
                mv_inter = interpolation(m,m.mV,n_vial,t_meas_scaled[i])
                mv_pred.append(mv_inter)
                
                if not np.isnan(mv_meas[i]):
                    res_m = mv_pred[i]-mv_meas[i]            
                    count_m += 1
                    res_m_.append(res_m/0.01)
                    obj_mi += (res_m/0.01)**2 # 0.01g error             
                    ob_mi += res_m**2
            if n_vial >= data_stru['data_config']['n_v0']: 
                Count_m += count_m
                res_m_assemble.extend(res_m_)
                obj_m += obj_mi#/count_m/collect_vial
                ob_m += ob_mi
    
            # permeate residual squared / uncertainty / # of measurements
            if type(data_stru['data_raw'][n_vial-1]['cV_avg']) != list:
                if n_vial > data_stru['data_config']['n_extra']:
                    res_cp = m.cV[n_vial,m.tau.last()]-data_stru['data_raw'][n_vial-1]['cV_avg']
                    res_cp_assemble.append(res_cp/(0.03*data_stru['data_raw'][n_vial-1]['cV_avg']))
                    obj_cp += (res_cp/(0.03*data_stru['data_raw'][n_vial-1]['cV_avg']))**2 # 3% error
                    ob_cp += (res_cp/(data_stru['data_raw'][n_vial-1]['cV_avg']))**2
                    Count_cp = collect_vial
            else:
                cp_pred=[]
                obj_cpi = 0
                ob_cpi = 0
                count_cp = 0
                if not all(np.isnan(cp_meas)):
                    for i in range(0,len(t_meas_scaled)):
                        cp_inter = interpolation(m,m.cV,n_vial,t_meas_scaled[i])
                        cp_pred.append(cp_inter)
                        if not np.isnan(cp_meas[i]):   
                            res_cp = cp_pred[i]-cp_meas[i]
                            count_cp += 1
                            res_cp_assemble.append(res_cp/(0.03*cp_meas[i]))
                            obj_cpi += (res_cp/(0.03*cp_meas[i]))**2 # 3% error
                            ob_cpi += (res_cp/(cp_meas[i]))**2
                    Count_cp = collect_vial
                    obj_cp += obj_cpi#/count_cp/collect_vial
                    ob_cp += ob_cpi
         
            # retentate residual squared / uncertainty / # of measurements
            cf_pred=[]
            obj_cfi = 0
            ob_cfi = 0
            count_cf = 0
            if not all(np.isnan(cf_meas)):
                for i in range(0,len(t_meas_scaled)):
                    cf_inter = interpolation(m,m.cF,n_vial,t_meas_scaled[i])
                    cf_pred.append(cf_inter)
                    if not np.isnan(cf_meas[i]):   
                        res_cf = cf_pred[i]-cf_meas[i]
                        count_cf += 1
                        res_cf_assemble.append(res_cf/(0.003*cf_meas[i]))
                        obj_cfi += (res_cf/(0.003*cf_meas[i]))**2 # 0.3% error
                        ob_cfi += (res_cf/(cf_meas[i]))**2
                if n_vial >= data_stru['data_config']['n_v0']:
                    Count_cf += count_cf
                    obj_cf += obj_cfi#/count_cf/collect_vial
                    ob_cf += ob_cfi
                else: 
                    count_cf0 += count_cf
                    obj_cf0 += obj_cfi#/count_cf/collect_vial
                    ob_cf0 += ob_cfi
        
        m.count_m = Count_m
        m.count_cv = Count_cp
        m.count_cr = Count_cf
        m.count_cr0 = count_cf0
        m.count = Count_m+Count_cp+Count_cf+count_cf0
        
        m.res_m = res_m_assemble
        m.res_cp = res_cp_assemble
        m.res_cf = res_cf_assemble
        
        m.obj_m = obj_m/Count_m
        m.obj_cv = obj_cp/Count_cp
        m.obj_cr = (obj_cf0+obj_cf)/(count_cf0+Count_cf)
        m.obj_cr_tru = obj_cf/Count_cf
        m.obj_tru = 1e4*(obj_m/Count_m+obj_cp/Count_cp+obj_cf/Count_cf)
        m.llh1 = m.count_m*log(m.obj_m) + m.count_cv*log(m.obj_cv) + (m.count_cr0+m.count_cr)*log(m.obj_cr)#m.count * log((obj_m+obj_cp+obj_cf0+obj_cf)/m.count)
        m.llh2 = Count_m*log(ob_m/Count_m) + Count_cp*log(ob_cp/Count_cp) + (count_cf0+Count_cf)*log((ob_cf0+ob_cf)/(count_cf0+Count_cf))
        
        return 1e4*(obj_m/Count_m + obj_cp/Count_cp + (obj_cf0+obj_cf)/(count_cf0+Count_cf))
    
    # pyomo model instance
    model_builder = model_builder or model_construct_inter
    instance = model_builder(data_stru, mode, theta, sim_opt, B_form)
    #instance.pprint()
    if sim_opt:
        instance.Obj_1 = Objective(expr = 1)
        instance.Obj = Expression(rule=obj_rule)
    else:
        instance.Obj = Objective(rule=obj_rule, sense=minimize)

    #Try initialize
    try:
        #Simulate the model using scipy
        sim = Simulator(instance, package='casadi') 
        tsim, profiles = sim.simulate(numpoints=300, integrator='idas')
        #Discretize the model using finite_difference
        TransformationFactory('dae.finite_difference').apply_to(instance, nfe=300, scheme='BACKWARD')
        #Initialize the discretized model using the simulator profiles
        sim.initialize_model()
    except Exception as e:
        print(f"Initialization failed: {e}. Applying discretization without initialization.")
        TransformationFactory('dae.finite_difference').apply_to(instance, nfe=300, scheme='BACKWARD')

    solver = SolverFactory('ipopt')
    solver.options["linear_solver"] = "ma97"
    solver.options["max_iter"] = 3000
    #solver.options["halt_on_ampl_error"] = "yes" # option 1
    #solver.options["print_level"] = 1 # option 2
    try:
        results = solver.solve(instance,tee=True)
        if results.solver.termination_condition == TerminationCondition.optimal:
            results.write()
        else:
            print(f"Solver termination: {results.solver.termination_condition}, resolving...")
            results = solver.solve(instance,tee=True)
            assert results.solver.termination_condition == TerminationCondition.optimal, (
                    f"Solver failed: Non-optimal termination condition "
                    f"{results.solver.termination_condition}. ")
    except:
        print("Retrying with adjusted solver settings...")
        solver.options["linear_solver"] = "ma57"
        solver.options['max_iter']=5000       
        results = solver.solve(instance,tee=True)
        assert results.solver.termination_condition == TerminationCondition.optimal, (
                f"Solver failed again: Non-optimal termination condition "
                f"{results.solver.termination_condition}. "
                "Try alternative initialization or further debugging.")

    # Process results if optimal
    if results.solver.termination_condition == TerminationCondition.optimal:
        if mode !='DATA':
            fit_stru, sim_stru=save_model(instance, LO=True, B_form=B_form, LOUD=True)
        else:
            fit_stru, sim_stru=save_model(instance, LO=False, B_form=B_form, LOUD=True)
        sim_inter = inter_model(data_stru, sim_stru, fit_stru)
    else:
        print("Non-optimal solution after retries. Check initialization or parameter settings.")
        fit_stru, sim_stru, sim_inter = None, None, None
        return fit_stru, sim_stru, sim_inter

    if LOUD:
        time = []
        cIn = []
        cH = []
        Jw = []
        Js = []       
        for i in range(data_stru['data_config']['n']):   
            time.extend(sim_stru[i]['time'])
            cIn.extend(sim_stru[i]['cIn'])
            cH.extend(sim_stru[i]['cH'])
            Jw.extend(sim_stru[i]['Jw'])
            Js.extend(sim_stru[i]['Js'])
        
        sim_data = {'time': time,
                'cIn': cIn,
                'cH': cH,
                'Jw': Jw,
                'Js': Js}
        print(sim_data)
        fname = 'sim_data-dat'+str(data_stru['dataset'])       
        #create data frame from dictionary
        sim_datapd = pd.DataFrame(sim_data)
        
        #save dataframe to csv file
        sim_datapd.to_csv(fname+".csv", index=False)
        
        #validate the csv file by importing it
        #print(pd.read_csv(fname+".csv"))

    # Finish print statement
    print("###################################################################")

    return fit_stru, sim_stru, sim_inter


def nested_dict_values(dict_nest):
    ''' 
    This function accepts a nested dictionary as argument
        and iterates over all values of nested dictionaries
    '''
    dict_value=[]
    # Iterate over all values of given dictionary
    for value in dict_nest.values():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over all its values
            for v in  nested_dict_values(value):
                dict_value.append(v)
        else:
            # If value is not dict type then append the value
            dict_value.append(value)
    return dict_value

def nested_dict_keys(dict_nest):
    ''' 
    This function accepts a nested dictionary as argument
        and iterates over all keys of nested dictionaries
    '''   
    dict_key=[]
    # Iterate over all values of given dictionary
    for key in dict_nest.keys():
        # Check if value is of dict type
        if isinstance(dict_nest[key], dict):           
            # If value is dict then iterate over all its keys
            for v in nested_dict_keys(dict_nest[key]):
                dict_key.append(key+'-'+str(v))
        else:
            # If value is not dict type then yield the value
            dict_key.append(key)
    return dict_key

def nested_dict_update(dict_nest,new_value,i=0):
    ''' 
    This function accepts a nested dictionary and a new list of values as arguments
        and updates all values of nested dictionaries
    '''    
    # Iterate over all keys of given dictionary
    for key in dict_nest.keys():
        # Check if value is of dict type
        if isinstance(dict_nest[key], dict):
            # If value is dict then iterate over all its keys
            for v in nested_dict_update(dict_nest[key],new_value,i=i):
                v = new_value[i]
                i += 1
        else:
            # If value is not dict type then update the value
            dict_nest[key] = new_value[i]
            i += 1
    return dict_nest


def calc_FIM(data_stru, mode, theta=None, step=1e-8, formula='backward', B_form='single', model_builder=None):
    """
    Calculate FIM
    
    Arguments:
        data_stru: dict, experimental data dictionary
        mode: str, experiment mode, {DATA, lag, overflow}
        theta: dict, preset parameter values
        step: float, relative step change for parameters
        formula: str, finite difference scheme option, {backward, forward, central}
        B_form: float/str, different solute permeability coefficient formula
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration 
    
    Returns:
        doe_stru: dict, FIM results
    """
    sim_opt = True
    if theta is None:
        sim_opt = False
        fit_stru_p, sim_stru_p, sim_inter_p = solve_model(data_stru, mode, theta, sim_opt, B_form, model_builder=model_builder)
        theta = fit_stru_p['parameters']
        sim_opt = True
    theta_p = copy.deepcopy(theta)
    theta_p_v=nested_dict_values(theta_p)
    theta_p1_v = []
    theta_p2_v = []

    # Step changes
    for i in theta_p_v: 
        
        if formula == 'central':
            if i!=0:
                theta_p1_v.append(i * (1-step))
                theta_p2_v.append(i * (1+step))
            else:
                theta_p1_v.append(-step)
                theta_p2_v.append(step)
        elif formula == 'backward':
            if i!=0:
                theta_p1_v.append(i * (1-step))
            else:
                theta_p1_v.append(-step)
            theta_p2_v.append(i)
        else: #forward
            theta_p1_v.append(i)
            if i!=0:
                theta_p2_v.append(i * (1+step))
            else:
                theta_p2_v.append(step)
    
    # Prepare prediction covariance matrix        
    if sim_opt == True:
        fit_stru_p, sim_stru_p, sim_inter_p = solve_model(data_stru, mode, theta_p, sim_opt, B_form, model_builder=model_builder)
    var_pred=[]
    for n_vial in range(data_stru['data_config']['n']):
        var_pred = np.append(var_pred,0.01 ** 2 * np.ones(len(sim_inter_p[n_vial]['mV'])))
        var_pred = np.append(var_pred,(0.03 * sim_inter_p[n_vial]['cV'])**2)
        var_pred = np.append(var_pred,(0.003 * sim_inter_p[n_vial]['cF'])**2)
    cov_pred = np.diag(var_pred)
    
    # Prepare Jacobian
    Jac = []
    for i in range(len(theta_p_v)):
        theta_pk1_v = theta_p_v.copy()
        theta_pk2_v = theta_p_v.copy()
        theta_pk1_v[i] = theta_p1_v[i]    
        theta_pk2_v[i] = theta_p2_v[i]       
        theta_pk1 = copy.deepcopy(theta_p)
        theta_pk2 = copy.deepcopy(theta_p)
        theta_pk1 = nested_dict_update(theta_pk1,theta_pk1_v)
        theta_pk2 = nested_dict_update(theta_pk2,theta_pk2_v)
        print(theta_pk1)
        print(theta_pk2)
        
        sim_opt = True
        fit_stru_pk1, sim_stru_pk1, sim_inter_pk1 = solve_model(data_stru, mode, theta_pk1, sim_opt, B_form, model_builder=model_builder)
        fit_stru_pk2, sim_stru_pk2, sim_inter_pk2 = solve_model(data_stru, mode, theta_pk2, sim_opt, B_form, model_builder=model_builder)
        jac=[]
        for n_vial in range(data_stru['data_config']['n']):
            jac = np.append(jac,sim_inter_pk2[n_vial]['mV']-sim_inter_pk1[n_vial]['mV'])
            jac = np.append(jac,sim_inter_pk2[n_vial]['cV']-sim_inter_pk1[n_vial]['cV'])
            jac = np.append(jac,sim_inter_pk2[n_vial]['cF']-sim_inter_pk1[n_vial]['cF'])
        # Predictions in row
        if theta_p_v[i]!=0:
            delta = step*theta_p_v[i]
        else:
            delta = step
        if len(Jac) == 0:
            Jac = jac/delta
        else:
            Jac = np.vstack([Jac, jac/delta])

    doe_stru = dict()
    if formula == 'central':
        doe_stru['Jac'] = Jac/2
    else:
        doe_stru['Jac'] = Jac
    FIM = doe_stru['Jac'] @ np.linalg.inv(cov_pred) @ doe_stru['Jac'].T
    doe_stru['Jac'] = doe_stru['Jac'].tolist()
    doe_stru['FIM'] = FIM.tolist()

    # Compute eigenvalues of FIM
    w, v = np.linalg.eigh(FIM)
    doe_stru['eig_val'] = w.tolist()
    doe_stru['eig_vec'] = v.tolist() # in columns
    theta_p_name=nested_dict_keys(theta_p)
    maxind = np.argmax(abs(v), axis=0)
    doe_stru['eig_dir']=[theta_p_name[i] for i in (maxind)]
    doe_stru['trace'] = np.trace(FIM).tolist()
    doe_stru['det'] = np.linalg.det(FIM).tolist()
    doe_stru['min_eig'] = min(w)
    doe_stru['cond'] = max(w) / min(w)
    try:
        doe_stru['V'] = np.linalg.inv(FIM).tolist()
        doe_stru['std'] = np.sqrt(np.diag(doe_stru['V'])).tolist()
    except:
        print('No inverse of FIM')

    print(doe_stru)
    return doe_stru

def correlation_from_covariance(covariance):
    ''' 
    This function calculate correlation matrix from covariance matrix.
    '''    
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    correlation[covariance == 0] = 0
    print(correlation)
    return correlation

def store_json(file_name,structure):
    with open(file_name, "w") as json_file:
        json.dump(structure, json_file)
        
        
def solve_model_B_fix(data_stru, mode, theta=None, sim_opt=False, B_form=1, sigma_fixed=True, LOUD=False, model_builder=None):
    """
    Solve pyomo model
    
    Arguments:
        data_stru: dict, experimental data dictionary
        mode: str, experiment mode, {DATA, lag, overflow}
        theta: dict, preset parameter values
        sim_opt: boolean, if run simulation with fixed parameter values
        B_form: float/str, different solute permeability coefficient formula
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration
        sigma_fixed: boolean, if fix sigma value
        LOUD: boolean, if print out parameter results and store model predictions
    
    Returns:
        fit_stru: dict, parameter fit results
        sim_stru: dict, model predictions
        sim_inter: dict, model predictions for experimental measurements 
    """

    print("\n\n###################################################################")
    print("Creating and solving the Pyomo model with the following settings: ")
    print("mode =", mode)
    print("theta =",theta)
    print("sim_opt =", sim_opt)
    print("B_form =", B_form)
    print("sigma_fixed =",sigma_fixed)
    print(" ")
    
    # interpolation function
    def interpolation(m,var,n_vial,t):
        """
        Interpolation function for pyomo model
        """
        tp = list(m.tau)
        for j in range(0,len(tp)-1):
            if tp[j+1]>=t and tp[j]<=t:
                var_inter = (var[n_vial,tp[j+1]]-var[n_vial,tp[j]]) * (t - tp[j])/(tp[j+1]-tp[j]) + var[n_vial,tp[j]]
        return var_inter
    
    # objective function
    def obj_rule(m):
        """
        pyomo Objective function (residual calculated using interpolation from model)
        """
        obj_m = 0
        obj_cp = 0
        obj_cf0 = 0
        obj_cf = 0
        ob_m = 0
        ob_cp = 0
        ob_cf0 = 0
        ob_cf = 0
        res_m_assemble = []
        res_cp_assemble = []
        res_cf_assemble = []
        
        collect_vial = data_stru['data_config']['n'] - data_stru['data_config']['n_extra']
        Count_m = 0
        Count_cp = 0
        Count_cf = 0
        count_cf0 = 0
        t_delay = data_stru['data_raw'][0]['time'][0]
        TF_list = [data_stru['data_raw'][i]['time'][-1]-t_delay for i in range(data_stru['data_config']['n'])]
        TF_dict = dict(zip(m.n_vial,TF_list)) # unscaled time elapse for each vial 
        TI_list = [data_stru['data_raw'][i]['time'][0]-t_delay for i in range(data_stru['data_config']['n'])]
        TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time for each vial
    
        for n_vial in m.n_vial:
            t_meas = data_stru['data_raw'][n_vial-1]['time'] - t_delay
            mv_meas = data_stru['data_raw'][n_vial-1]['mass']
            cp_meas = data_stru['data_raw'][n_vial-1]['cV_avg']
            cf_meas = data_stru['data_raw'][n_vial-1]['cF_exp']
    
            t_meas_scaled = [(t-TI_dict[n_vial])/(TF_dict[n_vial]-TI_dict[n_vial]) for t in t_meas]
            
            mv_pred=[]
            res_m_=[]
            obj_mi = 0
            ob_mi = 0
            count_m = 0
            for i in range(0,len(t_meas_scaled)):
                mv_inter = interpolation(m,m.mV,n_vial,t_meas_scaled[i])
                mv_pred.append(mv_inter)
                
                if not np.isnan(mv_meas[i]):
                    res_m = mv_pred[i]-mv_meas[i]            
                    count_m += 1
                    res_m_.append(res_m/0.01)
                    obj_mi += (res_m/0.01)**2 # 0.01g error             
                    ob_mi += res_m**2
            if n_vial >= data_stru['data_config']['n_v0']: 
                Count_m += count_m
                res_m_assemble.extend(res_m_)
                obj_m += obj_mi#/count_m/collect_vial
                ob_m += ob_mi
    
            # permeate residual squared / uncertainty / # of measurements
            if type(data_stru['data_raw'][n_vial-1]['cV_avg']) != list:
                if n_vial > data_stru['data_config']['n_extra']:
                    res_cp = m.cV[n_vial,m.tau.last()]-data_stru['data_raw'][n_vial-1]['cV_avg']
                    res_cp_assemble.append(res_cp/(0.03*data_stru['data_raw'][n_vial-1]['cV_avg']))
                    obj_cp += (res_cp/(0.03*data_stru['data_raw'][n_vial-1]['cV_avg']))**2 # 3% error
                    ob_cp += (res_cp/(data_stru['data_raw'][n_vial-1]['cV_avg']))**2
                    Count_cp = collect_vial
            else:
                cp_pred=[]
                obj_cpi = 0
                ob_cpi = 0
                count_cp = 0
                if not all(np.isnan(cp_meas)):
                    for i in range(0,len(t_meas_scaled)):
                        cp_inter = interpolation(m,m.cV,n_vial,t_meas_scaled[i])
                        cp_pred.append(cp_inter)
                        if not np.isnan(cp_meas[i]):   
                            res_cp = cp_pred[i]-cp_meas[i]
                            count_cp += 1
                            res_cp_assemble.append(res_cp/(0.03*cp_meas[i]))
                            obj_cpi += (res_cp/(0.03*cp_meas[i]))**2 # 3% error
                            ob_cpi += (res_cp/(cp_meas[i]))**2
                    Count_cp = collect_vial
                    obj_cp += obj_cpi#/count_cp/collect_vial
                    ob_cp += ob_cpi
         
            # retentate residual squared / uncertainty / # of measurements
            cf_pred=[]
            obj_cfi = 0
            ob_cfi = 0
            count_cf = 0
            if not all(np.isnan(cf_meas)):
                for i in range(0,len(t_meas_scaled)):
                    cf_inter = interpolation(m,m.cF,n_vial,t_meas_scaled[i])
                    cf_pred.append(cf_inter)
                    if not np.isnan(cf_meas[i]):   
                        res_cf = cf_pred[i]-cf_meas[i]
                        count_cf += 1
                        res_cf_assemble.append(res_cf/(0.003*cf_meas[i]))
                        obj_cfi += (res_cf/(0.003*cf_meas[i]))**2 # 0.3% error
                        ob_cfi += (res_cf/(cf_meas[i]))**2
                if n_vial >= data_stru['data_config']['n_v0']:
                    Count_cf += count_cf
                    obj_cf += obj_cfi#/count_cf/collect_vial
                    ob_cf += ob_cfi
                else: 
                    count_cf0 += count_cf
                    obj_cf0 += obj_cfi#/count_cf/collect_vial
                    ob_cf0 += ob_cfi
        
        m.count_m = Count_m
        m.count_cv = Count_cp
        m.count_cr = Count_cf
        m.count_cr0 = count_cf0
        m.count = Count_m+Count_cp+Count_cf+count_cf0
        
        m.res_m = res_m_assemble
        m.res_cp = res_cp_assemble
        m.res_cf = res_cf_assemble
        
        m.obj_m = obj_m/Count_m
        m.obj_cv = obj_cp/Count_cp
        m.obj_cr = (obj_cf0+obj_cf)/(count_cf0+Count_cf)
        m.obj_cr_tru = obj_cf/Count_cf
        m.obj_tru = 1e4*(obj_m/Count_m+obj_cp/Count_cp+obj_cf/Count_cf)
        m.llh1 = m.count_m*log(m.obj_m) + m.count_cv*log(m.obj_cv) + (m.count_cr0+m.count_cr)*log(m.obj_cr)#m.count * log((obj_m+obj_cp+obj_cf0+obj_cf)/m.count)
        m.llh2 = Count_m*log(ob_m/Count_m) + Count_cp*log(ob_cp/Count_cp) + (count_cf0+Count_cf)*log((ob_cf0+ob_cf)/(count_cf0+Count_cf))
        
        return 1e4*(obj_m/Count_m + obj_cp/Count_cp + (obj_cf0+obj_cf)/(count_cf0+Count_cf))
    
    # pyomo model instance
    model_builder = model_builder or model_construct_inter
    instance = model_builder(data_stru, mode, theta, sim_opt, B_form)
    if B_form=='single':
        instance.B.fixed=True
    else:
        instance.beta_0.fixed=True
        instance.beta_1.fixed=True
        if B_form > 1:
            instance.beta_2.fixed=True
    if sigma_fixed:
        instance.sigma.fixed=True
    if sim_opt:
        instance.Obj_1 = Objective(expr = 1)
        instance.Obj = Expression(rule=obj_rule)
    else:
        instance.Obj = Objective(rule=obj_rule, sense=minimize)

    #Try initialize
    try:
        #Simulate the model using scipy
        sim = Simulator(instance, package='casadi') 
        tsim, profiles = sim.simulate(numpoints=300, integrator='idas')
        #Discretize the model using finite difference
        TransformationFactory('dae.finite_difference').apply_to(instance, nfe=300, scheme='BACKWARD')
        #Initialize the discretized model using the simulator profiles
        sim.initialize_model()
    except:
        TransformationFactory('dae.finite_difference').apply_to(instance, nfe=300, scheme='BACKWARD')

    solver = SolverFactory('ipopt')
    solver.options["linear_solver"] = "ma97"

    results = solver.solve(instance,tee=True)
    assert results.solver.termination_condition == TerminationCondition.optimal, "Optimization fail"
    results.write()
    #solver.options["halt_on_ampl_error"] = "yes" # option 1
    #solver.options["print_level"] = 1 # option 2
    #solver.solve(instance,tee=True).write()
    #instance.display()

    if mode !='DATA':
        fit_stru, sim_stru=save_model(instance, LO=True, B_form=B_form, LOUD=True)
    else:
        fit_stru, sim_stru=save_model(instance, LO=False, B_form=B_form, LOUD=True)
    sim_inter = inter_model(data_stru, sim_stru, fit_stru)

    if LOUD:
        time = []
        cIn = []
        cH = []
        Jw = []
        Js = []       
        for i in range(data_stru['data_config']['n']):   
            time.extend(sim_stru[i]['time'])
            cIn.extend(sim_stru[i]['cIn'])
            cH.extend(sim_stru[i]['cH'])
            Jw.extend(sim_stru[i]['Jw'])
            Js.extend(sim_stru[i]['Js'])
        
        sim_data = {'time': time,
                'cIn': cIn,
                'cH': cH,
                'Jw': Jw,
                'Js': Js}
        print(sim_data)
        fname = 'sim_data-dat'+str(data_stru['dataset'])       
        #create data frame from dictionary
        sim_datapd = pd.DataFrame(sim_data)
        
        #save dataframe to csv file
        sim_datapd.to_csv(fname+".csv", index=False)
        
        #validate the csv file by importing it
        #print(pd.read_csv(fname+".csv"))

    print("###################################################################")

    return fit_stru, sim_stru, sim_inter
