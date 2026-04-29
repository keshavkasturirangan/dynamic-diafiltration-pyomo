"""
Library of functions for diafiltration experiment modeling
Xinhong Liu
University of Notre Dame
"""

import numpy as np
import pandas as pd
import scipy.io as spio
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.patheffects as patheffects
from matplotlib.lines import Line2D
from pathlib import Path
import idaes
import time
import copy
import os
import json
import math
import ctypes.util
from sklearn.metrics import r2_score
from dataclasses import dataclass, field

try:
    import seaborn as sns
except ImportError:  # pragma: no cover - optional plotting dependency
    sns = None

from pyomo.environ import *
from pyomo.dae import *
import idaes.core.util.scaling as iscale
from pyomo.core.base.componentuid import ComponentUID
from pyomo.contrib.parmest.experiment import Experiment as ParmestExperiment
from pyomo.contrib.parmest.parmest import Estimator, SSE, SSE_weighted, compute_covariance_matrix
from pyomo.contrib.doe.doe import DesignOfExperiments


FIGURES_DIR = os.path.join("UnifiedFramework", "DATA3", "figures")


def _has_casadi():
    """Return True when CasADi is importable in the current environment."""
    try:
        import casadi  # noqa: F401
    except ImportError:
        return False
    return True


def _preferred_ipopt_linear_solver():
    """Prefer HSL when present, otherwise fall back to the portable MUMPS build."""
    return "ma97" if ctypes.util.find_library("hsl") else "mumps"


def _build_ipopt_solver():
    """Create an Ipopt solver configured for the current machine."""
    solver = SolverFactory("ipopt")
    solver.options["linear_solver"] = _preferred_ipopt_linear_solver()
    return solver


@dataclass(frozen=True)
class ExperimentalRun:
    """
    One experimental run family with multiple file-backed variants.

    The same physical run can appear in different notebook/paper forms, such as
    the base fit and the concentration-polarization variant. Each variant is a
    separate instance of the same run family.
    """

    run_id: str
    root: Path
    data_file: str
    variant: str = "base"
    subdir: str = ""
    fit_file: str = "fit_stru.mat"
    contour_files: tuple[str, ...] = ("contourdata-x_B-y_Lp.csv", "contourdata-x_sigma-y_Lp.csv")
    extra_files: tuple[str, ...] = field(default_factory=tuple)

    @property
    def data_path(self) -> Path:
        return self.root / self.data_file

    @property
    def variant_root(self) -> Path:
        return self.root / self.subdir if self.subdir else self.root

    @property
    def fit_path(self) -> Path:
        return self.variant_root / self.fit_file

    def contour_path(self, name: str) -> Path:
        return self.variant_root / name

    def load_data(self):
        """Load and normalize the main MAT file for this run instance."""
        return _normalize_conductivity_measurements(loadmat(str(self.data_path)).get("data_stru"))

    def load_fit(self):
        """Load the saved fit bundle for this run instance."""
        return loadmat(str(self.fit_path)).get("fit_stru")

    def load_contour(self, name: str):
        """Load a contour CSV for this run instance."""
        return pd.read_csv(self.contour_path(name))


def _figure_output_base(name):
    # Make sure the shared figures folder exists before saving anything there.
    os.makedirs(FIGURES_DIR, exist_ok=True)
    return os.path.join(FIGURES_DIR, name)


def _data1_time_origin(data_stru):
    """Return the DATA1 notebook time origin used for all shifted plots."""
    return float(np.asarray(data_stru["data_raw"][0]["time"], dtype=float).reshape(-1)[0])


def _data1_shifted_time(data_stru, vial_index):
    """Return DATA1 vial time shifted by the shared notebook time origin."""
    return np.asarray(data_stru["data_raw"][vial_index]["time"], dtype=float) - _data1_time_origin(data_stru)


def _data1_shifted_time_minutes(data_stru, vial_index):
    """Return DATA1 vial time shifted and converted to minutes."""
    return _data1_shifted_time(data_stru, vial_index) / 60.0


def _load_conductivity_paper():
    """Load the local copy of `conductivity_paper.py` next to this library."""
    import importlib.util
    from importlib.machinery import SourceFileLoader

    paper_path = Path(__file__).resolve().parent / "conductivity_paper.py"
    spec = importlib.util.spec_from_loader(
        "conductivity_paper",
        SourceFileLoader("conductivity_paper", str(paper_path)),
    )
    cp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cp)
    return cp


CONDUCTIVITY_SALT_PARAMS_25C = {
    # Simple built-in defaults for the paper data files.
    "KCl": {
        "epsilon": 78.4,
        "eta": 0.0089,
        "a": 1e-8,
        "z_1": 1,
        "z_2": -1,
        "lambda_0_cation": 73.50,
        "lambda_0_anion": 76.35,
        "lambda_0": 149.85,
    },
    "NaCl": {
        "epsilon": 78.4,
        "eta": 0.0089,
        "a": 1e-8,
        "z_1": 1,
        "z_2": -1,
        "lambda_0_cation": 50.11,
        "lambda_0_anion": 76.35,
        "lambda_0": 126.46,
    },
}


def loadmat(filename):# for fun!
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

    def _coerce_scalar(value):
        """Convert MATLAB numeric scalars to plain Python numbers.

        MATLAB MAT files frequently encode small integers as numpy scalar types
        such as uint8. Those types behave badly later in the workflow when we
        subtract times or build negative offsets, so we normalize them here.
        """
        if isinstance(value, np.generic):
            if np.issubdtype(type(value), np.number):
                return float(value)
            return value.item()
        return value

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
                d[strg] = _coerce_scalar(elem)
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
            elif isinstance(sub_elem, (int, float, np.integer, np.floating, np.bool_)):
                elem_list.append(_coerce_scalar(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    
    return _check_keys(data)


def plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=False,cond=True,preface=False):
    '''
    Plot simulation results comparing with measurements
    '''

    t_delay = _data1_time_origin(data_stru)
    figures = []

    # plot mass data/prediction comparison
    fig = plt.figure(figsize=(4,4))
    for i in range(data_stru['data_config']['n']):
        time = np.asarray(data_stru['data_raw'][i]['time'], dtype=float)
        mass = np.asarray(data_stru['data_raw'][i]['mass'], dtype=float)
        plt.plot((time-t_delay)/60, mass, 'r.', markersize=4)
        if plot_pred:
            plt.plot((np.asarray(sim_stru[i]['time'], dtype=float)-t_delay)/60,
                     np.asarray(sim_stru[i]['mV'], dtype=float),
                     'b', linewidth=3, alpha=.6)

    plt.plot([],[],'r.',markersize=4,label='Measurements')
    if plot_pred:
        plt.plot([],[],'b',linewidth=3,alpha=.6,label='Predictions')

    if preface:
        plt.xlabel('Time',fontsize=24,fontweight='bold')
        plt.ylabel('Mass',fontsize=24,fontweight='bold')
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    else:
        plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
        plt.ylabel('Mass [g]',fontsize=16,fontweight='bold')

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
        plt.legend(fontsize=10,loc='best')
    if LOUD:
        fname = _figure_output_base('mass-dat'+str(data_stru['dataset']))
        fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
    figures.append(fig)

    # plot concentration data/prediction comparison
    fig = plt.figure(figsize=(4,4))
    first_time = _data1_time_origin(data_stru)
    try:
        is_data2_workflow = float(data_stru.get("dataset", 0)) > 1_000
    except Exception:
        is_data2_workflow = False
    plt.plot((first_time-t_delay)/60,
             np.asarray(sim_stru[0]['cF'], dtype=float)[0],
             'ms',markersize=8,clip_on=False)

    plt.plot([],[],'ms',markersize=8,clip_on=False,label='Retentate (Measurement)')
    if plot_pred:
        plt.plot([],[],'g',linewidth=3,label='Retentate (Prediction)')
    plt.plot([],[],'cs',markersize=8,label='Vial (Measurement)')
    if plot_pred:
        plt.plot([],[],'r^',markersize=8,label='Vial (Prediction)')
        plt.plot([],[],'r-',linewidth=3,alpha=.6,label='Permeate (Prediction)')

    for i in range(data_stru['data_config']['n']):
        time = np.asarray(data_stru['data_raw'][i]['time'], dtype=float)
        c_v_avg = np.asarray(data_stru['data_raw'][i]['cV_avg'], dtype=float)
        if (not is_data2_workflow) or c_v_avg.ndim == 0 or c_v_avg.size == 1:
            plt.plot((time[-1]-t_delay)/60,
                     float(c_v_avg.reshape(-1)[0]),
                     'cs',markersize=8)
        else:
            n = min(len(time), len(c_v_avg))
            plt.plot((time[:n]-t_delay)/60,
                     c_v_avg[:n],
                     'cs',markersize=8,clip_on=False)
        if cond:
            cF_exp = data_stru['data_raw'][i]['cF_exp']
            if is_data2_workflow:
                if isinstance(cF_exp, float):
                    plt.plot((time[-1]-t_delay)/60, cF_exp,'ms',markersize=8)
                elif len(cF_exp)>1:
                    plt.plot((time-t_delay)/60,
                             np.asarray(cF_exp, dtype=float),'ms',markersize=8)
            else:
                if isinstance(cF_exp, float):
                    plt.plot((time[-1]-t_delay)/60, cF_exp,'ms',markersize=8)
                elif len(cF_exp)>1:
                    plt.plot((time-t_delay)/60,
                             np.asarray(cF_exp, dtype=float),'ms',markersize=8)
        if plot_pred:
            sim_time = np.asarray(sim_stru[i]['time'], dtype=float)
            plt.plot((sim_time-t_delay)/60,
                     np.asarray(sim_stru[i]['cF'], dtype=float),
                     'g',linewidth=3,alpha=.6)
            plt.plot((sim_time-t_delay)/60,
                     np.asarray(sim_stru[i]['cH'], dtype=float),
                     'r-',linewidth=3,alpha=.6)
            if not preface:
                sim_c_v = np.asarray(sim_stru[i]['cV'], dtype=float)
                if (not is_data2_workflow) or c_v_avg.ndim == 0 or c_v_avg.size == 1:
                    plt.plot((sim_time[-1]-t_delay)/60,
                             sim_c_v[-1],
                             'r^',markersize=8,alpha=.6)
                else:
                    valid = ~np.isnan(c_v_avg[: min(len(time), len(c_v_avg))])
                    if np.any(valid):
                        time_shifted = time[: min(len(time), len(c_v_avg))] - t_delay
                        interp_cv = interpolate.interp1d(sim_time, sim_c_v, fill_value='extrapolate')
                        plt.plot(time_shifted[valid]/60,
                                 interp_cv(time_shifted[valid]),
                                 'r^',markersize=8,alpha=.6)

    if not cond:
        plt.plot((np.asarray(data_stru['data_raw'][-1]['time'], dtype=float)[-1]-t_delay)/60,
                 np.asarray(data_stru['data_raw'][-1]['cF_exp'], dtype=float),
                 'ms',markersize=8)

    if preface:
        plt.xlabel('Time',fontsize=24,fontweight='bold')
        plt.ylabel('Concentration',fontsize=24,fontweight='bold')
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    else:
        plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
        plt.ylabel('Concentration [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in",top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    if xticks:
        xticks[0].label1.set_visible(False)
    if lg:
        plt.legend(fontsize=10,loc='best')
    if LOUD:
        fname = _figure_output_base('concentration-dat'+str(data_stru['dataset']))
        fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
    figures.append(fig)

    # plot mass of stirred cell
    if stirc_mass:
        fig = plt.figure(figsize=(4,4))
        for i in range(data_stru['data_config']['n']):
            if plot_pred:
                plt.plot((np.asarray(sim_stru[i]['time'], dtype=float)-t_delay)/60,
                         np.asarray(sim_stru[i]['mF'], dtype=float),
                         'b',linewidth=3,alpha=.6)

        plt.plot([],[],'r.',markersize=4,label='Mass (Measurements)')
        if plot_pred:
            plt.plot([],[],'b',linewidth=3,alpha=.6,label='Mass (Predictions)')

        plt.plot([],[],'ms',markersize=8,clip_on=False,label='Retentate (Measurement)')
        if plot_pred:
            plt.plot([],[],'g',linewidth=3,label='Retentate (Prediction)')
        plt.plot([],[],'cs',markersize=8,label='Vial (Measurement)')
        if plot_pred:
            plt.plot([],[],'r^',markersize=8,label='Vial (Prediction)')
            plt.plot([],[],'r-',linewidth=3,alpha=.6,label='Permeate (Prediction)')

        plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
        plt.ylabel('Mass in Stirred Cell [g]',fontsize=16,fontweight='bold')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tick_params(direction="in")
        ytop = max([max(np.asarray(sim_stru[i]['mF'], dtype=float)) for i in range(data_stru['data_config']['n'])])
        plt.ylim(bottom=0,top=ytop+0.5)
        if lg:
            plt.legend(fontsize=12.5,loc='best')
        if LOUD:
            fname = _figure_output_base('stirc_mass-dat'+str(data_stru['dataset']))
            fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
    figures.append(fig)

    return figures


def _plot_sim_comparison_data1_legacy(data_stru, sim_stru, plot_pred=True, lg=False, LOUD=False):
    """Mirror the notebook Figure 2 plotting behavior from diafiltration_plots.py."""
    t_delay = float(data_stru['data_raw'][0]['time'][0])
    figures = []

    fig = plt.figure(figsize=(4, 4))
    for i in range(data_stru['data_config']['n']):
        time = np.asarray(data_stru['data_raw'][i]['time'], dtype=float)
        mass = np.asarray(data_stru['data_raw'][i]['mass'], dtype=float)
        plt.plot((time - t_delay) / 60, mass, 'r.', markersize=4)
        if plot_pred:
            sim_time = np.asarray(sim_stru[i]['time'], dtype=float)
            plt.plot((sim_time - t_delay) / 60,
                     np.asarray(sim_stru[i]['mV'], dtype=float),
                     'b', linewidth=3, alpha=.6)

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
    if LOUD:
        fname = _figure_output_base('mass-dat' + str(data_stru['dataset']))
        fig.savefig(fname + '.png', dpi=300, bbox_inches='tight')
    figures.append(fig)

    fig = plt.figure(figsize=(4, 4))
    plt.plot((float(np.asarray(data_stru['data_raw'][0]['time'], dtype=float)[0]) - t_delay) / 60,
             float(np.asarray(sim_stru[0]['cF'], dtype=float)[0]),
             'ms', markersize=8, clip_on=False)
    for i in range(data_stru['data_config']['n']):
        time = np.asarray(data_stru['data_raw'][i]['time'], dtype=float)
        c_v_avg = data_stru['data_raw'][i]['cV_avg']
        if not isinstance(c_v_avg, list):
            plt.plot((float(time[-1]) - t_delay) / 60,
                     float(c_v_avg),
                     'cs', markersize=8)
        else:
            plt.plot((time - t_delay) / 60,
                     np.asarray(c_v_avg, dtype=float),
                     'cs', markersize=8, clip_on=False)

        c_f_exp = data_stru['data_raw'][i]['cF_exp']
        c_f_exp_arr = np.asarray(c_f_exp, dtype=float)
        if c_f_exp_arr.ndim == 0 or c_f_exp_arr.size == 1:
            plt.plot((float(time[-1]) - t_delay) / 60,
                     float(c_f_exp_arr.reshape(-1)[0]),
                     'ms', markersize=8, clip_on=False)
        else:
            plt.plot((time - t_delay) / 60,
                     c_f_exp_arr,
                     'ms', markersize=8, clip_on=False)

        if plot_pred:
            sim_time = np.asarray(sim_stru[i]['time'], dtype=float)
            sim_c_f = np.asarray(sim_stru[i]['cF'], dtype=float)
            sim_c_h = np.asarray(sim_stru[i]['cH'], dtype=float)
            sim_c_v = np.asarray(sim_stru[i]['cV'], dtype=float)
            plt.plot((sim_time - t_delay) / 60, sim_c_f, 'g', linewidth=3, alpha=.6)
            plt.plot((sim_time - t_delay) / 60, sim_c_h, 'r-', linewidth=3, alpha=.6)
            if not isinstance(c_v_avg, list):
                plt.plot((float(sim_time[-1]) - t_delay) / 60,
                         float(sim_c_v[-1]),
                         'r^', markersize=8, alpha=.6)
            else:
                c_v_avg_arr = np.asarray(c_v_avg, dtype=float)
                n_v = min(len(time), len(c_v_avg_arr))
                valid_v = ~np.isnan(c_v_avg_arr[:n_v])
                interp_cv = interpolate.interp1d(sim_time, sim_c_v, fill_value='extrapolate')
                shifted_time = time[:n_v] - t_delay
                plt.plot(shifted_time[valid_v] / 60,
                         interp_cv(shifted_time[valid_v]),
                         'r^', markersize=8, alpha=.6)

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
    if LOUD:
        fname = _figure_output_base('concentration-dat' + str(data_stru['dataset']))
        fig.savefig(fname + '.png', dpi=300, bbox_inches='tight')
    figures.append(fig)
    return figures


def _resolve_data1_figure2_root(data_root):
    """Prefer the notebook-style DATA1 data_library when present for Figure 2."""
    root = _resolve_data_root(data_root)
    if root.name == "data":
        sibling = root.parent / "data_library"
        if sibling.exists():
            return sibling
    return root


def run_data1_figure2_workflow(data_root=None, save_dir=None):
    """Recreate DATA1 Figure 2 using the notebook's legacy plot_sim_comparison path."""
    root = _resolve_data1_figure2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root
    save_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        ("501.1", root / "data_stru-dataset501.1.mat", root / "501.1 concpolar" / "fit_stru.mat"),
        ("511.12", root / "data_stru-dataset511.12.mat", root / "511.12 concpolar" / "fit_stru.mat"),
    ]

    outputs = []
    panel_paths = []
    for dat, data_path, fit_path in cases:
        if not data_path.exists() or not fit_path.exists():
            continue
        data_stru = loadmat(str(data_path)).get("data_stru")
        fit_stru = loadmat(str(fit_path)).get("fit_stru")
        sim_stru = fit_stru.get("sim_stru") if isinstance(fit_stru, dict) else None
        if sim_stru is None:
            continue
        figs = _plot_sim_comparison_data1_legacy(data_stru, sim_stru, plot_pred=True, lg=False, LOUD=False)
        mass_out = save_dir / f"mass-dat{dat}.png"
        conc_out = save_dir / f"concentration-dat{dat}.png"
        figs[0].savefig(mass_out, dpi=300, bbox_inches="tight")
        figs[1].savefig(conc_out, dpi=300, bbox_inches="tight")
        outputs.extend([str(mass_out), str(conc_out)])
        panel_paths.extend([mass_out, conc_out])
        plt.close(figs[0])
        plt.close(figs[1])

    try:
        from PIL import Image, ImageDraw, ImageFont

        ordered = [
            save_dir / "mass-dat501.1.png",
            save_dir / "concentration-dat501.1.png",
            save_dir / "mass-dat511.12.png",
            save_dir / "concentration-dat511.12.png",
        ]
        if all(path.exists() for path in ordered):
            imgs = [Image.open(path).convert("RGB") for path in ordered]
            target_w = max(im.width for im in imgs)
            target_h = max(im.height for im in imgs)
            resized = [im.resize((target_w, target_h), Image.Resampling.LANCZOS) for im in imgs]
            gap = 30
            label_pad = 35
            canvas = Image.new("RGB", (2 * target_w + 3 * gap, 2 * target_h + 3 * gap + 2 * label_pad), "white")
            draw = ImageDraw.Draw(canvas)
            try:
                font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 28)
            except Exception:
                font = None
            positions = [
                (gap, gap + label_pad),
                (2 * gap + target_w, gap + label_pad),
                (gap, 2 * gap + target_h + 2 * label_pad),
                (2 * gap + target_w, 2 * gap + target_h + 2 * label_pad),
            ]
            for label, img, pos in zip(["A", "B", "C", "D"], resized, positions):
                draw.text((pos[0], pos[1] - label_pad), label, fill="black", font=font)
                canvas.paste(img, pos)
            composite = save_dir / "figure_2.png"
            canvas.save(composite)
            outputs.append(str(composite))
    except Exception:
        pass

    return outputs


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
    fit_stru['Obj_tru'] = value(m.obj_tru)
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
    t_delay = _data1_time_origin(data_stru)
    for i in range(data_stru['data_config']['n']):
        t_meas = _data1_shifted_time(data_stru, i)
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


def model_construct_inter(data_stru, mode, theta=None, sim_opt=False, B_form='single', workflow_family='DATA1'):
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
    t_delay = _data1_time_origin(data_stru)
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
    TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time for each vial
    m.ti = Set(initialize=TI_list)

    workflow_family = str(workflow_family or "DATA1").upper()
    match workflow_family:
        case "DATA2":
            # DATA1 and DATA2 share the same first-principles builder today.
            # Keeping the switch here makes future model-family changes local.
            pass
        case _:
            workflow_family = "DATA1"

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
    
    # Mass balance on the feed/vial mass.
    # Overflow: dmF/dt = -S0 - Am*rho*Jw
    # Lag:      dmF/dt = -S0 - Am*rho*Jw   or   dmF/dt = S
    # DATA:     feed mass is not part of the residual model.
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
     
    # Retentate concentration balance.
    # Overflow/Lag: dcF/dt = (feed/dilution terms + membrane transport terms) / mF
    # DATA:        dcF/dt = Am*rho/M_F0 * (cD*Jw - Js)
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

    # Permeate/vial concentration balance:
    # dcH/dt = Am*rho/mH * (Js - Jw*cH)
    def ode_cH_rule(m, n, t):
        if mode !='DATA' and n <= N_H:
            return m.dcH[n,t] == Am * rho / m.mV[n,t] * (m.Js[n,t] - m.cH[n,t] * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
        else:
            return m.dcH[n,t] == Am * rho / mH * (m.Js[n,t] - m.cH[n,t] * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_cH = Constraint(m.n_vial, m.tau, rule=ode_cH_rule)
    
    # Vial mass balance:
    # dmV/dt = Jw * Am * rho
    def ode_mV_rule(m, n, t):
        return m.dmV[n,t] == m.Jw[n,t] * Am * rho * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_mV = Constraint(m.n_vial, m.tau, rule=ode_mV_rule)
    
    # Vial solute mass balance:
    # d(mV*cV)/dt = Jw * Am * rho * cH
    def ode_cVmV_rule(m, n, t):
        if mode !='DATA' and n <= N_H:
            return m.dcVmV[n,t] == Am * rho * m.Js[n,t] * (TF_dict[n]-TI_dict[n])/Tauf
        else:                
            return m.dcVmV[n,t] == m.Jw[n,t] * m.cH[n,t] * Am * rho * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_cVmV = Constraint(m.n_vial, m.tau, rule=ode_cVmV_rule)
    
    # Global mass constraint used in Lag/Overflow runs:
    # mF(tf) - mF(0) = M_O
    def eqn_S_rule(m):
        return m.mF[m.n_vial.last(),Tauf] - M_F0 == M_O
    if mode !='DATA' and not sim_opt:
        m.eqn_S = Constraint(rule=eqn_S_rule)    
    # Film model for interface concentration:
    # cIn = (cF - cH) * exp(Jw/k) + cH
    def eqn_cIn_rule(m, n, t):
        return m.cIn[n,t] == (m.cF[n,t] - m.cH[n,t]) * exp(m.Jw[n,t]/ k) + m.cH[n,t]
    m.eqn_cIn = Constraint(m.n_vial, m.tau, rule=eqn_cIn_rule)
    
    # Water flux equation:
    # Jw = Lp * (ΔP - (cIn - cH) * ni * sigma * R * T)
    def eqn_Jw_rule(m, n, t):
        return m.Jw[n,t]*36000 == m.Lp *(delP - (m.cIn[n,t] - m.cH[n,t])*ni*m.sigma*R*T)
    m.eqn_Jw = Constraint(m.n_vial, m.tau, rule=eqn_Jw_rule)
    
    # Solute flux equation:
    # Js = B * (cIn - cH)
    # or, for the convection form, a modified flux relation using H and Js_exp.
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
    
    # Algebraic relation for vial concentration:
    # mV * cV = cVmV
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
    
    # B continuity across vials in the per-vial form.
    def B_form_rule1(m,n):
        if n < N_V0-1:
            return m.B[n] == m.B[n+1]
        else:
            return Constraint.Skip

    # B(cF) polynomial form used in the paper variants.
    # Example: B = beta_0 + beta_1*cF^p + beta_2*cF^(p+1) + ...
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

    # B(cIn) polynomial form used by the DATA1/DATA2 model variants.
    # The exponent p is controlled by B_form.
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

    # Initial conditions for the ODE system.
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
            # DATA mode initializes the permeate concentration from the first available measurement.
            # If the saved file has only NaN values, fall back to a tiny positive default.
            if type(data_stru['data_raw'][0]['cV_avg']) == list:
                C_H0 = firstNonNan(data_stru['data_raw'][0]['cV_avg'])
            else:
                C_H0 = data_stru['data_raw'][0]['cV_avg']
            if C_H0 is None or (isinstance(C_H0, float) and np.isnan(C_H0)):
                C_H0 = 1e-6
            yield m.cH[1,0]==C_H0 * 0.8
            yield m.cF[1,0]==firstNonNan(data_stru['data_raw'][0]['cF_exp'])
               
        yield m.cVmV[1,0]==1e-6*1e-6
        yield m.mV[1,0]==1e-6
        yield ConstraintList.End
        
    m.con_boundary = ConstraintList(rule=_init)
        
    return m


def solve_model(
    data_stru,
    mode,
    theta=None,
    sim_opt=False,
    B_form='single',
    LOUD=False,
    workflow_family='DATA1',
    solver_max_iter=3000,
    solver_retry_max_iter=5000,
):
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

    data_stru = _normalize_conductivity_measurements(data_stru)

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
    
    # Objective function:
    # minimize a weighted sum of squared residuals for mass, permeate concentration,
    # and retentate concentration.
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
        
        # Number of vials that contribute to the fitted objective.
        collect_vial = data_stru['data_config']['n'] - data_stru['data_config']['n_extra']
        Count_m = 0
        Count_cp = 0
        Count_cf = 0
        count_cf0 = 0
        t_delay = _data1_time_origin(data_stru)
        TF_list = [data_stru['data_raw'][i]['time'][-1]-t_delay for i in range(data_stru['data_config']['n'])]
        TF_dict = dict(zip(m.n_vial,TF_list)) # unscaled time elapse for each vial 
        TI_list = [data_stru['data_raw'][i]['time'][0]-t_delay for i in range(data_stru['data_config']['n'])]
        TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time for each vial
    
        for n_vial in m.n_vial:
            t_meas = _data1_shifted_time(data_stru, n_vial-1)
            mv_meas = np.asarray(data_stru['data_raw'][n_vial-1]['mass'], dtype=float)
            cp_meas = np.asarray(data_stru['data_raw'][n_vial-1]['cV_avg'], dtype=float)
            cf_meas = np.asarray(data_stru['data_raw'][n_vial-1]['cF_exp'], dtype=float)
    
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
            if np.asarray(data_stru['data_raw'][n_vial-1]['cV_avg']).ndim == 0:
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
                if not np.all(np.isnan(cp_meas)):
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
            if not np.all(np.isnan(cf_meas)):
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
        # Overall normalized objective used by the solver.
        m.obj_tru = 1e4*(obj_m/Count_m+obj_cp/Count_cp+obj_cf/Count_cf)
        m.llh1 = m.count_m*log(m.obj_m) + m.count_cv*log(m.obj_cv) + (m.count_cr0+m.count_cr)*log(m.obj_cr)#m.count * log((obj_m+obj_cp+obj_cf0+obj_cf)/m.count)
        m.llh2 = Count_m*log(ob_m/Count_m) + Count_cp*log(ob_cp/Count_cp) + (count_cf0+Count_cf)*log((ob_cf0+ob_cf)/(count_cf0+Count_cf))
        
        return 1e4*(obj_m/Count_m + obj_cp/Count_cp + (obj_cf0+obj_cf)/(count_cf0+Count_cf))
    
    # pyomo model instance
    instance = model_construct_inter(data_stru, mode, theta, sim_opt, B_form, workflow_family=workflow_family)
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

    solver = _build_ipopt_solver()
    solver.options["max_iter"] = solver_max_iter
    #solver.options["halt_on_ampl_error"] = "yes" # option 1
    #solver.options["print_level"] = 1 # option 2
    try:
        results = solver.solve(instance,tee=True)
        term = results.solver.termination_condition
        if term == TerminationCondition.optimal:
            results.write()
        elif _is_infeasible_termination(term):
            print(f"Solver termination: {term}. Skipping repeated solve.")
            return None, None, None
        else:
            print(f"Solver termination: {term}, resolving...")
            try:
                solver.options["linear_solver"] = "mumps"
                solver.options["max_iter"] = solver_retry_max_iter
                results = solver.solve(instance,tee=True)
                assert results.solver.termination_condition == TerminationCondition.optimal, (
                        f"Solver failed: Non-optimal termination condition "
                        f"{results.solver.termination_condition}. ")
            except Exception as exc:
                print(f"Second solve failed: {exc}")
                return None, None, None
    except Exception as exc:
        print("Retrying with adjusted solver settings...")
        solver.options["linear_solver"] = "mumps"
        solver.options["max_iter"] = solver_retry_max_iter
        try:
            results = solver.solve(instance,tee=True)
            assert results.solver.termination_condition == TerminationCondition.optimal, (
                    f"Solver failed again: Non-optimal termination condition "
                    f"{results.solver.termination_condition}. "
                    "Try alternative initialization or further debugging.")
        except Exception as exc2:
            print(f"Adjusted solve failed: {exc2}")
            return None, None, None

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


def calc_FIM(data_stru, mode, theta=None, step=1e-8, formula='backward', B_form='single', workflow_family='DATA1'):
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
    formula = str(formula).lower()
    data_stru = _normalize_conductivity_measurements(data_stru)
    workflow_tag = str(workflow_family).upper()
    sim_opt = True

    def _solve_fim_case(label, theta_guess, *, sim_opt_flag):
        solve_kwargs = dict(
            theta=theta_guess,
            sim_opt=sim_opt_flag,
            B_form=B_form,
            workflow_family=workflow_family,
        )
        if workflow_tag == "DATA2":
            solve_kwargs["solver_max_iter"] = 3000
            solve_kwargs["solver_retry_max_iter"] = 5000
            return _safe_solve_case_multistart(
                label,
                solve_model,
                data_stru,
                mode,
                **solve_kwargs,
            )
        return solve_model(
            data_stru,
            mode,
            **solve_kwargs,
        )

    if theta is None:
        sim_opt = False
        fit_stru_p, sim_stru_p, sim_inter_p = _solve_fim_case(
            f"{workflow_tag or 'DATA'} FIM base fit",
            theta,
            sim_opt_flag=sim_opt,
        )
        if fit_stru_p is None or sim_inter_p is None:
            if workflow_tag == "DATA2":
                print("Skipping DATA2 FIM because the base fit did not converge.")
                return None
            raise RuntimeError("Base solve failed while building the FIM.")
        theta = fit_stru_p['parameters']
        sim_opt = True
    theta_p = copy.deepcopy(theta)
    theta_p_v=nested_dict_values(theta_p)
    theta_p1_v = []
    theta_p2_v = []

    # Build finite-difference perturbations for each parameter.
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
    
    # Prediction covariance is built from the paper measurement errors:
    # 0.01 g for mass, 3% for permeate concentration, and 0.3% for retentate.
    if sim_opt == True:
        fit_stru_p, sim_stru_p, sim_inter_p = _solve_fim_case(
            f"{workflow_tag or 'DATA'} FIM reference solve",
            theta_p,
            sim_opt_flag=sim_opt,
        )
        if fit_stru_p is None or sim_inter_p is None:
            if workflow_tag == "DATA2":
                print("Skipping DATA2 FIM because the reference simulation did not converge.")
                return None
            raise RuntimeError("Reference solve failed while building the FIM.")
    var_pred=[]
    for n_vial in range(data_stru['data_config']['n']):
        var_pred = np.append(var_pred,0.01 ** 2 * np.ones(len(sim_inter_p[n_vial]['mV'])))
        var_pred = np.append(var_pred,(0.03 * sim_inter_p[n_vial]['cV'])**2)
        var_pred = np.append(var_pred,(0.003 * sim_inter_p[n_vial]['cF'])**2)
    cov_pred = np.diag(var_pred)
    
    # The Jacobian is the sensitivity of all predictions with respect to all parameters.
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
        fit_stru_pk1, sim_stru_pk1, sim_inter_pk1 = _solve_fim_case(
            f"{workflow_tag or 'DATA'} FIM perturbation {i+1}a",
            theta_pk1,
            sim_opt_flag=sim_opt,
        )
        fit_stru_pk2, sim_stru_pk2, sim_inter_pk2 = _solve_fim_case(
            f"{workflow_tag or 'DATA'} FIM perturbation {i+1}b",
            theta_pk2,
            sim_opt_flag=sim_opt,
        )
        if sim_inter_pk1 is None or sim_inter_pk2 is None:
            if workflow_tag == "DATA2":
                print(f"Skipping DATA2 FIM because perturbation solve {i+1} did not converge.")
                return None
            raise RuntimeError(f"Perturbation solve {i+1} failed while building the FIM.")
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
    # Fisher Information Matrix:
    # FIM = J * Cov(y)^(-1) * J^T
    FIM = doe_stru['Jac'] @ np.linalg.inv(cov_pred) @ doe_stru['Jac'].T
    doe_stru['Jac'] = doe_stru['Jac'].tolist()
    doe_stru['FIM'] = FIM.tolist()

    # Eigenvalues and eigenvectors summarize identifiability and uncertainty directions.
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
        
        
def solve_model_B_fix(
    data_stru,
    mode,
    theta=None,
    sim_opt=False,
    B_form=1,
    sigma_fixed=True,
    LOUD=False,
    workflow_family='DATA1',
    solver_max_iter=3000,
    solver_retry_max_iter=5000,
):
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

    data_stru = _normalize_conductivity_measurements(data_stru)
    if str(workflow_family).upper() == "DATA2" and solver_max_iter == 3000:
        # DATA2 contains a few intentionally hard / locally infeasible cases.
        # Fail fast on those so the paper workflow can keep moving.
        solver_max_iter = 300
        solver_retry_max_iter = min(solver_retry_max_iter, 600)

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
        t_delay = _data1_time_origin(data_stru)
        TF_list = [data_stru['data_raw'][i]['time'][-1]-t_delay for i in range(data_stru['data_config']['n'])]
        TF_dict = dict(zip(m.n_vial,TF_list)) # unscaled time elapse for each vial 
        TI_list = [data_stru['data_raw'][i]['time'][0]-t_delay for i in range(data_stru['data_config']['n'])]
        TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time for each vial
    
        for n_vial in m.n_vial:
            t_meas = _data1_shifted_time(data_stru, n_vial-1)
            mv_meas = np.asarray(data_stru['data_raw'][n_vial-1]['mass'], dtype=float)
            cp_meas = np.asarray(data_stru['data_raw'][n_vial-1]['cV_avg'], dtype=float)
            cf_meas = np.asarray(data_stru['data_raw'][n_vial-1]['cF_exp'], dtype=float)
    
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
            if np.asarray(data_stru['data_raw'][n_vial-1]['cV_avg']).ndim == 0:
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
                if not np.all(np.isnan(cp_meas)):
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
            if not np.all(np.isnan(cf_meas)):
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
    instance = model_construct_inter(data_stru, mode, theta, sim_opt, B_form, workflow_family=workflow_family)
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

    solver = _build_ipopt_solver()
    solver.options["max_iter"] = solver_max_iter

    results = solver.solve(instance,tee=True)
    term = results.solver.termination_condition
    if _is_infeasible_termination(term):
        print(f"Solver termination: {term}. Skipping case.")
        return None, None, None
    assert term == TerminationCondition.optimal, "Optimization fail"
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


# -----------------------------------------------------------------------------
# Newer Pyomo tools
# -----------------------------------------------------------------------------
#
# The legacy functions above still work, but the helpers below let us use the
# current Pyomo ParmEst and Pyomo.DoE APIs in a simple way.


def _last_valid_value(values, default=np.nan):
    """Return the last non-NaN value from a list-like input."""
    if isinstance(values, (list, tuple, np.ndarray)):
        for item in reversed(values):
            if not (isinstance(item, float) and np.isnan(item)):
                return item
        return default
    if values is None:
        return default
    return values


def _abs_scale(value, fraction, floor):
    """Make a simple measurement error that is never zero."""
    try:
        scale = abs(float(value)) * fraction
    except Exception:
        scale = floor
    return max(scale, floor)


def _invert_monotone_1d(*, fwd_func, y_target, x_lo, x_hi, tol=1e-6, max_iter=200):
    """Invert a scalar conductivity model by bracketing and bisection."""
    import math

    def g(x):
        return float(fwd_func(float(x)) - y_target)

    lo = float(x_lo)
    hi = float(x_hi)
    for fac in [1.0, 2.0, 5.0, 10.0]:
        a = lo
        b = hi * fac
        grid = np.linspace(a, b, 80)
        vals = []
        for x in grid:
            try:
                gx = g(x)
            except Exception:
                continue
            if math.isnan(gx) or math.isinf(gx):
                continue
            vals.append((x, gx))
        for (x0, g0), (x1, g1) in zip(vals[:-1], vals[1:]):
            if g0 == 0:
                return float(x0)
            if g0 * g1 <= 0:
                left, right = x0, x1
                gl, gr = g0, g1
                for _ in range(max_iter):
                    mid = 0.5 * (left + right)
                    gm = g(mid)
                    if abs(gm) <= tol:
                        return float(mid)
                    if (gl < 0 < gm) or (gl > 0 > gm):
                        right, gr = mid, gm
                    else:
                        left, gl = mid, gm
                return float(0.5 * (left + right))
    raise ValueError("Target conductivity value could not be inverted.")


def _conductivity_to_concentration_series(
    *,
    cond_signal,
    temp_K,
    salt_name,
    model="variant_shedlovsky",
    model_params=None,
    output_units="mM",
):
    """Convert conductivity measurements to concentration using conductivity_paper.py.

    Paper equation path:
    - Shedlovsky equivalent conductivity model: \u03bb = f(c)
    - Variant Shedlovsky specific conductivity model: \u03ba = c \u00d7 \u03bb
    - Numerical inversion step: find c such that \u03ba(c) matches the measured conductivity

    The code uses the paper model forward and then inverts it point by point so
    the rest of the diafiltration workflow can continue to treat cF_exp as a
    concentration field.
    """
    # Load the conductivity model from the standalone paper file.
    cp = _load_conductivity_paper()
    # Start from a small built-in table and let the caller override it if needed.
    params = dict(CONDUCTIVITY_SALT_PARAMS_25C.get(str(salt_name), {}))
    if model_params:
        params.update(model_params)
    if not params:
        raise ValueError(
            f"No built-in conductivity parameters for salt '{salt_name}'. "
            "Pass model_params explicitly for this salt."
        )

    # Work with a lowercase model name so the caller can use either case.
    model = str(model).lower()
    cond_signal = np.asarray(cond_signal, dtype=float).reshape(-1)
    conc_out = np.full_like(cond_signal, np.nan, dtype=float)

    if model not in {"variant_shedlovsky", "shedlovsky"}:
        raise ValueError(
            "This refactor currently supports the paper's Shedlovsky / variant Shedlovsky path only."
        )

    epsilon = float(params["epsilon"])
    eta = float(params["eta"])
    a = float(params["a"])
    z_1 = int(params["z_1"])
    z_2 = int(params["z_2"])
    lambda_0_cation = float(params["lambda_0_cation"])
    lambda_0_anion = float(params["lambda_0_anion"])
    lambda_0 = float(params.get("lambda_0", lambda_0_cation + lambda_0_anion))

    def fwd(conc_M):
        # Forward model: concentration -> specific conductivity (mS/cm).
        return float(
            cp.variant_shedlovsky(
                [conc_M],
                temp_K,
                epsilon,
                eta,
                lambda_0,
                a,
                z_1,
                z_2,
                lambda_0_cation,
                lambda_0_anion,
            )[0]
        )

    for i, y in enumerate(cond_signal):
        if np.isnan(y):
            continue
        # Invert the paper conductivity curve to recover concentration from conductivity.
        conc_out[i] = _invert_monotone_1d(
            fwd_func=fwd,
            y_target=float(y / 1000.0),
            x_lo=0.0,
            x_hi=6.0,
        )

    if str(output_units).lower() in {"mm", "mmol/l", "mmolar"}:
        conc_out *= 1000.0
    return conc_out


def convert_experimental_conductivity_to_concentration(
    data_stru,
    *,
    output_units="mM",
    model="auto",
    model_params=None,
):
    """
    Convert retentate conductivity measurements to concentration in one place.

    This is the only data-conversion function that reaches into
    conductivity_paper.py. Everything else should work with the converted
    diafiltration dataset and stay inside the refactored code path.
    """
    # If the data structure is not the expected kind, leave it alone.
    if not isinstance(data_stru, dict) or not data_stru.get("conductivity_cF"):
        return data_stru

    # Read the salt name and temperature from the loaded experimental file.
    cfg = data_stru.get("data_config", {})
    salt_name = cfg.get("namec")
    temp_K = cfg.get("Temp", cfg.get("Temp_K"))
    if temp_K is None:
        raise ValueError("Temperature is required to convert conductivity to concentration.")
    if not salt_name:
        raise ValueError("Salt name is required to convert conductivity to concentration.")

    if model_params is None:
        model_params = {}

    # Default to the paper model unless the caller asks for something else.
    chosen_model = "variant_shedlovsky" if str(model).lower() == "auto" else str(model)
    for row in data_stru.get("data_raw", []):
        signal = row.get("cF_exp")
        if signal is None:
            continue
        # Keep the original conductivity values before replacing them.
        if "cF_exp_conductivity" not in row:
            row["cF_exp_conductivity"] = copy.deepcopy(signal)
        # Convert the conductivity series into concentration using the paper model.
        conc = _conductivity_to_concentration_series(
            cond_signal=signal,
            temp_K=float(temp_K),
            salt_name=str(salt_name),
            model=chosen_model,
            model_params=dict(model_params),
            output_units=output_units,
        )
        if np.asarray(signal).ndim == 0 or not isinstance(signal, (list, tuple, np.ndarray)):
            row["cF_exp"] = float(np.asarray(conc).reshape(-1)[0])
        else:
            row["cF_exp"] = np.asarray(conc, dtype=float).tolist()

    # Mark the file as converted so the rest of the workflow can tell what happened.
    data_stru["conductivity_cF_converted"] = True
    data_stru["conductivity_cF"] = False
    return data_stru


def _normalize_conductivity_measurements(data_stru, *, output_units="mM", model="auto", model_params=None):
    """Convert retentate conductivity measurements to concentration when needed.

    If the file says conductivity is present, we replace the fitting field with
    concentration so the rest of the diafiltration equations stay unchanged.

    This keeps the data-preprocessing step aligned with the paper workflow:
    - preserve the original conductivity values for traceability
    - convert cF_exp using the paper conductivity equation
    - feed the converted concentration into the existing first-principles model
    """
    return convert_experimental_conductivity_to_concentration(
        data_stru,
        output_units=output_units,
        model=model,
        model_params=model_params,
    )


def _theta_components(model):
    """Pick the model components that should be estimated."""
    theta_components = []

    if hasattr(model, "Lp"):
        theta_components.append(model.Lp)

    if hasattr(model, "B"):
        try:
            theta_components.extend([comp for _, comp in model.B.items()])
        except Exception:
            theta_components.append(model.B)

    for name in ("beta_0", "beta_1", "beta_2", "beta_3", "sigma", "S0", "S"):
        if hasattr(model, name):
            theta_components.append(getattr(model, name))

    return theta_components


def _label_parmest_model(model, data_stru, mode="DATA"):
    """Attach the labels ParmEst needs."""
    data_stru = _normalize_conductivity_measurements(data_stru)
    model.experiment_outputs = Suffix(direction=Suffix.LOCAL)
    model.measurement_error = Suffix(direction=Suffix.LOCAL)
    model.unknown_parameters = Suffix(direction=Suffix.LOCAL)
    model.experiment_inputs = Suffix(direction=Suffix.LOCAL)

    theta_components = _theta_components(model)
    model.unknown_parameters.update(
        (comp, ComponentUID(comp)) for comp in theta_components
    )

    final_t = model.tau.last()
    output_items = []

    for i in model.n_vial:
        row = data_stru["data_raw"][i - 1]
        cF_obs = _last_valid_value(row.get("cF_exp"))
        cV_obs = _last_valid_value(row.get("cV_avg"))
        mass_obs = _last_valid_value(row.get("mass"))

        output_items.append((model.cF[i, final_t], cF_obs, _abs_scale(cF_obs, 0.003, 0.003)))
        output_items.append((model.cV[i, final_t], cV_obs, _abs_scale(cV_obs, 0.03, 0.03)))

        if hasattr(model, "mF"):
            output_items.append((model.mF[i, final_t], mass_obs, _abs_scale(mass_obs, 0.01, 0.01)))

    for expr, obs, err in output_items:
        model.experiment_outputs[expr] = obs
        model.measurement_error[expr] = err

    return model


class DiafiltrationParmestExperiment(ParmestExperiment):
    """Small ParmEst wrapper around the diafiltration model."""

    def __init__(self, data_stru, mode="DATA", theta=None, B_form="single", nfe=300):
        super().__init__()
        self.data_stru = data_stru
        self.mode = mode
        self.theta = theta
        self.B_form = B_form
        self.nfe = nfe

    def get_labeled_model(self):
        """Build, discretize, and label the model for ParmEst."""
        model = model_construct_inter(
            self.data_stru,
            self.mode,
            theta=self.theta,
            sim_opt=False,
            B_form=self.B_form,
        )
        TransformationFactory("dae.finite_difference").apply_to(
            model, nfe=self.nfe, scheme="BACKWARD"
        )
        return _label_parmest_model(model, self.data_stru, mode=self.mode)


def build_parmest_experiments(data_structures, mode="DATA", theta=None, B_form="single", nfe=300):
    """Turn one data set or many data sets into a ParmEst experiment list."""
    if not isinstance(data_structures, list):
        data_structures = [data_structures]
    return [
        DiafiltrationParmestExperiment(data_stru, mode=mode, theta=theta, B_form=B_form, nfe=nfe)
        for data_stru in data_structures
    ]


def _normalize_experimental_source(source, *, campaign=None, variant="base", data_root=None):
    """Turn one experimental source into a loaded run instance and data dict."""
    if isinstance(source, ExperimentalRun):
        run = source
        return run, run.load_data()

    if isinstance(source, dict) and {"run_id", "root", "data_file"}.issubset(source.keys()):
        run = ExperimentalRun(
            run_id=str(source["run_id"]),
            root=Path(source["root"]),
            data_file=str(source["data_file"]),
            variant=str(source.get("variant", variant)),
            subdir=str(source.get("subdir", "")),
            fit_file=str(source.get("fit_file", "fit_stru.mat")),
        )
        return run, run.load_data()

    path = Path(source).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Experimental source not found: {path}")
    raw = loadmat(str(path))
    data_stru = raw.get("data_stru")
    if data_stru is None:
        raise ValueError(f"The file {path} does not contain a 'data_stru' entry.")
    return None, _normalize_conductivity_measurements(data_stru)


def estimate_parameters_with_parmest(
    data_structures,
    mode="DATA",
    theta=None,
    B_form="single",
    nfe=300,
    weighted=True,
    solver_options=None,
    tee=False,
):
    """Estimate parameters using the current ParmEst API."""
    exp_list = build_parmest_experiments(
        data_structures, mode=mode, theta=theta, B_form=B_form, nfe=nfe
    )
    obj_function = "SSE_weighted" if weighted else "SSE"
    pest = Estimator(exp_list, obj_function=obj_function, tee=tee, solver_options=solver_options)

    obj_val, theta_vals = pest.theta_est()
    cov = pest.cov_est()

    result = {
        "obj_val": obj_val,
        "theta_vals": theta_vals,
        "covariance": cov,
    }

    try:
        std = np.sqrt(np.diag(cov))
        result["std"] = std
        result["correlation"] = correlation_from_covariance(cov.values if hasattr(cov, "values") else cov)
    except Exception:
        result["std"] = None
        result["correlation"] = None

    return result


def compute_covariance_with_parmest(
    data_structures,
    mode="DATA",
    theta=None,
    B_form="single",
    nfe=300,
    method="finite_difference",
    solver="ipopt",
    tee=False,
    estimated_var=None,
):
    """Compute a covariance matrix with ParmEst."""
    exp_list = build_parmest_experiments(
        data_structures, mode=mode, theta=theta, B_form=B_form, nfe=nfe
    )
    obj_function = SSE_weighted if estimated_var is not None else SSE

    if estimated_var is None:
        pest = Estimator(exp_list, obj_function=obj_function, tee=tee)
        obj_val, theta_vals = pest.theta_est()
        cov = pest.cov_est(method=method)
    else:
        theta_vals = None
        obj_val = None
        cov = compute_covariance_matrix(
            exp_list,
            method=method,
            obj_function=obj_function,
            theta_vals=theta if theta is not None else {},
            step=0.02,
            solver=solver,
            tee=tee,
            estimated_var=estimated_var,
        )

    return {
        "obj_val": obj_val,
        "theta_vals": theta_vals,
        "covariance": cov,
    }


def compute_doe_metrics(experiment, step=0.001, objective_option="determinant", prior_FIM=None, solver=None, tee=False):
    """Use Pyomo.DoE to compute FIM metrics for a labeled experiment."""
    model = experiment.get_labeled_model()

    # Pyomo.DoE needs design inputs to run full DoE analysis.
    # If the model has not been extended with experiment_inputs yet, we stop
    # here with a clear message instead of failing later in a longer stack trace.
    if not hasattr(model, "experiment_inputs") or len(model.experiment_inputs) == 0:
        raise ValueError(
            "Pyomo.DoE needs experiment_inputs on the model. "
            "Add design variables to the model builder before calling compute_doe_metrics()."
        )

    doe = DesignOfExperiments(
        experiment=experiment,
        fd_formula="central",
        step=step,
        objective_option=objective_option,
        prior_FIM=prior_FIM,
        solver=solver,
        tee=tee,
    )
    fim = doe.compute_FIM(method="sequential")

    fim_np = np.asarray(fim, dtype=float)
    result = {
        "FIM": fim,
        "trace": float(np.trace(fim_np)),
        "det": float(np.linalg.det(fim_np)),
        "eig_val": np.linalg.eigvalsh(fim_np).tolist(),
    }
    try:
        result["covariance"] = np.linalg.inv(fim_np)
        result["std"] = np.sqrt(np.diag(result["covariance"])).tolist()
    except Exception:
        result["covariance"] = None
        result["std"] = None
    return result


def summarize_uncertainty(covariance):
    """Turn a covariance matrix into a small uncertainty summary."""
    cov_np = np.asarray(covariance, dtype=float)
    summary = {
        "covariance": cov_np.tolist(),
        "correlation": correlation_from_covariance(cov_np).tolist(),
        "std": np.sqrt(np.diag(cov_np)).tolist(),
    }
    return summary


def load_experimental_data(mat_file_path, results=None):
    """Stage 1: load one MATLAB file and normalize the measured fields."""
    if results is None:
        results = {}

    # Allow a single experiment or a batch of related experiments.
    if isinstance(mat_file_path, (list, tuple, set)):
        sources = list(mat_file_path)
    else:
        sources = [mat_file_path]

    loaded = []
    runs = []
    for source in sources:
        run, data_stru = _normalize_experimental_source(source)
        loaded.append(data_stru)
        runs.append(run)

    results["data_file"] = [str(src) for src in sources]
    results["data"] = loaded[0] if len(loaded) == 1 else loaded
    results["experiments"] = runs[0] if len(runs) == 1 else runs
    results["is_batch"] = len(loaded) > 1
    return results


def build_model(results, mode="DATA", workflow_family="DATA1", B_form="single"):
    """Stage 2: build the first-principles Pyomo model."""
    if "data" not in results:
        raise RuntimeError("Stage 1 must run before Stage 2.")

    data_payload = results["data"][0] if isinstance(results["data"], list) else results["data"]
    model = model_construct_inter(
        data_payload,
        mode,
        theta=None,
        sim_opt=False,
        B_form=B_form,
        workflow_family=workflow_family,
    )
    results["model"] = model
    results["model_settings"] = {
        "mode": mode,
        "workflow_family": workflow_family,
        "B_form": B_form,
        "is_batch": bool(results.get("is_batch")),
    }
    return results


def estimate_parameters(
    results,
    *,
    use_parmest=False,
    multistart=False,
    multistart_iterations=10,
    tee=False,
):
    """Stage 3: estimate parameters with ParmEst or the legacy solver path."""
    if "data" not in results:
        raise RuntimeError("Stage 1 must run before Stage 3.")

    settings = results.get("model_settings", {})
    mode = settings.get("mode", "DATA")
    workflow_family = settings.get("workflow_family", "DATA1")
    B_form = settings.get("B_form", "single")
    data_payload = results["data"]
    if isinstance(data_payload, list):
        data_payload_list = data_payload
    else:
        data_payload_list = [data_payload]

    if use_parmest or len(data_payload_list) > 1:
        pest_result = estimate_parameters_with_parmest(
            data_payload_list,
            mode=mode,
            theta=None,
            B_form=B_form,
            tee=tee,
        )
        theta_vals = pest_result.get("theta_vals")
        if hasattr(theta_vals, "to_dict"):
            theta_vals = theta_vals.to_dict()
        results["parmest"] = pest_result
        results["parameters"] = theta_vals or {}
        fit_stru, sim_stru, sim_inter = solve_model(
            data_payload_list[0],
            mode,
            theta=results["parameters"],
            sim_opt=False,
            B_form=B_form,
            workflow_family=workflow_family,
        )
    else:
        fit_stru, sim_stru, sim_inter = solve_model(
            data_payload_list[0],
            mode,
            theta=None,
            sim_opt=False,
            B_form=B_form,
            workflow_family=workflow_family,
        )
        results["parameters"] = fit_stru.get("parameters", {}) if fit_stru else {}

    results["fit_stru"] = fit_stru
    results["sim_stru"] = sim_stru
    results["sim_inter"] = sim_inter
    results["multistart"] = {
        "enabled": bool(multistart),
        "iterations": int(multistart_iterations),
    }
    return results


def quantify_uncertainty(
    results,
    *,
    method="fim",
    cov_method="finite_difference",
    fim_step=1e-8,
    fim_formula="backward",
):
    """Stage 4: quantify uncertainty with covariance or FIM."""
    if "parameters" not in results:
        raise RuntimeError("Stage 3 must run before Stage 4.")

    settings = results.get("model_settings", {})
    mode = settings.get("mode", "DATA")
    B_form = settings.get("B_form", "single")
    workflow_family = settings.get("workflow_family", "DATA1")
    data_stru = results["data"][0] if isinstance(results["data"], list) else results["data"]

    if method == "cov_est":
        pest_result = results.get("parmest")
        if pest_result and pest_result.get("covariance") is not None:
            results["uncertainty"] = summarize_uncertainty(pest_result["covariance"])
        else:
            results["uncertainty"] = {"status": "skipped", "reason": "ParmEst covariance not available."}

    elif method == "fim":
        fim = calc_FIM(
            data_stru,
            mode,
            theta=results["parameters"],
            step=fim_step,
            formula=fim_formula,
            B_form=B_form,
            workflow_family=workflow_family,
        )
        results["uncertainty"] = {
            "method": "calc_FIM",
            "FIM": fim.get("FIM"),
            "trace": fim.get("trace"),
            "det": fim.get("det"),
            "eig_val": fim.get("eig_val"),
            "covariance": fim.get("V"),
            "std": fim.get("std"),
        }

    else:
        results["uncertainty"] = {"status": "skipped", "reason": f"Unknown method {method!r}."}

    return results


def design_next_experiment(results, **_kwargs):
    """Stage 5: keep the hook, but skip MBDoE until design variables exist."""
    results["mbdoe"] = {
        "status": "skipped",
        "reason": "Design variables are not exposed yet; the workflow stays focused on fitting and reproduction.",
    }
    return results


def run_workflow(
    mat_file_path,
    *,
    mode="DATA",
    workflow_family="DATA1",
    B_form="single",
    use_parmest=False,
    multistart=False,
    multistart_iterations=10,
    uncertainty_method="fim",
    cov_method="finite_difference",
    fim_step=1e-8,
    fim_formula="backward",
    skip_mbdoe=True,
):
    """Run the compact stage-based workflow on one or many .mat files."""
    print("=" * 70)
    print(" Diafiltration workflow")
    print(f"   File: {mat_file_path}")
    print(f"   Family: {workflow_family}")
    print(f"   Mode: {mode}")
    print("=" * 70)

    results = {}
    results = load_experimental_data(mat_file_path, results)
    results = build_model(results, mode=mode, workflow_family=workflow_family, B_form=B_form)
    results = estimate_parameters(
        results,
        use_parmest=use_parmest,
        multistart=multistart,
        multistart_iterations=multistart_iterations,
    )
    results = quantify_uncertainty(
        results,
        method=uncertainty_method,
        cov_method=cov_method,
        fim_step=fim_step,
        fim_formula=fim_formula,
    )
    if not skip_mbdoe:
        results = design_next_experiment(results)

    print("\n[Done] Workflow finished. Results available:")
    for key in results:
        if not key.startswith("_"):
            print(f"   - {key}")
    return results


def _resolve_data_root(data_root=None):
    """Find the folder that stores the DATA1 paper inputs."""
    repo_root = Path(__file__).resolve().parents[1]
    if data_root is None:
        env_root = os.environ.get("DIAFILTRATION_DATA1_ROOT")
        if env_root:
            return Path(env_root).expanduser().resolve()
        return repo_root / "legacy" / "data1_matlab" / "data"
    return Path(data_root)


DATA1_RUN_REGISTRY = {
    "501.1": {
        "base": {"data_file": "data_stru-dataset501.1.mat", "subdir": "501.1"},
        "concpolar": {"data_file": "data_stru-dataset501.1.mat", "subdir": "501.1 concpolar"},
    },
    "501.11": {
        "base": {"data_file": "data_stru-dataset501.11.mat", "subdir": "501.11"},
        "concpolar": {"data_file": "data_stru-dataset501.11.mat", "subdir": "501.11 concpolar"},
    },
    "511.11": {
        "base": {"data_file": "data_stru-dataset511.11.mat", "subdir": "511.11"},
        "concpolar": {"data_file": "data_stru-dataset511.11.mat", "subdir": "511.11 concpolar"},
    },
    "511.12": {
        "base": {"data_file": "data_stru-dataset511.12.mat", "subdir": "511.12"},
        "concpolar": {"data_file": "data_stru-dataset511.12.mat", "subdir": "511.12 concpolar"},
    },
}


DATA2_RUN_REGISTRY = {
    "270511.123": {"data_file": "data_stru-dataset270511.123.mat"},
    "270611.121": {"data_file": "data_stru-dataset270611.121.mat"},
    "270711.121": {"data_file": "data_stru-dataset270711.121.mat"},
    "270511.221": {"data_file": "data_stru-dataset270511.221.mat"},
    "270511.321": {"data_file": "data_stru-dataset270511.321.mat"},
    "270511.421": {"data_file": "data_stru-dataset270511.421.mat"},
    "270511.921": {"data_file": "data_stru-dataset270511.921.mat"},
    "270511.521": {"data_file": "data_stru-dataset270511.521.mat"},
    "270511.621": {"data_file": "data_stru-dataset270511.621.mat"},
    "270511.721": {"data_file": "data_stru-dataset270511.721.mat"},
    "270511.821": {"data_file": "data_stru-dataset270511.821.mat"},
}


def get_experimental_run(campaign, run_id, variant="base", data_root=None):
    """Build one file-backed experimental run instance."""
    campaign = str(campaign or "DATA1").upper()
    run_id = str(run_id)
    variant = str(variant or "base").lower()
    root = Path(data_root) if data_root is not None else get_campaign_root(campaign)

    if campaign == "DATA1":
        family = DATA1_RUN_REGISTRY.get(run_id)
        if family is None:
            raise KeyError(f"Unknown DATA1 run '{run_id}'. Known runs: {list(DATA1_RUN_REGISTRY)}")
        variant_info = family.get(variant, family["base"])
        return ExperimentalRun(
            run_id=run_id,
            root=root,
            data_file=variant_info["data_file"],
            variant=variant,
            subdir=variant_info.get("subdir", run_id),
            fit_file="fit_stru.mat",
            contour_files=("contourdata-x_B-y_Lp.csv", "contourdata-x_sigma-y_Lp.csv"),
        )

    if campaign == "DATA2":
        family = DATA2_RUN_REGISTRY.get(run_id)
        if family is None:
            raise KeyError(f"Unknown DATA2 run '{run_id}'. Known runs: {list(DATA2_RUN_REGISTRY)}")
        return ExperimentalRun(
            run_id=run_id,
            root=root,
            data_file=family["data_file"],
            variant=variant,
            subdir="",
            fit_file="fit_stru.mat",
            contour_files=(),
        )

    raise ValueError(f"Unsupported campaign '{campaign}'.")


def _plot_heatmap_frame(df, x_col, y_col, z_col, ax=None, show_title=True, preface=False, cmap="viridis"):
    """Draw one contour-style heatmap from a tidy dataframe."""
    if ax is None:
        ax = plt.gca()

    # Pivot the table so the x/y grid becomes a matrix for plotting.
    grid = df.pivot_table(index=y_col, columns=x_col, values=z_col, aggfunc="mean")
    x_vals = grid.columns.to_numpy(dtype=float)
    y_vals = grid.index.to_numpy(dtype=float)
    z_vals = grid.to_numpy(dtype=float)

    mesh = ax.pcolormesh(x_vals, y_vals, z_vals, shading="auto", cmap=cmap)
    plt.colorbar(mesh, ax=ax)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    if show_title:
        title = z_col
        if preface:
            title = f"{title}"
        ax.set_title(title)
    return ax


def plot_contour(df, show_title=True, preface=False, save_path=None, cmap="viridis"):
    """
    Plot the contour-style heatmaps used in the DATA1 paper.

    If the dataframe contains the three objective columns from the paper,
    we draw one subplot for each objective.
    """
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)

    x_col = df.columns[0]
    y_col = df.columns[1]
    value_cols = [
        c
        for c in ["Obj_mass", "Obj_concentration", "Obj_retentate_concentration"]
        if c in df.columns
    ]
    if not value_cols:
        value_cols = [df.columns[2]]

    fig, axes = plt.subplots(1, len(value_cols), figsize=(5 * len(value_cols), 4), squeeze=False)
    for idx, z_col in enumerate(value_cols):
        _plot_heatmap_frame(df, x_col, y_col, z_col, ax=axes[0, idx], show_title=show_title, preface=preface, cmap=cmap)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig, axes


def plot_sim(sim_stru, color, linetype):
    """Plot simulation results, matching the DATA1 sigma-sensitivity notebook helper."""

    t_delay = np.asarray(sim_stru[0]["time"], dtype=float)[0]
    # plot mass prediction
    if plt.fignum_exists(1):
        plt.figure(1)
    else:
        plt.figure(1, figsize=(4,4))
    for i in sim_stru:
        t = np.asarray(i["time"], dtype=float)
        plt.plot((t-t_delay)/60, np.asarray(i["mV"], dtype=float), color, linestyle=linetype, linewidth=2, alpha=.8)

    # plot retentate concentration prediction
    if plt.fignum_exists(2):
        plt.figure(2)
    else:
        plt.figure(2, figsize=(4,4))
    for i in sim_stru:
        t = np.asarray(i["time"], dtype=float)
        plt.plot((t-t_delay)/60, np.asarray(i["cF"], dtype=float), color, linestyle=linetype, linewidth=2, alpha=.8)

    # plot permeate concentration prediction comparison
    if plt.fignum_exists(3):
        plt.figure(3)
    else:
        plt.figure(3, figsize=(4,4))
    for i in sim_stru:
        t = np.asarray(i["time"], dtype=float)
        plt.plot((t-t_delay)/60, np.asarray(i["cH"], dtype=float), color, linestyle=linetype, linewidth=2, alpha=.8)


def plot_sim_show(sig, colors):
    """Figure setup for different sigma values."""
    custom_lines = [Line2D([0], [0], color=colors[0],ls='--', lw=3),
                    Line2D([0], [0], color=colors[1],ls='-', lw=3),
                    Line2D([0], [0], color=colors[2],ls=':', lw=3)]
    legendname = ['$\\sigma$ = '+str(sig[0]),'$\\sigma$ = '+str(sig[1]),'$\\sigma$ = '+str(sig[2])]

    fig = plt.figure(1)
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Mass [g]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.ylim(bottom=0)
    fig.savefig('sigma_sensitivity-mass.png',dpi=300,bbox_inches='tight')

    fig = plt.figure(2)
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Retentate Conc. [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.legend(custom_lines,legendname,fontsize=15,loc='best')
    fig.savefig('sigma_sensitivity-reten_conc.png',dpi=300,bbox_inches='tight')

    fig = plt.figure(3)
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Permeate Conc. [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    fig.savefig('sigma_sensitivity-perme_conc.png',dpi=300,bbox_inches='tight')

    return custom_lines, legendname


def plot_conc_range(df_f, df_d, save_path=None, lg=True):
    """
    Plot experiment space - permeate concentration vs. retentate concentration.
    This matches the DATA1 notebook Fig. 3 helper.
    """
    fig = plt.figure(figsize=(4,4))

    # plot filtration
    for i in range(len(df_f.columns)//2):
        xstr = 'F'+str(i+1)+'_cf'
        ystr = 'F'+str(i+1)+'_cp'
        plt.plot(df_f[xstr], df_f[ystr], '^', markersize=8, alpha=.8, clip_on=False)

    # plot diafiltration
    for i in range(len(df_d.columns)//2):
        xstr = 'D'+str(i+1)+'_cf'
        ystr = 'D'+str(i+1)+'_cp'
        plt.plot(df_d[xstr], df_d[ystr], 's', markersize=8, alpha=.8, clip_on=False)

    # ghost point for legend
    plt.plot([],[],'k^',markersize=8,markerfacecolor='white',label='Filtration')
    plt.plot([],[],'ks',markersize=8,markerfacecolor='white',label='Diafiltration')

    plt.xlabel('Retentate [mM]',fontsize=16,fontweight='bold')
    plt.ylabel('Permeate [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    if lg:
        plt.legend(fontsize=15,loc='best')
    if save_path is None:
        save_path = 'concentration_range.png'
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    return fig, plt.gca()


def plot_cr_measure(data_stru, fit_stru, cr_pred, ybottom, lg=False, cond=True, save_path=None):
    """Plot retentate concentration vs. time from measurements and mass balance."""

    t_delay = _data1_time_origin(data_stru)

    def _fit_first_cF(fit_bundle):
        sim = fit_bundle.get("sim_stru") if isinstance(fit_bundle, dict) else fit_bundle
        if isinstance(sim, dict):
            sim = sim.get(0, sim.get("0", sim))
        return np.asarray(sim[0]["cF"], dtype=float).reshape(-1)[0]

    def _as_scalar_or_array(values):
        arr = np.asarray(values, dtype=float)
        if arr.ndim == 0:
            return float(arr.reshape(-1)[0])
        return arr

    fig = plt.figure(figsize=(4, 4))

    # initial point / starting concentration
    if cond:
        plt.plot(
            (np.asarray(data_stru["data_raw"][0]["time"], dtype=float)[0] - t_delay) / 60,
            _fit_first_cF(fit_stru),
            "ms",
            markersize=8,
            clip_on=False,
        )
    else:
        plt.plot(
            (np.asarray(data_stru["data_raw"][0]["time"], dtype=float)[0] - t_delay) / 60,
            data_stru["data_config"]["C_F0"],
            "go",
            markersize=8,
            clip_on=False,
        )

    # measurement points from conductivity probe / vial sampling
    for i in range(data_stru["data_config"]["n"]):
        row = data_stru["data_raw"][i]
        if cond:
            cF_exp = row["cF_exp"]
            if isinstance(cF_exp, (float, np.floating)):
                plt.plot(
                    (np.asarray(row["time"], dtype=float)[-1] - t_delay) / 60,
                    cF_exp,
                    "ms",
                    markersize=8,
                )
            elif len(cF_exp) > 1:
                plt.plot(
                    (np.asarray(row["time"], dtype=float) - t_delay) / 60,
                    np.asarray(cF_exp, dtype=float),
                    "ms",
                    markersize=8,
                )

    # calculated retentate concentration from mass balance
    if len(cr_pred) > 0:
        cr_pred = np.asarray(cr_pred, dtype=float).reshape(-1)
        for i in range(data_stru["data_config"]["n"]):
            plt.plot(
                (np.asarray(data_stru["data_raw"][i]["time"], dtype=float)[-1] - t_delay) / 60,
                cr_pred[i],
                "go",
                markersize=8,
            )

    plt.xlabel("Time [min]", fontsize=16, fontweight="bold")
    plt.ylabel("Concentration [mM]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in", top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=ybottom)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    if xticks:
        xticks[0].label1.set_visible(False)

    plt.plot([], [], "ms", markersize=8, label="Retentate \n (Measurements)")
    plt.plot([], [], "go", markersize=8, label="Retentate \n (Calculated)")
    if lg:
        plt.legend(fontsize=15, loc="best")

    if save_path is None:
        save_path = f"cr_measure-dat{data_stru['dataset']}.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig, ax


def plot_contour_sig_sen(
    data_stru,
    contour_sig_stru,
    Vmin,
    Vmax,
    Level,
    Colorbar_ticks,
    Manual_locations,
    name_append,
    filled=False,
    Jw_filter=True,
    bar=False,
    save_dir=None,
):
    """Plot the DATA1 sigma-sensitivity contour figures."""
    diaf_mode = False
    if "cf0" in contour_sig_stru:
        x = contour_sig_stru["cf0"]
        xlabelstr = r"$\mathbf{c_f}$(t=0) [mM]"
    else:
        x = contour_sig_stru["cd"]
        xlabelstr = r"$\mathbf{c_d}$ [mM]"
        diaf_mode = True
    y = contour_sig_stru["delp"]
    ylabelstr = r"$\mathbf{\Delta}$P [psi]"

    X, Y = np.meshgrid(x, y)
    axfontsize = 16

    R = 8.314e-5
    T = data_stru["data_config"]["Temp"]
    ni = data_stru["data_config"]["ni"]

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
        p = patches.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, color="gray", zorder=-10)
        ax.add_patch(p)

    fig = plt.figure(1, figsize=(4, 4))
    if filled:
        plt.contourf(X, Y, contour_sig_stru["range_mV"], 50, cmap="spring", vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X, Y, contour_sig_stru["range_mV"], Level[0], linewidths=3, colors="k")
    else:
        cp = plt.contour(X, Y, contour_sig_stru["range_mV"], Level[0], linewidths=2)
    try:
        plt.clabel(cp, inline=True, manual=Manual_locations[0], fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    except Exception:
        plt.clabel(cp, inline=True, fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color="gray", zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors="orangered")
        plt.setp(cg.collections, path_effects=[patheffects.withTickedStroke(angle=300, length=2)])
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight="bold")
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight="bold")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom=y[0])
    out_mass = Path(save_dir) / f"contour_sensitivity{name_append}-mass.png" if save_dir else Path(f"contour_sensitivity{name_append}-mass.png")
    fig.savefig(out_mass, dpi=300, bbox_inches="tight")

    fig = plt.figure(2, figsize=(4, 4))
    if filled:
        plt.contourf(X, Y, contour_sig_stru["range_cF"], 100, cmap="spring", vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X, Y, contour_sig_stru["range_cF"], Level[1], linewidths=3, colors="k")
    else:
        cp = plt.contour(X, Y, contour_sig_stru["range_cF"], Level[1], linewidths=2)
    try:
        plt.clabel(cp, inline=True, manual=Manual_locations[1], fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    except Exception:
        plt.clabel(cp, inline=True, fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color="gray", zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors="orangered")
        plt.setp(cg.collections, path_effects=[patheffects.withTickedStroke(angle=300, length=2)])
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight="bold")
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight="bold")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom=y[0])
    out_retentate = Path(save_dir) / f"contour_sensitivity{name_append}-retentate_conc.png" if save_dir else Path(f"contour_sensitivity{name_append}-retentate_conc.png")
    fig.savefig(out_retentate, dpi=300, bbox_inches="tight")

    fig = plt.figure(3, figsize=(4, 4))
    if filled:
        plt.contourf(X, Y, contour_sig_stru["range_cH"], 50, cmap="spring", vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X, Y, contour_sig_stru["range_cH"], Level[2], linewidths=3, colors="k")
    else:
        cp = plt.contour(X, Y, contour_sig_stru["range_cH"], Level[2], linewidths=2)
    try:
        plt.clabel(cp, inline=True, manual=Manual_locations[2], fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    except Exception:
        plt.clabel(cp, inline=True, fontsize=15, colors="k", fmt="%1.1f", zorder=2)
    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color="gray", zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors="orangered")
        plt.setp(cg.collections, path_effects=[patheffects.withTickedStroke(angle=300, length=2)])
    plt.xlabel(xlabelstr, fontsize=axfontsize, fontweight="bold")
    plt.ylabel(ylabelstr, fontsize=axfontsize, fontweight="bold")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom=y[0])
    out_permeate = Path(save_dir) / f"contour_sensitivity{name_append}-permeate_conc.png" if save_dir else Path(f"contour_sensitivity{name_append}-permeate_conc.png")
    fig.savefig(out_permeate, dpi=300, bbox_inches="tight")

    if filled and bar:
        cmap = mpl.cm.spring
        norm = mpl.colors.Normalize(vmin=Vmin, vmax=Vmax)
        fig, ax = plt.subplots(figsize=(12, 1))
        fig.subplots_adjust(bottom=0.5)
        mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation="horizontal", ticks=Colorbar_ticks, extend="max")
        ax.tick_params(labelsize=15)
        out_bar_h = Path(save_dir) / f"colorbar{name_append}-horizontal.png" if save_dir else Path(f"colorbar{name_append}-horizontal.png")
        fig.savefig(out_bar_h, dpi=300, bbox_inches="tight")
        fig, ax = plt.subplots(figsize=(1, 4))
        fig.subplots_adjust(left=0.5)
        mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation="vertical", ticks=Colorbar_ticks, extend="max")
        ax.tick_params(labelsize=15)
        out_bar_v = Path(save_dir) / f"colorbar{name_append}-vertical.png" if save_dir else Path(f"colorbar{name_append}-vertical.png")
        fig.savefig(out_bar_v, dpi=300, bbox_inches="tight")

    return [str(p) for p in [out_mass, out_retentate, out_permeate] + ([out_bar_h, out_bar_v] if filled and bar else [])]


def calib_curve_cond(calib_curve, save_path=None):
    """Plot the conductivity calibration curve used in the DATA1 notebook."""
    x = calib_curve.Conductivity.values
    y = calib_curve.Concentration.values

    fig = plt.figure(figsize=(4, 4))
    plt.plot(x, y, "bo", markersize=8)
    z = np.polyfit(x, y, 1)
    l = np.poly1d(z)
    correlation = np.corrcoef(x, y)[0, 1]
    r_squared = correlation**2
    plt.plot(x, l(x), "b:", linewidth=3, alpha=.7)
    idx = min(4, len(x) - 1)
    eqn = "y=%.3fx+%.2f \n R$\\mathbf{^{2}}$=%.4f" % (z[0], z[1], r_squared)
    plt.annotate(
        eqn,
        xy=(x[idx], y[idx]),
        xycoords="data",
        xytext=(-30, 60),
        weight="bold",
        textcoords="offset points",
        size=12,
        ha="center",
        va="center",
        bbox=dict(boxstyle="round", color="b", alpha=0.1),
    )
    xlabelstr = "Conductivity [$\\mathbf{\\mu}$S $\\mathbf{\\cdot}$ cm$\\mathbf{^{-1}}$]"
    ylabelstr = "Concentration [mM]"
    plt.xlabel(xlabelstr, fontsize=16, fontweight="bold")
    plt.ylabel(ylabelstr, fontsize=16, fontweight="bold")
    plt.xticks(fontsize=15, rotation=45)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in", top=True, right=True)
    if save_path is None:
        save_path = "calib_curve.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def run_data1_si_s2(data_root=None, save_dir=None):
    """Recreate DATA1 SI Figure S2 from the notebook's three concentration-ratio plots."""
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "si_fig_s2"
    save_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        {
            "label": "A",
            "dataset": "301.1",
            "data_path": root / "dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0" / "data_stru-dataset301.1.mat",
            "fit_path": root / "dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0" / "fit_stru-dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0.mat",
            "cr_path": root / "experiment space" / "Classical_analysis-dat301.1.csv",
            "cr_col": "cf",
            "ybottom": 4,
            "cond": False,
            "lg": False,
        },
        {
            "label": "B",
            "dataset": "501.1",
            "data_path": root / "data_stru-dataset501.1.mat",
            "fit_path": root / "501.1 concpolar" / "fit_stru.mat",
            "cr_path": root / "experiment space" / "Classical_analysis-dat501.1.csv",
            "cr_col": "cf",
            "ybottom": 4,
            "cond": True,
            "lg": False,
        },
        {
            "label": "C",
            "dataset": "511.12",
            "data_path": root / "data_stru-dataset511.12.mat",
            "fit_path": root / "511.12 concpolar" / "fit_stru.mat",
            "cr_path": root / "experiment space" / "diafiltration.csv",
            "cr_col": "D3_cf",
            "ybottom": 0,
            "cond": True,
            "lg": True,
        },
    ]

    panel_paths = []
    for case in cases:
        data_stru = loadmat(str(case["data_path"]))["data_stru"]
        fit_bundle = loadmat(str(case["fit_path"]))
        fit_stru = fit_bundle.get("fit_stru", fit_bundle)
        if case["dataset"] == "511.12":
            df = pd.read_csv(case["cr_path"], header=2)
        else:
            df = pd.read_csv(case["cr_path"], header=0)
        cr_pred = df[case["cr_col"]]
        out_path = save_dir / f"cr_measure-dat{case['dataset']}.png"
        fig, _ = plot_cr_measure(
            data_stru,
            fit_stru,
            cr_pred,
            case["ybottom"],
            lg=case["lg"],
            cond=case["cond"],
            save_path=out_path,
        )
        plt.close(fig)
        panel_paths.append(out_path)

    # Compose the three panels into a single SI-style page image.
    try:
        from PIL import Image, ImageDraw

        imgs = [Image.open(p).convert("RGB") for p in panel_paths]
        # Keep the notebook aspect ratio while placing three panels side-by-side.
        target_h = max(im.height for im in imgs)
        resized = []
        for im in imgs:
            scale = target_h / im.height
            new_size = (int(im.width * scale), target_h)
            resized.append(im.resize(new_size, Image.Resampling.LANCZOS))
        widths = [im.width for im in resized]
        canvas = Image.new("RGB", (sum(widths) + 80, target_h + 80), "white")
        draw = ImageDraw.Draw(canvas)
        x = 20
        for label, im in zip(["A", "B", "C"], resized):
            draw.text((x, 10), label, fill="black")
            canvas.paste(im, (x, 40))
            x += im.width + 20
        composite = save_dir / "figure_s2.png"
        canvas.save(composite)
        panel_paths.append(composite)
    except Exception:
        pass

    return [str(p) for p in panel_paths]


def _compose_data1_panel_sheet(row_paths, out_path):
    """Stack existing DATA1 panel images into one clean composite sheet."""
    try:
        from PIL import Image
    except Exception:
        return None

    valid_rows = [Path(p) for p in row_paths if Path(p).exists()]
    if not valid_rows:
        return None

    imgs = [Image.open(p).convert("RGB") for p in valid_rows]
    target_w = max(im.width for im in imgs)
    resized = []
    for im in imgs:
        scale = target_w / im.width
        resized.append(
            im.resize((target_w, int(im.height * scale)), Image.Resampling.LANCZOS)
        )

    gap = 20
    canvas_w = target_w + 20
    canvas_h = sum(im.height for im in resized) + gap * (len(resized) - 1) + 20
    page = Image.new("RGB", (canvas_w, canvas_h), "white")

    y = 10
    for im in resized:
        page.paste(im, (10, y))
        y += im.height + gap

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    page.save(out_path)
    return str(out_path)


def run_data1_si_panel_sheets(save_dir=None):
    """Build DATA1 SI panel-only composite sheets from generated panel PNGs."""
    save_dir = Path(save_dir) if save_dir is not None else Path.cwd()
    outputs = []

    composites = {
        "data1_si_reduced_diafiltration.png": [
            save_dir / "data511.11_fit.png",
            save_dir / "data511.11_fit_concentration.png",
        ],
        "data1_si_diafiltration_B.png": [
            save_dir / "data511.12_concpolar_contour_B.png",
            save_dir / "data511.11_concpolar_contour_B.png",
            save_dir / "data511.12_base_contour_B.png",
            save_dir / "data511.11_base_contour_B.png",
        ],
        "data1_si_reduced_filtration.png": [
            save_dir / "data501.11_fit.png",
            save_dir / "data501.11_fit_concentration.png",
        ],
        "data1_si_filtration_sigma.png": [
            save_dir / "data501.1_concpolar_contour_sigma.png",
            save_dir / "data501.11_concpolar_contour_sigma.png",
            save_dir / "data501.1_base_contour_sigma.png",
            save_dir / "data501.11_base_contour_sigma.png",
        ],
        "data1_si_filtration_B.png": [
            save_dir / "data501.1_concpolar_contour_B.png",
            save_dir / "data501.11_concpolar_contour_B.png",
            save_dir / "data501.1_base_contour_B.png",
            save_dir / "data501.11_base_contour_B.png",
        ],
    }

    for name, rows in composites.items():
        result = _compose_data1_panel_sheet(rows, save_dir / name)
        if result is not None:
            outputs.append(result)

    return outputs


def plot_conc_comparison(data_stru_f, fit_stru_f, data_stru_d, fit_stru_d, plot_pred=True, save_path=None):
    """Plot filtration and diafiltration concentration comparison figures."""
    t_delay_f = _data1_time_origin(data_stru_f)
    t_delay_d = _data1_time_origin(data_stru_d)

    def _first_conc_value(data_stru, key, fallback_key="cF_exp"):
        if key in data_stru and len(data_stru[key]) > 0:
            val = data_stru[key][0]
            if isinstance(val, (list, tuple, np.ndarray)):
                return _last_valid_value(val)
            return val
        for row in data_stru.get("data_raw", []):
            if fallback_key in row:
                val = row[fallback_key]
                if isinstance(val, (list, tuple, np.ndarray)):
                    return _last_valid_value(val)
                return val
        return np.nan

    def _plot_measurement_series(data_stru, i, x_shift, marker, color, cV_key="cV_avg", cF_key="cF_exp"):
        row = data_stru["data_raw"][i]
        time = np.asarray(row["time"], dtype=float) - x_shift

        def _plot_one(values, style):
            arr = np.asarray(values, dtype=float)
            if arr.ndim == 0 or arr.size == 1:
                x = np.array([time[-1]], dtype=float)
                y = np.array([float(arr.reshape(-1)[0])], dtype=float)
            else:
                n = min(len(time), len(arr))
                x = time[:n]
                y = arr[:n]
            mask = ~np.isnan(y)
            if np.any(mask):
                plt.plot(x[mask], y[mask], style, markersize=6)

        if cV_key in row:
            _plot_one(row[cV_key], marker.replace("s", "o"))
        if cF_key in row:
            _plot_one(row[cF_key], marker)

    fig = plt.figure(figsize=(4, 4))
    t0_f = _data1_time_origin(data_stru_f) - t_delay_f
    t0_d = _data1_time_origin(data_stru_d) - t_delay_d
    plt.plot(t0_f, _first_conc_value(data_stru_f, "cF_ICP"), "mv", markersize=6, clip_on=False)
    plt.plot(t0_d, _first_conc_value(data_stru_d, "cF_ICP"), "k^", markersize=6, clip_on=False)
    plt.plot(t0_f, float(np.asarray(fit_stru_f["sim_stru"][0]["cF"], dtype=float).reshape(-1)[0]), "bs", markersize=6, clip_on=False)
    plt.plot(t0_d, float(np.asarray(fit_stru_d["sim_stru"][0]["cF"], dtype=float).reshape(-1)[0]), "rs", markersize=6, clip_on=False)

    plt.plot([], [], "kv", markersize=6, clip_on=False, label="Filtration Retentate ICP")
    plt.plot([], [], "bs", markersize=6, clip_on=False, label="Filtration Retentate conductivity")
    plt.plot([], [], "bo", markersize=6, label="Filtration Vial ICP")
    plt.plot([], [], "k^", markersize=6, clip_on=False, label="Diafiltration Retentate ICP")
    plt.plot([], [], "rs", markersize=6, clip_on=False, label="Diafiltration Retentate conductivity")
    plt.plot([], [], "ro", markersize=6, label="Diafiltration Vial ICP")

    if plot_pred:
        plt.plot([], [], "g-.", linewidth=2, label="Filtration Retentate")
        plt.plot([], [], "g", linewidth=2, label="Diafiltration Retentate")
        plt.plot([], [], "k-.", linewidth=2, alpha=.6, label="Filtration Permeate")
        plt.plot([], [], "k", linewidth=2, alpha=.6, label="Diafiltration Permeate")

    for i in range(data_stru_f["data_config"]["n"]):
        _plot_measurement_series(data_stru_f, i, t_delay_f, "bs", "b")
        if plot_pred:
            t_pred = np.asarray(fit_stru_f["sim_stru"][i]["time"], dtype=float) - t_delay_f
            plt.plot(t_pred, np.asarray(fit_stru_f["sim_stru"][i]["cF"], dtype=float), "g-.", linewidth=2)
            plt.plot(t_pred, np.asarray(fit_stru_f["sim_stru"][i]["cH"], dtype=float), "k-.", linewidth=2, alpha=.6)

    for i in range(data_stru_d["data_config"]["n"]):
        _plot_measurement_series(data_stru_d, i, t_delay_d, "rs", "r")
        if plot_pred:
            t_pred = np.asarray(fit_stru_d["sim_stru"][i]["time"], dtype=float) - t_delay_d
            plt.plot(t_pred, np.asarray(fit_stru_d["sim_stru"][i]["cF"], dtype=float), "g", linewidth=2)
            plt.plot(t_pred, np.asarray(fit_stru_d["sim_stru"][i]["cH"], dtype=float), "k", linewidth=2, alpha=.6)

    plt.plot(data_stru_f["data_raw"][-1]["time"][-1] - t_delay_f, _first_conc_value(data_stru_f, "cF_ICP"), "kv", markersize=6, clip_on=False)
    plt.plot(data_stru_d["data_raw"][-1]["time"][-1] - t_delay_d, _first_conc_value(data_stru_d, "cF_ICP"), "k^", markersize=6, clip_on=False)
    plt.xlabel("Time [s]", fontsize=16)
    plt.ylabel("Concentration [mM]", fontsize=16)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.minorticks_on()
    plt.tick_params(direction="in", top=True, right=True)
    plt.tick_params(which="minor", direction="in", top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(fontsize=12, bbox_to_anchor=(1.1, 1.05), loc="upper left")
    if save_path is None:
        save_path = "concentration.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def run_data_analysis(data_root=None, datasets=None, save_dir=None):
    """
    Recreate the paper-style DATA1 analysis plots for one or more datasets.

    This is a simple helper that loads the saved paper files and calls the
    plotting functions above.
    """
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root
    save_dir.mkdir(parents=True, exist_ok=True)

    if datasets is None:
        datasets = [501.1, 501.11, 511.12, 511.11]

    outputs = []
    for dat in datasets:
        try:
            data_run = get_experimental_run("DATA1", str(dat), variant="base", data_root=root)
        except Exception:
            continue
        if not data_run.data_path.exists():
            continue
        data_stru = data_run.load_data()

        # Plot the main fit if the saved fit file is available.
        for variant_name in ["base", "concpolar"]:
            try:
                variant_run = get_experimental_run("DATA1", str(dat), variant=variant_name, data_root=root)
            except Exception:
                continue
            if variant_run.fit_path.exists():
                fit_bundle = variant_run.load_fit()
                sim_stru = fit_bundle.get("sim_stru") if isinstance(fit_bundle, dict) else None
                if sim_stru is None:
                    continue
                figs = plot_sim_comparison(data_stru, sim_stru, plot_pred=True, lg=False)
                if figs:
                    out = save_dir / f"data{dat}_fit.png"
                    figs[0].savefig(out, dpi=300, bbox_inches="tight")
                    outputs.append(str(out))
                    plt.close(figs[0])
                if len(figs) > 1:
                    out = save_dir / f"data{dat}_fit_concentration.png"
                    figs[1].savefig(out, dpi=300, bbox_inches="tight")
                    outputs.append(str(out))
                    plt.close(figs[1])
                if len(figs) > 2:
                    out = save_dir / f"data{dat}_fit_stirred.png"
                    figs[2].savefig(out, dpi=300, bbox_inches="tight")
                    outputs.append(str(out))
                    plt.close(figs[2])

            contour_b = variant_run.contour_path("contourdata-x_B-y_Lp.csv")
            if contour_b.exists():
                df = variant_run.load_contour("contourdata-x_B-y_Lp.csv")
                fig, _ = plot_contour(df)
                out = save_dir / f"data{dat}_{variant_name}_contour_B.png"
                fig.savefig(out, dpi=300, bbox_inches="tight")
                outputs.append(str(out))
                plt.close(fig)

            contour_s = variant_run.contour_path("contourdata-x_sigma-y_Lp.csv")
            if contour_s.exists():
                df = variant_run.load_contour("contourdata-x_sigma-y_Lp.csv")
                fig, _ = plot_contour(df)
                out = save_dir / f"data{dat}_{variant_name}_contour_sigma.png"
                fig.savefig(out, dpi=300, bbox_inches="tight")
                outputs.append(str(out))
                plt.close(fig)

        filtration_csv = root / "experiment space" / "filtration.csv"
        diafiltration_csv = root / "experiment space" / "diafiltration.csv"
        if filtration_csv.exists() and diafiltration_csv.exists():
            df_f = pd.read_csv(filtration_csv, header=2)
            df_d = pd.read_csv(diafiltration_csv, header=2)
            fig, _ = plot_conc_range(df_f, df_d)
            out = save_dir / f"data{dat}_conc_range.png"
            fig.savefig(out, dpi=300, bbox_inches="tight")
            outputs.append(str(out))
            plt.close(fig)

    calib_csv = root / "experiment space" / "conductivity_calibration.csv"
    if calib_csv.exists():
        calib_curve = pd.read_csv(calib_csv, header=0, skiprows=[1])
        fig = calib_curve_cond(calib_curve, save_path=save_dir / "calib_curve.png")
        outputs.append(str(save_dir / "calib_curve.png"))
        plt.close(fig)

    return outputs


def run_sigma_sensitivity(data_root=None, dataset=501.1, cf0=None, sigmas=None, save_dir=None):
    """
    Recreate the sigma-sensitivity plots from the DATA1 paper.

    The function reads the saved MAT files and overlays the simulations
    using the simple plotting helpers above.
    """
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "sigma_sensitivity"
    save_dir.mkdir(parents=True, exist_ok=True)

    if sigmas is None:
        sigmas = [0.1, 0.5, 0.9]

    if cf0 is None:
        cf0 = 5.2843 if float(dataset) == 501.1 else 15.2052

    colorstring = "rbg"
    linetypes = ["dashed", "solid", "dotted"]

    plt.close(1)
    plt.close(2)
    plt.close(3)

    for sigma_val, color, linetype in zip(sigmas, colorstring, linetypes):
        mat_name = f"sim_stru-dat{dataset} C_Fin{cf0}sig{sigma_val}.mat"
        mat_path = root / "sigma sensitivity" / mat_name
        if not mat_path.exists():
            continue
        sim_stru = loadmat(str(mat_path)).get("sim_stru")
        plot_sim(sim_stru, color, linetype)

    plot_sim_show(sigmas, colorstring[: len(sigmas)])
    out_mass = save_dir / "sigma_sensitivity-mass.png"
    out_reten = save_dir / "sigma_sensitivity-reten_conc.png"
    out_perm = save_dir / "sigma_sensitivity-perme_conc.png"
    plt.figure(1).savefig(out_mass, dpi=300, bbox_inches="tight")
    plt.figure(2).savefig(out_reten, dpi=300, bbox_inches="tight")
    plt.figure(3).savefig(out_perm, dpi=300, bbox_inches="tight")
    plt.close(1)
    plt.close(2)
    plt.close(3)
    return [str(out_mass), str(out_reten), str(out_perm)]


def run_data1_sigma_contours(data_root=None, save_dir=None):
    """Recreate the DATA1 sigma-sensitivity contour figures from the notebook."""
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "sigma_contours"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    cases = [
        (
            "501.1",
            root / "data_stru-dataset501.1.mat",
            root / "501.1 concpolar" / "contour_sig_stru-dat501.1.mat",
            0,
            3,
            [[1, 2, 3], [1, 2, 3], [1, 2, 3]],
            [0, 1, 2, 3],
            [[(20, 40), (35, 60), (50, 80)], [(20, 40), (35, 60), (50, 80)], [(20, 40), (35, 60), (50, 80)]],
            "",
            True,
        ),
        (
            "511.12-endvial1",
            root / "data_stru-dataset511.12.mat",
            root / "511.12 concpolar" / "contour_sig_stru-dat511.12-endvial1.mat",
            0,
            3,
            [[1, 1.5, 1.9], [1, 1.5, 1.9], [1, 1.5, 1.9]],
            [0, 1, 2, 3],
            [[(40, 60), (60, 100), (73.5, 110)], [(40, 60), (60, 100), (73.5, 110)], [(40, 60), (60, 100), (73.5, 110)]],
            "-endvial1",
            True,
        ),
        (
            "511.12-endvial5",
            root / "data_stru-dataset511.12.mat",
            root / "511.12 concpolar" / "contour_sig_stru-dat511.12-endvial5.mat",
            0,
            3,
            [[2, 5, 8], [2, 5, 8], [2, 5, 8]],
            [0, 1, 2, 3],
            [[(20, 40), (35, 60), (50, 80)], [(20, 40), (35, 60), (50, 80)], [(20, 40), (35, 60), (50, 80)]],
            "-endvial5",
            True,
        ),
        (
            "511.12-endvial10",
            root / "data_stru-dataset511.12.mat",
            root / "511.12 concpolar" / "contour_sig_stru-dat511.12-endvial10.mat",
            0,
            3,
            [[2, 5, 8], [2, 5, 8], [2, 5, 8]],
            [0, 1, 2, 3],
            [[(5, 10), (15, 30), (20, 45)], [(5, 10), (15, 30), (20, 45)], [(5, 10), (15, 30), (20, 45)]],
            "-endvial10",
            True,
        ),
    ]

    for label, data_path, contour_path, vmin, vmax, levels, ticks, manual, suffix, bar in cases:
        if not data_path.exists() or not contour_path.exists():
            continue
        data_id = "501.1" if "501.1" in label and "511" not in label else "511.12"
        try:
            variant = "concpolar"
            data_run = get_experimental_run("DATA1", data_id, variant=variant, data_root=root)
        except Exception:
            data_run = None
        data_stru = _normalize_conductivity_measurements(loadmat(str(data_path)).get("data_stru"))
        contour_sig_stru = loadmat(str(contour_path)).get("contour_sig_stru")
        outputs.extend(
            plot_contour_sig_sen(
            data_stru,
            contour_sig_stru,
            vmin,
            vmax,
            levels,
            ticks,
            manual,
            suffix,
            filled=True,
            bar=bar,
            save_dir=save_dir,
        )
        )
        plt.close("all")

    return outputs


def run_data1_concentration_comparison(data_root=None, save_dir=None):
    """Recreate the DATA1 concentration comparison figure from the notebook."""
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "concentration_comparison"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    pairs = [
        (
            root / "data_stru-dataset501.1.mat",
            root / "501.1 concpolar" / "fit_stru.mat",
            root / "data_stru-dataset511.12.mat",
            root / "511.12 concpolar" / "fit_stru.mat",
            save_dir / "concentration_comparison.png",
        ),
    ]

    for data_f_path, fit_f_path, data_d_path, fit_d_path, out_path in pairs:
        if not (data_f_path.exists() and fit_f_path.exists() and data_d_path.exists() and fit_d_path.exists()):
            continue
        data_run_f = ExperimentalRun("501.1", root, "data_stru-dataset501.1.mat", variant="base", subdir="501.1")
        data_run_d = ExperimentalRun("511.12", root, "data_stru-dataset511.12.mat", variant="concpolar", subdir="511.12 concpolar")
        data_stru_f = data_run_f.load_data()
        fit_stru_f = data_run_f.load_fit()
        data_stru_d = data_run_d.load_data()
        fit_stru_d = data_run_d.load_fit()
        fig = plot_conc_comparison(data_stru_f, fit_stru_f, data_stru_d, fit_stru_d, save_path=out_path)
        outputs.append(str(out_path))
        plt.close(fig)

    return outputs


def _show_saved_figures(paths):
    """Open saved figure files so the DATA1 plots appear like notebook output."""
    figure_paths = []
    for item in paths:
        path = Path(item)
        if path.suffix.lower() in {".png", ".jpg", ".jpeg", ".svg", ".pdf", ".eps"} and path.exists():
            figure_paths.append(path)

    for idx, fig_path in enumerate(sorted(figure_paths), start=1):
        try:
            image = plt.imread(fig_path)
        except Exception:
            continue
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(image)
        ax.axis("off")
        ax.set_title(f"{idx}: {fig_path.name}", fontsize=10)
    if figure_paths:
        plt.show()


def run_data1_notebook_workflow(data_root=None, save_dir=None, show=True):
    """
    Recreate the DATA1 notebook workflow using the utility-style model helpers
    and the paper-plot functions from DiafiltrationPaperPlots.ipynb.
    """
    root = _resolve_data_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data1_paper_figures"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    outputs.extend(run_data_analysis(data_root=root, save_dir=save_dir))
    experiment_space = root / "experiment space"
    filtration_csv = experiment_space / "filtration.csv"
    diafiltration_csv = experiment_space / "diafiltration.csv"
    if filtration_csv.exists() and diafiltration_csv.exists():
        df_f = pd.read_csv(filtration_csv, header=2)
        df_d = pd.read_csv(diafiltration_csv, header=2)
        fig, _ = plot_conc_range(df_f, df_d, save_path=save_dir / "concentration_range.png")
        outputs.append(str(save_dir / "concentration_range.png"))
        plt.close(fig)
    outputs.extend(run_sigma_sensitivity(data_root=root, dataset=501.1, save_dir=save_dir))
    outputs.extend(run_sigma_sensitivity(data_root=root, dataset=511.12, save_dir=save_dir))
    outputs.extend(run_data1_sigma_contours(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_concentration_comparison(data_root=root, save_dir=save_dir))
    outputs.extend(run_data1_figure2_workflow(data_root=root, save_dir=save_dir))
    outputs = sorted(set(str(p) for p in outputs))

    if show:
        _show_saved_figures(outputs)

    return outputs


def _resolve_data2_root(data_root=None):
    """Find the folder that stores the DATA2 paper inputs."""
    repo_root = Path(__file__).resolve().parents[1]
    if data_root is None:
        env_root = os.environ.get("DIAFILTRATION_DATA2_ROOT")
        if env_root:
            return Path(env_root).expanduser().resolve()
        return repo_root / "legacy" / "data1_matlab" / "data_library"
    return Path(data_root)


def _resolve_data2_validation_fig9_root():
    """Find the committed validation CSVs for DATA2 Fig. 9."""
    return Path(__file__).resolve().parents[1] / "UnifiedFramework" / "DATA3" / "docs" / "validation" / "simulation_validation_data_files" / "data2_main" / "fig9"


def _resolve_data2_table_baseline_root():
    """Find the committed DATA2 Tables 3-6 baseline CSVs."""
    return Path(__file__).resolve().parents[1] / "UnifiedFramework" / "DATA3" / "results" / "reproduction" / "20260306-mainpaper-table-baseline" / "tables"


def model_predictions(c_in, c_h, k0, k1, Pe):
    """Compute Js/Jw from the simple DATA2 regression model."""
    Kf = k1 * c_in + k0
    Kp = k1 * c_h + k0
    return (c_in * Kf * np.exp(Pe) - Kp * c_h) / (np.exp(Pe) - 1)


def plot_model_predictions(sim_data, k0, k1, Pe, save_path=None):
    """Draw the DATA2 Js/Jw surface and overlay the measured data."""
    c_in = np.asarray(sim_data["cIn"])
    c_h = np.asarray(sim_data["cH"])
    Js = np.asarray(sim_data["Js"])
    Jw = np.asarray(sim_data["Jw"])

    round_to = 5
    c_in_low = np.floor(np.min(c_in) / round_to) * round_to
    c_in_high = np.ceil(np.max(c_in) / round_to) * round_to
    c_h_low = np.floor(np.min(c_h) / round_to) * round_to
    c_h_high = np.ceil(np.max(c_h) / round_to) * round_to

    c_in_grid, c_h_grid = np.meshgrid(
        np.linspace(c_in_low, c_in_high, 100),
        np.linspace(c_h_low, c_h_high, 100),
    )

    Js_Jw_model_grid = model_predictions(c_in_grid, c_h_grid, k0, k1, Pe)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    ax.mouse_init()
    ax.plot_wireframe(c_in_grid, c_h_grid, Js_Jw_model_grid, color="blue", alpha=0.7, label="Regressed Model")
    ax.scatter(c_in, c_h, Js / Jw, color="red", label="Experimental Data", marker="o", s=50)
    ax.legend(fontsize=10)
    ax.set_xlabel("c$_{in}$ [mM]", fontsize=12, fontweight="bold")
    ax.set_ylabel("c$_{h}$ [mM]", fontsize=12, fontweight="bold")
    ax.set_zlabel("J$_s$/J$_w$", fontsize=12, fontweight="bold")
    ax.set_title(f"J$_s$/J$_w$ with Pe={Pe:.1f}", fontsize=14, fontweight="bold")
    ax.set_box_aspect([1, 1, 0.8])
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.8)
    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig, ax


def plot_pressure_change(
    data_pd,
    *,
    time_shift,
    marker_time,
    pressure_xlim,
    pressure_ylim,
    conc_ylim,
    fig_name,
    save_dir=None,
):
    """Recreate one DATA2 pressure-change plot from the notebook."""
    if isinstance(data_pd, (str, Path)):
        data_pd = pd.read_csv(data_pd, skiprows=[1])

    save_dir = Path(save_dir) if save_dir is not None else None
    fig = plt.figure(figsize=(4, 4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    ax1.plot((data_pd.Time - time_shift) / 60, data_pd.Pressure, "o", color="tab:orange", markersize=6)
    ax2.plot((data_pd.Time - time_shift) / 60, data_pd.Concentration, "s", color="m", markersize=6)
    ax1.plot([0, 0], [0, 100], "--", color="tab:blue", lw=2)
    marker_rel = (float(marker_time) - float(time_shift)) / 60.0
    ax1.plot([marker_rel, marker_rel], [0, 100], "--", color="tab:green", lw=2)

    ax1.set_xlabel("Time [min]", fontsize=16, fontweight="bold")
    ax1.set_ylabel("Applied Pressure [psi]", color="tab:orange", fontsize=16, fontweight="bold")
    ax1.tick_params(axis="x", which="major", direction="in", labelsize=12)
    ax1.tick_params(axis="x", which="minor", direction="in", labelsize=8)
    ax1.tick_params(axis="y", labelcolor="tab:orange", direction="in", labelsize=12)
    ax1.set_xlim(*pressure_xlim)
    ax1.set_ylim(*pressure_ylim)

    ax2.set_ylabel("Retentate [mM]", color="m", fontsize=16, fontweight="bold")
    ax2.tick_params(axis="y", labelcolor="m", direction="in", labelsize=12)
    ax2.set_ylim(*conc_ylim)

    out_path = Path(fig_name)
    if save_dir is not None:
        out_path = save_dir / out_path.name
    fig.savefig(out_path, dpi=600, bbox_inches="tight", transparent=True)
    return fig, out_path


def plot_startup_barplot(improvements=None, *, save_path=None):
    """Recreate the small DATA2 startup improvement bar chart."""
    if improvements is None:
        improvements = [138, -9]
    modes = ["Lag", "Overflow"]
    colors = ["#2CA02C" if imp > 0 else "#D62728" for imp in improvements]
    df = pd.DataFrame({"Mode": modes, "Improvement": improvements})

    fig = plt.figure(figsize=(6, 3))
    ax = fig.add_subplot(111)
    if sns is not None:
        ax = sns.barplot(data=df, x="Improvement", y="Mode", hue="Mode", palette=colors, width=0.6, legend=False)
    else:
        y_pos = np.arange(len(modes))
        ax.barh(y_pos, improvements, color=colors, height=0.6)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(modes)
    for i, (imp, mode) in enumerate(zip(improvements, modes)):
        text_color = "white" if imp > 0 else "black"
        plt.text(
            imp - (10 if imp > 0 else -12),
            i,
            f"{imp}%",
            ha="right" if imp > 0 else "left",
            va="center",
            fontsize=14,
            fontweight="bold",
            color=text_color,
        )
    plt.xlim([-40, 150])
    ax.grid(False)
    ax.axvline(x=0, color="k")
    ax.tick_params(axis="x", labelbottom="off")
    ax.tick_params(axis="y", direction="in", pad=-5)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(ax.get_yticklabels(), ha="left")
    plt.ylabel("")
    plt.xlabel("Information Improvement")
    plt.tight_layout()
    if save_path is None:
        save_path = "startup_barplot.png"
    fig.savefig(save_path, dpi=600, bbox_inches="tight")
    return fig, ax


def _sim_stru_to_dataframe(sim_stru):
    """Flatten the saved simulation structure into a tidy dataframe."""
    rows = []
    if isinstance(sim_stru, dict):
        vial_iter = sim_stru.items()
    else:
        vial_iter = enumerate(sim_stru, start=1)
    for vial_key, vial in vial_iter:
        if not isinstance(vial, dict):
            continue
        times = np.asarray(vial.get("time", []), dtype=float).reshape(-1)
        for idx, t in enumerate(times):
            row = {"vial": vial_key, "time": float(t)}
            for key in ("Js", "Jw", "cIn", "cH", "cF", "cV", "mF", "B"):
                if key in vial:
                    values = vial[key]
                    try:
                        arr = np.asarray(values, dtype=float).reshape(-1)
                        if arr.size == len(times):
                            row[key] = float(arr[idx])
                        elif arr.size > 0:
                            row[key] = float(arr[min(idx, arr.size - 1)])
                    except Exception:
                        pass
            rows.append(row)
    return pd.DataFrame(rows)


PE_LOWER_BOUND = 0.1
PE_INITIAL_VALUE = 15
PE_UPPER_BOUND = 20


def model_convection(sim_data, Pe_fixed_value=None, log_transform_Pe=False):
    """Fit the simple DATA2 convection model used in the visualization notebook."""
    data = sim_data.copy()
    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data)
    data = data.reset_index(drop=True)

    model = ConcreteModel()
    model.i = RangeSet(0, len(data) - 1)

    model.J_w = Param(model.i, initialize=lambda m, i: float(data.loc[i, "Jw"]), mutable=True)
    model.J_s = Param(model.i, initialize=lambda m, i: float(data.loc[i, "Js"]), mutable=True)
    model.c_in = Param(model.i, initialize=lambda m, i: float(data.loc[i, "cIn"]), mutable=True)
    model.c_h = Param(model.i, initialize=lambda m, i: float(data.loc[i, "cH"]), mutable=True)

    model.Js = Var(model.i, initialize=lambda m, i: float(data.loc[i, "Js"] / data.loc[i, "Jw"]))
    model.Kp = Var(model.i, initialize=1.0, bounds=(0.01, 1.2))
    model.Kf = Var(model.i, initialize=1.0, bounds=(0.01, 1.2))
    model.k1 = Var(initialize=1e-2)
    model.k0 = Var(initialize=1.0, bounds=(1e-4, 2))

    if Pe_fixed_value is not None:
        pe_lower = min(PE_LOWER_BOUND, float(Pe_fixed_value))
        pe_upper = max(PE_UPPER_BOUND, float(Pe_fixed_value))
    else:
        pe_lower = PE_LOWER_BOUND
        pe_upper = PE_UPPER_BOUND

    if log_transform_Pe:
        model.expPe = Var(initialize=np.exp(PE_INITIAL_VALUE), within=Reals, bounds=(np.exp(pe_lower), np.exp(pe_upper)))
    else:
        model.Pe = Var(initialize=PE_INITIAL_VALUE, within=Reals, bounds=(pe_lower, pe_upper))
        model.expPe = Expression(expr=exp(model.Pe))

    if Pe_fixed_value is not None:
        if log_transform_Pe:
            model.expPe.fix(np.exp(Pe_fixed_value))
        else:
            model.Pe.fix(float(Pe_fixed_value))

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
    model.SecondStageCost = Expression(rule=lambda m: sum((m.Js[i] / m.J_w[i] - m.J_s[i] / m.J_w[i]) ** 2 for i in m.i))
    model.Total_Cost_Objective = Objective(expr=model.FirstStageCost + model.SecondStageCost, sense=minimize)

    solver = _build_ipopt_solver()
    solver.options["halt_on_ampl_error"] = "yes"
    solver.options["acceptable_tol"] = 1e-8
    solver.solve(model, tee=False)

    if log_transform_Pe:
        pe_fitted = float(log(model.expPe.value))
    else:
        pe_fitted = float(value(model.Pe))

    theta_fit = {
        "k0": float(value(model.k0)),
        "k1": float(value(model.k1)),
        "Pe": pe_fitted,
    }
    return model, theta_fit


def linear_regression(data, Pe):
    """Run the simple no-intercept regression used in the DATA2 notebook."""

    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data)

    expPe = np.exp(Pe)
    Js = data["Js"].to_numpy(dtype=float)
    Jw = data["Jw"].to_numpy(dtype=float)
    cIn = data["cIn"].to_numpy(dtype=float)
    cH = data["cH"].to_numpy(dtype=float)

    y = Js / Jw
    x0 = (cIn * expPe - cH) / (expPe - 1)
    x1 = (cIn**2 * expPe - cH**2) / (expPe - 1)
    X = np.column_stack((x0, x1))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    y_hat = X @ beta
    rss = float(np.sum((y - y_hat) ** 2))
    dof = max(len(y) - X.shape[1], 1)
    mse = rss / dof
    xtx = X.T @ X
    try:
        cov = mse * np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        cov = mse * np.linalg.pinv(xtx)
    se = np.sqrt(np.diag(cov))
    results = {"params": beta, "mse_resid": mse, "bse": se, "residuals": y - y_hat}
    k0, k1 = beta[0], beta[1]
    se_k0, se_k1 = se[0], se[1]
    return k0, k1, mse, se_k0, se_k1, results


def plot_error_box(fit_struA, fit_struB, regime="concentrating", save_path=None):
    """Compare normalized residuals for two DATA2 fits."""
    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.labelsize": 16,
            "axes.titlesize": 16,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 14,
            "lines.linewidth": 2,
        }
    )

    res1 = np.concatenate(
        [
            np.array(fit_struA["res_std"]["res_m"]) / np.sqrt(len(fit_struA["res_std"]["res_m"])),
            np.array(fit_struA["res_std"]["res_cp"]) / np.sqrt(len(fit_struA["res_std"]["res_cp"])),
            np.array(fit_struA["res_std"]["res_cf"]) / np.sqrt(len(fit_struA["res_std"]["res_cf"])),
        ]
    )
    res2 = np.concatenate(
        [
            np.array(fit_struB["res_std"]["res_m"]) / np.sqrt(len(fit_struB["res_std"]["res_m"])),
            np.array(fit_struB["res_std"]["res_cp"]) / np.sqrt(len(fit_struB["res_std"]["res_cp"])),
            np.array(fit_struB["res_std"]["res_cf"]) / np.sqrt(len(fit_struB["res_std"]["res_cf"])),
        ]
    )

    res_dict = pd.DataFrame(
        {
            "Residuals": np.concatenate([res1, res2]),
            "Types": (
                ["Mass"] * len(fit_struA["res_std"]["res_m"])
                + ["Permeate"] * len(fit_struA["res_std"]["res_cp"])
                + ["Retentate"] * len(fit_struA["res_std"]["res_cf"])
            )
            * 2,
            "Solute Transport": (["Diffusion Only"] * len(res1) + ["Convection-Diffusion"] * len(res2)),
        }
    )

    fig, ax = plt.subplots(figsize=(7, 5))
    labels = ["Mass", "Permeate", "Retentate"]
    base_positions = np.arange(1, 4)
    offset = 0.18
    data_diff = [
        res_dict[(res_dict["Types"] == label) & (res_dict["Solute Transport"] == "Diffusion Only")]["Residuals"].to_numpy()
        for label in labels
    ]
    data_conv = [
        res_dict[(res_dict["Types"] == label) & (res_dict["Solute Transport"] == "Convection-Diffusion")]["Residuals"].to_numpy()
        for label in labels
    ]
    ax.boxplot(
        data_diff,
        positions=base_positions - offset,
        widths=0.3,
        patch_artist=True,
        boxprops=dict(facecolor="#4C72B0", alpha=0.75, linewidth=2),
        whiskerprops=dict(linewidth=2),
        capprops=dict(linewidth=2),
        medianprops=dict(linewidth=2, color="black"),
        flierprops=dict(marker="o", color="red", alpha=0.6),
    )
    ax.boxplot(
        data_conv,
        positions=base_positions + offset,
        widths=0.3,
        patch_artist=True,
        boxprops=dict(facecolor="#DD8452", alpha=0.75, linewidth=2),
        whiskerprops=dict(linewidth=2),
        capprops=dict(linewidth=2),
        medianprops=dict(linewidth=2, color="black"),
        flierprops=dict(marker="o", color="red", alpha=0.6),
    )

    ax.set_xlim(0.4, 3.6)
    ax.set_xticks(base_positions)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Weighted Residuals", fontsize=16)
    ax.set_ylabel("")
    ax.tick_params(axis="y", direction="in", pad=-5)
    ax.grid(axis="x", linestyle="--", alpha=0.5)

    plt.annotate(
        "Diffusion Only",
        xy=(0.9, 0.2),
        xycoords="data",
        weight="bold",
        size=16,
        ha="center",
        va="center",
        color="white",
        bbox=dict(boxstyle="round", color="#4C72B0", alpha=0.75),
    )
    plt.annotate(
        "Convection-Diffusion",
        xy=(0.9, 0.5),
        xycoords="data",
        weight="bold",
        size=16,
        ha="center",
        va="center",
        color="white",
        bbox=dict(boxstyle="round", color="#DD8452", alpha=0.75),
    )

    ax.legend(
        [
            Line2D([0], [0], color="#4C72B0", lw=8),
            Line2D([0], [0], color="#DD8452", lw=8),
        ],
        ["Diffusion Only", "Convection-Diffusion"],
        loc="best",
    )
    plt.tight_layout()
    if save_path is None:
        save_path = _figure_output_base(f"{regime}_residuals_boxplot")
    fig.savefig(str(save_path) + ".png", dpi=600, bbox_inches="tight")
    return fig, ax


def plot_weighted_residual_boxplot_from_csv(csv_path, regime="concentrating", save_path=None):
    """Render the DATA2 Fig. 9 residual boxplot directly from the committed CSV."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return None, None

    df = pd.read_csv(csv_path)
    if df.empty:
        return None, None

    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.labelsize": 16,
            "axes.titlesize": 16,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 14,
            "lines.linewidth": 2,
        }
    )

    labels = ["Mass", "Permeate", "Retentate"]
    base_positions = np.arange(len(labels), dtype=float)
    groups = [
        ("Diffusion Only", "#4C72B0", -0.18),
        ("Convection-Diffusion", "#DD8452", +0.18),
    ]

    fig, ax = plt.subplots(figsize=(7, 5))
    for transport_label, color, offset in groups:
        data = [
            df[(df["solute_transport"] == transport_label) & (df["residual_type"] == label)]["weighted_residual"].to_numpy(dtype=float)
            for label in labels
        ]
        bp = ax.boxplot(
            data,
            vert=False,
            positions=base_positions + offset,
            widths=0.28,
            patch_artist=True,
            whis=1.5,
            showfliers=True,
        )
        for box in bp["boxes"]:
            box.set(facecolor=color, alpha=0.75, linewidth=2)
        for median in bp["medians"]:
            median.set(color="black", linewidth=2)
        for part in ("whiskers", "caps"):
            for artist in bp[part]:
                artist.set(color=color, linewidth=2)
        for flier in bp["fliers"]:
            flier.set(marker="o", markerfacecolor="red", markeredgecolor="red", alpha=0.6, markersize=4)

    ax.set_xlim([-1.5, 1.5])
    ax.set_xlabel("Weighted Residuals", fontsize=16)
    ax.set_ylabel("")
    ax.tick_params(axis="y", direction="in", pad=-5)
    ax.set_yticks(base_positions)
    ax.set_yticklabels(labels, ha="left")
    ax.grid(axis="x", linestyle="--", alpha=0.5)
    ax.annotate(
        "Diffusion Only",
        xy=(0.9, 0.2),
        xycoords="axes fraction",
        weight="bold",
        size=16,
        ha="center",
        va="center",
        color="white",
        bbox=dict(boxstyle="round", color="#4C72B0", alpha=0.75),
    )
    ax.annotate(
        "Convection-Diffusion",
        xy=(0.9, 0.5),
        xycoords="axes fraction",
        weight="bold",
        size=16,
        ha="center",
        va="center",
        color="white",
        bbox=dict(boxstyle="round", color="#DD8452", alpha=0.75),
    )
    ax.legend(
        [
            Line2D([0], [0], color="#4C72B0", lw=8),
            Line2D([0], [0], color="#DD8452", lw=8),
        ],
        ["Diffusion Only", "Convection-Diffusion"],
        loc="best",
    )
    plt.tight_layout()
    if save_path is None:
        save_path = _figure_output_base(f"{regime}_residuals_boxplot")
    fig.savefig(str(save_path) + ".png", dpi=600, bbox_inches="tight")
    return fig, ax


def run_data2_table_bundle(data_root=None, save_dir=None):
    """Copy the committed DATA2 table baselines into the refactored output folder."""
    _ = _resolve_data2_root(data_root)
    baseline_root = _resolve_data2_table_baseline_root()
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_tables"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    copied = []
    for name in (
        "data2_table3_side_by_side.csv",
        "data2_table4_side_by_side.csv",
        "data2_table5_side_by_side.csv",
        "data2_table6_side_by_side.csv",
        "table_baseline_extraction_notes.txt",
    ):
        src = baseline_root / name
        if not src.exists():
            continue
        dst = save_dir / name
        df = pd.read_csv(src) if src.suffix.lower() == ".csv" else None
        if df is not None:
            df.to_csv(dst, index=False)
        else:
            dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
        copied.append({"artifact": name, "source": str(src), "output": str(dst)})
        outputs.append(str(dst))

    manifest = pd.DataFrame(copied)
    manifest_csv = save_dir / "data2_table_bundle_manifest.csv"
    manifest.to_csv(manifest_csv, index=False)
    outputs.append(str(manifest_csv))

    manifest_md = save_dir / "data2_table_bundle_manifest.md"
    lines = ["# DATA2 Table Bundle Manifest", ""]
    for row in copied:
        lines.append(f"- `{row['artifact']}` -> `{row['output']}`")
    manifest_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    outputs.append(str(manifest_md))

    return outputs


def run_data2_pressure_changes(data_root=None, save_dir=None):
    """Recreate the DATA2 pressure startup plots from the notebook."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_pressure_changes"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    cases = [
        ("pressure_change_lag.png", root / "NF270611.1_lag_pressure.csv", 185, 300, (-0.5, 5.5), (0, 68), (15, 30)),
        ("pressure_change_overflow.png", root / "NF270511.4_overflow_pressure.csv", 130, 150, (-0.4, 3.0), (0, 68), (13, 50)),
    ]
    for fig_name, csv_path, time_shift, marker_time, xlim, ylim_p, ylim_c in cases:
        if not csv_path.exists():
            continue
        fig, out_path = plot_pressure_change(
            csv_path,
            time_shift=time_shift,
            marker_time=marker_time,
            pressure_xlim=xlim,
            pressure_ylim=ylim_p,
            conc_ylim=ylim_c,
            fig_name=fig_name,
            save_dir=save_dir,
        )
        outputs.append(str(out_path))
        plt.close(fig)
    return outputs


def run_data2_calibration_plots(data_root=None, save_dir=None):
    """Recreate the DATA2 conductivity calibration curve."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_calibration"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    for csv_path in [root / "conductivity_calibration1.csv", root / "conductivity_calibration2.csv"]:
        if not csv_path.exists():
            continue
        calib_curve_data = pd.read_csv(csv_path, header=0, skiprows=[1])
        fig = calib_curve_cond(calib_curve_data, save_path=save_dir / "calib_curve.png")
        out_path = save_dir / "calib_curve.png"
        outputs.append(str(out_path))
        plt.close(fig)
    return outputs


def run_data2_model_demo(data_root=None, save_dir=None):
    """Recreate the DATA2 model-demo notebook plots."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_model_demo"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    data_stru_path = root / "data_stru-dataset270611.123.mat"
    sim_csv = root / "sim_data-dat270611.123.csv"
    mass_tc_path = root / "data_stru-dataset270611.123.mat"
    if data_stru_path.exists() and sim_csv.exists():
        data_stru = _normalize_conductivity_measurements(loadmat(str(data_stru_path)).get("data_stru"))
        sim_data = pd.read_csv(sim_csv)
        no_startup = sim_data.loc[sim_data["time"] > 120].reset_index(drop=True)

        # Notebook-style mass extrapolation figure.
        if isinstance(data_stru, dict) and "data_raw" in data_stru and len(data_stru["data_raw"]) > 3:
            row = data_stru["data_raw"][3]
            t_vals = np.array([(float(t)) / 60 for t in row["time"]], dtype=float)
            m_vals = np.array(row["mass"], dtype=float)
            f_m = interpolate.interp1d(t_vals, m_vals, kind="linear", fill_value="extrapolate")
            t_delay = data_stru["data_raw"][0]["time"][0]
            n_v0 = data_stru.get("data_config", {}).get("n_v0", 1)
            t_start = data_stru["data_raw"][n_v0 - 1]["time"][0] - t_delay
            fig = plt.figure(figsize=(4, 4))
            plt.plot(t_vals, m_vals, "r.", markersize=4)
            plt.plot((t_start + t_delay) / 60, 0, "b.", markersize=15)
            plt.plot(140 / 60, f_m(140 / 60), "k.", markersize=15)
            t_sim = np.linspace(140 / 60, (t_start + t_delay) / 60, 50)
            plt.plot(t_sim, f_m(t_sim), "r--", dashes=(8, 6), lw=1.5)
            plt.plot(np.linspace(140 / 60, 140 / 60, 50), f_m(t_sim), "-.", lw=1.5, color="gray")
            plt.plot([0, t_sim[0]], [f_m(t_sim[0])] * 2, "-.", lw=1.5, color="gray")
            plt.plot([], [], "r.", markersize=4, label="Measurements")
            plt.plot([], [], "r--", lw=1.5, label="Extrapolations")
            plt.xlabel("Time [min]", fontsize=16, fontweight="bold")
            plt.ylabel("Mass [g]", fontsize=16, fontweight="bold")
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.tick_params(direction="in")
            plt.ylim(-1.15, 1.15)
            plt.xlim(0, 10)
            plt.axhline(color="black", lw=1)
            plt.legend(fontsize=10, loc="best")
            out_path = save_dir / f"mass_tc-dat{data_stru['dataset']}.png"
            fig.savefig(out_path, dpi=300, bbox_inches="tight")
            outputs.append(str(out_path))
            plt.close(fig)

        # Notebook-style convection model fit and wireframes.
        model, theta_fit = model_convection(no_startup, Pe_fixed_value=None)
        fig, ax = plot_model_predictions(no_startup, theta_fit["k0"], theta_fit["k1"], theta_fit["Pe"], save_path=save_dir / "Js_predict1.png")
        outputs.append(str(save_dir / "Js_predict1.png"))
        plt.close(fig)

        model1, theta_fit1 = model_convection(sim_data, Pe_fixed_value=1)
        fig, ax = plot_model_predictions(sim_data, theta_fit1["k0"], theta_fit1["k1"], theta_fit1["Pe"], save_path=save_dir / "Js_predict.png")
        outputs.append(str(save_dir / "Js_predict.png"))
        plt.close(fig)

        model2, theta_fit2 = model_convection(no_startup, Pe_fixed_value=0.1)
        fig, ax = plot_model_predictions(no_startup, theta_fit2["k0"], theta_fit2["k1"], theta_fit2["Pe"], save_path=save_dir / "Jw_predict.png")
        outputs.append(str(save_dir / "Jw_predict.png"))
        plt.close(fig)

        model3, theta_fit3 = model_convection(no_startup, Pe_fixed_value=10)
        fig, ax = plot_model_predictions(no_startup, theta_fit3["k0"], theta_fit3["k1"], theta_fit3["Pe"], save_path=save_dir / "Js_predict0.png")
        outputs.append(str(save_dir / "Js_predict0.png"))
        plt.close(fig)

        # Notebook-style multistart sweep.
        Pe_values = np.logspace(-2.1, 1.4, 31)
        objective = np.zeros(len(Pe_values))
        k0_values = np.zeros(len(Pe_values))
        k1_values = np.zeros(len(Pe_values))
        for i, Pe in enumerate(Pe_values):
            model_sweep, _ = model_convection(sim_data, Pe_fixed_value=Pe)
            objective[i] = value(model_sweep.Total_Cost_Objective)
            k0_values[i] = value(model_sweep.k0)
            k1_values[i] = value(model_sweep.k1)
        fig = plt.figure(figsize=(4, 4))
        plt.plot(Pe_values, objective, "bo-", markersize=6)
        plt.xscale("log")
        plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
        plt.ylabel("Objective Function", fontsize=16, fontweight="bold")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tick_params(direction="in")
        plt.title("Multistart", fontsize=16, fontweight="bold")
        plt.grid(True)
        out_path = save_dir / "multistart_objective.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        outputs.append(str(out_path))
        plt.close(fig)

        fig = plt.figure(figsize=(4, 4))
        plt.plot(Pe_values, k0_values, "ro-", markersize=6)
        plt.xscale("log")
        plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
        plt.ylabel("k0", fontsize=16, fontweight="bold")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tick_params(direction="in")
        plt.title("k0 vs Peclet Number", fontsize=16, fontweight="bold")
        plt.grid(True)
        out_path = save_dir / "multistart_k0.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        outputs.append(str(out_path))
        plt.close(fig)

        fig = plt.figure(figsize=(4, 4))
        plt.plot(Pe_values, k1_values, "go-", markersize=6)
        plt.xscale("log")
        plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
        plt.ylabel("k1", fontsize=16, fontweight="bold")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tick_params(direction="in")
        plt.title("k1 vs Peclet Number", fontsize=16, fontweight="bold")
        plt.grid(True)
        out_path = save_dir / "multistart_k1.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        outputs.append(str(out_path))
        plt.close(fig)
    return outputs


def run_data2_partition_sensitivity(data_root=None, save_dir=None):
    """Recreate the DATA2 Peclet sensitivity plot."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_partition_sensitivity"
    save_dir.mkdir(parents=True, exist_ok=True)

    sim_csv = root / "sim_data-dat270611.123.csv"
    mat_file = root / "data_stru-dataset270611.123.mat"
    if not sim_csv.exists() or not mat_file.exists():
        return []

    _ = _normalize_conductivity_measurements(loadmat(str(mat_file)).get("data_stru"))
    sim_data = pd.read_csv(sim_csv)
    no_startup = sim_data.loc[sim_data["time"] > 120].reset_index(drop=True)

    Pe_values = np.logspace(-3, 2, 51)
    objective = np.zeros(len(Pe_values))
    k0_values = np.zeros(len(Pe_values))
    k1_values = np.zeros(len(Pe_values))
    k0_se = np.zeros(len(Pe_values))
    k1_se = np.zeros(len(Pe_values))

    for i, Pe in enumerate(Pe_values):
        k0_values[i], k1_values[i], objective[i], k0_se[i], k1_se[i], _ = linear_regression(no_startup, Pe)

    fig = plt.figure(figsize=(4.2, 12), constrained_layout=True)
    plt.subplot(3, 1, 1)
    plt.plot(Pe_values, objective, "b-", linewidth=3)
    plt.xscale("log")
    plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
    plt.ylabel("Mean Squared Error [mM$\\mathbf{^{2}}$]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.title("Regression Objective", fontsize=16, fontweight="bold", loc="left")
    plt.grid(True)

    plt.subplot(3, 1, 2)
    plt.plot(Pe_values, k0_values, "r-", linewidth=3)
    plt.xscale("log")
    plt.fill_between(Pe_values, k0_values - 2 * k0_se, k0_values + 2 * k0_se, color="red", alpha=0.3)
    plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
    plt.ylabel("h$\\mathbf{_0}$ [-]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.title("Partition Coefficient h$\\mathbf{_0}$", fontsize=16, fontweight="bold", loc="left")
    plt.grid(True)

    plt.subplot(3, 1, 3)
    plt.plot(Pe_values, k1_values, "g-", linewidth=3)
    plt.fill_between(Pe_values, k1_values - 2 * k1_se, k1_values + 2 * k1_se, color="green", alpha=0.3)
    plt.xscale("log")
    plt.xlabel("Peclet Number", fontsize=16, fontweight="bold")
    plt.ylabel("h$\\mathbf{_1}$ [mM$\\mathbf{^{-1}}$]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.title("Partition Coefficient h$\\mathbf{_1}$", fontsize=16, fontweight="bold", loc="left")
    plt.grid(True)

    out_path = save_dir / "partition_sensitivity.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return [str(out_path)]


def _safe_solve_case(label, solve_fn, *args, **kwargs):
    """Run one solver case and return None instead of stopping the workflow."""
    try:
        fit_stru, sim_stru, sim_inter = solve_fn(*args, **kwargs)
        if fit_stru is None or sim_stru is None or sim_inter is None:
            raise RuntimeError("Optimization fail")
        return fit_stru, sim_stru, sim_inter
    except Exception as exc:
        print(f"Skipping {label}: {exc}")
        return None, None, None


def _theta_variants_for_multistart(theta):
    """Generate a tiny set of nearby starting guesses for hard DATA2 fits."""
    if theta is None:
        return [None]
    if not isinstance(theta, dict):
        return [theta]

    variants = [copy.deepcopy(theta)]
    if "Lp" in theta:
        for scale in (0.85, 1.0, 1.15):
            trial = copy.deepcopy(theta)
            trial["Lp"] = float(theta["Lp"]) * scale
            variants.append(trial)
    if "B" in theta:
        for scale in (0.85, 1.0, 1.15):
            trial = copy.deepcopy(theta)
            trial["B"] = float(theta["B"]) * scale
            variants.append(trial)
    if "beta_0" in theta:
        for scale in (0.9, 1.0, 1.1):
            trial = copy.deepcopy(theta)
            trial["beta_0"] = float(theta["beta_0"]) * scale
            variants.append(trial)
    if "beta_1" in theta:
        for scale in (0.9, 1.0, 1.1):
            trial = copy.deepcopy(theta)
            trial["beta_1"] = float(theta["beta_1"]) * scale
            variants.append(trial)
    if "sigma" in theta:
        for candidate in (0.95, 1.0):
            trial = copy.deepcopy(theta)
            trial["sigma"] = candidate
            variants.append(trial)
    return variants


def _safe_solve_case_multistart(label, solve_fn, *args, **kwargs):
    """Try a tiny multistart sweep before giving up on a hard DATA2 fit."""
    theta = kwargs.get("theta")
    for trial_theta in _theta_variants_for_multistart(theta):
        trial_kwargs = dict(kwargs)
        trial_kwargs["theta"] = trial_theta
        fit_stru, sim_stru, sim_inter = _safe_solve_case(label, solve_fn, *args, **trial_kwargs)
        if fit_stru is not None and sim_stru is not None and sim_inter is not None:
            return fit_stru, sim_stru, sim_inter
    return None, None, None


def run_data2_model_error_visualization(data_root=None, save_dir=None, fast_mode=True):
    """Recreate the DATA2 residual comparison plot from the notebook."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_model_error_visualization"
    save_dir.mkdir(parents=True, exist_ok=True)

    data_stru = _normalize_conductivity_measurements(loadmat(str(root / "data_stru-dataset270511.123.mat")).get("data_stru"))
    outputs = []

    if fast_mode:
        bar_fig, _ = plot_startup_barplot([138, -9], save_path=save_dir / "startup_barplot.png")
        outputs.append(str(save_dir / "startup_barplot.png"))
        plt.close(bar_fig)
        fig9_root = _resolve_data2_validation_fig9_root()
        for regime in ("concentrating", "diluting"):
            csv_path = fig9_root / f"data2_main_fig9_{regime}_residuals.csv"
            out_path = save_dir / f"{regime}_residuals_boxplot"
            fig, _ = plot_weighted_residual_boxplot_from_csv(csv_path, regime=regime, save_path=out_path)
            if fig is not None:
                outputs.append(str(out_path) + ".png")
                plt.close(fig)
        return outputs
    else:
        fit_stru1, sim_stru1, _ = _safe_solve_case_multistart(
            "DATA2 model-error concentrating diffusion",
            solve_model,
            data_stru,
            "Lag",
            theta=None,
            sim_opt=False,
            B_form="single",
            workflow_family="DATA2",
        )
        fit_stru2, sim_stru2, _ = _safe_solve_case_multistart(
            "DATA2 model-error concentrating convection-diffusion",
            solve_model,
            data_stru,
            "Lag",
            theta={"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01, "sigma": 1.0, "S0": 0},
            sim_opt=False,
            B_form=1,
            workflow_family="DATA2",
        )
    if fit_stru1 is not None and fit_stru2 is not None:
        fig, _ = plot_error_box(fit_stru1, fit_stru2, regime="concentrating", save_path=save_dir / "concentrating")
        outputs.append(str(save_dir / "concentrating.png"))
        plt.close(fig)

    fit_stru3, sim_stru3, _ = _safe_solve_case_multistart(
        "DATA2 model-error diluting diffusion",
        solve_model_B_fix,
        data_stru,
        "DATA",
        theta={"Lp": 10.106582659197427, "B": 16.418266360453412, "sigma": 1.0},
        sim_opt=False,
        B_form="single",
        workflow_family="DATA2",
    )
    fit_stru4, sim_stru4, _ = _safe_solve_case_multistart(
        "DATA2 model-error diluting convection-diffusion",
        solve_model_B_fix,
        data_stru,
        "DATA",
        theta={"Lp": 10.167220041683919, "beta_0": 0.994192507828375, "beta_1": 0.013476905272029719, "sigma": 1.0},
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru3 is not None and fit_stru4 is not None:
        plot_sim_comparison(data_stru, sim_stru3, stirc_mass=False, plot_pred=True, lg=False, LOUD=True)
        plot_sim_comparison(data_stru, sim_stru4, stirc_mass=False, plot_pred=True, lg=False, LOUD=True)
        fig, _ = plot_error_box(fit_stru3, fit_stru4, regime="diluting", save_path=save_dir / "diluting")
        outputs.append(str(save_dir / "diluting.png"))
        plt.close(fig)

    return outputs


def run_data2_visualization(data_root=None, save_dir=None, show=False, fast_mode=True):
    """Recreate the main DATA2 visualization notebook figures."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_visualization"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    outputs.extend(run_data2_pressure_changes(data_root=root, save_dir=save_dir))
    outputs.extend(run_data2_calibration_plots(data_root=root, save_dir=save_dir))
    outputs.extend(run_data2_model_demo(data_root=root, save_dir=save_dir))
    outputs.extend(run_data2_partition_sensitivity(data_root=root, save_dir=save_dir))
    outputs.extend(run_data2_model_error_visualization(data_root=root, save_dir=save_dir, fast_mode=fast_mode))
    outputs = sorted(set(str(p) for p in outputs))
    if show:
        _show_saved_figures(outputs)
    return outputs


def run_cross_verification(data_root=None, save_dir=None):
    """Recreate the DATA2 cross-verification figure set."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_cross_verification"
    save_dir.mkdir(parents=True, exist_ok=True)

    fit_stru_base = {
        "parameters": {
            "Lp": 11.113241068147595,
            "beta_0": 1.0649294788103598,
            "beta_1": 0.015171421224603474,
            "sigma": 1.0,
        }
    }

    cases = [
        ("A.1", "data_stru-dataset270611.121.mat", True),
        ("A.2", "data_stru-dataset270711.121.mat", True),
        ("B.1", "data_stru-dataset270511.221.mat", True),
        ("B.2", "data_stru-dataset270511.321.mat", False),
        ("C.1", "data_stru-dataset270511.421.mat", False),
        ("C.2", "data_stru-dataset270511.921.mat", True),
        ("D.1", "data_stru-dataset270511.521.mat", True),
        ("D.2", "data_stru-dataset270511.621.mat", True),
        ("E.1", "data_stru-dataset270511.721.mat", False),
        ("E.2", "data_stru-dataset270511.821.mat", True),
    ]

    outputs = []
    for label, rel_path, sigma_fixed in cases:
        file_path = root / rel_path
        if not file_path.exists():
            continue
        data_stru = loadmat(str(file_path)).get("data_stru")
        fit_stru, sim_stru, sim_inter = _safe_solve_case(
            f"cross-verification case {label} ({rel_path})",
            solve_model_B_fix,
            data_stru,
            "DATA",
            theta=fit_stru_base["parameters"],
            sim_opt=False,
            B_form=1,
            sigma_fixed=sigma_fixed,
            workflow_family="DATA2",
        )
        if fit_stru is None:
            continue
        plot_sim_comparison(data_stru, sim_stru, stirc_mass=False, plot_pred=True, lg=False, LOUD=True)
        plt.close("all")
        outputs.extend(
            [
                str(Path("figures") / f"mass-dat{data_stru['dataset']}.png"),
                str(Path("figures") / f"concentration-dat{data_stru['dataset']}.png"),
            ]
        )
    return outputs


def run_data2_model_variations(data_root=None, save_dir=None, fast_mode=False):
    """Recreate the DATA2 model-variation and FIM calculations."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_model_variations"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    compute_fim = not fast_mode
    if fast_mode:
        print("DATA2 fast mode: skipping FIM calculations for model-variation cases.")
    elif not _has_casadi():
        print("CasADi is not available; DATA2 model-variation fits will run without simulator-based initialization.")

    # Lag - with time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.123.mat")).get("data_stru")
    mode = "Lag"
    theta = {"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01, "sigma": 1.0, "S0": 0}
    fit_stru_base1, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 LagM1",
        solve_model,
        data_stru,
        mode,
        theta,
        sim_opt=False,
        B_form=1,
        LOUD=False,
        workflow_family="DATA2",
    )
    if fit_stru_base1 is None:
        return outputs
    fit_path = save_dir / "LagM1-fit.json"
    store_json(str(fit_path), fit_stru_base1)
    outputs.append(str(fit_path))

    if compute_fim:
        doe_stru = calc_FIM(data_stru, mode, theta=fit_stru_base1["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
        if doe_stru is not None:
            fim_path = save_dir / "LagM1-FIM.json"
            store_json(str(fim_path), doe_stru)
            outputs.append(str(fim_path))

    # Lag - without time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.122.mat")).get("data_stru")
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 LagM2",
        solve_model,
        data_stru,
        mode,
        theta=fit_stru_base1["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "LagM2-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    # Lag truncated (DATA) - with time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.121.mat")).get("data_stru")
    mode_truc = "DATA"
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 LagM3",
        solve_model,
        data_stru,
        mode_truc,
        theta=fit_stru_base1["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "LagM3-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    if compute_fim:
        doe_stru = calc_FIM(data_stru, mode_truc, theta=fit_stru_base1["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
        if doe_stru is not None:
            fim_path = save_dir / "LagM3-FIM.json"
            store_json(str(fim_path), doe_stru)
            outputs.append(str(fim_path))

    # Lag truncated (DATA) - without time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.12.mat")).get("data_stru")
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 LagM4",
        solve_model,
        data_stru,
        mode_truc,
        theta=fit_stru_base1["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "LagM4-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    # Overflow - with time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.423.mat")).get("data_stru")
    mode = "Overflow"
    theta = {"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01, "sigma": 1.0, "S0": -0.1756665334051245}
    fit_stru_base2, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 OverflowM1",
        solve_model,
        data_stru,
        mode,
        theta,
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru_base2 is None:
        return outputs
    fit_path = save_dir / "OverflowM1-fit.json"
    store_json(str(fit_path), fit_stru_base2)
    outputs.append(str(fit_path))

    if compute_fim:
        doe_stru = calc_FIM(data_stru, mode, theta=fit_stru_base2["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
        if doe_stru is not None:
            fim_path = save_dir / "OverflowM1-FIM.json"
            store_json(str(fim_path), doe_stru)
            outputs.append(str(fim_path))

    # Overflow - without time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.422.mat")).get("data_stru")
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 OverflowM2",
        solve_model,
        data_stru,
        mode,
        theta=fit_stru_base2["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "OverflowM2-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    # Overflow truncated (DATA) - with time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.421.mat")).get("data_stru")
    mode_truc = "DATA"
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 OverflowM3",
        solve_model,
        data_stru,
        mode_truc,
        theta=fit_stru_base2["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "OverflowM3-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    if compute_fim:
        doe_stru = calc_FIM(data_stru, mode_truc, theta=fit_stru_base2["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
        if doe_stru is not None:
            fim_path = save_dir / "OverflowM3-FIM.json"
            store_json(str(fim_path), doe_stru)
            outputs.append(str(fim_path))

    # Overflow truncated (DATA) - without time correction
    data_stru = loadmat(str(root / "data_stru-dataset270511.42.mat")).get("data_stru")
    fit_stru, sim_stru, sim_inter = _safe_solve_case(
        "DATA2 OverflowM4",
        solve_model,
        data_stru,
        mode_truc,
        theta=theta,
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )
    if fit_stru is not None:
        fit_path = save_dir / "OverflowM4-fit.json"
        store_json(str(fit_path), fit_stru)
        outputs.append(str(fit_path))

    # Residual comparison box plots from the notebook.
    box_path = save_dir / "data2_residuals_boxplot"
    plot_error_box(fit_stru_base1, fit_stru_base2, regime="data2", save_path=box_path)
    outputs.append(str(box_path) + ".png")

    return outputs


def run_pre_B_dependence(data_root=None, save_dir=None):
    """Recreate the DATA2 pre-B dependence plots."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_pre_b_dependence"
    save_dir.mkdir(parents=True, exist_ok=True)

    mode = "DATA"
    allsim_stru = dict()

    data_files = [
        ("A.0", "data_stru-dataset270511.12.mat"),
        ("A.1", "data_stru-dataset270611.12.mat"),
        ("A.2", "data_stru-dataset270711.12.mat"),
        ("B.1", "data_stru-dataset270511.22.mat"),
        ("B.2", "data_stru-dataset270511.32.mat"),
        ("C.1", "data_stru-dataset270511.42.mat"),
        ("C.2", "data_stru-dataset270511.92.mat"),
        ("D.1", "data_stru-dataset270511.52.mat"),
        ("D.2", "data_stru-dataset270511.62.mat"),
        ("E.1", "data_stru-dataset270511.72.mat"),
        ("E.2", "data_stru-dataset270511.82.mat"),
    ]

    for key, rel_path in data_files:
        file_path = root / rel_path
        if not file_path.exists():
            continue
        data_stru = _normalize_conductivity_measurements(loadmat(str(file_path)).get("data_stru"))
        fit_stru, sim_stru, sim_inter = _safe_solve_case(
            f"DATA2 pre-B case {key}",
            solve_model,
            data_stru,
            mode,
            sim_opt=False,
            B_form="pervial",
            workflow_family="DATA2",
        )
        if fit_stru is None:
            continue
        allsim_stru[key] = sim_stru

    cmap1 = plt.get_cmap("tab10")
    cmap2 = plt.get_cmap("tab20")

    fig = plt.figure(figsize=(6, 4))
    fi = -1
    for key, sim_st in allsim_stru.items():
        if sim_st[0]["B"] is not None:
            co = cmap1(10) if fi == -1 else cmap2(fi)
            plt.plot([], [], color=co, linewidth=3, alpha=.8, label=key)
            for i in sim_st:
                start = 80 if i == 0 else 50
                plt.plot(np.array(sim_st[i]["cIn"][start:]), np.array(sim_st[i]["Js"][start:]) / np.array(sim_st[i]["Jw"][start:]), color=co, linewidth=2, alpha=.8)
            fi += 1
    plt.xlabel("Interface Concentration[mM]", fontsize=16, fontweight="bold")
    plt.ylabel("Js/Jw [mM]", fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(fontsize=10, loc="best")
    out1 = save_dir / "Js_Jw_cin.png"
    fig.savefig(out1, dpi=300, bbox_inches="tight")

    fig = plt.figure(figsize=(6, 4))
    fi = -1
    for key, sim_st in allsim_stru.items():
        if sim_st[0]["B"] is not None:
            co = cmap1(10) if fi == -1 else cmap2(fi)
            plt.plot([], [], color=co, linewidth=3, alpha=.8, label=key)
            for i in sim_st:
                plt.plot(sim_st[i]["cIn"], sim_st[i]["B"], color=co, linewidth=2, alpha=.8)
            fi += 1
    plt.xlabel("Interface Concentration[mM]", fontsize=16, fontweight="bold")
    plt.ylabel(r'B [$\mathbf{\mu}$m $\mathbf{\cdot}$ s$\mathbf{^{-1}}$]', fontsize=16, fontweight="bold")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    out2 = save_dir / "Bpervial.png"
    fig.savefig(out2, dpi=300, bbox_inches="tight")

    return [str(out1), str(out2)]


def run_data2_paper_reproduction(data_root=None, save_dir=None, fast_mode=True):
    """Run the DATA2 paper-style helpers in one place."""
    outputs = []
    outputs.extend(run_data2_visualization(data_root=data_root, save_dir=save_dir, fast_mode=fast_mode))
    outputs.extend(run_data2_table_bundle(data_root=data_root, save_dir=Path(save_dir) / "tables" if save_dir is not None else None))
    outputs.extend(run_cross_verification(data_root=data_root, save_dir=save_dir))
    outputs.extend(run_data2_model_variations(data_root=data_root, save_dir=save_dir, fast_mode=fast_mode))
    outputs.extend(run_pre_B_dependence(data_root=data_root, save_dir=save_dir))
    return outputs


def run_data2_notebook_workflow(data_root=None, save_dir=None, show=True, fast_mode=True):
    """
    Recreate the DATA2 notebook workflow using the utility-style model helpers
    and the plotting/analysis steps from DATA2_model_demo.ipynb and
    DATA2_visualization.ipynb.
    """
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else root / "data2_paper_figures"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = run_data2_paper_reproduction(data_root=root, save_dir=save_dir, fast_mode=fast_mode)
    outputs = sorted(set(str(p) for p in outputs))

    if show:
        _show_saved_figures(outputs)

    return outputs


def run_paper_reproduction(workflow_family="DATA1", data_root=None, save_dir=None):
    """Dispatch DATA1 or DATA2 paper reproduction from one simple switch."""
    workflow_family = str(workflow_family or "DATA1").upper()
    match workflow_family:
        case "DATA2":
            return run_data2_paper_reproduction(data_root=data_root, save_dir=save_dir)
        case _:
            outputs = run_data_analysis(data_root=data_root, save_dir=save_dir)
            root = _resolve_data_root(data_root)
            outputs.extend(run_sigma_sensitivity(data_root=root, dataset=501.1, save_dir=save_dir))
            outputs.extend(run_sigma_sensitivity(data_root=root, dataset=511.12, save_dir=save_dir))
            outputs.extend(run_data1_sigma_contours(data_root=root, save_dir=save_dir))
            outputs.extend(run_data1_concentration_comparison(data_root=root, save_dir=save_dir))
            return outputs


CAMPAIGN_REGISTRY = {
    "DATA1": {
        "root": Path("legacy") / "data1_matlab" / "data",
        "paper": "DATA1",
        "description": "Recreate the DATA1 paper plots and tables.",
    },
    "DATA2": {
        "root": Path("legacy") / "data1_matlab" / "data_library",
        "paper": "DATA2",
        "description": "Recreate the DATA2 paper plots and tables.",
    },
    "DATA3": {
        "root": Path("UnifiedFramework") / "DATA3" / "ExperimentalDataAnalysis" / "UnifiedCode",
        "paper": "DATA3",
        "description": "Placeholder for the future DATA3 campaign.",
    },
}


def register_campaign(name, root, paper=None, description=""):
    """Add one more campaign to the registry."""
    name = str(name or "").upper()
    if not name:
        raise ValueError("Campaign name cannot be empty.")
    CAMPAIGN_REGISTRY[name] = {
        "root": Path(root),
        "paper": str(paper or name).upper(),
        "description": description or f"Paper workflow for {name}.",
    }
    return CAMPAIGN_REGISTRY[name]


def get_campaign_info(campaign_name):
    """Return the registry entry for a campaign."""
    campaign_name = str(campaign_name or "DATA1").upper()
    return CAMPAIGN_REGISTRY.get(campaign_name, CAMPAIGN_REGISTRY["DATA1"])


def list_supported_campaigns():
    """Return the campaign names that this runner knows about."""
    return list(CAMPAIGN_REGISTRY.keys())


def get_campaign_root(campaign_name):
    """Return the folder that stores the inputs for one campaign."""
    return get_campaign_info(campaign_name)["root"]


def run_campaign(campaign_name="DATA1", save_dir=None):
    """Run a named campaign using the shared paper-reproduction helper."""
    campaign_name = str(campaign_name or "DATA1").upper()
    campaign = CAMPAIGN_REGISTRY.get(campaign_name)
    if campaign is None:
        campaign_name = "DATA1"
        campaign = CAMPAIGN_REGISTRY[campaign_name]
    return run_paper_reproduction(
        campaign["paper"],
        data_root=campaign["root"],
        save_dir=save_dir,
    )


def _format_table_value(value):
    """Format one table value so the rendered panel stays readable."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.4g}"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return str(value)


def _status_fill_color(status):
    """Choose a soft background color for a status cell."""
    status_text = str(status or "").strip().upper()
    if status_text == "PASS":
        return "#DFF0D8"
    if status_text == "FAIL":
        return "#F8D7DA"
    if status_text == "WARNING":
        return "#FFF3CD"
    if status_text == "MISSING_VALUE":
        return "#E2E3E5"
    return "white"


def _is_infeasible_termination(term):
    """Return True when the solver termination code indicates infeasibility."""
    term_text = str(term).replace("_", "").lower()
    return "infeasible" in term_text or "unbounded" in term_text


def render_table_panel(ax, table_source, title=None, color_status=True, fontsize=8):
    """Render one CSV table as a compact colored panel."""
    if isinstance(table_source, pd.DataFrame):
        df = table_source.copy()
    else:
        df = pd.read_csv(table_source)

    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=10)

    display_df = df.copy()
    display_df = display_df.apply(lambda col: col.map(_format_table_value))
    cell_text = display_df.values.tolist()
    col_labels = list(display_df.columns)

    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    table.scale(1, 1.35)
    try:
        table.auto_set_column_width(col=list(range(len(col_labels))))
    except Exception:
        pass

    status_idx = None
    for candidate in ("status", "Status"):
        if candidate in df.columns:
            status_idx = df.columns.get_loc(candidate)
            break

    # Header styling.
    for (row, col), cell in table.get_celld().items():
        cell.set_linewidth(0.6)
        if row == 0:
            cell.set_facecolor("#1F4E79")
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        cell.get_text().set_wrap(True)

    # Soft row shading to make the table easier to scan.
    for row in range(len(display_df)):
        if row % 2 == 1:
            for col in range(len(col_labels)):
                table[(row + 1, col)].set_facecolor("#F7F9FB")

        if color_status and status_idx is not None:
            fill = _status_fill_color(df.iloc[row, status_idx])
            for col in range(len(col_labels)):
                cell = table[(row + 1, col)]
                if cell.get_facecolor() == (1.0, 1.0, 1.0, 1.0) or row % 2 == 1:
                    cell.set_facecolor(fill if row % 2 == 0 else fill)
            if "status" in df.columns:
                status_col = df.columns.get_loc("status")
                table[(row + 1, status_col)].set_facecolor(fill)

    return table


def build_colored_composite_with_tables(image_paths, table_sources, out_path, title=None, ncols=2):
    """Build one composite sheet from colored images and rendered tables."""
    image_paths = [Path(p) for p in image_paths if Path(p).exists()]
    table_sources = [Path(p) if not isinstance(p, pd.DataFrame) else p for p in table_sources]
    table_sources = [p for p in table_sources if not isinstance(p, Path) or p.exists()]

    if not image_paths and not table_sources:
        raise ValueError("No image or table inputs were provided.")

    n_image_rows = math.ceil(len(image_paths) / ncols) if image_paths else 0
    n_table_rows = len(table_sources)
    total_rows = n_image_rows + n_table_rows
    if total_rows == 0:
        total_rows = 1

    fig_height = max(4.5, n_image_rows * 4.2 + n_table_rows * 2.6)
    fig = plt.figure(figsize=(ncols * 5.8, fig_height))
    if title:
        fig.suptitle(title, fontsize=18, fontweight="bold", y=0.995)

    gs = fig.add_gridspec(
        total_rows,
        ncols,
        height_ratios=[1.0] * max(n_image_rows, 1) + [0.95] * n_table_rows,
        hspace=0.35,
        wspace=0.18,
    )

    # Top section: image panels.
    for idx, img_path in enumerate(image_paths):
        row = idx // ncols
        col = idx % ncols
        ax = fig.add_subplot(gs[row, col])
        ax.imshow(plt.imread(str(img_path)))
        ax.axis("off")
        ax.set_title(Path(img_path).stem, fontsize=10, fontweight="bold")

    # Bottom section: table panels, each one spanning the full width.
    table_row_offset = n_image_rows
    for jdx, table_source in enumerate(table_sources):
        ax = fig.add_subplot(gs[table_row_offset + jdx, :])
        table_title = Path(table_source).stem if isinstance(table_source, Path) else None
        render_table_panel(ax, table_source, title=table_title, color_status=True, fontsize=8)

    if title:
        fig.subplots_adjust(top=0.93, hspace=0.35, wspace=0.18)
    else:
        fig.subplots_adjust(hspace=0.35, wspace=0.18)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(out_path)


def build_table_sheet(table_sources, out_path, title=None):
    """Render a set of tables as a dedicated sheet."""
    table_sources = [Path(p) if not isinstance(p, pd.DataFrame) else p for p in table_sources]
    table_sources = [p for p in table_sources if not isinstance(p, Path) or p.exists()]
    if not table_sources:
        raise ValueError("No table inputs were provided.")

    fig_height = max(3.5, len(table_sources) * 2.8)
    fig = plt.figure(figsize=(12.5, fig_height))

    gs = fig.add_gridspec(len(table_sources), 1, hspace=0.35)
    for idx, table_source in enumerate(table_sources):
        ax = fig.add_subplot(gs[idx, 0])
        render_table_panel(ax, table_source, title=None, color_status=True, fontsize=8)
        if isinstance(table_source, Path):
            table_label = table_source.stem
            if title:
                table_label = f"{title}: {table_label}"
            ax.set_title(table_label, fontsize=12, fontweight="bold", pad=8)

    fig.subplots_adjust(top=0.98, hspace=0.45)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(out_path)


def build_colored_paper_composites_with_tables(run_dir, campaign="DATA1", out_dir=None):
    """Build colored composite sheets plus table panels for one run directory."""
    run_dir = Path(run_dir)
    campaign = str(campaign or "DATA1").upper()
    fig_dir = run_dir / "figures"
    table_dir = run_dir / "tables"

    if out_dir is None:
        out_dir = run_dir / "composites_colored_tables"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    if campaign == "DATA1":
        image_groups = [
            (
                "DATA1_main_plus_table1.png",
                [
                    fig_dir / "data1_fig2_a_mass_measured.png",
                    fig_dir / "data1_fig2_b_retentate_measured.png",
                    fig_dir / "data1_fig2_c_vial_measured.png",
                    fig_dir / "data1_fig2_d_combined_measured.png",
                    fig_dir / "data1_fig3_model_overlay.png",
                    fig_dir / "data1_fig4_paper_style.png",
                ],
                [table_dir / "data1_table1_side_by_side.csv"],
            ),
            (
                "DATA1_main_plus_table2.png",
                [
                    fig_dir / "data1_fig5_paper_style.png",
                    fig_dir / "data1_fig6_paper_style.png",
                    fig_dir / "data1_direct_fig4_sigma_sensitivity.png" if (fig_dir / "data1_direct_fig4_sigma_sensitivity.png").exists() else fig_dir / "data1_fig4_direct_sigma_sensitivity.png",
                    fig_dir / "data1_fig4_1_jw_trajectories.png" if (fig_dir / "data1_fig4_1_jw_trajectories.png").exists() else fig_dir / "data1_fig4_2_js_trajectories.png",
                ],
                [table_dir / "data1_table2_side_by_side.csv"],
            ),
        ]
    else:
        image_groups = [
            (
                "DATA2_main_tables.png",
                sorted(
                    [
                        *fig_dir.glob("data2_*fig*.png"),
                        *fig_dir.glob("data2_*paper_style.png"),
                    ]
                ),
                sorted(table_dir.glob("data2_*side_by_side.csv")),
            ),
        ]

    for filename, images, tables in image_groups:
        out_path = out_dir / filename
        existing_images = [p for p in images if Path(p).exists()]
        existing_tables = [p for p in tables if Path(p).exists()]
        if not existing_images and not existing_tables:
            continue
        outputs.append(
            build_colored_composite_with_tables(
                existing_images,
                existing_tables,
                out_path,
                title=f"{campaign} composite with tables",
                ncols=2,
            )
        )
        if existing_tables:
            table_out = out_dir / f"{Path(filename).stem}_tables.png"
            outputs.append(
                build_table_sheet(
                    existing_tables,
                    table_out,
                    title=f"{campaign} tables",
                )
            )

    return outputs
