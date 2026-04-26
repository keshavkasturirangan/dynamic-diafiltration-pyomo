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
from matplotlib.lines import Line2D
from pathlib import Path
import idaes
import time
import copy
import os
import json
from sklearn.metrics import r2_score

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


def _figure_output_base(name):
    # Make sure the shared figures folder exists before saving anything there.
    os.makedirs(FIGURES_DIR, exist_ok=True)
    return os.path.join(FIGURES_DIR, name)


def _load_conductivity_paper():
    """Load the paper conductivity model directly from the repo."""
    try:
        # If the module is already importable, use it directly.
        import conductivity_paper as cp
        return cp
    except ModuleNotFoundError:
        # Otherwise, load it from the DATA3 folder by file path.
        import importlib.util
        from importlib.machinery import SourceFileLoader

        paper_path = (
            Path(__file__).resolve().parents[1]
            / "UnifiedFramework"
            / "DATA3"
            / "ExperimentalDataAnalysis"
            / "UnifiedCode"
            / "conductivity_paper.py"
        )
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


def plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=False):
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
    n_v0 = data_stru.get('data_config', {}).get('n_v0', 1)
    vial1 = n_v0 - 1
    t_start = 0

    def _is_sequence(value):
        return isinstance(value, (list, tuple, np.ndarray)) and not isinstance(value, (str, bytes))

    def _plot_measurement(x_vals, y_vals, style, **kwargs):
        """Plot a full series when lengths match, otherwise fall back to the last point."""
        x_list = list(x_vals)
        if not x_list:
            return
        if _is_sequence(y_vals):
            y_list = list(y_vals)
            if len(y_list) == 0:
                return
            if len(y_list) == len(x_list):
                plt.plot(x_list, y_list, style, **kwargs)
            else:
                plt.plot(x_list[-1], _last_valid_value(y_list), style, **kwargs)
        else:
            plt.plot(x_list[-1], y_vals, style, **kwargs)

    figures = []

    # plot mass data/prediction comparison
    fig = plt.figure(figsize=(4,4))
    for i in range(data_stru['data_config']['n']):
        _plot_measurement(
            [(float(t)-t_delay-t_start)/60 for t in data_stru['data_raw'][i]['time']],
            data_stru['data_raw'][i]['mass'],
            'r.',
            markersize=4,
        )
        if plot_pred and i >= vial1:
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
        fname = _figure_output_base('mass-dat'+str(data_stru['dataset']))
        fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
    figures.append(fig)

    # plot concentration data/prediction comparison
    fig = plt.figure(figsize=(4,4))

    for i in range(data_stru['data_config']['n']):
        _plot_measurement(
            [(float(t)-t_delay-t_start)/60 for t in data_stru['data_raw'][i]['time']],
            data_stru['data_raw'][i]['cV_avg'],
            'cs',
            markersize=6,
            clip_on=False,
        )
        _plot_measurement(
            [(float(t)-t_delay-t_start)/60 for t in data_stru['data_raw'][i]['time']],
            data_stru['data_raw'][i]['cF_exp'],
            'ms',
            markersize=6,
            clip_on=False,
        )
        if plot_pred:
            plt.plot([(time-t_start)/60 for time in sim_stru[i]['time']],
                     sim_stru[i]['cF'],
                     'g',linewidth=2)
            plt.plot([(time-t_start)/60 for time in sim_stru[i]['time']],
                     sim_stru[i]['cH'],
                     'r-',linewidth=2,alpha=.6)
            if not _is_sequence(data_stru['data_raw'][i]['cV_avg']):
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
        fname = _figure_output_base('concentration-dat'+str(data_stru['dataset']))
        fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
    figures.append(fig)
        
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
            fname = _figure_output_base('stirc_mass-dat'+str(data_stru['dataset']))
            fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
        figures.append(fig)

    return figures


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


def solve_model(data_stru, mode, theta=None, sim_opt=False, B_form='single', LOUD=False, workflow_family='DATA1'):
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
    sim_opt = True
    if theta is None:
        sim_opt = False
        fit_stru_p, sim_stru_p, sim_inter_p = solve_model(data_stru, mode, theta, sim_opt, B_form, workflow_family=workflow_family)
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
        fit_stru_p, sim_stru_p, sim_inter_p = solve_model(data_stru, mode, theta_p, sim_opt, B_form, workflow_family=workflow_family)
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
        fit_stru_pk1, sim_stru_pk1, sim_inter_pk1 = solve_model(data_stru, mode, theta_pk1, sim_opt, B_form, workflow_family=workflow_family)
        fit_stru_pk2, sim_stru_pk2, sim_inter_pk2 = solve_model(data_stru, mode, theta_pk2, sim_opt, B_form, workflow_family=workflow_family)
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
        
        
def solve_model_B_fix(data_stru, mode, theta=None, sim_opt=False, B_form=1, sigma_fixed=True, LOUD=False, workflow_family='DATA1'):
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

    The conversion is done by solving the paper model:
    conductivity = f(concentration)
    then numerically inverting it point by point.
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
        raise ValueError("This simple refactor currently supports variant_shedlovsky only.")

    epsilon = float(params["epsilon"])
    eta = float(params["eta"])
    a = float(params["a"])
    z_1 = int(params["z_1"])
    z_2 = int(params["z_2"])
    lambda_0_cation = float(params["lambda_0_cation"])
    lambda_0_anion = float(params["lambda_0_anion"])
    lambda_0 = float(params.get("lambda_0", lambda_0_cation + lambda_0_anion))

    def fwd(conc_M):
        # Forward model: c -> kappa.
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
        # Invert kappa = f(c) to recover concentration c.
        conc_out[i] = _invert_monotone_1d(
            fwd_func=fwd,
            y_target=float(y / 1000.0),
            x_lo=0.0,
            x_hi=6.0,
        )

    if str(output_units).lower() in {"mm", "mmol/l", "mmolar"}:
        conc_out *= 1000.0
    return conc_out


def _normalize_conductivity_measurements(data_stru, *, output_units="mM", model="auto", model_params=None):
    """Convert retentate conductivity measurements to concentration when needed.

    If the file says conductivity is present, we replace the fitting field with
    concentration so the rest of the diafiltration equations stay unchanged.
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
        # Convert the conductivity series into concentration.
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


def _resolve_data_root(data_root=None):
    """Find the folder that stores the DATA1 paper inputs."""
    if data_root is None:
        return Path("legacy") / "data1_matlab" / "data"
    return Path(data_root)


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


def plot_sim(sim_stru, colorstring="r", linestyle="solid", ax=None, show_permeate=False):
    """
    Plot a sigma-sensitivity simulation curve.

    The MATLAB scripts use this helper to overlay several simulation runs.
    """
    if ax is None:
        ax = plt.gca()

    # The saved MAT file may come back as a dict or a list of vial records.
    if isinstance(sim_stru, dict):
        vial_items = sim_stru.values()
    else:
        vial_items = sim_stru

    first = True
    for vial in vial_items:
        if not isinstance(vial, dict):
            continue
        t = np.asarray(vial.get("time", []), dtype=float)
        if t.size == 0:
            continue
        if "cF" in vial:
            ax.plot(t, vial["cF"], color=colorstring, linestyle=linestyle, linewidth=2, alpha=0.85)
        if show_permeate and "cH" in vial:
            ax.plot(t, vial["cH"], color=colorstring, linestyle="--", linewidth=1.5, alpha=0.6)
        if first:
            first = False

    ax.set_xlabel("Time")
    ax.set_ylabel("Concentration")
    return ax


def plot_sim_show(sigma, colorstring):
    """Create the small legend used in the sigma-sensitivity figures."""
    handles = []
    labels = []
    for s, c in zip(sigma, colorstring):
        handles.append(Line2D([0], [0], color=c, linewidth=2))
        labels.append(f"sigma = {s}")
    plt.legend(handles, labels, loc="best")
    return handles, labels


def plot_conc_range(df_f, df_d, save_path=None):
    """
    Plot the concentration range summary used in the paper notebook.

    This is intentionally lightweight: it just compares the available
    concentration columns from the filtration and diafiltration tables.
    """
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    # Use the columns if they exist, otherwise fall back to the second column.
    f_col = "cf" if "cf" in df_f.columns else df_f.columns[1]
    d_col = "cf" if "cf" in df_d.columns else df_d.columns[1]

    ax.plot(df_f.index, df_f[f_col], label=f"filtration: {f_col}", linewidth=2)
    ax.plot(df_d.index, df_d[d_col], label=f"diafiltration: {d_col}", linewidth=2)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Concentration")
    ax.legend(loc="best")
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig, ax


def plot_cr_measure(data_stru, fit_stru, cr_pred, ybottom, lg=False, cond=True, save_path=None):
    """
    Plot the concentration-ratio measurement figure from the notebook.

    The inputs are kept flexible so the function can work with the raw
    paper tables or with fitted model outputs.
    """
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    # Plot the predicted curve from the paper table.
    cr_pred = np.asarray(cr_pred, dtype=float)
    ax.plot(np.arange(len(cr_pred)), cr_pred, "k-", linewidth=2, label="paper prediction")

    # Overlay the last available retentate concentration for each vial.
    if isinstance(data_stru, dict) and "data_raw" in data_stru:
        y_meas = []
        for row in data_stru["data_raw"]:
            val = row.get("cF_exp", np.nan)
            if isinstance(val, (list, np.ndarray)):
                val = _last_valid_value(val)
            y_meas.append(val)
        ax.plot(np.arange(len(y_meas)), y_meas, "ms", markersize=5, label="measurements")

    ax.set_ylim(bottom=ybottom)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Concentration ratio / concentration")
    if lg:
        ax.legend(loc="best")
    if cond:
        ax.set_title("Concentration ratio comparison")
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig, ax


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
        data_stru_path = root / f"data_stru-dataset{dat}.mat"
        if not data_stru_path.exists():
            continue
        data_stru = _normalize_conductivity_measurements(loadmat(str(data_stru_path)).get("data_stru"))

        # Plot the main fit if the saved fit file is available.
        for variant in [f"{dat}", f"{dat} concpolar"]:
            fit_path = root / variant / "fit_stru.mat"
            if fit_path.exists():
                fit_bundle = loadmat(str(fit_path)).get("fit_stru")
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

            contour_b = root / variant / "contourdata-x_B-y_Lp.csv"
            if contour_b.exists():
                df = pd.read_csv(contour_b)
                fig, _ = plot_contour(df)
                out = save_dir / f"data{dat}_{variant.replace(' ', '_')}_contour_B.png"
                fig.savefig(out, dpi=300, bbox_inches="tight")
                outputs.append(str(out))
                plt.close(fig)

            contour_s = root / variant / "contourdata-x_sigma-y_Lp.csv"
            if contour_s.exists():
                df = pd.read_csv(contour_s)
                fig, _ = plot_contour(df)
                out = save_dir / f"data{dat}_{variant.replace(' ', '_')}_contour_sigma.png"
                fig.savefig(out, dpi=300, bbox_inches="tight")
                outputs.append(str(out))
                plt.close(fig)

        classical = root / "experiment space" / f"Classical_analysis-dat{dat}.csv"
        if classical.exists():
            df = pd.read_csv(classical)
            fig, _ = plot_conc_range(df, df)
            out = save_dir / f"data{dat}_conc_range.png"
            fig.savefig(out, dpi=300, bbox_inches="tight")
            outputs.append(str(out))
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
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    for sigma_val, color in zip(sigmas, colorstring):
        mat_name = f"sim_stru-dat{dataset} C_Fin{cf0}sig{sigma_val}.mat"
        mat_path = root / "sigma sensitivity" / mat_name
        if not mat_path.exists():
            continue
        sim_stru = loadmat(str(mat_path)).get("sim_stru")
        plot_sim(sim_stru, color, "solid", ax=ax)

    plot_sim_show(sigmas, colorstring[: len(sigmas)])
    out = save_dir / f"sigma_sensitivity-dat{dataset}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(out)


def _resolve_data2_root(data_root=None):
    """Find the folder that stores the DATA2 paper inputs."""
    if data_root is None:
        return Path("legacy") / "data1_matlab" / "data_library"
    return Path(data_root)


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


def run_data2_model_error_visualization(data_root=None, save_dir=None):
    """Recreate the DATA2 residual comparison plot from the notebook."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_model_error_visualization"
    save_dir.mkdir(parents=True, exist_ok=True)

    data_stru = _normalize_conductivity_measurements(loadmat(str(root / "data_stru-dataset270511.123.mat")).get("data_stru"))
    mode = "DATA"

    fit_stru_base = {
        "parameters": {
            "Lp": 10.106582659197427,
            "B": 16.418266360453412,
            "sigma": 1.0,
        }
    }
    fit_stru3, sim_stru3, sim_inter3 = solve_model_B_fix(
        data_stru,
        mode,
        theta=fit_stru_base["parameters"],
        sim_opt=False,
        B_form="single",
        workflow_family="DATA2",
    )

    fit_stru_base = {
        "parameters": {
            "Lp": 10.167220041683919,
            "beta_0": 0.994192507828375,
            "beta_1": 0.013476905272029719,
            "sigma": 1.0,
        }
    }
    fit_stru4, sim_stru4, sim_inter4 = solve_model_B_fix(
        data_stru,
        mode,
        theta=fit_stru_base["parameters"],
        sim_opt=False,
        B_form=1,
        workflow_family="DATA2",
    )

    plot_error_box(fit_stru3, fit_stru4, regime="diluting", save_path=save_dir / "diluting")
    return [str(save_dir / "diluting.png")]


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
        data_stru = _normalize_conductivity_measurements(loadmat(str(file_path)).get("data_stru"))
        fit_stru, sim_stru, sim_inter = solve_model_B_fix(
            data_stru,
            "DATA",
            theta=fit_stru_base["parameters"],
            sim_opt=False,
            B_form=1,
            sigma_fixed=sigma_fixed,
            workflow_family="DATA2",
        )
        plot_sim_comparison(data_stru, sim_stru, stirc_mass=False, plot_pred=True, lg=False, LOUD=True)
        plt.close("all")
        outputs.extend(
            [
                str(Path("figures") / f"mass-dat{data_stru['dataset']}.png"),
                str(Path("figures") / f"concentration-dat{data_stru['dataset']}.png"),
            ]
        )
    return outputs


def run_data2_model_variations(data_root=None, save_dir=None):
    """Recreate the DATA2 model-variation and FIM calculations."""
    root = _resolve_data2_root(data_root)
    save_dir = Path(save_dir) if save_dir is not None else Path(FIGURES_DIR) / "data2_model_variations"
    save_dir.mkdir(parents=True, exist_ok=True)

    outputs = []

    # Lag - with time correction
    data_stru = _normalize_conductivity_measurements(loadmat(str(root / "data_stru-dataset270511.123.mat")).get("data_stru"))
    mode = "Lag"
    theta = {"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01, "sigma": 1.0, "S0": 0}
    fit_stru_base1, sim_stru, sim_inter = solve_model(
        data_stru, mode, theta, sim_opt=False, B_form=1, LOUD=False, workflow_family="DATA2"
    )
    fit_path = save_dir / "LagM1-fit.json"
    store_json(str(fit_path), fit_stru_base1)
    outputs.append(str(fit_path))

    doe_stru = calc_FIM(data_stru, mode, theta=fit_stru_base1["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
    fim_path = save_dir / "LagM1-FIM.json"
    store_json(str(fim_path), doe_stru)
    outputs.append(str(fim_path))

    # Lag - without time correction
    data_stru = _normalize_conductivity_measurements(loadmat(str(root / "data_stru-dataset270511.122.mat")).get("data_stru"))
    fit_stru, sim_stru, sim_inter = solve_model(
        data_stru, mode, theta=fit_stru_base1["parameters"], sim_opt=False, B_form=1, workflow_family="DATA2"
    )
    fit_path = save_dir / "LagM2-fit.json"
    store_json(str(fit_path), fit_stru)
    outputs.append(str(fit_path))

    # Lag truncated (DATA) - with time correction
    data_stru = _normalize_conductivity_measurements(loadmat(str(root / "data_stru-dataset270511.121.mat")).get("data_stru"))
    mode_truc = "DATA"
    fit_stru, sim_stru, sim_inter = solve_model(
        data_stru, mode_truc, theta=fit_stru_base1["parameters"], sim_opt=False, B_form=1, workflow_family="DATA2"
    )
    fit_path = save_dir / "LagM3-fit.json"
    store_json(str(fit_path), fit_stru)
    outputs.append(str(fit_path))

    doe_stru = calc_FIM(data_stru, mode_truc, theta=fit_stru_base1["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
    fim_path = save_dir / "LagM3-FIM.json"
    store_json(str(fim_path), doe_stru)
    outputs.append(str(fim_path))

    # Lag truncated (DATA) - without time correction
    data_stru = _normalize_conductivity_measurements(loadmat(str(root / "data_stru-dataset270511.12.mat")).get("data_stru"))
    fit_stru, sim_stru, sim_inter = solve_model(
        data_stru, mode_truc, theta=fit_stru_base1["parameters"], sim_opt=False, B_form=1, workflow_family="DATA2"
    )
    fit_path = save_dir / "LagM4-fit.json"
    store_json(str(fit_path), fit_stru)
    outputs.append(str(fit_path))

    # Overflow - with time correction
    data_stru = _normalize_conductivity_measurements(loadmat(str(root / "data_stru-dataset270511.423.mat")).get("data_stru"))
    mode = "Overflow"
    theta = {"Lp": 11, "beta_c": 15, "beta_0": 1, "beta_1": 0.01, "sigma": 1.0, "S0": -0.1756665334051245}
    fit_stru_base2, sim_stru, sim_inter = solve_model(
        data_stru, mode, theta, sim_opt=False, B_form=1, workflow_family="DATA2"
    )
    fit_path = save_dir / "OverflowM1-fit.json"
    store_json(str(fit_path), fit_stru_base2)
    outputs.append(str(fit_path))

    doe_stru = calc_FIM(data_stru, mode, theta=fit_stru_base2["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
    fim_path = save_dir / "OverflowM1-FIM.json"
    store_json(str(fim_path), doe_stru)
    outputs.append(str(fim_path))

    # Overflow - without time correction
    data_stru = _normalize_conductivity_measurements(loadmat(str(root / "data_stru-dataset270511.422.mat")).get("data_stru"))
    fit_stru, sim_stru, sim_inter = solve_model(
        data_stru, mode, theta=fit_stru_base2["parameters"], sim_opt=False, B_form=1, workflow_family="DATA2"
    )
    fit_path = save_dir / "OverflowM2-fit.json"
    store_json(str(fit_path), fit_stru)
    outputs.append(str(fit_path))

    # Overflow truncated (DATA) - with time correction
    data_stru = _normalize_conductivity_measurements(loadmat(str(root / "data_stru-dataset270511.421.mat")).get("data_stru"))
    mode_truc = "DATA"
    fit_stru, sim_stru, sim_inter = solve_model(
        data_stru, mode_truc, theta=fit_stru_base2["parameters"], sim_opt=False, B_form=1, workflow_family="DATA2"
    )
    fit_path = save_dir / "OverflowM3-fit.json"
    store_json(str(fit_path), fit_stru)
    outputs.append(str(fit_path))

    doe_stru = calc_FIM(data_stru, mode_truc, theta=fit_stru_base2["parameters"], step=1e-8, formula="Backward", B_form=1, workflow_family="DATA2")
    fim_path = save_dir / "OverflowM3-FIM.json"
    store_json(str(fim_path), doe_stru)
    outputs.append(str(fim_path))

    # Overflow truncated (DATA) - without time correction
    data_stru = _normalize_conductivity_measurements(loadmat(str(root / "data_stru-dataset270511.42.mat")).get("data_stru"))
    fit_stru, sim_stru, sim_inter = solve_model(
        data_stru, mode_truc, theta=theta, sim_opt=False, B_form=1, workflow_family="DATA2"
    )
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
        fit_stru, sim_stru, sim_inter = solve_model(
            data_stru, mode, sim_opt=False, B_form="pervial", workflow_family="DATA2"
        )
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


def run_data2_paper_reproduction(data_root=None, save_dir=None):
    """Run the DATA2 paper-style helpers in one place."""
    outputs = []
    outputs.extend(run_data2_model_error_visualization(data_root=data_root, save_dir=save_dir))
    outputs.extend(run_cross_verification(data_root=data_root, save_dir=save_dir))
    outputs.extend(run_data2_model_variations(data_root=data_root, save_dir=save_dir))
    outputs.extend(run_pre_B_dependence(data_root=data_root, save_dir=save_dir))
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
            outputs.append(run_sigma_sensitivity(data_root=root, dataset=501.1, save_dir=save_dir))
            outputs.append(run_sigma_sensitivity(data_root=root, dataset=511.12, save_dir=save_dir))
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
