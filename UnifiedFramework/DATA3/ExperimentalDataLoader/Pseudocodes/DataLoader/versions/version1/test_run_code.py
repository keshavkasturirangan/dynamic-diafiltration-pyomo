#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 11:30:00 2025

@authors: Keshav Kasturi Rangan, Alex Dowling
"""

from experiment_read_xlsx_functions import read_excel_workbook
from experiment_read_mat_functions import read_mat_experiment

# Excel workbook: each sheet = experiment
exps = read_excel_workbook("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/ExperimentalDataFiles/NF270_MC3.xlsx")
print("Excel experiments:", len(exps), "first sheet:", exps[0].run_id, "vials:", exps[0].n_vials)

# MAT file: one dataset
exp_m = read_mat_experiment("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset270511.121.mat")
print("MAT run:", exp_m.run_id, "vials:", exp_m.n_vials)
