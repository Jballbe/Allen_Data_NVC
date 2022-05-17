#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:47:59 2022

@author: julienballbe
"""
#%% Import
import allensdk
import json
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.cell_types_api import CellTypesApi
from allensdk.internal.model.biophysical import ephys_utils
import webbrowser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ast import literal_eval
import time
from plotnine import ggplot, geom_line, aes, geom_abline, geom_point, geom_text, labels,geom_histogram
import scipy
from scipy.stats import linregress
from scipy import optimize
from scipy.optimize import curve_fit
import random
import plotly.express as plotly

from allensdk.core.swc import Marker
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
import warnings
from allensdk.core.cell_types_cache import CellTypesCache
import pandas
from allensdk.api.queries.cell_types_api import CellTypesApi
import os
from lmfit.models import LinearModel, StepModel, ExpressionModel, Model,ExponentialModel,ConstantModel
from lmfit import Parameters, Minimizer,fit_report
from plotnine.scales import scale_y_continuous,ylim,xlim,scale_color_manual
from plotnine.labels import xlab
from sklearn.metrics import mean_squared_error
import seaborn as sns

ctc= CellTypesCache(manifest_file="/Users/julienballbe/My_Work/Allen_Data/Common_Script/Full_analysis_cell_types/manifest.json")
full_analysis_cells=get_cells("/Users/julienballbe/My_Work/Allen_Data/Common_Script/Full_analysis_cell_types/manifest","Mouse")
mouse_dict,human_dict=dict_specimen("/Users/julienballbe/My_Work/Allen_Data/Common_Script/Full_analysis_cell_types/manifest")
mouse_sweep_stim_table=create_species_sweeps_stim_table(mouse_dict)
mouse_id_list=mouse_sweep_stim_table['specimen_id']
mouse_id_list=pd.Series(mouse_id_list)

#%%Randomization
id_list=random_id_list(mouse_id_list,40)
test_id=id_list[7]

#%%
mycol=['Cell_id','Starting_frequency_A','Adapt_cst_B','Steady_state_frequency_C']
units=['--',
        'Hz_rel',
        'Spike_index',
        'Hz_rel']
new_fit_comparison_table=pd.DataFrame(columns=mycol)
first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])
mylist=mouse_id_list
number_done=0
total_time=0
for cell_id in mylist:
    start_time=time.time()
    A,B,C,plot=new_fit_exponential_decay(cell_id,mouse_sweep_stim_table,per_time=True,first_x_ms=1000,do_plot=False)
    new_line=pd.Series([str(cell_id),A,B,C],
                   index=['cell_id','Starting_frequency_A','Adapt_cst_B', 'Steady_state_frequency_C'])
    new_fit_comparison_table=new_fit_comparison_table.append(new_line,ignore_index=True)
    end_time=time.time()
    number_done+=1
    number_remaining=len(mylist)-number_done
    total_time+=end_time-start_time
    remaining_time=number_remaining*total_time/number_done
    total_remaining_minutes=remaining_time//60
    remaining_hours=total_remaining_minutes//60
    remaining_minutes=total_remaining_minutes%60
    remaining_seconds=remaining_time%60
    print('Remaining time: ',remaining_hours,'h ',remaining_minutes,'min ',round(remaining_seconds,0),"s")

adaptation_table=pd.concat([first_two_lines,new_fit_comparison_table])
adaptation_table=adaptation_table.iloc[1:,:]
adaptation_table.to_csv(path_or_buf="/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_04_14/Adaptation_computation_1000ms.csv",na_rep="nan",index=False)
print("Total time= ",total_time)


    
    
#%%5ms

mycol=['Cell_id',"fit","best_single_QNRMSE","best_single_amplitude","best_single_center","best_single_sigma","best_compo_QNRMSE","best_step_amplitude","best_step_center","best_step_sigma","best_sigmoid_amplitude","best_sigmoid_center","best_sigmoid_sigma"]
units=['--',
       '--',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       'Hz',
       'pA',
       'Hz/pA']
        


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])
start_time=time.time()
total_time=0
number_done=0
mylist=mouse_id_list

new_fit_table=pd.DataFrame(columns=mycol)
for cell_id in mylist:
    start_time=time.time()
    new_line=pd.Series([str(cell_id),*new_fit_sigmoid(cell_id,mouse_sweep_stim_table,per_time=True,first_x_ms=5,do_plot=False)],
                   index=mycol)
    new_fit_table=new_fit_table.append(new_line,ignore_index=True)
    end_time=time.time()
    number_done+=1
    number_remaining=len(mylist)-number_done
    total_time+=end_time-start_time
    remaining_time=number_remaining*total_time/number_done
    total_remaining_minutes=remaining_time//60
    remaining_hours=total_remaining_minutes//60
    remaining_minutes=total_remaining_minutes%60
    remaining_seconds=remaining_time%60
    print('Remaining time: ',remaining_hours,'h ',remaining_minutes,'min ',round(remaining_seconds,0),"s")

print("Total time new_fit_sigmoid= ",total_time)


fit_table=pd.concat([first_two_lines,new_fit_table])
fit_table=fit_table.iloc[1:,:]
fit_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/fit_table_5ms.csv',na_rep='nan',index=False)
    
units=['--','--','Hz/pA','pA','Hz']
units=pd.Series(units,index=["Cell_id","Fit","Gain","Threshold","Saturation"])
main_start_time=time.time()

f_I_table=compute_f_I_params(new_fit_table,mouse_sweep_stim_table,per_time=True,first_x_ms=5,per_nth_spike=False,first_nth_spike=0,do_plot=False)
f_I_table=pd.concat([units,f_I_table])
f_I_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/f_I_table_5ms.csv',na_rep='nan',index=False)
end_time=time.time()

print("time for f_I_table=",end_time-start_time,"s")


#%% 10ms

mycol=['Cell_id',"fit","best_single_QNRMSE","best_single_amplitude","best_single_center","best_single_sigma","best_compo_QNRMSE","best_step_amplitude","best_step_center","best_step_sigma","best_sigmoid_amplitude","best_sigmoid_center","best_sigmoid_sigma"]
units=['--',
       '--',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       'Hz',
       'pA',
       'Hz/pA']
        


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])
start_time=time.time()
total_time=0
number_done=0
mylist=mouse_id_list

new_fit_table=pd.DataFrame(columns=mycol)
for cell_id in mylist:
    start_time=time.time()
    new_line=pd.Series([str(cell_id),*new_fit_sigmoid(cell_id,mouse_sweep_stim_table,per_time=True,first_x_ms=10,do_plot=False)],
                   index=mycol)
    new_fit_table=new_fit_table.append(new_line,ignore_index=True)
    end_time=time.time()
    number_done+=1
    number_remaining=len(mylist)-number_done
    total_time+=end_time-start_time
    remaining_time=number_remaining*total_time/number_done
    total_remaining_minutes=remaining_time//60
    remaining_hours=total_remaining_minutes//60
    remaining_minutes=total_remaining_minutes%60
    remaining_seconds=remaining_time%60
    print('Remaining time: ',remaining_hours,'h ',remaining_minutes,'min ',round(remaining_seconds,0),"s")

print("Total time new_fit_sigmoid= ",total_time)


fit_table=pd.concat([first_two_lines,new_fit_table])
fit_table=fit_table.iloc[1:,:]
fit_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/fit_table_10ms.csv',na_rep='nan',index=False)
    
units=['--','--','Hz/pA','pA','Hz']
units=pd.Series(units,index=["Cell_id","Fit","Gain","Threshold","Saturation"])
main_start_time=time.time()

f_I_table=compute_f_I_params(new_fit_table,mouse_sweep_stim_table,per_time=True,first_x_ms=10,per_nth_spike=False,first_nth_spike=0,do_plot=False)
f_I_table=pd.concat([units,f_I_table])
f_I_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/f_I_table_10ms.csv',na_rep='nan',index=False)
end_time=time.time()

print("time for f_I_table=",end_time-start_time,"s")
#%% 25ms

mycol=['Cell_id',"fit","best_single_QNRMSE","best_single_amplitude","best_single_center","best_single_sigma","best_compo_QNRMSE","best_step_amplitude","best_step_center","best_step_sigma","best_sigmoid_amplitude","best_sigmoid_center","best_sigmoid_sigma"]
units=['--',
       '--',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       'Hz',
       'pA',
       'Hz/pA']
        


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])
start_time=time.time()
total_time=0
number_done=0
mylist=mouse_id_list

new_fit_table=pd.DataFrame(columns=mycol)
for cell_id in mylist:
    start_time=time.time()
    new_line=pd.Series([str(cell_id),*new_fit_sigmoid(cell_id,mouse_sweep_stim_table,per_time=True,first_x_ms=25,do_plot=False)],
                   index=mycol)
    new_fit_table=new_fit_table.append(new_line,ignore_index=True)
    end_time=time.time()
    number_done+=1
    number_remaining=len(mylist)-number_done
    total_time+=end_time-start_time
    remaining_time=number_remaining*total_time/number_done
    total_remaining_minutes=remaining_time//60
    remaining_hours=total_remaining_minutes//60
    remaining_minutes=total_remaining_minutes%60
    remaining_seconds=remaining_time%60
    print('Remaining time: ',remaining_hours,'h ',remaining_minutes,'min ',round(remaining_seconds,0),"s")

print("Total time new_fit_sigmoid= ",total_time)


fit_table=pd.concat([first_two_lines,new_fit_table])
fit_table=fit_table.iloc[1:,:]
fit_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/fit_table_25ms.csv',na_rep='nan',index=False)
    
units=['--','--','Hz/pA','pA','Hz']
units=pd.Series(units,index=["Cell_id","Fit","Gain","Threshold","Saturation"])
main_start_time=time.time()

f_I_table=compute_f_I_params(new_fit_table,mouse_sweep_stim_table,per_time=True,first_x_ms=25,per_nth_spike=False,first_nth_spike=0,do_plot=False)
f_I_table=pd.concat([units,f_I_table])
f_I_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/f_I_table_25ms.csv',na_rep='nan',index=False)
end_time=time.time()

print("time for f_I_table=",end_time-start_time,"s")
#%% 50ms


mycol=['Cell_id',"fit","best_single_QNRMSE","best_single_amplitude","best_single_center","best_single_sigma","best_compo_QNRMSE","best_step_amplitude","best_step_center","best_step_sigma","best_sigmoid_amplitude","best_sigmoid_center","best_sigmoid_sigma"]
units=['--',
       '--',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       'Hz',
       'pA',
       'Hz/pA']
        


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])
start_time=time.time()
total_time=0
number_done=0
mylist=mouse_id_list

new_fit_table=pd.DataFrame(columns=mycol)
for cell_id in mylist:
    start_time=time.time()
    new_line=pd.Series([str(cell_id),*new_fit_sigmoid(cell_id,mouse_sweep_stim_table,per_time=True,first_x_ms=50,do_plot=False)],
                   index=mycol)
    new_fit_table=new_fit_table.append(new_line,ignore_index=True)
    end_time=time.time()
    number_done+=1
    number_remaining=len(mylist)-number_done
    total_time+=end_time-start_time
    remaining_time=number_remaining*total_time/number_done
    total_remaining_minutes=remaining_time//60
    remaining_hours=total_remaining_minutes//60
    remaining_minutes=total_remaining_minutes%60
    remaining_seconds=remaining_time%60
    print('Remaining time: ',remaining_hours,'h ',remaining_minutes,'min ',round(remaining_seconds,0),"s")

print("Total time new_fit_sigmoid= ",total_time)


fit_table=pd.concat([first_two_lines,new_fit_table])
fit_table=fit_table.iloc[1:,:]
fit_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/fit_table_50ms.csv',na_rep='nan',index=False)
    
units=['--','--','Hz/pA','pA','Hz']
units=pd.Series(units,index=["Cell_id","Fit","Gain","Threshold","Saturation"])
main_start_time=time.time()

f_I_table=compute_f_I_params(new_fit_table,mouse_sweep_stim_table,per_time=True,first_x_ms=50,per_nth_spike=False,first_nth_spike=0,do_plot=False)
f_I_table=pd.concat([units,f_I_table])
f_I_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/f_I_table_50ms.csv',na_rep='nan',index=False)
end_time=time.time()

print("time for f_I_table=",end_time-start_time,"s")
#%% 100ms


mycol=['Cell_id',"fit","best_single_QNRMSE","best_single_amplitude","best_single_center","best_single_sigma","best_compo_QNRMSE","best_step_amplitude","best_step_center","best_step_sigma","best_sigmoid_amplitude","best_sigmoid_center","best_sigmoid_sigma"]
units=['--',
       '--',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       'Hz',
       'pA',
       'Hz/pA']
        


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])
start_time=time.time()
total_time=0
number_done=0
mylist=mouse_id_list

new_fit_table=pd.DataFrame(columns=mycol)
for cell_id in mylist:
    start_time=time.time()
    new_line=pd.Series([str(cell_id),*new_fit_sigmoid(cell_id,mouse_sweep_stim_table,per_time=True,first_x_ms=100,do_plot=False)],
                   index=mycol)
    new_fit_table=new_fit_table.append(new_line,ignore_index=True)
    end_time=time.time()
    number_done+=1
    number_remaining=len(mylist)-number_done
    total_time+=end_time-start_time
    remaining_time=number_remaining*total_time/number_done
    total_remaining_minutes=remaining_time//60
    remaining_hours=total_remaining_minutes//60
    remaining_minutes=total_remaining_minutes%60
    remaining_seconds=remaining_time%60
    print('Remaining time: ',remaining_hours,'h ',remaining_minutes,'min ',round(remaining_seconds,0),"s")

print("Total time new_fit_sigmoid= ",total_time)


fit_table=pd.concat([first_two_lines,new_fit_table])
fit_table=fit_table.iloc[1:,:]
fit_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/fit_table_100ms.csv',na_rep='nan',index=False)
    
units=['--','--','Hz/pA','pA','Hz']
units=pd.Series(units,index=["Cell_id","Fit","Gain","Threshold","Saturation"])
main_start_time=time.time()

f_I_table=compute_f_I_params(new_fit_table,mouse_sweep_stim_table,per_time=True,first_x_ms=100,per_nth_spike=False,first_nth_spike=0,do_plot=False)
f_I_table=pd.concat([units,f_I_table])
f_I_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/f_I_table_100ms.csv',na_rep='nan',index=False)
end_time=time.time()

print("time for f_I_table=",end_time-start_time,"s")
#%% 250ms

mycol=['Cell_id',"fit","best_single_QNRMSE","best_single_amplitude","best_single_center","best_single_sigma","best_compo_QNRMSE","best_step_amplitude","best_step_center","best_step_sigma","best_sigmoid_amplitude","best_sigmoid_center","best_sigmoid_sigma"]
units=['--',
       '--',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       'Hz',
       'pA',
       'Hz/pA']
        


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])
start_time=time.time()
total_time=0
number_done=0
mylist=mouse_id_list

new_fit_table=pd.DataFrame(columns=mycol)
for cell_id in mylist:
    start_time=time.time()
    new_line=pd.Series([str(cell_id),*new_fit_sigmoid(cell_id,mouse_sweep_stim_table,per_time=True,first_x_ms=250,do_plot=False)],
                   index=mycol)
    new_fit_table=new_fit_table.append(new_line,ignore_index=True)
    end_time=time.time()
    number_done+=1
    number_remaining=len(mylist)-number_done
    total_time+=end_time-start_time
    remaining_time=number_remaining*total_time/number_done
    total_remaining_minutes=remaining_time//60
    remaining_hours=total_remaining_minutes//60
    remaining_minutes=total_remaining_minutes%60
    remaining_seconds=remaining_time%60
    print('Remaining time: ',remaining_hours,'h ',remaining_minutes,'min ',round(remaining_seconds,0),"s")

print("Total time new_fit_sigmoid= ",total_time)


fit_table=pd.concat([first_two_lines,new_fit_table])
fit_table=fit_table.iloc[1:,:]
fit_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/fit_table_250ms.csv',na_rep='nan',index=False)
    
units=['--','--','Hz/pA','pA','Hz']
units=pd.Series(units,index=["Cell_id","Fit","Gain","Threshold","Saturation"])
main_start_time=time.time()

f_I_table=compute_f_I_params(new_fit_table,mouse_sweep_stim_table,per_time=True,first_x_ms=250,per_nth_spike=False,first_nth_spike=0,do_plot=False)
f_I_table=pd.concat([units,f_I_table])
f_I_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/f_I_table_250ms.csv',na_rep='nan',index=False)
end_time=time.time()

print("time for f_I_table=",end_time-start_time,"s")
#%% 500ms


mycol=['Cell_id',"fit","best_single_QNRMSE","best_single_amplitude","best_single_center","best_single_sigma","best_compo_QNRMSE","best_step_amplitude","best_step_center","best_step_sigma","best_sigmoid_amplitude","best_sigmoid_center","best_sigmoid_sigma"]
units=['--',
       '--',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       'Hz',
       'pA',
       'Hz/pA']
        


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])
start_time=time.time()
total_time=0
number_done=0
mylist=mouse_id_list

new_fit_table=pd.DataFrame(columns=mycol)
for cell_id in mylist:
    start_time=time.time()
    new_line=pd.Series([str(cell_id),*new_fit_sigmoid(cell_id,mouse_sweep_stim_table,per_time=True,first_x_ms=500,do_plot=False)],
                   index=mycol)
    new_fit_table=new_fit_table.append(new_line,ignore_index=True)
    end_time=time.time()
    number_done+=1
    number_remaining=len(mylist)-number_done
    total_time+=end_time-start_time
    remaining_time=number_remaining*total_time/number_done
    total_remaining_minutes=remaining_time//60
    remaining_hours=total_remaining_minutes//60
    remaining_minutes=total_remaining_minutes%60
    remaining_seconds=remaining_time%60
    print('Remaining time: ',remaining_hours,'h ',remaining_minutes,'min ',round(remaining_seconds,0),"s")

print("Total time new_fit_sigmoid= ",total_time)


fit_table=pd.concat([first_two_lines,new_fit_table])
fit_table=fit_table.iloc[1:,:]
fit_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/fit_table_500ms.csv',na_rep='nan',index=False)
    
units=['--','--','Hz/pA','pA','Hz']
units=pd.Series(units,index=["Cell_id","Fit","Gain","Threshold","Saturation"])
main_start_time=time.time()

f_I_table=compute_f_I_params(new_fit_table,mouse_sweep_stim_table,per_time=True,first_x_ms=500,per_nth_spike=False,first_nth_spike=0,do_plot=False)
f_I_table=pd.concat([units,f_I_table])
f_I_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/f_I_table_500ms.csv',na_rep='nan',index=False)
end_time=time.time()

print("time for f_I_table=",end_time-start_time,"s")
#%% 1000ms


mycol=['Cell_id',"fit","best_single_QNRMSE","best_single_amplitude","best_single_center","best_single_sigma","best_compo_QNRMSE","best_step_amplitude","best_step_center","best_step_sigma","best_sigmoid_amplitude","best_sigmoid_center","best_sigmoid_sigma"]
units=['--',
       '--',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       '--',
       'Hz',
       'pA',
       'Hz/pA',
       'Hz',
       'pA',
       'Hz/pA']
        


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])
start_time=time.time()
total_time=0
number_done=0
mylist=mouse_id_list

new_fit_table=pd.DataFrame(columns=mycol)
for cell_id in mylist:
    start_time=time.time()
    new_line=pd.Series([str(cell_id),*new_fit_sigmoid(cell_id,mouse_sweep_stim_table,per_time=True,first_x_ms=1000,do_plot=False)],
                   index=mycol)
    new_fit_table=new_fit_table.append(new_line,ignore_index=True)
    end_time=time.time()
    number_done+=1
    number_remaining=len(mylist)-number_done
    total_time+=end_time-start_time
    remaining_time=number_remaining*total_time/number_done
    total_remaining_minutes=remaining_time//60
    remaining_hours=total_remaining_minutes//60
    remaining_minutes=total_remaining_minutes%60
    remaining_seconds=remaining_time%60
    print('Remaining time: ',remaining_hours,'h ',remaining_minutes,'min ',round(remaining_seconds,0),"s")

print("Total time new_fit_sigmoid= ",total_time)


fit_table=pd.concat([first_two_lines,new_fit_table])
fit_table=fit_table.iloc[1:,:]
fit_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/fit_table_1000ms.csv',na_rep='nan',index=False)
    
units=['--','--','Hz/pA','pA','Hz']
units=pd.Series(units,index=["Cell_id","Fit","Gain","Threshold","Saturation"])
main_start_time=time.time()

f_I_table=compute_f_I_params(new_fit_table,mouse_sweep_stim_table,per_time=True,first_x_ms=1000,per_nth_spike=False,first_nth_spike=0,do_plot=False)
f_I_table=pd.concat([units,f_I_table])
f_I_table.to_csv(path_or_buf='/Users/julienballbe/My_Work/Allen_Data/Feature_computation/2022_05_04/f_I_table_1000ms.csv',na_rep='nan',index=False)
end_time=time.time()

print("time for f_I_table=",end_time-start_time,"s")
#%% 4spikes

mycol=['Cell_id','Gain','Threshold','Saturation','Neuron_type','NRMSE_sigmoid','NRMSE_composite','Parameter_amplitude','Parameter_center','Parameter_sigma','Starting_frequency_A','Adapt_cst_B','Steady_state_frequency_C','Normalized_starting_freq_Anorm','Normalized_ss_freq_Cnorm','Input_resistance_MOhm']
units=['--',
        'Hz/pA',
        'pA',
        'Hz',
        '--',
        'WU',
        'WU',
        'Hz',
        'pA',
        'Hz/pA',
        'WU',
        'spike_index',
        'WU',
        'WU',
        'WU',
        'MOhm']


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])





start_time=time.time()
features_4spikes=compute_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=4)
features_4spikes=features_4spikes.reindex(columns=mycol)
File_4spikes=pd.concat([first_two_lines,features_4spikes])
File_4spikes=File_4spikes.iloc[1:,:]
File_4spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/features_4spikes.csv"),na_rep="nan",index=False)
end_4spike=time.time()
print("time for 4spikes=",end_4spike-start_time,"s")

#%% 5spikes
mycol=['Cell_id','Gain','Threshold','Saturation','Neuron_type','NRMSE_sigmoid','NRMSE_composite','Parameter_amplitude','Parameter_center','Parameter_sigma','Starting_frequency_A','Adapt_cst_B','Steady_state_frequency_C','Normalized_starting_freq_Anorm','Normalized_ss_freq_Cnorm','Input_resistance_MOhm']
units=['--',
        'Hz/pA',
        'pA',
        'Hz',
        '--',
        'WU',
        'WU',
        'Hz',
        'pA',
        'Hz/pA',
        'WU',
        'spike_index',
        'WU',
        'WU',
        'WU',
        'MOhm']


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])





start_time=time.time()
features_5spikes=compute_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=5)
features_5spikes=features_5spikes.reindex(columns=mycol)
File_5spikes=pd.concat([first_two_lines,features_5spikes])
File_5spikes=File_5spikes.iloc[1:,:]
File_5spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/features_5spikes.csv"),na_rep="nan",index=False)
end_5spike=time.time()
print("time for 5spikes=",end_5spike-start_time,"s")

#%% 6spikes

mycol=['Cell_id','Gain','Threshold','Saturation','Neuron_type','NRMSE_sigmoid','NRMSE_composite','Parameter_amplitude','Parameter_center','Parameter_sigma','Starting_frequency_A','Adapt_cst_B','Steady_state_frequency_C','Normalized_starting_freq_Anorm','Normalized_ss_freq_Cnorm','Input_resistance_MOhm']
units=['--',
        'Hz/pA',
        'pA',
        'Hz',
        '--',
        'WU',
        'WU',
        'Hz',
        'pA',
        'Hz/pA',
        'WU',
        'spike_index',
        'WU',
        'WU',
        'WU',
        'MOhm']


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])




start_time=time.time()
features_6spikes=compute_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=6)
features_6spikes=features_6spikes.reindex(columns=mycol)
File_6spikes=pd.concat([first_two_lines,features_6spikes])
File_6spikes=File_6spikes.iloc[1:,:]
File_6spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/features_6spikes.csv"),na_rep="nan",index=False)
end_6spike=time.time()
print("time for 6spikes=",end_6spike-start_time,"s")

#%% 7spikes
mycol=['Cell_id','Gain','Threshold','Saturation','Neuron_type','NRMSE_sigmoid','NRMSE_composite','Parameter_amplitude','Parameter_center','Parameter_sigma','Starting_frequency_A','Adapt_cst_B','Steady_state_frequency_C','Normalized_starting_freq_Anorm','Normalized_ss_freq_Cnorm','Input_resistance_MOhm']
units=['--',
        'Hz/pA',
        'pA',
        'Hz',
        '--',
        'WU',
        'WU',
        'Hz',
        'pA',
        'Hz/pA',
        'WU',
        'spike_index',
        'WU',
        'WU',
        'WU',
        'MOhm']


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])





start_time=time.time()
features_7spikes=compute_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=7)
features_7spikes=features_7spikes.reindex(columns=mycol)
File_7spikes=pd.concat([first_two_lines,features_7spikes])
File_7spikes=File_7spikes.iloc[1:,:]
File_7spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/features_7spikes.csv"),na_rep="nan",index=False)
end_7spike=time.time()
print("time for 7spikes=",end_7spike-start_time,"s")

#%% 8spikes

mycol=['Cell_id','Gain','Threshold','Saturation','Neuron_type','NRMSE_sigmoid','NRMSE_composite','Parameter_amplitude','Parameter_center','Parameter_sigma','Starting_frequency_A','Adapt_cst_B','Steady_state_frequency_C','Normalized_starting_freq_Anorm','Normalized_ss_freq_Cnorm','Input_resistance_MOhm']
units=['--',
        'Hz/pA',
        'pA',
        'Hz',
        '--',
        'WU',
        'WU',
        'Hz',
        'pA',
        'Hz/pA',
        'WU',
        'spike_index',
        'WU',
        'WU',
        'WU',
        'MOhm']


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])



start_time=time.time()
features_8spikes=compute_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=8)
features_8spikes=features_8spikes.reindex(columns=mycol)
File_8spikes=pd.concat([first_two_lines,features_8spikes])
File_8spikes=File_8spikes.iloc[1:,:]
File_8spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/features_8spikes.csv"),na_rep="nan",index=False)
end_8spike=time.time()
print("time for 8spikes=",end_8spike-start_time,"s")

#%% 9spikes
mycol=['Cell_id','Gain','Threshold','Saturation','Neuron_type','NRMSE_sigmoid','NRMSE_composite','Parameter_amplitude','Parameter_center','Parameter_sigma','Starting_frequency_A','Adapt_cst_B','Steady_state_frequency_C','Normalized_starting_freq_Anorm','Normalized_ss_freq_Cnorm','Input_resistance_MOhm']
units=['--',
        'Hz/pA',
        'pA',
        'Hz',
        '--',
        'WU',
        'WU',
        'Hz',
        'pA',
        'Hz/pA',
        'WU',
        'spike_index',
        'WU',
        'WU',
        'WU',
        'MOhm']


first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])




start_time=time.time()
features_9spikes=compute_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=9)
features_9spikes=features_9spikes.reindex(columns=mycol)
File_9spikes=pd.concat([first_two_lines,features_9spikes])
File_9spikes=File_9spikes.iloc[1:,:]
File_9spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/features_9spikes.csv"),na_rep="nan",index=False)
end_9spike=time.time()
print("time for 9spikes=",end_9spike-start_time,"s")
#%% 10spikes
mycol=['Cell_id','Gain','Threshold','Saturation','Neuron_type','NRMSE_sigmoid','NRMSE_composite','Parameter_amplitude','Parameter_center','Parameter_sigma','Starting_frequency_A','Adapt_cst_B','Steady_state_frequency_C','Normalized_starting_freq_Anorm','Normalized_ss_freq_Cnorm','Input_resistance_MOhm']
units=['--',
        'Hz/pA',
        'pA',
        'Hz',
        '--',
        'WU',
        'WU',
        'Hz',
        'pA',
        'Hz/pA',
        'WU',
        'spike_index',
        'WU',
        'WU',
        'WU',
        'MOhm']

first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])






start_time=time.time()
features_10spikes=compute_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=10)
features_10spikes=features_10spikes.reindex(columns=mycol)
File_10spikes=pd.concat([first_two_lines,features_10spikes])
File_10spikes=File_10spikes.iloc[1:,:]
File_10spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/features_10spikes.csv"),na_rep="nan",index=False)
end_10spike=time.time()
print("time for 10spikes=",end_10spike-start_time,"s")
start_time=time.time()

#%%Functions
def fit_specimen_fi_slope(stim_amps, avg_rates):
    """
    Fit the rate and stimulus amplitude to a line and return the slope of the fit.

    Parameters
    ----------
    stim_amps: array of sweeps amplitude in mA
    avg_rates: array of sweeps avergae firing rate in Hz
    Returns
    -------
    m: f-I curve slope for the specimen
    c:f-I curve intercept for the specimen

    """

    x = stim_amps
    y = avg_rates

    A = np.vstack([x, np.ones_like(x)]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]

    return m, c
def get_cells(file_name, species, reconstruction=False, morphology=False):
    if species == "Mouse":
        my_species = [CellTypesApi.MOUSE]
        ctc.get_cells(file_name=file_name, require_reconstruction=reconstruction, require_morphology=morphology,
                      species=my_species)
    elif species == "Human":
        my_species = [CellTypesApi.HUMAN]
        ctc.get_cells(file_name=file_name, require_reconstruction=reconstruction, require_morphology=morphology,
                      species=my_species)
    elif species == "All":
        ctc.get_cells(file_name=file_name, require_reconstruction=reconstruction, require_morphology=morphology)

    return ("Data saved in file: " + str(file_name))

def stimulus_type(cell_id):
    '''
    Returns for the different stimuli the sweeps using those stimuli
    Parameters
    ----------
    cell_id : str
        cell_id found in the id_list from take_id function.
    Returns
    -------
    dict_stimuli : dictionary
        dictionary with the keys being the different stimulus types and the values are the lists of sweep numbers.
    '''
    all_sweeps = ctc.get_ephys_sweeps(cell_id)
    dict_stimuli = dict()
    for i in all_sweeps:
        number = i["sweep_number"]
        if i["stimulus_name"] not in dict_stimuli:
            dict_stimuli[i["stimulus_name"]] = [number]
        else:
            dict_stimuli[i["stimulus_name"]].append(number)
    return (dict_stimuli)

def ephys_web_page(cell_id):
    '''
    The electrophysiological web page of the cell id in the Allen Brain Atlas is opened
    Parameters
    ----------
    cell_id : str
        one str value of the id_list from take_id function
    '''

    link = "http://celltypes.brain-map.org/experiment/electrophysiology/" + str(cell_id)
    webbrowser.open(link)

def dict_specimen(name_file):
    '''
    Generate lists of dictionary for all specimen per species
    Parameters
    ----------
    name_file : str
        File name entered in get_cells function as string.
    Returns
    -------
    d_mouse : list
        List of Dictionnary for all mouse specimen
    d_human : list
        List of dictionnary for all human specimen.
    '''

    f = open(name_file)
    data = json.load(f)
    d_mouse = []
    d_human = []
    for i in data:
        if i["donor__species"] == "Mus musculus":
            d_mouse.append(i)
        else:
            d_human.append(i)
    return (d_mouse, d_human)

def get_structure_name(species_dict, speciment_index):
    '''
    Return structure substructure and layer for create_species_table function
    Parameters
    ----------
    species_dict : Dictionnary
        From dict_specimen function.
    speciment_index : int
        Index of the specimen of interest.
    Returns
    -------
    name : str
    substructure : str
    layer : str
    '''

    full_name = species_dict[speciment_index]['structure__name']  # full_name=dict_mouse[i]["structure__name"]
    j = 0
    for e in full_name:
        if e == ",":
            j = j + 1
    if j == 0:
        name = full_name.replace('"', '')
        substructure = "No substructure"
        layer = "No layer"
    elif j == 1:
        name1, layer1 = full_name.split(",")
        name = name1.replace('"', '')
        layer = layer1.replace('"', '')
        substructure = "No substructure"

    elif j == 2:
        name1, substructure, layer1 = full_name.split(",")
        name = name1.replace('"', '')
        layer = layer1.replace('"', '')
    return name, substructure, layer

def create_species_sweeps_stim_table(species_dict, from_csv=False):
    '''
    Either from species dict (default) or from already saved table(from_csv=True)
    Create a dataframe containing for a given species in each row
    the specimen_id,
    the structure acronyme,
    the structure name,
    the substructure name,
    the layer,
    the dendrite_type,
    the total number of sweeps and
    the sweeps_id per stimulus kind
    Parameters
    ----------
    species_dict : Dictionnary or str
        From dict_specimen function
        Or path to csv file
    from_csv : Bool
        Default False: Create table from a csv file
    Returns
    -------
    full_dataframe : DataFrame

    '''
    if from_csv == True:
        full_dataframe = pd.read_csv(str(species_dict))
        full_dataframe = full_dataframe.iloc[:, 1:]

        for row in range(full_dataframe.shape[0]):
            for col in range(7, 20):
                if type(full_dataframe.iloc[row, col]) == str:
                    full_dataframe.iat[row, col] = literal_eval(full_dataframe.iat[row, col])

                else:
                    full_dataframe.iloc[row, col] = -1
        return full_dataframe

    start = time.time()
    full_dataframe = pd.DataFrame()

    for current_specimen in range(len(species_dict)):
        name, substructure, layer = get_structure_name(species_dict, current_specimen)
        number_of_sweeps = 0

        specimen_stimulus_dict = stimulus_type(species_dict[current_specimen]["specimen__id"])

        for key in stimulus_type(species_dict[current_specimen]["specimen__id"]).keys():
            number_of_sweeps += len(stimulus_type(species_dict[current_specimen]["specimen__id"])[key])

        current_row = pd.Series([species_dict[current_specimen]["specimen__id"],
                                 species_dict[current_specimen]["structure__acronym"],
                                 name,
                                 substructure,
                                 layer,
                                 species_dict[current_specimen]["tag__dendrite_type"],
                                 number_of_sweeps],
                                index=['specimen_id', 'structure_acronyme', 'structure_name', 'structure_substructure',
                                       'layer', 'dendrite_type', 'number_of_sweeps'])
        current_row = pd.DataFrame(current_row).T

        stim_dict = pd.Series(specimen_stimulus_dict)
        stim_dict = pd.DataFrame(stim_dict).T

        full_row = pd.concat([current_row, stim_dict], axis=1)

        full_dataframe = pd.concat([full_dataframe, full_row], ignore_index=True, axis=0, join='outer')
    end = time.time()
    print('Running time= ', end - start)
    return full_dataframe


def compute_feature(specimen_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_nth_spike=False,first_nth_spike=0):
    '''
    

    Parameters
    ----------
    specimen_id : list
        list of specimen id.
    species_sweep_stim_table : DataFrame
        DataFrame coming from create_species_sweeps_stim_table function.
    per_time : Boolean, optional
        Extract feature per time . The default is False.
    first_x_ms : int, optional
        If per_time==True , indicate the time to take into account to extract features in ms (x ms after the start of the stimulus). The default is 0.
    per_nth_spike : Boolean, optional
        Extract feature per spike number. The default is False.
    first_nth_spike : int, optional
        If per_nth_spike==True, indicate the n first spike to take into account for feature extraction. The default is 0.

    Returns
    -------
    full_table : DataFrame
        DataFrame containing the different feature values for each cell .

    '''
    mycolumns=['Cell_id','Gain','Threshold','Saturation','Neuron_type','NRMSE_sigmoid','NRMSE_composite','Parameter_amplitude','Parameter_center','Parameter_sigma','Starting_frequency_A','Adapt_cst_B','Steady_state_frequency_C','Normalized_starting_freq_Anorm','Normalized_ss_freq_Cnorm','Input_resistance_MOhm']
    full_table=pd.DataFrame(columns=mycolumns)
    ephys_features_table=ctc.get_ephys_features()
    ephys_features_table=pd.DataFrame(ephys_features_table)
    for current_specimen_id in specimen_id:
        print(current_specimen_id)
        estimated_slope,estimated_threshold,estimated_saturation,Neuron_type,my_plot,NRMSE_sigmoid,NRMSE_composite,parameter_amplitude,parameter_center,parameter_sigma=fit_sigmoid(current_specimen_id,
                                                                                                                                                                                    species_sweep_stim_table,
                                                                                                                                                                                    per_time=per_time,
                                                                                                                                                                                    first_x_ms=first_x_ms,
                                                                                                                                                                                    per_nth_spike=per_nth_spike,
                                                                                                                                                                                    first_nth_spike=first_nth_spike,
                                                                                                                                                                                    do_plot=False)
        A,B,C,my_plot=fit_exponential_decay(current_specimen_id,
                                            species_sweep_stim_table,
                                            per_time=per_time,
                                            first_x_ms=first_x_ms,
                                            per_spike_nb=per_nth_spike,
                                            first_nth_spike=first_nth_spike,
                                            do_plot=False)
       
        A_norm=A/(A+C)
        C_norm=C/(A+C)
        R_in=ephys_features_table.loc[ephys_features_table["specimen_id"]==current_specimen_id]['input_resistance_mohm'].values[0]
        new_line=pd.Series([str(current_specimen_id),
                            round(estimated_slope,3),
                            round(estimated_threshold,3),
                            round(estimated_saturation,3),
                            Neuron_type,
                            round(NRMSE_sigmoid,3),
                            round(NRMSE_composite,3),
                            round(parameter_amplitude,5),
                            round(parameter_center,5),
                            round(parameter_sigma,5),
                            round(A,3),
                            round(B,3),
                            round(C,3),
                            round(A_norm,3),
                            round(C_norm,3),
                            round(R_in,3)],
                           index=mycolumns)

        full_table=full_table.append(new_line, ignore_index=True)

    return full_table

def extract_stim_freq(specimen_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_nth_spike=False,first_nth_spike=0):
    '''
    Function to extract for each specified specimen_id and the corresponding stimulus the frequency of the response
    Frequency is defined as the number of spikes divided by the time between the stimulus start and the time of the specified index
    Parameters
    ----------
    specimen_id : int

    
    Returns
    -------
    f_I_table : DataFrame
        DataFrame with a column "specimen_id"(factor),the sweep number (int),the stimulus amplitude in pA(float),and the computed frequency of the response (float).
    '''
    
    f_I_table = pd.DataFrame(columns=['specimen', 'sweep_number', 'Stim_amp_pA', 'Frequency_Hz'])
    index_stim = species_sweep_stim_table.columns.get_loc('Long Square')
    index_specimen = species_sweep_stim_table.index[species_sweep_stim_table["specimen_id"] == specimen_id][0]
    
    my_specimen_data = ctc.get_ephys_data(specimen_id)
    sweep_numbers = species_sweep_stim_table.iloc[index_specimen, index_stim]
    
    
    for current_sweep in sweep_numbers:
        
        index_range=my_specimen_data.get_sweep(current_sweep)["index_range"]
        
        sampling_rate=my_specimen_data.get_sweep(current_sweep)["sampling_rate"]
        current_stim_array=(my_specimen_data.get_sweep(current_sweep)["stimulus"][0:index_range[1]+1])* 1e12 #to pA
    
        
        stim_start_index=index_range[0]+next(x for x, val in enumerate(current_stim_array[index_range[0]:]) if val != 0 )
        current_time_array=np.arange(0, len(current_stim_array)) * (1.0 / sampling_rate)
        
        stim_start_time=current_time_array[stim_start_index]
       
       
        
        if len(my_specimen_data.get_spike_times(current_sweep)) <1:
            freq = 0

        else :
            if per_nth_spike==True:
                reshaped_spike_times=my_specimen_data.get_spike_times(current_sweep)[:first_nth_spike]
    
                #nb_spike = len(reshaped_spike_times)
                t_last_spike = reshaped_spike_times[-1]
                #freq = nb_spike / (t_last_spike - stim_start_time)
                freq=len(reshaped_spike_times)/((t_last_spike - stim_start_time))

            elif per_time==True:
                end_time=stim_start_time+(first_x_ms*1e-3)

                reshaped_spike_times=my_specimen_data.get_spike_times(current_sweep)[my_specimen_data.get_spike_times(current_sweep) <= end_time ]
                nb_spike = len(reshaped_spike_times)

                if nb_spike !=0:
                    
                    freq=nb_spike/(first_x_ms*1e-3)
                else:
                    freq=0
        new_line = pd.Series([int(specimen_id), current_sweep,
                              my_specimen_data.get_sweep_metadata(current_sweep)['aibs_stimulus_amplitude_pa'],
                              freq],
                             index=['specimen', 'sweep_number', 'Stim_amp_pA', 'Frequency_Hz'])
        f_I_table = f_I_table.append(new_line, ignore_index=True)
    
    f_I_table = f_I_table.sort_values(by=["specimen", 'Stim_amp_pA'])
    f_I_table['specimen'] = pd.Categorical(f_I_table['specimen'])
    f_I_table['sweep_number']=np.int64(f_I_table['sweep_number'])
    f_I_table=pd.merge(f_I_table,coarse_or_fine_sweep(specimen_id),on='sweep_number')
    return f_I_table
    

def mysigmoid(x,maxi,x0,slope):
    y=maxi/(1+np.exp((x0-x)/slope))
    return y

def old_fit_sigmoid(f_I_table):
    '''
    Fit a sigmoid curve to the stim/amplitude data points to extract several I/O metrcis : the threshold, the saturation and the gain

    Parameters
    ----------
    f_I_table : DataFrame
        Stimulus_frequency table for one cell.

    Returns
    -------
    estimated_gain : float
        estimated gain of the I/O.
    estimated_threshold : float
        estimated neuron threshold.
    estimated_saturation : float
        estimated neuron saturation firing rate.
    my_plot : ggplot
        plot of the data point with the sigmoid fit and the linear fit to the linear part of the sigmoid.
    pcov: 2-D array
        The estimated covariance of popt
    popt: 1D array
        Estimated parameters of function fit
    
    '''
    try:
        x_data=f_I_table.iloc[:,2]
        y_data=f_I_table.iloc[:,3]
        
        
        ##Get the initial estimate for the fit of sigmoid
        #Get the maximum firing rate of the data
        maxi=max(x_data)

        #Get the index corresponding to the median non-zero firing rate
        

        without_zero_index=next(x for x, val in enumerate(y_data) if val >0 )
        median_firing_rate_index=next(x for x, val in enumerate(y_data) if val >= np.median(y_data.iloc[without_zero_index:]))
        #Get the stimulus amplitude correspondingto the median non-zero firing rate
    
        x0=x_data.iloc[median_firing_rate_index]

        #Get the slope from the linear fit of the firing rate
        slope=fit_specimen_fi_slope(x_data,y_data)[0]
    
    
        initial_estimate=[maxi,x0,slope]
        parameters_boundaries=([0,0,0],[np.inf,np.inf,np.inf])
        
        
        popt,pcov=curve_fit(mysigmoid,x_data,y_data,p0=initial_estimate,bounds=parameters_boundaries,check_finite=False)
        
        new_x_data=pd.Series(np.arange(min(x_data),max(x_data),1))
        new_y_data=pd.Series(mysigmoid(new_x_data,*popt))
        new_data=pd.concat([new_x_data,new_y_data],axis=1,ignore_index=True)
        new_data.columns=["stim_amplitude_pA","frequence_Hz"]
        #slope_confidence_threshold=5
        
        
        #get index 25% and 75% of max firing rate
        twentyfive_index=next(x for x, val in enumerate(new_y_data) if val >(0.25*popt[0]))
        seventyfive_index=next(x for x, val in enumerate(new_y_data) if val >(0.75*popt[0]))
        #fit linear line to linear sigmoid portion
        linear_estimated_slope,linear_estimated_intercept=fit_specimen_fi_slope(new_x_data.iloc[twentyfive_index:seventyfive_index],mysigmoid(new_x_data.iloc[twentyfive_index:seventyfive_index],*popt))
        estimated_threshold=(0-linear_estimated_intercept)/linear_estimated_slope
        my_derivative=np.array(derivative(mysigmoid,new_x_data,dx=1e-1,args=(popt[0],popt[1],popt[2])))
       
        if my_derivative[-1]<0.001:
            
            estimated_saturation=popt[0]
        else:
            estimated_saturation=np.nan
        estimated_gain=linear_estimated_slope
        my_plot=np.nan
        my_plot=ggplot(f_I_table,aes(x=f_I_table.columns[2],y=f_I_table.columns[3]))+geom_point()+geom_line(new_data,aes(x=new_data.columns[0],y=new_data.columns[1]),color='blue')+geom_abline(aes(intercept=linear_estimated_intercept,slope=linear_estimated_slope))
        my_plot+=geom_text(x=10,y=estimated_saturation,label="popt2="+str(round(popt[2],2))+
                              'id='+str(f_I_table.iloc[0,0])
                              ,size=10,color="black")
        #print(my_plot)

       
        return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt
    except (StopIteration):
        print("Stopped Iteration")

        estimated_gain=np.nan
        estimated_saturation=np.nan
        estimated_threshold=np.nan
        my_plot=np.nan
        pcov=np.nan
        popt=np.nan
        return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt
    except (ValueError):
        print("stopped_valueError")

        estimated_gain=np.nan
        estimated_saturation=np.nan
        estimated_threshold=np.nan
        my_plot=np.nan
        pcov=np.nan
        popt=np.nan
        return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt
    except (RuntimeError):
        print("Can't fit sigmoid, least-square optimization failed")

        estimated_gain=np.nan
        estimated_saturation=np.nan
        estimated_threshold=np.nan
        my_plot=np.nan
        pcov=np.nan
        popt=np.nan
        return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt

def random_id_list(id_list, n):
    max_value = len(id_list)
    new_id_list = []
    for elt in range(n):
        index = random.randint(1, max_value - 1)

        new_id_list.append(id_list[index])
    return new_id_list    
from scipy.misc import derivative
def trust_sigmoid(new_x_data,maxi,x0,slope,slope_confidence_threshold):
    
    my_derivative=np.array(derivative(mysigmoid,new_x_data,dx=1e-1,args=(maxi,x0,slope)))
   
    if max(my_derivative)>slope_confidence_threshold:
        print("slope too high: slope=",str(max(my_derivative)))
       
        return False
    elif my_derivative[-1]==0:
        print("end slope too low: endslope=",str(my_derivative[-1]))
        return False
    
    
    return True

def extract_inst_freq_table(specimen_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_spike_nb=False,first_nth_spikes=0):
    '''
    Compute the instananous frequency in each interspike interval per sweep for a cell

    Parameters
    ----------
    specimen_id : int
        specimencell id.
    species_sweep_stim_table : DataFrame
        Coming from create_species_sweeps_stim_table function.

    Returns
    -------
    inst_freq_table: DataFrame
        Table containing for a given cell for each sweep the stimulus amplitude and the instantanous frequency per interspike interval.

    '''
    index_stim = species_sweep_stim_table.columns.get_loc('Long Square')
    index_specimen = species_sweep_stim_table.index[species_sweep_stim_table["specimen_id"] == specimen_id][0]
    
    my_specimen_data = ctc.get_ephys_data(specimen_id)
    sweep_numbers = species_sweep_stim_table.iloc[index_specimen, index_stim]
    maximum_nb_interval =0
 
    for current_sweep in sweep_numbers:
        if per_time==True:
            index_range=my_specimen_data.get_sweep(current_sweep)["index_range"]
            sampling_rate=my_specimen_data.get_sweep(current_sweep)["sampling_rate"]
            current_stim_array=(my_specimen_data.get_sweep(current_sweep)["stimulus"][0:index_range[1]+1])* 1e12 #to pA
            stim_start_index=index_range[0]+next(x for x, val in enumerate(current_stim_array[index_range[0]:]) if val != 0 )
            current_time_array=np.arange(0, len(current_stim_array)) * (1.0 / sampling_rate)
            stim_start_time=current_time_array[stim_start_index]
            end_time=stim_start_time+(first_x_ms*1e-3)
            
            spike_times=my_specimen_data.get_spike_times(current_sweep)[my_specimen_data.get_spike_times(current_sweep) <= end_time ]
            if len(spike_times)>maximum_nb_interval:
                maximum_nb_interval=len(spike_times)
        elif per_spike_nb==True:
            spike_times=my_specimen_data.get_spike_times(current_sweep)[:first_nth_spikes]
            if len(spike_times)>maximum_nb_interval:
                maximum_nb_interval=len(spike_times)
        else:
            spike_times=my_specimen_data.get_spike_times(current_sweep)
            if len(spike_times)>maximum_nb_interval:
                maximum_nb_interval=len(spike_times)
       
    mycolumns=["specimen","sweep","stim_amplitude_pA"]+["interval_"+str(i) for i in range(1,(maximum_nb_interval))]
    inst_freq_table=pd.DataFrame(index=np.arange(len(sweep_numbers)),columns=mycolumns)
    for col in range(inst_freq_table.shape[1]):
        inst_freq_table.iloc[:,col]=np.nan
        
    for line in range(len(sweep_numbers)):
        current_sweep=sweep_numbers[line]
        stim_amplitude=my_specimen_data.get_sweep_metadata(current_sweep)['aibs_stimulus_amplitude_pa']
        if per_time==True:
            index_range=my_specimen_data.get_sweep(current_sweep)["index_range"]
            sampling_rate=my_specimen_data.get_sweep(current_sweep)["sampling_rate"]
            current_stim_array=(my_specimen_data.get_sweep(current_sweep)["stimulus"][0:index_range[1]+1])* 1e12 #to pA
            stim_start_index=index_range[0]+next(x for x, val in enumerate(current_stim_array[index_range[0]:]) if val != 0 )
            current_time_array=np.arange(0, len(current_stim_array)) * (1.0 / sampling_rate)
            stim_start_time=current_time_array[stim_start_index]
            end_time=stim_start_time+(first_x_ms*1e-3)
            
            spike_times=my_specimen_data.get_spike_times(current_sweep)[my_specimen_data.get_spike_times(current_sweep) <= end_time ]
           
        elif per_spike_nb==True:
            spike_times=my_specimen_data.get_spike_times(current_sweep)[:first_nth_spikes]
            
        else:
            spike_times=my_specimen_data.get_spike_times(current_sweep)
            
        
        inst_freq_table.iloc[line,0]=specimen_id
        inst_freq_table.iloc[line,1]=current_sweep   
        inst_freq_table.iloc[line,2]=stim_amplitude
        # Put a minimum number of spikes to compute adaptation
        if len(spike_times) >2:
            for current_spike_time_index in range(1,len(spike_times)):
                current_inst_frequency=1/(spike_times[current_spike_time_index]-spike_times[current_spike_time_index-1])
                
                inst_freq_table.iloc[line,(current_spike_time_index+2)]=current_inst_frequency
        
            inst_freq_table.iloc[line,3:]/=inst_freq_table.iloc[line,3]
    inst_freq_table = inst_freq_table.sort_values(by=["specimen", 'stim_amplitude_pA'])
    inst_freq_table['specimen']=pd.Categorical(inst_freq_table['specimen'])
    
    interval_freq_table=pd.DataFrame(columns=['interval','inst_frequency','stimulus_amplitude'])
    isnull_table=inst_freq_table.isnull()
    for col in range(3,(inst_freq_table.shape[1])):
        for line in range(inst_freq_table.shape[0]):
            if isnull_table.iloc[line,col] == False:
                new_line=pd.Series([int(col-2),inst_freq_table.iloc[line,col],np.float64(inst_freq_table.iloc[line,2])],
                                   index=['interval','inst_frequency','stimulus_amplitude'])
                interval_freq_table=interval_freq_table.append(new_line,ignore_index=True)
   
    specimen=pd.Series(np.array([inst_freq_table.iloc[0,0]]*interval_freq_table.shape[0]))
    interval_freq_table=pd.concat([specimen,interval_freq_table],axis=1)
    interval_freq_table.columns=["specimen",'interval','inst_frequency','stimulus_amplitude']
    interval_freq_table['specimen']=pd.Categorical(interval_freq_table['specimen'])

    return interval_freq_table
    # return(inst_freq_table)

def table_to_fit(inst_freq_table):
    '''
    Create table of interspike index-instantanous frequency- stimulus amplitude

    Parameters
    ----------
    inst_freq_table : DataFrame
        Coming from extract_inst_freq_table.

    Returns
    -------
    interval_freq_table : DataFrame
        Reorganized table to fit exponential function in fit_exponential_decay function.

    '''
    interval_freq_table=pd.DataFrame(columns=['interval','inst_frequency','stimulus_amplitude'])
    isnull_table=inst_freq_table.isnull()
    for col in range(3,(inst_freq_table.shape[1])):
        for line in range(inst_freq_table.shape[0]):
            if isnull_table.iloc[line,col] == False:
                new_line=pd.Series([int(col-2),inst_freq_table.iloc[line,col],np.float64(inst_freq_table.iloc[line,2])],
                                   index=['interval','inst_frequency','stimulus_amplitude'])
                interval_freq_table=interval_freq_table.append(new_line,ignore_index=True)
   
    specimen=pd.Series(np.array([inst_freq_table.iloc[0,0]]*interval_freq_table.shape[0]))
    interval_freq_table=pd.concat([specimen,interval_freq_table],axis=1)
    interval_freq_table.columns=["specimen",'interval','inst_frequency','stimulus_amplitude']
    interval_freq_table['specimen']=pd.Categorical(interval_freq_table['specimen'])

    return interval_freq_table


def my_exponential_decay(x,A,B,C):
    '''
    Parameters
    ----------
    x : Array
        interspike interval index array.
    A: flt
        initial instantanous frequency .
    B : flt
        Adaptation index constant.
    C : flt
        intantaneous frequency limit.

    Returns
    -------
    y : array
        Modelled instantanous frequency.

    '''
    y=A*np.exp(-(x-1)/B)+C
    
    return y


def expo_to_minimize(params,x_data,data):
    A=params['amplitude']
    B=params['decay']
    C=params['steady_state']
    model=A*np.exp(-(x_data-1)/B)+C
    return model-data
    
    
    
def new_fit_exponential_decay(cell_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_spike_nb=False,first_nth_spike=0,do_plot=False):
    '''
    Parameters
    ----------
    interval_freq_table : DataFrame
        Comming from table_to_fit function.

    Returns
    -------
    my_plot : ggplot
        
    starting_freq : flt
        estimated initial instantanous frequency.
    adapt_cst : flt
        Adaptation index constant.
    limit_freq : flt
        intantaneous frequency limit.
    pcov_overall : 2-D array
        The estimated covariance of popt

    '''
    
    try:
        
        if type(cell_id)!=np.float64:
            str_cell_id=str(cell_id)
            cell_id=np.float64(cell_id)
        else:
            str_cell_id=str(cell_id)
        interval_freq_table=extract_inst_freq_table(cell_id,species_sweep_stim_table,per_time=per_time,first_x_ms=first_x_ms,per_spike_nb=per_spike_nb,first_nth_spikes=first_nth_spike)
        if interval_freq_table.shape[0]==0:
            print('not enough spike')
            my_plot=np.nan
            A=np.nan
            B=np.nan
            C=np.nan

            return A,B,C,my_plot
        x_data=interval_freq_table.iloc[:,1]
        
        median_table=interval_freq_table.groupby(by=["interval"],dropna=True).median()
        median_table["count_weights"]=pd.DataFrame(interval_freq_table.groupby(by=["interval"],dropna=True).count()).iloc[:,-1]
        median_table["interval"]=np.arange(1,(median_table.shape[0]+1))
        median_table["interval"]=np.float64(median_table["interval"])  
        

        try:
            lower_index=next(x for x, val in enumerate(median_table["inst_frequency"]) if val<((min(median_table["inst_frequency"])+median_table["inst_frequency"][1])/2))
        except(StopIteration):
            lower_index=np.inf
        try:
        
            higher_index= next(x for x, val in enumerate(median_table["inst_frequency"]) if val>((max(median_table["inst_frequency"])+median_table["inst_frequency"][1])/2))
        except(StopIteration):
            higher_index=np.inf
        #med_index=next(x for x, val in enumerate(median_table["inst_frequency"]) if val<0.5*max(median_table['inst_frequency']))

        if lower_index<higher_index:
            med_index=lower_index
            initial_amplitude=1
        else:
            med_index=higher_index
            initial_amplitude=-1
        initial_decay_value=np.mean(median_table["interval"][med_index:med_index+1])

        

        decay_step=15
        while max(x_data)%decay_step==0:
            decay_step+=1
        params=Parameters()
        params.add('amplitude',value=initial_amplitude)
        params.add('decay',value=initial_decay_value,min=0.001)
        params.add('steady_state',value=0.1,min=0,expr="1-amplitude")
        
        params['amplitude'].set(brute_step=5)
        params['decay'].set(brute_step=decay_step)
        params['steady_state'].set(brute_step=0.1,min=0)
        
        
           
        #data=my_exponential_decay(x=median_table["interval"], A=0.5, B=2, C=0.1)
        fitter=Minimizer(expo_to_minimize,params,fcn_args=(median_table["interval"],median_table["inst_frequency"]))
        result_brute=fitter.minimize(method='brute',Ns=40,keep=40)


        
        
        # plt.plot(median_table["interval"],median_table["inst_frequency"],'o')
        # plt.plot(median_table["interval"],median_table["inst_frequency"]+expo_to_minimize(result_brute.params,median_table["interval"],median_table["inst_frequency"]),'--')
        # plt.plot(median_table["interval"],my_exponential_decay(median_table["interval"],1,-34,0.3),'x')
        # plt.show()
        
        # myparams=result_brute.candidates[1].params
        # myamplitude=myparams['amplitude'].value
        # mydecay=myparams['decay'].value
        # mysteady_state=myparams['steady_state'].value
        # print(myamplitude,mydecay,mysteady_state)
        
        best_chi=None
        
        for current_result in result_brute.candidates:
            current_amplitude=current_result.params['amplitude'].value
            current_decay=current_result.params['decay'].value
            current_steady_state=current_result.params['steady_state'].value

            
            exponential_model=ExponentialModel(prefix="expo_")
            expo_params=exponential_model.make_params()
            expo_params["expo_amplitude"].set(value=current_amplitude)
            expo_params["expo_decay"].set(value=current_decay,min=0.1)
            
            constant_model=ConstantModel(prefix="const_")
            expo_params.update(constant_model.make_params())
            expo_params["const_c"].set(value=current_steady_state,min=0,expr='1-expo_amplitude')
            
            full_expo_model=exponential_model+constant_model
            full_expo_out=full_expo_model.fit(median_table["inst_frequency"],expo_params,x=median_table["interval"],weights=median_table["count_weights"])
            
            
            if best_chi==None:
                best_chi=full_expo_out.chisqr
                best_A=full_expo_out.best_values["expo_amplitude"]
                best_B=full_expo_out.best_values["expo_decay"]
                best_C=full_expo_out.best_values["const_c"]
            elif best_chi>full_expo_out.chisqr:
                best_chi=full_expo_out.chisqr
                best_A=full_expo_out.best_values["expo_amplitude"]
                best_B=full_expo_out.best_values["expo_decay"]
                best_C=full_expo_out.best_values["const_c"]
        
        A_norm=best_A/(best_A+best_C)
        C_norm=best_C/(best_A+best_C)
        interval_range=np.arange(1,max(median_table["interval"])+1,1)

        simulation=my_exponential_decay(interval_range,best_A,best_B,best_C)
        norm_simulation=my_exponential_decay(interval_range,A_norm,best_B,C_norm)
        sim_table=pd.DataFrame(np.column_stack((interval_range,simulation)),columns=["Interval","Normalized_freq"])
        norm_sim_table=pd.DataFrame(np.column_stack((interval_range,norm_simulation)),columns=["Interval","Normalized_freq"])
        
        
        
        my_plot=np.nan
        if do_plot==True:
            plot_results_brute(result_brute,best_vals=True,varlabels=None)
            my_plot=ggplot(interval_freq_table,aes(x=interval_freq_table["interval"],y=interval_freq_table["inst_frequency"]))+geom_point(aes(color=interval_freq_table["stimulus_amplitude"]))
            
            my_plot=my_plot+geom_point(median_table,aes(x='interval',y='inst_frequency',size=median_table["count_weights"]),shape='s',color='red')
            my_plot=my_plot+geom_line(sim_table,aes(x='Interval',y='Normalized_freq'),color='black')
            my_plot=my_plot+geom_line(norm_sim_table,aes(x='Interval',y='Normalized_freq'),color="green")
            my_plot+=xlab(str("ISI index; cell_id: "+str_cell_id))
            print(my_plot)


        
        return best_A,best_B,best_C,my_plot
    except (StopIteration):
        print("Stopped Iteration")
        my_plot=np.nan
        A=np.nan
        B=np.nan
        C=np.nan

        return A,B,C,my_plot
    except (ValueError):
        print("stopped_valueError")
        my_plot=np.nan
        A=np.nan
        B=np.nan
        C=np.nan
        return A,B,C,my_plot
    except(RuntimeError):
        print("Can't fit exponential, least-square optimization failed")
        my_plot=np.nan
        A=np.nan
        B=np.nan
        C=np.nan
        return A,B,C,my_plot
    except(TypeError):
        print("Stopped TypeError")
        my_plot=np.nan
        A=np.nan
        B=np.nan
        C=np.nan
        return A,B,C,my_plot

def plot_results_brute(result, best_vals=True, varlabels=None,
                       output=None):
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square value per parameter and contour
    plots for all combination of two parameters.

    Inspired by the `corner` package (https://github.com/dfm/corner.py).

    Parameters
    ----------
    result : :class:`~lmfit.minimizer.MinimizerResult`
        Contains the results from the :meth:`brute` method.

    best_vals : bool, optional
        Whether to show the best values from the grid search (default is True).

    varlabels : list, optional
        If None (default), use `result.var_names` as axis labels, otherwise
        use the names specified in `varlabels`.

    output : str, optional
        Name of the output PDF file (default is 'None')
    """
    npars = len(result.var_names)
    _fig, axes = plt.subplots(npars, npars)

    if not varlabels:
        varlabels = result.var_names
    if best_vals and isinstance(best_vals, bool):
        best_vals = result.params

    for i, par1 in enumerate(result.var_names):
        for j, par2 in enumerate(result.var_names):

            # parameter vs chi2 in case of only one parameter
            if npars == 1:
                axes.plot(result.brute_grid, result.brute_Jout, 'o', ms=3)
                axes.set_ylabel(r'$\chi^{2}$')
                axes.set_xlabel(varlabels[i])
                if best_vals:
                    axes.axvline(best_vals[par1].value, ls='dashed', color='r')

            # parameter vs chi2 profile on top
            elif i == j and j < npars-1:
                if i == 0:
                    axes[0, 0].axis('off')
                ax = axes[i, j+1]
                red_axis = tuple(a for a in range(npars) if a != i)
                ax.plot(np.unique(result.brute_grid[i]),
                        np.minimum.reduce(result.brute_Jout, axis=red_axis),
                        'o', ms=3)
                ax.set_ylabel(r'$\chi^{2}$')
                ax.yaxis.set_label_position("right")
                ax.yaxis.set_ticks_position('right')
                ax.set_xticks([])
                if best_vals:
                    ax.axvline(best_vals[par1].value, ls='dashed', color='r')

            # parameter vs chi2 profile on the left
            elif j == 0 and i > 0:
                ax = axes[i, j]
                red_axis = tuple(a for a in range(npars) if a != i)
                ax.plot(np.minimum.reduce(result.brute_Jout, axis=red_axis),
                        np.unique(result.brute_grid[i]), 'o', ms=3)
                ax.invert_xaxis()
                ax.set_ylabel(varlabels[i],rotation="horizontal")

                if i != npars-1:
                    ax.set_xticks([])
                else:
                    ax.set_xlabel(r'$\chi^{2}$')
                if best_vals:
                    ax.axhline(best_vals[par1].value, ls='dashed', color='r')

            # contour plots for all combinations of two parameters
            elif j > i:
                ax = axes[j, i+1]
                red_axis = tuple(a for a in range(npars) if a not in (i, j))
                X, Y = np.meshgrid(np.unique(result.brute_grid[i]),
                                   np.unique(result.brute_grid[j]))
                lvls1 = np.linspace(result.brute_Jout.min(),
                                    np.median(result.brute_Jout)/2.0, 7, dtype='int')
                lvls2 = np.linspace(np.median(result.brute_Jout)/2.0,
                                    np.median(result.brute_Jout), 3, dtype='int')
                lvls = np.unique(np.concatenate((lvls1, lvls2)))
                ax.contourf(X.T, Y.T, np.minimum.reduce(result.brute_Jout, axis=red_axis),
                            lvls, norm=LogNorm())
                ax.set_yticks([])
                if best_vals:
                    ax.axvline(best_vals[par1].value, ls='dashed', color='r')
                    ax.axhline(best_vals[par2].value, ls='dashed', color='r')
                    ax.plot(best_vals[par1].value, best_vals[par2].value, 'rs', ms=3)
                if j != npars-1:
                    ax.set_xticks([])
                else:
                    ax.set_xlabel(varlabels[i],rotation="vertical")
                if j - i >= 2:
                    axes[i, j].axis('off')

    if output is not None:
        plt.savefig(output)


def fit_exponential_decay(cell_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_spike_nb=False,first_nth_spike=0,do_plot=False):
    '''
    Parameters
    ----------
    interval_freq_table : DataFrame
        Comming from table_to_fit function.

    Returns
    -------
    my_plot : ggplot
        
    starting_freq : flt
        estimated initial instantanous frequency.
    adapt_cst : flt
        Adaptation index constant.
    limit_freq : flt
        intantaneous frequency limit.
    pcov_overall : 2-D array
        The estimated covariance of popt

    '''
    
    try:
        if type(cell_id)!=np.float64:
            str_cell_id=str(cell_id)
            cell_id=np.float64(cell_id)
        else:
            str_cell_id=str(cell_id)
        interval_freq_table=extract_inst_freq_table(cell_id,species_sweep_stim_table,per_time=per_time,first_x_ms=first_x_ms,per_spike_nb=per_spike_nb,first_nth_spikes=first_nth_spike)
        x_data=interval_freq_table.iloc[:,1]
        y_data=interval_freq_table.iloc[:,2]
        initial_steady_state=np.mean(interval_freq_table[interval_freq_table['interval']==max(interval_freq_table['interval'])]['inst_frequency'])
        median_table=interval_freq_table.groupby(by=["interval"],dropna=True).median()
        median_table["count_weights"]=pd.DataFrame(interval_freq_table.groupby(by=["interval"],dropna=True).count()).iloc[:,-1]
        median_table["interval"]=np.arange(1,(median_table.shape[0]+1))
        median_table["interval"]=np.float64(median_table["interval"])      
        
        exponential_model=ExponentialModel(prefix="expo_")
        expo_params=exponential_model.make_params()
        expo_params["expo_amplitude"].set(value=-1)
        expo_params["expo_decay"].set(value=15)
        
        constant_model=ConstantModel(prefix="const_")
        expo_params.update(constant_model.make_params())
        expo_params["const_c"].set(value=initial_steady_state,min=0,expr='1-expo_amplitude')
        
        full_expo_model=exponential_model+constant_model
        full_expo_out=full_expo_model.fit(median_table["inst_frequency"],expo_params,x=median_table["interval"],weights=median_table["count_weights"])
        
        A=full_expo_out.best_values["expo_amplitude"]
        B=full_expo_out.best_values["expo_decay"]
        C=full_expo_out.best_values["const_c"]
        A_norm=A/(A+C)
        C_norm=C/(A+C)
        interval_range=np.arange(1,max(median_table["interval"])+1,1)

        simulation=my_exponential_decay(interval_range,A,B,C)
        norm_simulation=my_exponential_decay(interval_range,A_norm,B,C_norm)
        sim_table=pd.DataFrame(np.column_stack((interval_range,simulation)),columns=["Interval","Normalized_freq"])
        norm_sim_table=pd.DataFrame(np.column_stack((interval_range,norm_simulation)),columns=["Interval","Normalized_freq"])
        
        print('chi_square=',full_expo_out.chisqr)
        
        my_plot=np.nan
        if do_plot==True:
            my_plot=ggplot(interval_freq_table,aes(x=interval_freq_table["interval"],y=interval_freq_table["inst_frequency"]))+geom_point(aes(color=interval_freq_table["stimulus_amplitude"]))
            
            my_plot=my_plot+geom_point(median_table,aes(x='interval',y='inst_frequency',size=median_table["count_weights"]),shape='s',color='red')
            my_plot=my_plot+geom_line(sim_table,aes(x='Interval',y='Normalized_freq'),color='black')
            my_plot=my_plot+geom_line(norm_sim_table,aes(x='Interval',y='Normalized_freq'),color="green")
            my_plot+=xlab(str("ISI index; cell_id: "+str_cell_id))
            print(my_plot)
        return A,B,C,my_plot
    except (StopIteration):
        print("Stopped Iteration")
        my_plot=np.nan
        A=np.nan
        B=np.nan
        C=np.nan
        return A,B,C,my_plot
    except (ValueError):
        print("stopped_valueError")
        my_plot=np.nan
        A=np.nan
        B=np.nan
        C=np.nan
        return A,B,C,my_plot
    except(RuntimeError):
        print("Can't fit exponential, least-square optimization failed")
        my_plot=np.nan
        A=np.nan
        B=np.nan
        C=np.nan
        return A,B,C,my_plot
    except(TypeError):
        print("Stopped TypeError")
        my_plot=np.nan
        A=np.nan
        B=np.nan
        C=np.nan
        return A,B,C,my_plot


def fit_sigmoid (cell_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_nth_spike=False,first_nth_spike=0,do_plot=False):
    try:
         if type(cell_id)!=np.float64:
             str_cell_id=str(cell_id)
             cell_id=np.float64(cell_id)
         else:
             str_cell_id=str(cell_id)
         # extract f_I table for the specimen
         f_I_table=extract_stim_freq(float(cell_id), species_sweep_stim_table,per_time=per_time,first_x_ms=first_x_ms,per_nth_spike=per_nth_spike,first_nth_spike=first_nth_spike)
         x_data=f_I_table.iloc[:,2]
         y_data=f_I_table.iloc[:,3]
         
         #get initial estimate of parameters
         without_zero_index=next(x for x, val in enumerate(y_data) if val >0 )
         last_zero_index=without_zero_index-1
         median_firing_rate_index=next(x for x, val in enumerate(y_data) if val >= np.median(y_data.iloc[without_zero_index:]))
         #Get the stimulus amplitude correspondingto the median non-zero firing rate
         x0=x_data.iloc[median_firing_rate_index]
         #Get the slope from the linear fit of the firing rate
         slope,intercept=fit_specimen_fi_slope(x_data,y_data)
         last_zero_x=x_data.iloc[last_zero_index]
         first_non_zero_x=x_data.iloc[without_zero_index]
         Q3_index=next(x for x, val in enumerate(y_data) if val >0.75*max(y_data) )
         Q3=y_data.iloc[Q3_index]
         Q1_index=next(x for x, val in enumerate(y_data) if val >0.25*max(y_data) )
         Q1=y_data.iloc[Q1_index]
         # Define sigmoid function model
         sigmoid_mod=StepModel(form='logistic',prefix='sigmoid_')
         

         #params=sigmoid_mod.guess(y_data,x=x_data,center=x0)
         params=sigmoid_mod.make_params()
         params['sigmoid_amplitude'].set(value=max(y_data),min=0.1)
         params['sigmoid_center'].set(value=x0,min=last_zero_x)
         params['sigmoid_sigma'].set(value=100,min=0)
         out=sigmoid_mod.fit(y_data,params,x=x_data)
         general_table=pd.DataFrame(np.column_stack((x_data,y_data)),columns=["Stim_amp_pA","Frequency"])
         
         #Get parameters best estimations
         amplitude=out.best_values['sigmoid_amplitude']
         center=out.best_values['sigmoid_center']
         sigma=out.best_values['sigmoid_sigma']

        
         #Create more continuous x_data in the same range, compute corresponding y_data according to computed parameters
         new_x_data=pd.Series(np.arange(min(x_data),max(x_data),1))
         NRMSE_sigmoid=mean_squared_error(y_data, sigmoid_function(x_data,amplitude,center,sigma),squared=(False))/(Q3-Q1)

         
         
         #Make composite model fit by multipying the sigmoid with a step function
         new_sigmoid_model=StepModel(form='logistic',prefix='new_sigmoid_')
         new_sigmoid_params=new_sigmoid_model.make_params()
         new_sigmoid_params['new_sigmoid_amplitude'].set(value=max(y_data),min=max(y_data))
         new_sigmoid_params['new_sigmoid_center'].set(value=first_non_zero_x,min=first_non_zero_x)
         new_sigmoid_params['new_sigmoid_sigma'].set(value=100,min=0.0)     

         step_mod=StepModel(form='logistic',prefix="step_")
         new_sigmoid_params.update(step_mod.make_params())
         new_sigmoid_params['step_amplitude'].set(value=1,vary=False)
         new_sigmoid_params['step_center'].set(value=first_non_zero_x)
         new_sigmoid_params['step_sigma'].set(value=0.0001,vary=False) 
         
         compo_model=new_sigmoid_model*step_mod
         compo_out=compo_model.fit(y_data,new_sigmoid_params,x=x_data)
         # Get parameters best estimations
         new_sigmoid_amplitude=compo_out.best_values['new_sigmoid_amplitude']
         new_sigmoid_center=compo_out.best_values['new_sigmoid_center']
         new_sigmoid_sigma=compo_out.best_values['new_sigmoid_sigma']
         step_jump_x=compo_out.best_values['step_center']
         step_amplitude=1
         NRMSE_composite=mean_squared_error(y_data,pd.Series(compo(x_data,step_jump_x,step_amplitude,new_sigmoid_amplitude,new_sigmoid_center,new_sigmoid_sigma)),squared=False )/(Q3-Q1)
       
         if NRMSE_composite<=NRMSE_sigmoid:
             # in case where step function + sigmoid function fits better (Type II)
             computed_y_data=pd.Series(compo(x_data,step_jump_x,step_amplitude,new_sigmoid_amplitude,new_sigmoid_center,new_sigmoid_sigma))
            
             #get index 25% and 75% of max firing rate

             twentyfive_index=next(x for x, val in enumerate(computed_y_data) if val >(0.25*max(computed_y_data)))
             seventyfive_index=next(x for x, val in enumerate(computed_y_data) if val >(0.75*max(computed_y_data)))
             #fit linear line to linear sigmoid portion
             estimated_slope,estimated_intercept=fit_specimen_fi_slope(new_x_data.iloc[twentyfive_index:seventyfive_index],compo(new_x_data.iloc[twentyfive_index:seventyfive_index],step_jump_x,step_amplitude,new_sigmoid_amplitude,new_sigmoid_center,new_sigmoid_sigma))
             estimated_threshold=step_jump_x
             my_derivative=np.array(derivative(compo,new_x_data,dx=1e-1,args=(step_jump_x,step_amplitude,new_sigmoid_amplitude,new_sigmoid_center,new_sigmoid_sigma)))
             end_slope=np.mean(my_derivative[-10:])
             parameter_amplitude=new_sigmoid_amplitude
             parameter_center=new_sigmoid_center
             parameter_sigma=new_sigmoid_sigma
             if end_slope <=0.001:
                 estimated_saturation=np.mean(computed_y_data[-10:])

             else:
                 estimated_saturation=np.nan
             
             color_plot='green'
             Neuron_type='Type_II'
             linetype_composite='solid'
             linetype_sigmoid="dashed"
         else:
             # In case where simple sigmoid fits better (Type I)
             computed_y_data=sigmoid_function(new_x_data, amplitude, center, sigma)
             
             twentyfive_index=next(x for x, val in enumerate(computed_y_data) if val >(0.25*max(computed_y_data)))
             seventyfive_index=next(x for x, val in enumerate(computed_y_data) if val >(0.75*max(computed_y_data)))
             #fit linear line to linear sigmoid portion
             estimated_slope,estimated_intercept=fit_specimen_fi_slope(new_x_data.iloc[twentyfive_index:seventyfive_index],sigmoid_function(new_x_data.iloc[twentyfive_index:seventyfive_index],amplitude,center,sigma))
             estimated_threshold=(0-estimated_intercept)/estimated_slope
             my_derivative=np.array(derivative(sigmoid_function,new_x_data,dx=1e-1,args=(amplitude,center,sigma)))
             end_slope=np.mean(my_derivative[-10:])
             parameter_amplitude=amplitude
             parameter_center=center
             parameter_sigma=sigma
             if end_slope <=0.001:
                 estimated_saturation=np.mean(computed_y_data[-10:])

             else:
                 estimated_saturation=np.nan
             Neuron_type='Type_I'
             color_plot='blue'
             linetype_composite='dashed'
             linetype_sigmoid="solid"
         #if you want to see both curve (incompatible with do_plot=True)
         # computed_y_data=pd.Series(compo(new_x_data,step_jump_x,new_sigmoid_amplitude,new_sigmoid_center,new_sigmoid_sigma))
         # secondcomputed_y_data=sigmoid_function(new_x_data, amplitude, center, sigma)
         # my_model_table=pd.DataFrame(np.column_stack((new_x_data,computed_y_data)),columns=["Stim_amp_pA","Frequency"])
         # my_second_model_table=pd.DataFrame(np.column_stack((new_x_data,secondcomputed_y_data)),columns=["Stim_amp_pA","Frequency"])
         # my_new_plot=ggplot(general_table,aes(x=general_table["Stim_amp_pA"],y=general_table["Frequency"]))+geom_point()
         # my_new_plot+=geom_line(my_model_table,aes(x=my_model_table["Stim_amp_pA"],y=my_model_table['Frequency']),color='green',linetype=linetype_composite)
         # my_new_plot+=geom_line(my_second_model_table,aes(x=my_second_model_table["Stim_amp_pA"],y=my_second_model_table['Frequency']),color='blue',linetype=linetype_sigmoid)
         # my_new_plot+=xlab(str("Stim_amp_pA_id: "+str_cell_id))
         # print(my_new_plot)
         my_plot=np.nan
         if do_plot == True:
             model_table=pd.DataFrame(np.column_stack((new_x_data,computed_y_data)),columns=["Stim_amp_pA","Frequency"])
             my_plot=ggplot(general_table,aes(x=general_table["Stim_amp_pA"],y=general_table["Frequency"]))+geom_point()
             my_plot+=geom_line(model_table,aes(x=model_table["Stim_amp_pA"],y=model_table['Frequency']),color=color_plot)
             my_plot+=geom_abline(aes(intercept=estimated_intercept,slope=estimated_slope))
             my_plot+=xlab(str("Stim_amp_pA_id: "+str_cell_id))
             print(my_plot)
             
         
         #return (amplitude,center,sigma,end_slope,NRMSE_sigmoid,NRMSE_composite,step_jump_x,new_sigmoid_amplitude,new_sigmoid_center,new_sigmoid_sigma)
         return estimated_slope,estimated_threshold,estimated_saturation,Neuron_type,my_plot,NRMSE_sigmoid,NRMSE_composite,parameter_amplitude,parameter_center,parameter_sigma
    except(StopIteration):
         print("Stop Iteration")
         #return(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
         return np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
    except (ValueError):
         print("stopped_valueError")
         #return(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
         return np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
        
    except (RuntimeError):
         print("Can't fit sigmoid, least-square optimization failed")
         #return(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
         return np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
    
def sigmoid_to_minimize(params,x_data,data):
    sigmoid_amplitude=params['sigmoid_amplitude']
    sigmoid_center=params["sigmoid_center"]
    sigmoid_sigma=params["sigmoid_sigma"]
    step_amplitude=params['step_amplitude']
    step_center=params['step_center']
    step_sigma=0.00001
    model=(sigmoid_amplitude*(1-(1/(1+np.exp((x_data-sigmoid_center)/sigmoid_sigma)))))*(step_amplitude*(1-(1/(1+np.exp((x_data-step_center)/step_sigma)))))
    return model-data
    

def single_sigmoid_to_minimize(params,x_data,data):
    single_sigmoid_amplitude=params['single_sigmoid_amplitude']
    single_sigmoid_center=params['single_sigmoid_center']
    single_sigmoid_sigma=params['single_sigmoid_sigma']

    model=single_sigmoid_amplitude*(1-(1/(1+np.exp((x_data-single_sigmoid_center)/single_sigmoid_sigma))))

    return model-data

def new_fit_sigmoid (cell_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_nth_spike=False,first_nth_spike=0,do_plot=False):
    try:

         
         if type(cell_id)!=np.float64:
             str_cell_id=str(cell_id)
             cell_id=np.float64(cell_id)
         else:
             str_cell_id=str(cell_id)
         # extract f_I table for the specimen and use only the "coarse" annotated sweeps
         f_I_table=extract_stim_freq(float(cell_id), species_sweep_stim_table,per_time=per_time,first_x_ms=first_x_ms,per_nth_spike=per_nth_spike,first_nth_spike=first_nth_spike)
         coarse_f_I_table=f_I_table[f_I_table['stimulus_description']=='COARSE']
         if max(coarse_f_I_table.loc[:,'Frequency_Hz'])==0:
             fit='No_COARSE_Response'
             return (fit,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
         x_data=coarse_f_I_table.loc[:,'Stim_amp_pA']
         y_data=coarse_f_I_table.loc[:,"Frequency_Hz"]
         
         #get initial estimate of parameters for single sigmoid fit
         without_zero_index=next(x for x, val in enumerate(y_data) if val >0 )

         
         median_firing_rate_index=next(x for x, val in enumerate(y_data) if val >= np.median(y_data.iloc[without_zero_index:]))
         #Get the stimulus amplitude correspondingto the median non-zero firing rate
         x0=x_data.iloc[median_firing_rate_index]
         #Get the slope from the linear fit of the firing rate
         slope,intercept=fit_specimen_fi_slope(x_data,y_data)
         
         first_non_zero_x=x_data.iloc[without_zero_index]
         new_x_data=pd.Series(np.arange(min(x_data),max(x_data),1))
         first_non_zero_extended_x_index=next(x for x, val in enumerate(new_x_data) if val >=first_non_zero_x )
         best_single_amplitude=np.nan
         best_single_center=np.nan
         best_single_sigma=np.nan
         best_compo_QNRMSE=None
         best_step_amplitude=np.nan
         best_step_center=np.nan
         best_step_sigma=np.nan
         best_sigmoid_amplitude=np.nan
         best_sigmoid_center=np.nan
         best_sigmoid_sigma=np.nan

         ##First, try to fit a single sigmoid
         params_single_sigmoid=Parameters()
         params_single_sigmoid.add('single_sigmoid_amplitude',value=max(y_data),min=0)
         params_single_sigmoid.add('single_sigmoid_center',value=x0)
         params_single_sigmoid.add('single_sigmoid_sigma',value=500,min=0.1)
         params_single_sigmoid['single_sigmoid_amplitude'].set(brute_step=20)
         params_single_sigmoid['single_sigmoid_center'].set(brute_step=21)
         params_single_sigmoid['single_sigmoid_sigma'].set(brute_step=31)

         single_sigmoid_fitter=Minimizer(single_sigmoid_to_minimize,params_single_sigmoid,fcn_args=(x_data,y_data))

         single_result_brute=single_sigmoid_fitter.minimize(method='brute',Ns=10,keep=10)

         best_single_QNRMSE=None
         fit='Rejected'
         typeII_tested=False
         #plot_results_brute(single_result_brute,best_vals=True,varlabels=None)
         for current_single_result in single_result_brute.candidates:

             current_single_sigmoid_amplitude=current_single_result.params["single_sigmoid_amplitude"].value
             current_single_sigmoid_center=current_single_result.params["single_sigmoid_center"].value
             current_single_sigmoid_sigma=current_single_result.params["single_sigmoid_sigma"].value

             single_sigmoid_mod=StepModel(form='logistic',prefix='single_sigmoid_')
             single_sigmoid_mod_params=single_sigmoid_mod.make_params()
             single_sigmoid_mod_params['single_sigmoid_amplitude'].set(value=current_single_sigmoid_amplitude)
             single_sigmoid_mod_params['single_sigmoid_center'].set(value=current_single_sigmoid_center)
             single_sigmoid_mod_params['single_sigmoid_sigma'].set(value=current_single_sigmoid_sigma)
             
             single_sigmoid_out=single_sigmoid_mod.fit(y_data,single_sigmoid_mod_params,x=x_data)
             current_best_single_sigmoid_amplitude=single_sigmoid_out.best_values['single_sigmoid_amplitude']
             current_best_single_sigmoid_center=single_sigmoid_out.best_values['single_sigmoid_center']
             current_best_single_sigmoid_sigma=single_sigmoid_out.best_values['single_sigmoid_sigma']
             true=y_data.iloc[without_zero_index:]
             pred=pd.Series(sigmoid_function(x_data.iloc[without_zero_index:],current_best_single_sigmoid_amplitude,
                                                                                    current_best_single_sigmoid_center,
                                                                                    current_best_single_sigmoid_sigma))
             pred_extended=pd.Series(sigmoid_function(new_x_data.loc[first_non_zero_extended_x_index:],current_best_single_sigmoid_amplitude,
                                                                                    current_best_single_sigmoid_center,
                                                                                    current_best_single_sigmoid_sigma))

             if best_single_QNRMSE==None or best_single_QNRMSE>normalized_root_mean_squared_error(true,pred,pred_extended):
                 
                 best_single_amplitude=current_best_single_sigmoid_amplitude
                 best_single_center=current_best_single_sigmoid_center
                 best_single_sigma=current_best_single_sigmoid_sigma
                 true=y_data.iloc[without_zero_index:]
                 pred=pd.Series(sigmoid_function(x_data.iloc[without_zero_index:],current_best_single_sigmoid_amplitude,
                                                                                        current_best_single_sigmoid_center,
                                                                                        current_best_single_sigmoid_sigma))
                 pred_extended=pd.Series(sigmoid_function(new_x_data.loc[first_non_zero_extended_x_index:],current_best_single_sigmoid_amplitude,
                                                                                        current_best_single_sigmoid_center,
                                                                                        current_best_single_sigmoid_sigma))

                 best_single_QNRMSE=normalized_root_mean_squared_error(true,pred,pred_extended)

        
         single_sigmoid_y_data=pd.Series(sigmoid_function(new_x_data,best_single_amplitude,
                                                                                best_single_center,
                                                                                best_single_sigma))



         if best_single_QNRMSE<0.5:
             fit='TypeI'
         ##Define condition to test double sigmoid fit
         if best_single_QNRMSE<1e-3 or best_single_QNRMSE>0.5 or best_single_sigma<1:

             params=Parameters()
             params.add('sigmoid_amplitude',value=max(y_data),min=0)
             params.add('sigmoid_center',value=x0)
             params.add("sigmoid_sigma",value=500,min=0.1)
             params.add('step_amplitude',value=1,vary=False)
             params.add('step_center',value=first_non_zero_x,min=0)

             
             params['sigmoid_amplitude'].set(brute_step=80)
             params["sigmoid_center"].set(brute_step=60)
             params["sigmoid_sigma"].set(brute_step=120)
             params['step_amplitude'].set(brute_step=0.1)
             params['step_center'].set(brute_step=(max(x_data//9)))


             fitter=Minimizer(sigmoid_to_minimize,params,fcn_args=(x_data,y_data))

             result_brute=fitter.minimize(method='brute',Ns=15,keep=15)
             
            
    
    
             for current_results in result_brute.candidates:
                  current_sigmoid_amplitude=current_results.params['sigmoid_amplitude'].value
                  current_sigmoid_center=current_results.params['sigmoid_center'].value
                  current_sigmoid_sigma=current_results.params['sigmoid_sigma'].value
                  current_step_amplitude=current_results.params['step_amplitude'].value
                  current_step_center=current_results.params['step_center'].value
                 
                  new_sigmoid_model=StepModel(form='logistic',prefix='sigmoid_')
                  new_sigmoid_params=new_sigmoid_model.make_params()
                  new_sigmoid_params['sigmoid_amplitude'].set(value=current_sigmoid_amplitude)
                  new_sigmoid_params['sigmoid_center'].set(value=current_sigmoid_center)
                  new_sigmoid_params['sigmoid_sigma'].set(value=current_sigmoid_sigma)     
                 
                  step_mod=StepModel(form='logistic',prefix="step_")
                  new_sigmoid_params.update(step_mod.make_params())
                  new_sigmoid_params['step_amplitude'].set(value=current_step_amplitude)
                  new_sigmoid_params['step_center'].set(value=current_step_center)
                  new_sigmoid_params['step_sigma'].set(value=0.00001,vary=False) 
                 
                  compo_model=new_sigmoid_model*step_mod
    
                  compo_out=compo_model.fit(y_data,new_sigmoid_params,x=x_data)
    
                  # Get parameters best estimations
                  new_sigmoid_amplitude=compo_out.best_values['sigmoid_amplitude']
                  new_sigmoid_center=compo_out.best_values['sigmoid_center']
                  new_sigmoid_sigma=compo_out.best_values['sigmoid_sigma']
                  true=y_data.iloc[without_zero_index:]
                  pred=pd.Series(compo(x_data.iloc[without_zero_index:],compo_out.best_values['step_center'],compo_out.best_values['step_amplitude'],compo_out.best_values["sigmoid_amplitude"],compo_out.best_values["sigmoid_center"],compo_out.best_values["sigmoid_sigma"]))
                  pred_extended=pd.Series(compo(new_x_data.loc[first_non_zero_extended_x_index:],compo_out.best_values['step_center'],compo_out.best_values['step_amplitude'],compo_out.best_values["sigmoid_amplitude"],compo_out.best_values["sigmoid_center"],compo_out.best_values["sigmoid_sigma"]))
                 

                  if best_compo_QNRMSE==None or best_compo_QNRMSE>normalized_root_mean_squared_error(true,pred,pred_extended):
                     
                      true=y_data.iloc[without_zero_index:]
                      pred=pd.Series(compo(x_data.iloc[without_zero_index:],compo_out.best_values['step_center'],compo_out.best_values['step_amplitude'],compo_out.best_values["sigmoid_amplitude"],compo_out.best_values["sigmoid_center"],compo_out.best_values["sigmoid_sigma"]))
                      pred_extended=pd.Series(compo(new_x_data.loc[first_non_zero_extended_x_index:],compo_out.best_values['step_center'],compo_out.best_values['step_amplitude'],compo_out.best_values["sigmoid_amplitude"],compo_out.best_values["sigmoid_center"],compo_out.best_values["sigmoid_sigma"]))
                      best_compo_QNRMSE=normalized_root_mean_squared_error(true,pred,pred_extended)
                      
                      best_sigmoid_amplitude=compo_out.best_values["sigmoid_amplitude"]
                      best_sigmoid_center=compo_out.best_values["sigmoid_center"]
                      best_sigmoid_sigma=compo_out.best_values["sigmoid_sigma"]
                     
                      best_step_amplitude=compo_out.best_values['step_amplitude']
                      best_step_center=compo_out.best_values['step_center']
                      best_step_sigma=compo_out.best_values['step_sigma']
                      
             
             computed_y_data=pd.Series(compo(new_x_data,best_step_center,best_step_amplitude,best_sigmoid_amplitude,best_sigmoid_center,best_sigmoid_sigma))
             model_table=pd.DataFrame(np.column_stack((new_x_data,computed_y_data)),columns=["Stim_amp_pA","Frequency_Hz"])
             typeII_tested=True

             if 2*best_compo_QNRMSE<=best_single_QNRMSE and best_compo_QNRMSE<0.1 and best_compo_QNRMSE>1e-4:
                 fit= 'TypeII'
             elif best_single_QNRMSE<0.5 and best_single_sigma>1:
                 fit='TypeI'
             else:
                 fit='Rejected'
                 
         color_dict={"COARSE":"green",
                     "FINEST":'red'}

         if best_compo_QNRMSE==None:
             best_compo_QNRMSE=np.nan

         

         if do_plot == True:
              single_sigmoid_table=pd.DataFrame(np.column_stack((new_x_data,single_sigmoid_y_data)),columns=["Stim_amp_pA","Frequency_Hz"])
              my_plot=ggplot(f_I_table,aes(x=f_I_table["Stim_amp_pA"],y=f_I_table["Frequency_Hz"],color=f_I_table["stimulus_description"]))+geom_point()+scale_color_manual(values=color_dict)
              if fit=='TypeII':
                   compo_line='solid'
                   single_line='dashed'
                   my_plot+=geom_line(model_table,aes(x=model_table["Stim_amp_pA"],y=model_table['Frequency_Hz']),color='red',linetype=compo_line)
              elif fit=='TypeI':
                   compo_line='dashed'
                   single_line='solid'

                   if typeII_tested==True:
                       print('jijij')
                       my_plot+=geom_line(model_table,aes(x=model_table["Stim_amp_pA"],y=model_table['Frequency_Hz']),color='red',linetype=compo_line)

              else:
                  single_line="dashed"
              
              my_plot+=geom_line(single_sigmoid_table,aes(x=single_sigmoid_table["Stim_amp_pA"],y=single_sigmoid_table['Frequency_Hz']),color='blue',linetype=single_line)

              my_plot+=xlab(str("Stim_amp_pA_id: "+str_cell_id))
              print(my_plot)
              
         end_time=time.time()
         
         #print("Total time= ",end_time-start_time,"s ; Time compo brute=",time_1-start_time,"s ; Time single brute=",time_2-time_1)
         return fit,best_single_QNRMSE,best_single_amplitude,best_single_center,best_single_sigma,best_compo_QNRMSE,best_step_amplitude,best_step_center,best_step_sigma,best_sigmoid_amplitude,best_sigmoid_center,best_sigmoid_sigma
         
    except(StopIteration):
         print("Stop Iteration")
         return("Failed",np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
         
    except (ValueError):
         print("stopped_valueError")
         
         return("Failed",np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
        
    except (RuntimeError):
         print("Can't fit sigmoid, least-square optimization failed")
        
         return('Failed',np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
    except (TypeError):
         print("Stop Type Error")
         return('Failed',np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
     
        
 
def heaviside_fit_sigmoid (cell_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_nth_spike=False,first_nth_spike=0,do_plot=False):
    try:

         
         if type(cell_id)!=np.float64:
             str_cell_id=str(cell_id)
             cell_id=np.float64(cell_id)
         else:
             str_cell_id=str(cell_id)
         # extract f_I table for the specimen and use only the "coarse" annotated sweeps
         f_I_table=extract_stim_freq(float(cell_id), species_sweep_stim_table,per_time=per_time,first_x_ms=first_x_ms,per_nth_spike=per_nth_spike,first_nth_spike=first_nth_spike)
         coarse_f_I_table=f_I_table[f_I_table['stimulus_description']=='COARSE']
         if max(coarse_f_I_table.loc[:,'Frequency_Hz'])==0:
             fit='No_COARSE_Response'
             return (fit,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
         x_data=coarse_f_I_table.loc[:,'Stim_amp_pA']
         y_data=coarse_f_I_table.loc[:,"Frequency_Hz"]
         
         #get initial estimate of parameters for single sigmoid fit
         without_zero_index=next(x for x, val in enumerate(y_data) if val >0 )

         
         median_firing_rate_index=next(x for x, val in enumerate(y_data) if val >= np.median(y_data.iloc[without_zero_index:]))
         #Get the stimulus amplitude correspondingto the median non-zero firing rate
         x0=x_data.iloc[median_firing_rate_index]
         #Get the slope from the linear fit of the firing rate
         slope,intercept=fit_specimen_fi_slope(x_data,y_data)
         
         first_non_zero_x=x_data.iloc[without_zero_index]
         new_x_data=pd.Series(np.arange(min(x_data),max(x_data),1))
         first_non_zero_extended_x_index=next(x for x, val in enumerate(new_x_data) if val >=first_non_zero_x )
         best_single_amplitude=np.nan
         best_single_center=np.nan
         best_single_sigma=np.nan
         best_compo_QNRMSE=None
         best_step_amplitude=np.nan
         best_step_center=np.nan
         best_step_sigma=np.nan
         best_sigmoid_amplitude=np.nan
         best_sigmoid_center=np.nan
         best_sigmoid_sigma=np.nan

         ##First, try to fit a single sigmoid
         params_single_sigmoid=Parameters()
         params_single_sigmoid.add('single_sigmoid_amplitude',value=max(y_data),min=0)
         params_single_sigmoid.add('single_sigmoid_center',value=x0)
         params_single_sigmoid.add('single_sigmoid_sigma',value=500,min=0.1)
         params_single_sigmoid['single_sigmoid_amplitude'].set(brute_step=20)
         params_single_sigmoid['single_sigmoid_center'].set(brute_step=21)
         params_single_sigmoid['single_sigmoid_sigma'].set(brute_step=31)

         single_sigmoid_fitter=Minimizer(single_sigmoid_to_minimize,params_single_sigmoid,fcn_args=(x_data,y_data))

         single_result_brute=single_sigmoid_fitter.minimize(method='brute',Ns=10,keep=10)

         best_single_QNRMSE=None
         fit='Rejected'
         typeII_tested=False

         #plot_results_brute(single_result_brute,best_vals=True,varlabels=None)
         for current_single_result in single_result_brute.candidates:

             current_single_sigmoid_amplitude=current_single_result.params["single_sigmoid_amplitude"].value
             current_single_sigmoid_center=current_single_result.params["single_sigmoid_center"].value
             current_single_sigmoid_sigma=current_single_result.params["single_sigmoid_sigma"].value

             single_sigmoid_mod=StepModel(form='logistic',prefix='single_sigmoid_')
             single_sigmoid_mod_params=single_sigmoid_mod.make_params()
             single_sigmoid_mod_params['single_sigmoid_amplitude'].set(value=current_single_sigmoid_amplitude)
             single_sigmoid_mod_params['single_sigmoid_center'].set(value=current_single_sigmoid_center)
             single_sigmoid_mod_params['single_sigmoid_sigma'].set(value=current_single_sigmoid_sigma)
             
             single_sigmoid_out=single_sigmoid_mod.fit(y_data,single_sigmoid_mod_params,x=x_data)
             current_best_single_sigmoid_amplitude=single_sigmoid_out.best_values['single_sigmoid_amplitude']
             current_best_single_sigmoid_center=single_sigmoid_out.best_values['single_sigmoid_center']
             current_best_single_sigmoid_sigma=single_sigmoid_out.best_values['single_sigmoid_sigma']
             true=y_data.iloc[without_zero_index:]
             pred=pd.Series(sigmoid_function(x_data.iloc[without_zero_index:],current_best_single_sigmoid_amplitude,
                                                                                    current_best_single_sigmoid_center,
                                                                                    current_best_single_sigmoid_sigma))
             pred_extended=pd.Series(sigmoid_function(new_x_data.loc[first_non_zero_extended_x_index:],current_best_single_sigmoid_amplitude,
                                                                                    current_best_single_sigmoid_center,
                                                                                    current_best_single_sigmoid_sigma))

             if best_single_QNRMSE==None or best_single_QNRMSE>normalized_root_mean_squared_error(true,pred,pred_extended):
                 
                 best_single_amplitude=current_best_single_sigmoid_amplitude
                 best_single_center=current_best_single_sigmoid_center
                 best_single_sigma=current_best_single_sigmoid_sigma
                 true=y_data.iloc[without_zero_index:]
                 pred=pd.Series(sigmoid_function(x_data.iloc[without_zero_index:],current_best_single_sigmoid_amplitude,
                                                                                        current_best_single_sigmoid_center,
                                                                                        current_best_single_sigmoid_sigma))
                 pred_extended=pd.Series(sigmoid_function(new_x_data.loc[first_non_zero_extended_x_index:],current_best_single_sigmoid_amplitude,
                                                                                        current_best_single_sigmoid_center,
                                                                                        current_best_single_sigmoid_sigma))

                 best_single_QNRMSE=normalized_root_mean_squared_error(true,pred,pred_extended)

        
         single_sigmoid_y_data=pd.Series(sigmoid_function(new_x_data,best_single_amplitude,
                                                                                best_single_center,
                                                                                best_single_sigma))



         if best_single_QNRMSE<0.5:
             fit='TypeI'
             

         ##Define condition to test double sigmoid fit
         if best_single_QNRMSE<1e-3 or best_single_QNRMSE>0.5 or best_single_sigma<1:

             params=Parameters()
             params.add('sigmoid_amplitude',value=max(y_data),min=0)
             params.add('sigmoid_center',value=x0)
             params.add("sigmoid_sigma",value=104,min=40)
             #params.add('heaviside_amplitude',value=1,min=0.1)
             params.add('heaviside_step',value=first_non_zero_x)

             
             params['sigmoid_amplitude'].set(brute_step=80)
             params["sigmoid_center"].set(brute_step=30)
             params["sigmoid_sigma"].set(brute_step=20)
             #params['heaviside_amplitude'].set(brute_step=0.1)
             params['heaviside_step'].set(brute_step=5)


             fitter=Minimizer(sigmoid_heaviside_to_minimize,params,fcn_args=(x_data,y_data))

             result_brute=fitter.minimize(method='brute',Ns=20,keep=20)
             

    
    
             for current_results in result_brute.candidates:
                  current_sigmoid_amplitude=current_results.params['sigmoid_amplitude'].value
                  current_sigmoid_center=current_results.params['sigmoid_center'].value
                  current_sigmoid_sigma=current_results.params['sigmoid_sigma'].value
                  #current_heaviside_amplitude=current_results.params['heaviside_amplitude'].value
                  current_heaviside_step=current_results.params['heaviside_step'].value
                 
                  composite_model=Model(Heaviside_function)*Model(sigmoid_function)
                  pars=composite_model.make_params(amplitude=current_sigmoid_amplitude,
                                                   center=current_sigmoid_center,
                                                   sigma=current_sigmoid_sigma,
                                                   mid=current_heaviside_step)
                 
    
                  compo_out=composite_model.fit(y_data,pars,x=x_data)
    
                  # Get parameters best estimations
                  
                  true=y_data.iloc[without_zero_index:]

                  pred=pd.Series(sigmoid_heaviside(x_data.iloc[without_zero_index:],compo_out.best_values['amplitude'],compo_out.best_values['center'],compo_out.best_values["sigma"],compo_out.best_values["mid"]))

                  pred_extended=pd.Series(sigmoid_heaviside(new_x_data.loc[first_non_zero_extended_x_index:],compo_out.best_values['amplitude'],compo_out.best_values['center'],compo_out.best_values["sigma"],compo_out.best_values["mid"]))


                  if best_compo_QNRMSE==None or best_compo_QNRMSE>normalized_root_mean_squared_error(true,pred,pred_extended):

                      true=y_data.iloc[without_zero_index:]
                      pred=pd.Series(sigmoid_heaviside(x_data.iloc[without_zero_index:],compo_out.best_values['amplitude'],compo_out.best_values['center'],compo_out.best_values["sigma"],compo_out.best_values["mid"]))
                      pred_extended=pd.Series(sigmoid_heaviside(new_x_data.loc[first_non_zero_extended_x_index:],compo_out.best_values['amplitude'],compo_out.best_values['center'],compo_out.best_values["sigma"],compo_out.best_values["mid"]))
                     
                      best_compo_QNRMSE=normalized_root_mean_squared_error(true,pred,pred_extended)
                      
                      best_sigmoid_amplitude=compo_out.best_values["amplitude"]
                      best_sigmoid_center=compo_out.best_values["center"]
                      best_sigmoid_sigma=compo_out.best_values["sigma"]
                      best_heaviside_step=compo_out.best_values['mid']
                      
                      

             computed_y_data=pd.Series(pd.Series(sigmoid_heaviside(new_x_data,best_sigmoid_amplitude,best_sigmoid_center,best_sigmoid_sigma,best_heaviside_step)))
            
             model_table=pd.DataFrame(np.column_stack((new_x_data,computed_y_data)),columns=["Stim_amp_pA","Frequency_Hz"])

             typeII_tested=True

             if 2*best_compo_QNRMSE<=best_single_QNRMSE and best_compo_QNRMSE<0.1 and best_compo_QNRMSE>1e-4:
                 fit= 'TypeII'
             elif best_single_QNRMSE<0.5 and best_single_sigma>1:
                 fit='TypeI'
             else:
                 fit='Rejected'
                 
         color_dict={"COARSE":"green",
                     "FINEST":'red'}

         if best_compo_QNRMSE==None:
             best_compo_QNRMSE=np.nan

         

         if do_plot == True:
              single_sigmoid_table=pd.DataFrame(np.column_stack((new_x_data,single_sigmoid_y_data)),columns=["Stim_amp_pA","Frequency_Hz"])
              my_plot=ggplot(f_I_table,aes(x=f_I_table["Stim_amp_pA"],y=f_I_table["Frequency_Hz"],color=f_I_table["stimulus_description"]))+geom_point()+scale_color_manual(values=color_dict)
              if fit=='TypeII':
                   compo_line='solid'
                   single_line='dashed'
                   my_plot+=geom_line(model_table,aes(x=model_table["Stim_amp_pA"],y=model_table['Frequency_Hz']),color='red',linetype=compo_line)
              elif fit=='TypeI':
                   compo_line='dashed'
                   single_line='solid'

                   if typeII_tested==True:

                       my_plot+=geom_line(model_table,aes(x=model_table["Stim_amp_pA"],y=model_table['Frequency_Hz']),color='red',linetype=compo_line)

              else:
                  single_line="dashed"
              
              my_plot+=geom_line(single_sigmoid_table,aes(x=single_sigmoid_table["Stim_amp_pA"],y=single_sigmoid_table['Frequency_Hz']),color='blue',linetype=single_line)

              my_plot+=xlab(str("Stim_amp_pA_id: "+str_cell_id))
              print(my_plot)
              
         end_time=time.time()
         
         #print("Total time= ",end_time-start_time,"s ; Time compo brute=",time_1-start_time,"s ; Time single brute=",time_2-time_1)
         return fit,best_single_QNRMSE,best_single_amplitude,best_single_center,best_single_sigma,best_compo_QNRMSE,best_heaviside_step,best_sigmoid_amplitude,best_sigmoid_center,best_sigmoid_sigma
         
    except(StopIteration):
         print("Stop Iteration")
         return("Failed",np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
         
    except (ValueError):
         print("stopped_valueError")
         
         return("Failed",np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
        
    except (RuntimeError):
         print("Can't fit sigmoid, least-square optimization failed")
        
         return('Failed',np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
    except (TypeError):
         print("Stop Type Error")
         return('Failed',np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
     

def heaviside_compute_f_I_params(fit_table,species_sweep_stim_table,per_time=False,first_x_ms=0,per_nth_spike=False,first_nth_spike=0,do_plot=False):
    mycol=["Cell_id","Fit","Gain","Threshold","Saturation"]
    f_I_params_table=pd.DataFrame(columns=mycol)
    for current_cell_id in fit_table.loc[:,"Cell_id"]:
        current_fit_value=fit_table[fit_table['Cell_id']==str(current_cell_id)].fit.values[0]
        
        current_f_I_table=extract_stim_freq(float(current_cell_id), species_sweep_stim_table,per_time=per_time,first_x_ms=first_x_ms,per_nth_spike=per_nth_spike,first_nth_spike=first_nth_spike)
        coarse_f_I_table=current_f_I_table[current_f_I_table['stimulus_description']=='COARSE']
        
        x_data=coarse_f_I_table.loc[:,'Stim_amp_pA']
        new_x_data=pd.Series(np.arange(min(x_data),max(x_data),0.1))
        if current_fit_value == 'TypeI':
            
            single_amplitude_value=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_single_amplitude.values[0]
            single_center_value=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_single_center.values[0]
            single_sigma_value=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_single_sigma.values[0]
            single_sigmoid_y_data=pd.Series(sigmoid_function(new_x_data,single_amplitude_value,single_center_value,single_sigma_value))
                                                                                   
            
            twentyfive_index=next(x for x, val in enumerate(single_sigmoid_y_data) if val >(0.25*max(single_sigmoid_y_data)))
            seventyfive_index=next(x for x, val in enumerate(single_sigmoid_y_data) if val >(0.75*max(single_sigmoid_y_data)))
            #fit linear line to linear sigmoid portion
            Gain,Intercept=fit_specimen_fi_slope(new_x_data.iloc[twentyfive_index:seventyfive_index],sigmoid_function(new_x_data.iloc[twentyfive_index:seventyfive_index],single_amplitude_value,single_center_value,single_sigma_value))
            Threshold=(0-Intercept)/Gain
            extended_f_I_table=pd.DataFrame(np.column_stack((new_x_data,single_sigmoid_y_data)),columns=["Stim_amp_pA","Frequency_Hz"])
            my_derivative=np.array(derivative(sigmoid_function,new_x_data,dx=1e-1,args=(single_amplitude_value,single_center_value,single_sigma_value)))
            end_slope=np.mean(my_derivative[-10:])
            Saturation=np.nan
            if end_slope <=0.001:
                Saturation=np.mean(single_sigmoid_y_data[-10:])
                
            if do_plot==True:
                
                myplot=ggplot(coarse_f_I_table,aes(x=coarse_f_I_table["Stim_amp_pA"],y=coarse_f_I_table["Frequency_Hz"]))+geom_point()
                myplot+=geom_line(extended_f_I_table,aes(x=extended_f_I_table["Stim_amp_pA"],y=extended_f_I_table["Frequency_Hz"]),color='blue')
                myplot+=geom_abline(aes(intercept=Intercept,slope=Gain))
                Threshold_table=pd.DataFrame({'Stim_amp_pA':[Threshold],'Frequency_Hz':[0]})
                myplot+=geom_point(Threshold_table,aes(x=Threshold_table["Stim_amp_pA"],y=Threshold_table["Frequency_Hz"]),color='green')
                if Saturation!=np.nan:
                    myplot+=geom_abline(aes(intercept=Saturation,slope=0))
                
                
                
        elif current_fit_value == 'TypeII':
            compo_step_amplitude=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_step_amplitude.values[0]
            compo_step_center=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_step_center.values[0]
            compo_sigmoid_amplitude=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_sigmoid_amplitude.values[0]
            compo_sigmoid_center=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_sigmoid_center.values[0]
            compo_sigmoid_sigma=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_sigmoid_sigma.values[0]
            new_y_data=pd.Series(compo(new_x_data,compo_step_center,compo_step_amplitude,compo_sigmoid_amplitude,compo_sigmoid_center,compo_sigmoid_sigma))
         
            extended_f_I_table=pd.DataFrame(np.column_stack((new_x_data,new_y_data)),columns=["Stim_amp_pA","Frequency_Hz"])
           
            second_derivative=new_y_data.diff().diff()
            linear_portion_start_index=next(x for x, val in enumerate(second_derivative) if val <0)
            linear_portion_start_index+=1
            linear_portion_y_data=new_y_data.iloc[linear_portion_start_index:]
            linear_portion_x_data=new_x_data.iloc[linear_portion_start_index:]
            
            twentyfive_index=next(x for x, val in enumerate(linear_portion_y_data) if val >((min(linear_portion_y_data)+0.25*(max(linear_portion_y_data)-min(linear_portion_y_data)))))
            seventyfive_index=next(x for x, val in enumerate(linear_portion_y_data) if val >((min(linear_portion_y_data)+0.75*(max(linear_portion_y_data)-min(linear_portion_y_data)))))
            
            Gain,Intercept=fit_specimen_fi_slope(linear_portion_x_data.iloc[twentyfive_index:seventyfive_index],linear_portion_y_data.iloc[twentyfive_index:seventyfive_index])
            first_derivative=new_y_data.diff()
            Threshold=compo_step_center
            Saturation=np.nan
            if np.mean(first_derivative[-100:]) <=0.001:
                Saturation=np.mean(first_derivative[-100:])
            
            
        else:
            Gain,Threshold,Saturation=np.nan,np.nan,np.nan
            
        if do_plot==True and current_fit_value=='TypeI' or do_plot==True and current_fit_value=='TypeII':
            
            myplot=ggplot(coarse_f_I_table,aes(x=coarse_f_I_table["Stim_amp_pA"],y=coarse_f_I_table["Frequency_Hz"]))+geom_point()
            myplot+=geom_line(extended_f_I_table,aes(x=extended_f_I_table["Stim_amp_pA"],y=extended_f_I_table["Frequency_Hz"]),color='blue')
            myplot+=geom_abline(aes(intercept=Intercept,slope=Gain))
            Threshold_table=pd.DataFrame({'Stim_amp_pA':[Threshold],'Frequency_Hz':[0]})
            myplot+=geom_point(Threshold_table,aes(x=Threshold_table["Stim_amp_pA"],y=Threshold_table["Frequency_Hz"]),color='green')
            if Saturation!=np.nan:
                myplot+=geom_abline(aes(intercept=Saturation,slope=0))
            myplot+=xlab(str("Stim_amp_pA_id: "+str(current_cell_id)))
            print(myplot)
        new_line=pd.Series([str(current_cell_id),current_fit_value,Gain,Threshold,Saturation],
                       index=mycol)
        f_I_params_table=f_I_params_table.append(new_line,ignore_index=True)
        
        
    return(f_I_params_table)
        



def compute_f_I_params(fit_table,species_sweep_stim_table,per_time=False,first_x_ms=0,per_nth_spike=False,first_nth_spike=0,do_plot=False):
    mycol=["Cell_id","Fit","Gain","Threshold","Saturation"]
    f_I_params_table=pd.DataFrame(columns=mycol)
    for current_cell_id in fit_table.loc[:,"Cell_id"]:
        current_fit_value=fit_table[fit_table['Cell_id']==str(current_cell_id)].fit.values[0]
        
        current_f_I_table=extract_stim_freq(float(current_cell_id), species_sweep_stim_table,per_time=per_time,first_x_ms=first_x_ms,per_nth_spike=per_nth_spike,first_nth_spike=first_nth_spike)
        coarse_f_I_table=current_f_I_table[current_f_I_table['stimulus_description']=='COARSE']
        
        x_data=coarse_f_I_table.loc[:,'Stim_amp_pA']
        new_x_data=pd.Series(np.arange(min(x_data),max(x_data),0.1))
        if current_fit_value == 'TypeI':
            
            single_amplitude_value=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_single_amplitude.values[0]
            single_center_value=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_single_center.values[0]
            single_sigma_value=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_single_sigma.values[0]
            single_sigmoid_y_data=pd.Series(sigmoid_function(new_x_data,single_amplitude_value,single_center_value,single_sigma_value))
                                                                                   
            
            twentyfive_index=next(x for x, val in enumerate(single_sigmoid_y_data) if val >(0.25*max(single_sigmoid_y_data)))
            seventyfive_index=next(x for x, val in enumerate(single_sigmoid_y_data) if val >(0.75*max(single_sigmoid_y_data)))
            #fit linear line to linear sigmoid portion
            Gain,Intercept=fit_specimen_fi_slope(new_x_data.iloc[twentyfive_index:seventyfive_index],sigmoid_function(new_x_data.iloc[twentyfive_index:seventyfive_index],single_amplitude_value,single_center_value,single_sigma_value))
            Threshold=(0-Intercept)/Gain
            extended_f_I_table=pd.DataFrame(np.column_stack((new_x_data,single_sigmoid_y_data)),columns=["Stim_amp_pA","Frequency_Hz"])
            my_derivative=np.array(derivative(sigmoid_function,new_x_data,dx=1e-1,args=(single_amplitude_value,single_center_value,single_sigma_value)))
            end_slope=np.mean(my_derivative[-10:])
            Saturation=np.nan
            if end_slope <=0.001:
                Saturation=np.mean(single_sigmoid_y_data[-10:])
                
            if do_plot==True:
                
                myplot=ggplot(coarse_f_I_table,aes(x=coarse_f_I_table["Stim_amp_pA"],y=coarse_f_I_table["Frequency_Hz"]))+geom_point()
                myplot+=geom_line(extended_f_I_table,aes(x=extended_f_I_table["Stim_amp_pA"],y=extended_f_I_table["Frequency_Hz"]),color='blue')
                myplot+=geom_abline(aes(intercept=Intercept,slope=Gain))
                Threshold_table=pd.DataFrame({'Stim_amp_pA':[Threshold],'Frequency_Hz':[0]})
                myplot+=geom_point(Threshold_table,aes(x=Threshold_table["Stim_amp_pA"],y=Threshold_table["Frequency_Hz"]),color='green')
                if Saturation!=np.nan:
                    myplot+=geom_abline(aes(intercept=Saturation,slope=0))
                
                
                
        elif current_fit_value == 'TypeII':
            compo_step_amplitude=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_step_amplitude.values[0]
            compo_step_center=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_step_center.values[0]
            compo_sigmoid_amplitude=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_sigmoid_amplitude.values[0]
            compo_sigmoid_center=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_sigmoid_center.values[0]
            compo_sigmoid_sigma=fit_table[fit_table['Cell_id']==str(current_cell_id)].best_sigmoid_sigma.values[0]
            new_y_data=pd.Series(compo(new_x_data,compo_step_center,compo_step_amplitude,compo_sigmoid_amplitude,compo_sigmoid_center,compo_sigmoid_sigma))
         
            extended_f_I_table=pd.DataFrame(np.column_stack((new_x_data,new_y_data)),columns=["Stim_amp_pA","Frequency_Hz"])
           
            second_derivative=new_y_data.diff().diff()
            linear_portion_start_index=next(x for x, val in enumerate(second_derivative) if val <0)
            linear_portion_start_index+=1
            linear_portion_y_data=new_y_data.iloc[linear_portion_start_index:]
            linear_portion_x_data=new_x_data.iloc[linear_portion_start_index:]
            
            twentyfive_index=next(x for x, val in enumerate(linear_portion_y_data) if val >((min(linear_portion_y_data)+0.25*(max(linear_portion_y_data)-min(linear_portion_y_data)))))
            seventyfive_index=next(x for x, val in enumerate(linear_portion_y_data) if val >((min(linear_portion_y_data)+0.75*(max(linear_portion_y_data)-min(linear_portion_y_data)))))
            
            Gain,Intercept=fit_specimen_fi_slope(linear_portion_x_data.iloc[twentyfive_index:seventyfive_index],linear_portion_y_data.iloc[twentyfive_index:seventyfive_index])
            first_derivative=new_y_data.diff()
            Threshold=compo_step_center
            Saturation=np.nan
            if np.mean(first_derivative[-100:]) <=0.001:
                Saturation=np.mean(first_derivative[-100:])
            
            
        else:
            Gain,Threshold,Saturation=np.nan,np.nan,np.nan
            
        if do_plot==True and current_fit_value=='TypeI' or do_plot==True and current_fit_value=='TypeII':
            
            myplot=ggplot(coarse_f_I_table,aes(x=coarse_f_I_table["Stim_amp_pA"],y=coarse_f_I_table["Frequency_Hz"]))+geom_point()
            myplot+=geom_line(extended_f_I_table,aes(x=extended_f_I_table["Stim_amp_pA"],y=extended_f_I_table["Frequency_Hz"]),color='blue')
            myplot+=geom_abline(aes(intercept=Intercept,slope=Gain))
            Threshold_table=pd.DataFrame({'Stim_amp_pA':[Threshold],'Frequency_Hz':[0]})
            myplot+=geom_point(Threshold_table,aes(x=Threshold_table["Stim_amp_pA"],y=Threshold_table["Frequency_Hz"]),color='green')
            if Saturation!=np.nan:
                myplot+=geom_abline(aes(intercept=Saturation,slope=0))
            myplot+=xlab(str("Stim_amp_pA_id: "+str(current_cell_id)))
            print(myplot)
        new_line=pd.Series([str(current_cell_id),current_fit_value,Gain,Threshold,Saturation],
                       index=mycol)
        f_I_params_table=f_I_params_table.append(new_line,ignore_index=True)
        
        
    return(f_I_params_table)
        
def composite_function(x,sigmoid_amplitude,sigmoid_center,sigmoid_sigma,step_amplitude,step_center,step_sigma):
    y=(sigmoid_amplitude*(1-(1/(1+np.exp((x-sigmoid_center)/sigmoid_sigma)))))*(step_amplitude*(1-(1/(1+np.exp((x-step_center)/step_sigma)))))
    return y

def sigmoid_function(x,amplitude,center,sigma):
     y=amplitude*(1-(1/(1+np.exp((x-center)/sigma))))
     return y   
def my_expression(x,step_center,step_amplitude):
    y=step_amplitude*(1-(1/(1+np.exp((x-step_center)/0.00001))))
    return y
def compo(x,step_center,step_amplitude,amplitude,center,sigma):
    y=sigmoid_function(x,amplitude,center,sigma)*my_expression(x,step_center,step_amplitude)
    return y

def sigmoid_heaviside(x,sigmoid_amplitude,sigmoid_center,sigmoid_sigma,heaviside_step):
    heaviside=Heaviside_function(x,heaviside_step)
    y=(sigmoid_amplitude*(1-(1/(1+np.exp((x-sigmoid_center)/sigmoid_sigma)))))*heaviside
    return y
def Heaviside_function(x, mid):
    """Heaviside step function."""
    
    if mid<=min(x):
        o=np.ones(x.size)
        return o
    elif mid>=max(x):
        o=np.zeros(x.size)
        return o
    else:
        o = np.zeros(x.size)
        
        imid = max(np.where(x < mid)[0])
        
        o[imid:] = 1
        return o

def sigmoid_heaviside_to_minimize(params,x_data,data):
    sigmoid_amplitude=params['sigmoid_amplitude']
    sigmoid_center=params["sigmoid_center"]
    sigmoid_sigma=params["sigmoid_sigma"]
    heaviside_step=params["heaviside_step"]
   # heaviside_amplitude=params["heaviside_amplitude"]
    
   

    heaviside=Heaviside_function(x_data,heaviside_step)
    
    model=(sigmoid_amplitude*(1-(1/(1+np.exp((x_data-sigmoid_center)/sigmoid_sigma)))))*heaviside
    return model

def relative_squared_error(true, pred):
    true_mean = np.mean(true)
    squared_error_num = np.sum(np.square(true - pred))
    squared_error_den = np.sum(np.square(true - true_mean))
    rse_loss = squared_error_num / squared_error_den
    return rse_loss

def normalized_root_mean_squared_error(true, pred,pred_extended):
    #Normalization by the interquartile range
    squared_error = np.square((true - pred))
    sum_squared_error = np.sum(squared_error)
    rmse = np.sqrt(sum_squared_error / true.size)
    Q1=np.percentile(pred_extended,25)
    Q3=np.percentile(pred_extended,75)
    #print('Q1=',Q1,", Q3=",Q3)
    nrmse_loss = rmse/(Q3-Q1)
    return nrmse_loss

def relative_root_mean_squared_error(true, pred):
    num = np.sum(np.square(true - pred))
    den = np.sum(np.square(pred))
    squared_error = num/den
    rrmse_loss = np.sqrt(squared_error)
    return rrmse_loss

def coarse_or_fine_sweep(cell_id):
    all_sweeps_table=pd.DataFrame(ctc.get_ephys_sweeps(cell_id))
    all_sweeps_table=all_sweeps_table[all_sweeps_table['stimulus_name']=="Long Square"]
    all_sweeps_table=all_sweeps_table[['sweep_number','stimulus_description']]
    for elt in range(all_sweeps_table.shape[0]):
        if 'COARSE' in all_sweeps_table.iloc[elt,1]:
            all_sweeps_table.iloc[elt,1]='COARSE'
        elif 'FINEST' in all_sweeps_table.iloc[elt,1]:
            all_sweeps_table.iloc[elt,1]='FINEST'    
    return(all_sweeps_table)
    
def take_upper(string):
    upper = ''
    for char in string:
    #check uppercase characters
        if char.isupper():
            upper += char
# print uppercase characters
    return(upper)
def take_until_number(string):
    upper = ''
    for char in string:
    #check uppercase characters
        if char.isnumeric():
            return upper
        else:
            upper += char
# print uppercase characters
    return(upper)

#%%
start_time=time.time()
total_time=0
number_done=0
mylist=mouse_id_list
mycol=['cell_id',"fit","best_single_QNRMSE","best_single_amplitude","best_single_center","best_single_sigma","best_compo_QNRMSE","best_step_amplitude","best_step_center","best_step_sigma","best_sigmoid_amplitude","best_sigmoid_center","best_sigmoid_sigma"]
new_fit_table=pd.DataFrame(columns=mycol)
for cell_id in mylist:
    start_time=time.time()
    new_line=pd.Series([str(cell_id),*new_fit_sigmoid(cell_id,mouse_sweep_stim_table,per_time=True,first_x_ms=1000,do_plot=False)],
                   index=mycol)
    new_fit_table=new_fit_table.append(new_line,ignore_index=True)
    end_time=time.time()
    number_done+=1
    number_remaining=len(mylist)-number_done
    total_time+=end_time-start_time
    remaining_time=number_remaining*total_time/number_done
    total_remaining_minutes=remaining_time//60
    remaining_hours=total_remaining_minutes//60
    remaining_minutes=total_remaining_minutes%60
    remaining_seconds=remaining_time%60
    print('Remaining time: ',remaining_hours,'h ',remaining_minutes,'min ',round(remaining_seconds,0),"s")

print("Total time= ",total_time)
f_I_table=compute_f_I_params(new_fit_table,mouse_sweep_stim_table,per_time=True,first_x_ms=1000,per_nth_spike=False,first_nth_spike=0,do_plot=False)
#%%
#plt.hist(new_fit_table['bestNRMSE'],bins=600)
plt.scatter(new_fit_table['bestNRMSE'],new_fit_table['best_sigmoid_sigma'],alpha=0.05)
plt.ylim([-50,1500])
#plt.xlim([-10,1000])
plt.show()
#%% in case of a decreasing phase - to be fixed
 
         # max_y_index=next(x for x, val in enumerate(y_data) if val == max(y_data))
         # max_y_x=x_data.iloc[max_y_index]
         # #Make a composite model to account for an possible decrease of Firing rate
         # increase_sigmoid_model=StepModel(form='logistic',prefix='increase_phase_')
         # increase_sigmoid_params=increase_sigmoid_model.make_params()
         # increase_sigmoid_params['increase_phase_amplitude'].set(value=max(y_data),min=max(y_data))
         # increase_sigmoid_params['increase_phase_center'].set(value=first_non_zero_x,min=first_non_zero_x)
         # increase_sigmoid_params['increase_phase_sigma'].set(value=100,min=0.0)
         
         
         # decrease_sigmoid_model=StepModel(form='logistic',prefix='decrease_phase_')
         # increase_sigmoid_params.update(decrease_sigmoid_model.make_params())
         # increase_sigmoid_params['decrease_phase_amplitude'].set(value=max(y_data),expr='0-increase_phase_amplitude')
         # increase_sigmoid_params['decrease_phase_center'].set(value=max(x_data), min=max(x_data))
         # increase_sigmoid_params['decrease_phase_sigma'].set(value=100)
         
         # double_phase_model=increase_sigmoid_model+decrease_sigmoid_model
         # double_phase_out=double_phase_model.fit(y_data,increase_sigmoid_params,x=x_data)
         # increase_amplitude=double_phase_out.best_values['increase_phase_amplitude']
         # increase_center=double_phase_out.best_values['increase_phase_center']
         # increase_sigma=double_phase_out.best_values['increase_phase_sigma']
         
         # decrease_amplitude=double_phase_out.best_values['decrease_phase_amplitude']
         # decrease_center=double_phase_out.best_values['decrease_phase_center']
         # decrease_sigma=double_phase_out.best_values['decrease_phase_sigma']
         
         # double_y_data=pd.Series(sigmoid_function(new_x_data,increase_amplitude,increase_center,increase_sigma)-sigmoid_function(new_x_data,decrease_amplitude,decrease_center,decrease_sigma))
         # double_phase_table=pd.concat([new_x_data,double_y_data],axis=1,ignore_index=True)
         # double_phase_table.columns=["Stim_amp_pA","Frequency"]
         # my_plot+=geom_line(double_phase_table,aes(x='Stim_amp_pA',y='Frequency'))
         # print(double_phase_table)
#%% My modified function

def modified_fit_sigmoid(f_I_table):
    '''
    Fit a sigmoid curve to the stim/amplitude data points to extract several I/O metrcis : the threshold, the saturation and the gain

    Parameters
    ----------
    f_I_table : DataFrame
        Stimulus_frequency table for one cell.

    Returns
    -------
    estimated_gain : float
        estimated gain of the I/O.
    estimated_threshold : float
        estimated neuron threshold.
    estimated_saturation : float
        estimated neuron saturation firing rate.
    my_plot : ggplot
        plot of the data point with the sigmoid fit and the linear fit to the linear part of the sigmoid.
    pcov: 2-D array
        The estimated covariance of popt
    popt: 1D array
        Estimated parameters of function fit
    
    '''
    try:
        x_data=f_I_table.iloc[:,2]
        y_data=f_I_table.iloc[:,3]
        
        
        ##Get the initial estimate for the fit of sigmoid
        #Get the maximum firing rate of the data
        maxi=max(x_data)
        
        #Get the index corresponding to the median non-zero firing rate
        

        without_zero_index=next(x for x, val in enumerate(y_data) if val >0 )
        median_firing_rate_index=next(x for x, val in enumerate(y_data) if val >= np.median(y_data.iloc[without_zero_index:]))
        #Get the stimulus amplitude correspondingto the median non-zero firing rate
    
        x0=x_data.iloc[median_firing_rate_index]
    
        #Get the slope from the linear fit of the firing rate
        slope=fit_specimen_fi_slope(x_data,y_data)[0]
    
        
        initial_estimate=[maxi,x0,slope]
        parameters_boundaries=([0,0,0],[np.inf,np.inf,np.inf])
        

        popt,pcov=curve_fit(mysigmoid,x_data,y_data,p0=initial_estimate,bounds=parameters_boundaries,check_finite=False)
        
        new_x_data=pd.Series(np.arange(min(x_data),max(x_data),1))
        new_y_data=pd.Series(mysigmoid(new_x_data,*popt))
        new_data=pd.concat([new_x_data,new_y_data],axis=1,ignore_index=True)
        new_data.columns=["stim_amplitude_pA","frequence_Hz"]
        estimated_saturation=popt[0]
        twentyfive_index=next(x for x, val in enumerate(new_y_data) if val >(0.25*popt[0]))
        seventyfive_index=next(x for x, val in enumerate(new_y_data) if val >(0.75*popt[0]))
        #fit linear line to linear sigmoid portion
        linear_estimated_slope,linear_estimated_intercept=fit_specimen_fi_slope(new_x_data.iloc[twentyfive_index:seventyfive_index],mysigmoid(new_x_data.iloc[twentyfive_index:seventyfive_index],*popt))
        estimated_threshold=(0-linear_estimated_intercept)/linear_estimated_slope
        
        my_plot=ggplot(f_I_table,aes(x=f_I_table.columns[2],y=f_I_table.columns[3]))+geom_point()+geom_line(new_data,aes(x=new_data.columns[0],y=new_data.columns[1]),color='blue')+geom_abline(aes(intercept=linear_estimated_intercept,slope=linear_estimated_slope))
        my_plot+=geom_text(x=10,y=10,label="popt2="+str(round(popt[2],2))+
                              'id='+str(f_I_table.iloc[0,0])
                              ,size=10,color="black")
        my_derivative=np.array(derivative(mysigmoid,new_x_data,dx=1e-1,args=(popt[0],popt[1],popt[2])))
        end_slope=my_derivative[-1]
        print(my_plot)
        return popt[0],popt[1],popt[2],end_slope
        #slope_confidence_threshold=5
        
        # if trust_sigmoid(new_x_data, *popt, slope_confidence_threshold)==False:
        #     print("failed to fit a sigmoid")
        #     estimated_gain=np.nan

        #     #estimated_saturation=np.nan
        #     #estimated_threshold=np.nan
        #     #my_plot=np.nan
        #     #pcov=np.nan
        #     #popt=np.nan
        #     ##########
        #     estimated_saturation=popt[0]
        #     twentyfive_index=next(x for x, val in enumerate(new_y_data) if val >(0.25*popt[0]))
        #     seventyfive_index=next(x for x, val in enumerate(new_y_data) if val >(0.75*popt[0]))
        #     #fit linear line to linear sigmoid portion
        #     linear_estimated_slope,linear_estimated_intercept=fit_specimen_fi_slope(new_x_data.iloc[twentyfive_index:seventyfive_index],mysigmoid(new_x_data.iloc[twentyfive_index:seventyfive_index],*popt))
        #     estimated_threshold=(0-linear_estimated_intercept)/linear_estimated_slope
        #     my_derivative=np.array(derivative(mysigmoid,new_x_data,dx=1e-1,args=(popt[0],popt[1],popt[2])))
        #     my_plot=ggplot(f_I_table,aes(x=f_I_table.columns[2],y=f_I_table.columns[3]))+geom_point()+geom_line(new_data,aes(x=new_data.columns[0],y=new_data.columns[1]),color='blue')+geom_abline(aes(intercept=linear_estimated_intercept,slope=linear_estimated_slope))
        #     my_plot+=geom_text(x=10,y=estimated_saturation,label="popt2="+str(round(popt[2],2))+
        #                          'id='+str(f_I_table.iloc[0,0])
        #                          ,size=10,color="black")
        #     print('non accepted',popt[2])
        #     print(my_plot)
        #     ######
        #     return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt
        # else:

        #     #get index 25% and 75% of max firing rate
        #     twentyfive_index=next(x for x, val in enumerate(new_y_data) if val >(0.25*popt[0]))
        #     seventyfive_index=next(x for x, val in enumerate(new_y_data) if val >(0.75*popt[0]))
        #     #fit linear line to linear sigmoid portion
        #     linear_estimated_slope,linear_estimated_intercept=fit_specimen_fi_slope(new_x_data.iloc[twentyfive_index:seventyfive_index],mysigmoid(new_x_data.iloc[twentyfive_index:seventyfive_index],*popt))
        #     estimated_threshold=(0-linear_estimated_intercept)/linear_estimated_slope
        #     my_derivative=np.array(derivative(mysigmoid,new_x_data,dx=1e-1,args=(popt[0],popt[1],popt[2])))
           
        #     if my_derivative[-1]<0.001:
                
        #         estimated_saturation=popt[0]
        #     else:
        #         estimated_saturation=np.nan
        #     estimated_gain=linear_estimated_slope
        #     my_plot=np.nan
        #     my_plot=ggplot(f_I_table,aes(x=f_I_table.columns[2],y=f_I_table.columns[3]))+geom_point()+geom_line(new_data,aes(x=new_data.columns[0],y=new_data.columns[1]),color='blue')+geom_abline(aes(intercept=linear_estimated_intercept,slope=linear_estimated_slope))
        #     my_plot+=geom_text(x=10,y=estimated_saturation,label="popt2="+str(round(popt[2],2))+
        #                          'id='+str(f_I_table.iloc[0,0])
        #                          ,size=10,color="black")
            
        #     print("accepted",popt[2])
        #     print(my_plot)
        #     return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt
    except (StopIteration):
        print("Stopped Iteration")

        estimated_gain=np.nan
        estimated_saturation=np.nan
        estimated_threshold=np.nan
        my_plot=np.nan
        pcov=np.nan
        popt=np.nan
        return np.nan,np.nan,np.nan,np.nan
    except (ValueError):
        print("stopped_valueError")

        estimated_gain=np.nan
        estimated_saturation=np.nan
        estimated_threshold=np.nan
        my_plot=np.nan
        pcov=np.nan
        popt=np.nan
        return np.nan,np.nan,np.nan,np.nan
    except (RuntimeError):
        print("Can't fit sigmoid, least-square optimization failed")

        estimated_gain=np.nan
        estimated_saturation=np.nan
        estimated_threshold=np.nan
        my_plot=np.nan
        pcov=np.nan
        popt=np.nan
        return np.nan,np.nan,np.nan,np.nan
    
#%%Mytestcell
mycol=['Cell_id','Amlitude','center','sigma','end_slope']
units=['--','Hz','pA','Hz/pA','Hz/pA']

start_time=time.time()
first_line=pd.Series(mycol,index=mycol)
second_line=pd.Series(units,index=mycol)
first_two_lines=pd.DataFrame([first_line,second_line])
mytable=pd.DataFrame(columns=mycol)
for current_specimen_id in mouse_id_list:
    
    Maxi,x0,slope,end_slope=mytest_curve_fitting(current_specimen_id,mouse_sweep_stim_table)

    new_line=pd.Series([str(current_specimen_id),
                        round(Maxi,2),
                        round(x0,2),
                        round(slope,2),
                        round(end_slope,2)],
                        
                        index=mycol)
    #print(new_line)
    
    mytable=mytable.append(new_line, ignore_index=True)

#mytable=compute_feature(mouse_id_list,mouse_sweep_stim_table,per_time=True,first_x_ms=1000)
mytable=mytable.reindex(columns=mycol)
end_1000ms=time.time()
print("time for 1000ms=",end_1000ms-start_time,"s")
mytable=pd.concat([first_two_lines,mytable])

#%%treat_Lyle_Data
def lyle_extract_inst_freq(spike_time_table):
    
    response_number=spike_time_table.shape[0]
    interval_number=spike_time_table.shape[1]
    mycolumns=["specimen","response",'stimulus']+["inst_freq_"+str(i) for i in range(1,(interval_number+1))]

    
    
    inst_freq_table=pd.DataFrame(index=np.arange(response_number),columns=mycolumns)
    for col in range(3,inst_freq_table.shape[1]):
        inst_freq_table.iloc[:,col]=np.nan
        
    for line in range(response_number):

        current_response=spike_time_table.iloc[line,:]
        try:
            non_nan_index=next(x for x, val in enumerate(current_response) if str(val)=='nan')
            current_response=current_response.iloc[:non_nan_index]
        except(StopIteration):
            print(' no nan')



        if len(current_response)>2:

            for col in range(len(current_response)):
                inst_freq_table.iloc[line,col+3]=current_response.iloc[0]/current_response.iloc[col]
        
            
        
        inst_freq_table.iloc[line,0]='lyle_data'
        inst_freq_table.iloc[line,1]=str(line)   
        inst_freq_table.iloc[line,2]=str(line) 
        # Put a minimum number of spikes to compute adaptation
    

    inst_freq_table['specimen']=pd.Categorical(inst_freq_table['specimen'])
    
    return(inst_freq_table)


def lyle_table_to_fit(inst_freq_table):
    '''
    Create table of interspike index-instantanous frequency- stimulus amplitude

    Parameters
    ----------
    inst_freq_table : DataFrame
        Coming from extract_inst_freq_table.

    Returns
    -------
    interval_freq_table : DataFrame
        Reorganized table to fit exponential function in fit_exponential_decay function.

    '''
    interval_freq_table=pd.DataFrame(columns=['interval','inst_frequency'])
    
    for col in range(3,(inst_freq_table.shape[1])):
        for line in range(inst_freq_table.shape[0]):
            if str(inst_freq_table.iloc[line,col]) != 'nan':
                new_line=pd.Series([int(col-2),inst_freq_table.iloc[line,col]],
                                   index=['interval','inst_frequency'])
                interval_freq_table=interval_freq_table.append(new_line,ignore_index=True)
   
    specimen=pd.Series(np.array(['lyle_data']*interval_freq_table.shape[0]))
    interval_freq_table=pd.concat([specimen,interval_freq_table],axis=1)
    interval_freq_table.columns=["specimen",'interval','inst_frequency']
   # interval_freq_table['specimen']=pd.Categorical(interval_freq_table['specimen'])

    return interval_freq_table
