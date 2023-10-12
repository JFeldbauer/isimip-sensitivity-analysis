# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 09:46:31 2023

@author: tkan & jf
"""


import subprocess
import sys
import mikeio
import SALib
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt, matplotlib.dates as mdates
import numpy as np
from SALib.analyze import delta

# load results
folder = "../data/"
results_df = pd.read_csv(folder+'results_lhc.csv')
results_head = results_df.head(10)

## remove rmse larger than 20
rmse_flter = results_df['rmse'] < 20
results_df = results_df[rmse_flter]

## remove r that are inf or nas
r_flter = np.invert(results_df['r'].isna())
results_df = results_df[r_flter]

## remove nmae that are inf or nas
nmae_flter = results_df['nmae'] == np.inf
results_df['nmae'][nmae_flter] = 1e2

# get all lake names 
lakes = results_df['lake'].unique()
models = results_df['model'].unique()

# add random dummy variable
results_df['dummy'] = np.random.rand(results_df.shape[0], 1)

# define parameters for analyses
pars_gotm = ['wind_speed', 'swr', 'turb_param.k_min', 'bottom.h0b',
             'turb_param.const_num', 'Kw', 'dummy']
pars_glm = ['wind_speed', 'swr', 'mixing.coef_mix_hyp', 'mixing.coef_mix_conv',
            'mixing.coef_mix_turb', 'Kw', 'dummy']
pars_flake = ['wind_speed', 'swr', 'c_relax_C', 'fetch_lk',
              'depth_bs_lk', 'Kw', 'dummy']
pars_simstrat = ['wind_speed', 'swr', 'a_seiche', 'hgeo',
                 'cd', 'Kw', 'dummy']

# model specific bounds
bound_gotm = [[0.25, 1.5],
              [0.7, 1.3],
              [1.4e-07, 10.0e-06],
              [0.025, 0.75],
              [0.000250, 0.000750],
              [0.2310, 0.43],
              [0, 1]]
              
bound_glm = [[0.25, 1.5],
             [0.7, 1.3],
             [0.1, 2],
             [0.1, 3],
             [0.35, 0.65],
             [0.2310, 0.43],
             [0, 1]]
              
bound_flake = [[0.25, 1.5],
               [0.7, 1.3],
               [0.0001, 0.01],
               [500, 3000],
               [2, 8],
               [0.2310, 0.43],
               [0, 1]]
              
bound_simstrat = [[0.25, 1.5],
                  [0.7, 1.3],
                  [0.0008, 0.003],
                  [0.0, 0.5],
                  [0.00075, 0.000325],
                  [0.2310, 0.43],
                  [0, 1]]

#for par in pars:
#    print(result_test[par].describe())

# define stat metric
stat = ['rmse', 'r', 'nse', 'bias', 'mae', 'nmae']


# data frame for the results
res = pd.DataFrame()

# for loop over all lakes
for l in lakes:
  
  print('started lake: ' + l)
  # loop over the four models
  for m in models:
    
    # set parameters to model specific parameters
    if m == 'GLM':
      pars = pars_glm
      bound = bound_glm
    
    if m == 'GOTM':
      pars = pars_gotm
      bound = bound_gotm
    
    if m == 'Simstrat':
      pars = pars_simstrat
      bound = bound_simstrat
    
    if m == 'FLake':
      pars = pars_flake
      bound = bound_flake
    
    lake_flter = results_df['lake'] == l
    model_flter = results_df['model'] == m
    result_test_ = results_df[lake_flter & model_flter]

    result_test = result_test_.dropna(how='all', axis=1)
    
    # update range of Kw for selected lake
    bound[5] = [result_test['Kw'].min(), result_test['Kw'].max()]
    
    # SA
    problem = {
        'num_vars': len(pars),
        'names': pars,
        'bounds': bound
        }
    
    # loop over the different statistic metrics
    for s in stat:
      
      X = result_test[pars].to_numpy()
    
      # select stat metrics to compare against
      Y = result_test[s].to_numpy()
      
      Si = delta.analyze(problem, X , Y)
      
      si_df = pd.DataFrame(Si)
      
      delta_sum = si_df['delta'].sum()
      
      si_df['var'] = s
      si_df['lake'] = l
      si_df['model'] = m
      
      
      si_df = si_df.sort_values('delta', ascending=False)
      
      res = pd.concat([res, si_df])
    

# write to file
res.to_csv('res_sens.csv', index=False)





# filter out dummy parameter and calculates top conf interval of max conf of par. less than dummy
#dummy_filt = si_df[si_df['names'] == 'dummy']
#dummy_val = dummy_filt['delta'].iloc[-1] + dummy_filt['delta_conf'].iloc[-1]
#less_dummy = si_df[si_df['delta'] < dummy_val]
#sens_threshold = dummy_filt['delta'].iloc[-1] + less_dummy['delta_conf'].max()
#si_df['sensitive'] = si_df['delta'] >= sens_threshold

### plot

# def plot_delta(si_df, SA_var, var, site, output_folder=None):
#         si_df = si_df.sort_values(SA_var, ascending=False)
# 
#         plt.bar(si_df['names'], si_df[SA_var])
#         plt.errorbar(si_df['names'], si_df[SA_var], yerr=si_df[SA_var+'_conf'], fmt=".", color="r", alpha=0.5)
#         
#         """
#         if SA_var == 'delta':
#             plt.plot(si_df['names'], [dummy_val]*len(si_df['names']), '--', color='gray')
#             plt.plot(si_df['names'], [sens_threshold]*len(si_df['names']), '-', color='gray')      
#         """
#         
#         plt.xticks(rotation=75)
#         plt.title(var)
#         plt.ylabel(SA_var)
#         
#         if output_folder:
#             plt.savefig(output_folder+SA_var+'_'+site+'_'+var+'.png', dpi=400)
#         
#         plt.show()
# 
# plot_folder = folder
# 
# plot_delta(si_df, 'delta', stat, lake, output_folder=plot_folder)
# plot_delta(si_df, 'S1', stat, lake, output_folder=plot_folder)


