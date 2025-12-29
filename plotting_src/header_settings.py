import sys
import os

import pandas as pd
import copy

import numpy as np 
import scipy.stats
from scipy.optimize import curve_fit
from scipy import special

from inspect import Parameter

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
from matplotlib.legend_handler import HandlerTuple
import matplotlib.colors as mcolors
import matplotlib.cm as cmcolors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

import matplotlib.lines as mlines

TINY_SIZE = 10

SMALL_SIZE = 12
MEDIUM_SIZE = 16

FULL_SIZE = 13 # in inches
HALF_SIZE = FULL_SIZE/2 
THIRD_SIZE = FULL_SIZE/3

FULL_HEIGHT = 4.
HALF_HEIGHT = FULL_HEIGHT/2.
THIRD_HEIGHT = FULL_HEIGHT/3.

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

import builtins
if hasattr(builtins, "NAMING"):
    NAMING = builtins.NAMING

else:
    NAMING = "old"
print("header_setting.py: using", NAMING, "naming")


colors  = ['tab:blue',
           'tab:orange',
           'tab:green', 
           'tab:red', 
           'tab:purple', 
           'tab:brown', 
           'tab:pink', 
           'tab:olive']

single_color_cmaps = ['Blues', "Oranges", "Greens", "Reds"]

# colors = [  '#58341cff',
#             '#3a6a5eff',
#             '#526586ff',
#             '#ae8859ff',
#             '#a7aebaff',
#             '#d69f96ff']


markers = ['o', '^', 's',  'p', '<', '>']

    
if NAMING == "old":
    ## CONDITIONS
    conditions = ["acetate", "glycerol", "glucose", "glucoseaa"]
    marker_by_condition = {cond: m for cond,m  in zip(conditions, markers)}
    color_by_condition = {cond: col for cond,col  in zip(conditions, colors)}
    cmap_by_condition = {cond: col for cond,col  in zip(conditions, single_color_cmaps)}

    ## PROMOTER
    promoters = ['hi1', 'hi3', 'med2', 'med3', 'rrnB', 'rpsB', 'rpmB', 'rplN']
    ribosomal = ['rrnB', 'rpsB', 'rpmB', 'rplN']
    synthetic = ['hi1', 'hi3', 'med2', 'med3']

    color_by_promoter = {p:c for p,c in zip(promoters, colors)}


elif NAMING == "new":
    ## CONDITIONS
    conditions = ["acetate005", "glycerol040", "glucose020", "glucoseaa020"]
    conditions_labels = ["acetate", "glycerol", "glucose", "gluc.+aa"]
    label_by_condition = {c:l for l,c in zip(conditions_labels, conditions)}
    marker_by_condition = {cond: m for cond,m  in zip(conditions, markers)}
    # color_by_condition = {cond: col for cond,col  in zip(conditions, colors)}
    cmap_by_condition = {cond: col for cond,col  in zip(conditions, single_color_cmaps)}
    beta_by_condition = {"acetate005": 0.0008, "glycerol040": 0.0016, "glucose020": 0.0032, "glucoseaa020":0.0064}

    ## PROMOTER
    promoters = ['hi1', 'hi3', 'med2', 'med3', 'rplN', 'rrnB']
    ribosomal = promoters[-2:]
    synthetic = ['hi1', 'hi3', 'med2', 'med3']
    
    import matplotlib.colors

    def categorical_cmap(nc, nsc, cmap="tab10", continuous=False):
        if nc > plt.get_cmap(cmap).N:
            raise ValueError("Too many categories for colormap.")
        if continuous:
            ccolors = plt.get_cmap(cmap)(np.linspace(0,1,nc))
        else:
            ccolors = plt.get_cmap(cmap)(np.arange(nc, dtype=int))
        cols = np.zeros((nc*nsc, 3))
        for i, c in enumerate(ccolors):
            chsv = matplotlib.colors.rgb_to_hsv(c[:3])
            arhsv = np.tile(chsv,nsc).reshape(nsc,3)
            arhsv[:,1] = np.linspace(chsv[1],0.25,nsc)
            arhsv[:,2] = np.linspace(chsv[2],1,nsc)
            rgb = matplotlib.colors.hsv_to_rgb(arhsv)
            cols[i*nsc:(i+1)*nsc,:] = rgb       
        cmap = matplotlib.colors.ListedColormap(cols)
        return cmap, cols

    def rgb_array_to_hex(rgbs):
        hexs = []
        for rgb in rgbs:
            hexs.append(mpl.colors.rgb2hex(rgb))
        return hexs

    _, cols = categorical_cmap(10, 4, cmap="tab10")
    hexs = rgb_array_to_hex(cols)

    # mix_tab = (tab20c[:16]+tab20b[8:16])[::-1] ## old version
    mix_tab = (hexs[:20]+hexs[24:28])
    mix_tab = np.ravel([mix_tab[i*4: (i+1)*4][::1] for i in range(6)])


    color_by_condition = {cond: col for cond,col  in zip(conditions, mix_tab)}

    # color_by_condition = {cond: col for cond,col  in zip(conditions, colors_conditions)}
    color_by_promoter = {p:c for p,c in zip(promoters, mix_tab[0::4])}
    color_by_both = {p: {cond: col for cond,col in zip(conditions, mix_tab[i*4:])} for i,p in enumerate(promoters)}
    grey_by_condition = {c:g for c,g in zip(conditions, ["#636363" ,"#969696" ,"#bdbdbd" ,"#d9d9d9"][::1])}
    ls_by_condition = {c:g for c,g in zip(conditions, ["-" ,"--" ,":" ,"-."][::1])}

    # for i, m in enumerate(hexs):
    #     plt.scatter(i,0,c=m)

    # for i, m in enumerate(mix_tab):
    #     plt.scatter(i,1,c=m)

    # for i, m in enumerate(grey_by_condition.values()):
    #     plt.scatter(i,2,c=m)



    marker_by_promoter = {p: m for p,m  in zip(promoters, markers)}
    # color_by_promoter = {p:c for p,c in zip(promoters, colors)}

    dates = ['20190515',
            '20190529',
            '20220921',
            '20221111',
            '20230223',
            '20230228',
            '20230302',
            '20230413',
            '20230420',
            '20230428']
    years = ['2019', '2022', '2023']
    marker_by_year = {y: m for y,m  in zip(years, markers)}

    marker_by_replicate = {"1": "o", "2":"^"}





LATEX_LABEL = {'mean_lambda': r'$\bar{\lambda}$',
            'gamma_lambda': r'$\gamma_\lambda$',
            'var_lambda': r'$\sigma^2_\lambda$',
            'mean_q': r'$\bar{q}$',
            'gamma_q': r'$\gamma_q$',
            'var_q': r'$\sigma^2_q$',
            'beta': r'$\beta$',
            'var_x': r'$\sigma^2_x$',
            'var_g': r'$\sigma^2_g$',
            'var_dx': r'$\sigma^2_{dx}$',
            'var_dg': r'$\sigma^2_{dg}$'} 

parameter_names = {'mean_lambda': r'$\bar{\lambda}\ \mathrm{(1/min)}$',
                    'gamma_lambda': r'$\gamma_\lambda\ \mathrm{(1/min)}$',
                    'var_lambda': r'$\sigma_\lambda^2\ \mathrm{(1/min^2)}$',
                    'mean_q': r'$\bar{q}\ \mathrm{(1/min)}$',
                    'gamma_q': r'$\gamma_q\ \mathrm{(1/min)}$',
                    'var_q': r'$\sigma_q^2\ \mathrm{(1/min^2)}$',
                    'beta': r'$\beta\ \mathrm{(1/min)}$',
                    'var_x': r'$\sigma_x^2\ \mathrm{(\mu m^2)}$',
                    'var_g':r'$\sigma_g^2\ \mathrm{(a.u.)}$',
                    'var_dx': r'$\sigma_{dx}^2$',
                    'var_dg':r'$\sigma_{dg}^2$',
                    'std_lambda':r'$\sqrt{\frac{\sigma_\lambda^2}{2\gamma_\lambda}}\ \mathrm{(1/min)}$',
                    'std_q':r'$\sqrt{\frac{\sigma_q^2}{2\gamma_q}}\ \mathrm{(1/min)}$', 
                    'sigma_x': r'$\sigma_x\ \mathrm{\mu m}$',
                    'sigma_g': r'$\sigma_g\ \mathrm{a.u.}$'}
                    

parameter_names_unit_free = {'mean_lambda': r'$\bar{\lambda}}$',
                            'gamma_lambda': r'$\gamma_\lambda}$',
                            'tau_lambda': r'$1/\gamma_\lambda$',
                            'var_lambda': r'$\sigma_\lambda^2}$',
                            'sigma_lambda': r'$\sigma_\lambda}$',
                            'cv_lambda': r'$\frac{\sqrt{\sigma_\lambda^2 /(2 \gamma_\lambda)}}{\bar{\lambda}}$',
                            'mean_q': r'$\bar{q}}$',
                            'gamma_q': r'$\gamma_q}$',
                            'tau_q': r'$1/\gamma_q$',
                            'var_q': r'$\sigma_q^2}$',
                            'sigma_q': r'$\sigma_q}$',
                            'cv_q': r'$\frac{\sqrt{\sigma_q^2 /(2 \gamma_q)}}{\bar{q}}$',
                            'beta': r'$\beta}$',
                            'var_x': r'$\sigma_x^2}$',
                            'var_g':r'$\sigma_g^2}$',
                            'var_dx': r'$\sigma_{dx}^2$',
                            'var_dg':r'$\sigma_{dg}^2$',
                            'sigma_dx': r'$\sigma_{dx}$',
                            'sigma_dg':r'$\sigma_{dg}$',
                            # 'std_lambda':r'$\sqrt{\frac{\sigma_\lambda^2}{2\gamma_\lambda}}$',
                            'std_lambda':r's.d.($\lambda$)',
                            # 'std_q':r'$\sqrt{\frac{\sigma_q^2}{2\gamma_q}}$',
                            'std_q':r's.d.($q$)',
                            'sigma_x': r'$\sigma_x$',
                            'sigma_g': r'$\sigma_g$'}

