import sys
import os


from ggp_notebooks.plotting_src.header_settings import *

from ggp_notebooks.plotting_src.read_ggp_run import * # you might need to change these two line, such that you can use read ggp_run
from ggp_notebooks.simulations.cells_simulation import * # you might need to change these two line, such that you can use read ggp_run

from ggp_notebooks.plotting_src.cells_utils import *
from ggp_notebooks.plotting_src.read_utils import *


import matplotlib.cm as cm

def return_zero(data):    
    return 0

def cv(data):
    return np.std(data)/np.mean(data)

def strat_cells(cells, code=[1, 2,3]):
    cell_buckets = {}
    for cell in cells:
        bucket_code = '.'.join([cell.cell_id.split('.')[i] for i in code])
        if bucket_code in cell_buckets:
            cell_buckets[bucket_code].append(cell)
        else:
            cell_buckets[bucket_code] = [cell]

    return cell_buckets


def sort_buckets(buckets, order=['bt' ,'pos', 'lane']):
    split_b = [tuple(buck.split('.')) for buck in buckets]

    dtype = [('pos', int), ('bt', "U10"), ('lane', int)]

    a = np.array(split_b, dtype=dtype)
    a = np.sort(a, order=order)  

    return [".".join([str(aii) for aii in ai])  for ai in a]
    
    
def calc_conc(cells):
    for i, _ in enumerate(cells):
        cells[i].concentration = cells[i].gfp/np.exp(cells[i].log_length)
    return cells

def get_conc_from_cells(cells):
    conc = np.array([])
    for cell in cells:
        conc = np.append(conc, cell.concentration)
    return conc

def get_log_length_from_cells(cells):
    ll = np.array([])
    for cell in cells:
        ll = np.append(ll, cell.log_length)
    return ll

def get_length_from_cells(cells):
    ll = np.array([])
    for cell in cells:
        ll = np.append(ll, np.exp(cell.log_length))
    return ll

def get_gfp_from_cells(cells):
    gfp = np.array([])
    for cell in cells:
        gfp = np.append(gfp, cell.gfp)
    return gfp


def calc_cv(arr):
    return np.std(arr)/np.mean(arr)

    
def group_by_position(buckets, code=[0,1]):
    positions = {}
    lanes = {}
    for k in buckets.keys():
        pos = '.'.join([k.split('.')[i] for i in code])
        if pos in positions.keys():
            positions[pos].append(buckets[k][0])
            lanes[pos].append(int(k.split('.')[-1]))
        else:
            positions[pos] = [buckets[k][0]]
            lanes[pos] = [int(k.split('.')[-1])]
    return positions, lanes


def pairwise_difference(a):
    diff = []
    for i in range(len(a)-1):
        for j in range(i+1, len(a)):
            diff.append(a[j] - a[i])
    return np.array(diff)
    
    
        
def get_lane_distance_of_position(lanes):
    distance = {}
    for pos in lanes.keys():
        if "b" in pos:
            distance[pos] = np.array(lanes[pos]) - lanes[pos][0]
        else:
            dist = np.array(lanes[pos]) * (-1)
            distance[pos] = dist - dist[0] 
    return distance
    
from scipy.optimize import curve_fit  
def lin_func(x, a, b):
    return a*x + b
    
    
def gfp_per_biomass(cells):
    gfp = np.array([])
    biomass = np.array([])
    for cell in cells:
        gfp = np.append(gfp, cell.gfp)
        biomass = np.append(biomass, np.exp(cell.log_length))
    return [np.sum(gfp)/np.sum(biomass)]

    

##### plot functions
def growth_lane_stat(prediction_files, 
                     raw_input_dir, 
                     get_func=get_conc_from_cells,
                     filter_func=None, 
                     stat_funcs=[np.mean, np.std], 
                     ylabel=None,
                     cell_id_col=None,
                     cell_rank_col=None,
                     title=None):

    # plot layout
    n_files = len(prediction_files)
    fig, axes = plt.subplots(1, n_files, figsize=(n_files*6,4))
    np.ravel(axes)[0].set_ylabel(ylabel)   
    
    for a in np.ravel(axes):
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

        a.spines['right'].set_color('none')
        a.yaxis.tick_left()

        a.spines['top'].set_color('none')
        a.xaxis.tick_bottom()
    
    for ax, prediction_file in zip(np.ravel(axes), prediction_files):
        prediction = pd.read_csv(prediction_file, 
                                 skiprows=header_lines(prediction_file), 
                                 dtype={"parent_id":str, "cell_id":str})
        cells = ggp_df2cells(prediction)
        raw_file = fetch_raw_file(prediction_file, raw_input_dir)
        
        raw_data_pd = pd.read_csv(raw_file)
        if cell_id_col!=None and cell_rank_col!=None:
            cells = assign_cell_rank_to_cells(cells, raw_data_pd, cell_id_col, cell_rank_col)
        
        if filter_func!=None:
            cells = filter_func(cells)

        cells = calc_concentration(cells)
        cell_buckets = strat_cells(cells)
        bucket_keys = sort_buckets(cell_buckets.keys())

        bucket_stat = {}
        
        for k in bucket_keys:
            quantity_arr = get_func(cell_buckets[k])
            bucket_stat[k] = [stat_func(quantity_arr) for stat_func in stat_funcs]
                

        ax.set_title(' '.join(prediction_file.split('/')[-1].split('_')[:2]))    

        ax.set_xlabel("positions")

        color_norm = mcolors.Normalize(vmin=0, vmax=9.9)

        color_i = 0
#         last = bucket_keys[0]
        last = ''
#         x_position = [0]
        x_position = []
        
#         x_ticks = [0]
        x_ticks = []
#         x_ticks_labels = [".".join(last.split('.')[:2])]
        x_ticks_labels = []
        
        for i, k in enumerate(bucket_keys):
            if ".".join(k.split('.')[:2]) != ".".join(last.split('.')[:2]):
                if len(x_position)>0:
                    color_i+=1
                    x_position.append(x_position[-1]+2)
                    ax.axvline(x_position[-1]-1, ls='--' , color = "gray")
                else:
                    x_position = [0]
                x_ticks_labels.append(".".join(k.split('.')[:2]))
            else:
                dx = int(k.split('.')[-1]) - int(last.split('.')[-1]) 
                x_position.append(x_position[-1]+dx)
                x_ticks_labels.append('')
                
            x_ticks.append(x_position[-1])

            color = cm.tab10(color_norm(color_i))
            if "b" in k:
                marker = "o" 
            else:
                marker ="v"

            ax.errorbar(x_position[-1], bucket_stat[k][0], yerr=bucket_stat[k][1], c=color, marker=marker)
            
            last = k

        xlim = ax.get_xlim()
        ax.set_xticks(x_ticks)
        
        ax.set_xticklabels(x_ticks_labels)
        
        ax.set_xticks(np.arange(100))
        ax.set_xlim(xlim)
        
        ax.tick_params(axis='x', rotation=70)
    plt.tight_layout()
    plt.show()
    
    
    


def growth_lane_xy(prediction_files, 
                     raw_input_dir, 
                     get_func_x=get_conc_from_cells,
                     get_func_y=get_conc_from_cells,
                     filter_func=None, 
                     stat_funcs=[np.mean, return_zero], 
                     xlabel=None, 
                     ylabel=None, 
                     title=None):

    # plot layout
    n_files = len(prediction_files)
    fig, axes = plt.subplots(1, n_files, figsize=(n_files*6,4))
    np.ravel(axes)[0].set_ylabel(ylabel)   
    
    for a in np.ravel(axes):
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

        a.spines['right'].set_color('none')
        a.yaxis.tick_left()

        a.spines['top'].set_color('none')
        a.xaxis.tick_bottom()
    
    for ax, prediction_file in zip(np.ravel(axes), prediction_files):
        prediction = pd.read_csv(prediction_file, 
                                 skiprows=header_lines(prediction_file), 
                                 dtype={"parent_id":str, "cell_id":str})
        cells = ggp_df2cells(prediction)
        raw_file = fetch_raw_file(prediction_file, raw_input_dir)
        
        raw_data_pd = pd.read_csv(raw_file)
        cells = assign_cell_rank_to_cells(cells, raw_data_pd)
        
        if filter_func!=None:
            cells = filter_func(cells)

        cells = calc_concentration(cells)
        cell_buckets = strat_cells(cells)
        bucket_keys = sort_buckets(cell_buckets.keys())

        bucket_stat = {}
        bucket_stat2 = {}
        
        for k in bucket_keys:
            quantity_arr = get_func_x(cell_buckets[k])
            bucket_stat[k] = [stat_func(quantity_arr) for stat_func in stat_funcs]
            
            quantity_arr2 = get_func_y(cell_buckets[k])
            bucket_stat2[k] = [stat_func(quantity_arr2) for stat_func in stat_funcs]
            
                
        ax.set_title(' '.join(prediction_file.split('/')[-1].split('_')[:2]))    

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        color_norm = mcolors.Normalize(vmin=0, vmax=9.9)
        color_i = -1
#         last = bucket_keys[0]
        last = ''

        for i, k in enumerate(bucket_keys):
            if ".".join(k.split('.')[:2]) != ".".join(last.split('.')[:2]):
                color_i+=1

            color = cm.tab10(color_norm(color_i))
            if "b" in k:
                marker = "o" 
            else:
                marker ="v"
            ax.errorbar(bucket_stat[k][0], bucket_stat2[k][0], 
                        xerr=bucket_stat[k][1], yerr=bucket_stat2[k][1], 
                        c=color, marker=marker)
            last = k

    plt.tight_layout()
    plt.show()
    
    

def growth_lane_stat_pairs(prediction_files, 
                     raw_input_dir, 
                     get_func=get_conc_from_cells,
                     filter_func=None, 
                     stat_funcs=[np.mean, np.std],
                    xlabel=None,
                     ylabel=None, 
                     title=None):

    # plot layout
    n_files = len(prediction_files)
    fig, axes = plt.subplots(1, n_files, figsize=(n_files*6,5))
    np.ravel(axes)[0].set_ylabel(ylabel)   
    
    for a in np.ravel(axes):
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

        a.spines['right'].set_color('none')
        a.yaxis.tick_left()

        a.spines['top'].set_color('none')
        a.xaxis.tick_bottom()
    
    for ax, prediction_file in zip(np.ravel(axes), prediction_files):
        prediction = pd.read_csv(prediction_file, 
                                 skiprows=header_lines(prediction_file), 
                                 dtype={"parent_id":str, "cell_id":str})
        cells = ggp_df2cells(prediction)
        raw_file = fetch_raw_file(prediction_file, raw_input_dir)
        
        raw_data_pd = pd.read_csv(raw_file)
        cells = assign_cell_rank_to_cells(cells, raw_data_pd)
        
        if filter_func!=None:
            cells = filter_func(cells)

        cells = calc_concentration(cells)
        cell_buckets = strat_cells(cells)
        bucket_keys = sort_buckets(cell_buckets.keys())

        bucket_stat = {}
        
        for k in bucket_keys:
            quantity_arr = get_func(cell_buckets[k])            
            bucket_stat[k] = [stat_func(quantity_arr) for stat_func in stat_funcs]
                
        positions, lanes = group_by_position(bucket_stat)


        ax.set_title(' '.join(prediction_file.split('/')[-1].split('_')[:2]))    
        ax.set_xlabel(xlabel)

        color_norm = mcolors.Normalize(vmin=0, vmax=9.9)

        color_i = 0

        
        dx_coll = np.array([])
        dstat_coll = np.array([])
        for i, k in enumerate(positions.keys()):
            dx = pairwise_difference(lanes[k])        
            dstat = pairwise_difference(positions[k])
            
            if "b" in k:
                marker = "o" 
            else:
                marker ="v"
                dx *= -1
            
            color = cm.tab10(color_norm(i))
            if len(dx)>0:
                ax.scatter(dx, dstat, c=[color], marker=marker, label=k)
                dx_coll = np.append(dx_coll, dx)
                dstat_coll = np.append(dstat_coll, dstat)
        
        ax.legend(ncol=3, bbox_to_anchor=(0,1.2), loc="lower left")

        corr = np.corrcoef(dx_coll, dstat_coll)[0,1]
        ax.set_title(' '.join(prediction_file.split('/')[-1].split('_')[:2]) + ", C=" +str(np.around(corr, decimals=2))) 
        
        
        
    plt.tight_layout()
    plt.show()
    



def growth_lane_regression(prediction_files, 
                            raw_input_dir, 
                            get_func=get_conc_from_cells,
                            filter_func=None, 
                            stat_funcs=[np.mean, np.std],
                            xlabel=None,
                            ylabel=None, 
                            title=None):

    # plot layout
    n_files = len(prediction_files)
    fig, axes = plt.subplots(int(n_files/4), 4, figsize=(4*6,int(n_files/4)*5))
#     plt.figure(facecolor='white')
    np.ravel(axes)[0].set_ylabel(ylabel)   
    
    for a in np.ravel(axes):
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

        a.spines['right'].set_color('none')
        a.yaxis.tick_left()

        a.spines['top'].set_color('none')
        a.xaxis.tick_bottom()
    
    for ax, prediction_file in zip(np.ravel(axes), prediction_files):
        prediction = pd.read_csv(prediction_file, 
                                 skiprows=header_lines(prediction_file), 
                                 dtype={"parent_id":str, "cell_id":str})
        cells = ggp_df2cells(prediction)
        raw_file = fetch_raw_file(prediction_file, raw_input_dir)
        
        raw_data_pd = pd.read_csv(raw_file)
        cells = assign_cell_rank_to_cells(cells, raw_data_pd)
        
        if filter_func!=None:
            cells = filter_func(cells)

        cells = calc_concentration(cells)
        cell_buckets = strat_cells(cells)
        bucket_keys = sort_buckets(cell_buckets.keys())
        bucket_stat = {}
        
        for k in bucket_keys:
            quantity_arr = get_func(cell_buckets[k])  
            if len(quantity_arr)>1:
                bucket_stat[k] = [np.mean(quantity_arr)]
            else:
                bucket_stat[k] = quantity_arr
                
        positions, lanes = group_by_position(bucket_stat)
        lanes = get_lane_distance_of_position(lanes)

        ax.set_title(' '.join(prediction_file.split('/')[-1].split('_')[:2]))    
        ax.set_xlabel(xlabel)

        color_norm = mcolors.Normalize(vmin=0, vmax=9.9)
        color_i = 0   
        ### PLOTS
        for i, k in enumerate(positions.keys()):
            dx = lanes[k]        
            dstat = positions[k]
            
            if "b" in k:
                marker = "o" 
            else:
                marker ="v"
            if len(dx)>2:
                color = cm.tab10(color_norm(i))
                ax.scatter(dx, dstat, c=[color], marker=marker, label=k)
                popt, pcov = curve_fit(lin_func, dx, dstat)
                ax.plot(dx, lin_func(dx, *popt), color=color)        
                ax.fill_between(dx, lin_func(dx, popt[0]+np.sqrt(pcov[0,0]), popt[1]), 
                                    lin_func(dx, popt[0]-np.sqrt(pcov[0,0]), popt[1]), color=color, alpha=0.15)        
        ### Final layout 
        ax.axvline(0, color="gray")
        ax.legend(ncol=3, bbox_to_anchor=(0,1.2), loc="lower left")
        ax.set_title(' '.join(prediction_file.split('/')[-1].split('_')[:2]) ) 
        
    plt.tight_layout()
    fig.savefig('growth_lane_regression.pdf', facecolor="white", edgecolor='none')
    plt.show()
    
    
  