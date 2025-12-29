
import copy
import numpy as np

from RealTrace_paper_analysis.plotting_src.header_settings import *
from RealTrace_paper_analysis.plotting_src.cells_utils import *


def deep_flatten(arr):
    new = np.array([], dtype=object)
    for a in arr:
        new = np.append(new, np.array(a))
    return np.array(new, dtype=float)

def index_nested(arr, idx, flat=False):
    new = []
    for i in idx:
        new.append(arr[i])
    if flat:
        return deep_flatten(new)
    return new

def distribute_over_cells(cells, time, rate):
    t =[]
    r =[]
    for cell in cells:
        t0 = time>=cell.time[0]-1e-5
        t1 = time<=cell.time[-1]+1e-5
        t.append(time[t0*t1])
        r.append(rate[t0*t1])
    return t, r


def sliding_mean(x, n=1):
    if n==1:
        return x
    idx = np.arange(n) + np.arange(len(x)-n+1)[:,None]
    mean = np.nanmean(x[idx],axis=1)
    return mean

def naive_growth_rate(cells, window_size=1, beta=None):
    new_cells = copy.deepcopy(cells)
    time = np.array([])
    log_length = np.array([])
    for idx, _ in enumerate(cells):
        time = np.append(time, cells[idx].time)
        log_length = np.append(log_length, cells[idx].log_length + idx*np.log(2)) 

    growth_rate =  np.diff(log_length)/np.diff(time)

    time_avg = sliding_mean((time[:-1]+np.diff(time)/2), window_size)
    rate_avg = sliding_mean(growth_rate,  window_size)
    t, r = distribute_over_cells(cells, time_avg,  rate_avg)
    for idx, _ in enumerate(cells):
        new_cells[idx].naive_rate = r[idx]
        new_cells[idx].time_naive_rate = t[idx]
    return new_cells


def naive_growth_rate_skip_div(cells, window_size=1, beta=None):
    new_cells = copy.deepcopy(cells)
    time = []
    growth_rate = []
    
    for idx, _ in enumerate(cells):
        temp = cells[idx].time[:-1]+np.diff(cells[idx].time)/2
        # gr = np.append(np.diff(cells[idx].log_length)/np.diff(cells[idx].time), np.nan)
        
        if idx < len(cells)-1:
            t = np.append(temp, cells[idx].time[-1] + np.diff(cells[idx].time)[-1]/2)
            gr_last = (cells[idx+1].log_length[0]+np.log(2) - cells[idx].log_length[-1])/(cells[idx+1].time[0] - cells[idx].time[-1])
            gr = np.append(np.diff(cells[idx].log_length)/np.diff(cells[idx].time), gr_last)
        else:
            t = temp
            gr = np.diff(cells[idx].log_length)/np.diff(cells[idx].time)
            
        time.append(t)
        growth_rate.append(gr)

    time_avg = sliding_mean(deep_flatten(time), window_size)
    rate_avg = sliding_mean(deep_flatten(growth_rate),  window_size)
    
    t, r = distribute_over_cells(cells, time_avg,  rate_avg)
    for idx, _ in enumerate(cells):
        new_cells[idx].naive_rate = r[idx]
        new_cells[idx].time_naive_rate = t[idx]
    return new_cells


def naive_production_rate(cells, window_size=1, beta=0):
    new_cells = copy.deepcopy(cells)
    time = np.array([])
    gfp = np.array([])
    log_length = np.array([])

    for idx, _ in enumerate(cells):
        time = np.append(time, cells[idx].time)
        gfp = np.append(gfp, np.array(cells[idx].gfp) * 2**idx) 
        log_length = np.append(log_length, cells[idx].log_length + idx*np.log(2)) 

    log_length_interpl = np.exp(log_length)[:-1] + np.diff(np.exp(log_length))/.2
    gfp_interpl = gfp[:-1] + np.diff(gfp)/2.
    production_rate = np.diff(gfp)/np.diff(time)/log_length_interpl + gfp_interpl*beta/log_length_interpl

    time_avg = sliding_mean((time[:-1]+np.diff(time)/2), window_size)
    rate_avg = sliding_mean(production_rate,  window_size)
    t, r = distribute_over_cells(cells, time_avg, rate_avg)
    for idx, _ in enumerate(cells):
        new_cells[idx].naive_rate = r[idx]
        new_cells[idx].time_naive_rate = t[idx]
    return new_cells


def naive_production_rate_skip_div(cells, window_size=1, beta=0):
    new_cells = copy.deepcopy(cells)
    
    time = []
    production_rate = []

    for idx, _ in enumerate(cells):
        temp = cells[idx].time[:-1]+np.diff(cells[idx].time)/2
        log_length_interpl = np.exp(cells[idx].log_length[:-1] + np.diff(cells[idx].log_length)/2.) 
        gfp_interpl = cells[idx].gfp[:-1] + np.diff(cells[idx].gfp)/2.

        if idx < len(cells)-1:
            log_length_interpl_last = np.exp(cells[idx].log_length[-1] + (cells[idx+1].log_length[0]-cells[idx].log_length[-1])/2.) 
            pr_last = (cells[idx+1].gfp[0]*2-cells[idx].gfp[-1])/(cells[idx+1].time[0]- cells[idx].time[-1])/log_length_interpl_last

            t = np.append(temp, cells[idx].time[-1] + np.diff(cells[idx].time)[-1]/2)
            pr = np.append(np.diff(cells[idx].gfp)/np.diff(cells[idx].time)/log_length_interpl + gfp_interpl*beta/log_length_interpl, pr_last)
        else:
            t = temp
            pr = np.diff(cells[idx].gfp)/np.diff(cells[idx].time)/log_length_interpl + gfp_interpl*beta/log_length_interpl
            
        time.append(t)
        production_rate.append(pr)

#     print(len(production_rate[0]), len(cells[0].time))
    time_avg = sliding_mean(deep_flatten(time), window_size)
    rate_avg = sliding_mean(deep_flatten(production_rate),  window_size)
    t, r = distribute_over_cells(cells, time_avg,  rate_avg)
        
    for idx, _ in enumerate(cells):
        new_cells[idx].naive_rate = r[idx]
        new_cells[idx].time_naive_rate = t[idx]
    return new_cells


def naive_rate_for_all_cells(cells, window_size, func, beta=0):
    new_cells = []
    for cell in cells:
        _, parent_cell = get_parent(cell, cells)
        daughter_cells = get_daughters(cell, cells)

        cell_conc = [cell]
        take_cell = 0
        if parent_cell !=None:
            cell_conc = [parent_cell] + cell_conc
            take_cell = 1

        if daughter_cells[0] != None and daughter_cells[1] != None:
            cella = func(cell_conc + [daughter_cells[0]], window_size, beta)[take_cell]
            cellb = func(cell_conc + [daughter_cells[1]], window_size, beta)[take_cell]
            cell.naive_rate = (cella.naive_rate + cellb.naive_rate)/2
            cell.time_naive_rate = cella.time_naive_rate

        elif daughter_cells[0] != None:
            cell = func(cell_conc + [daughter_cells[0]], window_size, beta)[take_cell]

        else:
            cell = func(cell_conc, window_size, beta)[take_cell]
        new_cells.append(cell)
    return new_cells

###########


def sliding_growth_rate(cells, window_size=1):
    new_cells = copy.deepcopy(cells)
    time = np.array([])
    growth_rate = np.array([])

    for idx, _ in enumerate(cells):
        time = np.append(time, cells[idx].time)
        growth_rate = np.append(growth_rate, cells[idx].lt) 

    time_avg = sliding_mean(time, window_size)
    rate_avg = sliding_mean(growth_rate,  window_size)
    t, r = distribute_over_cells(cells, time_avg,  rate_avg)
    for idx, _ in enumerate(cells):
        new_cells[idx].sliding_rate = r[idx]
        new_cells[idx].time_sliding_rate = t[idx]
    return new_cells

def sliding_production_rate(cells, window_size=1):
    new_cells = copy.deepcopy(cells)
    time = np.array([])
    production_rate = np.array([])

    for idx, _ in enumerate(cells):
        time = np.append(time, cells[idx].time)
        production_rate = np.append(production_rate, np.array(cells[idx].qt)) 

    time_avg = sliding_mean(time, window_size)
    rate_avg = sliding_mean(production_rate,  window_size)
    t, r = distribute_over_cells(cells, time_avg,  rate_avg)
    for idx, _ in enumerate(cells):
        new_cells[idx].sliding_rate = r[idx]
        new_cells[idx].time_sliding_rate = t[idx]
    return new_cells


def sliding_rate_for_all_cells(cells, window_size, func=sliding_growth_rate):
    new_cells = []
    for cell in cells:
        _, parent_cell = get_parent(cell, cells)
        daughter_cells = get_daughters(cell, cells)

        cell_conc = [cell]
        take_cell = 0
        if parent_cell !=None:
            cell_conc = [parent_cell] + cell_conc
            take_cell = 1

        if daughter_cells[0] != None and daughter_cells[1] != None:
            cella = func(cell_conc + [daughter_cells[0]], window_size)[take_cell]
            cellb = func(cell_conc + [daughter_cells[1]], window_size)[take_cell]
            cell.sliding_rate = (cella.sliding_rate + cellb.sliding_rate)/2
            cell.time_sliding_rate = cella.time_sliding_rate

        elif daughter_cells[0] != None:
            cell = func(cell_conc + [daughter_cells[0]], window_size)[take_cell]

        else:
            cell = func(cell_conc, window_size)[take_cell]
        new_cells.append(cell)
    return new_cells

