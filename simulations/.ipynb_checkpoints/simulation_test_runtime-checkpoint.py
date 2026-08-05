import sys

import sys
from cells_simulation import *

import sys
sys.path.insert(0, "../../")

from RealTrace_paper_analysis.plotting_src.read_ggp_run import * 
from RealTrace_paper_analysis.plotting_src.naive_rates import * 

import numpy as np
import subprocess
from multiprocessing import Pool, get_context
import copy
import time

def add_noise_in_log(x, scale):
    return x * np.random.lognormal(0, sigma=scale)
    
def write_param_file_var_noise(filename, 
                               parameters, 
                               rel_step, 
                               mean_lambda_init,
                               mean_q_init):
    with open(filename, "w") as fin:
        
        parameters_init = {}
        n_level2 = 1
        fin.write("# Generated config file for simulated data with noise\n")
        
        range_lim=20
        
        mean_lambda =  mean_lambda_init
        fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format("mean_lambda", 
                                                   mean_lambda,
                                                   mean_lambda*0.001,
                                                   mean_lambda/10,
                                                   mean_lambda*10))
        
        gamma_lambda = add_noise_in_log(parameters["gamma_lambda"], n_level2)
        fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format("gamma_lambda", 
                                                   gamma_lambda,
                                                   gamma_lambda*rel_step,
                                                   gamma_lambda/range_lim,
                                                   gamma_lambda*range_lim))
        
        var_lambda = add_noise_in_log(parameters["var_lambda"], n_level2)
        fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format("var_lambda", 
                                                   var_lambda,
                                                   var_lambda*rel_step,
                                                   var_lambda/range_lim,
                                                   var_lambda*range_lim))
        
        mean_q = mean_q_init
        fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format("mean_q", 
                                                   mean_q,
                                                   mean_q*0.001,
                                                   mean_q/10,
                                                   mean_q*10))
        
        gamma_q = add_noise_in_log(parameters["gamma_q"], n_level2)
        fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format("gamma_q", 
                                                   gamma_q,
                                                   gamma_q*rel_step,
                                                   gamma_q/range_lim,
                                                   gamma_q*range_lim))
        
        fin.write("{:s} = {:.2E}\n".format("beta", parameters["beta"]))
        
        var_q = add_noise_in_log(parameters["var_q"], n_level2)
        fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format("var_q", 
                                                   var_q,
                                                   var_q*rel_step,
                                                   var_q/range_lim,
                                                   var_q*range_lim))
        
        var_x = add_noise_in_log(parameters["var_x"], n_level2)
        fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format("var_x", 
                                                   var_x,
                                                   var_x*rel_step,
                                                   var_x/range_lim,
                                                   var_x*range_lim))
        
        var_g = add_noise_in_log(parameters["var_g"], n_level2)
        fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format("var_g", 
                                                   var_g,
                                                   var_g*rel_step,
                                                   var_g/range_lim,
                                                   var_g*range_lim))
        
        var_dx = add_noise_in_log(parameters["var_dx"], n_level2)
        fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format("var_dx", 
                                                   var_dx,
                                                   var_dx*rel_step,
                                                   var_dx/range_lim,
                                                   var_dx*range_lim))
        
        var_dg = add_noise_in_log(parameters["var_dg"], n_level2)
        fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format("var_dg", 
                                                   var_dg,
                                                   var_dg*rel_step,
                                                   var_dg/range_lim,
                                                   var_dg*range_lim))
        


    return filename


# ========== Simulation parameters ========== #
save_dataset = True    # saves data set as csv
run_ggp = True         # runs the ggp code (located in "../bin/" relative to this nb)

dt = 1e-3
dt_measument = 1 
discard_cells = 5
lanes = 10

div_mode = "adder"

division_log_length = 1+np.log(2)   # for sizer: division, when log_length hits division_log_length
division_time = 1        # for timer: division, when cell cycle time hits division_time
division_addition = 1    # for adder: divsion, when division_addition in log_length was added in cell cycle

n=10000
base_parameter_set  ={ "mean_lambda": 0.0114137, 
                    "gamma_lambda": 0.166994,
                    "var_lambda": 4.9481e-06,
                    "mean_q": 118.524,
                    "gamma_q": 0.0595599,
                    "var_q": 430.055,
                    "beta": 0, #0.0032,
                    "var_x": 0.000136701,
                    "var_g": 12.255,
                    "var_dx": 0.00468653,
                    "var_dg": 1.53195}


def main_sim(sim_idx):
    np.random.seed(0)
    out_dir = mk_missing_dir('../../data_realtrace_paper/test_runtime')

    # np.array([1.5**n for n in np.arange(-5,5)])

    n_cells = [10, 50, 100, 500, 1000][sim_idx]
    # print(n_cells)
    
    parameter_set = copy.deepcopy(base_parameter_set)
    
    # dataset = pd.DataFrame()
    # for i in np.arange(lanes):
    #     cells_simulated = simulate_cells(dt, n_cells, parameter_set, div_mode,
    #                                 division_log_length, 
    #                                 division_time, 
    #                                 division_addition, 
    #                                 tree=True, 
    #                                 discard_cells=discard_cells)

    #     gfp_smaller_zero = False
    #     for cell in cells_simulated:
    #         if np.any(np.array(cell.gfp)<0):
    #             gfp_smaller_zero=True

    #     if gfp_smaller_zero:
    #         continue

    #     temp_dataset = build_data_set_scale_gfp_noise(cells_simulated, 
    #                                             parameter_set['var_x'], 
    #                                             parameter_set['var_g'], 
    #                                             dt_measument,
    #                                             atol=dt/10)
    #     temp_dataset['lane'] = i+1
    #     # add the last "lane" to the data frame
    #     dataset = pd.concat([dataset, temp_dataset], ignore_index=True)
    
    # print("SAVE", sim_idx)
    # ----------- SAVE ----------- #
    directory, filename = get_next_file_name(out_dir,force_index=sim_idx)
    # # write_param_file(os.path.join(directory, "parameters_simulation.txt"), parameter_set) 

    # mean_lambda_init = np.mean(np.concatenate([cell.naive_rate for cell in naive_rate_for_all_cells(df2cells(dataset), 1, naive_growth_rate_skip_div)]))
    
    # mean_q_init = np.mean(np.concatenate([cell.naive_rate for cell in naive_rate_for_all_cells(df2cells(dataset), 1, naive_production_rate_skip_div,beta=parameter_set['beta'])]))
    
    param_file_temp = "parameters_noise{:d}.txt".format(0)
    # param_file = write_param_file_var_noise(os.path.join(directory, param_file_temp), 
    #                                         parameter_set, 
    #                                         0.1, 
    #                                         mean_lambda_init,
    #                                         mean_q_init)
    param_file = os.path.join(directory, param_file_temp)

    out_file = param_file[:-4] + "_out"

    csv_config_file = write_csv_config(os.path.join(directory, "csv_config.txt"), lane="lane")
    # dataset.to_csv(filename)

    print(time.process_time())

    start = time.perf_counter()
    # ----------- RUN ----------- #
    run_cmd = "./../../RealTrace/bin/RealTrace" + \
                " -c " + csv_config_file + \
                " -b " + param_file + \
                " -t " + str(1e-40) + \
                " -i " + filename + \
                " -space log -m -p -j -div binomial -noise scaled" + \
                " -o " + out_file
    s = subprocess.run(run_cmd, capture_output=True, shell=True) 
    elapsed_time = time.perf_counter() - start
    print("\n", n_cells, elapsed_time, "\n")
        
        
        
def main():
    for i in range(5):
        main_sim(i)
    # n_cores = os.environ['SLURM_JOB_CPUS_PER_NODE']
    # print("slurm: number of cpus {}".format(n_cores))
    # print(n_cores)
    # with get_context("spawn").Pool(int(n_cores)) as p:
    #     print("Start multiprocessing")
    #     p.map(main_sim, np.arange(10))
    #     p.close()
    #     p.join()
    #     print("Join")
        
if __name__ == "__main__":
    main()
   