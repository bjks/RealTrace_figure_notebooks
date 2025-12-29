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
n_cells = 100
discard_cells = 5
lanes = 10

div_mode = "adder"

division_log_length = 1+np.log(2)   # for sizer: division, when log_length hits division_log_length
division_time = 1        # for timer: division, when cell cycle time hits division_time
division_addition = 1    # for adder: divsion, when division_addition in log_length was added in cell cycle

n=10000
parameter_ranges = {"mean_lambda":  np.geomspace(0.005, 0.05, n), 
                    "gamma_lambda": np.geomspace(0.05,  0.5, n), 
                    # "var_lambda":   np.geomspace(1e-8,  1e-4, n), 
                    "mean_q":       np.geomspace(100,     1000, n), 
                    "gamma_q":      np.geomspace(0.05,  0.5, n), 
                    # "var_q":        np.geomspace(0.01,  1, n),
                    "beta":         np.geomspace(1e-4,  1e-3, n),
                    "var_x":        np.geomspace(2e-5,  2e-4, n),
                    "var_g":        np.geomspace(1.5,15 , n),
                    "var_dx":       np.geomspace(1e-4,  1e-2, n),
                    "var_dg":       np.geomspace(1,  100, n)}

cv_range = np.linspace(0.2, 0.8, n)


def main_sim(seed):
    np.random.seed(seed)
    out_dir = mk_missing_dir('../../data_sim_manuscript/random_parameters_simulations')
    # ========== Model parameters ========== #
    parameter_set = {   "mean_lambda":np.random.choice(parameter_ranges["mean_lambda"]) ,        # 0.01, 
                        "gamma_lambda": np.random.choice(parameter_ranges["gamma_lambda"]),       # 0.01,
                        # "var_lambda": np.random.choice(parameter_ranges["var_lambda"]),         # 1e-7,
                        "mean_q": np.random.choice(parameter_ranges["mean_q"]),               # 10,
                        "gamma_q": np.random.choice(parameter_ranges["gamma_q"]),            # 0.01,
                        # "var_q": np.random.choice(parameter_ranges["var_q"]),               # 0.1,
                        "beta": np.random.choice(parameter_ranges["beta"]),               # 0.001, 
                        "var_x": np.random.choice(parameter_ranges["var_x"]),              # 1e-3
                        "var_g": np.random.choice(parameter_ranges["var_g"]),              # 5000
                        "var_dx": np.random.choice(parameter_ranges["var_dx"]),
                        "var_dg": np.random.choice(parameter_ranges["var_dg"])}

    parameter_set["var_lambda"] = 2.* parameter_set["gamma_lambda"] * (np.random.choice(cv_range) * parameter_set["mean_lambda"])**2 
    parameter_set["var_q"] = 2.* parameter_set["gamma_q"] * (np.random.choice(cv_range) * parameter_set["mean_q"])**2 

    dataset = pd.DataFrame()
    for i in np.arange(lanes):
        cells_simulated = simulate_cells(dt, n_cells, parameter_set, div_mode,
                                    division_log_length, 
                                    division_time, 
                                    division_addition, 
                                    tree=True, 
                                    discard_cells=discard_cells)

        gfp_smaller_zero = False
        for cell in cells_simulated:
            if np.any(np.array(cell.gfp)<0):
                gfp_smaller_zero=True

        if gfp_smaller_zero:
            continue

        temp_dataset = build_data_set_scale_gfp_noise(cells_simulated, 
                                                parameter_set['var_x'], 
                                                parameter_set['var_g'], 
                                                dt_measument,
                                                atol=dt/10)
        temp_dataset['lane'] = i+1
        # add the last "lane" to the data frame
        dataset = pd.concat([dataset, temp_dataset], ignore_index=True)
    
    print("SAVE", seed)
    # ----------- SAVE ----------- #
    directory, filename = get_next_file_name(out_dir)
    write_param_file(os.path.join(directory, "parameters_simulation.txt"), parameter_set) 

    mean_lambda_init = np.mean(np.concatenate([cell.naive_rate for cell in naive_rate_for_all_cells(df2cells(dataset), 1, naive_growth_rate_skip_div)]))
    
    mean_q_init = np.mean(np.concatenate([cell.naive_rate for cell in naive_rate_for_all_cells(df2cells(dataset), 1, naive_production_rate_skip_div,beta=parameter_set['beta'])]))
    
    param_file_temp = "parameters_noise{:d}.txt".format(0)
    param_file = write_param_file_var_noise(os.path.join(directory, param_file_temp), 
                                            parameter_set, 
                                            0.1, 
                                            mean_lambda_init,
                                            mean_q_init)


    out_file = param_file[:-4] + "_out"

    csv_config_file = write_csv_config(os.path.join(directory, "csv_config.txt"), lane="lane")
    dataset.to_csv(filename)
    
    print("RUN", seed)

    # ----------- RUN ----------- #
    run_cmd = "./../../RealTrace/bin/RealTrace" + \
                " -c " + csv_config_file + \
                " -b " + param_file + \
                " -t " + str(1e-40) + \
                " -i " + filename + \
                " -space log -m -p -div binomial -noise scaled" + \
                " -o " + out_file
    s = subprocess.run(run_cmd, capture_output=True, shell=True) 

        
        
        
def main():
    # main_sim(1)
    n_cores = os.environ['SLURM_JOB_CPUS_PER_NODE']
    print("slurm: number of cpus {}".format(n_cores))
    print(n_cores)
    with get_context("spawn").Pool(int(n_cores)) as p:
        print("Start multiprocessing")
        p.map(main_sim, np.arange(100))
        p.close()
        p.join()
        print("Join")
        
if __name__ == "__main__":
    main()
   