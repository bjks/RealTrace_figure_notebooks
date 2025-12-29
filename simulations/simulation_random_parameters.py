import sys

import sys
from cells_simulation import *

import sys
sys.path.insert(0, "../../")

from RealTrace_paper_analysis.plotting_src.read_ggp_run import * # you might need to change these two line, such that you can use read ggp_run

import numpy as np
import subprocess
from multiprocessing import Pool, get_context

np.random.seed(1)

# ========== Simulation parameters ========== #
save_dataset = True    # saves data set as csv
run_ggp = True         # runs the ggp code (located in "../bin/" relative to this nb)

dt = 1e-2
dt_measument = 1 # in minutes
n_cells = 100 
discard_cells = 5
lanes = 5

div_mode = "adder"

division_log_length = 1+np.log(2)   # for sizer: division, when log_length hits division_log_length
division_time = 3000 - 1e-10          # for timer: division, when cell cycle time hits division_time
division_addition = np.log(2)       # for adder: divsion, when division_addition in log_length was added in cell cycle

n=10000
parameter_ranges = {"mean_lambda":  np.geomspace(0.001, 0.01, n), 
                    "gamma_lambda": np.geomspace(0.01,  0.1, n), 
                    # "var_lambda":   np.geomspace(1e-8,  1e-4, n), 
                    "mean_q":       np.geomspace(10,     100, n), 
                    "gamma_q":      np.geomspace(0.01,  0.1, n), 
                    # "var_q":        np.geomspace(0.01,  1, n),
                    "beta":         np.geomspace(1e-4,  1e-3, n),
                    "var_x":        np.geomspace(1e-4,  1e-3, n),
                    "var_g":        np.geomspace(0.3,   3, n),
                    "var_dx":       np.geomspace(1e-4,  1e-2, n),
                    "var_dg":       np.geomspace(1,  100, n)}

cv_range = np.linspace(0.2, 0.8, n)

def main_sim():
    out_dir = mk_missing_dir('../../data_sim_manuscript/random_parameters_control_cv_0.3_var_dg'+ '_n{:d}_{:s}'.format(n_cells*lanes, div_mode))
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
    
    # ----------- SAVE ----------- #
    directory, filename = get_next_file_name(out_dir)
    write_param_file(os.path.join(directory, "parameters_simulation.txt"), parameter_set) 

    for j in range(1):
        param_file_temp = "parameters_noise{:d}.txt".format(j)
        param_file = write_param_file_noise(os.path.join(directory, param_file_temp), 
                                                parameter_set, 
                                                rel_step=0.1, 
                                                log_noise=0.3,
                                                fixed = "beta")


        out_file = param_file[:-4] + "_out"

        csv_config_file = write_csv_config(os.path.join(directory, "csv_config.txt"), lane="lane")
        dataset.to_csv(filename)
        
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
    n_cores = os.environ['SLURM_JOB_CPUS_PER_NODE']
    print("slurm: number of cpus {}".format(n_cores))
    print(n_cores)
    with get_context("spawn").Pool(int(n_cores)) as p:
        print("Start multiprocessing")
        p.starmap(main_sim, np.arange(100))
        p.close()
        p.join()
        print("Join")
        
if __name__ == "__main__":
    main()
   