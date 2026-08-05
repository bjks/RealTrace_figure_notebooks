#!/bin/bash

#SBATCH --job-name=runtime
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=1G

#SBATCH --time=23:50:00
#SBATCH --qos=1day

#SBATCH --output=std.out
#SBATCH --mail-type=END,FAIL,TIME_LIMIT

module load matplotlib
python3 simulation_test_runtime.py 