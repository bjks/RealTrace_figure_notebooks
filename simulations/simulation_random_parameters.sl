#!/bin/bash

#SBATCH --job-name=ggp
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

#SBATCH --time=6-23:00:00
#SBATCH --qos=1week

#SBATCH --output=std.out
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=bjoern.kscheschinski@unibas.sh

module load matplotlib
python3 simulation_random_parameters.py

