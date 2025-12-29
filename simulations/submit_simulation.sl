#!/bin/bash

#SBATCH --job-name=auto_corr
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=1G

#SBATCH --time=23:50:00
#SBATCH --qos=1day

#SBATCH --output=std.out
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=bjoern.kscheschinski@unibas.ch

module load matplotlib
python3 cells_simulation.py --o /scicore/home/nimwegen/ksches0000/simulation_data/inference_output_1/ --n 20 --l 50 --d 10 --i /scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/constitexpr/20$
