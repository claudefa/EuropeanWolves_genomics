#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu 100G
#SBATCH --time=200:00:00
#SBATCH --output=/projects/psg/people/qvw641/ResQ_wolf/FEEMS/out/Feems_dense_ASIA_%A_%a.log
#SBATCH --job-name FEEMS

conda activate feems_e

python run_feems_claudia.py


