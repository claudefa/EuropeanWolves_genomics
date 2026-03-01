#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 10G
#SBATCH --time=100:00:00
#SBATCH --array=1-15 
#SBATCH --output=/projects/psg/people/qvw641/ResQ_wolf/HapNe/out/HapNe_LD_all_%A_%a.log
#SBATCH --job-name HapNe

conda activate hapne-snakemake

dir="/projects/psg/people/qvw641/ResQ_wolf/HapNe/"
group=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/ResQ_wolf/Groups/Groupings | awk '{print $1}')

for i in {1..10};
do
        echo "Try$i"
	mkdir -p ${dir}/results_LD_${group}
	mkdir -p ${dir}/results_LD_${group}/steps/
	mkdir -p ${dir}/results_LD_${group}/steps/logs/

	snakemake -c8 --configfile ${dir}/ConfigFiles/config_LD_${group}.yaml
        mv ${dir}/results_LD_${group} ${dir}/${i}Try/
        echo "moved directory"
done
