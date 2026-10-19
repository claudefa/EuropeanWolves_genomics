#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 3G
#SBATCH --time=24:00:00
#SBATCH --array=1-15 # local populations
#SBATCH --array=1-2 # metapopulations
#SBATCH --output=/projects/psg/people/qvw641/ResQ_wolf/CurrentNe/out/currentNe2_%A_%a.log
#SBATCH --job-name currentNe2

dir="/projects/psg/people/qvw641/ResQ_wolf/CurrentNe/"

#group=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/ResQ_wolf/Groups/Groupings | awk '{print $1}')
group=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/ResQ_wolf/Groups/Groupings_metapop | awk '{print $1}')

peddir=${dir}${group}

echo $group

for i in {1..10}
do
#	/home/qvw641/bin/currentNe2/currentne2 -k 0 ${peddir}/${group}_subset_${i}.ped 38 -o ${dir}/CurrNe_${group}_rep${i} # only for local Ne
	/home/qvw641/bin/currentNe2/currentne2 -x ${peddir}/${group}_subset_${i}.ped 38 -o ${dir}/CurrNe2_${group}_rep${i}
done
