#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu 2G
#SBATCH --time=12:00:00
##SBATCH --array=1-15 # for local populations
#SBATCH --array=1-2 # for metapopulations
#SBATCH --output=/projects/psg/people/qvw641/ResQ_wolf/CurrentNe/out/Plink_ped_%A_%a.log
#SBATCH --job-name Ped

module unload plink
module load plink/1.9.0

vcf="/projects/psg/people/qvw641/ResQ_wolf/VCF/combined.europeanwolves.norelated.maf01.vcf.gz"
dir="/projects/psg/people/qvw641/ResQ_wolf/CurrentNe/"

# local Ne estimates
#group=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/ResQ_wolf/Groups/Groupings | awk '{print $1}')

# metapopulation for broader Northern and Central Europe
group=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/ResQ_wolf/Groups/Groupings_metapop | awk '{print $1}')

outdir=${dir}${group}
mkdir -p ${outdir}

echo $group

vcftools --gzvcf ${vcf} --keep /home/qvw641/ResQ_wolf/Groups/${group} --plink --out ${outdir}/$group

for i in {1..10}
do
	awk '{print $2}' ${outdir}/${group}.map | shuf -n 200000 > ${outdir}/random_snps_${i}.txt
	plink --file ${outdir}/${group}  --cm-map /projects/psg/people/qvw641/ResQ_wolf/HapNe/VCFs/@.average_canFam3.1.txt  --extract ${
outdir}/random_snps_${i}.txt --recode --out ${outdir}/${group}_subset_${i} --dog
done
