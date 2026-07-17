#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu 2G
#SBATCH --time=12:00:00
#SBATCH --array=1-671%30
#SBATCH --output=/projects/psg/people/qvw641/ResQ_wolf/Heterozygosity_noDup/out/Het_%A_%a.log
#SBATCH --job-name Het
module load vcflib

# Script to estimate relative heterozygosity from imputed  VCF file 

scripts="/home/qvw641/ResQ_wolf/"
dir="/projects/psg/people/qvw641/ResQ_wolf/"

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${scripts}European_671Samples | awk '{print $1}')

vcf_dir="/projects/psg/people/qvw641/ResQ_wolf/VCF/combined_allchr_targetpanel.vcf.gz"
sites="29544701"

het_dir="${dir}/Heterozygosity_noDup/"
mkdir -p $het_dir
vcftools --bcf $vcf_dir  --indv $sample --recode --recode-INFO-all --stdout | vcfhetcount | tail -n 1 |  awk -v s=$sample -v g=$sites '{print s"\t"$1"\t"$1/g}'  > ${het_dir}/${sample}.het