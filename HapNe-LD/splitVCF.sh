#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu 2G
#SBATCH --time=12:00:00
#SBATCH --array=1-15
#SBATCH --output=/projects/psg/people/qvw641/ResQ_wolf/HapNe/out/SplitVCF%A_%a.log
#SBATCH --job-name SplitVCF

module unload plink
module load plink/1.9.0

vcf="/projects/psg/people/qvw641/ResQ_wolf/VCF/combined.europeanwolves.norelated.maf01.vcf.gz"
dir="/projects/psg/people/qvw641/ResQ_wolf/HapNe/"

group=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/ResQ_wolf/Groups/Groupings | awk '{print $1}')
#group=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/ResQ_wolf/Groups_Pure80/Groupings | awk '{print $1}')

outdir=${dir}"/Groups/"${group}
mkdir -p ${outdir}

echo $group

for chr in {1..38}; do
  echo $chr
   bcftools view -S /home/qvw641/ResQ_wolf/Groups/${group} -r chr"$chr" -Oz -o "${outdir}/${chr}.vcf.gz" $vcf
   bcftools index "${outdir}/${chr}.vcf.gz"
   plink --vcf ${outdir}/${chr}.vcf.gz --cm-map ${dir}/VCFs/${chr}.average_canFam3.1.txt ${chr} --recode --out ${outdir}/${chr}.shapeit --double-id --dog
  sed 's/^/chr/' ${outdir}/${chr}.shapeit.map > ${outdir}/${chr}.chr.shapeit.map
  echo "done"
done
