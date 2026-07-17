#!/bin/bash
#SBATCH -c 10
#SBATCH --mem-per-cpu 20G
#SBATCH --time=48:00:00
##SBATCH --array=1-1117%50
#SBATCH --output=/projects/psg/people/qvw641/ResQ_wolf/Heterozygosity_ANGSD_notrans/out/Heterozygosity_angsd_notrans_%A_%a.log
#SBATCH --job-name Het_angsd

module load angsd/0.940
module load winsfs

dir="/projects/psg/people/qvw641/ResQ_wolf/Heterozygosity_ANGSD_notrans/"
ref="/projects/psg/data/genomes/dog2/canFam31/canFam31.fasta"

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/ResQ_wolf/list_angsd | awk '{print $1}')
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/ResQ_wolf/list_angsd | awk '{print $2}')
cov=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/ResQ_wolf/list_angsd | awk '{print $3}')


#min coverage set at 3 but max coverage set at dp10 for low coverage (below 3x), and 3 times average coverage.

#angsd  -i ${sample} -ref $ref -anc $ref -out ${dir}${name} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1  -C 50 -baq 1 -minMapQ 30 -minQ 20  -setMinDepth 3 -setMaxDepth $cov -doCounts 1 -nThreads 10 -GL 2 -doSaf 1 -checkBamHeaders 0
#winsfs ${dir}${name}.saf.idx > ${dir}${name}.ml

# transversion
angsd  -i $sample -ref $ref -anc $ref -out ${dir}${name}_rmtrans -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -rmTrans 1 -C 50 -baq 1 -minMapQ 30 -minQ 20 -setMinDepth 3 -setMaxDepth $cov -doCounts 1 -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -nThreads 10 -checkBamHeaders 0


zcat ${dir}${name}_rmtrans.mafs.gz  | awk '{print $1 , $2}' | sed 1d > ${dir}${name}_rmtrans.list
angsd sites index ${dir}${name}_rmtrans.list
echo "index done"

angsd  -i ${sample} -ref $ref -sites ${dir}/${name}_rmtrans.list -anc $ref -out ${dir}${name}_filtList_notrans -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1  -noTrans 1 -C 50 -baq 1 -minMapQ 30 -minQ 20  -setMinDepth 3 -setMaxDepth $cov  -doCounts 1 -nThreads 10 -GL 2 -doSaf 1 -checkBamHeaders 0

winsfs ${dir}${name}_filtList_notrans.saf.idx > ${dir}${name}_notrans_list.ml

