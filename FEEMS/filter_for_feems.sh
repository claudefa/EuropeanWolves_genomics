conda activate feems_e
module load plink/2.0.0

# Convert to Plink
# removing one sample Novembre_Chinese_Wolf that has no coordinates in the metadata

plink2 --vcf  /projects/psg/people/qvw641/ResQ_wolf/VCF/phased.european.info08.name.norelated.maf01.wolves.vcf.gz --keep sampleslist --make-pgen vzs --out  modernWolf --double-id --dog --threads 8

# unique IDs
plink2 --pfile modernWolf vzs  --set-all-var-ids @:#:\$r:\$a  --make-pgen vzs --out modernWolf_uniq --double-id --dog

# LD-prune
plink2 --pfile modernWolf_uniq vzs --indep-pairwise 500 50 0.1 --out modernWolf --double-id --dog

plink2 --pfile modernWolf_uniq vzs --extract modernWolf.prune.in --make-bed --out modernWolf_pruned --threads 8  --double-id --dog


# Run feems
conda activate feems_e
sbatch feems_bash.sh
