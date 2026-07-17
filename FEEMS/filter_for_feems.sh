conda activate feems_e
module load plink/2.0.0

# Convert to Plink, without captive
plink2 --vcf   /projects/psg/people/pkb156/MW_sub2/popgen/combined.europeanwolves.norelated.maf01.vcf.gz --keep sampleslist_nocaptive --make-pgen vzs --out modernWolf_nocapt --double-id --dog --threads 8
plink2 --pfile modernWolf_nocapt vzs  --set-all-var-ids @:#:\$r:\$a  --make-pgen vzs --out modernWolf_nocapt_uniq --double-id --dog
plink2 --pfile modernWolf_nocapt_uniq vzs --indep-pairwise 500 50 0.1 --out modernWolf_nocapt --double-id --dog
plink2 --pfile modernWolf_nocapt_uniq vzs --extract modernWolf_nocapt.prune.in --make-bed --out modernWolf_nocapt_pruned --threads 8  --double-id --dog

# Run feems
conda activate feems_e
sbatch feems_bash.sh
# Run feems
conda activate feems_e
sbatch feems_bash.sh
