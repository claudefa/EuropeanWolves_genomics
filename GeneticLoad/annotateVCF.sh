conda activate j21

# Annotate VCFs
vcf="/projects/psg/people/pkb156/MW_imputation/impute/combine/combined_allchr_targetpanel.vcf.gz"
outvcf="/projects/psg/people/qvw641/ResQ_wolf/GenLoad/SNPeff/Annotation/combined_allchr_targetpanel.ann.vcf.gz"

java -Xmx8g -jar /home/qvw641/bin/snpEff/snpEff.jar CanFam3.1.99 $vcf | bgzip > $outvcf
tabix -p vcf  $outvcf

# Extract different profiles
zgrep -E "#|HIGH" $outvcf | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT]\n' > high.txt
zgrep -E "#|MODERATE" $outvcf | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT]\n' > moderate.txt
zgrep -E "#|LOW" $outvcf | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT]\n' > low.txt
