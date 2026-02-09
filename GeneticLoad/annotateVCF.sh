conda activate j21

# Annotate VCFs
vcf="/projects/psg/people/pkb156/MW/popgen/phased.all.info08.name.norelated.vcf.gz"
outvcf="/projects/psg/people/qvw641/ResQ_wolf/GenLoad/SNPeff/Annotation/phased.all.info08.name.norelated.ann.vcf.gz"
java -Xmx8g -jar /home/qvw641/bin/snpEff/snpEff.jar CanFam3.1.99 $vcf | bgzip > $outvcf
tabix -p vcf  $outvcf

# Extract deleterious variants, from each category
zgrep -E "#|HIGH" phased.all.info08.name.norelated.ann.vcf.gz | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT]\n' > high.txt
zgrep -E "#|MODERATE" phased.all.info08.name.norelated.ann.vcf.gz | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT]\n' > moderate.txt
zgrep -E "#|LOW" phased.all.info08.name.norelated.ann.vcf.gz | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT]\n' > low.txt
