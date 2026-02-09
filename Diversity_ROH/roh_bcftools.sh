# Estimate ROH using bcftools with dog recomb maps
bcftools roh -G 30  --estimate-AF European_650Samples  /projects/psg/people/qvw641/ResQ_wolf/VCF/phased.european.info08.name.norelated.maf01.wolves.vcf.gz -m genetic_map/{CHROM}_average_canFam3.1_forROH.txt | grep -E '^RG|^#' > modernWolf_AF_map.txt
