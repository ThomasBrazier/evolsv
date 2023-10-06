#!/bin/bash

# Variables
sra=$1
min_sv_size=$2
min_coverage=$3
svim_quality=$4

# SV calling
svim alignment $sra $sra/${sra}_sorted.bam $sra/ref_$sra.fna.gz --insertion_sequences --read_names --min_sv_size $min_sv_size --minimum_depth $min_coverage
# SVIM does not filter SV itself and outputs all variants
bcftools view -i "QUAL >= $svim_quality" ${sra}/variants.vcf > $sra/${sra}_svim.vcf
cat ${sra}/${sra}_svim.vcf | grep -v svim.BND > $sra/${sra}_svim_noBND.vcf

#mv ${sra}/variants.vcf $sra/${sra}_svim.vcf


