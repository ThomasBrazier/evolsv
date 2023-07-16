#!/bin/bash

# Variables
sra=$1
min_sv_size=$2
min_coverage=$3

# SV calling
svim alignment $sra $sra/${sra}_sorted.bam $sra/ref_$sra.fna.gz --insertion_sequences --min_sv_size $min_sv_size --minimum_depth $min_coverage
mv ${sra}/variants.vcf $sra/${sra}_svim.vcf
