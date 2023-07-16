#!/bin/bash

# Variables
sra=$1
min_sv_size=$2
min_coverage=$3

# SV calling
sniffles --input $sra/${sra}_sorted.bam --vcf $sra/${sra}_sniffles.vcf --reference $sra/ref_$sra.fna.gz --threads 4 --allow-overwrite --minsvlen $min_sv_size --qc-coverage $min_coverage

