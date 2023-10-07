#!/bin/bash
sra=$1


# Genotyping
svjedi -v $sra/${sra}_merged.vcf -r $sra/ref_${sra}.fna -i $sra/${sra}.fastq.gz -o $sra/${sra}_genotypes.vcf
