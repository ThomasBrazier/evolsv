#!/bin/bash

# Variables
sra=$1
min_sv_size=$2
min_coverage=$3

# SV calling
mkdir ${sra}/cutesv
gunzip -f --keep ${sra}/ref_${sra}.fna.gz
rm -rf ${sra}/cutesv/*
cuteSV --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --max_size -1 --min_support $min_coverage --min_size $min_sv_size $sra/${sra}_sorted.bam $sra/ref_$sra.fna $sra/${sra}_cutesv.vcf ${sra}/cutesv/

