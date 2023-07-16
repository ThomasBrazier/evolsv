#!/bin/bash

# Variables
sra=$1

# Merging
jasmine file_list=$sra/${sra}_vcf_list.txt out_file=$sra/${sra}_merged.vcf genome_file=$sra/ref_$sra.fna.gz out_dir=$sra/ bam_list=$sra/${sra}_bam_list.txt  --output_genotypes --run_iris

