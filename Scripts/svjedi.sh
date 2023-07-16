#!/bin/bash

# Genotyping
svjedi -v Vcf/merged.vcf -r Ref/ref_v_cardui_ERR6608653.fna -i Fastq/v_cardui_ERR6608653.fastq.gz -o Vcf/genotypes.vcf
