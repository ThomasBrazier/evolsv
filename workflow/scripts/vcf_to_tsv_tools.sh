#! /bin/bash
# Convert the full vcf to tabular data frame for R scripts
# Remove REF/ALT sequences for lighter storage

input_vcf=$1
output_tsv=$2

cat $input_vcf | grep '#CHROM' | sed 's/#//' | awk -v OFS='\t' '{ print $1, $2, $3, $6, $7, $8, $9, "genotype" }' > $output_tsv
cat $input_vcf | grep -v '#' | awk -v OFS='\t' '{ print $1, $2, $3, $6, $7, $8, $9, $10 }' >> $output_tsv