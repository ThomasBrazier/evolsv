#! /bin/bash
# Convert the full vcf to tabular data frame for R scripts
# Remove REF/ALT sequences for lighter storage

input_vcf=$1
output_tsv=$2

cat $input_vcf | grep '#CHROM' | sed 's/#//' | awk -v OFS='\t' '{ print $1, $2, $3, $6, $7, $8, $9, "minimap2_sniffles", "minimap2_svim", "minimap2_cutesv", "minimap2_debreak", "ngmlr_sniffles", "ngmlr_svim", "ngmlr_cutesv", "ngmlr_debreak", $18 }' > $output_tsv
cat $input_vcf | grep -v '#' | awk -v OFS='\t' '{ print $1, $2, $3, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18 }' >> $output_tsv