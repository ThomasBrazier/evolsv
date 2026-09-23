#! /bin/bash
# Convert the full vcf to tabular data frame for R scripts
# Remove REF/ALT sequences for lighter storage
#
# The VCF holds one sample column per aligner + SV caller callset, then one last sample:
# the SVjedi-graph genotypes merged in by rule all_samples_vcf. The callsets are named on
# the command line rather than hardcoded here, because how many there are depends on the
# `aligners` config key, and their ORDER is the order of Jasmine's SUPP_VEC bits (see
# `callsets` in workflow/Snakefile). workflow/scripts/merging_qc.R and finalQC.Rmd read
# both the shape of the run and that bit order out of the header written here.
#
# usage: vcf_to_tsv.sh <input.vcf> <output.tsv> <callset1,callset2,...>

input_vcf=$1
output_tsv=$2
callsets=$3

if [ -z "$callsets" ]; then
    echo "vcf_to_tsv.sh: no callset names given (third argument)" >&2
    exit 1
fi

# Number of callset columns, i.e. of VCF sample columns before the genotyped one.
n=$(awk -F, '{ print NF }' <<< "$callsets")

# Header: the 7 fixed fields, one column per callset, then the genotyped sample, whose
# name is taken from the VCF itself.
grep '#CHROM' "$input_vcf" | sed 's/#//' |
    awk -v OFS='\t' -v callsets="$callsets" -v n="$n" '{
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $6, $7, $8, $9
        split(callsets, name, ",")
        for (i = 1; i <= n; i++) printf "\t%s", name[i]
        printf "\t%s\n", $(10 + n)
    }' > "$output_tsv"

# Data: the same columns, dropping REF ($4) and ALT ($5).
grep -v '#' "$input_vcf" |
    awk -v OFS='\t' -v n="$n" '{
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $6, $7, $8, $9
        for (i = 10; i <= 9 + n; i++) printf "\t%s", $i
        printf "\t%s\n", $(10 + n)
    }' >> "$output_tsv"
