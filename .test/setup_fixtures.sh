#!/usr/bin/env bash
# Generate the placeholder fixtures used by the CI dry-run jobs.
#
# These are NOT test data. `snakemake -n` never reads the content of an input file,
# and the start_from_bam validation in workflow/rules/common.smk:17-56 only checks
# that each declared BAM/BAI/FASTQ exists and is non-empty. A few bytes per file is
# therefore enough to build the BAM-mode DAG. Running the pipeline for real on these
# files would fail immediately, by design.
#
# Fixtures are generated rather than committed because .gitignore excludes *.tsv,
# *.fa and *.txt repo-wide.
set -euo pipefail

fixtures="$(dirname "$0")/fixtures"
mkdir -p "$fixtures"

# Genome accession must match config/samples.tsv so both dry-run modes agree.
genome="GCA_947247005.1"
sample="SAMEA8724893"

printf '>chr1\nACGTACGTACGT\n' > "$fixtures/ref.fna"

for aligner in minimap2 ngmlr; do
    printf 'placeholder\n' > "$fixtures/${sample}_${aligner}.bam"
    printf 'placeholder\n' > "$fixtures/${sample}_${aligner}.bam.bai"
done
printf 'placeholder\n' > "$fixtures/${sample}_reads.fastq.gz"

# The sra column is intentionally empty: no SRA download happens in start_from_bam mode.
{
    printf 'sample_name\tsra\tgenome\tbam_minimap2\tbam_ngmlr\tfastq\n'
    printf '%s\t\t%s\t%s\t%s\t%s\n' \
        "$sample" \
        "$genome" \
        "$fixtures/${sample}_minimap2.bam" \
        "$fixtures/${sample}_ngmlr.bam" \
        "$fixtures/${sample}_reads.fastq.gz"
} > "$fixtures/samples_bam.tsv"

echo "Wrote fixtures to $fixtures"
