#!/usr/bin/env bash
# Generate the placeholder fixtures used by the CI dry-run tests (.test/test_dryrun.py).
#
# These are NOT test data. `snakemake -n` never reads the content of an input file,
# and the start_from_bam validation in workflow/rules/common.smk:17-56 only checks
# that each declared BAM/BAI/FASTQ exists and is non-empty. A few bytes per file is
# therefore enough to build the BAM-mode DAG. Running the pipeline for real on these
# files would fail immediately, by design.
#
# Fixtures are generated rather than committed because .gitignore excludes *.tsv,
# *.fa and *.txt repo-wide, and .test/fixtures/ outright.
#
# The files under fixtures/bad/ are deliberately malformed: each one trips exactly one
# of the WorkflowError branches in workflow/rules/common.smk so the negative tests can
# assert that the pipeline still refuses bad input at DAG-construction time.
set -euo pipefail

fixtures="$(dirname "$0")/fixtures"
bad="$fixtures/bad"
mkdir -p "$fixtures" "$bad"

# Genome accession must match config/samples.tsv so both dry-run modes agree.
genome="GCA_947247005.1"
sample="SAMEA8724893"

# The first of the two SRA runs listed in config/samples.tsv. samples_single_run.tsv
# keeps only this one, so the tests can prove the sheet drives the download_sra count.
sra="ERR10287556"

printf '>chr1\nACGTACGTACGT\n' > "$fixtures/ref.fna"

for aligner in minimap2 ngmlr; do
    printf 'placeholder\n' > "$fixtures/${sample}_${aligner}.bam"
    printf 'placeholder\n' > "$fixtures/${sample}_${aligner}.bam.bai"
done
printf 'placeholder\n' > "$fixtures/${sample}_reads.fastq.gz"

good_minimap2="$fixtures/${sample}_minimap2.bam"
good_ngmlr="$fixtures/${sample}_ngmlr.bam"
good_fastq="$fixtures/${sample}_reads.fastq.gz"

# ---------------------------------------------------------------------------
# Sample sheets for the FASTQ entry point
# ---------------------------------------------------------------------------

# A single SRA run. config/samples.tsv has two, so `download_sra` must drop to 1 here.
{
    printf 'sample_name\tsra\tgenome\n'
    printf '%s\t%s\t%s\n' "$sample" "$sra" "$genome"
} > "$fixtures/samples_single_run.tsv"

# ---------------------------------------------------------------------------
# Sample sheets for the pre-aligned BAM entry point
# ---------------------------------------------------------------------------

# Writes a start_from_bam sample sheet. Any argument may be the empty string, which is
# how the "column left blank" failure modes are expressed.
# usage: write_bam_sheet <path> <bam_minimap2> <bam_ngmlr> <fastq>
write_bam_sheet() {
    local path="$1" bam_minimap2="$2" bam_ngmlr="$3" fastq="$4"
    {
        printf 'sample_name\tsra\tgenome\tbam_minimap2\tbam_ngmlr\tfastq\n'
        # The sra column is intentionally empty: no SRA download happens in start_from_bam mode.
        printf '%s\t\t%s\t%s\t%s\t%s\n' \
            "$sample" "$genome" "$bam_minimap2" "$bam_ngmlr" "$fastq"
    } > "$path"
}

write_bam_sheet "$fixtures/samples_bam.tsv" "$good_minimap2" "$good_ngmlr" "$good_fastq"

# --- failure modes -----------------------------------------------------------

# check_readable_file: path does not exist.
write_bam_sheet "$bad/samples_missing_bam.tsv" \
    "$bad/does_not_exist_minimap2.bam" "$good_ngmlr" "$good_fastq"

# check_readable_file: path exists but is zero bytes.
: > "$bad/empty_minimap2.bam"
printf 'placeholder\n' > "$bad/empty_minimap2.bam.bai"
write_bam_sheet "$bad/samples_empty_bam.tsv" \
    "$bad/empty_minimap2.bam" "$good_ngmlr" "$good_fastq"

# resolve_bam_index: a readable BAM with neither <bam>.bai nor <stem>.bai next to it.
printf 'placeholder\n' > "$bad/noindex_minimap2.bam"
rm -f "$bad/noindex_minimap2.bam.bai" "$bad/noindex_minimap2.bai"
write_bam_sheet "$bad/samples_no_index.tsv" \
    "$bad/noindex_minimap2.bam" "$good_ngmlr" "$good_fastq"

# The fastq column is blank: SVJedi-graph has no reads to genotype with.
write_bam_sheet "$bad/samples_no_fastq.tsv" "$good_minimap2" "$good_ngmlr" ""

# The bam_ngmlr column is blank: the ensemble needs one BAM per aligner.
write_bam_sheet "$bad/samples_no_ngmlr.tsv" "$good_minimap2" "" "$good_fastq"

# ---------------------------------------------------------------------------
# Several individuals in one sample sheet
# ---------------------------------------------------------------------------

# A second individual with one SRA run. The accession is a placeholder: a dry run never
# downloads it. The first individual keeps its two runs from config/samples.tsv, so the
# tests can check that runs are merged per individual and never across individuals.
sample2="SAMEA0000002"
sra2="ERR00000001"

{
    printf 'sample_name\tsra\tgenome\n'
    printf '%s\t%s\t%s\n' "$sample" "ERR10287556" "$genome"
    printf '%s\t%s\t%s\n' "$sample" "ERR10287555" "$genome"
    printf '%s\t%s\t%s\n' "$sample2" "$sra2" "$genome"
} > "$fixtures/samples_multi.tsv"

for aligner in minimap2 ngmlr; do
    printf 'placeholder\n' > "$fixtures/${sample2}_${aligner}.bam"
    printf 'placeholder\n' > "$fixtures/${sample2}_${aligner}.bam.bai"
done
printf 'placeholder\n' > "$fixtures/${sample2}_reads.fastq.gz"

{
    printf 'sample_name\tsra\tgenome\tbam_minimap2\tbam_ngmlr\tfastq\n'
    printf '%s\t\t%s\t%s\t%s\t%s\n' "$sample" "$genome" "$good_minimap2" "$good_ngmlr" "$good_fastq"
    printf '%s\t\t%s\t%s\t%s\t%s\n' "$sample2" "$genome" \
        "$fixtures/${sample2}_minimap2.bam" "$fixtures/${sample2}_ngmlr.bam" \
        "$fixtures/${sample2}_reads.fastq.gz"
} > "$fixtures/samples_bam_multi.tsv"

# --- failure modes -----------------------------------------------------------

# Two individuals on two different reference genomes.
{
    printf 'sample_name\tsra\tgenome\n'
    printf '%s\t%s\t%s\n' "$sample" "$sra" "$genome"
    printf '%s\t%s\t%s\n' "$sample2" "$sra2" "GCA_000000001.1"
} > "$bad/samples_two_genomes.tsv"

# One individual declares two bam_minimap2 files on two rows.
{
    printf 'sample_name\tsra\tgenome\tbam_minimap2\tbam_ngmlr\tfastq\n'
    printf '%s\t\t%s\t%s\t%s\t%s\n' "$sample" "$genome" "$good_minimap2" "$good_ngmlr" "$good_fastq"
    printf '%s\t\t%s\t%s\t\t\n' "$sample" "$genome" "$fixtures/${sample2}_minimap2.bam"
} > "$bad/samples_duplicate_bam.tsv"

# A sample_name that cannot be a directory name.
{
    printf 'sample_name\tsra\tgenome\n'
    printf '%s\t%s\t%s\n' "bad/name" "$sra" "$genome"
} > "$bad/samples_bad_name.tsv"

echo "Wrote fixtures to $fixtures"
