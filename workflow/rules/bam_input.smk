"""
Entry point for datasets that are already aligned (config `start_from_bam: true`).

This file replaces rules/data_qc.smk and rules/mapping.smk: the SRA download, the
read-level QC (FastQC, NanoPlot), the chopper filtering and both alignments are
skipped. Instead, the BAM files declared in the sample sheet are exposed under the
canonical names the rest of the workflow expects, so no downstream rule changes.

One BAM per aligner is required. Reusing a single alignment under both aligner names
would make rule jasmine count the same evidence twice (it writes both BAM paths into
its bam_list) and would inflate the ensemble consensus.

Alignment-level QC (rules/mapping_qc.smk, rules/mappability.smk) runs as usual.
"""


rule stage_bam:
    """
    Validate a user-supplied pre-aligned BAM and expose it at the canonical path.

    The BAM is symlinked rather than copied: these files routinely reach tens of GB
    and a copy would double the storage cost of the run for no benefit.
    Validation runs first, so a rejected BAM never leaves a canonical path behind.
    """
    input:
        bam = lambda wildcards: input_bams[wildcards.aligner],
        bai = lambda wildcards: input_bais[wildcards.aligner],
        fai = "{wdir}/genome/{genome}.fna.fai"
    output:
        bam = "{wdir}/bam/{genome}_{aligner}_sorted.bam",
        bai = "{wdir}/bam/{genome}_{aligner}_sorted.bam.bai",
        report = "{wdir}/bam/{genome}_{aligner}_bam_check.txt"
    conda:
        "../envs/bamcheck.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_stage_bam.log"
    shell:
        """
        mkdir --parents {wdir}/bam

        python workflow/scripts/check_bam_reference.py \
        --bam {input.bam} --fai {input.fai} \
        --sample-id {sample_id} --aligner {wildcards.aligner} \
        --report {output.report} 2> {log}

        ln -sf $(realpath {input.bam}) {output.bam}
        ln -sf $(realpath {input.bai}) {output.bai}
        """


rule stage_fastq:
    """
    Expose the reads declared in the sample sheet where the genotyping rules expect
    them.

    SVJedi-graph maps reads onto a variation graph, so the five genotyping rules need
    the reads themselves and cannot work from the BAM. The reads are used as supplied:
    unlike the FASTQ entry point, no chopper filtering is applied (see the caveat in
    each bam_check.txt report). The `_filtered` name is kept so that the genotyping
    rules are identical in both input modes.

    A single file is symlinked; several are concatenated, which is valid for gzip and
    matches what rule merge_fastq does in FASTQ mode.
    """
    input:
        fastq = input_fastqs
    output:
        merged_fastq = "{wdir}/fastq/{genome}_filtered.fastq.gz"
    log:
        "{wdir}/logs/{genome}_stage_fastq.log"
    shell:
        """
        mkdir --parents {wdir}/fastq

        nfastq=$(echo {input.fastq} | wc -w)
        if [ "$nfastq" -eq 1 ]
        then
        ln -sf $(realpath {input.fastq}) {output.merged_fastq} 2> {log}
        else
        cat {input.fastq} > {output.merged_fastq} 2> {log}
        fi
        """
