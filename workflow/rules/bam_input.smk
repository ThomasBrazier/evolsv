"""
Entry point for datasets that are already aligned (config `start_from_bam: true`).

This file replaces rules/data_qc.smk and rules/mapping.smk: the SRA download, the
read-level QC (FastQC, NanoPlot), the chopper filtering and both alignments are
skipped. Instead, the BAM files declared in the sample sheet are exposed under the
canonical names the rest of the workflow expects, so no downstream rule changes.

One BAM per aligner is required. Reusing a single alignment under both aligner names
would make rule jasmine count the same evidence twice (it writes both BAM paths into
its bam_list) and would inflate the ensemble consensus.

The reads are also needed, for SVJedi-graph. Two rules can supply them, one per
individual: rule stage_fastq when the sample sheet declares a `fastq` file, and rule
bam_to_fastq when it does not. Both write the same canonical path, so their `sample`
wildcard domains are made disjoint (see fastq_individuals and bam_fastq_individuals in
rules/common.smk).

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
        bam=lambda wildcards: input_bams[(wildcards.sample, wildcards.aligner)],
        bai=lambda wildcards: input_bais[(wildcards.sample, wildcards.aligner)],
        fai="{wdir}/genome/{genome}.fna.fai",
    output:
        bam="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam.bai",
        report="{wdir}/{sample}/bam/{genome}_{aligner}_bam_check.txt",
    params:
        # Where the reads used for genotyping come from, recorded in the report so the
        # provenance of a run is readable from its output alone.
        reads_source=lambda wildcards: (
            "sample-sheet" if input_fastqs[wildcards.sample] else "bam-derived"
        ),
        reads_aligner=reads_aligner,
    conda:
        "../envs/bamcheck.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_{aligner}_stage_bam.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_{aligner}.stage_bam.tsv"
    shell:
        """
        mkdir --parents {wdir}/{wildcards.sample}/bam

        python workflow/scripts/check_bam_reference.py \
        --bam {input.bam} --fai {input.fai} \
        --sample-id {wildcards.sample} --aligner {wildcards.aligner} \
        --reads-source {params.reads_source} \
        --reads-aligner {params.reads_aligner} \
        --report {output.report} 2> {log}

        ln -sf $(realpath {input.bam}) {output.bam}
        ln -sf $(realpath {input.bai}) {output.bai}
        """


if fastq_individuals:

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
        wildcard_constraints:
            # Only the individuals that declare a `fastq` file: rule bam_to_fastq writes
            # the same path for the others.
            sample="|".join(re.escape(i) for i in fastq_individuals),
        input:
            fastq=lambda wildcards: input_fastqs[wildcards.sample],
        output:
            merged_fastq="{wdir}/{sample}/fastq/{genome}_filtered.fastq.gz",
        log:
            "{wdir}/{sample}/logs/{genome}_stage_fastq.log",
        benchmark:
            "{wdir}/{sample}/benchmarks/{genome}.stage_fastq.tsv"
        shell:
            """
            mkdir --parents {wdir}/{wildcards.sample}/fastq

            nfastq=$(echo {input.fastq} | wc -w)
            if [ "$nfastq" -eq 1 ]
            then
            ln -sf $(realpath {input.fastq}) {output.merged_fastq} 2> {log}
            else
            cat {input.fastq} > {output.merged_fastq} 2> {log}
            fi
            """


if bam_fastq_individuals:

    rule bam_to_fastq:
        """
        Recover the reads from the alignment, for the individuals whose sample sheet
        leaves the `fastq` column blank.

        SVJedi-graph cannot genotype from a linear BAM, but the reads themselves are
        still in it, so requiring the original FASTQ alongside a BAM is unnecessary.
        The first selected aligner's BAM is used (minimap2 whenever it is selected, see
        `aligners` in workflow/Snakefile); all the alignments of one individual are
        assumed to come from the same read set (see README), so the choice only matters
        through what each aligner kept.

        The extracted reads are NOT identical to the sequenced reads:

        * -F 0x900 drops secondary (0x100) and supplementary (0x800) records. Without it
          every clipped segment of a split long-read alignment would be emitted as its
          own read, inflating both the coverage the callers see and the read support
          SVJedi-graph counts.
        * unmapped records (0x4) are kept, but reads the aligner never wrote cannot be
          recovered. A BAM produced with, for instance, minimap2 --sam-hit-only (which
          rule minimap2 itself uses) or filtered to mapped reads yields fewer reads than
          the original FASTQ.
        * samtools fastq restores reverse-strand reads to their original orientation,
          but bases removed by *hard* clipping on a primary alignment are gone.
        * no chopper filtering is applied, exactly as for rule stage_fastq. The
          `_filtered` name is kept only so the downstream rules are identical in both
          input modes.

        The staged BAM is the input rather than the sample-sheet path, so extraction
        only happens once rule stage_bam has accepted the BAM.
        """
        wildcard_constraints:
            # Only the individuals that declare no `fastq` file; the others are handled
            # by rule stage_fastq, which writes the same path.
            sample="|".join(re.escape(i) for i in bam_fastq_individuals),
        input:
            bam="{{wdir}}/{{sample}}/bam/{{genome}}_{aligner}_sorted.bam".format(
                aligner=reads_aligner
            ),
        output:
            fastq=temp("{wdir}/{sample}/fastq/{genome}_filtered.fastq.gz"),
        conda:
            "../envs/samtools_fastq.yaml"
        log:
            "{wdir}/{sample}/logs/{genome}_bam_to_fastq.log",
        benchmark:
            "{wdir}/{sample}/benchmarks/{genome}.bam_to_fastq.tsv"
        shell:
            """
            mkdir --parents {wdir}/{wildcards.sample}/fastq

            samtools fastq -F 0x900 --threads {resources.cpus_per_task} {input.bam} \
            2> {log} | bgzip --threads {resources.cpus_per_task} > {output.fastq}
            """
