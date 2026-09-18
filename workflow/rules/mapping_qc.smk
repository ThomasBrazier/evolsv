"""
Alignment-level QC.

These rules only need a sorted BAM and its index, so they run in both input modes:
after mapping (rules/mapping.smk) or on user-supplied BAMs (rules/bam_input.smk).
Their outputs feed rule final_report, which reads
mapping_QC/{genome}_{aligner}_mapping.stats.tsv for both aligners.
"""


rule samtools_stats:
    """
    Mapping QC
    """
    input:
        bam="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam.bai",
    output:
        stats="{wdir}/{sample}/mapping_QC/{genome}_{aligner}_mapping.stats",
        stattsv="{wdir}/{sample}/mapping_QC/{genome}_{aligner}_mapping.stats.tsv",
        plot="{wdir}/{sample}/mapping_QC/{genome}_{aligner}_mapping_plot.html",
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_{aligner}.samtools_stats.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_{aligner}.samtools_stats.tsv"
    shell:
        """
        mkdir -p {wdir}/{wildcards.sample}/mapping_QC
        samtools stats {input.bam} > {output.stats}
        cat {output.stats} | grep ^SN | cut -f 2- > {output.stattsv}
        # QC visualization
        plot-bamstats -p {wdir}/{wildcards.sample}/mapping_QC/{genome}_{wildcards.aligner}_mapping_plot {output.stats}
        """


rule samtools_coverage:
    """
    Mapping coverage along the genome
    """
    input:
        bam="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam.bai",
    output:
        coverage="{wdir}/{sample}/mapping_QC/{genome}_{aligner}_coverage.tsv",
        coverage_hist="{wdir}/{sample}/mapping_QC/{genome}_{aligner}_coverage_hist.txt",
        coverage_depth="{wdir}/{sample}/mapping_QC/{genome}_{aligner}_coverage_depthplot.txt",
    conda:
        "../envs/samtools_coverage.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_{aligner}.samtools_coverage.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_{aligner}.samtools_coverage.tsv"
    shell:
        """
        # store coverage in a tsv/tab-separated file
        samtools coverage {input.bam} > {output.coverage}

        # get coverage as a histogram
        samtools coverage {input.bam} --histogram > {output.coverage_hist}

        # get coverage as depth plot
        samtools coverage {input.bam} --plot-depth > {output.coverage_depth}
        """
