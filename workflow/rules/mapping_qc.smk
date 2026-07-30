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
        bam = "{wdir}/bam/{genome}_{aligner}_sorted.bam",
        bai = "{wdir}/bam/{genome}_{aligner}_sorted.bam.bai"
    output:
        stats = "{wdir}/mapping_QC/{genome}_{aligner}_mapping.stats",
        stattsv = "{wdir}/mapping_QC/{genome}_{aligner}_mapping.stats.tsv",
        plot = "{wdir}/mapping_QC/{genome}_{aligner}_mapping_plot.html"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p {wdir}/mapping_QC
        samtools stats {input.bam} > {output.stats}
        cat {output.stats} | grep ^SN | cut -f 2- > {output.stattsv}
        # QC visualization
        plot-bamstats -p {wdir}/mapping_QC/{genome}_{wildcards.aligner}_mapping_plot {output.stats}
        """


rule samtools_coverage:
    """
    Mapping coverage along the genome
    """
    input:
        bam = "{wdir}/bam/{genome}_{aligner}_sorted.bam",
        bai = "{wdir}/bam/{genome}_{aligner}_sorted.bam.bai"
    output:
        coverage = "{wdir}/mapping_QC/{genome}_{aligner}_coverage.tsv",
        coverage_hist = "{wdir}/mapping_QC/{genome}_{aligner}_coverage_hist.txt",
        coverage_depth = "{wdir}/mapping_QC/{genome}_{aligner}_coverage_depthplot.txt"
    conda:
        "../envs/samtools_coverage.yaml"
    shell:
        """
        # store coverage in a tsv/tab-separated file
        samtools coverage {input.bam} > {output.coverage}

        # get coverage as a histogram
        samtools coverage {input.bam} --histogram > {output.coverage_hist}

        # get coverage as depth plot
        samtools coverage {input.bam} --plot-depth > {output.coverage_depth}
        """
