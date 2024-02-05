rule fastqc:
    """
    Report data quality for long reads
    """
    input:
        "{wdir}/{sra}.fastq.gz"
    output:
        "{wdir}/fastqc/{sra}_fastqc.html",
        "{wdir}/fastqc/{sra}_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p {wdir}/fastqc
        fastqc --outdir {wdir}/fastqc/ {input}
        """


rule nanoplot:
    """
    Quality control of raw data
    """
    input:
        fastq = "{wdir}/{sra}.fastq.gz",
        html = "{wdir}/fastqc/{sra}_fastqc.html",
        zip = "{wdir}/fastqc/{sra}_fastqc.zip"
    output:
        "{wdir}/nanoplot/{sra}_NanoStats.txt",
        "{wdir}/nanoplot/{sra}_LengthvsQualityScatterPlot_dot.html",
        "{wdir}/nanoplot/{sra}_LengthvsQualityScatterPlot_dot.png",
        "{wdir}/nanoplot/{sra}_LengthvsQualityScatterPlot_kde.html",
        "{wdir}/nanoplot/{sra}_LengthvsQualityScatterPlot_kde.png",
        "{wdir}/nanoplot/{sra}_NanoPlot-report.html",
        "{wdir}/nanoplot/{sra}_Non_weightedHistogramReadlength.html",
        "{wdir}/nanoplot/{sra}_Non_weightedHistogramReadlength.png",
        "{wdir}/nanoplot/{sra}_Non_weightedLogTransformed_HistogramReadlength.html",
        "{wdir}/nanoplot/{sra}_Non_weightedLogTransformed_HistogramReadlength.png",
        "{wdir}/nanoplot/{sra}_WeightedHistogramReadlength.html",
        "{wdir}/nanoplot/{sra}_WeightedHistogramReadlength.png",
        "{wdir}/nanoplot/{sra}_WeightedLogTransformed_HistogramReadlength.html",
        "{wdir}/nanoplot/{sra}_WeightedLogTransformed_HistogramReadlength.png",
        "{wdir}/nanoplot/{sra}_Yield_By_Length.html",
        "{wdir}/nanoplot/{sra}_Yield_By_Length.png"
    threads: workflow.cores
    conda:
        "../envs/nanoplot.yaml"
    shell:
        """
        NanoPlot --fastq {input.fastq} -t {threads} --tsv_stats --outdir {wdir}/nanoplot/ --prefix '{sra}_' --N50 --verbose --title {sra}
        """

