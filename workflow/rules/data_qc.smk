rule fastqc:
    """
    Report data quality for long reads
    """
    input:
        expand("{wdir}/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"])
    output:
        expand("{wdir}/fastqc/{sample}.fastqc.html", wdir=wdir, sample=samples["sra"]),
        expand("{wdir}/fastqc/{sample}.fastqc.zip", wdir=wdir, sample=samples["sra"])
    conda:
        "../envs/fastqc.yaml"
    log:
        expand("{wdir}/{sample}.fastqc.log", wdir=wdir, sample=samples["sra"])
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
        fastq = "{wdir}/{sample}.fastq.gz",
        html = "{wdir}/fastqc/{sample}.fastqc.html",
        qczip = "{wdir}/fastqc/{sample}.fastqc.zip",
    output:
        "{wdir}/nanoplot/{sample}_NanoStats.txt",
        "{wdir}/nanoplot/{sample}_LengthvsQualityScatterPlot_dot.html",
        "{wdir}/nanoplot/{sample}_LengthvsQualityScatterPlot_dot.png",
        "{wdir}/nanoplot/{sample}_LengthvsQualityScatterPlot_kde.html",
        "{wdir}/nanoplot/{sample}_LengthvsQualityScatterPlot_kde.png",
        "{wdir}/nanoplot/{sample}_NanoPlot-report.html",
        "{wdir}/nanoplot/{sample}_Non_weightedHistogramReadlength.html",
        "{wdir}/nanoplot/{sample}_Non_weightedHistogramReadlength.png",
        "{wdir}/nanoplot/{sample}_Non_weightedLogTransformed_HistogramReadlength.html",
        "{wdir}/nanoplot/{sample}_Non_weightedLogTransformed_HistogramReadlength.png",
        "{wdir}/nanoplot/{sample}_WeightedHistogramReadlength.html",
        "{wdir}/nanoplot/{sample}_WeightedHistogramReadlength.png",
        "{wdir}/nanoplot/{sample}_WeightedLogTransformed_HistogramReadlength.html",
        "{wdir}/nanoplot/{sample}_WeightedLogTransformed_HistogramReadlength.png",
        "{wdir}/nanoplot/{sample}_Yield_By_Length.html",
        "{wdir}/nanoplot/{sample}_Yield_By_Length.png"
    threads: workflow.cores
    conda:
        "../envs/nanoplot.yaml"
    log:
        "{wdir}/logs/{sample}_nanoplot.log"
    shell:
        """
        NanoPlot --fastq {input.fastq} -t {threads} --tsv_stats --outdir {wdir}/nanoplot/ --prefix '{wildcards.sample}_' --N50 --verbose --title {wildcards.sample}
        """

