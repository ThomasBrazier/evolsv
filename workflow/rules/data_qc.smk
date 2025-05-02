rule fastqc:
    """
    Report data quality for long reads
    """
    input:
        expand("{wdir}/fastq/{sample}_sra.fastq.gz", wdir=wdir, sample=samples["sra"])
    output:
        expand("{wdir}/fastqc/{sample}_sra_fastqc.html", wdir=wdir, sample=samples["sra"]),
        expand("{wdir}/fastqc/{sample}_sra_fastqc.zip", wdir=wdir, sample=samples["sra"])
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
        fastq = "{wdir}/fastq/{sample}_sra.fastq.gz",
        html = "{wdir}/fastqc/{sample}_sra_fastqc.html",
        qczip = "{wdir}/fastqc/{sample}_sra_fastqc.zip",
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


rule filter_reads_chopper:
    """
    Filter long reads with chopper
    --headcrop      Trim N nucleotides from the start of a read
    --maxlength     Sets a maximum read length
    -l, --minlength     Sets a minimum read length
    -q, --quality       Sets a minimum Phred average quality score
    --tailcrop      Trim N nucleotides from the end of a read
    """
    input:
        reads = "{wdir}/fastq/{genome}.fastq.gz"
    output:
        filtered_reads = "{wdir}/fastq/{genome}_filtered.fastq.gz"
    conda:
        "../envs/chopper.yaml"
    shell:
        """
        chopper -q {config[chopper_quality]} \
        -l {config[chopper_minlength]} \
        --maxlength {config[chopper_maxlength]} \
        --headcrop {config[chopper_headcrop]} \
        --tailcrop {config[chopper_tailcrop]} \
        --threads {resources.cpus_per_task} \
        -i {input.reads} | gzip > {output.filtered_reads}
        """


rule nanoplot_after_filtering:
    """
    Quality control after filtering long reads
    """
    input:
        fastq = "{wdir}/fastq/{genome}_filtered.fastq.gz"
    output:
        "{wdir}/nanoplot_filtered/{genome}_NanoStats.txt",
        "{wdir}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_dot.html",
        "{wdir}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_dot.png",
        "{wdir}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_kde.html",
        "{wdir}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_kde.png",
        "{wdir}/nanoplot_filtered/{genome}_NanoPlot-report.html",
        "{wdir}/nanoplot_filtered/{genome}_Non_weightedHistogramReadlength.html",
        "{wdir}/nanoplot_filtered/{genome}_Non_weightedHistogramReadlength.png",
        "{wdir}/nanoplot_filtered/{genome}_Non_weightedLogTransformed_HistogramReadlength.html",
        "{wdir}/nanoplot_filtered/{genome}_Non_weightedLogTransformed_HistogramReadlength.png",
        "{wdir}/nanoplot_filtered/{genome}_WeightedHistogramReadlength.html",
        "{wdir}/nanoplot_filtered/{genome}_WeightedHistogramReadlength.png",
        "{wdir}/nanoplot_filtered/{genome}_WeightedLogTransformed_HistogramReadlength.html",
        "{wdir}/nanoplot_filtered/{genome}_WeightedLogTransformed_HistogramReadlength.png",
        "{wdir}/nanoplot_filtered/{genome}_Yield_By_Length.html",
        "{wdir}/nanoplot_filtered/{genome}_Yield_By_Length.png"
    threads: workflow.cores
    conda:
        "../envs/nanoplot.yaml"
    log:
        "{wdir}/logs/{genome}_nanoplot_filtered.log"
    shell:
        """
        NanoPlot --fastq {input.fastq} -t {threads} --tsv_stats --outdir {wdir}/nanoplot_filtered/ --prefix '{wildcards.genome}_' --N50 --verbose --title {wildcards.genome}
        """
