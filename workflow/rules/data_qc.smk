rule fastqc:
    """
    Report data quality for long reads, one SRA run at a time
    """
    input:
        "{wdir}/{sample}/fastq/{run}_sra.fastq.gz",
    output:
        html="{wdir}/{sample}/fastqc/{run}_sra_fastqc.html",
        qczip="{wdir}/{sample}/fastqc/{run}_sra_fastqc.zip",
    threads: workflow.cores
    conda:
        "../envs/fastqc.yaml"
    log:
        "{wdir}/{sample}/logs/{run}.fastqc.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.html),
    shell:
        """
        mkdir -p {params.outdir}
        fastqc --threads {resources.cpus_per_task} --outdir {params.outdir}/ {input} &> {log}
        """


rule hifiadapterfilt:
    """
    Remove PacBio HiFi reads that contain adapter sequences, one SRA run at a time.
    Only used when sequencing_technology is hifi (see rule merge_fastq).

    HiFiAdapterFilt: Sim et al. (2022) BMC Genomics 23:157,
    https://doi.org/10.1186/s12864-022-08375-1
    Reads with a BLAST match to the PacBio adapter or primer sequences are removed
    as whole reads; no bases are trimmed from the remaining reads.
    -l      Minimum adapter match length to remove a read [default: 44]
    -m      Minimum adapter match percentage to remove a read [default: 97]
    -t      Number of threads for blastn

    hifiadapterfilt.sh looks for its input as ${{prefix}}*.f*q* in the current
    directory. It runs in a private work directory, so that no other file can match.
    """
    input:
        "{wdir}/{sample}/fastq/{run}_sra.fastq.gz",
    output:
        filt_fastq=temp("{wdir}/{sample}/fastq/{run}_sra.filt.fastq.gz"),
        stats="{wdir}/{sample}/hifiadapterfilt/{run}.stats",
        blocklist="{wdir}/{sample}/hifiadapterfilt/{run}.blocklist",
        blastout="{wdir}/{sample}/hifiadapterfilt/{run}.contaminant.blastout",
    threads: workflow.cores
    conda:
        "../envs/hifiadapterfilt.yaml"
    log:
        "{wdir}/{sample}/logs/{run}.hifiadapterfilt.log",
    params:
        workdir=lambda wildcards, output: os.path.join(
            os.path.dirname(output.stats), f"{wildcards.run}_work"
        ),
    shell:
        """
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        ln -s "$(realpath {input})" {params.workdir}/{wildcards.run}.fastq.gz
        log="$(realpath {log})"

        ( cd {params.workdir} && hifiadapterfilt.sh \
        -p {wildcards.run} \
        -l {config[hifiadapterfilt_min_length]} \
        -m {config[hifiadapterfilt_min_match]} \
        -t {resources.cpus_per_task} \
        -o . ) &> "$log"

        mv {params.workdir}/{wildcards.run}.filt.fastq.gz {output.filt_fastq}
        mv {params.workdir}/{wildcards.run}.stats {output.stats}
        mv {params.workdir}/{wildcards.run}.blocklist {output.blocklist}
        mv {params.workdir}/{wildcards.run}.contaminant.blastout {output.blastout}
        rm -rf {params.workdir}
        """


rule nanoplot:
    """
    Quality control of raw data
    """
    input:
        fastq="{wdir}/{sample}/fastq/{run}_sra.fastq.gz",
        html="{wdir}/{sample}/fastqc/{run}_sra_fastqc.html",
        qczip="{wdir}/{sample}/fastqc/{run}_sra_fastqc.zip",
    output:
        "{wdir}/{sample}/nanoplot/{run}_NanoStats.txt",
        # "{wdir}/{sample}/nanoplot/{run}_LengthvsQualityScatterPlot_dot.html",
        # "{wdir}/{sample}/nanoplot/{run}_LengthvsQualityScatterPlot_dot.png",
        # "{wdir}/{sample}/nanoplot/{run}_LengthvsQualityScatterPlot_kde.html",
        # "{wdir}/{sample}/nanoplot/{run}_LengthvsQualityScatterPlot_kde.png",
        "{wdir}/{sample}/nanoplot/{run}_NanoPlot-report.html",
        "{wdir}/{sample}/nanoplot/{run}_Non_weightedHistogramReadlength.html",
        # "{wdir}/{sample}/nanoplot/{run}_Non_weightedHistogramReadlength.png",
        "{wdir}/{sample}/nanoplot/{run}_Non_weightedLogTransformed_HistogramReadlength.html",
        # "{wdir}/{sample}/nanoplot/{run}_Non_weightedLogTransformed_HistogramReadlength.png",
        "{wdir}/{sample}/nanoplot/{run}_WeightedHistogramReadlength.html",
        # "{wdir}/{sample}/nanoplot/{run}_WeightedHistogramReadlength.png",
        "{wdir}/{sample}/nanoplot/{run}_WeightedLogTransformed_HistogramReadlength.html",
        # "{wdir}/{sample}/nanoplot/{run}_WeightedLogTransformed_HistogramReadlength.png",
        "{wdir}/{sample}/nanoplot/{run}_Yield_By_Length.html",
        # "{wdir}/{sample}/nanoplot/{run}_Yield_By_Length.png"
    threads: workflow.cores
    conda:
        "../envs/nanoplot.yaml"
    log:
        "{wdir}/{sample}/logs/{run}_nanoplot.log",
    shell:
        """
        NanoPlot --fastq {input.fastq} -t {resources.cpus_per_task} --tsv_stats --outdir {wildcards.wdir}/{wildcards.sample}/nanoplot/ --prefix '{wildcards.run}_' --N50 --no_static --verbose --title {wildcards.run}
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
        reads="{wdir}/{sample}/fastq/{genome}.fastq.gz",
    output:
        filtered_reads=temp("{wdir}/{sample}/fastq/{genome}_filtered.fastq.gz"),
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
        fastq="{wdir}/{sample}/fastq/{genome}_filtered.fastq.gz",
    output:
        "{wdir}/{sample}/nanoplot_filtered/{genome}_NanoStats.txt",
        "{wdir}/{sample}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_dot.html",
        # "{wdir}/{sample}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_dot.png",
        "{wdir}/{sample}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_kde.html",
        # "{wdir}/{sample}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_kde.png",
        "{wdir}/{sample}/nanoplot_filtered/{genome}_NanoPlot-report.html",
        "{wdir}/{sample}/nanoplot_filtered/{genome}_Non_weightedHistogramReadlength.html",
        # "{wdir}/{sample}/nanoplot_filtered/{genome}_Non_weightedHistogramReadlength.png",
        "{wdir}/{sample}/nanoplot_filtered/{genome}_Non_weightedLogTransformed_HistogramReadlength.html",
        # "{wdir}/{sample}/nanoplot_filtered/{genome}_Non_weightedLogTransformed_HistogramReadlength.png",
        "{wdir}/{sample}/nanoplot_filtered/{genome}_WeightedHistogramReadlength.html",
        # "{wdir}/{sample}/nanoplot_filtered/{genome}_WeightedHistogramReadlength.png",
        "{wdir}/{sample}/nanoplot_filtered/{genome}_WeightedLogTransformed_HistogramReadlength.html",
        # "{wdir}/{sample}/nanoplot_filtered/{genome}_WeightedLogTransformed_HistogramReadlength.png",
        "{wdir}/{sample}/nanoplot_filtered/{genome}_Yield_By_Length.html",
        # "{wdir}/{sample}/nanoplot_filtered/{genome}_Yield_By_Length.png"
    conda:
        "../envs/nanoplot.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_nanoplot_filtered.log",
    shell:
        """
        NanoPlot --fastq {input.fastq} -t {resources.cpus_per_task} --tsv_stats --outdir {wildcards.wdir}/{wildcards.sample}/nanoplot_filtered/ --prefix '{wildcards.genome}_' --N50 --no_static --verbose --title {wildcards.sample}
        """
