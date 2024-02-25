rule mapping:
    """
    Map reads to the reference genome with Minimap2
    """
    input:
        fastq = expand("{wdir}/fastq/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"]),
        fasta = "{wdir}/{genome}.fna",
        html = expand("{wdir}/fastqc/{sample}_fastqc.html", wdir=wdir, sample=samples["sra"]),
        qczip = expand("{wdir}/fastqc/{sample}_fastqc.zip", wdir=wdir, sample=samples["sra"]),
        nanoplot = expand("{wdir}/nanoplot/{sample}_NanoStats.txt", wdir=wdir, sample=samples["sra"])
    output:
        "{wdir}/{genome}.sam"
    threads: workflow.cores
    conda:
        "../envs/minimap2.yaml"
    log:
        "{wdir}/logs/{genome}_mapping.log"
    shell:
        """
        minimap2 -ax {config[minimap_ax]} --MD -2 --seed {config[minimap_seed]} --eqx -t {threads} --sam-hit-only {input.fasta} {input.fastq} > {output}
        """


rule samtools_view:
    """
    Transform the sam file to a bam file
    """
    input:
        "{wdir}/{genome}.sam"
    output:
        "{wdir}/{genome}.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -S -b {input} > {output}
        """


rule samtools_sort:
    """
    Sort the bam file
    """
    input:
        "{wdir}/{genome}.bam"
    output:
        "{wdir}/{genome}_sorted.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort {input} -o {output}
        """


rule samtools_index:
    """
    Create an index related file of the sorted bam file
    """
    input:
        "{wdir}/{genome}_sorted.bam"
    output:
        "{wdir}/{genome}_sorted.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input}
        """


rule samtools_stats:
    """
    Mapping QC
    """
    input:
        bam = "{wdir}/{genome}_sorted.bam",
        bai = "{wdir}/{genome}_sorted.bam.bai"
    output:
        stats = "{wdir}/mapping/{genome}_mapping.stats",
        stattsv = "{wdir}/mapping/{genome}_mapping.stats.tsv",
        plot = "{wdir}/mapping/{genome}_mapping_plot.html"
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/logs/{genome}_samtoolsstats.log"
    shell:
        """
        mkdir -p {wdir}/mapping
        samtools stats {input.bam} > {output.stats}
        cat {output.stats} | grep ^SN | cut -f 2- > {output.stattsv}
        # QC visualization
        plot-bamstats -p {wdir}/mapping/{genome}_mapping_plot {output.stats}
        """
