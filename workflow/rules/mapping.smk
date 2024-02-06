rule mapping:
    """
    Map reads to the reference genome with Minimap2
    """
    input:
        fastq = "{wdir}/{sra}.fastq.gz",
        fasta = "{wdir}/{sra}.fna",
        html = "{wdir}/fastqc/{sra}_fastqc.html",
        zip = "{wdir}/fastqc/{sra}_fastqc.zip",
        nanoplot = "{wdir}/nanoplot/{sra}_NanoStats.txt"
    output:
        "{wdir}/{sra}.sam"
    threads: workflow.cores
    conda:
        "../envs/minimap2.yaml"
    log:
        "{wdir}/logs/{sra}_mapping.log"
    shell:
        """
        minimap2 -ax {config[minimap_ax]} --MD -2 --seed {config[minimap_seed]} --eqx -t {threads} --sam-hit-only {input.fasta} {input.fastq} > {output}
        """


rule samtools_view:
    """
    Transform the sam file to a bam file
    """
    input:
        "{wdir}/{sra}.sam"
    output:
        "{wdir}/{sra}.bam"
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
        "{wdir}/{sra}.bam"
    output:
        "{wdir}/{sra}_sorted.bam"
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
        "{wdir}/{sra}_sorted.bam"
    output:
        "{wdir}/{sra}_sorted.bam.bai"
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
        bam = "{wdir}/{sra}_sorted.bam",
        bai = "{wdir}/{sra}_sorted.bam.bai"
    output:
        stats = "{wdir}/mapping/{sra}_mapping.stats",
        stattsv = "{wdir}/mapping/{sra}_mapping.stats.tsv",
        plot = "{wdir}/mapping/{sra}_mapping_plot.html"
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/logs/{sra}_samtoolsstats.log"
    shell:
        """
        mkdir -p {wdir}/mapping
        samtools stats {input.bam} > {output.stats}
        cat {output.stats} | grep ^SN | cut -f 2- > {output.stattsv}
        # QC visualization
        plot-bamstats -p {wdir}/mapping/{sra}_mapping_plot {output.stats}
        """
