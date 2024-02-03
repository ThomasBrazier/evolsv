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
    conda:
        "workflow/envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax {config[minimap_ax]} --MD --eqX -t {workflow.threads} --sam-hit-only {input.fasta} {input.fastq} > {output}
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
        "workflow/envs/samtools.yaml"
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
        "workflow/envs/samtools.yaml"
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
        "workflow/envs/samtools.yaml"
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
        plot = "{wdir}/mapping/{sra}_mapping_plot.html"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        """
        mkdir -p {wdir}/mapping
        samtools stats {input.bam} > {output.stats}
        # QC visualization
        plot-bamstats -p {output.plot} {output.stats}
        """