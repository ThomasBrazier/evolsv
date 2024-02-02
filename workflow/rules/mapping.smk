rule mapping:
    """
    Map reads to the reference genome
    """
    input:
        "{wdir}/{sra}.fastq.gz",
        "{wdir}/ref_{sra}.fna.gz"
    output:
        "{wdir}/{sra}.sam"
    conda:
        "Envs/minimap2.yaml"
    shell:
        "Scripts/minimap2.sh {sra}"

rule samtools_view:
    """
    Transform the sam file to a bam file
    """
    input:
        "{wdir}/{sra}.sam"
    output:
        "{wdir}/{sra}.bam"
    conda:
        "Envs/samtools.yaml"
    shell:
        "Scripts/samtools_view.sh {sra}"

rule samtools_sort:
    """
    Sort the bam file
    """
    input:
        "{wdir}/{sra}.bam"
    output:
        "{wdir}/{sra}_sorted.bam"
    conda:
        "Envs/samtools.yaml"
    shell:
        "Scripts/samtools_sort.sh {sra}"

rule samtools_index:
    """
    Create an index related file of the sorted bam file
    """
    input:
        "{wdir}/{sra}_sorted.bam"
    output:
        "{wdir}/{sra}_sorted.bam.bai"
    conda:
        "Envs/samtools.yaml"
    shell:
        "Scripts/samtools_index.sh {sra}"