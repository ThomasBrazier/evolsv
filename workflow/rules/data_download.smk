rule download_sra:
    """
    Download the long read fastq from the SRA archive
    and the genome from the NCBI genome assembly database
    """
    output:
        expand("{wdir}/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"])
    conda:
        "../envs/download.yaml"
    log:
        expand("{wdir}/logs/{sample}_downloadsra.log", wdir=wdir, sample=samples["sra"])
    shell:
        """
        bash workflow/scripts/download_sra.sh {output} {wdir}
        cp config/samples.tsv {wdir}/{genome}_samples.tsv
        """


rule download_genome:
    """
    Download the genome from the NCBI genome assembly database
    """
    output:
        "{wdir}/{genome}.fna",
        temporary("{wdir}/{genome}.zip"),
        "{wdir}/{genome}_config.yaml",
        "{wdir}/{genome}_assembly_data_report.jsonl",
        "{wdir}/{genome}_sequence_report.jsonl"
    conda:
        "../envs/download.yaml"
    log:
        "{wdir}/logs/{genome}_downloadgenome.log"
    shell:
        """
        bash workflow/scripts/download_genome.sh {genome} {wdir}
        cp config/config.yaml {wdir}/{genome}_config.yaml
        """
