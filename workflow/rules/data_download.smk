rule download:
    """
    Download the long read fastq from the SRA archive
    and the genome from the NCBI genome assembly database
    """
    output:
        "{wdir}/{sra}.fastq.gz",
        "{wdir}/{sra}.fna",
        temporary("{wdir}/{sra}/"),
        temporary("{wdir}/{sra}.zip"),
        "{wdir}/{sra}_config.yaml"
    conda:
        "workflow/envs/download.yaml"
    shell:
        """
        scripts/download.sh {sra} {genome} {wdir}
        cp config/config.yaml {wdir}/{sra}_config.yaml
        """