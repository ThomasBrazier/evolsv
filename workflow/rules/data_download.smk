rule download:
    """
    Download the long read fastq from the SRA archive
    and the genome from the NCBI genome assembly database
    """
    output:
        "{wdir}/{sra}.fastq.gz",
        "{wdir}/{sra}.fna.gz",
        temporary("{wdir}/{sra}/"),
        temporary("{wdir}/{sra}.zip")
    conda:
        "workflow/envs/download.yaml"
    shell:
        "worklow/scripts/download.sh {sra} {genome} {wdir}"