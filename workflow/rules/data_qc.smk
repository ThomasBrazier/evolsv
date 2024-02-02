rule download:
    """
    Download the long read fastq from the SRA archive
    and the genome from the NCBI genome assembly database
    """
    input:
        "{wdir}/{sra}.fastq.gz",
        "{wdir}/{sra}.fna.gz"
    output:
        "{wdir}/fastqc/{sra}"
    conda:
        "workflow/envs/fastqc.yaml"
    shell:
        """
        mkdir -p {wdir}/fastqc
        fastqc --noextract -outdir {wdir}/fastqc/ *.fastq.gz
        """
        