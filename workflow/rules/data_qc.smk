rule fastqc:
    """
    Report data quality for long reads
    """
    input:
        "{wdir}/{sra}.fastq.gz"
    output:
        "{wdir}/fastqc/{sra}"
    conda:
        "workflow/envs/fastqc.yaml"
    shell:
        """
        mkdir -p {wdir}/fastqc
        fastqc --noextract -outdir {wdir}/fastqc/ *.fastq.gz
        """
        