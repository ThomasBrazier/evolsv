rule download:
    """
    Download the long read fastq from the SRA archive
    and the genome from the NCBI genome assembly database
    """
    output:
        "{wdir}/{sra}.fastq.gz",
        "{wdir}/{sra}.fna",
        temporary("{wdir}/{sra}.zip"),
        "{wdir}/{sra}_config.yaml",
        "{wdir}/{sra}_assembly_data_report.jsonl",
        "{wdir}/{sra}_sequence_report.jsonl"
    conda:
        "../envs/download.yaml"
    log:
        "{wdir}/logs/{sra}_download.log"
    shell:
        """
        bash workflow/scripts/download.sh {sra} {genome} {wdir}
        cp config/config.yaml {wdir}/{sra}_config.yaml
        """
