rule download_sra:
    """
    Download the long read fastq from the SRA archive
    and the genome from the NCBI genome assembly database
    """
    output:
        "{wdir}/fastq/{sample}.fastq.gz"
    conda:
        "../envs/download.yaml"
    log:
        "{wdir}/logs/{sample}_download_sra.log"
    shell:
        """
        mkdir --parents {wdir}/fastq
        fastq-dump -v --gzip --outdir {wdir}/fastq/ {wildcards.sample}
        cp config/samples.tsv {wdir}/{genome}_samples.tsv
        """


rule download_genome:
    """
    Download the genome from the NCBI genome assembly database
    """
    input:
        expand("{wdir}/fastq/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"])
    output:
        "{wdir}/{genome}.fna",
        temporary("{wdir}/{genome}.zip"),
        # "{wdir}/{genome}.gff",
        "{wdir}/{genome}_config.yaml",
        "{wdir}/{genome}_assembly_data_report.jsonl",
        "{wdir}/{genome}_sequence_report.jsonl"
    conda:
        "../envs/download.yaml"
    log:
        "{wdir}/logs/{genome}_downloadgenome.log"
    shell:
        """
        datasets download genome accession {genome} --filename {wdir}/{genome}.zip --include genome,gff3,seq-report
	    unzip {wdir}/{genome}.zip -d {wdir}/
	    cp {wdir}/ncbi_dataset/data/{genome}/*_genomic.fna {wdir}/{genome}.fna
        if test -f {wdir}/ncbi_dataset/data/{genome}/genomic.gff
        then
        echo "GFF annotation exists."
        cp {wdir}/ncbi_dataset/data/{genome}/genomic.gff {wdir}/{genome}.gff
        fi
	    cp {wdir}/ncbi_dataset/data/assembly_data_report.jsonl {wdir}/{genome}_assembly_data_report.jsonl
	    cp {wdir}/ncbi_dataset/data/{genome}/sequence_report.jsonl {wdir}/{genome}_sequence_report.jsonl
        cp config/config.yaml {wdir}/{genome}_config.yaml
        """
