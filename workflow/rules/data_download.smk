rule download_sra:
    """
    Download the long read fastq from the SRA archive
    and the genome from the NCBI genome assembly database
    """
    output:
        "{wdir}/fastq/{sample}_sra.fastq.gz"
    conda:
        "../envs/download.yaml"
    shell:
        """
        mkdir --parents {wdir}/fastq
        fastq-dump -v --gzip --outdir {wdir}/fastq/ {wildcards.sample}
        mv "{wdir}/fastq/{wildcards.sample}.fastq.gz" "{wdir}/fastq/{wildcards.sample}_sra.fastq.gz"
        cp config/samples.tsv {wdir}/{genome}_samples.tsv
        """


rule download_genome:
    """
    Download the genome from the NCBI genome assembly database
    """
    input:
        expand("{wdir}/fastq/{sample}_sra.fastq.gz", wdir=wdir, sample=samples["sra"])
    output:
        "{wdir}/genome/{genome}.fna",
        temp("{wdir}/genome/{genome}.zip"),
        "{wdir}/genome/{genome}_config.yaml",
        "{wdir}/genome/{genome}_assembly_data_report.jsonl",
        "{wdir}/genome/{genome}_sequence_report.jsonl"
    conda:
        "../envs/download.yaml"
    shell:
        """
        SSL_CERT_FILE="ssl/cert.pem"
        echo $SSL_CERT_FILE

        datasets download genome accession {genome} --filename {wdir}/genome/{genome}.zip --include genome,gff3,seq-report
	    unzip -o {wdir}/genome/{genome}.zip -d {wdir}/genome/
	    cp {wdir}/genome/ncbi_dataset/data/{genome}/*_genomic.fna {wdir}/genome/{genome}.fna
        if test -f {wdir}/genome/ncbi_dataset/data/{genome}/genomic.gff
        then
        echo "GFF annotation exists."
        cp {wdir}/genome/ncbi_dataset/data/{genome}/genomic.gff {wdir}/genome/{genome}.gff
        fi
	    cp {wdir}/genome/ncbi_dataset/data/assembly_data_report.jsonl {wdir}/genome/{genome}_assembly_data_report.jsonl
	    cp {wdir}/genome/ncbi_dataset/data/{genome}/sequence_report.jsonl {wdir}/genome/{genome}_sequence_report.jsonl
        cp config/config.yaml {wdir}/genome/{genome}_config.yaml
        """


rule sample_ids:
    """
    Create a file with sample ids
    """
    input:
        expand("{wdir}/fastq/{sample}_sra.fastq.gz", wdir=wdir, sample=samples["sra"])
    output:
        sampleids = "{wdir}/{genome}.samples"
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/logs/{genome}_sample_ids.log"
    shell:
        """
        echo {sample_id} > {output.sampleids}
        """


rule merge_fastq:
    """
    Merge fastq files for mapping
    """
    input:
        fastq = expand("{wdir}/fastq/{sample}_sra.fastq.gz", wdir=wdir, sample=samples["sra"])
    output:
        merged_fastq = expand("{wdir}/fastq/{genome}.fastq.gz", wdir=wdir, genome=genome)
    conda:
        "../envs/samtools.yaml"
    params:
        fastqlist = " ".join(samplelist)
    shell:
        """
        cat {params.fastqlist} > {output.merged_fastq}
        """
