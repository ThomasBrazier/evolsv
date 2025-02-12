rule download_sra:
    """
    Download the long read fastq from the SRA archive
    and the genome from the NCBI genome assembly database
    """
    output:
        expand("{wdir}/fastq/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"])
    conda:
        "../envs/download.yaml"
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
        temp("{wdir}/{genome}.zip"),
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
	    unzip -o {wdir}/{genome}.zip -d {wdir}/
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


rule sample_ids:
    """
    Create a file with sample ids
    """
    input:
        expand("{wdir}/fastq/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"])
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
        fastq = expand("{wdir}/fastq/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"])
    output:
        merged_fastq = temp(expand("{wdir}/fastq/{genome}.fastq.gz", wdir=wdir, genome=genome))
    conda:
        "../envs/samtools.yaml"
    params:
        fastqlist = " ".join(samplelist)
    shell:
        """
        cat {params.fastqlist} > {output.merged_fastq}
        """
