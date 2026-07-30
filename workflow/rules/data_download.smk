if not bam_mode:

    rule download_sra:
        """
        Download the long read fastq from the SRA archive
        and the genome from the NCBI genome assembly database
        """
        output:
            temp("{wdir}/fastq/{sample}_sra.fastq.gz")
        conda:
            "../envs/download.yaml"
        shell:
            """
            mkdir --parents {wdir}/fastq
            fastq-dump -v --gzip --outdir {wdir}/fastq/ {wildcards.sample}
            mv "{wdir}/fastq/{wildcards.sample}.fastq.gz" "{wdir}/fastq/{wildcards.sample}_sra.fastq.gz"
            """


rule download_genome:
    """
    Download the reference genome and its metadata from the NCBI genome assembly
    database.

    The assembly metadata (assembly_data_report.jsonl, sequence_report.jsonl) is
    always downloaded, because rule autosomes_sexchromosomes and the final report
    read it. Only the FASTA itself can be replaced by a local copy, via the
    `reference_fasta` config key: a pre-aligned BAM must be checked against the exact
    FASTA it was aligned to, which is not necessarily the NCBI GenBank copy.
    """
    output:
        "{wdir}/genome/{genome}.fna",
        "{wdir}/genome/{genome}.fna.fai",
        temp("{wdir}/genome/{genome}.zip"),
        "{wdir}/genome/{genome}_config.yaml",
        "{wdir}/genome/{genome}_assembly_data_report.jsonl",
        "{wdir}/genome/{genome}_sequence_report.jsonl"
    conda:
        "../envs/download.yaml"
    params:
        local_fasta = config.get("reference_fasta", "")
    shell:
        """
        datasets download genome accession {genome} --filename {wdir}/genome/{genome}.zip --include genome,gff3,seq-report
	    unzip -o {wdir}/genome/{genome}.zip -d {wdir}/genome/

        if [ -n "{params.local_fasta}" ]
        then
        echo "Using the local reference FASTA {params.local_fasta} instead of the NCBI copy."
        ln -sf $(realpath {params.local_fasta}) {wdir}/genome/{genome}.fna
        else
	    cp {wdir}/genome/ncbi_dataset/data/{genome}/*_genomic.fna {wdir}/genome/{genome}.fna
        fi
        samtools faidx {wdir}/genome/{genome}.fna

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

    Also keeps a copy of the sample sheet next to the results for provenance. That
    copy used to be made by rule download_sra, which does not run when starting from
    pre-aligned BAM files.
    """
    output:
        sampleids = "{wdir}/{genome}.samples"
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/logs/{genome}_sample_ids.log"
    shell:
        """
        mkdir --parents {wdir}
        echo {sample_id} > {output.sampleids}
        cp {config[samples]} {wdir}/{genome}_samples.tsv
        """


if not bam_mode:

    rule merge_fastq:
        """
        Merge fastq files for mapping
        """
        input:
            fastq = expand("{wdir}/fastq/{sample}_sra.fastq.gz", wdir=wdir, sample=samples["sra"])
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
