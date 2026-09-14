if not bam_mode:

    rule download_sra:
        """
        Download the long read fastq from the SRA archive
        and the genome from the NCBI genome assembly database
        """
        output:
            fastq=temp("{wdir}/{sample}/fastq/{run}_sra.fastq.gz"),
        params:
            outdir=lambda wildcards, output: os.path.dirname(output.fastq),
        conda:
            "../envs/download.yaml"
        log:
            "{wdir}/{sample}/logs/download_sra/{run}.log",
        shell:
            """
            mkdir --parents {params.outdir}
            fastq-dump -v --gzip --outdir {params.outdir}/ {wildcards.run} &> {log}
            mv "{params.outdir}/{wildcards.run}.fastq.gz" "{output.fastq}" 2>> {log}
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
        "{wdir}/genome/{genome}_sequence_report.jsonl",
    conda:
        "../envs/download.yaml"
    params:
        local_fasta=config.get("reference_fasta", ""),
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
        sampleids="{wdir}/{sample}/{genome}.samples",
        sheet="{wdir}/{sample}/{genome}_samples.tsv",
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_sample_ids.log",
    shell:
        """
        echo {wildcards.sample} > {output.sampleids}
        cp {config[samples]} {output.sheet}
        """


if not bam_mode:

    rule merge_fastq:
        """
        Merge the fastq files of one individual for mapping

        For sequencing_technology hifi, the runs are merged after adapter filtering
        with rule hifiadapterfilt. For ont, the raw runs are merged.
        """
        input:
            fastq=lambda wildcards: expand(
                (
                    "{wdir}/{sample}/fastq/{run}_sra.filt.fastq.gz"
                    if sequencing_technology == "hifi"
                    else "{wdir}/{sample}/fastq/{run}_sra.fastq.gz"
                ),
                wdir=wildcards.wdir,
                sample=wildcards.sample,
                run=runs_of(wildcards.sample),
            ),
        output:
            merged_fastq=temp("{wdir}/{sample}/fastq/{genome}.fastq.gz"),
        conda:
            "../envs/samtools.yaml"
        shell:
            """
            cat {input.fastq} > {output.merged_fastq}
            """
