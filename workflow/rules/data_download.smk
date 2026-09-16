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


if reference_source == "local":

    rule stage_genome:
        """
        Stage a local reference genome and its local sequence report. Nothing is
        downloaded from NCBI (config keys reference_fasta and sequence_report, see
        reference_source in rules/common.smk).

        The FASTA is symlinked (it can be several GB). The sequence report is copied,
        so the run keeps the exact version it used. Rule check_reference_names then
        checks that the report and the FASTA name the same contigs.
        """
        input:
            fasta=reference_fasta,
            sequence_report=sequence_report,
        output:
            fasta="{wdir}/genome/{genome}.fna",
            fai="{wdir}/genome/{genome}.fna.fai",
            config="{wdir}/genome/{genome}_config.yaml",
            sequence_report="{wdir}/genome/{genome}_sequence_report.jsonl",
        conda:
            "../envs/samtools.yaml"
        shell:
            """
            echo "Using the local reference FASTA {input.fasta} instead of the NCBI copy."
            ln -sf $(realpath {input.fasta}) {output.fasta}
            samtools faidx {output.fasta}
            cp {input.sequence_report} {output.sequence_report}
            cp config/config.yaml {output.config}
            """

    if assembly_data_report:

        rule stage_assembly_data_report:
            """
            Copy the local assembly_data_report.jsonl (config key assembly_data_report).
            Only the final report reads it.
            """
            input:
                assembly_data_report,
            output:
                "{wdir}/genome/{genome}_assembly_data_report.jsonl",
            shell:
                """
                cp {input} {output}
                """

else:

    rule download_genome:
        """
        Download the reference genome and its metadata from the NCBI genome assembly
        database.

        The FASTA can be replaced by a local copy, via the `reference_fasta` config key:
        a pre-aligned BAM must be checked against the exact FASTA it was aligned to,
        which is not necessarily the NCBI GenBank copy. The metadata is then still
        downloaded; to use local metadata too, set sequence_report (rule stage_genome).
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
            local_fasta=reference_fasta,
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

        The runs are merged after adapter removal: rule hifiadapterfilt for
        sequencing_technology hifi, rule porechop_abi for ont.
        """
        input:
            fastq=lambda wildcards: expand(
                {
                    "hifi": "{wdir}/{sample}/fastq/{run}_sra.filt.fastq.gz",
                    "ont": "{wdir}/{sample}/fastq/{run}_sra.porechop.fastq.gz",
                }[sequencing_technology],
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
