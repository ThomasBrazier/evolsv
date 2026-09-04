rule minimap2:
    """
    Map reads to the reference genome with Minimap2
    """
    input:
        fastq = "{wdir}/fastq/{genome}_filtered.fastq.gz",
        fasta = "{wdir}/genome/{genome}.fna",
        html = expand("{wdir}/fastqc/{sample}_sra_fastqc.html", wdir=wdir, sample=samples["sra"]),
        qczip = expand("{wdir}/fastqc/{sample}_sra_fastqc.zip", wdir=wdir, sample=samples["sra"]),
        nanoplot = expand("{wdir}/nanoplot/{sample}_NanoStats.txt", wdir=wdir, sample=samples["sra"]),
        nanoplot_filtered = "{wdir}/nanoplot_filtered/{genome}_NanoStats.txt"
    output:
        sam = temp("{wdir}/bam/{genome}_minimap2.sam")
    conda:
        "../envs/minimap2.yaml"
    log:
        "{wdir}/logs/{genome}_minimap2.log"
    shell:
        """
        minimap2 -ax {config[minimap_ax]} --MD -2 \
        --seed {config[minimap_seed]} --eqx \
        -t {resources.cpus_per_task} \
        -R "@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:{config[read_group_platform]}" \
        --sam-hit-only {input.fasta} {input.fastq} > {output.sam}
        """



rule ngmlr:
    """
    Map reads to the reference genome with ngmlr
    --rg-id <string>
        Adds RG:Z:<string> to all alignments in SAM/BAM [none]
    --rg-sm <string>
        RG header: Sample [none]
    --rg-pl <string>
        RG header: Platform [none]
    -t <int>,  --threads <int>
        Number of threads [1]
    -x <pacbio, ont>,  --presets <pacbio, ont>
        Parameter presets for different sequencing technologies [pacbio]
        Set from config[ngmlr_preset], which the sequencing_technology preset in
        rules/common.smk fills in (hifi -> pacbio, ont -> ont).
    -i <0-1>,  --min-identity <0-1>
        Alignments with an identity lower than this threshold will be discarded [0.65]
    -R <int/float>,  --min-residues <int/float>
        Alignments containing less than <int> or (<float> * read length) residues will be discarded [0.25]
    --no-smallinv
        Don't detect small inversions [false]
    --no-lowqualitysplit
        Split alignments with poor quality [false]
    """
    input:
        fastq = "{wdir}/fastq/{genome}_filtered.fastq.gz",
        fasta = "{wdir}/genome/{genome}.fna",
        html = expand("{wdir}/fastqc/{sample}_sra_fastqc.html", wdir=wdir, sample=samples["sra"]),
        qczip = expand("{wdir}/fastqc/{sample}_sra_fastqc.zip", wdir=wdir, sample=samples["sra"]),
        nanoplot = expand("{wdir}/nanoplot/{sample}_NanoStats.txt", wdir=wdir, sample=samples["sra"]),
        nanoplot_filtered = "{wdir}/nanoplot_filtered/{genome}_NanoStats.txt"
    output:
        sam = temp("{wdir}/bam/{genome}_ngmlr.sam")
    conda:
        "../envs/ngmlr.yaml"
    log:
        "{wdir}/logs/{genome}_ngmlr.log"
    shell:
        """
        ngmlr -t {resources.cpus_per_task} \
        -r {input.fasta} -q {input.fastq} \
        --presets {config[ngmlr_preset]} \
        --min-identity {config[min-identity]} \
        --rg-id {sample_id} --rg-sm {sample_id} \
        --rg-pl {config[read_group_platform]} \
        -o {output.sam}
        """


rule samtools_view:
    """
    Transform the sam file to a bam file
    """
    input:
        sam_minimap2 = "{wdir}/bam/{genome}_minimap2.sam",
        sam_ngmlr = "{wdir}/bam/{genome}_ngmlr.sam"
    output:
        bam_minimap2 = temp("{wdir}/bam/{genome}_minimap2.bam"),
        bam_ngmlr = temp("{wdir}/bam/{genome}_ngmlr.bam")
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -S -b {input.sam_minimap2} > {output.bam_minimap2}
        samtools view -S -b {input.sam_ngmlr} > {output.bam_ngmlr}
        """


rule samtools_sort:
    """
    Sort the bam file
    """
    input:
        "{wdir}/bam/{genome}_{aligner}.bam"
    output:
        "{wdir}/bam/{genome}_{aligner}_sorted.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort {input} -o {output}
        """


rule samtools_index:
    """
    Create an index related file of the sorted bam file
    """
    input:
        "{wdir}/bam/{genome}_{aligner}_sorted.bam"
    output:
        "{wdir}/bam/{genome}_{aligner}_sorted.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input}
        """


#
# Alignment-level QC (samtools_stats, samtools_coverage) lives in
# rules/mapping_qc.smk, because it also applies to user-supplied BAM files.