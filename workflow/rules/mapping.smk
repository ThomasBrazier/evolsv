rule minimap2:
    """
    Map reads to the reference genome with Minimap2
    """
    input:
        fastq = "{wdir}/fastq/{genome}.fastq.gz",
        fasta = "{wdir}/{genome}.fna",
        html = expand("{wdir}/fastqc/{sample}_fastqc.html", wdir=wdir, sample=samples["sra"]),
        qczip = expand("{wdir}/fastqc/{sample}_fastqc.zip", wdir=wdir, sample=samples["sra"]),
        nanoplot = expand("{wdir}/nanoplot/{sample}_NanoStats.txt", wdir=wdir, sample=samples["sra"])
    output:
        sam = temp("{wdir}/{genome}_minimap2.sam")
    conda:
        "../envs/minimap2.yaml"
    log:
        "{wdir}/logs/{genome}_minimap2.log"
    shell:
        """
        minimap2 -ax {config[minimap_ax]} --MD -2 \
        --seed {config[minimap_seed]} --eqx \
        -t {resources.cpus_per_task} \
        -R "@RG\\tID:{sample_id}\\tSM:{sample_id}" \
        --sam-hit-only {input.fasta} {input.fastq} > {output.sam}
        """



rule nglmr:
    """
    Map reads to the reference genome with NGLMR
    --rg-id <string>
        Adds RG:Z:<string> to all alignments in SAM/BAM [none]
    --rg-sm <string>
        RG header: Sample [none]
    -t <int>,  --threads <int>
        Number of threads [1]
    -x <pacbio, ont>,  --presets <pacbio, ont>
        Parameter presets for different sequencing technologies [pacbio]
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
        fastq = "{wdir}/fastq/{genome}.fastq.gz",
        fasta = "{wdir}/{genome}.fna",
        html = expand("{wdir}/fastqc/{sample}_fastqc.html", wdir=wdir, sample=samples["sra"]),
        qczip = expand("{wdir}/fastqc/{sample}_fastqc.zip", wdir=wdir, sample=samples["sra"]),
        nanoplot = expand("{wdir}/nanoplot/{sample}_NanoStats.txt", wdir=wdir, sample=samples["sra"])
    output:
        sam = temp("{wdir}/{genome}_nglmr.sam")
    conda:
        "../envs/nglmr.yaml"
    log:
        "{wdir}/logs/{genome}_nglmr.log"
    shell:
        """
        ngmlr -t {resources.cpus_per_task} \
        -r {input.fasta} -q {input.fastq} \
        --presets pacbio --min-identity {config[min-identity]} \
        --min-residues {config[min-residues]} --no-smallinv {config[no-smallinv]} \
        --no-lowqualitysplit {config[no-lowqualitysplit]} \
        --rg-id {sample_id} --rg-sm {sample_id} \
        -o {output.sam}
        """


rule samtools_view:
    """
    Transform the sam file to a bam file
    """
    input:
        expand("{wdir}/{genome}_{aligner}.sam", wdir=wdir, genome=genome, aligner=aligner)
    output:
        expand("{wdir}/{genome}_{aligner}.bam", wdir=wdir, genome=genome, aligner=aligner)
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -S -b {input} > {output}
        """


rule samtools_sort:
    """
    Sort the bam file
    """
    input:
        expand("{wdir}/{genome}_{aligner}.bam", wdir=wdir, genome=genome, aligner=aligner)
    output:
        expand("{wdir}/{genome}_{aligner}_sorted.bam", wdir=wdir, genome=genome, aligner=aligner)
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
        expand("{wdir}/{genome}_{aligner}_sorted.bam", wdir=wdir, genome=genome, aligner=aligner)
    output:
        expand("{wdir}/{genome}_{aligner}_sorted.bam.bai", wdir=wdir, genome=genome, aligner=aligner)
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input}
        """


rule samtools_stats:
    """
    Mapping QC
    """
    input:
        bam = expand("{wdir}/{genome}_{aligner}_sorted.bam", wdir=wdir, genome=genome, aligner=aligner),
        bai = expand("{wdir}/{genome}_{aligner}_sorted.bam.bai", wdir=wdir, genome=genome, aligner=aligner)
    output:
        stats = expand("{wdir}/mapping/{genome}_{aligner}_mapping.stats", wdir=wdir, genome=genome, aligner=aligner),
        stattsv = expand("{wdir}/mapping/{genome}_{aligner}_mapping.stats.tsv", wdir=wdir, genome=genome, aligner=aligner),
        plot = expand("{wdir}/mapping/{genome}_{aligner}_mapping_plot.html", wdir=wdir, genome=genome, aligner=aligner)
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p {wdir}/mapping
        samtools stats {input.bam} > {output.stats}
        cat {output.stats} | grep ^SN | cut -f 2- > {output.stattsv}
        # QC visualization
        plot-bamstats -p {wdir}/mapping/{genome}_mapping_plot {output.stats}
        """
