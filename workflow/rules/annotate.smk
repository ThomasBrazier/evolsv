rule bwa_index:
    """
    BWA index of the shared reference, built once for all individuals.
    truvari anno grm remaps k-mers with BWA and needs this index next to the FASTA.
    """
    input:
        fasta="{wdir}/genome/{genome}.fna",
    output:
        amb=temp("{wdir}/genome/{genome}.fna.amb"),
        ann=temp("{wdir}/genome/{genome}.fna.ann"),
        bwt=temp("{wdir}/genome/{genome}.fna.bwt"),
        pac=temp("{wdir}/genome/{genome}.fna.pac"),
        sa=temp("{wdir}/genome/{genome}.fna.sa"),
    conda:
        "../envs/truvari.yaml"
    log:
        "{wdir}/logs/bwa_index/{genome}.txt",
    benchmark:
        "{wdir}/genome/benchmarks/{genome}.bwa_index.tsv"
    shell:
        """
        bwa index {input.fasta} 2> {log}
        """


rule truvari_grm:
    """
    truvari anno grm
    For every SV, we create a kmer over the the upstream and downstream reference and alternate breakpoints.
    We then remap that kmer to the reference genome and report alignment information.
    This does not alter the VCF traditional annotations,
    but instead will create a pandas DataFrame and save it to a joblib object.
    """
    input:
        vcf="{wdir}/{sample}/{genome}_final.vcf.gz",
        fasta="{wdir}/genome/{genome}.fna",
        bwa_index=multiext(
            "{wdir}/genome/{genome}.fna", ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    output:
        grm_pandas="{wdir}/{sample}/annotate_grm/{genome}_grm.jl",
        tabix="{wdir}/{sample}/{genome}_final.vcf.gz.tbi",
    conda:
        "../envs/truvari.yaml"
    params:
        kmersize=config["grm_kmersize"],
        min_sv_size=config["min_sv_size"],
    log:
        "{wdir}/{sample}/logs/truvari_grm/{genome}.txt",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.truvari_grm.tsv"
    shell:
        """
        tabix {input.vcf} 2> {log}
        singularity exec workflow/containers/truvari.sif truvari anno grm \
        -i {input.vcf} -r {input.fasta} -o {output.grm_pandas} \
        -k {params.kmersize} -m {params.min_sv_size} -t {resources.cpus_per_task} 2>> {log}
        """


# rule truvari_repeatmasker:
#     """
#     truvari anno repmask
#     Wrapper around RepeatMasker to annotate insertion sequences in a VCF.
#     """
#     input:
#         vcf = "{wdir}/{sample}/{genome}_final.vcf",
#         fasta = "{wdir}/{sample}/{genome}.fna"
#     output:
#         vcf = "{wdir}/{sample}/{genome}_repmasked.vcf"
#     conda:
#         "../envs/truvari.yaml"
#     shell:
#         """
#         truvari anno repmask -i {input.vcf} -o {output.vcf} [-e EXECUTABLE] [-m MIN_LENGTH] [-M MAX_LENGTH] [-t THRESHOLD] [-p PARAMS] [-T THREADS]
#         """
