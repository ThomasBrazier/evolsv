rule truvari_grm:
    """
    truvari anno grm
    For every SV, we create a kmer over the the upstream and downstream reference and alternate breakpoints.
    We then remap that kmer to the reference genome and report alignment information.
    This does not alter the VCF traditional annotations,
    but instead will create a pandas DataFrame and save it to a joblib object.
    """
    input:
        vcf = "{wdir}/{genome}_final.vcf.gz",
        fasta = "{wdir}/genome/{genome}.fna"
    output:
        grm_pandas = "{wdir}/annotate_grm/{genome}_grm.jl",
        tabix = "{wdir}/{genome}_final.vcf.gz.tbi",
        amb = temp("{wdir}/genome/{genome}.fna.amb"),
        ann = temp("{wdir}/genome/{genome}.fna.ann"),
        bwt = temp("{wdir}/genome/{genome}.fna.bwt"),
        pac = temp("{wdir}/genome/{genome}.fna.pac"),
        sa = temp("{wdir}/genome/{genome}.fna.sa")
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        bwa index {input.fasta} # Reference must be indexed
        tabix {input.vcf} # Tabix .ybi is required by truvari
        singularity exec workflow/containers/truvari.sif truvari anno grm \
        -i {input.vcf} -r {input.fasta} -o {output.grm_pandas} \
        -k {config[grm_kmersize]} -m {config[min_sv_size]} -t {resources.cpus_per_task}
        """


# rule truvari_repeatmasker:
#     """
#     truvari anno repmask
#     Wrapper around RepeatMasker to annotate insertion sequences in a VCF.
#     """
#     input:
#         vcf = "{wdir}/{genome}_final.vcf",
#         fasta = "{wdir}/{genome}.fna"
#     output:
#         vcf = "{wdir}/{genome}_repmasked.vcf"
#     conda:
#         "../envs/truvari.yaml"
#     shell:
#         """
#         truvari anno repmask -i {input.vcf} -o {output.vcf} [-e EXECUTABLE] [-m MIN_LENGTH] [-M MAX_LENGTH] [-t THRESHOLD] [-p PARAMS] [-T THREADS]
#         """