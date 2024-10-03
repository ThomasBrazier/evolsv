rule samplot_subset:
    """
    Randomly subset SVs for diagnostic plot
    """
    input:
        final = "{wdir}/{genome}_filtered.vcf"
    output:
        subset = "{wdir}/samplot/{genome}_subset.vcf"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        cat {input.final} | grep '^#' > {output.subset}
        cat {input.final} | grep -v '^#' | shuf -n {config[n_samplot]}  >> {output.subset}
        """

rule samplot_plot:
    """
    Plot a random subset of SVs
    For diagnostic purpose
    """
    input:
        subsetvcf = "{wdir}/samplot/{genome}_subset.vcf",
        fasta = "{wdir}/{genome}.fna",
        bam = "{wdir}/{genome}_sorted.bam"
    output:
        command_file = "{wdir}/samplot/{genome}_commands.sh"
    conda:
        "../envs/samplot.yaml"
    params:
        outdir = "{wdir}/samplot/"
    shell:
        """
        samplot vcf \
            --vcf {input.subsetvcf} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir} \
            -O jpg \
            -b {input.bam} \
            --command_file {output.command_file}
        """


rule finalreport:
    """
    Compute and print a summary report for assembly, mapping, SV calling, merging and genotyping
    """
    input:
        final = "{wdir}/{genome}_filtered.vcf",
        merged = "{wdir}/{genome}_merged_genotype.vcf",
        sniffles = "{wdir}/{genome}_sniffles_noBND.vcf",
        svim = "{wdir}/{genome}_svim_noBND.vcf",
        cutesv = "{wdir}/{genome}_cutesv_noBND.vcf",
        debreak = "{wdir}/{genome}_debreak_noBND.vcf",
        mappability = "{wdir}/callability/{genome}_callable_mappable.bed"
    output:
        "{wdir}/{genome}_finalQC.html"
    conda:
        "../envs/Renv.yaml"
    shell:
        """
        Rscript workflow/scripts/finalQC.R {wdir} {genome}
        """


# rule gzvcf:
#     """
#     BGzip final VCF
#     """
#     input:
#         vcf = "{wdir}/{genome}_filtered.vcf",
#         vcf_sexchr = "{wdir}/{genome}_filtered_sexchr.vcf",
#         html = "{wdir}/{genome}_finalQC.html"
#     output:
#         vcf = "{wdir}/{genome}_filtered.vcf.gz",
#         vcf_sexchr = "{wdir}/{genome}_filtered_sexchr.vcf.gz",
#         vcf_idx = "{wdir}/{genome}_filtered.vcf.gz.csi",
#         vcf_sexchr_idx = "{wdir}/{genome}_filtered_sexchr.vcf.gz.csi"
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         """
#         bcftools sort {input.vcf} -O v | bgzip > {output.vcf}
#         bcftools sort {input.vcf_sexchr} -O v | bgzip > {output.vcf_sexchr}
#         # bgzip --keep --force --threads {threads} {input.vcf}
#         # bgzip --keep --force --threads {threads} {input.vcf_sexchr}
#         tabix --csi {output.vcf}
#         tabix --csi {output.vcf_sexchr}
#         """