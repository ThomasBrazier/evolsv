rule samplot_subset_DUP:
    """
    Randomly subset SVs for diagnostic plot
    Subset N variants of each type DEL/DUP/INS/INV
    """
    input:
        final = "{wdir}/{genome}_filtered.vcf"
    output:
        subset_dup_tmp = temp("{wdir}/samplot/{genome}_samplot_DUP_tmp.vcf"),
        subset_dup = temp("{wdir}/samplot/{genome}_samplot_DUP.vcf")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools filter -i'INFO/SVTYPE="DUP"' {input.final} > {output.subset_dup_tmp}
        # Subset N random SNPs (default=100)
        N={config[n_samplot]}
        cat {output.subset_dup_tmp} | grep '^#' > {output.subset_dup}
        cat {output.subset_dup_tmp} | grep -v '^#' | shuf -n $N  >> {output.subset_dup}
        """

rule samplot_subset_INV:
    """
    Randomly subset SVs for diagnostic plot
    Subset N variants of each type DEL/DUP/INS/INV
    """
    input:
        final = "{wdir}/{genome}_filtered.vcf"
    output:
        subset_inv_tmp = temp("{wdir}/samplot/{genome}_samplot_INV_tmp.vcf"),
        subset_inv = temp("{wdir}/samplot/{genome}_samplot_INV.vcf")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools filter -i'INFO/SVTYPE="INV"' {input.final} > {output.subset_inv_tmp}
        # Subset N random SNPs (default=100)
        N={config[n_samplot]}
        cat {output.subset_inv_tmp} | grep '^#' > {output.subset_inv}
        cat {output.subset_inv_tmp} | grep -v '^#' | shuf -n $N  >> {output.subset_inv}
        """

# rule samplot_subset_INS:
#     """
#     Randomly subset SVs for diagnostic plot
#     Subset N variants of each type DEL/DUP/INS/INV
#     """
#     input:
#         final = "{wdir}/{genome}_filtered.vcf"
#     output:
#         subset_ins_tmp = temp("{wdir}/samplot/{genome}_samplot_INS_tmp.vcf"),
#         subset_ins = temp("{wdir}/samplot/{genome}_samplot_INS.vcf")
#     conda:
#         "../envs/bcftools.yaml"
#     shell:
#         """
#         bcftools filter -i'INFO/SVTYPE="INS"' {input.final} > {output.subset_ins_tmp}
#         # Subset N random SNPs (default=100)
#         N={config[n_samplot]}
#         cat {output.subset_ins_tmp} | grep '^#' > {output.subset_ins}
#         cat {output.subset_ins_tmp} | grep -v '^#' | shuf -n $N  >> {output.subset_ins}
#         """

rule samplot_subset_DEL:
    """
    Randomly subset SVs for diagnostic plot
    Subset N variants of each type DEL/DUP/INS/INV
    """
    input:
        final = "{wdir}/{genome}_filtered.vcf"
    output:
        subset_del_tmp = temp("{wdir}/samplot/{genome}_samplot_DEL_tmp.vcf"),
        subset_del = temp("{wdir}/samplot/{genome}_samplot_DEL.vcf")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools filter -i'INFO/SVTYPE="DEL"' {input.final} > {output.subset_del_tmp}
        # Subset N random SNPs (default=100)
        N={config[n_samplot]}
        cat {output.subset_del_tmp} | grep '^#' > {output.subset_del}
        cat {output.subset_del_tmp} | grep -v '^#' | shuf -n $N  >> {output.subset_del}
        """


rule samplot_plot:
    """
    Plot a random subset of SVs
    For diagnostic purpose
    """
    input:
        subset_DUP = expand("{wdir}/samplot/{genome}_samplot_DUP.vcf", wdir=wdir, genome=genome),
        subset_INV = expand("{wdir}/samplot/{genome}_samplot_INV.vcf", wdir=wdir, genome=genome),
        # subset_INS = expand("{wdir}/samplot/{genome}_samplot_INS.vcf", wdir=wdir, genome=genome),
        subset_DEL = expand("{wdir}/samplot/{genome}_samplot_DEL.vcf", wdir=wdir, genome=genome),
        fasta = expand("{wdir}/{genome}.fna", wdir=wdir, genome=genome),
        bam = expand("{wdir}/{genome}_sorted.bam", wdir=wdir, genome=genome)
    output:
        index_html_DUP = "{wdir}/samplot/DUP/index.html",
        index_html_INV = "{wdir}/samplot/INV/index.html",
        # index_html_INS = "{wdir}/samplot/INS/index.html",
        index_html_DEL = "{wdir}/samplot/DEL/index.html"
    conda:
        "../envs/samplot.yaml"
    params:
        outdir = "{wdir}/samplot"
    shell:
        """
        samplot vcf \
            --vcf {input.subset_DUP} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir}/DUP \
            -O jpg \
            --format GT,DP,AD,PL \
            -b {input.bam} \
            --sample_ids {sample_id} \
            --debug
        samplot vcf \
            --vcf {input.subset_INV} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir}/INV \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {sample_id} \
            -b {input.bam} \
            --debug
        # samplot vcf \
        #     --vcf {input.subset_INS} \
        #     --plot_all \
        #     --threads {threads} \
        #     -d {params.outdir}/INS \
        #     -O jpg \
        #     --format GT,DP,AD,PL \
        #     --sample_ids {sample_id} \
        #     -b {input.bam} \
        #     --debug
        samplot vcf \
            --vcf {input.subset_DEL} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir}/DEL \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {sample_id} \
            -b {input.bam} \
            --debug
        """


rule final_report:
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
        mappability = "{wdir}/callability/{genome}_callable_mappable.bed",
        index_html_DUP = "{wdir}/samplot/DUP/index.html",
        index_html_INV = "{wdir}/samplot/INV/index.html",
        # index_html_INS = "{wdir}/samplot/INS/index.html",
        index_html_DEL = "{wdir}/samplot/DEL/index.html"
    output:
        "{wdir}/{genome}_finalQC.html"
    conda:
        "../envs/Renv.yaml"
    shell:
        """
        Rscript workflow/scripts/finalQC.R {wdir} {genome}
        """


rule gzvcf:
    """
    BGzip final VCF
    """
    input:
        vcf = "{wdir}/{genome}_filtered.vcf",
        vcf_sexchr = "{wdir}/{genome}_filtered_sexchr.vcf",
        html = "{wdir}/{genome}_finalQC.html"
    output:
        tmp_vcf = temp("{wdir}/{genome}_filtered_newheader.vcf"),
        tmp_vcf_sexchr = temp("{wdir}/{genome}_filtered_sexchr_newheader.vcf"),
        vcf = "{wdir}/{genome}_filtered.vcf.gz",
        vcf_sexchr = "{wdir}/{genome}_filtered_sexchr.vcf.gz",
        vcf_idx = "{wdir}/{genome}_filtered.vcf.gz.csi",
        vcf_sexchr_idx = "{wdir}/{genome}_filtered_sexchr.vcf.gz.csi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools annotate --header-lines workflow/header/header.txt {input.vcf} > {output.tmp_vcf}
        bcftools annotate --header-lines workflow/header/header.txt {input.vcf_sexchr} > {output.tmp_vcf_sexchr}

        bcftools sort {output.tmp_vcf} -O v | bgzip > {output.vcf}
        bcftools sort {output.tmp_vcf_sexchr} -O v | bgzip > {output.vcf_sexchr}

        tabix --csi {output.vcf}
        tabix --csi {output.vcf_sexchr}
        """