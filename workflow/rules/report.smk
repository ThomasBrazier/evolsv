rule samplot_subset_DUP:
    """
    Randomly subset SVs for diagnostic plot
    Subset N variants of each type DEL/DUP/INS/INV
    """
    input:
        final = "{wdir}/{genome}_merged_genotype.vcf"
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
        final = "{wdir}/{genome}_merged_genotype.vcf"
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
#         final = "{wdir}/{genome}_final.vcf"
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
        final = "{wdir}/{genome}_merged_genotype.vcf"
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
        subset_DEL = expand("{wdir}/samplot/{genome}_samplot_DEL.vcf", wdir=wdir, genome=genome),
        fasta = expand("{wdir}/{genome}.fna", wdir=wdir, genome=genome),
        bam_minimap2 = "{wdir}/{genome}_minimap2_sorted.bam",
        bam_index_minimap2 = "{wdir}/{genome}_minimap2_sorted.bam.bai",
        bam_ngmlr = "{wdir}/{genome}_ngmlr_sorted.bam",
        bam_index_ngmlr = "{wdir}/{genome}_ngmlr_sorted.bam.bai"
    output:
        "{wdir}/{genome}_samplot_minimap2/DUP/index.html",
        "{wdir}/{genome}_samplot_minimap2/INV/index.html",
        "{wdir}/{genome}_samplot_minimap2/DEL/index.html",
        "{wdir}/{genome}_samplot_ngmlr/DUP/index.html",
        "{wdir}/{genome}_samplot_ngmlr/INV/index.html",
        "{wdir}/{genome}_samplot_ngmlr/DEL/index.html"
    conda:
        "../envs/samplot.yaml"
    params:
        outdir_minimap2 = "{wdir}/{genome}_samplot_minimap2",
        outdir_ngmlr = "{wdir}/{genome}_samplot_ngmlr"
    shell:
        """
        samplot vcf \
            --vcf {input.subset_DUP} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_minimap2}/DUP \
            -O jpg \
            --format GT,DP,AD,PL \
            -b {input.bam_minimap2} \
            --sample_ids {sample_id} \
            --debug
        samplot vcf \
            --vcf {input.subset_INV} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_minimap2}/INV \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {sample_id} \
            -b {input.bam_minimap2} \
            --debug
        samplot vcf \
            --vcf {input.subset_DEL} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_minimap2}/DEL \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {sample_id} \
            -b {input.bam_minimap2} \
            --debug

        samplot vcf \
            --vcf {input.subset_DUP} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_ngmlr}/DUP \
            -O jpg \
            --format GT,DP,AD,PL \
            -b {input.bam_ngmlr} \
            --sample_ids {sample_id} \
            --debug
        samplot vcf \
            --vcf {input.subset_INV} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_ngmlr}/INV \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {sample_id} \
            -b {input.bam_ngmlr} \
            --debug
        samplot vcf \
            --vcf {input.subset_DEL} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_ngmlr}/DEL \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {sample_id} \
            -b {input.bam_ngmlr} \
            --debug
        """


rule vcf_to_tsv:
    """
    Convert the full vcf to tabular data frame for R scripts
    """
    input:
        vcf = "{wdir}/{genome}_final.vcf"
    output:
        tsv = "{wdir}/{genome}_final.tsv"
    shell:
        """
        # All samples
        bash workflow/scripts/vcf_to_tsv.sh {input.vcf} {output.tsv}
        """


rule vcf_to_tsv_tools:
    """
    Convert the full vcf to tabular data frame for R scripts
    """
    input:
        minimap2_cutesv_vcf = "{wdir}/{genome}_minimap2_cutesv.vcf",
        minimap2_svim_vcf = "{wdir}/{genome}_minimap2_svim.vcf",
        minimap2_sniffles_vcf = "{wdir}/{genome}_minimap2_sniffles.vcf",
        minimap2_debreak_vcf = "{wdir}/{genome}_minimap2_debreak.vcf",
        ngmlr_cutesv_vcf = "{wdir}/{genome}_ngmlr_cutesv.vcf",
        ngmlr_svim_vcf = "{wdir}/{genome}_ngmlr_svim.vcf",
        ngmlr_sniffles_vcf = "{wdir}/{genome}_ngmlr_sniffles.vcf",
        ngmlr_debreak_vcf = "{wdir}/{genome}_ngmlr_debreak.vcf"
    output:
        minimap2_cutesv_tsv = "{wdir}/{genome}_minimap2_cutesv.tsv",
        minimap2_svim_tsv = "{wdir}/{genome}_minimap2_svim.tsv",
        minimap2_sniffles_tsv = "{wdir}/{genome}_minimap2_sniffles.tsv",
        minimap2_debreak_tsv = "{wdir}/{genome}_minimap2_debreak.tsv",
        ngmlr_cutesv_tsv = "{wdir}/{genome}_ngmlr_cutesv.tsv",
        ngmlr_svim_tsv = "{wdir}/{genome}_ngmlr_svim.tsv",
        ngmlr_sniffles_tsv = "{wdir}/{genome}_ngmlr_sniffles.tsv",
        ngmlr_debreak_tsv = "{wdir}/{genome}_ngmlr_debreak.tsv"
    shell:
        """
        # Tool specific output
        bash workflow/scripts/vcf_to_tsv_tools.sh {input.minimap2_cutesv_vcf} {output.minimap2_cutesv_tsv}
        bash workflow/scripts/vcf_to_tsv_tools.sh {input.minimap2_svim_vcf} {output.minimap2_svim_tsv}
        bash workflow/scripts/vcf_to_tsv_tools.sh {input.minimap2_sniffles_vcf} {output.minimap2_sniffles_tsv}
        bash workflow/scripts/vcf_to_tsv_tools.sh {input.minimap2_debreak_vcf} {output.minimap2_debreak_tsv}
        
        bash workflow/scripts/vcf_to_tsv_tools.sh {input.ngmlr_cutesv_vcf} {output.ngmlr_cutesv_tsv}
        bash workflow/scripts/vcf_to_tsv_tools.sh {input.ngmlr_svim_vcf} {output.ngmlr_svim_tsv}
        bash workflow/scripts/vcf_to_tsv_tools.sh {input.ngmlr_sniffles_vcf} {output.ngmlr_sniffles_tsv}
        bash workflow/scripts/vcf_to_tsv_tools.sh {input.ngmlr_debreak_vcf} {output.ngmlr_debreak_tsv}
        """


rule final_report:
    """
    Compute and print a summary report for assembly, mapping, SV calling, merging and genotyping
    """
    input:
        vcf = "{wdir}/{genome}_final.vcf",
        merged = "{wdir}/{genome}_merged_genotype.vcf",
        tsv = "{wdir}/{genome}_final.tsv",
        minimap2_cutesv_tsv = "{wdir}/{genome}_minimap2_cutesv.tsv",
        minimap2_svim_tsv = "{wdir}/{genome}_minimap2_svim.tsv",
        minimap2_sniffles_tsv = "{wdir}/{genome}_minimap2_sniffles.tsv",
        minimap2_debreak_tsv = "{wdir}/{genome}_minimap2_debreak.tsv",
        ngmlr_cutesv_tsv = "{wdir}/{genome}_ngmlr_cutesv.tsv",
        ngmlr_svim_tsv = "{wdir}/{genome}_ngmlr_svim.tsv",
        ngmlr_sniffles_tsv = "{wdir}/{genome}_ngmlr_sniffles.tsv",
        ngmlr_debreak_tsv = "{wdir}/{genome}_ngmlr_debreak.tsv",
        mapping_minimap2 = "{wdir}/mapping/{genome}_minimap2_mapping.stats.tsv",
        mapping_ngmlr = "{wdir}/mapping/{genome}_ngmlr_mapping.stats.tsv"
    output:
        "{wdir}/{genome}_finalQC.html",
        "{wdir}/performance/{genome}_performance.tsv"
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
        vcf = "{wdir}/{genome}_final.vcf",
        vcf_sexchr = "{wdir}/{genome}_final_sexchr.vcf",
        html = "{wdir}/{genome}_finalQC.html"
    output:
        tmp_vcf = temp("{wdir}/{genome}_final_newheader.vcf"),
        tmp_vcf_sexchr = temp("{wdir}/{genome}_final_sexchr_newheader.vcf"),
        vcf = "{wdir}/{genome}_final.vcf.gz",
        vcf_sexchr = "{wdir}/{genome}_final_sexchr.vcf.gz",
        vcf_idx = "{wdir}/{genome}_final.vcf.gz.csi",
        vcf_sexchr_idx = "{wdir}/{genome}_final_sexchr.vcf.gz.csi"
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