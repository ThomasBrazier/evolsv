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
        N={config["n_samplot"]}
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
        N={config["n_samplot"]}
        cat {output.subset_inv_tmp} | grep '^#' > {output.subset_inv}
        cat {output.subset_inv_tmp} | grep -v '^#' | shuf -n $N  >> {output.subset_inv}
        """

rule samplot_subset_INS:
    """
    Randomly subset SVs for diagnostic plot
    Subset N variants of each type DEL/DUP/INS/INV
    """
    input:
        final = "{wdir}/{genome}_filtered.vcf"
    output:
        subset_ins_tmp = temp("{wdir}/samplot/{genome}_samplot_INS_tmp.vcf"),
        subset_ins = temp("{wdir}/samplot/{genome}_samplot_INS.vcf")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools filter -i'INFO/SVTYPE="INS"' {input.final} > {output.subset_ins_tmp}
        # Subset N random SNPs (default=100)
        N={config["n_samplot"]}
        cat {output.subset_ins_tmp} | grep '^#' > {output.subset_ins}
        cat {output.subset_ins_tmp} | grep -v '^#' | shuf -n $N  >> {output.subset_ins}
        """

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
        N={config["n_samplot"]}
        cat {output.subset_del_tmp} | grep '^#' > {output.subset_del}
        cat {output.subset_del_tmp} | grep -v '^#' | shuf -n $N  >> {output.subset_del}
        """


rule samplot_plot:
    """
    Plot a random subset of SVs
    For diagnostic purpose
    """
    input:
        subset_DUP = "{wdir}/samplot/{genome}_samplot_DUP.vcf",
        subset_INV = "{wdir}/samplot/{genome}_samplot_INV.vcf",
        subset_INS = "{wdir}/samplot/{genome}_samplot_INS.vcf",
        subset_DEL = "{wdir}/samplot/{genome}_samplot_DEL.vcf",
        fasta = "{wdir}/{genome}.fna",
        bam = "{wdir}/{genome}_sorted.bam",
        gff3 = "{wdir}/{genome}.gff"
    output:
        command_file_DUP = "{wdir}/samplot/{genome}_commands_DUP.sh",command_file_INV = "{wdir}/samplot/{genome}_commands_INV.sh",
        command_file_INS = "{wdir}/samplot/{genome}_commands_INS.sh",
        command_file_DEL = "{wdir}/samplot/{genome}_commands_DEL.sh"
    conda:
        "../envs/samplot.yaml"
    params:
        outdir = "{wdir}/samplot/"
    shell:
        """
        samplot vcf \
            --vcf {input.subset_DUP} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir} \
            -O jpg \
            -b {input.bam} \
            --gff3 {gff3} \
            --command_file {output.command_file_DUP}
        samplot vcf \
            --vcf {input.subset_INV} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir} \
            -O jpg \
            -b {input.bam} \
            --gff3 {gff3} \
            --command_file {output.command_file_INV}
        samplot vcf \
            --vcf {input.subset_INS} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir} \
            -O jpg \
            -b {input.bam} \
            --gff3 {gff3} \
            --command_file {output.command_file_INS}
        samplot vcf \
            --vcf {input.subset_DEL} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir} \
            -O jpg \
            -b {input.bam} \
            --gff3 {gff3} \
            --command_file {output.command_file_DEL}
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
        command_file_DUP = "{wdir}/samplot/{genome}_commands_DUP.sh",command_file_INV = "{wdir}/samplot/{genome}_commands_INV.sh",
        command_file_INS = "{wdir}/samplot/{genome}_commands_INS.sh",
        command_file_DEL = "{wdir}/samplot/{genome}_commands_DEL.sh"
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
        vcf = "{wdir}/{genome}_filtered.vcf.gz",
        vcf_sexchr = "{wdir}/{genome}_filtered_sexchr.vcf.gz",
        vcf_idx = "{wdir}/{genome}_filtered.vcf.gz.csi",
        vcf_sexchr_idx = "{wdir}/{genome}_filtered_sexchr.vcf.gz.csi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools sort {input.vcf} -O v | bgzip > {output.vcf}
        bcftools sort {input.vcf_sexchr} -O v | bgzip > {output.vcf_sexchr}

        tabix --csi {output.vcf}
        tabix --csi {output.vcf_sexchr}
        """