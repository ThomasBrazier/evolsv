rule samplot_subset_DUP:
    """
    Randomly subset SVs for diagnostic plot
    Subset N variants of each type DEL/DUP/INS/INV
    """
    input:
        final="{wdir}/{sample}/genotype/{genome}_merged_genotype.vcf.gz",
    output:
        subset_dup_tmp=temp("{wdir}/{sample}/samplot/{genome}_samplot_DUP_tmp.vcf"),
        subset_dup=temp("{wdir}/{sample}/samplot/{genome}_samplot_DUP.vcf"),
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}.samplot_subset_DUP.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.samplot_subset_DUP.tsv"
    shell:
        """
        bcftools filter -i'INFO/OLDTYPE="DUP"' {input.final} > {output.subset_dup_tmp}
        # Subset N random SNPs (default=100)
        N={config[n_samplot]}
        cat {output.subset_dup_tmp} | grep '^#' > {output.subset_dup}
        cat {output.subset_dup_tmp} | grep -v '^#' | shuf -n $N  >> {output.subset_dup} || true
        """


rule samplot_subset_INV:
    """
    Randomly subset SVs for diagnostic plot
    Subset N variants of each type DEL/DUP/INS/INV
    """
    input:
        final="{wdir}/{sample}/genotype/{genome}_merged_genotype.vcf.gz",
    output:
        subset_inv_tmp=temp("{wdir}/{sample}/samplot/{genome}_samplot_INV_tmp.vcf"),
        subset_inv=temp("{wdir}/{sample}/samplot/{genome}_samplot_INV.vcf"),
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}.samplot_subset_INV.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.samplot_subset_INV.tsv"
    shell:
        """
        bcftools filter -i'INFO/SVTYPE="INV"' {input.final} > {output.subset_inv_tmp}
        # Subset N random SNPs (default=100)
        N={config[n_samplot]}
        cat {output.subset_inv_tmp} | grep '^#' > {output.subset_inv}
        cat {output.subset_inv_tmp} | grep -v '^#' | shuf -n $N  >> {output.subset_inv} || true
        """


# rule samplot_subset_INS:
#     """
#     Randomly subset SVs for diagnostic plot
#     Subset N variants of each type DEL/DUP/INS/INV
#     """
#     input:
#         final = "{wdir}/{sample}/{genome}_final.vcf.gz"
#     output:
#         subset_ins_tmp = temp("{wdir}/{sample}/samplot/{genome}_samplot_INS_tmp.vcf"),
#         subset_ins = temp("{wdir}/{sample}/samplot/{genome}_samplot_INS.vcf")
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
        final="{wdir}/{sample}/genotype/{genome}_merged_genotype.vcf.gz",
    output:
        subset_del_tmp=temp("{wdir}/{sample}/samplot/{genome}_samplot_DEL_tmp.vcf"),
        subset_del=temp("{wdir}/{sample}/samplot/{genome}_samplot_DEL.vcf"),
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}.samplot_subset_DEL.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.samplot_subset_DEL.tsv"
    shell:
        """
        bcftools filter -i'INFO/SVTYPE="DEL"' {input.final} > {output.subset_del_tmp}
        # Subset N random SNPs (default=100)
        N={config[n_samplot]}
        cat {output.subset_del_tmp} | grep '^#' > {output.subset_del}
        cat {output.subset_del_tmp} | grep -v '^#' | shuf -n $N  >> {output.subset_del} || true
        """


rule samplot_plot:
    """
    Plot a random subset of SVs
    For diagnostic purpose
    """
    input:
        subset_DUP="{wdir}/{sample}/samplot/{genome}_samplot_DUP.vcf",
        subset_INV="{wdir}/{sample}/samplot/{genome}_samplot_INV.vcf",
        subset_DEL="{wdir}/{sample}/samplot/{genome}_samplot_DEL.vcf",
        fasta="{wdir}/genome/{genome}.fna",
        bam_minimap2="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam",
        bam_index_minimap2="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai",
        bam_ngmlr="{wdir}/{sample}/bam/{genome}_ngmlr_sorted.bam",
        bam_index_ngmlr="{wdir}/{sample}/bam/{genome}_ngmlr_sorted.bam.bai",
    output:
        "{wdir}/{sample}/samplot/minimap2_{genome}/DUP/index.html",
        "{wdir}/{sample}/samplot/minimap2_{genome}/INV/index.html",
        "{wdir}/{sample}/samplot/minimap2_{genome}/DEL/index.html",
        "{wdir}/{sample}/samplot/ngmlr_{genome}/DUP/index.html",
        "{wdir}/{sample}/samplot/ngmlr_{genome}/INV/index.html",
        "{wdir}/{sample}/samplot/ngmlr_{genome}/DEL/index.html",
    conda:
        "../envs/samplot.yaml"
    params:
        outdir_minimap2="{wdir}/{sample}/samplot/minimap2_{genome}",
        outdir_ngmlr="{wdir}/{sample}/samplot/ngmlr_{genome}",
    log:
        "{wdir}/{sample}/logs/{genome}.samplot_plot.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.samplot_plot.tsv"
    shell:
        """
        if [[ $(cat {input.subset_DUP} | grep -v '#' | wc -l) -eq 0 ]]; then
            echo "No duplication found."
            touch {wdir}/{wildcards.sample}/samplot/minimap2_{genome}/DUP/index.html
        else
            samplot vcf \
            --vcf {input.subset_DUP} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_minimap2}/DUP \
            -O jpg \
            --format GT,DP,AD,PL \
            -b {input.bam_minimap2} \
            --sample_ids {wildcards.sample} \
            --debug
        fi

        if [[ $(cat {input.subset_INV} | grep -v '#' | wc -l) -eq 0 ]]; then
            echo "No inversion found."
            touch {wdir}/{wildcards.sample}/samplot/minimap2_{genome}/INV/index.html
        else
            samplot vcf \
            --vcf {input.subset_INV} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_minimap2}/INV \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {wildcards.sample} \
            -b {input.bam_minimap2} \
            --debug
        fi
        
        if [[ $(cat {input.subset_DEL} | grep -v '#' | wc -l) -eq 0 ]]; then
            echo "No inversion found."
            touch {wdir}/{wildcards.sample}/samplot/minimap2_{genome}/DEL/index.html
        else
            samplot vcf \
            --vcf {input.subset_DEL} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_minimap2}/DEL \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {wildcards.sample} \
            -b {input.bam_minimap2} \
            --debug
        fi

        if [[ $(cat {input.subset_DUP} | grep -v '#' | wc -l) -eq 0 ]]; then
            echo "No duplication found."
            touch {wdir}/{wildcards.sample}/samplot/ngmlr_{genome}/DUP/index.html
        else
            samplot vcf \
            --vcf {input.subset_DUP} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_ngmlr}/DUP \
            -O jpg \
            --format GT,DP,AD,PL \
            -b {input.bam_ngmlr} \
            --sample_ids {wildcards.sample} \
            --debug
        fi

        if [[ $(cat {input.subset_INV} | grep -v '#' | wc -l) -eq 0 ]]; then
            echo "No inversion found."
            touch {wdir}/{wildcards.sample}/samplot/ngmlr_{genome}/INV/index.html
        else
            samplot vcf \
            --vcf {input.subset_INV} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_ngmlr}/INV \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {wildcards.sample} \
            -b {input.bam_ngmlr} \
            --debug
        fi

        if [[ $(cat {input.subset_DEL} | grep -v '#' | wc -l) -eq 0 ]]; then
            echo "No inversion found."
            touch {wdir}/{wildcards.sample}/samplot/ngmlr_{genome}/DEL/index.html
        else
            samplot vcf \
            --vcf {input.subset_DEL} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir_ngmlr}/DEL \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {wildcards.sample} \
            -b {input.bam_ngmlr} \
            --debug
        fi
        """


rule vcf_to_tsv:
    """
    Convert the full vcf to tabular data frame for R scripts
    """
    input:
        vcf="{wdir}/{sample}/{genome}_final.vcf",
    output:
        tsv="{wdir}/{sample}/{genome}_final.tsv",
    log:
        "{wdir}/{sample}/logs/{genome}.vcf_to_tsv.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.vcf_to_tsv.tsv"
    shell:
        """
        # All samples
        bash workflow/scripts/vcf_to_tsv.sh {input.vcf} {output.tsv}
        """


rule merging_qc:
    """
    Run test on the merged VCF to find errors due to parsing or processing in Jasmine
    Or inconsistencies between callers
    """
    input:
        tsv="{wdir}/{sample}/{genome}_final.tsv",
    output:
        "{wdir}/{sample}/merging_QC/{genome}_svlen_equal_zero.tsv",
        "{wdir}/{sample}/merging_QC/{genome}_avglen_equal_zero.tsv",
        "{wdir}/{sample}/merging_QC/{genome}_no_avgend_field.tsv",
        "{wdir}/{sample}/merging_QC/{genome}_unmerged_sv.tsv",
    conda:
        "../envs/Renv.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}.merging_qc.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.merging_qc.tsv"
    shell:
        """
        Rscript workflow/scripts/merging_qc.R {wdir}/{wildcards.sample} {genome}
        """


rule vcf_to_tsv_tools:
    """
    Convert the full vcf to tabular data frame for R scripts
    """
    input:
        minimap2_cutesv_vcf="{wdir}/{sample}/calling/{genome}_minimap2_cutesv.vcf",
        minimap2_svim_vcf="{wdir}/{sample}/calling/{genome}_minimap2_svim.vcf",
        minimap2_sniffles_vcf="{wdir}/{sample}/calling/{genome}_minimap2_sniffles.vcf",
        minimap2_debreak_vcf="{wdir}/{sample}/calling/{genome}_minimap2_debreak.vcf",
        ngmlr_cutesv_vcf="{wdir}/{sample}/calling/{genome}_ngmlr_cutesv.vcf",
        ngmlr_svim_vcf="{wdir}/{sample}/calling/{genome}_ngmlr_svim.vcf",
        ngmlr_sniffles_vcf="{wdir}/{sample}/calling/{genome}_ngmlr_sniffles.vcf",
        ngmlr_debreak_vcf="{wdir}/{sample}/calling/{genome}_ngmlr_debreak.vcf",
    output:
        minimap2_cutesv_tsv="{wdir}/{sample}/calling/{genome}_minimap2_cutesv.tsv",
        minimap2_svim_tsv="{wdir}/{sample}/calling/{genome}_minimap2_svim.tsv",
        minimap2_sniffles_tsv="{wdir}/{sample}/calling/{genome}_minimap2_sniffles.tsv",
        minimap2_debreak_tsv="{wdir}/{sample}/calling/{genome}_minimap2_debreak.tsv",
        ngmlr_cutesv_tsv="{wdir}/{sample}/calling/{genome}_ngmlr_cutesv.tsv",
        ngmlr_svim_tsv="{wdir}/{sample}/calling/{genome}_ngmlr_svim.tsv",
        ngmlr_sniffles_tsv="{wdir}/{sample}/calling/{genome}_ngmlr_sniffles.tsv",
        ngmlr_debreak_tsv="{wdir}/{sample}/calling/{genome}_ngmlr_debreak.tsv",
    log:
        "{wdir}/{sample}/logs/{genome}.vcf_to_tsv_tools.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.vcf_to_tsv_tools.tsv"
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
        vcf="{wdir}/{sample}/{genome}_final.vcf",
        merged="{wdir}/{sample}/genotype/{genome}_merged_genotype.vcf",
        tsv="{wdir}/{sample}/{genome}_final.tsv",
        svlen_equal_zero="{wdir}/{sample}/merging_QC/{genome}_svlen_equal_zero.tsv",
        avglen_equal_zero="{wdir}/{sample}/merging_QC/{genome}_avglen_equal_zero.tsv",
        no_avgend_field="{wdir}/{sample}/merging_QC/{genome}_no_avgend_field.tsv",
        unmerged_sv="{wdir}/{sample}/merging_QC/{genome}_unmerged_sv.tsv",
        minimap2_cutesv_tsv="{wdir}/{sample}/calling/{genome}_minimap2_cutesv.tsv",
        minimap2_svim_tsv="{wdir}/{sample}/calling/{genome}_minimap2_svim.tsv",
        minimap2_sniffles_tsv="{wdir}/{sample}/calling/{genome}_minimap2_sniffles.tsv",
        minimap2_debreak_tsv="{wdir}/{sample}/calling/{genome}_minimap2_debreak.tsv",
        ngmlr_cutesv_tsv="{wdir}/{sample}/calling/{genome}_ngmlr_cutesv.tsv",
        ngmlr_svim_tsv="{wdir}/{sample}/calling/{genome}_ngmlr_svim.tsv",
        ngmlr_sniffles_tsv="{wdir}/{sample}/calling/{genome}_ngmlr_sniffles.tsv",
        ngmlr_debreak_tsv="{wdir}/{sample}/calling/{genome}_ngmlr_debreak.tsv",
        mapping_minimap2="{wdir}/{sample}/mapping_QC/{genome}_minimap2_mapping.stats.tsv",
        mapping_ngmlr="{wdir}/{sample}/mapping_QC/{genome}_ngmlr_mapping.stats.tsv",
        # Read by finalQC.Rmd. Absent only with a local sequence_report and no local
        # assembly_data_report; the report then skips its assembly section.
        assembly_report=(
            ["{wdir}/genome/{genome}_assembly_data_report.jsonl"]
            if reference_source != "local" or assembly_data_report
            else []
        ),
    output:
        html="{wdir}/{sample}/{genome}_finalQC.html",
        # "{wdir}/{sample}/{genome}_finalQC.pdf",
        performance="{wdir}/{sample}/performance/{genome}_performance.tsv",
    conda:
        "../envs/Renv.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}.final_report.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.final_report.tsv"
    shell:
        """
        if [[ $(cat {input.vcf} | grep -v '#' | wc -l) -eq 0 ]]; then
            echo "No SV found."
            touch {output.html}
            touch {output.performance}
        else
        Rscript workflow/scripts/finalQC.R {wdir} {genome} {wildcards.sample}
        fi
        """


rule light_vcf:
    """
    Make a light VCF for faster computation
    No sequences
    """
    input:
        vcf="{wdir}/{sample}/{genome}_final.vcf",
        vcf_sexchr="{wdir}/{sample}/{genome}_final_sexchr.vcf",
        html="{wdir}/{sample}/{genome}_finalQC.html",
    output:
        light_vcf="{wdir}/{sample}/{genome}_final_light.vcf",
    conda:
        "../envs/pysam_v2.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}.light_vcf.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.light_vcf.tsv"
    shell:
        """
        # Iterate over the VCF to add symbolic type to REF/ALT field
        # instead of sequence
        python workflow/scripts/light_vcf.py {input.vcf} {output.light_vcf}
        """


rule gzvcf:
    """
    BGzip final VCF
    """
    input:
        vcf="{wdir}/{sample}/{genome}_final.vcf",
        light_vcf="{wdir}/{sample}/{genome}_final_light.vcf",
        vcf_sexchr="{wdir}/{sample}/{genome}_final_sexchr.vcf",
        html="{wdir}/{sample}/{genome}_finalQC.html",
    output:
        tmp_vcf=temp("{wdir}/{sample}/{genome}_final_newheader.vcf"),
        tmp_vcf_sexchr=temp("{wdir}/{sample}/{genome}_final_sexchr_newheader.vcf"),
        vcf="{wdir}/{sample}/{genome}_final.vcf.gz",
        vcf_sexchr="{wdir}/{sample}/{genome}_final_sexchr.vcf.gz",
        vcf_idx="{wdir}/{sample}/{genome}_final.vcf.gz.csi",
        vcf_sexchr_idx="{wdir}/{sample}/{genome}_final_sexchr.vcf.gz.csi",
        light_tmp_vcf=temp("{wdir}/{sample}/{genome}_final_newheader_light.vcf"),
        light_vcf="{wdir}/{sample}/{genome}_final_light.vcf.gz",
        light_vcf_idx="{wdir}/{sample}/{genome}_final_light.vcf.gz.csi",
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}.gzvcf.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.gzvcf.tsv"
    shell:
        """
        bcftools annotate --header-lines workflow/header/header.txt {input.vcf} > {output.tmp_vcf}
        bcftools annotate --header-lines workflow/header/header.txt {input.vcf_sexchr} > {output.tmp_vcf_sexchr}
        bcftools annotate --header-lines workflow/header/header.txt {input.light_vcf} > {output.light_tmp_vcf}

        bcftools sort {output.tmp_vcf} -O v | bgzip > {output.vcf}
        bcftools sort {output.tmp_vcf_sexchr} -O v | bgzip > {output.vcf_sexchr}
        bcftools sort {output.light_tmp_vcf} -O v | bgzip > {output.light_vcf}

        tabix --csi {output.vcf}
        tabix --csi {output.vcf_sexchr}
        tabix --csi {output.light_vcf}
        """
