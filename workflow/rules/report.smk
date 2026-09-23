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

    One job per aligner and SV type. The six blocks this replaces were identical up to
    the aligner and the SV type, and enumerating the aligners made a single-aligner run
    (config key `aligners`) impossible.

    An empty subset VCF is not an error: rules samplot_subset_* draw at most
    config[n_samplot] calls of that type and a callset may hold none. The index.html is
    then created empty, so the rule still satisfies rule all.
    """
    wildcard_constraints:
        # Without this, {svtype} would also match the `DUP_tmp` subset files.
        svtype="DUP|INV|DEL",
    input:
        subset="{wdir}/{sample}/samplot/{genome}_samplot_{svtype}.vcf",
        fasta="{wdir}/genome/{genome}.fna",
        bam="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam",
        bam_index="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam.bai",
    output:
        index="{wdir}/{sample}/samplot/{aligner}_{genome}/{svtype}/index.html",
    conda:
        "../envs/samplot.yaml"
    params:
        outdir="{wdir}/{sample}/samplot/{aligner}_{genome}/{svtype}",
    log:
        "{wdir}/{sample}/logs/{genome}_{aligner}_{svtype}.samplot_plot.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_{aligner}_{svtype}.samplot_plot.tsv"
    shell:
        """
        if [[ $(cat {input.subset} | grep -v '#' | wc -l) -eq 0 ]]; then
            echo "No {wildcards.svtype} found." > {log}
            touch {output.index}
        else
            samplot vcf \
            --vcf {input.subset} \
            --plot_all \
            --threads {threads} \
            -d {params.outdir} \
            -O jpg \
            --format GT,DP,AD,PL \
            --sample_ids {wildcards.sample} \
            -b {input.bam} \
            --debug 2> {log}
        fi
        """


rule vcf_to_tsv:
    """
    Convert the full vcf to tabular data frame for R scripts

    The callset names are passed in rather than hardcoded in the script: they name the
    VCF sample columns, and their order is the order of Jasmine's SUPP_VEC bits (see
    `callsets` in workflow/Snakefile). merging_qc.R and finalQC.Rmd both read the shape
    of the run out of this header.
    """
    input:
        vcf="{wdir}/{sample}/{genome}_final.vcf",
    output:
        tsv="{wdir}/{sample}/{genome}_final.tsv",
    params:
        callsets=",".join(f"{aligner}_{caller}" for aligner, caller in callsets),
    log:
        "{wdir}/{sample}/logs/{genome}.vcf_to_tsv.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.vcf_to_tsv.tsv"
    shell:
        """
        # All samples
        bash workflow/scripts/vcf_to_tsv.sh {input.vcf} {output.tsv} {params.callsets}
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
        vcfs=expand(
            "{{wdir}}/{{sample}}/calling/{{genome}}_{aligner}_{caller}.vcf",
            zip,
            aligner=[aligner for aligner, _ in callsets],
            caller=[caller for _, caller in callsets],
        ),
    output:
        tsvs=expand(
            "{{wdir}}/{{sample}}/calling/{{genome}}_{aligner}_{caller}.tsv",
            zip,
            aligner=[aligner for aligner, _ in callsets],
            caller=[caller for _, caller in callsets],
        ),
    log:
        "{wdir}/{sample}/logs/{genome}.vcf_to_tsv_tools.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.vcf_to_tsv_tools.tsv"
    shell:
        """
        # Tool specific output. Both lists are built from `callsets`, so the nth VCF and
        # the nth TSV are the same callset.
        vcfs=({input.vcfs})
        tsvs=({output.tsvs})
        for i in "${{!vcfs[@]}}"; do
            bash workflow/scripts/vcf_to_tsv_tools.sh "${{vcfs[$i]}}" "${{tsvs[$i]}}" 2>> {log}
        done
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
        caller_tsvs=expand(
            "{{wdir}}/{{sample}}/calling/{{genome}}_{aligner}_{caller}.tsv",
            zip,
            aligner=[aligner for aligner, _ in callsets],
            caller=[caller for _, caller in callsets],
        ),
        mapping_stats=expand(
            "{{wdir}}/{{sample}}/mapping_QC/{{genome}}_{aligner}_mapping.stats.tsv",
            aligner=aligners,
        ),
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
