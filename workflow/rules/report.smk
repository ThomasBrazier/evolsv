rule finalreport:
    """
    Compute and print a summary report for assembly, mapping, SV calling, merging and genotyping
    """
    input:
        final = "{wdir}/{sra}_merged_genotype.vcf",
        sniffles = "{wdir}/{sra}_sniffles_noBND.vcf",
        svim = "{wdir}/{sra}_svim_noBND.vcf",
        cutesv = "{wdir}/{sra}_cutesv_noBND.vcf",
        debreak = "{wdir}/{sra}_debreak_noBND.vcf"
    output:
        "{wdir}/{sra}_finalQC.html"
    conda:
        "../envs/Renv.yaml"
    shell:
        """
        Rscript workflow/scripts/finalQC.R {wdir} {sra}
        """