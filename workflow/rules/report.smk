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
        debreak = "{wdir}/{genome}_debreak_noBND.vcf"
    output:
        "{wdir}/{genome}_finalQC.html"
    conda:
        "../envs/Renv.yaml"
    shell:
        """
        Rscript workflow/scripts/finalQC.R {wdir} {genome}
        """
