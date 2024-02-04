
rule sniffles2plot:
    """
    The sniffles2-lot package output a set of QC summary plots for a single VCF
    """
    input:
        vcf = "{wdir}/{sra}_merged_genotype.vcf"
    output:
        "{wdir}/vcf_QC/"
    conda:
        "../envs/sniffles.yaml"
    shell:
        """
        python3 -m sniffles2_plot -i {vcf} -o {wdir}/vcf_QC/
        """