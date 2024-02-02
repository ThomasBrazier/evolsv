rule svim:
    """
    Call SVs with the SVIM tool
    """
    input:
        "{wdir}/{sra}_sorted.bam.bai"
    output:
        "{wdir}/{sra}_svim.vcf"
    conda:
        "Envs/svim.yaml"
    shell:
        "Scripts/svim.sh {sra} {config[min_sv_size]} {config[min_coverage]} {config[svim_quality]}"


rule sniffles:
    """
    Call SVs with the Sniffles tool
    """
    input:
        "{wdir}/{sra}_svim.vcf"
    output:
        "{wdir}/{sra}_sniffles.vcf"
    conda:
        "Envs/sniffles.yaml"
    shell:
        "Scripts/sniffles.sh {sra} {config[min_sv_size]} {config[min_coverage]}"


rule cutesv:
    """
    Call SVs with the Cutesv tool
    """
    input:
        "{wdir}/{sra}_sniffles.vcf"
    output:
        "{wdir}/{sra}_cutesv.vcf"
    conda:
        "Envs/cutesv.yaml"
    shell:
        "Scripts/cutesv.sh {sra} {config[min_sv_size]} {config[min_coverage]}"