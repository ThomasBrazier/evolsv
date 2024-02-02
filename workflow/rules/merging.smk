rule vcf_list:
    """
    Create the txt file containing the list of the VCF files to merge
    """
    input:
        "{wdir}/{sra}_cutesv.vcf",
        "{wdir}/{sra}_sniffles.vcf",
        "{wdir}/{sra}_svim_noBND.vcf"
    output:
        "{wdir}/{sra}_vcf_list.txt"
    shell:
        """
        echo '{wdir}/{sra}_svim_noBND.vcf' > {wdir}/{sra}_vcf_list.txt
        echo '{wdir}/{sra}_sniffles.vcf' >> {wdir}/{sra}_vcf_list.txt
        echo '{wdir}/{sra}_cutesv.vcf' >> {wdir}/{sra}_vcf_list.txt
        """

rule bam_list:
    """
    Create the txt file containing the list of the bam files for IRIS
    """
    input :
        "{wdir}/{sra}_vcf_list.txt"
    output:
        "{wdir}/{sra}_bam_list.txt"
    shell:
        """
        echo '{wdir}/{sra}_sorted.bam' > {wdir}/{sra}_bam_list.txt
        echo '{wdir}/{sra}_sorted.bam' >> {wdir}/{sra}_bam_list.txt
        echo '{wdir}/{sra}_sorted.bam' >> {wdir}/{sra}_bam_list.txt
        """

rule jasmine:
    """
    Merge the VCF files obtained by the three SV callers
    """
    input:
        "{wdir}/{sra}_bam_list.txt"
    output:
        "{wdir}/{sra}_merged.vcf"
    conda:
        "Envs/jasminesv.yaml"
    shell:
        "Scripts/jasmine.sh {sra}"