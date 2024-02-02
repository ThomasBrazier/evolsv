rule vcf_list:
    """
    Create the txt file containing the list of the VCF files to merge
    """
    input:
        sniffles = "{wdir}/{sra}_sniffles.vcf",
        svim = "{wdir}/{sra}_svim.vcf",
        nobnd = "$sra/${sra}_svim_noBND.vcf",
        cutesv = "{wdir}/{sra}_cutesv.vcf"
    output:
        "{wdir}/{sra}_vcf_list.txt"
    shell:
        """
        echo '{input.nobnd}' > {output}
        echo '{input.sniffles}' >> {output}
        echo '{input.cutesv}' >> {output}
        """

rule bam_list:
    """
    Create the txt file containing the list of the bam files for IRIS
    """
    input :
        "{wdir}/{sra}_vcf_list.txt",
        bam = "{wdir}/{sra}_sorted.bam"
    output:
        "{wdir}/{sra}_bam_list.txt"
    shell:
        """
        echo '{input.bam}' > {output}
        echo '{input.bam}' >> {output}
        echo '{input.bam}' >> {output}
        """

rule jasmine:
    """
    Merge the VCF files obtained by the three SV callers
    """
    input:
        bamlist = "{wdir}/{sra}_bam_list.txt",
        vcflist = "{wdir}/{sra}_vcf_list.txt",
        sniffles = "{wdir}/{sra}_sniffles.vcf",
        svim = "{wdir}/{sra}_svim.vcf",
        nobnd = "$sra/${sra}_svim_noBND.vcf",
        cutesv = "{wdir}/{sra}_cutesv.vcf",
        fasta = "{wdir}/{sra}.fna"
    output:
        "{wdir}/{sra}_merged.vcf"
    conda:
        "workflow/envs/jasminesv.yaml"
    shell:
        """
        jasmine file_list={input.vcflist} out_file={output} genome_file={input.fasta} \
        out_dir={wdir} bam_list={input.bamlist} --output_genotypes --run_iris --dup_to_ins
        """