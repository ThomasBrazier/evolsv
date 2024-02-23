rule vcf_list:
    """
    Create the txt file containing the list of the VCF files to merge
    """
    input:
        sniffles = "{wdir}/{genome}_sniffles_noBND.vcf",
        svim = "{wdir}/{genome}_svim_noBND.vcf",
        cutesv = "{wdir}/{genome}_cutesv_noBND.vcf",
        debreak = "{wdir}/{genome}_debreak_noBND.vcf",
        snifflesQC = "{wdir}/{genome}_sniffles_QC/variant_count.jpg",
        svimQC = "{wdir}/{genome}_svim_QC/variant_count.jpg",
        cutesvQC = "{wdir}/{genome}_cutesv_QC/variant_count.jpg"
    output:
        "{wdir}/{genome}_vcf_list.txt"
    shell:
        """
        echo '{input.svim}' > {output}
        echo '{input.sniffles}' >> {output}
        echo '{input.cutesv}' >> {output}
        echo '{input.debreak}' >> {output}
        """

rule bam_list:
    """
    Create the txt file containing the list of the bam files for IRIS
    """
    input :
        vcflist = "{wdir}/{genome}_vcf_list.txt",
        bam = "{wdir}/{genome}_sorted.bam"
    output:
        "{wdir}/{genome}_bam_list.txt"
    shell:
        """
        echo '{input.bam}' > {output}
        echo '{input.bam}' >> {output}
        echo '{input.bam}' >> {output}
        echo '{input.bam}' >> {output}
        """


rule jasmine:
    """
    Merge the VCF files obtained by the three SV callers
    """
    input:
        bamlist = "{wdir}/{genome}_bam_list.txt",
        vcflist = "{wdir}/{genome}_vcf_list.txt",
        sniffles = "{wdir}/{genome}_sniffles_noBND.vcf",
        svim = "{wdir}/{genome}_svim_noBND.vcf",
        cutesv = "{wdir}/{genome}_cutesv_noBND.vcf",
        debreak = "{wdir}/{genome}_debreak_noBND.vcf",
        fasta = "{wdir}/{genome}.fna"
    output:
        "{wdir}/{genome}_merged.vcf"
    threads: workflow.cores
    conda:
        "../envs/jasminesv.yaml"
    log:
        "{wdir}/logs/{genome}_jasmine.log"
    shell:
        """
        jasmine file_list={input.vcflist} out_file={output} genome_file={input.fasta} \
        out_dir={wdir} bam_list={input.bamlist} --output_genotypes --run_iris --dup_to_ins
        """
