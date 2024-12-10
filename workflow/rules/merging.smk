rule vcf_list:
    """
    Create the txt file containing the list of the VCF files to merge
    """
    input:
        sniffles = "{wdir}/sniffles_genotype/{genome}_sniffles_genotype.vcf",
        svim = "{wdir}/svim_genotype/{genome}_svim_genotype.vcf",
        cutesv = "{wdir}/cutesv_genotype/{genome}_cutesv_genotype.vcf",
        debreak = "{wdir}/debreak_genotype/{genome}_debreak_genotype.vcf",
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
        sniffles = "{wdir}/sniffles_genotype/{genome}_sniffles_genotype.vcf",
        svim = "{wdir}/svim_genotype/{genome}_svim_genotype.vcf",
        cutesv = "{wdir}/cutesv_genotype/{genome}_cutesv_genotype.vcf",
        debreak = "{wdir}/debreak_genotype/{genome}_debreak_genotype.vcf",
        fasta = "{wdir}/{genome}.fna"
    output:
        "{wdir}/{genome}_merged.vcf",
        # "{wdir}/iris_{genome}/resultsstore.txt",
        # "{wdir}/iris_{genome}/results.csv"
    threads: workflow.cores
    conda:
        "../envs/jasminesv.yaml"
    log:
        "{wdir}/logs/{genome}_jasmine.log"
    shell:
        """
        jasmine file_list={input.vcflist} out_file={output} genome_file={input.fasta} \
        out_dir={wdir}/jasmine bam_list={input.bamlist} \
        --output_genotypes \
        # --run_iris iris_args="out_dir=iris_{genome}" \
        --dup_to_ins
        """
