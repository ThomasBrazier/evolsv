rule vcf_list:
    """
    Create the txt file containing the list of the VCF files to merge
    """
    input:
        sniffles_minimap2 = "{wdir}/sniffles_genotype/{genome}_minimap2_sniffles_genotype.vcf",
        svim_minimap2 = "{wdir}/svim_genotype/{genome}_minimap2_svim_genotype.vcf",
        cutesv_minimap2 = "{wdir}/cutesv_genotype/{genome}_minimap2_cutesv_genotype.vcf",
        debreak_minimap2 = "{wdir}/debreak_genotype/{genome}_minimap2_debreak_genotype.vcf",
        sniffles_nglmr = "{wdir}/sniffles_genotype/{genome}_nglmr_sniffles_genotype.vcf",
        svim_nglmr = "{wdir}/svim_genotype/{genome}_nglmr_svim_genotype.vcf",
        cutesv_nglmr = "{wdir}/cutesv_genotype/{genome}_nglmr_cutesv_genotype.vcf",
        debreak_nglmr = "{wdir}/debreak_genotype/{genome}_nglmr_debreak_genotype.vcf",
        snifflesQC = expand("{wdir}/{genome}_{aligner}_sniffles_QC/variant_count.jpg", wdir=wdir, genome=genome, aligner=aligner),
        svimQC = expand("{wdir}/{genome}_{aligner}_svim_QC/variant_count.jpg", wdir=wdir, genome=genome, aligner=aligner),
        cutesvQC = expand("{wdir}/{genome}_{aligner}_cutesv_QC/variant_count.jpg", wdir=wdir, genome=genome, aligner=aligner)
    output:
        vcflist = "{wdir}/{genome}_vcf_list.txt"
    shell:
        """
        echo '{input.svim_minimap2}' > {output.vcflist}
        echo '{input.sniffles_minimap2}' >> {output.vcflist}
        echo '{input.cutesv_minimap2}' >> {output.vcflist}
        echo '{input.debreak_minimap2}' >> {output.vcflist}
        echo '{input.svim_nglmr}' >> {output.vcflist}
        echo '{input.sniffles_nglmr}' >> {output.vcflist}
        echo '{input.cutesv_nglmr}' >> {output.vcflist}
        echo '{input.debreak_nglmr}' >> {output.vcflist}
        """

rule bam_list:
    """
    Create the txt file containing the list of the bam files for IRIS
    """
    input :
        bam_minimap2 = "{wdir}/{genome}_minimap2_sorted.bam",
        bam_nglmr = "{wdir}/{genome}_nglmr_sorted.bam"
    output:
        bamlist = "{wdir}/{genome}_bam_list.txt"
    shell:
        """
        echo '{input.bam_minimap2}' > {output.bamlist}
        echo '{input.bam_nglmr}' >> {output.bamlist}
        """


rule jasmine:
    """
    Merge the VCF files obtained by the three SV callers
    """
    input:
        bamlist = "{wdir}/{genome}_bam_list.txt",
        vcflist = "{wdir}/{genome}_vcf_list.txt",
        sniffles = expand("{wdir}/sniffles_genotype/{genome}_{aligner}_sniffles_genotype.vcf", wdir=wdir, genome=genome, aligner=aligner),
        svim = expand("{wdir}/svim_genotype/{genome}_{aligner}_svim_genotype.vcf", wdir=wdir, genome=genome, aligner=aligner),
        cutesv = expand("{wdir}/cutesv_genotype/{genome}_{aligner}_cutesv_genotype.vcf", wdir=wdir, genome=genome, aligner=aligner),
        debreak = expand("{wdir}/debreak_genotype/{genome}_{aligner}_debreak_genotype.vcf", wdir=wdir, genome=genome, aligner=aligner),
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
