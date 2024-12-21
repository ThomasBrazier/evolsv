# rule vcf_list:
#     """
#     Create the txt file containing the list of the VCF files to merge
#     """
#     input:
#         sniffles = "{wdir}/sniffles_genotype/{genome}_{aligner}_sniffles_genotype.vcf",
#         svim = "{wdir}/svim_genotype/{genome}_{aligner}_svim_genotype.vcf",
#         cutesv = "{wdir}/cutesv_genotype/{genome}_{aligner}_cutesv_genotype.vcf",
#         debreak = "{wdir}/debreak_genotype/{genome}_{aligner}_debreak_genotype.vcf"
#     output:
#         vcflist = "{wdir}/{genome}_vcf_list.txt",
#         bamlist = "{wdir}/{genome}_bam_list.txt"
#     shell:
#         """
#         echo "{wdir}/sniffles_genotype/{genome}_minimap2_sniffles_genotype.vcf" > {output.vcflist}
#         echo "{wdir}/sniffles_genotype/{genome}_minimap2_svim_genotype.vcf" >> {output.vcflist}
#         echo "{wdir}/sniffles_genotype/{genome}_minimap2_cutesv_genotype.vcf" >> {output.vcflist}
#         echo "{wdir}/sniffles_genotype/{genome}_minimap2_debreak_genotype.vcf" >> {output.vcflist}
#         echo "{wdir}/sniffles_genotype/{genome}_ngmlr_sniffles_genotype.vcf" >> {output.vcflist}
#         echo "{wdir}/sniffles_genotype/{genome}_ngmlr_svim_genotype.vcf" >> {output.vcflist}
#         echo "{wdir}/sniffles_genotype/{genome}_ngmlr_cutesv_genotype.vcf" >> {output.vcflist}
#         echo "{wdir}/sniffles_genotype/{genome}_ngmlr_debreak_genotype.vcf" >> {output.vcflist}

#         echo "{wdir}/{genome}_minimap2_sorted.bam" > {output.bamlist}
#         echo "{wdir}/{genome}_ngmlr_sorted.bam" >> {output.bamlist}
#         """

# rule bam_list:
#     """
#     Create the txt file containing the list of the bam files for IRIS
#     """
#     input :
#         bam1 = "{wdir}/{genome}_minimap2_sorted.bam",
#         bam2 = "{wdir}/{genome}_ngmlr_sorted.bam"
#     output:
#         bamlist = "{wdir}/{genome}_bam_list.txt"
#     shell:
#         """
#         echo "{wdir}/{genome}_minimap2_sorted.bam" > {output.bamlist}
#         echo "{wdir}/{genome}_ngmlr_sorted.bam" >> {output.bamlist}
#         """


rule jasmine:
    """
    Merge the VCF files obtained by the three SV callers
    """
    input:
        sniffles_minimap2 = "{wdir}/sniffles_genotype/{genome}_minimap2_sniffles_genotype.vcf",
        svim_minimap2 = "{wdir}/svim_genotype/{genome}_minimap2_svim_genotype.vcf",
        cutesv_minimap2 = "{wdir}/cutesv_genotype/{genome}_minimap2_cutesv_genotype.vcf",
        debreak_minimap2 = "{wdir}/debreak_genotype/{genome}_minimap2_debreak_genotype.vcf",
        sniffles_ngmlr = "{wdir}/sniffles_genotype/{genome}_ngmlr_sniffles_genotype.vcf",
        svim_ngmlr = "{wdir}/svim_genotype/{genome}_ngmlr_svim_genotype.vcf",
        cutesv_ngmlr = "{wdir}/cutesv_genotype/{genome}_ngmlr_cutesv_genotype.vcf",
        debreak_ngmlr = "{wdir}/debreak_genotype/{genome}_ngmlr_debreak_genotype.vcf",
        fasta = "{wdir}/{genome}.fna"
    output:
        vcf = "{wdir}/{genome}_merged.vcf",
        vcflist = "{wdir}/{genome}_vcf_list.txt",
        bamlist = "{wdir}/{genome}_bam_list.txt"
    conda:
        "../envs/jasminesv.yaml"
    log:
        "{wdir}/logs/{genome}_jasmine.log"
    shell:
        """
        echo "{wdir}/sniffles_genotype/{genome}_minimap2_sniffles_genotype.vcf" > {output.vcflist}
        echo "{wdir}/sniffles_genotype/{genome}_minimap2_svim_genotype.vcf" >> {output.vcflist}
        echo "{wdir}/sniffles_genotype/{genome}_minimap2_cutesv_genotype.vcf" >> {output.vcflist}
        echo "{wdir}/sniffles_genotype/{genome}_minimap2_debreak_genotype.vcf" >> {output.vcflist}
        echo "{wdir}/sniffles_genotype/{genome}_ngmlr_sniffles_genotype.vcf" >> {output.vcflist}
        echo "{wdir}/sniffles_genotype/{genome}_ngmlr_svim_genotype.vcf" >> {output.vcflist}
        echo "{wdir}/sniffles_genotype/{genome}_ngmlr_cutesv_genotype.vcf" >> {output.vcflist}
        echo "{wdir}/sniffles_genotype/{genome}_ngmlr_debreak_genotype.vcf" >> {output.vcflist}

        echo "{wdir}/{genome}_minimap2_sorted.bam" > {output.bamlist}
        echo "{wdir}/{genome}_ngmlr_sorted.bam" >> {output.bamlist}

        jasmine file_list={output.vcflist} out_file={output.vcf} genome_file={input.fasta} \
        out_dir={wdir}/jasmine bam_list={output.bamlist} \
        --output_genotypes \
        --dup_to_ins
        """
