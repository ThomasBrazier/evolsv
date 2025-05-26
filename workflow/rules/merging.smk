rule jasmine:
    """
    Merge the VCF files obtained by the three SV callers
    """
    input:
        sniffles_minimap2 = "{wdir}/filtered/{genome}_minimap2_sniffles_filtered.vcf",
        svim_minimap2 = "{wdir}/filtered/{genome}_minimap2_svim_filtered.vcf",
        cutesv_minimap2 = "{wdir}/filtered/{genome}_minimap2_cutesv_filtered.vcf",
        debreak_minimap2 = "{wdir}/filtered/{genome}_minimap2_debreak_filtered.vcf",
        sniffles_ngmlr = "{wdir}/filtered/{genome}_ngmlr_sniffles_filtered.vcf",
        svim_ngmlr = "{wdir}/filtered/{genome}_ngmlr_svim_filtered.vcf",
        cutesv_ngmlr = "{wdir}/filtered/{genome}_ngmlr_cutesv_filtered.vcf",
        debreak_ngmlr = "{wdir}/filtered/{genome}_ngmlr_debreak_filtered.vcf",
        fasta = "{wdir}/genome/{genome}.fna"
    output:
        tempvcf = "{wdir}/jasmine/{genome}_merged_noGenotypes.vcf",
        vcf = "{wdir}/merging/{genome}_merged.vcf",
        vcflist = "{wdir}/merging/{genome}_vcf_list.txt",
        bamlist = "{wdir}/merging/{genome}_bam_list.txt"
    conda:
        "../envs/jasminesv.yaml"
    shell:
        """
        # Make sure local decimal point is '.'
        LC_NUMERIC=C
        export LC_NUMERIC

        locale decimal_point

        LANG=en_US
        export LANG

        echo "{wdir}/filtered/{genome}_minimap2_sniffles_filtered.vcf" > {output.vcflist}
        echo "{wdir}/filtered/{genome}_minimap2_svim_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/filtered/{genome}_minimap2_cutesv_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/filtered/{genome}_minimap2_debreak_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/filtered/{genome}_ngmlr_sniffles_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/filtered/{genome}_ngmlr_svim_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/filtered/{genome}_ngmlr_cutesv_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/filtered/{genome}_ngmlr_debreak_filtered.vcf" >> {output.vcflist}

        echo "{wdir}/bam/{genome}_minimap2_sorted.bam" > {output.bamlist}
        echo "{wdir}/bam/{genome}_ngmlr_sorted.bam" >> {output.bamlist}

        jasmine file_list={output.vcflist} \
        out_file={output.vcf} genome_file={input.fasta} \
        out_dir={wdir}/jasmine bam_list={output.bamlist} \
        --ignore_strand --max_dist {config[jasmine_max_dist]} \
        --output_genotypes
        """
