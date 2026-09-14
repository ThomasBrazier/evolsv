rule jasmine:
    """
    Merge the VCF files obtained by the three SV callers
    """
    input:
        sniffles_minimap2="{wdir}/{sample}/filtered/{genome}_minimap2_sniffles_filtered.vcf",
        svim_minimap2="{wdir}/{sample}/filtered/{genome}_minimap2_svim_filtered.vcf",
        cutesv_minimap2="{wdir}/{sample}/filtered/{genome}_minimap2_cutesv_filtered.vcf",
        debreak_minimap2="{wdir}/{sample}/filtered/{genome}_minimap2_debreak_filtered.vcf",
        sniffles_ngmlr="{wdir}/{sample}/filtered/{genome}_ngmlr_sniffles_filtered.vcf",
        svim_ngmlr="{wdir}/{sample}/filtered/{genome}_ngmlr_svim_filtered.vcf",
        cutesv_ngmlr="{wdir}/{sample}/filtered/{genome}_ngmlr_cutesv_filtered.vcf",
        debreak_ngmlr="{wdir}/{sample}/filtered/{genome}_ngmlr_debreak_filtered.vcf",
        fasta="{wdir}/genome/{genome}.fna",
        fasta_fai="{wdir}/genome/{genome}.fna.fai",
    output:
        tempvcf=temp("{wdir}/{sample}/jasmine/{genome}_merged_noGenotypes.vcf"),
        vcf="{wdir}/{sample}/merging/{genome}_merged.vcf",
        vcf_annotated=temp("{wdir}/{sample}/merging/{genome}_annotated.vcf"),
        vcfgz="{wdir}/{sample}/merging/{genome}_merged.vcf.gz",
        vcftabix="{wdir}/{sample}/merging/{genome}_merged.vcf.gz.tbi",
        vcflist="{wdir}/{sample}/merging/{genome}_vcf_list.txt",
        bamlist="{wdir}/{sample}/merging/{genome}_bam_list.txt",
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

        echo "{wdir}/{wildcards.sample}/filtered/{genome}_minimap2_sniffles_filtered.vcf" > {output.vcflist}
        echo "{wdir}/{wildcards.sample}/filtered/{genome}_minimap2_svim_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/{wildcards.sample}/filtered/{genome}_minimap2_cutesv_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/{wildcards.sample}/filtered/{genome}_minimap2_debreak_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/{wildcards.sample}/filtered/{genome}_ngmlr_sniffles_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/{wildcards.sample}/filtered/{genome}_ngmlr_svim_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/{wildcards.sample}/filtered/{genome}_ngmlr_cutesv_filtered.vcf" >> {output.vcflist}
        echo "{wdir}/{wildcards.sample}/filtered/{genome}_ngmlr_debreak_filtered.vcf" >> {output.vcflist}

        echo "{wdir}/{wildcards.sample}/bam/{genome}_minimap2_sorted.bam" > {output.bamlist}
        echo "{wdir}/{wildcards.sample}/bam/{genome}_ngmlr_sorted.bam" >> {output.bamlist}

        # Modify header to prevent missing contig
        for aligner in minimap2 ngmlr; do
        for tool in sniffles svim cutesv debreak; do
        bcftools reheader --fai {wdir}/genome/{genome}.fna.fai -o {wdir}/{wildcards.sample}/filtered/{genome}_${{aligner}}_${{tool}}_filtered.reheadered.vcf {wdir}/{wildcards.sample}/filtered/{genome}_${{aligner}}_${{tool}}_filtered.vcf
        rm {wdir}/{wildcards.sample}/filtered/{genome}_${{aligner}}_${{tool}}_filtered.vcf
        mv {wdir}/{wildcards.sample}/filtered/{genome}_${{aligner}}_${{tool}}_filtered.reheadered.vcf {wdir}/{wildcards.sample}/filtered/{genome}_${{aligner}}_${{tool}}_filtered.vcf
        done 
        done 

        jasmine file_list={output.vcflist} \
        out_file={output.vcf} genome_file={input.fasta} \
        out_dir={wdir}/{wildcards.sample}/jasmine bam_list={output.bamlist} \
        threads={threads} \
        --ignore_strand --max_dist {config[jasmine_max_dist]} \
        --output_genotypes

        # Re-annotate header and fields, then sort
        bcftools annotate --header-lines workflow/header/header.txt -o {output.vcf_annotated} {output.vcf} 
        bcftools sort -o {output.vcf} {output.vcf_annotated}
        
        bgzip --force --keep {output.vcf}
        tabix {output.vcfgz}

        #  --allow_intrasample raises an error - no bugfix in JasmineSV
        # see https://github.com/mkirsche/Jasmine/issues/58
        """
