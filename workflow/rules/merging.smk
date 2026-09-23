rule jasmine:
    """
    Merge the VCF files obtained by the SV callers

    One callset per aligner + caller pair, in `callsets` order (see workflow/Snakefile):
    Jasmine numbers the SUPP_VEC bits of a merged call by the position of the callset in
    file_list, and workflow/scripts/vcf_to_tsv.sh names the VCF sample columns in that
    same order. Writing the list from {input.vcfs} keeps the two from drifting.

    With a single aligner this merges 4 callsets instead of 8, so SUPP counts are out of
    4 and the calls no longer carry cross-aligner agreement (see README).
    """
    input:
        vcfs=expand(
            "{{wdir}}/{{sample}}/filtered/{{genome}}_{aligner}_{caller}_filtered.vcf",
            zip,
            aligner=[aligner for aligner, _ in callsets],
            caller=[caller for _, caller in callsets],
        ),
        bams=expand(
            "{{wdir}}/{{sample}}/bam/{{genome}}_{aligner}_sorted.bam", aligner=aligners
        ),
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
    params:
        aligners=" ".join(aligners),
        callers=" ".join(CALLERS),
    log:
        "{wdir}/{sample}/logs/{genome}.jasmine.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.jasmine.tsv"
    shell:
        """
        # Make sure local decimal point is '.'
        LC_NUMERIC=C
        export LC_NUMERIC

        locale decimal_point

        LANG=en_US
        export LANG

        # One line per callset, in {input.vcfs} order: that order is what SUPP_VEC bits
        # refer to downstream.
        printf '%s\\n' {input.vcfs} > {output.vcflist}

        printf '%s\\n' {input.bams} > {output.bamlist}

        # Modify header to prevent missing contig
        for aligner in {params.aligners}; do
        for tool in {params.callers}; do
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
