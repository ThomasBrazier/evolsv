rule svjedigraph:
    """
    Use SVjedi-graph on the merged dataset to genotype SVs
    """
    input:
        merged = "{wdir}/merging/{genome}_merged.vcf",
        fasta = "{wdir}/genome/{genome}.fna",
        merged_fastq = "{wdir}/fastq/{genome}_filtered.fastq.gz",
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf = temp("{wdir}/genotype/{genome}_merged_genotype_tmp.vcf"),
        vcf_renamed = temp("{wdir}/genotype/{genome}_merged_genotype.vcf"),
        gfa = temp("{wdir}/genotype/{genome}_merged.gfa"),
        gaf = temp("{wdir}/genotype/{genome}_merged.gaf"),
        aln = "{wdir}/genotype/{genome}_merged_informative_aln.json"
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_svjedigraph.log"
    shell:
        """
        svjedi-graph.py -v {input.merged} -r {input.fasta} \
        -q {input.merged_fastq} -p {wdir}/genotype/{genome}_merged \
        -t {resources.cpus_per_task} --minsupport {config[svjedigraph_minsupport]}
        mv {output.vcf_renamed} {output.vcf}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf}
        """


rule autosomes_sexchromosomes:
    """
    Remove MT and Un chromosomes
    Get the names of sex chromosomes and autosomes to copy in separate vcf files
    """
    input: 
        seq = "{wdir}/genome/{genome}_sequence_report.jsonl"
    output: 
        sexchromosomes = "{wdir}/genome/{genome}.sexchromosomes",
        autosomes = "{wdir}/genome/{genome}.autosomes",
        chromosome_names = "{wdir}/genome/{genome}.chromosomes"
    conda:
        "../envs/Renv.yaml"
    params:
        seq = "{wdir}/genome/{genome}_sequence_report.jsonl",
        sexchromosomes = "{wdir}/genome/{genome}.sexchromosomes",
        autosomes = "{wdir}/genome/{genome}.autosomes",
        scaffolds_to_exclude = config["scaffolds_to_exclude"],
        chromosome_names = "{wdir}/genome/{genome}.chromosomes"
    script:
        "../scripts/autosomes_sexchromosomes.R"


rule all_samples_vcf:
    """
    Merge samples from Jasmine (_merged.vcf) and SVjedi-graph (_merged_genotype.vcf)
    """
    input:
        jasmine = "{wdir}/merging/{genome}_merged.vcf",
        svjedi = "{wdir}/genotype/{genome}_merged_genotype.vcf",
        genome_index="{wdir}/genome/{genome}.fna.fai"
    output:
        allsamples = temp("{wdir}/{genome}_allsamples.vcf"),
        jasmine_reheadered = temp("{wdir}/merging/{genome}_reheadered.vcf"),
        jasmine_gz = temp("{wdir}/merging/{genome}_merged.vcf.gz"),
        svjedi_reheadered = temp("{wdir}/genotype/{genome}_merged_genotype_reheadered.vcf"),
        svjedi_gz = temp("{wdir}/genotype/{genome}_merged_genotype.vcf.gz")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools reheader --fai {input.genome_index} {input.jasmine} > {output.jasmine_reheadered}

        bcftools annotate --header-lines workflow/header/header.txt -Ou {output.jasmine_reheadered} | \
        bcftools sort -Ou | \
        bcftools view -Oz -o {output.jasmine_gz}

        bcftools index {output.jasmine_gz}

        bcftools reheader --fai {input.genome_index} {input.svjedi} > {output.svjedi_reheadered}

        bcftools annotate --header-lines workflow/header/header.txt -Ou {output.svjedi_reheadered} | \
        bcftools sort -Ou | \
        bcftools view -Oz -o {output.svjedi_gz}

        bcftools index {output.svjedi_gz}

        bcftools merge --output {output.allsamples} {output.jasmine_gz} {output.svjedi_gz}
        """


rule final_vcf:
    """
    Do a last filtering step on the merged genotyped dataset
    Separate autosomes and sex chromosomes
    """
    input:
        vcf = "{wdir}/{genome}_allsamples.vcf",
        sexchromosomes = "{wdir}/genome/{genome}.sexchromosomes",
        autosomes = "{wdir}/genome/{genome}.autosomes"
    output:
        final_tmp = temp("{wdir}/{genome}_final_tmp.vcf"),
        final = "{wdir}/{genome}_final.vcf",
        final_sexchr = "{wdir}/{genome}_final_sexchr.vcf"   
    threads: workflow.cores
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/logs/{genome}_final_filtering.log"
    shell:
        """
        bcftools view -T {input.autosomes} -l 0 -o {output.final_tmp} {input.vcf}
        # Remove accessory INFO tags
        bcftools annotate --remove INFO/STRAND,INFO/STRANDS,INFO/AF,INFO/STDEV_POS,INFO/STDEV_LEN,INFO/COVERAGE,INFO/SUPPORT_INLINE,INFO/SUPPORT_LONG,INFO/CIPOS,INFO/CILEN,INFO/RE,INFO/PRECISE,INFO/IMPRECISE,INFO/REFINEDALT,INFO/MULTI,INFO/LARGEINS,INFO/STD_SPAN,INFO/STD_POS,INFO/STD_POS1,INFO/STD_POS2,INFO/CHR2,INFO/END,INFO/MAPQ,INFO/SUPPREAD,INFO/SUPPORT,INFO/READS,INFO/CUTPASTE,FILTER/MOSAIC_AF,FILTER/SVLEN_MIN,FILTER/STRAND,FILTER/ALN_NM,INFO/STARTVARIANCE,INFO/ENDVARIANCE,FILTER/COV_CHANGE_FRAC,FILTER/COV_CHANGE,FILTER/not_fully_covered,FILTER/hom_ref,INFO/SEQS,FILTER/q5,FILTER/STDEV_POS,FILTER/STDEV_LEN {output.final_tmp} > {output.final}

        # Filter sex chromosomes
        if [ $(wc -l < {input.sexchromosomes}) -eq 0 ]; then
        echo "No sex chr"
        cat {input.vcf} | grep "#" > {output.final_sexchr}
        else
        bcftools view -T {input.sexchromosomes} -l 0 -o {output.final_sexchr} {input.vcf}
        fi
        """



