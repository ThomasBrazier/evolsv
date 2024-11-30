rule svjedigraph:
    """
    Use SVjedi-graph on the merged dataset to genotype SVs
    """
    input:
        merged = "{wdir}/{genome}_merged.vcf",
        fasta = "{wdir}/{genome}.fna",
        fastq = expand("{wdir}/fastq/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"])
    output:
        "{wdir}/{genome}_merged_genotype.vcf",
        "{wdir}/{genome}_merged.gfa",
        "{wdir}/{genome}_merged.gaf",
        "{wdir}/{genome}_merged_informative_aln.json"
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_svjedigraph.log"
    shell:
        """
        svjedi-graph.py -v {input.merged} -r {input.fasta} -q {fqlist} -p {wdir}/{genome}_merged -t {resources.cpus_per_task} --minsupport {config[minsupport]}
        """


rule autosomes_sexchromosomes:
    """
    Remove MT and Un chromosomes
    Get the names of sex chromosomes and autosomes to copy in separate vcf files
    """
    input: 
        seq = "{wdir}/{genome}_sequence_report.jsonl"
    output: 
        sexchromosomes = "{wdir}/{genome}.sexchromosomes",
        autosomes = "{wdir}/{genome}.autosomes",
        chromosome_names = "{wdir}/{genome}.chromosomes"
    conda:
        "../envs/Renv.yaml"
    params:
        seq = "{wdir}/{genome}_sequence_report.jsonl",
        sexchromosomes = "{wdir}/{genome}.sexchromosomes",
        autosomes = "{wdir}/{genome}.autosomes",
        scaffolds_to_exclude = config["scaffolds_to_exclude"],
        chromosome_names = "{wdir}/{genome}.chromosomes"
    script:
        "../scripts/autosomes_sexchromosomes.R"


rule final_filtering:
    """
    Do a last filtering step on the merged genotyped dataset
    Separate autosomes and sex chromosomes
    """
    input:
        vcf = "{wdir}/{genome}_merged_genotype.vcf",
        sexchromosomes = "{wdir}/{genome}.sexchromosomes",
        autosomes = "{wdir}/{genome}.autosomes"
    output:
        filtered = "{wdir}/{genome}_filtered.vcf",
        filtered_sexchr = "{wdir}/{genome}_filtered_sexchr.vcf"   
    threads: workflow.cores
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/logs/{genome}_svjedigraph.log"
    shell:
        """
        bcftools view -T {input.autosomes} -l 0 -o {output.filtered} {input.vcf}
        bcftools view -T {input.sexchromosomes} -l 0 -o {output.filtered_sexchr} {input.vcf}
        """



