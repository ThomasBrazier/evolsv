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
    threads: workflow.cores
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_svjedigraph.log"
    shell:
        """
        svjedi-graph.py -v {input.merged} -r {input.fasta} -q {fqlist} -p {wdir}/{genome}_merged -t {threads} --minsupport {config[minsupport]}
        """


rule final_filtering:
    """
    Do a last filtering step on the merged genotyped dataset
    """
    input:
        merged = "{wdir}/{genome}_merged_genotype.vcf"
    output:
        filtered = "{wdir}/{genome}_filtered.vcf"
    threads: workflow.cores
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/logs/{genome}_svjedigraph.log"
    shell:
        """
        vcftools --vcf {input.merged} --remove-filtered-all --recode --stdout > {output.filtered}
        """
