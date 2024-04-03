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


rule filter_scaffolds:  
    input: 
        merged = "{wdir}/{genome}_merged_genotype.vcf"
    output: 
        filtered = temp("{wdir}/{genome}_merged_genotype_noscaffold.vcf")
    conda:
        "../envs/bcftools.yml"
    params:
        chr_ex = config["scaffolds_to_exclude"]
    shell:
        """
        if [ -z "{params.chr_ex}" ]
        then
            cp {input.vcf} {output.vcf}
        else
            bcftools view -t ^{params.chr_ex} \
            {input.vcf} -O u -o {output.vcf}
        fi
        """

rule final_filtering:
    """
    Do a last filtering step on the merged genotyped dataset
    """
    input:
        merged_noscaffold = "{wdir}/{genome}_merged_genotype_noscaffold.vcf"
    output:
        filtered = "{wdir}/{genome}_filtered.vcf"
    threads: workflow.cores
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/logs/{genome}_svjedigraph.log"
    shell:
        """
        vcftools --vcf {input.merged_noscaffold} --remove-filtered-all --recode --stdout > {output.filtered}
        """



