rule svjedi_jasmine:
    """
    Use SVjedi-graph on the merged dataset to genotype SVs
    """
    input:
        merged = "{wdir}/{sra}_merged.vcf",
        fasta = "{wdir}/{sra}.fna",
        fastq = "{wdir}/{sra}.fastq.gz"
    output:
        "{wdir}/{sra}_merged_genotype.vcf",
        "{wdir}/{sra}_merged.gfa",
        "{wdir}/{sra}_merged.gaf",
        "{wdir}/{sra}_merged_informative_aln.json"
    threads: workflow.cores
    conda:
        "../envs/svjedi-graph.yaml"
    shell:
        """
        svjedi-graph.py -v {input.merged} -r {input.fasta} -q {input.fastq} -p {wdir}/{sra}_merged -t {threads} --minsupport {config[minsupport]}
        """
