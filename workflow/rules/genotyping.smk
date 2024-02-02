rule svjedi_jasmine:
    """
    Use SVjedi-graph on the merged dataset to genotype SVs
    """
    input:
        "{wdir}/{sra}_merged.vcf"
    output:
        "{wdir}/{sra}_merged_genotype.vcf",
        "{wdir}/{sra}_merged.gfa",
        "{wdir}/{sra}_merged.gaf",
        "{wdir}/{sra}_merged_informative_aln.json"
    conda:
        "Envs/svjedi-graph.yaml"
    shell:
        """
        Scripts/svjedi.sh {sra} merged
        """