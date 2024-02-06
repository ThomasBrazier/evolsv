
rule vcf2tsv:
    """
    Convert the final SV callset vcf file into a tsv file
    """
    input:
        final = "{wdir}/{sra}_merged_genotype.vcf",
        sniffles = "{wdir}/{sra}_sniffles.vcf",
        svim = "{wdir}/{sra}_svim.vcf",
        cutesv = "{wdir}/{sra}_cutesv.vcf",
        debreak = "{wdir}/{sra}_debreak.vcf"
    output:
        final = "{wdir}/{sra}_merged_genotype.tsv",
        sniffles = "{wdir}/{sra}_sniffles.tsv",
        svim = "{wdir}/{sra}_svim.tsv",
        cutesv = "{wdir}/{sra}_cutesv.tsv",
        debreak = "{wdir}/{sra}_debreak.tsv"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # CuteSV
        echo -e "CHROM\tPOS\tEND\tTYPE\tLENGTH\tGENOTYPE\tGENOTYPE_QUALITY\tHQ_REFERENCE_READS\tHQ_VARIANTS_READS\tREADS_SUPPORT\tALLELE_FREQUENCY" > {output.cutesv}
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\t%INFO/RE\t%INFO/AF\n' {input.cutesv} >> {output.cutesv}

        # SVIM
        echo -e "CHROM\tPOS\tEND\tTYPE\tLENGTH\tGENOTYPE\tREADS_SUPPORT\tREADS_DEPTH\tRD_FOR_EACH_ALLELE" > {output.svim}
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t[%GT]\t%INFO/SUPPORT\t[%DP]\t[%AD]\n' {input.svim} >> {output.svim}

        # Sniffles
        echo -e "CHROM\tPOS\tEND\tTYPE\tLENGTH\tGENOTYPE\tGENOTYPE_QUALITY\tHQ_REFERENCE_READS\tHQ_VARIANTS_READS\tREADS_SUPPORT\tALLELE_FREQUENCY" > {output.sniffles}
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\t%INFO/SUPPORT\t%INFO/AF\n' {input.sniffles} >> {output.sniflles}

        # DeBreak
        echo -e "CHROM\tPOS\tEND\tTYPE\tLENGTH\tGENOTYPE\tGENOTYPE_QUALITY\tHQ_REFERENCE_READS\tHQ_VARIANTS_READS\tREADS_SUPPORT\tALLELE_FREQUENCY" > {output.debreak}
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\t%INFO/SUPPORT\t%INFO/AF\n' {input.debreak} >> {output.debreak}

        # Final callset
        echo -e "CHROM\tPOS\tEND\tTYPE\tLENGTH\tGENOTYPE\tHQ_REFERENCE_READS\tHQ_VARIANTS_READS\tREADS_SUPPORT\tTOOLS" > {output.final}
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t[ %GT]\t[ %DR]\t[ %DV]\t%INFO/SUPPORT\t%INFO/IDLIST\n' {input.final} >> {output.final}

        """



rule finalreport:
    """
    Compute and print a summary report for assembly, mapping, SV calling, merging and genotyping
    """
    input:
        vcf = "{wdir}/{sra}_merged_genotype.vcf",
        final = "{wdir}/{sra}_merged_genotype.tsv",
        sniffles = "{wdir}/{sra}_sniffles.tsv",
        svim = "{wdir}/{sra}_svim.tsv",
        cutesv = "{wdir}/{sra}_cutesv.tsv",
        debreak = "{wdir}/{sra}_debreak.tsv"
    output:
        "{wdir}/{sra}_finalQC.html"
    conda:
        "../envs/Renv.yaml"
    shell:
        """
        Rscript workflow/scripts/finalQC.R {wdir} {sra}
        """