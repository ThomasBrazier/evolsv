rule svim:
    """
    SV calling with SVIM
    """
    input:
        bam = "{wdir}/{sra}_sorted.bam",
        bai = "{wdir}/{sra}_sorted.bam.bai",
        fasta = "{wdir}/{sra}.fna",
        stats = "{wdir}/mapping/{sra}_mapping.stats",
        plot = "{wdir}/mapping/{sra}_mapping_plot"
    output:
        vcf = "{wdir}/{sra}_svim.vcf",
        nobnd = "$sra/${sra}_svim_noBND.vcf",
        temporary("${sra}/variants.vcf")
    conda:
        "workflow/envs/svim.yaml"
    shell:
        """
        svim alignment {sra} {input.bam} {input.fasta} --insertion_sequences --read_names --min_sv_size {config[min_sv_size]} --minimum_depth {config[min_coverage]}
        # SVIM does not filter SV itself and outputs all variants
        bcftools view -i "QUAL >= {config[svim_quality]}" ${sra}/variants.vcf > {output.vcf}
        cat {output.vcf} | grep -v svim.BND > {output.nobnd}
        """


rule sniffles:
    """
    SV calling with Sniffles
    """
    input:
        svim = "{wdir}/{sra}_svim.vcf",
        bam = "{wdir}/{sra}_sorted.bam",
        bai = "{wdir}/{sra}_sorted.bam.bai",
        fasta = "{wdir}/{sra}.fna"
    output:
        "{wdir}/{sra}_sniffles.vcf"
    conda:
        "workflow/envs/sniffles.yaml"
    shell:
        """
        sniffles --input {input.bam} --vcf {output} --reference {input.fasta} --threads {workflow.threads} --allow-overwrite --minsvlen {config[min_sv_size]} --qc-coverage {config[min_coverage]} --output-rnames
        """


rule cutesv:
    """
    SV calling with CuteSV
    """
    input:
        sniffles = "{wdir}/{sra}_sniffles.vcf",
        svim = "{wdir}/{sra}_svim.vcf",
        bam = "{wdir}/{sra}_sorted.bam",
        bai = "{wdir}/{sra}_sorted.bam.bai",
        fasta = "{wdir}/{sra}.fna"
    output:
        "{wdir}/{sra}_cutesv.vcf"
    conda:
        "workflow/envs/cutesv.yaml"
    shell:
        """
        mkdir -p {wdir}/cutesv
        rm -rf sra/cutesv/*
        cuteSV --max_cluster_bias_INS {config[max_cluster_bias_INS]} --diff_ratio_merging_INS {config[diff_ratio_merging_INS]} \
         --max_cluster_bias_DEL {config[max_cluster_bias_DEL]} --genotype --report_readid --diff_ratio_merging_DEL {config[diff_ratio_merging_DEL]} \
         --max_size [config[max_size]] --min_support {config[min_coverage]} --min_size {config[min_sv_size]} {input.bam} {input.fasta} {output} {sra}/cutesv/
        """