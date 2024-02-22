rule svim:
    """
    SV calling with SVIM
    """
    input:
        bam = "{wdir}/{sra}_sorted.bam",
        bai = "{wdir}/{sra}_sorted.bam.bai",
        fasta = "{wdir}/{sra}.fna",
        stats = "{wdir}/mapping/{sra}_mapping.stats",
        plot = "{wdir}/mapping/{sra}_mapping_plot.html"
    output:
        temporary("{wdir}/{sra}_svim/variants.vcf"),
        vcf = "{wdir}/{sra}_svim.vcf"
    threads: workflow.cores
    conda:
        "../envs/svim.yaml"
    log:
        "{wdir}/logs/{sra}_svim.log"
    shell:
        """
        svim alignment {wdir}/{sra}_svim {input.bam} {input.fasta} --insertion_sequences --read_names --min_sv_size {config[min_sv_size]} \
        --minimum_depth {config[min_coverage]} --segment_gap_tolerance {config[segment_gap_tolerance]} --segment_overlap_tolerance {config[segment_overlap_tolerance]}
        # SVIM does not filter SV itself and outputs all variants
        bcftools view -i "QUAL >= {config[svim_quality]}" {wdir}/{sra}_svim/variants.vcf > {output.vcf}
        # Correct the GT field for duplications (change DUP:INT ou DUP:TANDEM to DUP)
        sed -i 's/DUP:INT/DUP/g' {output.vcf}
        sed -i 's/DUP:TANDEM/DUP/g' {output.vcf}
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
    threads: workflow.cores
    conda:
        "../envs/sniffles.yaml"
    log:
        "{wdir}/logs/{sra}_sniffles.log"
    shell:
        """
        sniffles --input {input.bam} --vcf {output} --reference {input.fasta} --threads {threads} --allow-overwrite \
        --minsvlen {config[min_sv_size]} --qc-coverage {config[min_coverage]} --output-rnames
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
    threads: workflow.cores
    conda:
        "../envs/cutesv.yaml"
    log:
        "{wdir}/logs/{sra}_cutesv.log"
    shell:
        """
        mkdir -p {wdir}/cutesv
        cuteSV --max_cluster_bias_INS {config[max_cluster_bias_INS]} --diff_ratio_merging_INS {config[diff_ratio_merging_INS]} \
         --max_cluster_bias_DEL {config[max_cluster_bias_DEL]} --genotype --report_readid --diff_ratio_merging_DEL {config[diff_ratio_merging_DEL]} \
         --max_size {config[max_size]} --min_support {config[min_coverage]} --min_size {config[min_sv_size]} {input.bam} {input.fasta} {output} {wdir}/cutesv/
        """


rule debreak:
    """
    SV calling with DeBreak
    """
    input:
        sniffles = "{wdir}/{sra}_sniffles.vcf",
        svim = "{wdir}/{sra}_svim.vcf",
        cutesv = "{wdir}/{sra}_cutesv.vcf",
        bam = "{wdir}/{sra}_sorted.bam",
        bai = "{wdir}/{sra}_sorted.bam.bai",
        fasta = "{wdir}/{sra}.fna"
    output:
        "{wdir}/{sra}_debreak.vcf"
    threads: workflow.cores
    conda:
        "../envs/debreak.yaml"
    log:
        "{wdir}/logs/{sra}_debreak.log"
    shell:
        """
        debreak --bam {input.bam} --outpath {wdir}/debreak/ --rescue_large_ins --rescue_dup -t {threads} \
        --min_size {config[min_sv_size]} --min_support {config[min_coverage]} --poa --ref {input.fasta}
        cp {wdir}/debreak/debreak.vcf {wdir}/{sra}_debreak.vcf
        """


rule sniffles2plot:
    """
    Run sniffles2-plot for each SV caller
    The sniffles2-lot package output a set of QC summary plots for a single VCF
    """
    input:
        sniffles = "{wdir}/{sra}_sniffles.vcf",
        svim = "{wdir}/{sra}_svim.vcf",
        cutesv = "{wdir}/{sra}_cutesv.vcf",
        debreak = "{wdir}/{sra}_debreak.vcf"
    output:
        "{wdir}/{sra}_sniffles_QC/variant_count.jpg",
        "{wdir}/{sra}_svim_QC/variant_count.jpg",
        "{wdir}/{sra}_cutesv_QC/variant_count.jpg"
    threads: workflow.cores
    conda:
        "../envs/sniffles.yaml"
    log:
        "{wdir}/logs/{sra}_sniffles2plot.log"
    shell:
        """
        python3 -m sniffles2_plot -i {input.sniffles} -o {wdir}/{sra}_sniffles_QC/
        python3 -m sniffles2_plot -i {input.svim} -o {wdir}/{sra}_svim_QC/
        python3 -m sniffles2_plot -i {input.cutesv} -o {wdir}/{sra}_cutesv_QC/
        """
