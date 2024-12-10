rule svim:
    """
    SV calling with SVIM
    """
    input:
        bam = "{wdir}/{genome}_sorted.bam",
        bai = "{wdir}/{genome}_sorted.bam.bai",
        fasta = "{wdir}/{genome}.fna",
        stats = "{wdir}/mapping/{genome}_mapping.stats",
        plot = "{wdir}/mapping/{genome}_mapping_plot.html",
        sampleids = "{wdir}/{genome}.samples"
    output:
        svimvariants = temp("{wdir}/{genome}_svim/variants.vcf"),
        vcf = temp("{wdir}/{genome}_svim_tmp.vcf"),
        vcf_renamed = "{wdir}/{genome}_svim.vcf"
    resources:
        tmpdir = get_big_temp
    conda:
        "../envs/svim.yaml"
    log:
        "{wdir}/logs/{genome}_svim.log"
    shell:
        """
        svim alignment {wdir}/{genome}_svim {input.bam} {input.fasta} --insertion_sequences --read_names \
        --min_sv_size {config[min_sv_size]} --max_sv_size {config[max_sv_size]} \
        --minimum_depth {config[minimum_depth]} --min_mapq {config[min_mapq]} \
        --segment_gap_tolerance {config[segment_gap_tolerance]} --segment_overlap_tolerance {config[segment_overlap_tolerance]}
        # SVIM does not filter SV itself and outputs all variants
        bcftools filter -e "QUAL < {config[svim_quality]} || MIN(DP) < {config[svim_min_read_support]}" -o {output.vcf} -O v {wdir}/{genome}_svim/variants.vcf
        # Correct the GT field for duplications (change DUP:INT ou DUP:TANDEM to DUP)
        sed -i 's/DUP:INT/DUP/g' {output.vcf}
        sed -i 's/DUP:TANDEM/DUP/g' {output.vcf}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf}
        """


rule sniffles:
    """
    SV calling with Sniffles
    """
    input:
        bam = "{wdir}/{genome}_sorted.bam",
        bai = "{wdir}/{genome}_sorted.bam.bai",
        fasta = "{wdir}/{genome}.fna",
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf = temp("{wdir}/{genome}_sniffles_tmp.vcf"),
        vcf_renamed = "{wdir}/{genome}_sniffles.vcf"
    resources:
        tmpdir = get_big_temp
    conda:
        "../envs/sniffles.yaml"
    log:
        "{wdir}/logs/{genome}_sniffles.log"
    shell:
        """
        sniffles --input {input.bam} --vcf {output.vcf} --reference {input.fasta} \
        --threads {resources.cpus_per_task} --allow-overwrite \
        --minsvlen {config[min_sv_size]} --minsupport {config[minsupport]} \
        --minsvlen-screen-ratio {config[minsvlen-screen-ratio]} --mapq {config[mapq]} \
        --cluster-binsize {config[cluster-binsize]} --qc-coverage {config[min_coverage]} \
        --output-rnames
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf}
        """


rule cutesv:
    """
    SV calling with CuteSV
    """
    input:
        bam = "{wdir}/{genome}_sorted.bam",
        bai = "{wdir}/{genome}_sorted.bam.bai",
        fasta = "{wdir}/{genome}.fna",
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf = temp("{wdir}/{genome}_cutesv_tmp.vcf"),
        vcf_renamed = "{wdir}/{genome}_cutesv.vcf"
    resources:
        tmpdir = get_big_temp
    conda:
        "../envs/cutesv.yaml"
    log:
        "{wdir}/logs/{genome}_cutesv.log"
    shell:
        """
        mkdir -p {wdir}/cutesv
        cuteSV --max_cluster_bias_INS {config[max_cluster_bias_INS]} --diff_ratio_merging_INS {config[diff_ratio_merging_INS]} \
        --max_cluster_bias_DEL {config[max_cluster_bias_DEL]} --genotype --report_readid --diff_ratio_merging_DEL {config[diff_ratio_merging_DEL]} \
        --max_size {config[max_size]} --min_support {config[min_coverage]} --min_size {config[min_sv_size]} --min_siglength {config[min_siglength]} \
        {input.bam} {input.fasta} {output.vcf} {wdir}/cutesv/
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf}
        """


rule debreak:
    """
    SV calling with DeBreak
    """
    input:
        bam = "{wdir}/{genome}_sorted.bam",
        bai = "{wdir}/{genome}_sorted.bam.bai",
        fasta = "{wdir}/{genome}.fna",
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf = temp("{wdir}/{genome}_debreak_tmp.vcf"),
        vcf_renamed = "{wdir}/{genome}_debreak.vcf"
    conda:
        "../envs/debreak.yaml"
    # resources:
    #     tmpdir = get_big_temp
    log:
        "{wdir}/logs/{genome}_debreak.log"
    shell:
        """
        debreak --bam {input.bam} --outpath {wdir}/debreak/ --rescue_large_ins --rescue_dup -t {resources.cpus_per_task} \
        --min_size {config[min_sv_size]} --min_support {config[min_coverage]} --poa --ref {input.fasta}
        mv {wdir}/debreak/debreak.vcf {output.vcf}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf}
        """



rule removeBND:
    """
    Remove BND before merging - BND are difficult to treat in downstream analyses
    Remove TRANSLOCATION (TRA)
    """
    input:
        sniffles = "{wdir}/{genome}_sniffles.vcf",
        svim = "{wdir}/{genome}_svim.vcf",
        cutesv = "{wdir}/{genome}_cutesv.vcf",
        debreak = "{wdir}/{genome}_debreak.vcf"
    output:
        svim = "{wdir}/{genome}_svim_noBND.vcf",
        cutesv = "{wdir}/{genome}_cutesv_noBND.vcf",
        debreak = "{wdir}/{genome}_debreak_noBND.vcf",
        sniffles = "{wdir}/{genome}_sniffles_noBND.vcf"
    shell:
        """
        cat {input.svim} | grep -v '[a-zA-Z]*.BND' > {output.svim}
        cat {input.cutesv} | grep -v '[a-zA-Z]*.BND' > {output.cutesv}
        cat {input.debreak} | grep -v 'SVTYPE=BND' | grep -v 'SVTPE=TRA' > {output.debreak}
        cat {input.sniffles} | grep -v '[a-zA-Z]*.BND' > {output.sniffles}
        """


rule sniffles2plot:
    """
    Run sniffles2-plot for each SV caller
    The sniffles2-lot package output a set of QC summary plots for a single VCF
    """
    input:
        svim = "{wdir}/{genome}_svim_noBND.vcf",
        cutesv = "{wdir}/{genome}_cutesv_noBND.vcf",
        debreak = "{wdir}/{genome}_debreak_noBND.vcf",
        sniffles = "{wdir}/{genome}_sniffles_noBND.vcf"
    output:
        "{wdir}/{genome}_sniffles_QC/variant_count.jpg",
        "{wdir}/{genome}_svim_QC/variant_count.jpg",
        "{wdir}/{genome}_cutesv_QC/variant_count.jpg"
    threads: workflow.cores
    conda:
        "../envs/sniffles.yaml"
    log:
        "{wdir}/logs/{genome}_sniffles2plot.log"
    shell:
        """
        python3 -m sniffles2_plot -i {input.sniffles} -o {wdir}/{genome}_sniffles_QC/
        python3 -m sniffles2_plot -i {input.svim} -o {wdir}/{genome}_svim_QC/
        python3 -m sniffles2_plot -i {input.cutesv} -o {wdir}/{genome}_cutesv_QC/
        """


rule genotype_svim:
    """
    Use SVjedigraph to genotype one callset
    used downstream to estimate uncertainty with ensemble methods
    """
    input:
        vcf = "{wdir}/{genome}_svim_noBND.vcf",
        fasta = "{wdir}/{genome}.fna",
        fastq = expand("{wdir}/fastq/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"]),
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf = temp("{wdir}/svim_genotype/{genome}_svim_genotype_tmp.vcf"),
        vcf_renamed = "{wdir}/svim_genotype/{genome}_svim_genotype.vcf",
        gfa = "{wdir}/svim_genotype/{genome}_svim.gfa",
        gaf = "{wdir}/svim_genotype/{genome}_svim.gaf",
        aln = "{wdir}/svim_genotype/{genome}_svim_informative_aln.json"
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_genotype_svim.log"
    shell:
        """
        svjedi-graph.py -v {input.vcf} -r {input.fasta} \
        -q {fqlist} -p {wdir}/svim_genotype/{genome}_svim \
        -t {resources.cpus_per_task} \
        --minsupport 1
        mv {output.vcf_renamed} {output.vcf}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf}
        """

rule genotype_cutesv:
    """
    Use SVjedigraph to genotype one callset
    used downstream to estimate uncertainty with ensemble methods
    """
    input:
        vcf = "{wdir}/{genome}_cutesv_noBND.vcf",
        fasta = "{wdir}/{genome}.fna",
        fastq = expand("{wdir}/fastq/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"]),
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf = temp("{wdir}/cutesv_genotype/{genome}_cutesv_genotype_tmp.vcf"),
        vcf_renamed = "{wdir}/cutesv_genotype/{genome}_cutesv_genotype.vcf",
        gfa = "{wdir}/cutesv_genotype/{genome}_cutesv.gfa",
        gaf = "{wdir}/cutesv_genotype/{genome}_cutesv.gaf",
        aln = "{wdir}/cutesv_genotype/{genome}_cutesv_informative_aln.json"
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_genotype_cutesv.log"
    shell:
        """
        svjedi-graph.py -v {input.vcf} -r {input.fasta} \
        -q {fqlist} -p {wdir}/cutesv_genotype/{genome}_cutesv \
        -t {resources.cpus_per_task} \
        --minsupport 1
        mv {output.vcf_renamed} {output.vcf}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf}
        """

rule genotype_sniffles:
    """
    Use SVjedigraph to genotype one callset
    used downstream to estimate uncertainty with ensemble methods
    """
    input:
        vcf = "{wdir}/{genome}_sniffles_noBND.vcf",
        fasta = "{wdir}/{genome}.fna",
        fastq = expand("{wdir}/fastq/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"]),
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf = temp("{wdir}/sniffles_genotype/{genome}_sniffles_genotype_tmp.vcf"),
        vcf_renamed = "{wdir}/sniffles_genotype/{genome}_sniffles_genotype.vcf",
        gfa = "{wdir}/sniffles_genotype/{genome}_sniffles.gfa",
        gaf = "{wdir}/sniffles_genotype/{genome}_sniffles.gaf",
        aln = "{wdir}/sniffles_genotype/{genome}_sniffles_informative_aln.json"
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_genotype_sniffles.log"
    shell:
        """
        svjedi-graph.py -v {input.vcf} -r {input.fasta} \
        -q {fqlist} -p {wdir}/sniffles_genotype/{genome}_sniffles \
        -t {resources.cpus_per_task} \
        --minsupport 1
        mv {output.vcf_renamed} {output.vcf}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf}
        """


rule normalize_vcf_debreak:
    """
    Add INFO tags in Debreak output
    """
    input:
        vcf = "{wdir}/{genome}_debreak_noBND.vcf",
        fasta = "{wdir}/{genome}.fna",
        bam = "{wdir}/{genome}_sorted.bam"
    output:
        vcf = "{wdir}/{genome}_debreak_normalize.vcf",
        vcflist = "{wdir}/{genome}_vcf_list_debreak.txt",
        bamlist = "{wdir}/{genome}_bam_list_debreak.txt"
    threads: workflow.cores
    conda:
        "../envs/jasminesv.yaml"
    log:
        "{wdir}/logs/{genome}_normalize_debreak.log"
    shell:
        """
        echo "{input.vcf}" > {output.vcflist}
        echo "{input.bam}" > {output.bamlist}
        jasmine file_list={output.vcflist} out_file={output.vcf} \
        genome_file={input.fasta} \
        out_dir={wdir}/debreak_normalize bam_list={output.bamlist} \
        --preprocess-only --dup_to_ins
        """


rule genotype_debreak:
    """
    Use SVjedigraph to genotype one callset
    used downstream to estimate uncertainty with ensemble methods
    """
    input:
        vcf = "{wdir}/{genome}_debreak_normalize.vcf",
        fasta = "{wdir}/{genome}.fna",
        fastq = expand("{wdir}/fastq/{sample}.fastq.gz", wdir=wdir, sample=samples["sra"]),
        sampleids = "{wdir}/{genome}.samples"        
    output:
        vcf = temp("{wdir}/debreak_genotype/{genome}_debreak_genotype_tmp.vcf"),
        vcf_renamed = "{wdir}/debreak_genotype/{genome}_debreak_genotype.vcf",
        gfa = "{wdir}/debreak_genotype/{genome}_debreak.gfa",
        gaf = "{wdir}/debreak_genotype/{genome}_debreak.gaf",
        aln = "{wdir}/debreak_genotype/{genome}_debreak_informative_aln.json"
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_genotype_debreak.log"
    shell:
        """
        svjedi-graph.py -v {input.vcf} -r {input.fasta} \
        -q {input.fastq} -p {wdir}/debreak_genotype/{genome}_debreak \
        -t {resources.cpus_per_task} \
        --minsupport 1
        mv {output.vcf_renamed} {output.vcf}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf}
        """
