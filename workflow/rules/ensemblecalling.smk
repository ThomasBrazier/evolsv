rule svim:
    """
    SV calling with SVIM
    """
    input:
        bam = "{wdir}/bam/{genome}_{aligner}_sorted.bam",
        bai = "{wdir}/bam/{genome}_{aligner}_sorted.bam.bai",
        fasta = "{wdir}/genome/{genome}.fna",
        sampleids = "{wdir}/{genome}.samples"
    output:
        svimvariants = temp("{wdir}/{genome}_{aligner}_svim/variants.vcf"),
        vcf = temp("{wdir}/{genome}_{aligner}_svim_tmp.vcf"),
        vcf_raw = temp("{wdir}/{genome}_{aligner}_svim_raw.vcf"),
        vcf_renamed = "{wdir}/{genome}_{aligner}_svim.vcf"
    resources:
        tmpdir = get_big_temp
    conda:
        "../envs/svim.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_svim.log"
    shell:
        """
        svim alignment {wdir}/{genome}_{wildcards.aligner}_svim {input.bam} {input.fasta} \
        --insertion_sequences --read_names \
        --min_sv_size {config[min_sv_size]} \
        --max_sv_size {config[max_sv_size]} \
        --minimum_depth {config[minimum_depth]} \
        --min_mapq {config[min_mapq]} \
        --segment_gap_tolerance {config[segment_gap_tolerance]} \
        --segment_overlap_tolerance {config[segment_overlap_tolerance]}
        # SVIM does not filter SV itself and outputs all variants
        echo "Filter SVIM output"
        bcftools filter -e "QUAL < {config[svim_quality]} || MIN(DP) < {config[svim_min_read_support]}" -o {output.vcf} -O v {wdir}/{genome}_{wildcards.aligner}_svim/variants.vcf
        # Correct the GT field for duplications (change DUP:INT ou DUP:TANDEM to DUP)
        sed -i 's/DUP:INT/DUP/g' {output.vcf}
        sed -i 's/DUP:TANDEM/DUP/g' {output.vcf}
        # Consistent renaming of VCF header with sample id
        echo "Renaming VCF header with sample id"
        bcftools reheader --samples {input.sampleids} --output {output.vcf_raw} {output.vcf}
        bcftools view -f PASS --output-file {output.vcf_renamed} {output.vcf_raw}
        """

rule sniffles:
    """
    SV calling with Sniffles
    """
    input:
        bam = "{wdir}/bam/{genome}_{aligner}_sorted.bam",
        bai = "{wdir}/bam/{genome}_{aligner}_sorted.bam.bai",
        fasta = "{wdir}/genome/{genome}.fna",
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf = temp("{wdir}/{genome}_{aligner}_sniffles_tmp.vcf"),
        vcf_raw = temp("{wdir}/{genome}_{aligner}_sniffles_raw.vcf"),
        vcf_renamed = "{wdir}/{genome}_{aligner}_sniffles.vcf"
    resources:
        tmpdir = get_big_temp
    conda:
        "../envs/sniffles.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_sniffles.log"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --reference {input.fasta} \
        --threads {resources.cpus_per_task} \
        --allow-overwrite \
        --minsvlen {config[min_sv_size]} \
        --minsupport {config[minsupport]} \
        --minsvlen-screen-ratio {config[minsvlen-screen-ratio]} \
        --mapq {config[mapq]} \
        --cluster-binsize {config[cluster-binsize]} \
        --qc-coverage {config[min_coverage]} \
        --output-rnames
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_raw} {output.vcf}
        bcftools view -f PASS --output-file {output.vcf_renamed} {output.vcf_raw}
        """


rule cutesv:
    """
    SV calling with CuteSV
    """
    input:
        bam = "{wdir}/bam/{genome}_{aligner}_sorted.bam",
        bai = "{wdir}/bam/{genome}_{aligner}_sorted.bam.bai",
        fasta = "{wdir}/genome/{genome}.fna",
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf = temp("{wdir}/{genome}_{aligner}_cutesv_tmp.vcf"),
        vcf_raw = temp("{wdir}/{genome}_{aligner}_cutesv_raw.vcf"),
        vcf_renamed = "{wdir}/{genome}_{aligner}_cutesv.vcf"
    resources:
        tmpdir = get_big_temp
    conda:
        "../envs/cutesv.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_cutesv.log"
    shell:
        """
        mkdir -p {wdir}/cutesv_{wildcards.aligner}
        cuteSV --max_cluster_bias_INS {config[max_cluster_bias_INS]} \
        --diff_ratio_merging_INS {config[diff_ratio_merging_INS]} \
        --max_cluster_bias_DEL {config[max_cluster_bias_DEL]} \
        --genotype --report_readid \
        --diff_ratio_merging_DEL {config[diff_ratio_merging_DEL]} \
        --max_size {config[max_size]} \
        --min_support {config[min_coverage]} \
        --min_size {config[min_sv_size]} \
        --min_siglength {config[min_siglength]} \
        {input.bam} {input.fasta} {output.vcf} {wdir}/cutesv_{wildcards.aligner}/
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_raw} {output.vcf}
        bcftools view -f PASS --output-file {output.vcf_renamed} {output.vcf_raw}
        """


rule debreak:
    """
    SV calling with DeBreak
    """
    input:
        bam = "{wdir}/bam/{genome}_{aligner}_sorted.bam",
        bai = "{wdir}/bam/{genome}_{aligner}_sorted.bam.bai",
        fasta = "{wdir}/genome/{genome}.fna",
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf = temp("{wdir}/{genome}_{aligner}_debreak_tmp.vcf"),
        vcf_raw = temp("{wdir}/{genome}_{aligner}_debreak_raw.vcf"),
        vcf_renamed = "{wdir}/{genome}_{aligner}_debreak.vcf"
    conda:
        "../envs/debreak.yaml"
    resources:
        tmpdir = get_big_temp
    log:
        "{wdir}/logs/{genome}_{aligner}_debreak.log"
    shell:
        """
        debreak --bam {input.bam} \
        --outpath {wdir}/debreak_{wildcards.aligner}/ \
        --rescue_large_ins \
        --rescue_dup \
        -t {resources.cpus_per_task} \
        --min_size {config[min_sv_size]} \
        --min_support {config[min_coverage]} --poa \
        --ref {input.fasta}
        mv {wdir}/debreak_{wildcards.aligner}/debreak.vcf {output.vcf}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_raw} {output.vcf}
        bcftools view -f PASS --output-file {output.vcf_renamed} {output.vcf_raw}
        """


rule removeBND:
    """
    Remove BND before merging - BND are difficult to treat in downstream analyses
    Remove TRANSLOCATION (TRA)
    """
    input:
        sniffles_minimap2 = "{wdir}/{genome}_minimap2_sniffles.vcf",
        svim_minimap2 = "{wdir}/{genome}_minimap2_svim.vcf",
        cutesv_minimap2 = "{wdir}/{genome}_minimap2_cutesv.vcf",
        debreak_minimap2 = "{wdir}/{genome}_minimap2_debreak.vcf",
        sniffles_ngmlr = "{wdir}/{genome}_ngmlr_sniffles.vcf",
        svim_ngmlr = "{wdir}/{genome}_ngmlr_svim.vcf",
        cutesv_ngmlr = "{wdir}/{genome}_ngmlr_cutesv.vcf",
        debreak_ngmlr = "{wdir}/{genome}_ngmlr_debreak.vcf"
    output:
        svim_minimap2 = temp("{wdir}/{genome}_minimap2_svim_noBND.vcf"),
        cutesv_minimap2 = temp("{wdir}/{genome}_minimap2_cutesv_noBND.vcf"),
        debreak_minimap2 = temp("{wdir}/{genome}_minimap2_debreak_noBND.vcf"),
        sniffles_minimap2 = temp("{wdir}/{genome}_minimap2_sniffles_noBND.vcf"),
        svim_ngmlr = temp("{wdir}/{genome}_ngmlr_svim_noBND.vcf"),
        cutesv_ngmlr = temp("{wdir}/{genome}_ngmlr_cutesv_noBND.vcf"),
        debreak_ngmlr = temp("{wdir}/{genome}_ngmlr_debreak_noBND.vcf"),
        sniffles_ngmlr = temp("{wdir}/{genome}_ngmlr_sniffles_noBND.vcf")
    shell:
        """
        cat {input.svim_minimap2} | grep -v '[a-zA-Z]*.BND' > {output.svim_minimap2}
        cat {input.cutesv_minimap2} | grep -v '[a-zA-Z]*.BND' > {output.cutesv_minimap2}
        cat {input.debreak_minimap2} | grep -v 'SVTYPE=BND' | grep -v 'SVTYPE=TRA' > {output.debreak_minimap2}
        cat {input.sniffles_minimap2} | grep -v '[a-zA-Z]*.BND' > {output.sniffles_minimap2}

        cat {input.svim_ngmlr} | grep -v '[a-zA-Z]*.BND' > {output.svim_ngmlr}
        cat {input.cutesv_ngmlr} | grep -v '[a-zA-Z]*.BND' > {output.cutesv_ngmlr}
        cat {input.debreak_ngmlr} | grep -v 'SVTYPE=BND' | grep -v 'SVTYPE=TRA' > {output.debreak_ngmlr}
        cat {input.sniffles_ngmlr} | grep -v '[a-zA-Z]*.BND' > {output.sniffles_ngmlr}
        """

rule vcf_preprocess:
    """
    Preprocess VCF with Jasmine to convert DUP to INS and add consistent INFO tags
    """
    input:
        vcf = "{wdir}/{genome}_{aligner}_{caller}_noBND.vcf",
        fasta = "{wdir}/genome/{genome}.fna",
        bam = "{wdir}/bam/{genome}_{aligner}_sorted.bam"
    output:
        vcf_temp = temp("{wdir}/preprocess/{genome}_{aligner}_{caller}_preprocess_temp.vcf"),
        vcflist = temp("{wdir}/{genome}_{aligner}_{caller}_vcf_list.txt"),
        bamlist = temp("{wdir}/{genome}_{aligner}_{caller}_bam_list.txt")
    threads: workflow.cores
    conda:
        "../envs/jasminesv.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_{caller}_preprocess.log"
    shell:
        """
        echo "{input.vcf}" > {output.vcflist}
        echo "{input.bam}" > {output.bamlist}
        jasmine file_list={output.vcflist} out_file={output.vcf_temp} \
        genome_file={input.fasta} \
        out_dir={wdir}/preprocess bam_list={output.bamlist} \
        --preprocess-only --dup_to_ins
        """



rule dup_to_ins:
    """
    Convert DUP to INS
    - must be defined as an insertion event whith `CHR` and `POS` corresponding to the position of insertion of the novel copy
    - `INFO` field must contain `SVTYPE=INS`
    - `ALT` field must contain the sequence of the duplication
    """
    input:
        vcf_temp = "{wdir}/preprocess/{genome}_{aligner}_{caller}_preprocess_temp.vcf",
        fasta = "{wdir}/genome/{genome}.fna"
    output:
        vcf = temp("{wdir}/preprocess/{genome}_{aligner}_{caller}_preprocess.vcf")
    threads: workflow.cores
    conda:
        "../envs/pysam.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_{caller}_preprocess.log"
    shell:
        """
        # Add DUP sequence in ALT
        # ISSUE. This code require an older pysam version (= 0.10).
        # More recent versions remove the END INFO tag, hence errors in SVjedi-graph
        python workflow/scripts/add_dup_to_ins.py {input.vcf_temp} {output.vcf} {input.fasta}
        sed -i 's/SVTYPE=DUP/SVTYPE=INS/g' {output.vcf}
        """



rule sniffles2plot:
    """
    Run sniffles2-plot for each SV caller
    The sniffles2-plot package output a set of QC summary plots for a single VCF
    """
    input:
        svim_minimap2 = "{wdir}/{genome}_minimap2_svim_noBND.vcf",
        cutesv_minimap2 = "{wdir}/{genome}_minimap2_cutesv_noBND.vcf",
        sniffles_minimap2 = "{wdir}/{genome}_minimap2_sniffles_noBND.vcf",
        svim_ngmlr = "{wdir}/{genome}_ngmlr_svim_noBND.vcf",
        cutesv_ngmlr = "{wdir}/{genome}_ngmlr_cutesv_noBND.vcf",
        sniffles_ngmlr = "{wdir}/{genome}_ngmlr_sniffles_noBND.vcf"
    output:
        "{wdir}/minimap2_sniffles_QC_{genome}/variant_count.jpg",
        "{wdir}/minimap2_svim_QC_{genome}/variant_count.jpg",
        "{wdir}/minimap2_cutesv_QC_{genome}/variant_count.jpg",
        "{wdir}/ngmlr_sniffles_QC_{genome}/variant_count.jpg",
        "{wdir}/ngmlr_svim_QC_{genome}/variant_count.jpg",
        "{wdir}/ngmlr_cutesv_QC_{genome}/variant_count.jpg"
    conda:
        "../envs/sniffles.yaml"
    log:
        "{wdir}/logs/{genome}_sniffles2plot.log"
    shell:
        """
        python3 -m sniffles2_plot -i {input.sniffles_minimap2} -o {wdir}/minimap2_sniffles_QC_{genome}/
        python3 -m sniffles2_plot -i {input.svim_minimap2} -o {wdir}/minimap2_svim_QC_{genome}/
        python3 -m sniffles2_plot -i {input.cutesv_minimap2} -o {wdir}/minimap2_cutesv_QC_{genome}/

        python3 -m sniffles2_plot -i {input.sniffles_ngmlr} -o {wdir}/ngmlr_sniffles_QC_{genome}/
        python3 -m sniffles2_plot -i {input.svim_ngmlr} -o {wdir}/ngmlr_svim_QC_{genome}/
        python3 -m sniffles2_plot -i {input.cutesv_ngmlr} -o {wdir}/ngmlr_cutesv_QC_{genome}/
        """


rule genotype_svim:
    """
    Use SVjedigraph to genotype one callset
    used downstream to estimate uncertainty with ensemble methods
    """
    input:
        vcf = "{wdir}/preprocess/{genome}_{aligner}_svim_preprocess.vcf",
        fasta = "{wdir}/genome/{genome}.fna",
        merged_fastq = "{wdir}/fastq/{genome}_filtered.fastq.gz",
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf_temp = temp("{wdir}/svim_genotype/{genome}_{aligner}_svim_genotype_tmp.vcf"),
        vcf_renamed = "{wdir}/svim_genotype/{genome}_{aligner}_svim_genotype.vcf",
        gfa = "{wdir}/svim_genotype/{genome}_{aligner}_svim.gfa",
        gaf = "{wdir}/svim_genotype/{genome}_{aligner}_svim.gaf",
        aln = "{wdir}/svim_genotype/{genome}_{aligner}_svim_informative_aln.json"
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_genotype_svim.log"
    shell:
        """
        svjedi-graph.py -v {input.vcf} -r {input.fasta} \
        -q {input.merged_fastq} -p {wdir}/svim_genotype/{genome}_{wildcards.aligner}_svim \
        -t {resources.cpus_per_task} \
        --minsupport {config[minsupport]}
        mv {output.vcf_renamed} {output.vcf_temp}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf_temp}
        """

rule genotype_cutesv:
    """
    Use SVjedigraph to genotype one callset
    used downstream to estimate uncertainty with ensemble methods
    """
    input:
        vcf = "{wdir}/preprocess/{genome}_{aligner}_cutesv_preprocess.vcf",
        fasta = "{wdir}/genome/{genome}.fna",
        merged_fastq = "{wdir}/fastq/{genome}_filtered.fastq.gz",
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf_temp = temp("{wdir}/cutesv_genotype/{genome}_{aligner}_cutesv_genotype_tmp.vcf"),
        vcf_renamed = "{wdir}/cutesv_genotype/{genome}_{aligner}_cutesv_genotype.vcf",
        gfa = "{wdir}/cutesv_genotype/{genome}_{aligner}_cutesv.gfa",
        gaf = "{wdir}/cutesv_genotype/{genome}_{aligner}_cutesv.gaf",
        aln = "{wdir}/cutesv_genotype/{genome}_{aligner}_cutesv_informative_aln.json"
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_genotype_cutesv.log"
    shell:
        """
        svjedi-graph.py -v {input.vcf} -r {input.fasta} \
        -q {input.merged_fastq} -p {wdir}/cutesv_genotype/{genome}_{wildcards.aligner}_cutesv \
        -t {resources.cpus_per_task} \
        --minsupport {config[minsupport]}
        mv {output.vcf_renamed} {output.vcf_temp}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf_temp}
        """

rule genotype_sniffles:
    """
    Use SVjedigraph to genotype one callset
    used downstream to estimate uncertainty with ensemble methods
    """
    input:
        vcf = "{wdir}/preprocess/{genome}_{aligner}_sniffles_preprocess.vcf",
        fasta = "{wdir}/genome/{genome}.fna",
        merged_fastq = "{wdir}/fastq/{genome}_filtered.fastq.gz",
        sampleids = "{wdir}/{genome}.samples"
    output:
        vcf_temp = temp("{wdir}/sniffles_genotype/{genome}_{aligner}_sniffles_genotype_tmp.vcf"),
        vcf_renamed = "{wdir}/sniffles_genotype/{genome}_{aligner}_sniffles_genotype.vcf",
        gfa = "{wdir}/sniffles_genotype/{genome}_{aligner}_sniffles.gfa",
        gaf = "{wdir}/sniffles_genotype/{genome}_{aligner}_sniffles.gaf",
        aln = "{wdir}/sniffles_genotype/{genome}_{aligner}_sniffles_informative_aln.json"
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_genotype_sniffles.log"
    shell:
        """
        svjedi-graph.py -v {input.vcf} -r {input.fasta} \
        -q {input.merged_fastq} -p {wdir}/sniffles_genotype/{genome}_{wildcards.aligner}_sniffles \
        -t {resources.cpus_per_task} \
        --minsupport {config[minsupport]}
        mv {output.vcf_renamed} {output.vcf_temp}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf_temp}
        """


rule genotype_debreak:
    """
    Use SVjedigraph to genotype one callset
    used downstream to estimate uncertainty with ensemble methods
    """
    input:
        vcf = "{wdir}/preprocess/{genome}_{aligner}_debreak_preprocess.vcf",
        fasta = "{wdir}/genome/{genome}.fna",
        merged_fastq = "{wdir}/fastq/{genome}_filtered.fastq.gz",
        sampleids = "{wdir}/{genome}.samples"     
    output:
        vcf_temp = temp("{wdir}/debreak_genotype/{genome}_{aligner}_debreak_genotype_tmp.vcf"),
        vcf_renamed = "{wdir}/debreak_genotype/{genome}_{aligner}_debreak_genotype.vcf",
        gfa = "{wdir}/debreak_genotype/{genome}_{aligner}_debreak.gfa",
        gaf = "{wdir}/debreak_genotype/{genome}_{aligner}_debreak.gaf",
        aln = "{wdir}/debreak_genotype/{genome}_{aligner}_debreak_informative_aln.json"
    conda:
        "../envs/svjedi-graph.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_genotype_debreak.log"
    shell:
        """
        svjedi-graph.py -v {input.vcf} -r {input.fasta} \
        -q {input.merged_fastq} -p {wdir}/debreak_genotype/{genome}_{wildcards.aligner}_debreak \
        -t {resources.cpus_per_task} \
        --minsupport {config[minsupport]}
        mv {output.vcf_renamed} {output.vcf_temp}
        # Consistent renaming of VCF header with sample id
        bcftools reheader --samples {input.sampleids} --output {output.vcf_renamed} {output.vcf_temp}
        """


rule basic_filter:
    """
    After genotyping, apply basic filtering on the genotypes
    Filter DEPTH (min/amx DEPTH DP and min AD)
    Filter max SV size
    """
    input:
        vcf = "{wdir}/{caller}_genotype/{genome}_{aligner}_{caller}_genotype.vcf"
    output:
        vcf = "{wdir}/filtered/{genome}_{aligner}_{caller}_filtered.vcf"
    threads: workflow.cores
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/logs/{genome}_{aligner}_{caller}_basic_filter.log"
    shell:
        """
        bcftools filter -e "SVLEN > {config[max_sv_size]} || MIN(AD) < {config[min_alt_depth]} || MIN(DP) < {config[min_depth]} || MAX(DP) > {config[max_depth]}" -o {output.vcf} -O v {input.vcf}
        """
