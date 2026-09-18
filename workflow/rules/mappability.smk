rule mosdepth_summary:
    input:
        bam="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam.bai",
    output:
        dist="{wdir}/{sample}/callability/{genome}_{aligner}.mosdepth.global.dist.txt",
        summary="{wdir}/{sample}/callability/{genome}_{aligner}.mosdepth.summary.txt",
        coverage_windows="{wdir}/{sample}/callability/{genome}_{aligner}.regions.bed.gz",
    conda:
        "../envs/mosdepth.yaml"
    log:
        "{wdir}/{sample}/logs/mosdepth/{genome}_{aligner}.txt",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_{aligner}.mosdepth_summary.tsv"
    params:
        prefix="{wdir}/{sample}/callability/{genome}_{aligner}",
    shell:
        """
        mosdepth --no-per-base \
            -t {threads} \
            --by {config[mosdepth_windows_size]} \
            {params.prefix} \
            {input.bam}
        """


rule mosdepth_quantize:
    input:
        summary="{wdir}/{sample}/callability/{genome}_{aligner}.mosdepth.summary.txt",
        bam="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam.bai",
    output:
        quantized="{wdir}/{sample}/callability/{genome}_{aligner}.quantized.bed.gz",
        quantized_idx="{wdir}/{sample}/callability/{genome}_{aligner}.quantized.bed.gz.csi",
    conda:
        "../envs/mosdepth.yaml"
    log:
        "{wdir}/{sample}/logs/mosdepth_quantize/{genome}_{aligner}.txt",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_{aligner}.mosdepth_quantize.tsv"
    params:
        # prefix = "{wdir}/{sample}/callability/{genome}_{aligner}",
        lower=round(config["quantize_cov_threshold_lower"]),
        upper=round(config["quantize_cov_threshold_upper"]),
        sample_mean=lambda wildcards, input: get_mean_cov(input.summary),
        upper_threshold=lambda wildcards, input: round(
            config["quantize_cov_threshold_upper"] * get_mean_cov(input.summary)
        ),
    shell:
        """
        export MOSDEPTH_Q0=NO_COVERAGE
        export MOSDEPTH_Q1=LOW_COVERAGE
        export MOSDEPTH_Q2=CALLABLE
        export MOSDEPTH_Q3=HIGH_COVERAGE
        
        mosdepth --no-per-base -t {threads} \
        --quantize 0:1:{params.lower}:{params.upper_threshold}: \
        {wdir}/{wildcards.sample}/callability/{genome}_{wildcards.aligner} {input.bam}
        """


rule callable_bed:
    input:
        quantized="{wdir}/{sample}/callability/{genome}_{aligner}.quantized.bed.gz",
        quantized_idx="{wdir}/{sample}/callability/{genome}_{aligner}.quantized.bed.gz.csi",
    output:
        callable_bed="{wdir}/{sample}/callability/{genome}_{aligner}_callable.bed",
    conda:
        "../envs/mosdepth.yaml"
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_{aligner}.callable_bed.tsv"
    shell:
        """
        zcat {input.quantized} | grep CALLABLE | bedtools sort | bedtools merge > {output.callable_bed}
        # Check if file has any lines
        if [ $(wc -l < {output.callable_bed}) -eq 0 ]; then
            echo "File {output.callable_bed} is empty"
            exit 1
        fi
        """


rule genmap:
    input:
        ref="{wdir}/genome/{genome}.fna",
    output:
        bg=temp("{wdir}/genmap/{genome}.genmap.bedgraph"),
        sorted_bg="{wdir}/genmap/{genome}_sorted_mappability.bg",
    params:
        indir="{wdir}/genmap_index",
        outdir="{wdir}/genmap",
        kmer=config["mappability_k"],
    log:
        "{wdir}/logs/genmap/{genome}.txt",
    benchmark:
        "{wdir}/genome/benchmarks/{genome}.genmap.tsv"
    conda:
        "../envs/genmap.yaml"
    shell:
        # snakemake creates the output directory before the shell command, but genmap doesnt like this. so we remove the directory first.
        """
        rm -rf {params.indir} && genmap index -F {input.ref} -I {params.indir} &> {log}
        genmap map -K {params.kmer} -E 0 -I {params.indir} -O {params.outdir} -bg -T {threads} -v &> {log}
        sort -k1,1 -k2,2n {output.bg} > {output.sorted_bg} 2>> {log}
        """


rule mappability_bed:
    input:
        mappable="{wdir}/genmap/{genome}_sorted_mappability.bg",
    output:
        callable_sites="{wdir}/mappability/{genome}_mappable.bed",
        tmp_map=temp("{wdir}/mappability/{genome}_temp_map.bed"),
    conda:
        "../envs/genmap.yaml"
    params:
        merge=config["mappability_merge"],
        mappability=config["mappability_min"],
    log:
        "{wdir}/logs/mappability_bed/{genome}.txt",
    benchmark:
        "{wdir}/genome/benchmarks/{genome}.mappability_bed.tsv"
    shell:
        """
        awk 'BEGIN{{OFS="\\t";FS="\\t"}} {{ if($4>={params.mappability}) print $1,$2,$3 }}' {input.mappable} > {output.tmp_map} 2> {log}
        bedtools sort -i {output.tmp_map} | bedtools merge -d {params.merge} -i - > {output.callable_sites} 2>> {log}
        """


rule add_mappability:
    input:
        callable_bed_minimap2="{wdir}/{sample}/callability/{genome}_minimap2_callable.bed",
        callable_bed_ngmlr="{wdir}/{sample}/callability/{genome}_ngmlr_callable.bed",
        mappable_bed="{wdir}/mappability/{genome}_mappable.bed",
    output:
        callable="{wdir}/{sample}/callability/{genome}_callable.bed",
        callable_mappable="{wdir}/{sample}/callability/{genome}_callable_mappable.bed",
        callable_mappable_minimap2="{wdir}/{sample}/callability/{genome}_minimap2_callable_mappable.bed",
        callable_mappable_ngmlr="{wdir}/{sample}/callability/{genome}_ngmlr_callable_mappable.bed",
    conda:
        "../envs/mosdepth.yaml"
    log:
        "{wdir}/{sample}/logs/add_mappability/{genome}.txt",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.add_mappability.tsv"
    shell:
        """
        bedtools intersect -a {input.callable_bed_minimap2} -b {input.callable_bed_ngmlr} | bedtools sort | bedtools merge > {output.callable} 2> {log}
        bedtools intersect -a {output.callable} -b {input.mappable_bed} | bedtools sort | bedtools merge > {output.callable_mappable} 2>> {log}

        bedtools intersect -a {input.callable_bed_minimap2} -b {input.mappable_bed} | bedtools sort | bedtools merge > {output.callable_mappable_minimap2} 2>> {log}

        bedtools intersect -a {input.callable_bed_ngmlr} -b {input.mappable_bed} | bedtools sort | bedtools merge > {output.callable_mappable_ngmlr} 2>> {log}
        """
