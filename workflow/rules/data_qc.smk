rule fastqc:
    """
    Report data quality for long reads
    """
    input:
        "{wdir}/{sra}.fastq.gz"
    output:
        "{wdir}/fastqc/{sra}_fastqc.html",
        "{wdir}/fastqc/{sra}_fastqc.zip"
    conda:
        "workflow/envs/fastqc.yaml"
    shell:
        """
        mkdir -p {wdir}/fastqc
        fastqc --outdir {wdir}/fastqc/ *.fastq.gz
        """


rule nanoplot:
	"""
	Quality control of raw data
	"""
	input:
		fastq = "{wdir}/{sra}.fastq.gz",
        "{wdir}/fastqc/{sra}_fastqc.html",
        "{wdir}/fastqc/{sra}_fastqc.zip"
	output:
		"{wdir}/nanoplot/{sra}_NanoStats.txt"
	conda:
		"workflow/envs/nanoplot.yaml"
	shell:
		"""
        NanoPlot --fastq {input.fastq} -t {workflow.threads} --tsv_stats --outdir {wdir}/nanoplot/ --prefix {sra} --N50 --verbose --title {sra}
        """
