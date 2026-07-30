import sys
import os
import tempfile
import random
import string

from pathlib import Path

from snakemake.exceptions import WorkflowError

# Empty cells are possible when starting from pre-aligned BAM files, where no SRA
# accession is needed.
samplelist = [wdir + "/fastq/" + s + "_sra.fastq.gz" for s in samples["sra"] if s]
fqlist = ",".join(samplelist)


def check_readable_file(path, description):
    """Fail at DAG construction time rather than deep into an expensive run."""
    if not path:
        raise WorkflowError(
            "start_from_bam is set but no {} was declared in {}.".format(
                description, config["samples"]
            )
        )
    if not Path(path).is_file():
        raise WorkflowError(
            "The {} '{}' declared in {} does not exist.".format(
                description, path, config["samples"]
            )
        )
    if Path(path).stat().st_size == 0:
        raise WorkflowError(
            "The {} '{}' declared in {} is empty.".format(
                description, path, config["samples"]
            )
        )
    return path


def resolve_bam_index(bam):
    """Locate the index of a user-supplied BAM.

    Downstream callers require a BAI index (a CSI index is not accepted), and
    samtools writes it either as <bam>.bai or, with `samtools index -o`, as
    <bam without extension>.bai. Both layouts are accepted here.
    """
    candidates = [bam + ".bai", str(Path(bam).with_suffix(".bai"))]
    for candidate in candidates:
        if Path(candidate).is_file():
            return candidate
    raise WorkflowError(
        "No BAI index found for '{}'. Looked for {}. "
        "Index the BAM with `samtools index` before running in start_from_bam mode.".format(
            bam, " and ".join("'{}'".format(c) for c in candidates)
        )
    )


if bam_mode:
    # One pre-aligned BAM per aligner: the eight callsets of the ensemble stay
    # independent. Reusing a single alignment for both aligners would make Jasmine
    # count the same evidence twice (see rules/merging.smk).
    input_bams = {
        aligner: check_readable_file(
            samples["bam_" + aligner].iloc[0], "bam_" + aligner + " file"
        )
        for aligner in aligners
    }
    input_bais = {
        aligner: resolve_bam_index(bam) for aligner, bam in input_bams.items()
    }
    # Reads are still needed: SVJedi-graph genotypes by mapping reads onto a
    # variation graph, which a linear BAM cannot substitute for.
    input_fastqs = [
        check_readable_file(f, "fastq file") for f in samples["fastq"] if f
    ]
    if not input_fastqs:
        raise WorkflowError(
            "start_from_bam is set but no 'fastq' file was declared in {}. "
            "Reads are required to genotype the SVs with SVJedi-graph.".format(
                config["samples"]
            )
        )


def get_mean_cov(summary_file):

    if not Path(summary_file).exists():
        return -1

    with open(summary_file, "r") as f:
        for line in f:
            if line.startswith("total"):
                sample_mean = float(line.split("\t")[3])

    return sample_mean


def get_big_temp(wildcards):
    """Sets a temp dir for rules that need more temp space that is typical on some cluster environments. Defaults to system temp dir."""
    if config["bigtmp"]:
        if config["bigtmp"].endswith("/"):
            return (
                config["bigtmp"]
                + "".join(random.choices(string.ascii_uppercase, k=12))
                + "/"
            )
        else:
            return (
                config["bigtmp"]
                + "/"
                + "".join(random.choices(string.ascii_uppercase, k=12))
                + "/"
            )
    else:
        return tempfile.gettempdir()