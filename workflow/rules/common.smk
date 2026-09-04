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


# Sequencing-technology presets. Each entry is an ordinary config key: the preset only
# supplies a default (config.setdefault below), so a value set explicitly in the config
# file still wins.
#
# cuteSV clustering values are the author's per-technology recommendations,
# https://github.com/tjiangHIT/cuteSV#recommendation-parameters
# Sniffles2, SVIM, DeBreak and SVJedi-graph publish no technology presets, so nothing
# else in the workflow is technology-dependent. chopper is deliberately absent too:
# chopper_quality: 10 suits HiFi and modern ONT alike (see README).
TECH_PRESETS = {
    "hifi": {
        "minimap_ax": "map-hifi",
        "ngmlr_preset": "pacbio",
        "read_group_platform": "PACBIO",  # a SAM specification @RG PL value
        "max_cluster_bias_INS": 1000,
        "diff_ratio_merging_INS": 0.9,
        "max_cluster_bias_DEL": 1000,
        "diff_ratio_merging_DEL": 0.5,
    },
    "ont": {
        "minimap_ax": "map-ont",
        "ngmlr_preset": "ont",
        "read_group_platform": "ONT",
        "max_cluster_bias_INS": 100,
        "diff_ratio_merging_INS": 0.3,
        "max_cluster_bias_DEL": 100,
        "diff_ratio_merging_DEL": 0.3,
    },
}

sequencing_technology = str(config.get("sequencing_technology", "hifi")).strip().lower()
if sequencing_technology not in TECH_PRESETS:
    raise WorkflowError(
        "Unknown sequencing_technology '{}'. Supported values are: {}.".format(
            config.get("sequencing_technology"), ", ".join(sorted(TECH_PRESETS))
        )
    )
for preset_key, preset_value in TECH_PRESETS[sequencing_technology].items():
    config.setdefault(preset_key, preset_value)


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