import sys
import os
import tempfile
import random
import string

from pathlib import Path

from snakemake.exceptions import WorkflowError

def rows_of(sample):
    """Sample-sheet rows of one individual, always as a DataFrame.

    samples.loc[sample] is avoided on purpose: it returns a Series for an individual
    with one row and a DataFrame for an individual with several rows.
    """
    return samples[samples["sample_name"] == sample]


def runs_of(sample):
    """SRA runs of one individual. Empty cells are possible in start_from_bam mode."""
    return [run for run in rows_of(sample)["sra"] if run]


def run_qc_files(wildcards, pattern):
    """Per-run QC files of one individual, required before its reads are aligned."""
    return expand(
        f"{{wdir}}/{{sample}}/{pattern}",
        wdir=wildcards.wdir,
        sample=wildcards.sample,
        run=runs_of(wildcards.sample),
    )


def per_individual(pattern, **wildcards):
    """Expand a {wdir}/{sample}/ target for every individual in the sample sheet."""
    return expand(pattern, wdir=wdir, genome=genome, sample=individuals, **wildcards)


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


def config_flag(key):
    """Read a boolean config key.

    `--config key=false` gives the string "false", which Python treats as true, so
    strings are parsed explicitly and any other value is rejected.
    """
    value = config[key]
    if isinstance(value, bool):
        return value
    if str(value).strip().lower() in ("true", "false"):
        return str(value).strip().lower() == "true"
    raise WorkflowError(
        "Config key '{}' must be true or false. Found: '{}'.".format(key, value)
    )


def check_readable_file(path, description, declared_in=None):
    """Fail at DAG construction time rather than deep into an expensive run.

    declared_in names where the path comes from; it defaults to the sample sheet.
    """
    declared_in = declared_in or config["samples"]
    if not path:
        raise WorkflowError(
            "start_from_bam is set but no {} was declared in {}.".format(
                description, declared_in
            )
        )
    if not Path(path).is_file():
        raise WorkflowError(
            "The {} '{}' declared in {} does not exist.".format(
                description, path, declared_in
            )
        )
    if Path(path).stat().st_size == 0:
        raise WorkflowError(
            "The {} '{}' declared in {} is empty.".format(description, path, declared_in)
        )
    return path


# Source of the reference genome and of its assembly metadata. Only three combinations
# are valid, so the metadata always describes the FASTA it is used with:
#   ncbi         nothing local: FASTA and metadata downloaded from NCBI
#   local_fasta  reference_fasta only: local FASTA, metadata downloaded from NCBI
#   local        reference_fasta + sequence_report (+ optional assembly_data_report):
#                nothing is downloaded
reference_fasta = str(config.get("reference_fasta") or "")
sequence_report = str(config.get("sequence_report") or "")
assembly_data_report = str(config.get("assembly_data_report") or "")

if sequence_report and not reference_fasta:
    raise WorkflowError(
        "sequence_report is set but reference_fasta is not. A local sequence report "
        "must describe a local FASTA: set reference_fasta too, or remove sequence_report."
    )
if assembly_data_report and not sequence_report:
    raise WorkflowError(
        "assembly_data_report is set but sequence_report is not. Local assembly "
        "metadata needs both reference_fasta and sequence_report."
    )
for key, path in [
    ("reference_fasta", reference_fasta),
    ("sequence_report", sequence_report),
    ("assembly_data_report", assembly_data_report),
]:
    if path:
        check_readable_file(path, f"{key} file", declared_in="the config file")

if sequence_report:
    reference_source = "local"
elif reference_fasta:
    reference_source = "local_fasta"
else:
    reference_source = "ncbi"


def resolve_bam_index(bam):
    """Locate the index of a user-supplied BAM.

    Downstream callers require a BAI index (a CSI index is not accepted), and
    samtools writes it either as <bam>.bai or, with `samtools index -o`, as
    <bam without extension>.bai. Both layouts are accepted here.
    """
    candidates = [f"{bam}.bai", str(Path(bam).with_suffix(".bai"))]
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
    # Keyed by (individual, aligner). One pre-aligned BAM per aligner: the eight
    # callsets of the ensemble stay independent. Reusing a single alignment for both
    # aligners would make Jasmine count the same evidence twice (see rules/merging.smk).
    input_bams = {}
    input_bais = {}
    input_fastqs = {}
    for individual in individuals:
        rows = rows_of(individual)
        for aligner in aligners:
            declared = [bam for bam in rows[f"bam_{aligner}"] if bam]
            if len(declared) > 1:
                raise WorkflowError(
                    "Individual '{}' declares {} bam_{} files in {}. Declare exactly "
                    "one BAM per aligner per individual.".format(
                        individual, len(declared), aligner, config["samples"]
                    )
                )
            bam = check_readable_file(
                declared[0] if declared else "",
                f"bam_{aligner} file for individual '{individual}'",
            )
            input_bams[(individual, aligner)] = bam
            input_bais[(individual, aligner)] = resolve_bam_index(bam)

        # Reads are still needed: SVJedi-graph genotypes by mapping reads onto a
        # variation graph, which a linear BAM cannot substitute for.
        input_fastqs[individual] = [
            check_readable_file(f, f"fastq file for individual '{individual}'")
            for f in rows["fastq"]
            if f
        ]
        if not input_fastqs[individual]:
            raise WorkflowError(
                "start_from_bam is set but no 'fastq' file was declared for individual "
                "'{}' in {}. Reads are required to genotype the SVs with "
                "SVJedi-graph.".format(individual, config["samples"])
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
