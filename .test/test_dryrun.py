"""Dry-run tests for the EvolSV workflow.

Scope: these tests check that the DAG *builds* for a range of dataset shapes, and that
bad input is still rejected at DAG-construction time. They never execute a rule, so a
green run means "the workflow parses and the guards fire", not "the results are correct".

Run from the repo root with `pytest .test/`. CI runs exactly the same command.
"""

import re
import subprocess
from pathlib import Path

import pytest

TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parent

# Rules that only exist on the FASTQ entry point (rules/data_qc.smk + rules/mapping.smk).
# `filter_reads_chopper` is the real rule name; there is no rule called `chopper`.
FASTQ_ENTRY_RULES = ("download_sra", "merge_fastq", "filter_reads_chopper", "minimap2", "ngmlr")

# One terminal rule behind each group of `rule all` inputs in workflow/Snakefile.
# If any of these disappears from the DAG, a target has silently stopped being built.
RULE_ALL_TERMINALS = (
    "final_vcf",  # {genome}_final.vcf.gz
    "light_vcf",  # {genome}_final_light.vcf.gz
    "final_report",  # {genome}_finalQC.html
    "samplot_plot",  # samplot/*/{DUP,INV,DEL}/index.html
    "callable_bed",  # callability/{genome}_{aligner}_callable_mappable.bed
    "add_mappability",  # callability/{genome}_callable_mappable.bed
    "truvari_grm",  # annotate_grm/{genome}_grm.jl
    "samtools_coverage",  # mapping_QC/{genome}_{aligner}_coverage.tsv
)


def run_dryrun(configfile=None):
    """Run `snakemake -n` from the repo root and return the CompletedProcess.

    Always invoked from REPO_ROOT: workflow/Snakefile hardcodes
    `configfile: "config/config.yaml"` and several rules hardcode workflow/scripts/...,
    so `--directory` must not be used.

    profiles/ci supplies `default-resources: cpus_per_task`, which many rules
    interpolate into shell: as {resources.cpus_per_task}. Without it, --printshellcmds
    fails to format those rules.
    """
    command = [
        "snakemake",
        "-s", "workflow/Snakefile",
        "-n",
        "-p",
        "--profile", "profiles/ci",
    ]
    if configfile is not None:
        command += ["--configfile", configfile]

    return subprocess.run(
        command, cwd=REPO_ROOT, capture_output=True, text=True
    )


def parse_job_stats(output):
    """Parse snakemake's `Job stats:` table into {rule_name: count}.

    Reading the table rather than grepping the whole log means an assertion cannot be
    satisfied by the same word turning up in a file path or in the `Reasons:` section.

    Snakemake prints the table twice in a dry run (once up front, once at the end); the
    two are identical, so the first is used.
    """
    lines = output.splitlines()
    for index, line in enumerate(lines):
        if line.strip() == "Job stats:":
            break
    else:
        raise AssertionError("no 'Job stats:' table in snakemake output:\n" + output)

    counts = {}
    # Skip the "job    count" header and the dashed rule beneath it.
    for line in lines[index + 3:]:
        match = re.fullmatch(r"(\S+)\s+(\d+)", line.strip())
        if not match:
            break
        rule, count = match.groups()
        if rule == "total":
            break
        counts[rule] = int(count)

    assert counts, "parsed an empty 'Job stats:' table from:\n" + output
    return counts


def assert_succeeded(result):
    """Fail with the full snakemake output rather than a bare exit code."""
    assert result.returncode == 0, (
        "snakemake dry run failed with exit code {}\n"
        "--- stdout ---\n{}\n--- stderr ---\n{}".format(
            result.returncode, result.stdout, result.stderr
        )
    )


# ---------------------------------------------------------------------------
# Positive cases: different dataset shapes must all build a DAG
# ---------------------------------------------------------------------------

POSITIVE_CASES = [
    # (case name, configfile or None for the repo default)
    ("default", None),
    ("single-run", ".test/config_single_run.yaml"),
    ("bam-mode", ".test/config_bam.yaml"),
    ("local-reference", ".test/config_local_ref.yaml"),
    ("no-scaffold-exclusion", ".test/config_no_scaffold_exclusion.yaml"),
    ("bigtmp", ".test/config_bigtmp.yaml"),
    ("config-test", "config/config_test.yaml"),
]


@pytest.mark.parametrize("name,configfile", POSITIVE_CASES, ids=[c[0] for c in POSITIVE_CASES])
def test_dag_builds(name, configfile):
    """Every supported dataset shape produces a complete DAG ending in rule all."""
    result = run_dryrun(configfile)
    assert_succeeded(result)
    assert "all" in parse_job_stats(result.stdout)


def test_default_uses_the_fastq_entry_point():
    """The default config downloads reads and aligns them with both aligners."""
    counts = parse_job_stats(run_dryrun().stdout)

    for rule in FASTQ_ENTRY_RULES:
        assert rule in counts, "{} missing from the default DAG".format(rule)
    # config/samples.tsv declares two SRA runs for the one individual.
    assert counts["download_sra"] == 2
    # One alignment per aligner, then the eight caller x aligner callsets are merged.
    assert counts["minimap2"] == 1
    assert counts["ngmlr"] == 1
    assert counts["jasmine"] == 1


def test_sample_sheet_drives_the_download_count():
    """A one-run sheet yields one download_sra job, not a hardcoded two."""
    counts = parse_job_stats(run_dryrun(".test/config_single_run.yaml").stdout)
    assert counts["download_sra"] == 1


def test_bam_mode_stages_bams_and_skips_alignment():
    """start_from_bam replaces download/QC/alignment with staging, one BAM per aligner."""
    counts = parse_job_stats(run_dryrun(".test/config_bam.yaml").stdout)

    # One BAM per aligner: reusing a single alignment would make jasmine count the same
    # evidence twice (see rules/bam_input.smk).
    assert counts["stage_bam"] == 2
    assert counts["stage_fastq"] == 1
    for rule in FASTQ_ENTRY_RULES:
        assert rule not in counts, "{} should be skipped in BAM mode".format(rule)


def test_bam_mode_still_runs_alignment_qc_and_calling():
    """Skipping alignment must not skip anything downstream of it."""
    counts = parse_job_stats(run_dryrun(".test/config_bam.yaml").stdout)

    for rule in ("samtools_coverage", "mosdepth_summary", "svim", "sniffles", "cutesv", "debreak"):
        assert rule in counts, "{} missing from the BAM-mode DAG".format(rule)


def test_local_reference_reaches_the_shell_command():
    """reference_fasta only feeds a params:, so assert on the printed shell command."""
    result = run_dryrun(".test/config_local_ref.yaml")
    assert_succeeded(result)
    assert re.search(r"Using the local reference FASTA .*ref\.fna", result.stdout), (
        "the local reference FASTA never reached rule download_genome's shell command"
    )


@pytest.mark.parametrize("configfile", [None, ".test/config_bam.yaml"], ids=["default", "bam-mode"])
def test_rule_all_targets_are_reachable(configfile):
    """Every group of `rule all` inputs has its terminal rule scheduled, in both modes."""
    counts = parse_job_stats(run_dryrun(configfile).stdout)

    missing = [rule for rule in RULE_ALL_TERMINALS if rule not in counts]
    assert not missing, "rule all targets no longer reachable via: {}".format(missing)


# ---------------------------------------------------------------------------
# Negative cases: bad input must fail at DAG construction, with a usable message
# ---------------------------------------------------------------------------

# The point of the guards in workflow/rules/common.smk is to fail before an expensive
# run starts. Each case asserts the specific message too, so a test cannot pass by
# failing for an unrelated reason.
NEGATIVE_CASES = [
    ("missing-bam", ".test/config_bad_missing_bam.yaml", "does not exist"),
    ("empty-bam", ".test/config_bad_empty_bam.yaml", "is empty"),
    ("no-index", ".test/config_bad_no_index.yaml", "No BAI index found"),
    ("no-fastq", ".test/config_bad_no_fastq.yaml", "Reads are required to genotype the SVs"),
    ("no-ngmlr", ".test/config_bad_no_ngmlr.yaml", "no bam_ngmlr file was declared"),
]


@pytest.mark.parametrize(
    "name,configfile,message", NEGATIVE_CASES, ids=[c[0] for c in NEGATIVE_CASES]
)
def test_bad_input_is_rejected(name, configfile, message):
    """Malformed start_from_bam input fails the dry run with an actionable message."""
    result = run_dryrun(configfile)

    assert result.returncode != 0, (
        "expected the dry run to fail for '{}' but it succeeded:\n{}".format(name, result.stdout)
    )
    combined = result.stdout + result.stderr
    assert message in combined, (
        "expected {!r} in the failure output for '{}', got:\n{}".format(message, name, combined)
    )
