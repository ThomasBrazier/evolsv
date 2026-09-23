"""Dry-run tests for the EvolSV workflow.

Scope: these tests check that the DAG *builds* for a range of dataset shapes, and that
bad input is still rejected at DAG-construction time. They never execute a rule, so a
green run means "the workflow parses and the guards fire", not "the results are correct".

Run from the repo root with `pytest .test/`. CI runs exactly the same command.

Sibling suite: .test/test_profiles.py checks that the scheduler profiles under profiles/
and the rules defined in workflow/ agree about which rules exist. That one never invokes
snakemake; this one always does.
"""

import re
import subprocess
from pathlib import Path

import pytest

TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parent

# Rules that only exist on the FASTQ entry point (rules/data_qc.smk + rules/mapping.smk).
# `filter_reads_chopper` is the real rule name; there is no rule called `chopper`.
FASTQ_ENTRY_RULES = (
    "download_sra",
    "longqc",
    "merge_fastq",
    "filter_reads_chopper",
    "minimap2",
    "ngmlr",
)

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
    ("local-metadata", ".test/config_local_metadata.yaml"),
    ("local-metadata-no-assembly-report", ".test/config_local_metadata_no_asm.yaml"),
    ("no-scaffold-exclusion", ".test/config_no_scaffold_exclusion.yaml"),
    ("bigtmp", ".test/config_bigtmp.yaml"),
    ("config-test", "config/config_test.yaml"),
    ("ont", ".test/config_ont.yaml"),
    ("multi-individual", ".test/config_multi.yaml"),
    ("bam-mode-multi-individual", ".test/config_bam_multi.yaml"),
    ("bam-mode-no-fastq", ".test/config_bam_no_fastq.yaml"),
    ("bam-mode-mixed-fastq", ".test/config_bam_mixed_fastq.yaml"),
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


def test_bam_mode_extracts_reads_when_no_fastq_is_declared():
    """A blank `fastq` column is supported: the reads come back out of the minimap2 BAM."""
    result = run_dryrun(".test/config_bam_no_fastq.yaml")
    assert_succeeded(result)
    counts = parse_job_stats(result.stdout)

    assert counts["bam_to_fastq"] == 1
    assert "stage_fastq" not in counts, "nothing was declared to stage"
    # The reads exist for the reason the rule exists.
    assert counts["svjedigraph"] == 1
    # Secondary and supplementary records must be excluded: if they were kept, every
    # clipped segment of a split long-read alignment would enter the callers as an extra
    # read. Losing the flag changes the science without breaking the DAG, so it is
    # asserted on here.
    assert "samtools fastq -F 0x900" in result.stdout


def test_bam_mode_mixes_declared_and_extracted_reads():
    """Both producers of {genome}_filtered.fastq.gz can coexist in one DAG.

    stage_fastq and bam_to_fastq write the same path, so this is the case that raises
    AmbiguousRuleException if their wildcard_constraints stop being disjoint.
    """
    result = run_dryrun(".test/config_bam_mixed_fastq.yaml")
    assert_succeeded(result)
    counts = parse_job_stats(result.stdout)

    assert counts["stage_fastq"] == 1
    assert counts["bam_to_fastq"] == 1


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


def test_reference_names_are_checked_before_the_autosome_split():
    """The contig-name check runs in every reference mode, once for the shared genome."""
    for configfile in (None, ".test/config_local_ref.yaml", ".test/config_local_metadata.yaml"):
        result = run_dryrun(configfile)
        assert_succeeded(result)
        assert parse_job_stats(result.stdout)["check_reference_names"] == 1, configfile


def test_local_metadata_downloads_nothing_from_ncbi():
    """reference_fasta + sequence_report: rule stage_genome replaces download_genome."""
    result = run_dryrun(".test/config_local_metadata.yaml")
    assert_succeeded(result)
    counts = parse_job_stats(result.stdout)
    assert "download_genome" not in counts
    assert counts["stage_genome"] == 1
    assert counts["stage_assembly_data_report"] == 1
    assert "datasets download" not in result.stdout
    assert re.search(r"Using the local reference FASTA .*ref\.fna", result.stdout)

    # assembly_data_report is optional: no staging rule, and the report does not wait for it.
    result = run_dryrun(".test/config_local_metadata_no_asm.yaml")
    assert_succeeded(result)
    counts = parse_job_stats(result.stdout)
    assert counts["stage_genome"] == 1
    assert "stage_assembly_data_report" not in counts
    assert "download_genome" not in counts

    # A local FASTA alone still takes its metadata from NCBI.
    counts = parse_job_stats(run_dryrun(".test/config_local_ref.yaml").stdout)
    assert counts["download_genome"] == 1
    assert "stage_genome" not in counts


def test_default_technology_is_hifi():
    """The shipped default stays PacBio HiFi, so the preset table cannot flip it silently."""
    result = run_dryrun()
    assert_succeeded(result)
    for flag in ("-ax map-hifi", "--presets pacbio", "PL:PACBIO"):
        assert flag in result.stdout, "{!r} missing from the default shell commands".format(flag)


def test_adapter_filtering_matches_the_technology():
    """Each run goes through the adapter tool of its technology, then runs are merged."""
    # config/samples.tsv declares two SRA runs for the one individual.
    result = run_dryrun()
    assert_succeeded(result)
    counts = parse_job_stats(result.stdout)
    assert counts["hifiadapterfilt"] == 2
    assert "porechop_abi" not in counts
    assert "_sra.filt.fastq.gz" in result.stdout
    assert "_sra.porechop.fastq.gz" not in result.stdout

    result = run_dryrun(".test/config_ont.yaml")
    assert_succeeded(result)
    counts = parse_job_stats(result.stdout)
    assert counts["porechop_abi"] == 2
    assert "hifiadapterfilt" not in counts
    assert "_sra.porechop.fastq.gz" in result.stdout
    assert "_sra.filt.fastq.gz" not in result.stdout
    assert "--ab_initio" in result.stdout


def test_longqc_preset_follows_the_technology():
    """LongQC runs once per SRA run, with the -x preset of the sequencing technology."""
    result = run_dryrun()
    assert_succeeded(result)
    # config/samples.tsv declares two SRA runs for the one individual.
    assert parse_job_stats(result.stdout)["longqc"] == 2
    assert "-x pb-hifi" in result.stdout
    # profiles/ci has 1 core; LongQC exits below 4 CPUs, so -p must be raised to 4.
    assert "-p $(( 1 < 4 ? 4 : 1 ))" in result.stdout

    result = run_dryrun(".test/config_ont.yaml")
    assert_succeeded(result)
    assert "-x ont-ligation" in result.stdout
    assert "-x pb-hifi" not in result.stdout

    result = subprocess.run(
        [
            "snakemake", "-s", "workflow/Snakefile", "-n", "-p",
            "--profile", "profiles/ci",
            "--configfile", ".test/config_ont.yaml",
            "--config", "longqc_preset=ont-rapid",
        ],
        cwd=REPO_ROOT, capture_output=True, text=True, check=False,
    )
    assert_succeeded(result)
    assert "-x ont-rapid" in result.stdout
    assert "-x ont-ligation" not in result.stdout


def test_porechop_ab_initio_can_be_disabled():
    """porechop_ab_initio: false keeps Porechop_ABI on its adapter database only.

    `--config key=false` passes the string "false", which is truthy in Python, so both
    spellings are checked. Any other value must stop the run.
    """

    def run_with(value):
        return subprocess.run(
            [
                "snakemake", "-s", "workflow/Snakefile", "-n", "-p",
                "--profile", "profiles/ci",
                "--configfile", ".test/config_ont.yaml",
                "--config", f"porechop_ab_initio={value}",
            ],
            cwd=REPO_ROOT, capture_output=True, text=True, check=False,
        )

    for value in ("false", "False"):
        result = run_with(value)
        assert_succeeded(result)
        assert "porechop_abi -i" in result.stdout
        assert "--ab_initio" not in result.stdout, f"--ab_initio still set with {value!r}"

    result = run_with("maybe")
    assert result.returncode != 0
    assert "must be true or false" in result.stdout + result.stderr


def test_ont_presets_reach_the_shell_commands():
    """Besides hifiadapterfilt, sequencing_technology only sets config defaults.

    So assert on the commands.

    One flag per tool the preset owns: minimap2, ngmlr and cuteSV. cuteSV is the one
    that would otherwise silently keep HiFi clustering on an ONT run.
    """
    result = run_dryrun(".test/config_ont.yaml")
    assert_succeeded(result)

    for flag in ("-ax map-ont", "--presets ont", "PL:ONT", "--max_cluster_bias_INS 100"):
        assert flag in result.stdout, "{!r} missing from the ONT shell commands".format(flag)

    # The HiFi values must be gone, not merely joined by the ONT ones.
    for flag in ("-ax map-hifi", "--presets pacbio", "--max_cluster_bias_INS 1000"):
        assert flag not in result.stdout, "{!r} still present in an ONT run".format(flag)


def test_explicit_config_key_overrides_the_preset():
    """A key set by hand wins over the technology preset (config.setdefault, not assignment)."""
    result = subprocess.run(
        [
            "snakemake", "-s", "workflow/Snakefile", "-n", "-p",
            "--profile", "profiles/ci",
            "--configfile", ".test/config_ont.yaml",
            "--config", "minimap_ax=asm20",
        ],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert_succeeded(result)
    assert "-ax asm20" in result.stdout, "the explicit minimap_ax did not override the ont preset"
    # ngmlr is not overridden, so the rest of the preset must still apply.
    assert "--presets ont" in result.stdout


# Rules that run once per individual, and rules that run once for the shared reference.
PER_INDIVIDUAL_RULES = ("sample_ids", "jasmine", "svjedigraph", "final_vcf", "final_report", "truvari_grm")
SHARED_REFERENCE_RULES = (
    "download_genome",
    "genmap",
    "mappability_bed",
    "autosomes_sexchromosomes",
    "bwa_index",
)
MULTI_WDIR = ".test/fixtures/data/GCA_947247005.1"
INDIVIDUALS = ("SAMEA8724893", "SAMEA0000002")


def test_multi_individual_fastq_mode_runs_each_individual_once():
    """Two individuals: one result set each, one copy of the shared reference work."""
    counts = parse_job_stats(run_dryrun(".test/config_multi.yaml").stdout)

    for rule in PER_INDIVIDUAL_RULES + ("merge_fastq", "minimap2", "ngmlr"):
        assert counts.get(rule) == 2, f"{rule} should run once per individual"
    for rule in SHARED_REFERENCE_RULES:
        assert counts.get(rule) == 1, f"{rule} should run once for all individuals"
    # Three SRA runs in total: two for the first individual, one for the second.
    assert counts["download_sra"] == 3
    assert counts["fastqc"] == 3


def test_multi_individual_bam_mode_runs_each_individual_once():
    """BAM mode stages one BAM per aligner per individual."""
    counts = parse_job_stats(run_dryrun(".test/config_bam_multi.yaml").stdout)

    assert counts["stage_bam"] == 4
    assert counts["stage_fastq"] == 2
    for rule in PER_INDIVIDUAL_RULES:
        assert counts.get(rule) == 2, f"{rule} should run once per individual"
    for rule in SHARED_REFERENCE_RULES:
        assert counts.get(rule) == 1, f"{rule} should run once for all individuals"


def test_reads_are_never_merged_across_individuals():
    """Each individual's merged FASTQ is built from its own SRA runs only."""
    result = run_dryrun(".test/config_multi.yaml")
    assert_succeeded(result)

    merges = [
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip().startswith("cat ") and line.strip().endswith("/GCA_947247005.1.fastq.gz")
    ]
    first_dir = f"{MULTI_WDIR}/{INDIVIDUALS[0]}/fastq"
    second_dir = f"{MULTI_WDIR}/{INDIVIDUALS[1]}/fastq"
    # Runs are concatenated in sample-sheet order, after HiFiAdapterFilt (default hifi).
    first = (
        f"cat {first_dir}/ERR10287556_sra.filt.fastq.gz {first_dir}/ERR10287555_sra.filt.fastq.gz"
        f" > {first_dir}/GCA_947247005.1.fastq.gz"
    )
    second = (
        f"cat {second_dir}/ERR00000001_sra.filt.fastq.gz"
        f" > {second_dir}/GCA_947247005.1.fastq.gz"
    )
    assert first in merges, f"first individual's merge_fastq command not found in:\n{merges}"
    assert second in merges, f"second individual's merge_fastq command not found in:\n{merges}"


def test_each_individual_gets_its_own_read_group_and_final_vcf():
    """Read groups and final outputs are named after the individual, not the first row."""
    result = run_dryrun(".test/config_multi.yaml")
    assert_succeeded(result)

    for individual in INDIVIDUALS:
        assert f"SM:{individual}" in result.stdout, f"minimap2 read group missing for {individual}"
        assert f"--rg-sm {individual}" in result.stdout, f"ngmlr read group missing for {individual}"
        final = f"{MULTI_WDIR}/{individual}/GCA_947247005.1_final.vcf.gz"
        assert final in result.stdout, f"{final} is not scheduled"


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
    (
        "no-ngmlr",
        ".test/config_bad_no_ngmlr.yaml",
        "no bam_ngmlr file for individual 'SAMEA8724893' was declared",
    ),
    ("bad-technology", ".test/config_bad_tech.yaml", "Unknown sequencing_technology"),
    ("bad-longqc-preset", ".test/config_bad_longqc_preset.yaml", "Unknown longqc_preset 'pb-hifii'"),
    ("two-genomes", ".test/config_bad_two_genomes.yaml", "must use the same reference genome"),
    (
        "duplicate-bam",
        ".test/config_bad_duplicate_bam.yaml",
        "Individual 'SAMEA8724893' declares 2 bam_minimap2 files",
    ),
    ("bad-sample-name", ".test/config_bad_sample_name.yaml", "Invalid sample_name 'bad/name'"),
    (
        "report-without-fasta",
        ".test/config_bad_report_without_fasta.yaml",
        "sequence_report is set but reference_fasta is not",
    ),
    (
        "assembly-report-without-sequence-report",
        ".test/config_bad_asm_without_report.yaml",
        "assembly_data_report is set but sequence_report is not",
    ),
    (
        "missing-sequence-report",
        ".test/config_bad_missing_report.yaml",
        (
            "The sequence_report file '.test/fixtures/bad/does_not_exist_sequence_report.jsonl' "
            "declared in the config file does not exist"
        ),
    ),
]


@pytest.mark.parametrize(
    "name,configfile,message", NEGATIVE_CASES, ids=[c[0] for c in NEGATIVE_CASES]
)
def test_bad_input_is_rejected(name, configfile, message):
    """Malformed input fails the dry run with an actionable message."""
    result = run_dryrun(configfile)

    assert result.returncode != 0, (
        "expected the dry run to fail for '{}' but it succeeded:\n{}".format(name, result.stdout)
    )
    combined = result.stdout + result.stderr
    assert message in combined, (
        "expected {!r} in the failure output for '{}', got:\n{}".format(message, name, combined)
    )
