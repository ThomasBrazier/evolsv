"""Consistency tests for the scheduler profiles under profiles/.

Scope: these read the profile YAML and the rule definitions and check that the two
agree. They never invoke snakemake, so they are fast and do not need the PBS or SGE
executor plugins installed. Complementary to .test/test_dryrun.py, which builds the
actual DAG.

The invariants exist because every failure mode here is silent. Snakemake discards a
set-threads/set-resources key naming a rule that does not exist, and gives a rule that
is missing from the profile one thread -- neither produces a warning. Before these
tests, profiles/slurm/config.yaml had six keys for rules that had been renamed away,
including a `chopper: 16` whose real rule `filter_reads_chopper` was separately pinned
to 1, so read filtering had been running single-threaded.
"""

import os
import re
from pathlib import Path

import pytest
import yaml

import list_rules

TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parent
PROFILES_DIR = REPO_ROOT / "profiles"

# profiles/ci is deliberately minimal: it has no executor and no set-threads, because it
# only exists to make dry runs format. It is exempt from the sizing invariants.
SIZING_PROFILES = sorted(
    path
    for path in PROFILES_DIR.glob("*/config.yaml")
    if path.parent.name != "ci"
)

ALL_PROFILES = sorted(PROFILES_DIR.glob("*/config.yaml"))

# Resource keys that only mean something to one executor. A key from the wrong engine is
# not an error to snakemake -- it is passed through as an unused resource -- so
# copy-pasting between profiles fails silently. `qos` has no prefix but is SLURM-only.
ENGINE_KEY_PREFIXES = {
    "slurm": ("slurm_", "qos"),
    "pbs": ("pbs_",),
    "sge": ("sge_",),
}

# Rules that interpolate {resources.cpus_per_task} into shell:. Derived from the source
# so this cannot go stale, and used to explain *why* cpus_per_task must be defined.
CPUS_PER_TASK_INTERPOLATION = re.compile(r"\{resources\.cpus_per_task\}")


def load(path):
    return yaml.safe_load(path.read_text())


def profile_id(path):
    return path.parent.name


def resource_keys(profile):
    """Every resource name used anywhere in a profile's default- or set-resources."""
    keys = set(profile.get("default-resources") or {})
    for rule_resources in (profile.get("set-resources") or {}).values():
        keys.update(rule_resources or {})
    return keys


def rules_interpolating_cpus_per_task():
    sources = sorted(REPO_ROOT.glob("workflow/**/*.smk")) + [REPO_ROOT / "workflow" / "Snakefile"]
    hits = set()
    for source in sources:
        if not source.is_file():
            continue
        current = None
        for line in source.read_text().splitlines():
            match = re.match(r"^[ \t]*(?:rule|checkpoint)\s+(\w+)\s*:", line)
            if match:
                current = match.group(1)
            elif current and CPUS_PER_TASK_INTERPOLATION.search(line):
                hits.add(current)
    return hits


def test_rule_discovery_is_not_empty():
    """Guard against a regex change silently turning every other test into a no-op."""
    rules = list_rules.submittable_rules()
    assert len(rules) > 40, f"only found {len(rules)} rules; has the rule layout changed?"
    assert "all" not in rules, "rule all is local and must not require cluster resources"


@pytest.mark.parametrize("path", SIZING_PROFILES, ids=profile_id)
def test_thread_map_covers_every_rule(path):
    """Every submittable rule has an explicit thread count, and no key is stale."""
    declared = set(load(path).get("set-threads") or {})
    rules = list_rules.submittable_rules()

    missing = sorted(rules - declared)
    stale = sorted(declared - rules)
    assert not missing, (
        f"{profile_id(path)}: rules with no set-threads entry (they silently get 1 "
        f"thread): {missing}"
    )
    assert not stale, (
        f"{profile_id(path)}: set-threads names rules that do not exist (silently "
        f"ignored): {stale}"
    )


@pytest.mark.parametrize("path", ALL_PROFILES, ids=profile_id)
def test_no_stale_set_resources_keys(path):
    """set-resources is opt-in, so a subset -- but every key must name a real rule."""
    declared = set(load(path).get("set-resources") or {})
    stale = sorted(declared - list_rules.submittable_rules())
    assert not stale, (
        f"{profile_id(path)}: set-resources names rules that do not exist (silently "
        f"ignored): {stale}"
    )


@pytest.mark.parametrize("path", ALL_PROFILES, ids=profile_id)
def test_cpus_per_task_is_defined(path):
    """Without it, the rules that interpolate it into shell: fail to format."""
    defaults = load(path).get("default-resources") or {}
    consumers = sorted(rules_interpolating_cpus_per_task())
    assert "cpus_per_task" in defaults, (
        f"{profile_id(path)}: default-resources must define cpus_per_task; "
        f"{len(consumers)} rules interpolate it into shell: {consumers}"
    )


@pytest.mark.parametrize("path", SIZING_PROFILES, ids=profile_id)
def test_cpus_per_task_is_derived_from_threads(path):
    """set-threads must stay the single place a thread count is declared.

    A literal cpus_per_task escapes the --cores clamp while {threads} does not, so the
    two can disagree without anything failing -- which is how svim and cutesv ended up
    with set-threads: 8 and a commented-out cpus_per_task.
    """
    profile = load(path)
    assert (profile.get("default-resources") or {}).get("cpus_per_task") == "threads", (
        f"{profile_id(path)}: default-resources.cpus_per_task must be the expression "
        f"`threads`, so it tracks set-threads"
    )
    overrides = sorted(
        rule
        for rule, resources in (profile.get("set-resources") or {}).items()
        if "cpus_per_task" in (resources or {})
    )
    assert not overrides, (
        f"{profile_id(path)}: set cpus_per_task via set-threads, not set-resources, "
        f"or the two can disagree: {overrides}"
    )


@pytest.mark.parametrize("path", SIZING_PROFILES, ids=profile_id)
def test_cores_is_at_least_the_largest_thread_count(path):
    """`--cores N` clamps every rule's threads to N, even in cluster mode."""
    profile = load(path)
    threads = profile.get("set-threads") or {}
    largest = max(int(value) for value in threads.values())
    cores = profile.get("cores")
    assert cores is not None, f"{profile_id(path)}: must set cores:, or threads clamp to 1"
    assert int(cores) >= largest, (
        f"{profile_id(path)}: cores={cores} silently clamps set-threads (max {largest})"
    )


@pytest.mark.parametrize("path", SIZING_PROFILES, ids=profile_id)
def test_engine_specific_keys_do_not_leak(path):
    """A resource key from another engine is passed through unused, not rejected."""
    name = profile_id(path)
    assert name in ENGINE_KEY_PREFIXES, (
        f"new profile {name!r}: add its resource-key prefixes to ENGINE_KEY_PREFIXES"
    )
    foreign = tuple(
        prefix
        for engine, prefixes in ENGINE_KEY_PREFIXES.items()
        if engine != name
        for prefix in prefixes
    )
    leaked = sorted(key for key in resource_keys(load(path)) if key.startswith(foreign))
    assert not leaked, f"{name}: carries another engine's resource keys: {leaked}"


@pytest.mark.parametrize("path", SIZING_PROFILES, ids=profile_id)
def test_string_resources_are_quoted_for_eval(path):
    """Profile scalars are eval()'d, surviving as strings only via a NameError fallback.

    So a queue named `all`, `input` or `max` resolves to a Python builtin and aborts the
    run. The "'name'" form is immune.
    """
    profile = load(path)
    scopes = [("default-resources", profile.get("default-resources") or {})]
    for rule, resources in (profile.get("set-resources") or {}).items():
        scopes.append((f"set-resources.{rule}", resources or {}))

    for scope, resources in scopes:
        for key, value in resources.items():
            if not isinstance(value, str) or key == "cpus_per_task":
                continue
            # Arithmetic expressions (attempt * 20000) are meant to be eval'd.
            if re.fullmatch(r"[\w\s*/+()-]*[*/+()-][\w\s*/+()-]*", value):
                continue
            assert re.fullmatch(r"'.*'", value), (
                f"{profile_id(path)}: {scope}.{key} = {value!r} must be written "
                f'"\'{value}\'" so it is not eval()\'d as a Python name'
            )


@pytest.mark.parametrize("path", SIZING_PROFILES, ids=profile_id)
def test_cluster_generic_profiles_are_wired_up(path):
    """cluster-generic translates nothing itself, so the submit command carries it all."""
    profile = load(path)
    if profile.get("executor") != "cluster-generic":
        pytest.skip("not a cluster-generic profile")

    submit = profile.get("cluster-generic-submit-cmd")
    assert submit, f"{profile_id(path)}: cluster-generic needs cluster-generic-submit-cmd"
    assert "{threads}" in submit, (
        f"{profile_id(path)}: request {{threads}}, not {{resources.cpus_per_task}} -- "
        f"threads is the value after set-threads and after the --cores clamp"
    )
    assert "$" not in submit, (
        f"{profile_id(path)}: the profile parser runs os.path.expandvars() on scalars, "
        f"so $VAR expands on the submitting host at parse time"
    )

    # Helper scripts are referenced by bare filename; the parser rewrites them to the
    # in-profile path, which only works if the file is actually there and runnable.
    for key in ("cluster-generic-submit-cmd", "cluster-generic-status-cmd"):
        command = profile.get(key)
        if not command:
            continue
        script = path.parent / command.split()[0]
        if not script.exists():
            continue
        assert os.access(script, os.X_OK), f"{script.relative_to(REPO_ROOT)} is not executable"
