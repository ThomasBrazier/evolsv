"""Static checks on the conda environment files referenced by the workflow.

These run on every push and cost nothing. They are the cheap gate in front of the two
expensive CI jobs: `CondaSolve` (dependency resolution on every push) and `CondaInstall`
(a genuine install, weekly). Nothing here proves an environment is installable -- only
that the workflow points at a file that exists and is a well-formed env spec.

`snakemake -n` does not read conda: directives at all, so without this a rule could
reference a deleted or misspelled env file and every dry-run test would still pass.
"""

import re

import pytest
import yaml

from list_envs import REPO_ROOT, referenced_envs

ENVS = referenced_envs()

# Envs whose pip: section is only exercised by the weekly CondaInstall job, since
# `micromamba create --dry-run` never runs pip. Listed here so the set stays visible.
PIP_BEARING_ENVS = {"workflow/envs/nanoplot.yaml", "workflow/envs/sniffles.yaml"}


def load(env):
    return yaml.safe_load((REPO_ROOT / env).read_text())


def test_some_envs_are_referenced():
    """Guards against the discovery silently matching nothing after a refactor."""
    assert len(ENVS) > 15, "only found {} referenced envs: {}".format(len(ENVS), ENVS)


@pytest.mark.parametrize("env", ENVS)
def test_env_file_exists(env):
    """Every conda: directive resolves to a file on disk."""
    assert (REPO_ROOT / env).is_file(), "{} is referenced by a rule but does not exist".format(env)


@pytest.mark.parametrize("env", ENVS)
def test_env_file_is_a_valid_spec(env):
    """Every env file parses as YAML and declares channels and dependencies."""
    spec = load(env)
    assert isinstance(spec, dict), "{} is not a YAML mapping".format(env)
    for key in ("channels", "dependencies"):
        assert spec.get(key), "{} declares no {}".format(env, key)


@pytest.mark.parametrize("env", ENVS)
def test_env_declares_bioconda_channel_order(env):
    """conda-forge and bioconda are both present.

    Package resolution depends on both being available; a spec that lists only one
    resolves differently on a developer machine with extra channels configured than it
    does in CI.
    """
    channels = load(env)["channels"]
    for channel in ("conda-forge", "bioconda"):
        assert channel in channels, "{} does not list the {} channel".format(env, channel)


@pytest.mark.parametrize("env", sorted(PIP_BEARING_ENVS))
def test_pip_bearing_envs_declare_pip(env):
    """An env with a pip: block must also install pip itself.

    Without an explicit `pip` dependency the pip: section is resolved against whatever
    pip conda happens to pull in, which is how these blocks break silently.
    """
    dependencies = load(env)["dependencies"]
    pip_blocks = [d for d in dependencies if isinstance(d, dict) and "pip" in d]
    assert pip_blocks, "{} is listed as pip-bearing but has no pip: block".format(env)
    assert "pip" in dependencies, "{} has a pip: block but does not depend on pip".format(env)


def test_pip_bearing_env_list_is_complete():
    """PIP_BEARING_ENVS must cover every env with a pip: block.

    The weekly CondaInstall job smoke-tests exactly these, so a new pip: block that is
    not listed here would go unverified.
    """
    found = set()
    for env in ENVS:
        dependencies = load(env)["dependencies"]
        if any(isinstance(d, dict) and "pip" in d for d in dependencies):
            found.add(env)
    assert found == PIP_BEARING_ENVS, (
        "envs with a pip: block changed: {} (update PIP_BEARING_ENVS here and the smoke "
        "checks in the CondaInstall job)".format(found)
    )


@pytest.mark.parametrize("env", ENVS)
def test_env_has_no_duplicate_dependencies(env):
    """A package listed twice is a merge artefact and can carry conflicting pins."""
    names = []
    for dependency in load(env)["dependencies"]:
        if not isinstance(dependency, str):
            continue  # the pip: block
        # Strip any channel prefix and version constraint: bioconda::samtools=1.9 -> samtools
        name = dependency.split("::")[-1]
        name = re.split(r"[=<>!\s]", name, maxsplit=1)[0]
        names.append(name)

    duplicates = sorted({n for n in names if names.count(n) > 1})
    assert not duplicates, "{} lists {} more than once".format(env, duplicates)
