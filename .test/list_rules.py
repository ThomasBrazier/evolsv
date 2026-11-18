#!/usr/bin/env python3
"""List the rules defined by the workflow, and which of them need cluster resources.

Single source of truth for .test/test_profiles.py, so the scheduler profiles and the
workflow cannot drift apart. Rules are derived from the source rather than from a
hand-kept list: that is the whole point, since the drift this guards against
(profiles/slurm/config.yaml naming rules that no longer existed) came from a hand-kept
list going stale in silence.

Silence is the operative word. Snakemake looks up profile keys with
`if name in overwrite_threads`, so a set-threads or set-resources entry naming a rule
that does not exist is discarded without a warning, and a rule missing from the profile
falls back to one thread without a warning either.

Prints a JSON array of submittable rule names.
"""

import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Matches a rule or checkpoint definition, e.g.
#     rule minimap2:
# The leading [ \t]* is load-bearing: `download_sra` and `merge_fastq` are indented
# inside `if not bam_mode:` in workflow/rules/data_download.smk.
RULE_DIRECTIVE = re.compile(r"^[ \t]*(?:rule|checkpoint)\s+(\w+)\s*:", re.MULTILINE)

# Matches the names inside a localrules: block, e.g.
#     localrules:
#         all,
LOCALRULES_DIRECTIVE = re.compile(r"^localrules:\s*\n((?:\s+\w+,?\s*\n)+)", re.MULTILINE)


def _source_files():
    """The same file set .test/list_envs.py scans."""
    return sorted(REPO_ROOT.glob("workflow/**/*.smk")) + [REPO_ROOT / "workflow" / "Snakefile"]


def local_rules():
    """Return the rules declared local, which are never submitted to a scheduler."""
    names = set()
    for source in _source_files():
        if not source.is_file():
            continue
        for match in LOCALRULES_DIRECTIVE.finditer(source.read_text()):
            names.update(name.strip(" ,") for name in match.group(1).split() if name.strip(" ,"))
    return names


def all_rules():
    """Return every rule name defined anywhere in the workflow."""
    names = set()
    for source in _source_files():
        if not source.is_file():
            continue
        names.update(RULE_DIRECTIVE.findall(source.read_text()))
    return names


def submittable_rules():
    """Return the rules that a scheduler profile has to size.

    Rules defined inside `if not bam_mode:` are included: they run on the default
    (FASTQ) entry point and genuinely need resources there. Only local rules are
    excluded.
    """
    return all_rules() - local_rules()


def main():
    rules = sorted(submittable_rules())
    if not rules:
        print("no rule definitions found - has the rule layout changed?", file=sys.stderr)
        return 1
    print(json.dumps(rules))
    return 0


if __name__ == "__main__":
    sys.exit(main())
