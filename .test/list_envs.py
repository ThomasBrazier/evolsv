#!/usr/bin/env python3
"""List the conda environment files actually referenced by a workflow rule.

Single source of truth for both the CI env-matrix jobs and .test/test_conda_envs.py,
so the two cannot drift apart.

Envs are derived from the rules rather than globbed over workflow/envs/: new envs are
picked up automatically, and unreferenced legacy files (vapor, sratools, datasets,
gnuplot, pysam) do not fail the build.

Prints a JSON array of repo-relative paths, which is the form the GitHub Actions
`strategy.matrix` needs.
"""

import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Matches a conda: directive and its quoted path, e.g.
#     conda:
#         "../envs/bcftools.yaml"
CONDA_DIRECTIVE = re.compile(r"conda:\s*[\"']([^\"']+)[\"']", re.MULTILINE)


def referenced_envs():
    """Return the sorted, de-duplicated env files referenced by any rule."""
    envs = set()
    for smk in sorted(REPO_ROOT.glob("workflow/**/*.smk")) + [REPO_ROOT / "workflow" / "Snakefile"]:
        if not smk.is_file():
            continue
        text = smk.read_text()
        for match in CONDA_DIRECTIVE.finditer(text):
            # Paths in a conda: directive are relative to the file declaring the rule.
            resolved = (smk.parent / match.group(1)).resolve()
            envs.add(resolved.relative_to(REPO_ROOT).as_posix())
    return sorted(envs)


def main():
    envs = referenced_envs()
    if not envs:
        print("no conda: directives found - has the rule layout changed?", file=sys.stderr)
        return 1
    print(json.dumps(envs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
