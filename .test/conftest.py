"""Shared setup for the dry-run test suite.

The fixtures under .test/fixtures/ are git-ignored (.gitignore excludes *.tsv, *.fa and
*.txt repo-wide, plus .test/fixtures/ outright), so they are regenerated here rather
than committed. This makes `pytest .test/` work from a clean checkout, exactly as CI
does it.
"""

import subprocess
from pathlib import Path

import pytest

TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parent


@pytest.fixture(scope="session", autouse=True)
def fixtures():
    """Regenerate .test/fixtures/ once per session, before any dry run."""
    subprocess.run(
        ["bash", str(TEST_DIR / "setup_fixtures.sh")],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
