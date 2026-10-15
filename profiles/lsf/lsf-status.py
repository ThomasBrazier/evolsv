#!/usr/bin/env python3
"""Job-status probe for the cluster-generic executor on IBM Spectrum LSF.

Snakemake invokes this as `lsf-status.py <jobid>` and requires exactly one of
`running`, `success` or `failed` on stdout. Anything else, or a non-zero exit, makes
cluster-generic treat the job as failed -- so every path here prints one of the three
and exits 0.

Only the standard library is used: this runs under whichever interpreter is on the
submitting host's PATH, which is not necessarily the pipeline's conda environment.

`bjobs -noheader -o stat` reports PEND, PSUSP, RUN, USUSP, SSUSP, WAIT, PROV (live),
DONE (exit 0), EXIT (non-zero exit, or killed e.g. by a runtime/memory limit) and
UNKWN (sbatchd lost contact with the host; the job may still be running). Finished jobs
leave `bjobs` after CLEAN_PERIOD (default 1 h), so a job that `bjobs` no longer knows is
looked up with `bjobs -a`, then `bhist`. If it is in none of them, it is reported as
`success` and Snakemake's own output-file check decides -- a spurious `failed` would
abandon completed work, and a spurious `running` would hang the DAG forever.
"""

import subprocess
import sys
import time

ATTEMPTS = 5
BACKOFF_SECONDS = 3

LIVE = {"PEND", "PSUSP", "RUN", "USUSP", "SSUSP", "WAIT", "PROV", "UNKWN"}


def run(command):
    """Return stdout for `command`, or None if it could not be run or failed."""
    try:
        result = subprocess.run(
            command, capture_output=True, text=True, timeout=60, check=False
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0:
        return None
    return result.stdout


def state(jobid):
    """Return the LSF state of `jobid`, or None if LSF does not know it."""
    for command in (
        ["bjobs", "-noheader", "-o", "stat", jobid],
        ["bjobs", "-a", "-noheader", "-o", "stat", jobid],
    ):
        text = run(command)
        if text and text.split():
            return text.split()[0]
    return None


def status(jobid):
    for attempt in range(ATTEMPTS):
        current = state(jobid)
        if current in LIVE:
            return "running"
        if current == "DONE":
            return "success"
        if current == "EXIT":
            return "failed"
        if attempt < ATTEMPTS - 1:
            time.sleep(BACKOFF_SECONDS)
    # Unknown to bjobs: see the module docstring for why this is not "failed".
    return "success"


def main():
    if len(sys.argv) != 2:
        print("usage: lsf-status.py <jobid>", file=sys.stderr)
        print("failed")
        return 0
    print(status(sys.argv[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
