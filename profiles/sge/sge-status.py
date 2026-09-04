#!/usr/bin/env python3
"""Job-status probe for the cluster-generic executor on SGE / SoGE / UGE.

Snakemake invokes this as `sge-status.py <jobid>` and requires exactly one of
`running`, `success` or `failed` on stdout. Anything else, or a non-zero exit, makes
cluster-generic treat the job as failed -- so every path here prints one of the three
and exits 0.

Only the standard library is used: this runs under whichever interpreter is on the
submitting host's PATH, which is not necessarily the pipeline's conda environment.

SGE splits job state across two commands, and there is a window between them: `qstat`
forgets a job the moment it finishes, while `qacct` only learns about it once the
accounting file has been flushed, which on a busy cluster can take tens of seconds.
A job that is in neither is therefore polled several times before a verdict is given.
If it is still in neither after that, it is reported as `success` and Snakemake's own
output-file check decides -- a spurious `failed` would abandon completed work, and a
spurious `running` would hang the DAG forever.
"""

import subprocess
import sys
import time

# Polling budget for the qstat/qacct gap described above.
ATTEMPTS = 5
BACKOFF_SECONDS = 3


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


def parse_field(text, field):
    """Pull a `field   value` pair out of qacct -j output."""
    for line in text.splitlines():
        parts = line.split(None, 1)
        if len(parts) == 2 and parts[0] == field:
            return parts[1].strip()
    return None


def classify_accounting(text):
    """Map a qacct record to success/failed, or None if it is not decidable yet."""
    failed = parse_field(text, "failed")
    exit_status = parse_field(text, "exit_status")
    if failed is None and exit_status is None:
        return None
    # `failed` is reported as "0" or as "37 : qmaster enforced h_rt limit"; only the
    # leading code matters.
    failed_code = failed.split()[0] if failed else "0"
    if failed_code != "0":
        return "failed"
    return "success" if (exit_status or "0") == "0" else "failed"


def status(jobid):
    for attempt in range(ATTEMPTS):
        # qstat -j succeeds only while the job is still queued or running.
        if run(["qstat", "-j", jobid]) is not None:
            return "running"

        text = run(["qacct", "-j", jobid])
        if text:
            verdict = classify_accounting(text)
            if verdict is not None:
                return verdict

        if attempt < ATTEMPTS - 1:
            time.sleep(BACKOFF_SECONDS)

    # In neither qstat nor qacct. See the module docstring for why this is not "failed".
    return "success"


def main():
    if len(sys.argv) != 2:
        print("usage: sge-status.py <jobid>", file=sys.stderr)
        # Still print a verdict: a usage error must not be read as job output.
        print("failed")
        return 0
    print(status(sys.argv[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
