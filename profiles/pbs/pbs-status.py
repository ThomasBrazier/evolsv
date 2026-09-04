#!/usr/bin/env python3
"""Job-status probe for the cluster-generic executor on PBS Pro / OpenPBS / Torque.

Snakemake invokes this as `pbs-status.py <jobid>` and requires exactly one of
`running`, `success` or `failed` on stdout. Anything else, or a non-zero exit, makes
cluster-generic treat the job as failed -- so every path here prints one of the three
and exits 0.

Only the standard library is used: this runs under whichever interpreter is on the
submitting host's PATH, which is not necessarily the pipeline's conda environment.

The hard case is a job that PBS no longer knows about. `qstat -x` reports finished
jobs, but only for as long as the site's `job_history_duration` keeps them, and Torque
drops jobs from qstat immediately. Reporting `running` for an unknown job would hang
the DAG forever, so an unknown id is reported as `success` once the grace window below
has passed: Snakemake then checks the rule's output files itself and raises a clear
MissingOutputException if the job really did die. That is a better failure mode than a
silent hang, but it does mean a job killed by the memory enforcer after writing its
outputs can be recorded as successful.
"""

import subprocess
import sys
import time

# PBS job states that mean "not finished yet". F (finished) and C (completed, Torque)
# fall through to the exit-status check.
RUNNING_STATES = frozenset("QHRETWSB")

# An id can be unknown to qstat for a few seconds right after submission, before the
# server has registered it. Poll a few times before concluding the job is gone.
ATTEMPTS = 3
BACKOFF_SECONDS = 2


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
    """Pull `field = value` out of qstat -f output.

    qstat folds long lines by breaking them with a newline plus a leading tab, so the
    text is unfolded first.
    """
    unfolded = text.replace("\n\t", "")
    for line in unfolded.splitlines():
        key, separator, value = line.partition(" = ")
        if separator and key.strip() == field:
            return value.strip()
    return None


def classify(text):
    """Map one qstat -f record to running/success/failed, or None if undecidable."""
    state = parse_field(text, "job_state")
    if state is None:
        return None
    if state in RUNNING_STATES:
        return "running"
    exit_status = parse_field(text, "Exit_status")
    if exit_status is None:
        # Finished but the exit code has not been recorded yet.
        return "running"
    return "success" if exit_status == "0" else "failed"


def status(jobid):
    for attempt in range(ATTEMPTS):
        # -x includes finished jobs on PBS Pro. Torque has no -x, so the plain form is
        # tried as well; on Torque a finished job is simply absent from both.
        for command in (["qstat", "-f", "-x", jobid], ["qstat", "-f", jobid]):
            text = run(command)
            if text:
                verdict = classify(text)
                if verdict is not None:
                    return verdict
        if attempt < ATTEMPTS - 1:
            time.sleep(BACKOFF_SECONDS)

    # Torque only: tracejob can still see a job that has left qstat.
    text = run(["tracejob", "-n", "2", jobid.split(".")[0]])
    if text:
        for line in text.splitlines():
            if "Exit_status=" in line:
                code = line.split("Exit_status=")[1].split()[0]
                return "success" if code == "0" else "failed"

    # Unknown to PBS. See the module docstring for why this is not "running".
    return "success"


def main():
    if len(sys.argv) != 2:
        print("usage: pbs-status.py <jobid>", file=sys.stderr)
        # Still print a verdict: a usage error must not be read as job output.
        print("failed")
        return 0
    print(status(sys.argv[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
