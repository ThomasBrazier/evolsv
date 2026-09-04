#!/usr/bin/env bash
#
# qsub wrapper for the cluster-generic executor on SGE / SoGE / UGE.
#
# Invoked by snakemake as
#     sge-submit.sh <name> <threads> <mem_mb> <runtime_min> <queue> <pe> <jobscript>
# and must print the job id, and nothing else, on stdout: cluster-generic reads the
# first line of stdout verbatim as the id. Hence `qsub -terse`.
#
# The wrapper exists because SGE enforces h_vmem PER SLOT while the profile declares
# mem_mb as a total for the job. Snakemake resource expressions cannot reference each
# other, so the division cannot live in profiles/sge/config.yaml.
set -euo pipefail

if [ "$#" -ne 7 ]; then
    echo "usage: $(basename "$0") <name> <threads> <mem_mb> <runtime_min> <queue> <pe> <jobscript>" >&2
    exit 1
fi

name=$1
threads=$2
mem_mb=$3
runtime=$4
queue=$5
pe=$6
jobscript=$7

# Per-slot memory, rounded up so the job never gets less than mem_mb in total.
mem_per_slot=$(((mem_mb + threads - 1) / threads))
if [ "$mem_per_slot" -lt 1 ]; then
    mem_per_slot=1
fi

# h_rt takes [[hours:]minutes:]seconds and accepts a minutes field above 59, so
# "<minutes>:00" passes the profile's runtime through unchanged.
args=(
    -terse
    -N "$name"
    -cwd
    -V
    -j y
    -o logs/sge/
    -q "$queue"
    -l "h_vmem=${mem_per_slot}M"
    -l "h_rt=${runtime}:00"
)

# A parallel environment must not be requested for single-slot jobs: many sites
# configure PEs with an allocation rule that rejects a request of 1.
if [ "$threads" -gt 1 ]; then
    args+=(-pe "$pe" "$threads")
fi

exec qsub "${args[@]}" "$jobscript"
