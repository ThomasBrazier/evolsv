#!/usr/bin/env bash
#
# bsub wrapper for the cluster-generic executor on IBM Spectrum LSF.
#
# Invoked by snakemake as
#     lsf-submit.sh <name> <threads> <mem_mb> <runtime_min> <queue> <project> \
#                   <mem_unit> <mem_scope> <jobscript>
# and must print the job id, and nothing else, on stdout: cluster-generic reads the
# first line of stdout verbatim as the id. bsub prints "Job <123> is submitted to
# queue <normal>.", so the id is extracted from it.
set -euo pipefail

if [ "$#" -ne 9 ]; then
    echo "usage: $(basename "$0") <name> <threads> <mem_mb> <runtime_min> <queue> <project> <mem_unit> <mem_scope> <jobscript>" >&2
    exit 1
fi

name=$1
threads=$2
mem_mb=$3
runtime=$4
queue=$5
project=$6
mem_unit=$7
mem_scope=$8
jobscript=$9

# Memory for the request, rounded up so the job never gets less than mem_mb in total.
case "$mem_scope" in
job) mem=$mem_mb ;;
slot) mem=$(((mem_mb + threads - 1) / threads)) ;;
*)
    echo "lsf_mem_scope must be 'job' or 'slot', got '$mem_scope'" >&2
    exit 1
    ;;
esac

# Convert MB to the cluster's LSF_UNIT_FOR_LIMITS, rounding up.
case "$mem_unit" in
KB) mem=$((mem * 1024)) ;;
MB) ;;
GB) mem=$(((mem + 1023) / 1024)) ;;
*)
    echo "lsf_mem_unit must be KB, MB or GB, got '$mem_unit'" >&2
    exit 1
    ;;
esac
if [ "$mem" -lt 1 ]; then
    mem=1
fi

# -W takes [hours:]minutes, so plain minutes must be split.
walltime=$(printf '%d:%02d' $((runtime / 60)) $((runtime % 60)))

args=(
    -J "$name"
    -n "$threads"
    -q "$queue"
    -W "$walltime"
    -M "$mem"
    # span[hosts=1]: all slots on one node (threads, not MPI ranks).
    -R "span[hosts=1] rusage[mem=${mem}]"
    -o "logs/lsf/%J.out"
    -e "logs/lsf/%J.err"
)
if [ -n "$project" ]; then
    args+=(-P "$project")
fi

# bsub reads the job script on stdin (script-file arguments are not run as scripts).
out=$(bsub "${args[@]}" <"$jobscript")
if [[ $out =~ Job\ \<([0-9]+)\> ]]; then
    echo "${BASH_REMATCH[1]}"
else
    echo "could not parse job id from bsub output: $out" >&2
    exit 1
fi
