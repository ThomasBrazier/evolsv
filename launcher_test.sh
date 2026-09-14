#!/bin/bash
#
# Dry-run smoke test: checks that every supported input shape still builds a
# Snakemake DAG (snakemake -n), without running a single rule. Positive cases
# only -- for the deliberately malformed configs that must be rejected, see
# .test/test_dryrun.py (NEGATIVE_CASES).
#
# No --cores: it comes from the profile. Passing it here clamps every rule's
# threads to that value, even in cluster mode.
#SBATCH --mail-user=mail@mail.com
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=120GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=evolsv-test

source activate snakemake

cd "$(dirname "$0")" || exit 1

set -u
set -o pipefail

# Placeholder BAM/FASTQ/sample-sheet fixtures for the bam-mode, local-reference,
# single-run, no-scaffold-exclusion, bigtmp and ont cases. snakemake -n never reads
# file content, so a few bytes per file is enough to build the DAG.
bash .test/setup_fixtures.sh

# name:configfile pairs. An empty configfile means "no --configfile" (repo default,
# config/config.yaml). Keep this in sync with POSITIVE_CASES in .test/test_dryrun.py.
CASES=(
  "default:"
  "config-test:config/config_test.yaml"
  "single-run:.test/config_single_run.yaml"
  "bam-mode:.test/config_bam.yaml"
  "local-reference:.test/config_local_ref.yaml"
  "no-scaffold-exclusion:.test/config_no_scaffold_exclusion.yaml"
  "bigtmp:.test/config_bigtmp.yaml"
  "ont:.test/config_ont.yaml"
)

failed=()
for case in "${CASES[@]}"; do
  name="${case%%:*}"
  configfile="${case#*:}"

  args=(-s workflow/Snakefile -n -p --profile ./profiles/slurm)
  if [ -n "$configfile" ]; then
    args+=(--configfile "$configfile")
  fi

  echo "=== $name (${configfile:-repo default}) ==="
  if snakemake "${args[@]}"; then
    echo "PASS: $name"
  else
    echo "FAIL: $name"
    failed+=("$name")
  fi
  echo
done

if [ "${#failed[@]}" -ne 0 ]; then
  echo "Failed cases: ${failed[*]}"
  exit 1
fi

echo "All input cases built a DAG successfully."
