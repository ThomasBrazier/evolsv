#!/bin/bash
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

snakemake -s workflow/Snakefile --configfile config/config_test.yaml --use-conda --profile ./profiles/slurm
