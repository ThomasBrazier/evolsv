#!/bin/bash
#SBATCH --mail-user=mail@mail.com
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=120GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=evolsv

conda activate snakemake

snakemake -s workflow/Snakefile --config config/config.yaml --use-conda --cores 8
