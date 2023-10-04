#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=16
#SBATCH --mem=120GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=evolsv

. /local/env/envsnakemake-7.28.3.sh

snakemake -s pipeline.snake --use-conda --conda-frontend "mamba" --cores 8 --rerun-incomplete --stats statistiques.json
