#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=120GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=evolsv-test
#SBATCH --output=log/slurm-%A.out

species=$1

module load snakemake

echo "Snakemake version"
snakemake --version
echo "SLURM version"
srun --version

echo "Running Snakemake pipeline for species $species..."
snakemake -s workflow/Snakefile --configfile config/config.yaml \
--use-conda --profile ./profiles/slurm --cores 1 --rerun-incomplete \
--config samples="data/config/samples_$species.tsv"