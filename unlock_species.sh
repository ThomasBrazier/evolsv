#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=1-60:00:00
#SBATCH --job-name=evolsv-unlock
#SBATCH --output=log/slurm-%A.out

species=$1

module load snakemake/8.9.0

echo "Snakemake version"
snakemake --version
echo "SLURM version"
srun --version
echo "Conda version"
conda --version

git status

echo "Unlocking Snakemake pipeline for species $species..."

snakemake -s workflow/Snakefile --configfile config/config.yaml \
--use-conda --conda-frontend conda --unlock --profile ./profiles/slurm --cores 1 --rerun-incomplete \
--config samples="data/config/samples_$species.tsv"