#!/bin/bash
#SBATCH --mail-user=user@mail.com
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --time=20-60:00:00
#SBATCH --job-name=evolsv
#SBATCH --output=log/slurm-%A.out

species=$1

# module load snakemake/8.9.0
# module load snakemake/8.27.1
# module load snakemake/9.4.0
# . /local/env/envconda.sh
# source activate snakemake
module load conda
source activate snakemake_v2


echo "Snakemake version"
snakemake --version
echo "SLURM version"
srun --version
echo "Conda version"
conda --version

git status

echo "Running Snakemake pipeline for species $species..."

snakemake -s workflow/Snakefile --configfile data-lewontin/config/config.yaml \
--use-conda --conda-frontend conda --profile ./profiles/slurm \
--cores 1 --rerun-incomplete --printshellcmds --unlock \
--config samples="data-lewontin/config/samples_$species.tsv"

snakemake -s workflow/Snakefile --configfile data-lewontin/config/config.yaml \
--use-conda --conda-frontend conda --profile ./profiles/slurm \
--cores 1 --rerun-incomplete --printshellcmds --slurm-status-attempts 100000 \
--config samples="data-lewontin/config/samples_$species.tsv"