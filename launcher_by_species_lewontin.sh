#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=evolsv-lewontin
#SBATCH --output=log/slurm-%A.out

species=$1

# module load snakemake/8.9.0
. /local/env/envconda.sh
source activate snakemake_v2


echo "Snakemake version"
snakemake --version
echo "SLURM version"
srun --version
echo "Conda version"
conda --version

git status

echo "Running Snakemake pipeline for species $species..."

SSL_CERT_FILE="ssl/cert.pem"
echo $SSL_CERT_FILE
# cat $SSL_CERT_FILE

snakemake -s workflow/Snakefile --configfile data-lewontin/config/config.yaml \
--use-conda --conda-frontend conda --profile ./profiles/slurm \
--cores 1 --rerun-incomplete \
--config samples="data-lewontin/config/samples_$species.tsv"