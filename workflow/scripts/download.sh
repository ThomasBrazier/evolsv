#!/bin/bash

# Variables
sra=$1
genome=$2

# Téléchargement des reads
mkdir --parents $sra
fastq-dump -v --gzip --outdir $sra $sra

# Téléchargement du génome de référence
datasets download genome accession $genome --filename $sra/ref_$sra.zip
unzip $sra/ref_$sra.zip -d $sra/
mv $sra/ncbi_dataset/data/$genome/*_genomic.fna $sra/ref_$sra.fna
bgzip $sra/ref_$sra.fna
