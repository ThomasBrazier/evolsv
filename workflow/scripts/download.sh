#!/bin/bash

# Sample
sra=$1
genome=$2
wdir=$3

# Download long reads
mkdir --parents $wdir/$sra
fastq-dump -v --gzip --outdir $wdir/$sra $sra

# Download reference genome assembly
datasets download genome accession $genome --filename $wdir/$sra/$sra.zip
unzip $wdir/$sra/$sra.zip -d $wdir/$sra/
mv $wdir/$sra/ncbi_dataset/data/$genome/*_genomic.fna $wdir/$sra/$sra.fna
# bgzip $wdir/$sra/$sra.fna
