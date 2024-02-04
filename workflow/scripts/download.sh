#!/bin/bash

# Sample
sra=$1
genome=$2
wdir=$3

# Download long reads
mkdir --parents $wdir
fastq-dump -v --gzip --outdir $wdir $sra

# Download reference genome assembly
datasets download genome accession $genome --filename $wdir/$sra.zip
unzip $wdir/$sra.zip -d $wdir/
mv $wdir/ncbi_dataset/data/$genome/*_genomic.fna $wdir/$sra.fna
# bgzip $wdir/$sra/$sra.fna
