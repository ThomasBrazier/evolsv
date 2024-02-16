#!/bin/bash

# Sample
sra=$1
genome=$2
wdir=$3

# Download long reads
mkdir --parents $wdir
fastq-dump -v --gzip --outdir $wdir $sra

# Download reference genome assembly
datasets download genome accession $genome --filename $wdir/$sra.zip --include genome,seq-report
unzip $wdir/$sra.zip -d $wdir/
cp $wdir/ncbi_dataset/data/$genome/*_genomic.fna $wdir/$sra.fna
cp $wdir/ncbi_dataset/data/assembly_data_report.jsonl $wdir/${sra}_assembly_data_report.jsonl
cp $wdir/ncbi_dataset/data/$genome/sequence_report.jsonl $wdir/${sra}_sequence_report.jsonl
# bgzip $wdir/$sra/$sra.fna
