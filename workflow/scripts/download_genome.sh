#!/bin/bash

# Sample
genome=$1
wdir=$2

# Download reference genome assembly
datasets download genome accession $genome --filename $wdir/$sra.zip --include genome,seq-report
unzip $wdir/$sra.zip -d $wdir/
cp $wdir/ncbi_dataset/data/$genome/*_genomic.fna $wdir/$sra.fna
cp $wdir/ncbi_dataset/data/assembly_data_report.jsonl $wdir/${sra}_assembly_data_report.jsonl
cp $wdir/ncbi_dataset/data/$genome/sequence_report.jsonl $wdir/${sra}_sequence_report.jsonl
# bgzip $wdir/$sra/$sra.fna
