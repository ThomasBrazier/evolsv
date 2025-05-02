#!/bin/bash

# Sample
sra=$1
wdir=$2

# Download long reads
mkdir --parents $wdir
fastq-dump -v --gzip --outdir $wdir $sra