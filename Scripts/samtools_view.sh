#!/bin/bash

# Variables
sra=$1

# Conversion
samtools view -S -b $sra/$sra.sam > $sra/$sra.bam
