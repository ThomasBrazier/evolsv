#!/bin/bash

# Variables
sra=$1

# Contrôle qualité des reads
NanoPlot --fastq $sra/$sra.fastq.gz --N50 --verbose
mv NanoStats.txt $sra/${sra}_nanoplotstats.txt
