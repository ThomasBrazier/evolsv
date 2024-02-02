#!/bin/bash

# Variables
sra=$1

# Contrôle qualité du mapping
samtools stats $sra/${sra}_sorted.bam > $sra/${sra}_mapping.stats
