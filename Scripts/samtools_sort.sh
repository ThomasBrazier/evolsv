#!/bin/bash

# Variables
sra=$1

# Tri
samtools sort $sra/$sra.bam -o $sra/${sra}_sorted.bam
