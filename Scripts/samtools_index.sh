#!/bin/bash

# Variables
sra=$1

# Indexage
samtools index $sra/${sra}_sorted.bam
