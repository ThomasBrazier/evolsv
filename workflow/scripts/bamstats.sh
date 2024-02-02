#!/bin/bash

# Variables
sra=$1

# Visualisation du contrôle qualité du mapping
plot-bamstats -p $sra/${sra}_mapping_visu $sra/${sra}_mapping.stats
