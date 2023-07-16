#!/bin/bash

snakemake -s pipeline.snake --use-conda --conda-frontend "mamba" --cores 8 --rerun-incomplete --stats statistiques.json
