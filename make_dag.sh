#! /bin/bash

module load snakemake
module load graphviz

snakemake -s workflow/Snakefile --cores 1 --dag | dot -Tpdf > dag.pdf
