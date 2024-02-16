#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
wdir = args[1]
sra = args[2]

library(rmarkdown)
library(vcfR)
library(adegenet)
library(poppr)

rmarkdown::render('workflow/scripts/finalQC.Rmd',
                    output_file = paste0(sra, ".finalQC.html"),
                    output_dir = paste0(wdir),
                    params = list(wdir = wdir, sra = sra))
