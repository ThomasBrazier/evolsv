#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
wdir = args[1]
genome = args[2]
sample = args[3]

library(rmarkdown)
library(vcfR)
library(adegenet)
library(poppr)

rmarkdown::render('workflow/scripts/finalQC.Rmd',
                    output_file = paste0(genome, "_finalQC.html"),
                    output_dir = paste0(wdir, "/", sample),
                    params = list(wdir = wdir, genome = genome, sample = sample))


# rmarkdown::render('workflow/scripts/finalQC.Rmd',
#                     output_file = paste0(genome, "_finalQC.pdf"),
#                     output_dir = paste0(wdir, "/", sample),
#                     params = list(wdir = wdir, genome = genome, sample = sample))
