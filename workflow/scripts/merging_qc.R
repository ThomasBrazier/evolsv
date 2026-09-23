#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
wdir = args[1]
genome = args[2]


library(GenomicRanges)
library(ggplot2)
library(tidyverse)

# import a dataset
merged = read_tsv(paste0(wdir, "/", genome, "_final.tsv"))

# If no SV in the final vcf/tsv, just output empty files
if (nrow(merged) == 0) {
  write_tsv(merged, paste0(wdir, "/merging_QC/", genome, "_svlen_equal_zero.tsv"))
  write_tsv(merged, paste0(wdir, "/merging_QC/", genome, "_avglen_equal_zero.tsv"))
  write_tsv(merged, paste0(wdir, "/merging_QC/", genome, "_no_avgend_field.tsv"))
  write_tsv(merged, paste0(wdir, "/merging_QC/", genome, "_unmerged_sv.tsv"))
} else {
  #--------------------------------------------------------
  # TODO
  # Add automatic tests
  # - TODO SVLEN=0?
  # - TODO END field?
  # - TODO Check if variants not merged, as for example with start = start + 1
  #--------------------------------------------------------



  #--------------------------------------------------------
  # Parse
  #--------------------------------------------------------
  # Column layout of <genome>_final.tsv, written by workflow/scripts/vcf_to_tsv.sh:
  # 6 fixed fields, FORMAT, one column per aligner + caller callset, then the
  # SVjedi-graph genotyped sample last. The number of callsets follows the `aligners`
  # config key, so it is read from the header rather than assumed (4 with one aligner,
  # 8 with both). finalQC.Rmd derives the same thing the same way.
  FORMAT_COL = 7
  tool_cols = names(merged)[(FORMAT_COL + 1):(ncol(merged) - 1)]
  tools = data.frame(
    column = tool_cols,
    # Callset columns are named <aligner>_<caller>; aligner names carry no underscore.
    aligner = sub("_.*$", "", tool_cols),
    caller = sub("^[^_]*_", "", tool_cols),
    stringsAsFactors = FALSE
  )

  fix.merged = merged[,1:(FORMAT_COL - 1)]
  gt.merged = merged[,FORMAT_COL:ncol(merged)]

  fix.merged$sv_type = as.vector(gsub("SVTYPE=", "", str_match(fix.merged$INFO, "SVTYPE=[A-Z]+")))
  fix.merged$sv_length = as.numeric(as.vector(gsub("SVLEN=", "", str_match(fix.merged$INFO, "SVLEN=[-]*[0-9]+"))))

  fix.merged$supp_vec = as.vector(gsub("SUPP_VEC=", "", str_match(fix.merged$INFO, "SUPP_VEC=[0-1]+")))
  fix.merged$n_supp = as.numeric(as.vector(gsub("SUPP=", "", str_match(fix.merged$INFO, "SUPP=[0-9]+"))))
  fix.merged$jasmine = ifelse(!is.na(fix.merged$supp_vec), 1, 0)

  # Bit i of Jasmine's SUPP_VEC is the callset tools$column[i]: both orders come from the
  # order the VCFs are listed in workflow/rules/merging.smk. Decode once into a logical
  # matrix instead of indexing the string by hand, which assumed 8 callsets in a fixed
  # aligner order and silently mislabelled them with any other `aligners` setting.
  supp_mat = matrix(FALSE, nrow = nrow(fix.merged), ncol = nrow(tools),
                    dimnames = list(NULL, tools$column))
  supp_ok = !is.na(fix.merged$supp_vec) & nchar(fix.merged$supp_vec) == nrow(tools)
  if (any(!is.na(fix.merged$supp_vec) & !supp_ok)) {
    warning("SUPP_VEC length does not match the ", nrow(tools),
            " callset columns of the final TSV; those calls are treated as unsupported.")
  }
  if (any(supp_ok)) {
    supp_mat[supp_ok, ] = do.call(
      rbind, strsplit(fix.merged$supp_vec[supp_ok], "")) == "1"
  }

  # Per-aligner and per-caller support: a call is supported by an aligner (or a caller)
  # when any of its callsets supports it. Not read by the checks below; part of the parsed
  # frame, and kept in step with finalQC.Rmd, which reports on the same columns.
  for (aligner in unique(tools$aligner)) {
    fix.merged[[aligner]] = as.integer(
      rowSums(supp_mat[, tools$aligner == aligner, drop = FALSE]) > 0)
  }
  for (caller in unique(tools$caller)) {
    fix.merged[[caller]] = as.integer(
      rowSums(supp_mat[, tools$caller == caller, drop = FALSE]) > 0)
  }


  #--------------------------------------------------------
  # Check if SVLEN=0 (error during parsing or preprocessing)
  #--------------------------------------------------------
  svlen_func = function(x) {
    gsub("SVLEN=", "", strsplit(merged$INFO[x], ";")[[1]][grepl("SVLEN",strsplit(merged$INFO[x], ";")[[1]])])
  }
  svlen = unlist(lapply(1:nrow(merged), svlen_func))
  svlen = as.numeric(svlen)

  # summary(svlen)

  # sum(svlen == 0)

  if (sum(svlen == 0) > 0) {print("SVLEN = 0")}

  # which(svlen == 0)

  # View(merged[which(svlen == 0),])

  write_tsv(merged[which(svlen == 0),], paste0(wdir, "/merging_QC/", genome, "_svlen_equal_zero.tsv"))


  #--------------------------------------------------------
  # Even if SVLEN=0, check if AVG_LEN is more reliable
  #--------------------------------------------------------
  avglen_func = function(x) {
    gsub("AVG_LEN=", "", strsplit(merged$INFO[x], ";")[[1]][grepl("AVG_LEN",strsplit(merged$INFO[x], ";")[[1]])])
  }
  avglen = unlist(lapply(1:nrow(merged), avglen_func))
  avglen = as.numeric(avglen)
  # summary(avglen)

  # sum(avglen == 0)

  if (sum(avglen == 0) > 0) {print("AVG_LEN = 0")}

  # which(avglen == 0)

  # View(merged[which(avglen == 0),])

  write_tsv(merged[which(avglen == 0),], paste0(wdir, "/merging_QC/", genome, "_avglen_equal_zero.tsv"))



  #--------------------------------------------------------
  # END field exists
  #--------------------------------------------------------
  # grepl("END=",strsplit(merged$INFO[x], ";")[[1]])

  end_exists_func = function(x) {
    grepl(";AVG_END=", merged$INFO[x])
  }
  end_exists = unlist(lapply(1:nrow(merged), end_exists_func))

  # sum(!end_exists)

  if (sum(!end_exists) > 0) {print("END field (AVG_END) does not exist")}

  # View(merged[which(!end_exists),])


  write_tsv(merged[which(!end_exists),], paste0(wdir, "/merging_QC/", genome, "_no_avgend_field.tsv"))


  #--------------------------------------------------------
  # Check that end and length are in agreement
  #--------------------------------------------------------

  # WORK IN PROGRESS
  # Inconsistencies between callers and between intial calls and jasmine merging
  merged$start = merged$POS

  end_func = function(x) {
    end = gsub("AVG_END=", "", strsplit(merged$INFO[x], ";")[[1]][grepl("AVG_END=",strsplit(merged$INFO[x], ";")[[1]], perl = TRUE)], perl = TRUE)
    round(as.numeric(end), digits = 0)
  }

  end = unlist(lapply(1:nrow(merged), end_func))
  merged$end = end

  merged$svlen = as.numeric(svlen)

  merged$width = merged$end - merged$start

  avg_len_func = function(x) {
    end = gsub("AVG_LEN=", "", strsplit(merged$INFO[x], ";")[[1]][grepl("AVG_LEN=",strsplit(merged$INFO[x], ";")[[1]], perl = TRUE)], perl = TRUE)
    round(as.numeric(end), digits = 0)
  }

  avg_len = unlist(lapply(1:nrow(merged), avg_len_func))
  merged$avg_len = avg_len


  # abs(merged$svlen) - abs(merged$width)
  # which(abs(merged$svlen) != abs(merged$width))

  #--------------------------------------------------------
  # Check if variants not merged, as for example with start = start + 1
  #--------------------------------------------------------
  svtype_func = function(x) {
    gsub("SVTYPE=", "", strsplit(merged$INFO[x], ";")[[1]][grepl("SVTYPE=",strsplit(merged$INFO[x], ";")[[1]], perl = TRUE)], perl = TRUE)
  }

  svtype = unlist(lapply(1:nrow(merged), svtype_func))
  merged$svtype = svtype


  # Find variants with less than 25 bp differences in positions
  sv_ranges = merged %>%
    select(chrom = CHROM,
          old_start = start,
          old_end = end,
          avg_len = avg_len)

  sv_ranges$start = ifelse(sv_ranges$old_start > sv_ranges$old_end, sv_ranges$old_end, sv_ranges$old_start)
  sv_ranges$end = sv_ranges$start + abs(sv_ranges$avg_len)

  sv_ranges = makeGRangesFromDataFrame(sv_ranges)

  hits = findOverlaps(sv_ranges, sv_ranges, type = "equal", maxgap = 50)

  hits
  hits = hits[queryHits(hits) != subjectHits(hits)]
  hits

  # keep only one over two - duplicates
  # hits = hits[seq(1, length(hits), by = 2)]

  # length(hits)
  # length(hits) / nrow(merged)

  if (length(hits) > 0) {print(paste0( length(hits), " variants (", round(length(hits) / nrow(merged), digits = 3) * 100, "%) not correctly merged."))}


  # sv_ranges[1764]
  # sv_ranges[1765]
  # 
  # View(merged[c(1764, 1765),])
  # idx = sort(c(queryHits(hits), subjectHits(hits)))
  idx = queryHits(hits)
  # idx
  # View(merged[idx,])
  unmerged = merged[idx,]
  unmerged$query = queryHits(hits)
  unmerged$subject = subjectHits(hits)

  write_tsv(unmerged, paste0(wdir, "/merging_QC/", genome, "_unmerged_sv.tsv"))

  if (length(hits) / nrow(merged) > 0.01) {print("More than 1% of variants were not merged successfully")}
}




