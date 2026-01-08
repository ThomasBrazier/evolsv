#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
wdir = args[1]
genome = args[2]


library(GenomicRanges)
library(ggplot2)
library(tidyverse)

# import a dataset
merged = read_tsv(paste0(wdir, "/", genome, "_final.tsv"))

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
fix.merged = merged[,1:6]
gt.merged = merged[,7:16]

fix.merged$sv_type = as.vector(gsub("SVTYPE=", "", str_match(fix.merged$INFO, "SVTYPE=[A-Z]+")))
fix.merged$sv_length = as.numeric(as.vector(gsub("SVLEN=", "", str_match(fix.merged$INFO, "SVLEN=[-]*[0-9]+"))))

fix.merged$supp_vec = as.vector(gsub("SUPP_VEC=", "", str_match(fix.merged$INFO, "SUPP_VEC=[0-1]+")))
fix.merged$n_supp = as.numeric(as.vector(gsub("SUPP=", "", str_match(fix.merged$INFO, "SUPP=[0-9]+"))))
fix.merged$jasmine = ifelse(!is.na(fix.merged$supp_vec), 1, 0)

# Get aligner support
fix.merged$minimap2 = ifelse(lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 1) == "1" |
                               lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 2) == "1" |
                               lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 3) == "1" |
                               lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 4) == "1", 1, 0)
fix.merged$ngmlr = ifelse(lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 5) == "1" |
                            lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 6) == "1" |
                            lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 7) == "1" |
                            lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 8) == "1", 1, 0)

# Get SV caller support
fix.merged$sniffles = ifelse(lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 1) == "1" | lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 5) == "1", 1, 0)
fix.merged$svim = ifelse(lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 2) == "1" | lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 6) == "1", 1, 0)
fix.merged$cutesv = ifelse(lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 3) == "1" | lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 7) == "1", 1, 0)
fix.merged$debreak = ifelse(lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 4) == "1" | lapply(strsplit(fix.merged$supp_vec, ""), `[[`, 8) == "1", 1, 0)


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

