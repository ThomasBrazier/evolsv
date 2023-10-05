#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 1) {
  stop("<accession> must be passed in argument", call.=FALSE)
}
accession = args[1]
# vcf_svim = args[1]
# outfile = args[3]
# ins_bed = args[2]

vcf_svim = svim = paste0(accession, "/", accession, "_svim.vcf")
ins_bed = paste0(accession, "/candidates/candidates_novel_insertions.bed")
outfile = paste0(accession, "/", accession, "_svim_refined.vcf")
  
add_readnames_svim <- function(vcf_svim,ins_bed,outfile){ 

  library(tidyverse)
  library("glue")
  rn = read.table(ins_bed)
  
  #Extract the read names from the ins_bed
  readnames = function(x) {
    m = gregexpr(paste0(accession, ".[0-9]+"), x)
    j = unlist(regmatches(x, m))
    j = paste(as.character(j), collapse = ",")
    return(j)
  }
  
  # x = "[CAJMZP020000002.1|409|501|INS;cigar|ERR6608653.438229][CAJMZP020000002.1|580|685|INS;cigar|ERR6608653.352096]"
  # gregexpr(paste0(accession, ".[0-9]+"), x)
  
  l = lapply(as.list(rn$V7), function(x){readnames(x)})
  

  
  #grep_rnames <- function(x){  
    #k <- str_match_all(x, "./?(cigar|suppl)\\|\\s*(.*?)\\]")
    #y <- glue_collapse(k[[1]][,3], sep = ",")
    #return(y)
  #}
  
  #j <- sapply(rn$V6, grep_rnames)
  bind_cols(rn,l) -> rn
  colnames(rn) <- c("CHR","START","END","V4","V5","V6", "V8","RNAMES")
  
  
  #Split up vcf
  vcf <- as.data.frame(readLines(vcf_svim),stringsAsFactors = F)
  colnames(vcf) <- "V1"
  
  idx <- grepl("^#", vcf[,"V1"], ignore.case=TRUE)
  header <- vcf[idx,]
  data <- vcf[!idx,]
  
  rm(idx)
  rm(vcf)
  
  data <- as.data.frame(data)
  colnames(data) <- "V1"
  
  data1 <- data.frame(do.call('rbind', strsplit(as.character(data$V1),'\t',fixed=TRUE)))
  #data2 <- data.frame(do.call('bind_cols', strsplit(as.character(data1$X8),';',fixed=TRUE)))
  svtype = regexpr(paste0("SVTYPE=[A-Z]+"), data1$X8)
  svtype = regmatches(data1$X8, svtype)
  
  data1$X8[which(svtype == "SVTYPE=INS")] = paste0(data1$X8[which(svtype == "SVTYPE=INS")], ";RNAMES=", rn$RNAMES)
  
  data1$X8[which(svtype != "SVTYPE=INS")] = paste0(data1$X8[which(svtype != "SVTYPE=INS")], ";RNAMES=NULL")
  

  #endpos = regexpr(paste0("END=[0-9]+"), data1$X8)
  #endpos = ifelse(endpos == -1, NA ,regmatches(data1$X8, endpos))
  #svlen = regexpr(paste0("SVLEN=[0-9]+"), data1$X8)
  #svlen = ifelse(svlen == -1, NA ,regmatches(data1$X8, svlen))
  #support = regexpr(paste0("SUPPORT=[0-1]+"), data1$X8)
  #support = ifelse(support == -1, NA ,regmatches(data1$X8, support))
  #stdspan = regexpr(paste0("STD_SPAN=[A-Z0-9.]+"), data1$X8)
  #stdspan = ifelse(stdspan == -1, NA ,regmatches(data1$X8, stdspan))
  #stdpos = regexpr(paste0("STD_POS=[A-Z0-9.]+"), data1$X8)
  #stdpos = ifelse(stdpos == -1, NA ,regmatches(data1$X8, stdpos))
  #seqs = regexpr(paste0("SEQS=[A-Z]+"), data1$X8)
  #seqs = ifelse(seqs == -1, NA ,regmatches(data1$X8, seqs))


  #data <- bind_cols(data1[,1:7],data2,data1[,9:10])
  #data <- bind_cols(data1[,1:7],svtype, endpos, svlen, support, stdspan, stdpos, seqs,data1[,9:10])
  #colnames(data) <- c("CHR","START","a","b","c","d","e","SVTYPE","END","SVLEN","SUPPORT","STD_SPAN","STD_POS","seq","f","g")
  


  #if(data$SVTYPE == "SVTYPE=INS"){
  #  data$c <- data.frame(do.call('rbind', strsplit(as.character(data$seq),'=',fixed=TRUE)))[,2]
  #}
    
  
  # data %>% unite("pos",c("CHR","START"), sep =".",remove = F) -> data
  # rn %>% 
  #   #mutate(END2 = paste0("END=",as.character(END)))  %>% 
  #   unite("pos",c("CHR","START"), sep =".",remove = F) -> rn  
  # 
  # merge(data,rn[,c(1,8)], by = "pos",all.x = T) -> data
  # data$RNAMES %>% replace_na(".") -> data$RNAMES
  # 
  # data %>% mutate(READNAMES2 = paste0("RNAMES=",as.character(RNAMES))) -> data
  # 
  # #Remove text from info-field
  # data %>% unite(info, SVTYPE,READNAMES2,END ,SVLEN,sep = ";") %>% select(CHR,START,a,b,c,d,e,info,f,g) -> vcf
  # 
  # paste(vcf$CHR,vcf$START,vcf$a,vcf$b,vcf$c,vcf$d,vcf$e,vcf$info,vcf$f,vcf$g,sep = "\t" ) -> h
  # 
  # 
  as.data.frame(header) -> k
  as.data.frame(h) -> h

  colnames(k) <- "a"
  colnames(h) <- "a"

  dplyr::bind_rows(k,h) -> final
  
  write.table(final,outfile,quote = F,row.names = F, col.names = F)
}


add_readnames_svim(vcf_svim, ins_bed, outfile) 

