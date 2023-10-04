#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 2) {
  stop("<vcf sniffles>, <outfile> and <ins_bed> must be passed in argument", call.=FALSE)
}
vcf_svim = args[1]
outfile = args[2]
ins_bed = args[3]


add_readnames_svim <- function(svim_vcf,ins_bed,outfile){ 

  library(tidyverse)
  library("glue")
  rn <-read.table(ins_bed)
  
  #Extract the read names from the ins_bed
  grep_rnames <- function(x){  
    k <- str_match_all(x, "./?(cigar|suppl)\\|\\s*(.*?)\\]")
    y <- glue_collapse(k[[1]][,3], sep = ",")
    return(y)
  }
  
  j <- sapply(rn$V6, grep_rnames)
  bind_cols(rn,j) -> rn
  colnames(rn) <- c("CHR","START","END","V4","V5","V6","RNAMES")
  
  
  #Split up vcf
  vcf <- as.data.frame(readLines(svim_vcf),stringsAsFactors = F)
  colnames(vcf) <- "V1"
  
  idx <- grepl("^#", vcf[,"V1"], ignore.case=TRUE)
  header <- vcf[idx,]
  data <- vcf[!idx,]
  
  rm(idx)
  rm(vcf)
  
  data <- as.data.frame(data)
  colnames(data) <- "V1"
  
  data1 <- data.frame(do.call('rbind', strsplit(as.character(data$V1),'\t',fixed=TRUE)))
  data2 <- data.frame(do.call('rbind', strsplit(as.character(data1$X8),';',fixed=TRUE)))
  
  data <- bind_cols(data1[,1:7],data2,data1[,9:10])
  colnames(data) <- c("CHR","START","a","b","c","d","e","SVTYPE","END","SVLEN","SUPPORT","STD_SPAN","STD_POS","seq","f","g")
  

  
  if(data$SVTYPE == "SVTYPE=INS"){
    data$c <- data.frame(do.call('rbind', strsplit(as.character(data$seq),'=',fixed=TRUE)))[,2]
  }
    
  
  data %>% unite("pos",c("CHR","START"), sep =".",remove = F) -> data
  rn %>% 
    #mutate(END2 = paste0("END=",as.character(END)))  %>% 
    unite("pos",c("CHR","START"), sep =".",remove = F) -> rn  
  
  merge(data,rn[,c(1,8)], by = "pos",all.x = T) -> data
  data$RNAMES %>% replace_na(".") -> data$RNAMES
  
  data %>% mutate(READNAMES2 = paste0("RNAMES=",as.character(RNAMES))) -> data
  
  #Remove text from info-field
  data %>% unite(info, SVTYPE,READNAMES2,END ,SVLEN,sep = ";") %>% select(CHR,START,a,b,c,d,e,info,f,g) -> vcf
  
  paste(vcf$CHR,vcf$START,vcf$a,vcf$b,vcf$c,vcf$d,vcf$e,vcf$info,vcf$f,vcf$g,sep = "\t" ) -> h
  
  
  as.data.frame(header) -> k
  as.data.frame(h) -> h
  
  colnames(k) <- "a"
  colnames(h) <- "a"
  
  dplyr::bind_rows(k,h) -> final
  
  write.table(final,outfile,quote = F,row.names = F, col.names = F)
}


add_readnames_svim(vcf_svim,ins_bed,outfile) 

