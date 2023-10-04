#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 2) {
  stop("<vcf sniffles> and <outfile> must be passed in argument", call.=FALSE)
}
vcf_sniffles = args[1]
outfile = args[2]



reformat_sniffles <- function(vcf_sniffles,outfile){  
  library(tidyverse)
  
  vcf <- as.data.frame(readLines(vcf_sniffles),stringsAsFactors = F)
  colnames(vcf) <- "V1"
  
  idx <- grepl("^#", vcf[,"V1"], ignore.case=TRUE)
  header <- vcf[idx,]
  data <- vcf[!idx,]
  
  rm(idx)
  rm(vcf)
  
  #Split data
  data <- as.data.frame(data)
  colnames(data) <- "V1"
  
  data1 <- data.frame(do.call('rbind', strsplit(as.character(data$V1),'\t',fixed=TRUE)))
  data2 <- data.frame(do.call('rbind', strsplit(as.character(data1$X8),';',fixed=TRUE)))
  
  data <- bind_cols(data1[,1:7],data2,data1[,9:10])
  rm(data1)
  rm(data2)
  
  #Remove text from info-field
  data %>% unite(info, X9...16,X10...17,X4...11 ,X12,X15,sep = ";") %>% select(X1...1,
                                                                     X2...2,X3...3,X4...4,X5...5,X6...6,
                                                                     X7...7,info,X9...26,X10...27) -> data_new
  

 paste(data_new$X1...1,data_new$X2...2,data_new$X3...3,data_new$X4...4,data_new$X5...5,data_new$X6...6,
 data_new$X7...7,data_new$info,data_new$X9...26,data_new$X10...27,sep = "\t" ) -> h


  as.data.frame(header) -> k
  as.data.frame(h) -> h

  colnames(k) <- "a"
  colnames(h) <- "a"

  dplyr::bind_rows(k,h) -> final

  write.table(final,outfile,quote = F,row.names = F, col.names = F)
}


reformat_sniffles(vcf_sniffles,outfile)

