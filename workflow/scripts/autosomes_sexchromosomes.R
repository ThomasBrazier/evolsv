library(jsonlite)
library(dplyr)

seq = snakemake@params[['seq']]
sexchromosomes = snakemake@params[['sexchromosomes']]
autosomes = snakemake@params[['autosomes']]
scaffolds_to_exclude = snakemake@params[['scaffolds_to_exclude']]
chromosome_names = snakemake@params[['chromosome_names']]

# seq='data/GCA_905220545.2/GCA_905220545.2_sequence_report.jsonl'
# sexchromosomes='data/GCA_905220545.2/GCA_905220545.2.sexchromosomes'
# autosomes='data/GCA_905220545.2/GCA_905220545.2.autosomes'
# chromosome_names='data/GCA_905220545.2/GCA_905220545.2.chromosomes'
# scaffolds_to_exclude = "MT,Mt,Un"

cat(seq, "\n")

json_text = readLines(seq, warn = FALSE, encoding = "UTF-8")

df = lapply(json_text, fromJSON) |>
    bind_rows()

df = as.data.frame(df)

colnames(df)

df[,c("chrName", "genbankAccession", "length")]


df$start = 1
df$end = df$length


write.table(df[,c("chrName", "genbankAccession", "length")], chromosome_names, quote = F, row.names = F, col.names = T, sep = "\t")

formula = as.character(paste("--chr", as.vector(df$genbankAccession[which(df$chrName %in% c("W", "Z", "X", "Y") & df$role == "assembled-molecule")]), collapse = " "))
formula
formula = df[which(df$chrName %in% c("W", "Z", "X", "Y") & df$role == "assembled-molecule"), c("genbankAccession", "start", "end")]

write.table(formula, sexchromosomes, quote = F, row.names = F, col.names = F, sep = "\t")

exclude = c("W", "Z", "X", "Y", unlist(strsplit(scaffolds_to_exclude, ",")))

formula = as.character(paste("--chr", as.vector(df$genbankAccession[which(!(df$chrName %in% exclude) & df$role == "assembled-molecule")]), collapse = " "))
formula
formula = df[which(!(df$chrName %in% exclude) & df$role == "assembled-molecule"), c("genbankAccession", "start", "end")]

write.table(formula, autosomes, quote = F, row.names = F, col.names = F, sep = "\t")