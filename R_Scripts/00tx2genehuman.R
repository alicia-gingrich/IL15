if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")


library(tximport)
library(readr)
library(dplyr)
setwd ("C:/Users/Alicia/Desktop/il15")

#find and download tabular file for human from NCBI(feature_table.txt.gz)
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_feature_table.txt.gz",
              destfile = "inputs/GCF_000001405.39_GRCh38.p12_feature_table.txt.gz")
feat_table_hs <- read_tsv('inputs/GCF_000001405.39_GRCh38.p12_feature_table.txt.gz')

#select relevant columns: product accession and GeneID (in that order)
feat_table_hs <- select(feat_table_hs, product_accession, symbol, GeneID)
feat_table_hs <- feat_table_hs%>%
  unique()

#make sure table is capturing all salmon counts, should be 158690 obs in human
quant <- read_tsv("outputs/quant_hs_nonribo/Human_NK_Day_0_122118_12_21_18_quant.sf")

table(quant$Name %in% feat_table_hs$product_accession)
#should have 158689 = TRUE
write.table(feat_table_hs, "outputs/tx2gene/humantx2gene.tsv", quote = F, row.names = F, sep = "\t")
