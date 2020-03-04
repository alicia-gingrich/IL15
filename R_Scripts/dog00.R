library(tximport)
library(readr)
library(dplyr)
setwd ("C:/Users/Alicia/Desktop/il15")

#find and download tabular file for canine from NCBI (feature_table.txt.gz)
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_feature_table.txt.gz",
              destfile = "inputs/GCF_000002285.3_CanFam3.1_feature_table.txt.gz")
feat_table <- read_tsv('inputs/GCF_000002285.3_CanFam3.1_feature_table.txt.gz')

#select relevant columns: product accession and GeneID (in that order)
feat_table <- select(feat_table, product_accession, symbol, GeneID)
feat_table <- feat_table%>%
  unique()

#make sure table is capturing all salmon counts, should be 82430 obs in canine
quant <- read_tsv("outputs/quant_cf_nonribo/Canine_NK_Day_0_HVE5_2_15_19_quant.sf")

table(quant$Name %in% feat_table$product_accession)
write.table(feat_table, "outputs/tx2gene/dogtx2gene.tsv", quote = F, row.names = F, sep = "\t")
