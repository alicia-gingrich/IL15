
setwd ("C:/Users/Alicia/Desktop/il15")
install.packages("tidyverse")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
BiocManager::install("DESeq2")

library(readr)
library(dplyr)
library(clusterProfiler)
library(enrichplot)

res9 <- read.csv("C:/Users/Alicia/Desktop/il15/res9_treated_results.csv", stringsAsFactors = F)
feat_table_human <-read_tsv("C:/Users/Alicia/Desktop/il15/outputs/tx2gene/humantx2gene.tsv")

feat_table_human <- feat_table_human%>%
  select(symbol, GeneID)%>%
  distinct()


head(res9)
head(feat_table_human)

res9 <- left_join(res9, feat_table_human, by= c("X" = "symbol"))

res9up <- res9%>%
  filter(padj <0.05) %>% 
  filter(log2FoldChange >0)

head(res9up)

res9upenrich <- enrichKEGG(gene = unique(res9up$GeneID), organism = "hsa")

res9upenrich_results <- res9upenrich@result

write.csv(as.data.frame(res9upenrich_results), 
          file="res9upenrich.csv")

dotplot(res9upenrich, showCategory=10)
upsetplot(res9upenrich)
cnetplot(res9upenrich)

#barplot(res9upenrich)
#emapplot(res9upenrich)

