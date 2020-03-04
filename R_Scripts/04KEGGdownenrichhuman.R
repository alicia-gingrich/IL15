#install.packages("tidyverse")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
BiocManager::install("DESeq2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Cf.eg.db")

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

res9down <- res9%>%
  filter(padj <0.05) %>% 
  filter(log2FoldChange <0)

head(res9down)

res9downenrich <- enrichKEGG(gene = unique(res9down$GeneID), organism = "hsa")

res9downenrich_results <- res9downenrich@result

write.csv(as.data.frame(res9downenrich_results), 
          file="res9downenrich.csv")

dotplot(res9downenrich, showCategory=10)
upsetplot(res9downenrich)
cnetplot(res9downenrich)

#barplot(res9downenrich)
#emapplot(res9downenrich)