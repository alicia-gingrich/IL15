library(DESeq2)
setwd ("C:/Users/Alicia/Desktop/il15")
counts_hs <- read.csv("outputs/counts/human_nonribo_counts.csv",row.names = 1)
counts_hs <- apply(counts_hs, 1:2, round)
samples_hs <- read.csv("samples_human_nonribo.csv")


dds_hs<- DESeqDataSetFromMatrix(counts_hs,
                                colData = samples_hs,
                                design = ~ treatment)
dds_hs <- DESeq(dds_hs)

res7 <- results(dds_hs, contrast=c("treatment","resting","il15"))
res8 <- results(dds_hs, contrast=c("treatment","coculture","il15"))
res9 <- results(dds_hs, contrast=c("treatment","resting","coculture"))

#We can order our results table by the smallest p value:
resOrdered7 <- res7[order(res7$pvalue),]
resOrdered8 <- res8[order(res8$pvalue),]
resOrdered9 <- res9[order(res9$pvalue),]

plotMA(res7, ylim=c(-10,10))
plotMA(res8, ylim=c(-10,10))
plotMA(res9, ylim=c(-10,10))

idx <- identify(res1$baseMean, res1$log2FoldChange)
rownames(res1)[idx]

plotCounts(dds_hs, gene=which.min(res7$padj), intgroup="treatment")
plotCounts(dds_hs, gene=which(rownames(res7)=="CD247"), intgroup="treatment")

write.csv(as.data.frame(resOrdered7), 
          file="outputs/DESeq/res7_treated_results.csv")


vsd_hs <- vst(dds_hs, blind=FALSE)

library("pheatmap")
select_hs <- order(rowMeans(counts(dds_hs,normalized=TRUE)),
                   decreasing=TRUE)[1:20]
df_hs <- as.data.frame(colData(dds_hs)[,c("treatment","human")])
pheatmap(assay(vsd_hs)[select_hs,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df_hs)

plotPCA(vsd_hs, intgroup=c("treatment", "human"))
