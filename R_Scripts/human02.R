library(DESeq2)
setwd ("C:/Users/Alicia/Desktop/il15")
counts_hs <- read.csv("outputs/counts/human_nonribo_counts.csv",row.names = 1)
counts_hs <- apply(counts_hs, 1:2, round)
samples_hs <- read.csv("samples_human_nonribo.csv")


dds_hs<- DESeqDataSetFromMatrix(counts_hs,
                                colData = samples_hs,
                                design = ~ treatment)
dds_hs <- DESeq(dds_hs)

res4 <- results(dds_hs, contrast=c("treatment","resting","il15"))
res5 <- results(dds_hs, contrast=c("treatment","coculture","il15"))
res6 <- results(dds_hs, contrast=c("treatment","resting","coculture"))

#We can order our results table by the smallest p value:
resOrdered4 <- res4[order(res4$pvalue),]
resOrdered5 <- res5[order(res5$pvalue),]
resOrdered6 <- res6[order(res6$pvalue),]

plotMA(res4, ylim=c(-10,10))
plotMA(res5, ylim=c(-10,10))
plotMA(res6, ylim=c(-10,10))

idx <- identify(res1$baseMean, res1$log2FoldChange)
rownames(res1)[idx]

plotCounts(dds_hs, gene=which.min(res4$padj), intgroup="treatment")
plotCounts(dds_hs, gene=which(rownames(res4)=="TIGIT"), intgroup="treatment")

write.csv(as.data.frame(resOrdered4), 
          file="res4_treated_results.csv")

vsd_hs <- vst(dds_hs, blind=FALSE)

library("pheatmap")
select_hs <- order(rowMeans(counts(dds_hs,normalized=TRUE)),
                   decreasing=TRUE)[1:100]
df_hs <- as.data.frame(colData(dds_hs)[,c("treatment","human")])
pheatmap(assay(vsd_hs)[select_hs,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df_hs)

plotPCA(vsd_hs, intgroup=c("treatment", "human"))
