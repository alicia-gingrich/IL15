library(DESeq2)
#setwd ("C:/Users/Alicia/Desktop/il15")
counts <- read.csv("outputs/counts/dog_nonribo_counts.csv",row.names = 1)
counts <- apply(counts, 1:2, round)
samples <- read.csv("samples_dog_nonribo.csv")


dds<- DESeqDataSetFromMatrix(counts,
                             colData = samples,
                             design = ~ treatment)
dds <- DESeq(dds)

res1 <- results(dds, contrast=c("treatment","resting","il15"))
res2 <- results(dds, contrast=c("treatment","coculture","il15"))
res3 <- results(dds, contrast=c("treatment","resting","coculture"))

resOrdered1 <- res1[order(res1$pvalue),]
resOrdered2 <- res2[order(res2$pvalue),]
resOrdered3 <- res3[order(res3$pvalue),]

plotMA(res1, ylim=c(-10,10), alpha=0.05)
plotMA(res2, ylim=c(-15,15), alpha=0.05)
plotMA(res3, ylim=c(-10,10), alpha=0.05)

idx <- identify(res2$baseMean, res2$log2FoldChange)
rownames(res2)[idx]

plotCounts(dds, gene=which.min(res2$padj), intgroup="treatment")
plotCounts(dds, gene=which(rownames(res2)=="KLRF1"), intgroup="treatment")

write.csv(as.data.frame(resOrdered1), 
          file="outputs/DESeq/res2_treated_results.csv")

vsd <- vst(dds, blind=FALSE)

library(RColorBrewer)
pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(100)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("dog","treatment")])
heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
                    cluster_cols=TRUE, annotation_col=df, color = pal, legend = TRUE, annotation_colors = list(
                      dog = c(Dazey = "#FFFF00", Emma = "#FF0000", Iggy = "#FF00CC", Morgan = "#6600FF", Nick = "#9900CC", Santino = "#33FF00", Sessy= "#0000FF", HVE5 = "#FF9933", EKA5 = "#FF9933", IRS4 = "#FF6600", KPC6 = "#FF9900", VBE7 = "#FF6633"),
                      treatment = c(day0 = "#CCFFFF", day3 = "#99CCFF", day10 = "#6699FF", day17 = "#3366CC", day31 = "#003399", resting = "#00FFFF", il15 = "#00CCCC", coculture = "#009999")
                    ))

plotPCA(vsd, intgroup=c("dog", "treatment"))

library(ggplot2)

#customize PCA plot with ggplot function
pcaData <- plotPCA(vsd, intgroup=c("dog", "treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=dog, shape=treatment)) +
  geom_point(size=3) +
  scale_color_manual(values = c(Dazey = "#FFFF00", Emma = "#FF0000", Iggy = "#FF00CC", Morgan = "#6600FF", Nick = "#9900CC", Santino = "#33FF00", Sessy= "#0000FF", HVE5 = "#FF9933", EKA5 = "#FF9933", IRS4 = "#FF6600", KPC6 = "#FF9900", VBE7 = "#FF6633"))+
  scale_shape_manual(values = c(day0 = 3, day3 = 4, day10 = 7, day17 = 8, day31 = 10, resting = 15, il15 = 16, coculture = 17))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_light()
