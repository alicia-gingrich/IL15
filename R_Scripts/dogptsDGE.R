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
quant <- read_tsv("outputs/quant_cf_nonribo/Dazey_051719_dim_quant.sf")

table(quant$Name %in% feat_table$product_accession)
write.table(feat_table, "outputs/tx2gene/dogtx2gene.tsv", quote = F, row.names = F, sep = "\t")

# read in file names 
files <- read_csv("samples_dog_pts_only.csv")

# read in dog counts
tx2gene <- read_tsv("outputs/tx2gene/dogtx2gene.tsv")
canine <- tximport(files = files$files, type = "salmon", tx2gene = tx2gene)
caninecounts <- canine$counts

# change rownames for provenance
colnames(caninecounts) <- files$sample

# write to file
write.csv(caninecounts, "outputs/counts/dog_pts_counts.csv", quote = F, row.names = T)

library(DESeq2)
counts <- read.csv("outputs/counts/dog_pts_counts.csv",row.names = 1)
counts <- apply(counts, 1:2, round)
samples <- read.csv("samples_dog_pts_only.csv")


dds<- DESeqDataSetFromMatrix(counts,
                             colData = samples,
                             design = ~ treatment)
dds <- DESeq(dds)

#res1 <- results(dds, contrast=c("treatment","resting","il15"))
#res2 <- results(dds, contrast=c("treatment","coculture","il15"))
#res3 <- results(dds, contrast=c("treatment","resting","coculture"))

#resOrdered1 <- res1[order(res1$pvalue),]
#resOrdered2 <- res2[order(res2$pvalue),]
#resOrdered3 <- res3[order(res3$pvalue),]

#plotMA(res1, ylim=c(-10,10), alpha=0.05)
#plotMA(res2, ylim=c(-15,15), alpha=0.05)
#plotMA(res3, ylim=c(-10,10), alpha=0.05)

#idx <- identify(res2$baseMean, res2$log2FoldChange)
#rownames(res2)[idx]

#plotCounts(dds, gene=which.min(res2$padj), intgroup="treatment")
#plotCounts(dds, gene=which(rownames(res2)=="KLRF1"), intgroup="treatment")

#write.csv(as.data.frame(resOrdered1), 
#file="res2_treated_results.csv")

vsd <- vst(dds, blind=FALSE)

library(RColorBrewer)
pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(100)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("dog","treatment")])
heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, color = pal, legend = TRUE, annotation_colors = list(
           dog = c(Dazey = "#FFFF00", Emma = "#FF0000", Iggy = "#FF00CC", Morgan = "#6600FF", Nick = "#9900CC", Santino = "#33FF00", Sessy= "#0000FF"),
           treatment = c(day0 = "#CCFFFF", day3 = "#99CCFF", day10 = "#6699FF", day17 = "#3366CC", day31 = "#003399")
         ))

library(ggplot2)

plotPCA(vsd, intgroup=c("dog", "treatment"))

pcaData <- plotPCA(vsd, intgroup=c("dog", "treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=dog, shape=treatment)) +
  geom_point(size=3) +
  scale_color_manual(values = c(Dazey = "#FFFF00", Emma = "#FF0000", Iggy = "#FF00CC", Morgan = "#6600FF", Nick = "#9900CC", Santino = "#33FF00", Sessy= "#0000FF"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_light()

#stat_ellipse(type = "norm")

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$dog, vsd$treatment, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


library(vegan)
perm <- adonis(sampleDists ~ treatment*dog, 
               data = samples, 
               permutations = 10000, 
               method = "bray")
perm
