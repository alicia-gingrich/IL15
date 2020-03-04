###001_makecounts_human###
library(tximport)
library(readr)
library(dplyr)
setwd ("C:/Users/Alicia/Desktop/il15")
# read in file names 
files <- read_csv("samples_human_nonribo.csv")

# read in dog counts
tx2gene_hs <- read_tsv("outputs/tx2gene/humantx2gene.tsv")
human <- tximport(files = files$files, type = "salmon", tx2gene = tx2gene_hs)
humancounts <- human$counts

# change rownames for provenance
colnames(humancounts) <- files$sample

# write to file
write.csv(humancounts, "outputs/counts/human_nonribo_counts.csv", quote = F, row.names = T)
