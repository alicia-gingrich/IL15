library(tximport)
library(readr)
library(dplyr)
setwd ("C:/Users/Alicia/Desktop/il15")
# read in file names 
files <- read_csv("samples_dog_nonribo.csv")

# read in dog counts
tx2gene <- read_tsv("outputs/tx2gene/dogtx2gene.tsv")
canine <- tximport(files = files$files, type = "salmon", tx2gene = tx2gene)
caninecounts <- canine$counts

# change rownames for provenance
colnames(caninecounts) <- files$sample

# write to file
write.csv(caninecounts, "outputs/counts/dog_nonribo_counts.csv", quote = F, row.names = T)
