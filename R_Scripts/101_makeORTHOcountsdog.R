library(readr)
library(dplyr)


# ortho_refseq <- read_delim("C:/Users/Alicia/Desktop/il15/HSA_CFA_Ortho_RefSeq.txt",
#                            col_names = c("Human", "Dog", "Pair", "RefSeqID"),
#                            delim = " ")
#                            
# ortho_refseq$Human <- trimws(ortho_refseq$Human, which = c("both"))
# ortho_refseq$Dog <- trimws(ortho_refseq$Dog, which = c("both"))
# ortho_refseq$Pair <- trimws(ortho_refseq$Pair, which = c("both"))
# ortho_refseq$RefSeqID <- trimws(ortho_refseq$RefSeqID, which = c("both"))
# 
# head(ortho_refseq)

ortho_entrez <- read_tsv("C:/Users/Alicia/Desktop/il15/HSA_CFA_Ortho_Entrez.txt",
                           col_names = c("Human", "Dog", "Pair", "entrezID"))

head(ortho_entrez)

#filter ortho_entrez 1:1
ortho_entrez <- filter(ortho_entrez, Pair=="1:1")
table(is.na(ortho_entrez$entrezID))

feat_table_dog <-read_tsv("C:/Users/Alicia/Desktop/il15/outputs/tx2gene/dogtx2gene.tsv")
head(feat_table_dog)
table(ortho_entrez$Dog %in% feat_table_dog$GeneID)
table(feat_table_dog$GeneID %in% ortho_entrez$Dog)

head(ortho_entrez$Dog %in% feat_table_dog$GeneID, n=100)
head(feat_table_dog$GeneID %in% ortho_entrez$Dog, n=100)

ortho_entrez$Dog[7]
feat_table_dog[6,]
feat_table_dog$GeneID <- as.character(feat_table_dog$GeneID)

#join ortho_entrez and feat_table_dog
feat_table_dog_ortho <- left_join(feat_table_dog, ortho_entrez, by = c("GeneID"="Dog"))
head(feat_table_dog)
head(ortho_entrez)
head(feat_table_dog_ortho)

#filter to symbol, entrezID
feat_table_dog_ortho <- feat_table_dog_ortho%>%
  select(symbol, entrezID)%>%
  distinct()%>%
  filter(!is.na(entrezID))
head(feat_table_dog_ortho)

dog_counts <- read.csv("C:/Users/Alicia/Desktop/il15/outputs/counts/dog_nonribo_counts.csv")
head(dog_counts)

#join table with symbols and entrezID for all dog genes with human 1:1 orthos, to dog counts
dog_ortho_counts <- inner_join(dog_counts, feat_table_dog_ortho, by = c("X"="symbol"))
head(dog_ortho_counts)

write.csv(as.data.frame(dog_ortho_counts), 
          file="dog_ortho_counts.csv")

all_ortho_counts <- inner_join(dog_ortho_counts, human_ortho_counts, by = "entrezID")

head(all_ortho_counts)

write.csv(as.data.frame(all_ortho_counts), 
          file="all_ortho_counts.csv")
