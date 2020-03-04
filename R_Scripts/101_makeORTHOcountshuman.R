library(readr)
library(dplyr)

#read in Human and Dog Orthologous genes list from OMA with Entrez notation 
ortho_entrez <- read_tsv("C:/Users/Alicia/Desktop/il15/HSA_CFA_Ortho_Entrez.txt",
                         col_names = c("Human", "Dog", "Pair", "entrezID"))
head(ortho_entrez)

#filter ortho_entrez for 1:1 matches 
ortho_entrez <- filter(ortho_entrez, Pair=="1:1")
table(is.na(ortho_entrez$entrezID))
head(ortho_entrez)

#read in tx2gene human feature table
feat_table_human <-read_tsv("C:/Users/Alicia/Desktop/il15/outputs/tx2gene/humantx2gene.tsv")
head(feat_table_human)

#check tables for overlap
table(ortho_entrez$Human %in% feat_table_human$GeneID)

table(feat_table_human$GeneID %in% ortho_entrez$Human)


###use this if you need to interrogate individual cells for values
#head(ortho_entrez$Dog %in% feat_table_dog$GeneID, n=100)
#head(feat_table_dog$GeneID %in% ortho_entrez$Dog, n=100)
#ortho_entrez$Dog[7]
#feat_table_dog[6,]

#make GeneID column in feat_table_human a character instead of numeric 
feat_table_human$GeneID <- as.character(feat_table_human$GeneID)

#join ortho_entrez and feat_table_human to make a table with tx2gene and entrezID of 1:1 ortholog pairs
feat_table_human_ortho <- left_join(feat_table_human, ortho_entrez, by = c("GeneID"="Human"))
head(feat_table_human_ortho)

#filter this new table to symbol, entrezID
feat_table_human_ortho <- feat_table_human_ortho%>%
  select(symbol, entrezID)%>%
  distinct()%>%
  filter(!is.na(entrezID))
head(feat_table_human_ortho)
#should have 13716 obs. of 2 variables

#read in human counts from local computer
human_counts <- read.csv("C:/Users/Alicia/Desktop/il15/outputs/counts/human_nonribo_counts.csv")
head(human_counts)

#join table with symbols and entrezID for all human genes with dog 1:1 orthos, to human counts
human_ortho_counts <- inner_join(human_counts, feat_table_human_ortho, by = c("X"="symbol"))
head(human_ortho_counts)
#should have 13704 obs. of 14 variables

write.csv(as.data.frame(human_ortho_counts), 
          file="human_ortho_counts.csv")



###This section to use RefSeq gene annotation from OMA
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