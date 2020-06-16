library(DESeq2)
library(dplyr)
library(readr)
library(conflicted)
library(tidyr)
library(tibble)
library(pheatmap)

## Import counts file
counts <- readRDS("counts_cleaned.rds")

## Import metadata file
metadata <- read.csv("APOE-TR_RNAseq_metadata.csv", header=T)
mousemeta <- read.csv("APOE-TR_mouse_metadata.csv", header =T)
biospecimen <- read.csv("APOE-TR_biospecimen_metadata.csv", header=T)

##adapt biospecimen data to show only the samples present in the rnaseq metadata
meta <- biospecimen[biospecimen$specimenID %in% metadata$specimenID, ]
  ## the information on the rnaseq metadata is not important for the analysis, so I won't join the two dataframes
##adapt mousedata to show only the samples present in the rnaseq metadata
mousemeta1 <- mousemeta[mousemeta$individualID %in% meta$individualID, ]
## put all info together
meta <- cbind(meta, mousemeta1)
## Delete unnecessary columns in clinical data
colnames(meta)
meta <- subset(meta, select =c("individualID", "specimenID", "tissue", "species", "sex", "genotype", "genotypeBackground"))

## Column names have . instead of _
colnames(counts) <- gsub('[.]','_', colnames(counts))

## Specimen ID column names after row 48 have - instead of _
meta$specimenID <- gsub('-','_', meta$specimenID)

##Important: columns in counts are similar to specimenID on meta
 
## Check if sample IDs are in the same order
all(meta$specimenID == colnames(counts))
    ## check where it doesn't match
match(colnames(counts), meta$specimenID)
    ## reorder counts columns so the order is now the same as in specimenID
idx <- match(meta$specimenID, colnames(counts))
data <- counts[ ,idx]
all(meta$specimenID == colnames(data))


### Because I'm only comparing E4 vs E3, i'll remove E2 samples
meta <- meta[meta$genotype == "APOE4" | meta$genotype == "APOE3" , ]
levels(meta) <- c("APOE4", "APOE3")
data <- data[ , colnames(data) %in% meta$specimenID]
all(meta$specimenID == colnames(data))

## Make sure Genotype and Sex are factors
meta$genotype <- factor(meta$genotype)
meta$sex <- factor(meta$sex)

## Subsetting to age of 3 months##
data3M <- data[ , grep("3M", colnames(data))]
meta3M <- meta[grep("3M", meta$specimenID) , ]

## Create DESeq object
dds3M <- DESeq2::DESeqDataSetFromMatrix(countData = data3M,
                              colData = meta3M,
                              design = ~ genotype)

##Principal Component Analysis
vsd <- DESeq2::vst(dds3M, blind = FALSE)
plotPCA(vsd, intgroup = "genotype")

#### Compare E3 vs E4 only among males
    ## For that I'll subset meta3M and data3M one more time
meta3M <- meta3M[meta3M$sex == "male" , ]
data3M <- data3M[ , colnames(data3M) %in% meta3M$specimenID]
rownames(meta3M) <- meta3M$specimenID
identical(colnames(data3M), rownames(meta3M))

####  RERUN DESEQ OBJECT ###
dds3M <- DESeq2::DESeqDataSetFromMatrix(countData = data3M,
                                        colData = meta3M,
                                        design = ~ genotype)

##Principal Component Analysis
vsd <- DESeq2::vst(dds3M, blind = FALSE)
plotPCA(vsd, intgroup = "genotype")

## Run DE analysis
dds3MO <- DESeq2::DESeq(dds3M)        

## Results APOE4 vs 3 
results3M <- DESeq2::results(dds3MO, contrast = c("genotype", "APOE4", "APOE3"))
results_table3 <- data.frame(results3M)

## turn rownames into column with gene symbols
results_table3 <- tibble::rownames_to_column(results_table3, "Gene_Symbol")

## See DE genes
results_DE3 <- results_table3[which(results_table3$padj < 0.05), ]
results_DE3 <- results_DE3[order(results_DE3$padj),]

## predicted genes
predicted3M <- results_DE3[grep("Gm", results_DE3$Gene_Symbol) , ]

## Export processed data
write.csv(results_DE3,"3months_male_APOE4vs3.csv")

## Select up/down regulated genes
upregulated <- dplyr::filter(results_DE3, log2FoldChange>0)
downregulated <- dplyr::filter(results_DE3, log2FoldChange<0)

write.csv(upregulated,"upregulated_3months_male_APOE4vs3.csv")

## Heatmap to show expression of MiR genes
mirna <- grep("Mir", results_DE3$Gene_Symbol)
mirna <- data.frame(results_DE3[mirna, ])
# Subset normalized counts to significant genes
ddsnorm3 <- estimateSizeFactors(dds3M)
normalized_counts <- counts(ddsnorm3, normalized=TRUE)
mirna_norm_counts <- as.data.frame(normalized_counts[mirna$Gene_Symbol, ]) 
colnames(mirna_norm_counts) <- mirna$Gene_Symbol
mirna_norm_counts <- data.frame(t(mirna_norm_counts))

library(RColorBrewer) 
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(mirna_norm_counts, 
         color = heat_colors, 
         cluster_rows = F, 
         cluster_cols = T,
         show_rownames = T,
         annotation_col = select(meta3M, genotype), 
         scale = "row")

## 12 Months ###
  ## Subsetting to age of 12 months
meta12M <- meta[grep("12M", meta$individualID) , ]
data12M <- data[ , colnames(data) %in% meta12M$specimenID]
  ## Check if sample IDs are in the same order
all(meta12M$specimenID == colnames(data12M))

## Create DESeq object
dds12M <- DESeq2::DESeqDataSetFromMatrix(countData = data12M,
                                        colData = meta12M,
                                        design = ~ genotype)

##Principal Component Analysis
vsd <- DESeq2::vst(dds12M, blind = FALSE)
plotPCA(vsd, intgroup = "genotype") 

## Checking PCA for each sex, male samples are more distinguished among genotypes compared to females,so I'll

## run the analysis for 12 months with only males
meta12M <- meta12M[meta12M$sex == "male" , ]
data12M <- data12M[ , colnames(data12M) %in% meta12M$specimenID]
####  RERUN DESEQ OBJECT ###
dds12M <- DESeq2::DESeqDataSetFromMatrix(countData = data12M,
                                        colData = meta12M,
                                        design = ~ genotype)
##Principal Component Analysis
vsd <- DESeq2::vst(dds12M, blind = FALSE)
plotPCA(vsd, intgroup = "genotype")

## Run DE analysis
dds12MO <- DESeq2::DESeq(dds12M)        

## Results
DESeq2::resultsNames(dds12MO)
## 4 vs 3 
results12M <- DESeq2::results(dds12MO, contrast = c("genotype", "APOE4", "APOE3"))
head(results12MO, n=10)
results_table <- data.frame(results12M)

## turn rownames into column with gene symbols
results_table <- tibble::rownames_to_column(results_table, "Gene_Symbol")

## See DE genes
results_DE12 <- results_table[which(results_table$padj < 0.05), ]
results_DE12 <- results_DE12[order(results_DE12$padj),]

## Export processed data
write.csv(results_DE,"12months_male_APOE4vs3.csv")

## predicted genes
predicted12M <- results_DE12[grep("Gm", results_DE12$Gene_Symbol) , ]

## Select up/down regulated genes
upregulated <- dplyr::filter(results_DE12, log2FoldChange>0)
downregulated <- dplyr::filter(results_DE12, log2FoldChange<0)

## Heatmap to show expression of MiR genes
mirna <- grep("Mir", results_DE12$Gene_Symbol)
mirna <- data.frame(results_DE12[mirna, ])
# Subset normalized counts to significant genes
ddsnorm12 <- estimateSizeFactors(dds12M)
normalized_counts <- counts(ddsnorm12, normalized=TRUE)
mirna_norm_counts <- as.data.frame(normalized_counts[mirna$Gene_Symbol, ]) 

rownames(meta12M) <- meta12M$specimenID

library(RColorBrewer) 
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(mirna_norm_counts, 
         color = heat_colors, 
         cluster_rows = F, 
         cluster_cols = T,
         show_rownames = T,
         annotation_col = dplyr::select(meta12M, genotype), 
         scale = "row")

## 24 Months ###
  ## Subsetting to age of 24 months
meta24M <- meta[grep("24M", meta$individualID) , ]
data24M <- data[ , colnames(data) %in% meta24M$specimenID]
## Check if sample IDs are in the same order
all(meta24M$specimenID == colnames(data24M))

## Create DESeq object
dds24M <- DESeq2::DESeqDataSetFromMatrix(countData = data24M,
                                         colData = meta24M,
                                         design = ~ genotype)

##Principal Component Analysis
vsd <- DESeq2::vst(dds24M, blind = FALSE)
plotPCA(vsd, intgroup = "genotype") 

## Only males (could've been females, samples were clustered according to genotype, chose male to compare)
meta24M <- meta24M[meta24M$sex == "male" , ]
data24M <- data24M[ , colnames(data24M) %in% meta24M$specimenID]

## Create DESeq object
dds24M <- DESeq2::DESeqDataSetFromMatrix(countData = data24M,
                                         colData = meta24M,
                                         design = ~ genotype)

##Principal Component Analysis
vsd <- DESeq2::vst(dds24M, blind = FALSE)
plotPCA(vsd, intgroup = "genotype") 

## Run DE analysis
dds24MO <- DESeq2::DESeq(dds24M)        

## Results
DESeq2::resultsNames(dds24MO)
## 4 vs 3 
results24M <- DESeq2::results(dds24MO, contrast = c("genotype", "APOE4", "APOE3"))
results_table <- data.frame(results24M)

## turn rownames into column with gene symbols
results_table <- tibble::rownames_to_column(results_table, "Gene_Symbol")

## See DE genes
results_DE24 <- results_table[which(results_table$padj < 0.05), ]
results_DE24 <- results_DE24[order(results_DE24$padj),]

## Export processed data
write.csv(results_DE24,"24months_male_APOE4vs3.csv")

## predicted genes
predicted24M <- results_DE24[grep("Gm", results_DE24$Gene_Symbol) , ]

## Select up/down regulated genes
upregulated <- dplyr::filter(results_DE24, log2FoldChange>0)
downregulated <- dplyr::filter(results_DE24, log2FoldChange<0)

## Heatmap to show expression of MiR genes
mirna <- grep("Mir", results_DE24$Gene_Symbol)
mirna <- data.frame(results_DE24[mirna, ])
# Subset normalized counts to significant genes
ddsnorm24 <- estimateSizeFactors(dds24M)
normalized_counts <- counts(ddsnorm24, normalized=TRUE)
mirna_norm_counts <- as.data.frame(normalized_counts[mirna$Gene_Symbol, ]) 

rownames(meta24M) <- meta24M$specimenID

library(RColorBrewer) 
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(mirna_norm_counts, 
         color = heat_colors, 
         cluster_rows = F, 
         cluster_cols = T,
         show_rownames = T,
         annotation_col = dplyr::select(meta24M, genotype), 
         scale = "row")

## miRNA analysis ##
library(multiMiR)

## 3 months, male, APOE4 vs APOE3 ##
    ## mir676 is upregulated
        ##separate list with downregulated genes (upregulated miRNA)
downregulated3 <- dplyr::filter(results_DE3, log2FoldChange<0)
  ## Adapting the dataframes - follow miRNA notation + column with only gene names/IDs
mirna <- c("mmu-miR-676-3p", "mmu-miR-676-5p")
down3 <- as.data.frame( downregulated3[,1])
  ## retrieve interactions between a set of miRNAs and a set of genes
mirnagene <- get_multimir(org="mmu", mirna=mirna, target=down3, table="all", summary=TRUE, use.tibble = T)
summary(mirnagene@data)    
## To see how many predicted and validated interactions
table(mirnagene@data$type)
multimir <- data.frame(mirnagene@data)
write.csv(multimir,"multimir_3M_upmirna_downgene.csv")

## 12 months, male, APOE4 vs APOE3 ##
## mir128-2, MIR99ahg is upregulated
## mir700, mir3093, mir411, mir7227, mir190b, mir200c downregulated
##separate list with down/upregulated genes
downregulated12 <- dplyr::filter(results_DE12, log2FoldChange<0)
upregulated12 <- dplyr::filter(results_DE12, log2FoldChange>0)

## Adapting the dataframes - follow miRNA notation + column with only gene names/IDs
upmirna <- c("mmu-miR-128-3p", "mmu-miR-128-2-5p", "mmu-miR-99a-3p", "mmu-miR-99a-5p", "mmu-mir-125b-2-3p", "mmu-mir-125b-2-5p", "mmu-let7c-5p") 
downmirna <- c("mmu-miR-200c-3p","mmu-miR-200c-5p", "mmu-miR-700-5p", "mmu-miR-700-3p", "mmu-miR-3093-3p","mmu-miR-3093-5p", "mmu-miR-190b-3p", "mmu-miR-190b-5p", "mmu-miR-7227-5p", "mmu-miR-7227-3p","mmu-miR-411-5p", "mmu-miR-411-3p")
down12 <- as.data.frame( downregulated12[,1])
up12 <- as.data.frame( upregulated12[,1])

## retrieve interactions between a set of miRNAs and a set of genes
upmirna_downgene <- get_multimir(org="mmu", mirna=upmirna, target=down12, table="all", summary=TRUE, predicted.cutoff.type="p", predicted.cutoff=10, use.tibble = T)
summary(upmirna_downgene@data)    

downmirna_upgene <- get_multimir(org="mmu", mirna=downmirna, target=up12, table="all", summary=TRUE, predicted.cutoff.type="p", predicted.cutoff=10, use.tibble = T)
summary(downmirna_upgene@data)

## To see how many predicted and validated interactions
table(upmirna_downgene@data$type)
table(downmirna_upgene@data$type)

multimir_12M_downmirna_upgene <- data.frame(downmirna_upgene@data)
multimir_12M_upmirna_downgene <- data.frame(upmirna_downgene@data)

write.csv(multimir_12M_downmirna_upgene,"multimir_12M_downmirna_upgene.csv")
write.csv(multimir_12M_upmirna_downgene,"multimir_12M_upmirna_downgene.csv")


## 24 months, male, APOE4 vs APOE3 ##
## mir7b is upregulated
## mir7227, mir99ahg downregulated
##separate list with down/upregulated genes
downregulated24 <- dplyr::filter(results_DE24, log2FoldChange<0)
upregulated24 <- dplyr::filter(results_DE24, log2FoldChange>0)

## Adapting the dataframes - follow miRNA notation + column with only gene names/IDs
upmirna <- c("mmu-miR-7b-3p", "mmu-miR-7b-5p") 
downmirna <- c("mmu-miR-7227-5p", "mmu-miR-7227-3p", "mmu-miR-99a-3p", "mmu-miR-99a-5p", "mmu-mir-125b-2-3p", "mmu-mir-125b-2-5p", "mmu-let7c-5p")
down24 <- as.data.frame( downregulated24[,1])
up24 <- as.data.frame( upregulated24[,1])

## retrieve interactions between a set of miRNAs and a set of genes
upmirna_downgene <- get_multimir(org="mmu", mirna=upmirna, target=down24, table="all", summary=TRUE, predicted.cutoff.type="p", predicted.cutoff=10, use.tibble = T)
summary(upmirna_downgene@data)    

downmirna_upgene <- get_multimir(org="mmu", mirna=downmirna, target=up24, table="all", summary=TRUE, predicted.cutoff.type="p", predicted.cutoff=10, use.tibble = T)
summary(downmirna_upgene@data)

## To see how many predicted and validated interactions
table(upmirna_downgene@data$type)
table(downmirna_upgene@data$type)

multimir_24M_downmirna_upgene <- data.frame(downmirna_upgene@data)
multimir_24M_upmirna_downgene <- data.frame(upmirna_downgene@data)

write.csv(multimir_24M_downmirna_upgene,"multimir_24M_downmirna_upgene.csv")
write.csv(multimir_24M_upmirna_downgene,"multimir_24M_upmirna_downgene.csv")
