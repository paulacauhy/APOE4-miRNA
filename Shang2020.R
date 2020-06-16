library(GEOquery)
library(DESeq2)
library(dplyr)
library(readr)
library(conflicted)

## Import series matrix file from GEO
gse <- getGEO(filename="GSE140205_series_matrix.txt")
View(gse)
metadata <- data.frame(gse)

## Import counts 
data <- read.table("GSE140205_counts_tximport_4groups.txt", sep = ",", header=T)

## Rename rownames in metadata
rownames(metadata) <- metadata$title
## Subset metadata with relevant columns
metadata <- subset(metadata, select =c("title", "age.ch1", "genotype.ch1", "Sex.ch1", "strain.ch1", "tissue.ch1"))

## convert ensembl ids into rownames -- data
counts <- data[-1]
row.names(counts) <- data$X
head(counts)

## Need sample IDs to be in the same order in counts and metadata
## Check if sample IDs are in the same order
all(rownames(metadata) == colnames(counts))
   
## Make sure that genotype and sex are factors
metadata <- data.frame(Genotype = metadata$genotype.ch1, metadata)
metadata <- data.frame(Sex = metadata$Sex.ch1, metadata)
metadata$Genotype <- factor(metadata$Genotype)
metadata$Sex <- factor(metadata$Sex)

## Create DESeq object
    ## on the paper they put sex first on the design, I'll put the genotype
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Sex + Genotype + Sex:Genotype)

## Run DE analysis
dds1 <- DESeq(dds)

## Results 4 vs 3 
results <- results(dds1, contrast = c("Genotype", "APOE4", "APOE3"))
results_table <- data.frame(results)

  ## turn rownames into column with gene symbols
results_table <- tibble::rownames_to_column(results_table, "Gene_Symbol")

## See DE genes
results_DE <- results_table[which(results_table$padj < 0.05), ]
results_DE <- results_DE[order(results_DE$padj),]

## Turn unknown symbols into gene names
results_DE$Gene_Symbol <- gsub('UNKOWN_','', results_DE$Gene_Symbol)
library(biomaRt)
listMarts()
ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
attributes
unknown <- getBM(attributes=c("ensembl_transcript_id",'ensembl_gene_id', "mgi_symbol"), 
      filters = 'ensembl_transcript_id', 
      values = results_DE$Gene_Symbol, 
      mart = ensembl)

results_DE$Gene_Symbol <- gsub("ENSMUST00000215334", "Gm48348", results_DE$Gene_Symbol)

## Export processed data
write.csv(results_DE,"APOE4vs3.csv")

## Select upregulated genes
upregulated <- dplyr::filter(results_DE, log2FoldChange>0)
upregulated$Gene_Symbol <- gsub('UNKOWN_','', upregulated$Gene_Symbol)
write.csv(upregulated,"upregulated_Shang.csv")

## predicted genes
predicted <- results_DE[grep("Gm", results_DE$Gene_Symbol) , ]

## Heatmap to show expression of MiR genes
mirna <- grep("Mir", results_DE$Gene_Symbol)
mirna <- data.frame(results_DE[mirna, ])
rownames(mirna) <- mirna$Gene_Symbol
# Subset normalized counts to significant genes
ddsnorm <- estimateSizeFactors(dds)
normalized_counts <- counts(ddsnorm, normalized=TRUE)
mirna_norm_counts <- as.data.frame(normalized_counts[mirna$Gene_Symbol, ]) 
colnames(mirna_norm_counts) <- mirna$Gene_Symbol
mirna_norm_counts <- data.frame(t(mirna_norm_counts))

library(pheatmap)
library(RColorBrewer) 
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(mirna_norm_counts, 
         color = heat_colors, 
         cluster_rows = F, 
         cluster_cols = T,
         show_rownames = T,
         annotation_col = select(metadata, Genotype), 
         scale = "row")


#### Only male mice ####
metaM <- metadata[metadata$Sex == "Male", ]
countsM <- counts[ , colnames(counts) %in% rownames(metaM)]
identical(colnames(countsM), rownames(metaM))

metaM$Genotype <- factor(metaM$Genotype)

## Create DESeq object
ddsM <- DESeqDataSetFromMatrix(countData = countsM,
                              colData = metaM,
                              design = ~ Genotype)

## Run DE analysis
ddsMA <- DESeq(ddsM)

## Results
  ## 4 vs 3 
results <- results(ddsMA, contrast = c("Genotype", "APOE4", "APOE3"))
head(results, n=10)
results_table <- data.frame(results)

## turn rownames into column with gene symbols
results_table <- tibble::rownames_to_column(results_table, "Gene_Symbol")

## See DE genes
results_DE <- results_table[which(results_table$padj < 0.05), ]
results_DE <- results_DE[order(results_DE$padj),]

## predicted genes
predicted <- results_DE[grep("Gm", results_DE$Gene_Symbol) , ]


## multimiR
library(multiMiR)
  ## Mir99ahg is downregulated
    ##separate list with upregulated genes 
upregulated <- dplyr::filter(results_DE, log2FoldChange>0)
  ## Adapting the dataframes - follow miRNA notation + column with only gene names/IDs
mirna <- c("mmu-miR-99a-3p", "mmu-miR-99a-5p", "mmu-mir-125b-2-3p", "mmu-mir-125b-2-5p","mmu-let7c-5p")
up <- as.data.frame( upregulated[,1])

## retrieve interactions between a set of miRNAs and a set of genes
mirnagene <- get_multimir(org="mmu", mirna=mirna, target=up, table="all", summary=TRUE, use.tibble = T)
   
## To see how many predicted and validated interactions
table(mirnagene@data$type)
multimir <- data.frame(mirnagene@data)
write.csv(multimir,"multimir_shang_downmirna_upgene.csv")


