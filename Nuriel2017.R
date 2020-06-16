library(DESeq2)
library(dplyr)
library(tibble)
library(pheatmap)

## Import counts matrix
counts <- read.table("raw_counts_pseudofilter-3.txt", header=T)

## Import metadata
meta <- read.csv("design_updated-genotypes_minus-E44.csv", header=T)

## Make genes appear in rownames for count matrix/ sample IDs in rows for meta
rownames(counts) <- counts$geneSymbol
counts <- counts[ ,-1]
rownames(meta) <- meta$sample_id

## Use only EC samples
meta <- meta[meta$Region == "EC", ]
meta$Genotype <- gsub('/','_', meta$Genotype)
meta$Genotype <- factor(meta$Genotype)
levels(meta) <- c("E3/4", "E3/3")

counts <- counts[ ,colnames(counts) %in% rownames(meta)]
all(rownames(meta) == colnames(counts))

## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Genotype)

##Principal Component Analysis
vsd <- DESeq2::vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Genotype")

## Run DE analysis
dds1 <- DESeq(dds)

## Results- 4 vs 3
resultsNames(dds1)
results <- results(dds1, contrast = c("Genotype", "E3_4", "E3_3"))
results_table <- data.frame(results)

## See DE genes
results_DE <- results_table[which(results_table$padj < 0.05), ]
results_DE <- results_DE[order(results_DE$padj),]

results_DE <- tibble::rownames_to_column(results_DE, "Gene_Symbol")

## Export results
write.csv(results_DE,"APOE4vs3_fromcountmatrix.csv")

## Upregulated/Downregulated
upregulated <- dplyr::filter(results_DE, log2FoldChange>0)
downregulated <- dplyr::filter(results_DE, log2FoldChange<0)

## Heatmap to show expression of MiR genes
mirna <- grep("Mir", results_DE$Gene_Symbol)
mirna <- data.frame(results_DE[mirna, ])
# Subset normalized counts to significant genes
ddsnorm <- estimateSizeFactors(dds)
normalized_counts <- counts(ddsnorm, normalized=TRUE)
mirna_norm_counts <- as.data.frame(normalized_counts[mirna$Gene_Symbol, ]) 

library(RColorBrewer) 
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(mirna_norm_counts, 
         color = heat_colors, 
         cluster_rows = F, 
         cluster_cols = T,
         show_rownames = T,
         annotation_col = select(meta, Genotype), 
         scale = "row")

## predicted genes
predicted <- results_DE[grep("Gm", results_DE$Gene_Symbol) , ]

##multiMiR results 
library(multiMiR)
  ## Mir494 and Mir486 are downregulated ##
## Adapting the dataframes - follow miRNA notation + column with only gene names/IDs
mirna_list <- c("mmu-miR-486a-5p", "mmu-miR-486a-3p","mmu-miR-494-3p", "mmu-miR-494-5p")
upregulated_list <- as.data.frame( upregulated[,1])

  ## retrieve interactions between a set of miRNAs and a set of genes
mirnagene <- get_multimir(org="mmu", mirna=mirna_list, target=upregulated_list, table="all", summary=TRUE, use.tibble = T)

## To see how many predicted and validated interactions
table(mirnagene@data$type)

head(mirnagene@data)
output <- data.frame(mirnagene@data)

## Export dataframe as csv file
write.csv(output,"Nuriel_multiMiR_output.csv")

