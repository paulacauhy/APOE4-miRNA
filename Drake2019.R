##Meyer, 2019  -- APOE4/3 organoids
library(GEOquery)
library(DESeq2)
library(pheatmap)
library(dplyr)

## Import raw counts
counts <- read.table( "GSE117588_organoids-gene-counts.tsv", header=T, sep='\t' )

## Import metadata
gse <- getGEO(filename="GSE117588_series_matrix.txt")
metadata <- data.frame(gse)

## Adapt counts file so gene IDs are row names and not in a column
rownames(counts) <- counts$geneId
counts <- counts[, -1]

## Change rownames in metadata, so they match sample names on counts
rownames(metadata) <- metadata$title

## Exclude unnecessary columns on metadata
metadata <- subset(metadata, select= c("title", "geo_accession", "tissue.ch1", "apoe.genotype.ch1"))

## Check if samples are in the same order in the two files
all(rownames(metadata) == colnames(counts))

## Change genotype name to only E3E3 and E4E4- look better on PCA
metadata <- data.frame(Genotype = metadata$apoe.genotype.ch1, metadata)
metadata$Genotype <- gsub('_CTR','', metadata$Genotype)
metadata$Genotype <- gsub('_ISO','', metadata$Genotype)

## Make sure that genotype is a factor
metadata$Genotype <- factor(metadata$Genotype)

## Create DESeq object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Genotype)
## Run DE analysis
dds1 <- DESeq2::DESeq(dds)

## Results
DESeq2::resultsNames(dds1)
## 4 vs 3 
results <- DESeq2::results(dds1, contrast = c("Genotype", "E4E4", "E3E3"))
head(results, n=10)
results_table <- data.frame(results)

# convert ensembl ids to gene symbols
library(AnnotationDbi)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
annot <- AnnotationDbi ::select(org.Hs.eg.db, keys=rownames(results),
                                columns=c("SYMBOL"), keytype="ENSEMBL")
results_table <- merge(results_table, annot, by.x=0, by.y="ENSEMBL")
head(results_table)

## See DE genes
results_DE <- results_table[which(results_table$padj < 0.05), ]
results_DE <- results_DE[order(results_DE$padj),]

## Up/ Down regulated genes
upregulated <- dplyr::filter(results_DE, log2FoldChange>0)
downregulated <- dplyr::filter(results_DE, log2FoldChange<0)

## Export processed data
write.csv(results_DE,"APOE4vs3.csv")

##Principal Component Analysis
vsd <- DESeq2::vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Genotype")

## Heatmap to show expression of MiR genes
mirna <- grep("MIR", results_DE$SYMBOL)
mirna <- data.frame(results_DE[mirna, ])
  # Subset normalized counts to significant genes
ddsnorm <- estimateSizeFactors(dds)
normalized_counts <- counts(ddsnorm, normalized=TRUE)
mirna_norm_counts <- as.data.frame(normalized_counts[mirna$Row.names, ])

annot1 <- AnnotationDbi ::select(org.Hs.eg.db, keys=rownames(mirna_norm_counts),
                                columns=c("SYMBOL"), keytype="ENSEMBL")
mirna_norm_counts <- merge(mirna_norm_counts, annot1, by.x=0, by.y="ENSEMBL")
rownames(mirna_norm_counts) <- mirna_norm_counts$SYMBOL
mirna_norm_counts <-mirna_norm_counts[ ,-1]
mirna_norm_counts <-mirna_norm_counts[ ,-7]

library(RColorBrewer) 
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(mirna_norm_counts, 
         color = heat_colors, 
         cluster_rows = F, 
         cluster_cols = T,
         show_rownames = T,
         annotation_col = select(metadata, Genotype), 
         scale = "row")

## miRNA-mRNA interactions - multimiR
library(multiMiR)

## Adapting the dataframes - follow miRNA notation + column with only gene names/IDs
  ## Upregulated mirnas and down genes
upmirna <- c("hsa-miR-99a-5p", "hsa-miR-125b-2-3p", "hsa-let-7c-5p","hsa-miR-9-5p", "hsa-miR-124-3p","hsa-miR-124-5p", "hsa-miR-34a-3p","hsa-miR-34a-5p", "hsa-miR-3619-3p", "hsa-miR-3619-5p")
downgenes <- as.data.frame( downregulated[,8])

## downregulated mirnas and up genes
downmirna <- c("hsa-miR-3648", "hsa-miR-22-3p", "hsa-miR-22-5p")
upgenes <- as.data.frame( upregulated[,8])

## retrieve interactions between a set of miRNAs and a set of genes
  ## Up miRNAs and Down genes
mirnagene <- get_multimir(org="hsa", mirna=upmirna, target=downgenes, table="all", summary=TRUE, predicted.cutoff.type="p", predicted.cutoff=10, use.tibble = T)
table(mirnagene@data$type)
UPmirna_DOWNgene <- data.frame(mirnagene@data)

  ## Down miRNAs and Up genes
mirnagene1 <- get_multimir(org="hsa", mirna=downmirna, target=upgenes, table="all", summary=TRUE, predicted.cutoff.type="p", predicted.cutoff=10, use.tibble = T)
table(mirnagene1@data$type)
DOWNmirna_UPgene <- data.frame(mirnagene1@data)

write.csv(DOWNmirna_UPgene,"DOWNmirna_UPgene.csv")
write.csv(UPmirna_DOWNgene,"UPmirna_DOWNgene.csv")
