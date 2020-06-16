### Bar plot -- Pathway Enrichment Analysis ##
library(ggplot2)

## 12 months- down miRNAs, up genes
downmir_upgene_12 <- read.csv("gProfiler_12M_downmir_upgenes.csv", header=T)
downmir_upgene_12 <- subset(downmir_upgene_12, select =c("source", "term_name", "adjusted_p_value", "negative_log10_of_adjusted_p_value"))

## select top 10
downmir_upgene_12 <- downmir_upgene_12[1:10, ]

# Plot Graph
p <-ggplot(data=downmir_upgene_12, aes(x= reorder(term_name, negative_log10_of_adjusted_p_value), y=negative_log10_of_adjusted_p_value)) + geom_bar(stat="identity", width=0.5, fill="black") + xlab("- log10(adjusted p-value)") + theme(aspect.ratio = 5/1) + theme_bw() + theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank()) 
p <- p + coord_flip()
p <- p + theme(axis.title.y = element_blank())
p <- p + theme(axis.title.x = element_blank())
p

## 12 months- up miRNAs, down genes
upmir_downgene_12 <- read.csv("gProfiler_12M_upmir_downgenes.csv", header=T)
upmir_downgene_12 <- subset(upmir_downgene_12, select =c("source", "term_name", "adjusted_p_value", "negative_log10_of_adjusted_p_value"))

## select top 10
upmir_downgene_12 <- upmir_downgene_12[1:10, ]

# Plot Graph
p <-ggplot(data=upmir_downgene_12, aes(x= reorder(term_name, negative_log10_of_adjusted_p_value), y=negative_log10_of_adjusted_p_value)) + geom_bar(stat="identity", width=0.5, fill="black") + xlab("- log10(adjusted p-value)") + theme(aspect.ratio = 5/1) + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                                                                                                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                   panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                                                                   panel.border = element_blank(),
                                                                                                                                                                                                                                                                                   panel.background = element_blank()) 
p <- p + coord_flip()
p <- p + theme(axis.title.y = element_blank())
p <- p + theme(axis.title.x = element_blank())
p

## 24 months- down miRNAs, up genes
downmir_upgene_24 <- read.csv("gProfiler_24M_downmir_upgenes.csv", header=T)
downmir_upgene_24 <- subset(downmir_upgene_24, select =c("source", "term_name", "adjusted_p_value", "negative_log10_of_adjusted_p_value"))

# Plot Graph
p <-ggplot(data=downmir_upgene_24, aes(x= reorder(term_name, negative_log10_of_adjusted_p_value), y=negative_log10_of_adjusted_p_value)) + geom_bar(stat="identity", width=0.5, fill="black") + xlab("- log10(adjusted p-value)") + theme(aspect.ratio = 5/1) + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                                                                                                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                   panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                                                                   panel.border = element_blank(),
                                                                                                                                                                                                                                                                                   panel.background = element_blank()) 
p <- p + coord_flip()
p <- p + theme(axis.title.y = element_blank())
p <- p + theme(axis.title.x = element_blank())
p

## MEYER2019/ORGANOIDS- down miRNAs, up genes
downmir_upgene <- read.csv("gProfiler_downmirna_upgenes.csv", header=T)
downmir_upgene <- subset(downmir_upgene, select =c("source", "term_name", "adjusted_p_value", "negative_log10_of_adjusted_p_value"))

## select top 10
downmir_upgene <- downmir_upgene[1:10, ]

# Plot Graph
p <-ggplot(data=downmir_upgene, aes(x= reorder(term_name, negative_log10_of_adjusted_p_value), y=negative_log10_of_adjusted_p_value)) + geom_bar(stat="identity", width=0.5, fill="black") + xlab("- log10(adjusted p-value)") + theme(aspect.ratio = 5/1) + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                                                                                                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                   panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                                                                   panel.border = element_blank(),
                                                                                                                                                                                                                                                                                   panel.background = element_blank()) 
p <- p + coord_flip()
p <- p + theme(axis.title.y = element_blank())
p <- p + theme(axis.title.x = element_blank())
p

## MEYER2019/ORGANOIDS- up miRNAs, down genes
upmir_downgene <- read.csv("gProfiler_upmirna_downgenes.csv", header=T)
upmir_downgene <- subset(upmir_downgene, select =c("source", "term_name", "adjusted_p_value", "negative_log10_of_adjusted_p_value"))

## select top 10
upmir_downgene <- upmir_downgene[1:10, ]

# Plot Graph
p <-ggplot(data=upmir_downgene, aes(x= reorder(term_name, negative_log10_of_adjusted_p_value), y=negative_log10_of_adjusted_p_value)) + geom_bar(stat="identity", width=0.5, fill="black") + xlab("- log10(adjusted p-value)") + theme(aspect.ratio = 5/1) + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                                                                                                                                                                                                panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                                                                panel.border = element_blank(),
                                                                                                                                                                                                                                                                                panel.background = element_blank()) 
p <- p + coord_flip()
p <- p + theme(axis.title.y = element_blank())
p <- p + theme(axis.title.x = element_blank())
p

## NURIEL
downmir_upgene <- read.csv("gProfiler_downmirna_upgenes.csv", header=T)
downmir_upgene <- subset(downmir_upgene, select =c("source", "term_name", "adjusted_p_value", "negative_log10_of_adjusted_p_value"))

## select top 10
downmir_upgene <- downmir_upgene[1:10, ]

# Plot Graph
p <-ggplot(data=downmir_upgene, aes(x= reorder(term_name, negative_log10_of_adjusted_p_value), y=negative_log10_of_adjusted_p_value)) + geom_bar(stat="identity", width=0.5, fill="black") + xlab("- log10(adjusted p-value)") + theme(aspect.ratio = 5/1) + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                                                                                                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                   panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                                                                   panel.border = element_blank(),
                                                                                                                                                                                                                                                                                   panel.background = element_blank()) 
p <- p + coord_flip()
p <- p + theme(axis.title.y = element_blank())
p <- p + theme(axis.title.x = element_blank())
p
