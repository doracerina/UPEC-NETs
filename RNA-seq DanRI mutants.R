# Project: RNA-seq analysis of UPEC danRI mutants (Figure 4)
# Description: Differential Expression (DE) analysis of WT vs DEL3 mutants.
#              Sequences processed via Galaxy.eu.
# Author: Dora Cerina


# 1. Setup
rm(list = ls()) 

# Load required packages
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(readr)
library(dplyr)
library(tibble)
library(openxlsx)
library(rtracklayer)
library(tidyverse)

# 2. Data import

# Wild Type (WT)
WT_1 <- read.delim("data/WT_1_featureCounts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))
WT_2 <- read.delim("data/WT_2_featureCounts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))
WT_3 <- read.delim("data/WT_3_featureCounts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))

# Mutant DEL3 - DanRI
DEL3_1 <- read.delim("data/DEL3_1_featureCounts_Counts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))
DEL3_2 <- read.delim("data/DEL3_2_featureCounts_Counts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))
DEL3_3 <- read.delim("data/DEL3_3_featureCounts_Counts.tabular", header=T, col.names = c('GeneID', 'GeneCount')) 

# Mutant DEL31 - DanR
DEL31_1 <- read.delim("data/DEL31_1_featureCounts_Counts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))
DEL31_2 <- read.delim("data/DEL31_2_featureCounts_Counts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))
DEL31_3 <- read.delim("data/DEL31_3_featureCounts_Counts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))

# Mutant DEL32 - DanI
DEL32_1 <- read.delim("data/DEL32_1_featureCounts_Counts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))
DEL32_2 <- read.delim("data/DEL32_2_featureCounts_Counts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))
DEL32_3 <- read.delim("data/DEL32_3_featureCounts_Counts.tabular", header=T, col.names = c('GeneID', 'GeneCount'))

# Gene Lengths and Assembly Annotations
GeneID_Lenght <- read.delim("data/Feature_lengths_UTI89.tabular", header=T, col.names = c('GeneID', 'GeneLenght'))
UTI89_assembly <- rtracklayer::import("data/genomic.gtf") %>% base::as.data.frame()

# 3. Data clean

# Isolate count vectors by removing GeneID column
WT_1 <- WT_1[,-1]; WT_2 <- WT_2[,-1]; WT_3 <- WT_3[,-1]
DEL3_1 <- DEL3_1[,-1]; DEL3_2 <- DEL3_2[,-1]; DEL3_3 <- DEL3_3[,-1] 
DEL31_1 <- DEL31_1[,-1]; DEL31_2 <- DEL31_2[,-1]; DEL31_3 <- DEL31_3[,-1]
DEL32_1 <- DEL32_1[,-1]; DEL32_2 <- DEL32_2[,-1]; DEL32_3 <- DEL32_3[,-1]

GeneID <- GeneID_Lenght$GeneID

# Construct main Count Dataframe
df_UPEC <- data.frame(WT_1, WT_2, WT_3, DEL3_1, DEL3_2, DEL3_3, 
                      DEL31_1, DEL31_2, DEL31_3, DEL32_1, DEL32_2, DEL32_3, 
                      row.names = GeneID)

# Define Metadata
Mutant <- c(rep('WT', 3), rep('DEL3', 3), rep('DEL31', 3), rep('DEL32', 3))
md_UPEC <- data.frame(Mutant, row.names = colnames(df_UPEC))

# Refine Genome Annotations
UTI89_assembly <- tibble(UTI89_assembly)
UTI89_genes <- UTI89_assembly %>% filter(type == 'gene') %>% select(-score, -phase, -transcript_id)
UTI89_cds_products <- UTI89_assembly %>% filter(type == 'CDS')

# 4. DESeq2 analysis

# DESeqDataSet
dds_UPEC <- DESeqDataSetFromMatrix(countData = df_UPEC, colData = md_UPEC, design = ~Mutant)

# Set 'WT' as the reference factor level
dds_UPEC$Mutant <- relevel(dds_UPEC$Mutant, ref = 'WT')

# Normalization
dds_UPEC <- estimateSizeFactors(dds_UPEC)
normalized_counts_UPEC <- counts(dds_UPEC, normalized=TRUE)
vsd_UPEC <- vst(dds_UPEC, blind=TRUE)

# 5. Visualzation

# PCA Plot
pca_UPEC <- plotPCA(vsd_UPEC, intgroup=c('Mutant'), returnData=TRUE) 
percentVar <- round(100 * attr(pca_UPEC, "percentVar"))

ggplot(pca_UPEC, aes(PC1, PC2, color=Mutant)) +
  geom_point(size=4) +
  scale_color_manual(values = c("WT" = "black", "DEL3" = "firebrick", "DEL31" = "darkorchid", "DEL32" = "darkcyan")) +
  theme_bw() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# 6. Differental expression

dds_UPEC <- DESeq(dds_UPEC)

# Contrast: WT vs DEL3
de_WTvsDEL3 <- results(dds_UPEC, contrast = c('Mutant', 'DEL3', 'WT'), alpha = 0.05, lfcThreshold = 0.5)
de_WTvsDEL3 <- lfcShrink(dds_UPEC, coef = 2, res = de_WTvsDEL3, type = 'apeglm')

# Contrast: WT vs DEL31
de_WTvsDEL31 <- results(dds_UPEC, contrast = c('Mutant', 'DEL31', 'WT'), alpha = 0.05, lfcThreshold = 0.5)
de_WTvsDEL31 <- lfcShrink(dds_UPEC, coef = 3, res = de_WTvsDEL31, type = 'apeglm')

# Contrast: WT vs DEL32
de_WTvsDEL32 <- results(dds_UPEC, contrast = c('Mutant', 'DEL32', 'WT'), alpha = 0.05, lfcThreshold = 0.5)
de_WTvsDEL32 <- lfcShrink(dds_UPEC, coef = 4, res = de_WTvsDEL32, type = 'apeglm')

# Merge with gene annotations
UTI89_merge <- merge(UTI89_genes, UTI89_cds_products, by = 'gene_id', all.x = TRUE) %>% 
  distinct(gene_id, .keep_all = TRUE) %>% 
  select(-ends_with('.y'))

# Results tables
de_WTvsDEL3_all  <- data.frame(de_WTvsDEL3) %>% rownames_to_column(var = 'gene_id') %>% left_join(UTI89_merge, by = 'gene_id')
de_WTvsDEL31_all <- data.frame(de_WTvsDEL31) %>% rownames_to_column(var = 'gene_id') %>% left_join(UTI89_merge, by = 'gene_id')
de_WTvsDEL32_all <- data.frame(de_WTvsDEL32) %>% rownames_to_column(var = 'gene_id') %>% left_join(UTI89_merge, by = 'gene_id')

theme_bw()

# 9. Export results
write_xlsx(de_WTvsDEL3_all, "results/WT_vs_DEL3_results.xlsx")
write_xlsx(de_WTvsDEL31_all, "results/WT_vs_DEL31_results.xlsx")
write_xlsx(de_WTvsDEL32_all, "results/WT_vs_DEL32_results.xlsx")
write_xlsx(as.data.frame(normalized_counts_UPEC), "results/normalized_counts.xlsx")
