# Project: RNA-seq analysis of UPEC with Neutrophil Extracellular Traps (NETs) for Figure 1
# Description: Differential Expression (DE) analysis using DESeq2. 
#              Raw counts generated via Galaxy.eu.
# Author: Dora Cerina

# 1. Setup
rm(list = ls()) 

# Load required packages
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(DESeq2)
library(rtracklayer)
library(ggplot2)
library(ggrepel)
library(stringr) 
library(arsenal)
library(VennDiagram)
library(writexl)

# 2. Data import

# Unstimulated UPEC
UPEC_1_2_1h <- read.delim('data/1_2_1h_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_2_2_1h <- read.delim('data/2_2_1h_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_3_2_1h <- read.delim('data/3_2_1h_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_4_2_1h <- read.delim('data/4_2_1h_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_5_2_1h <- read.delim('data/5_2_1h_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))

# UPEC + PMA
UPEC_PMA_1_3_1h <- read.delim('data/1_3_1h_UPEC_PMA_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_PMA_2_3_1h <- read.delim('data/2_3_1h_UPEC_PMA_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_PMA_3_3_1h <- read.delim('data/3_3_1h_UPEC_PMA_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_PMA_4_3_1h <- read.delim('data/4_3_1h_UPEC_PMA_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_PMA_5_3_1h <- read.delim('data/5_3_1h_UPEC_PMA_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))

# UPEC + DNase
UPEC_DNase_1_4_1h <- read.delim('data/1_4_1h_UPEC_DNase_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_DNase_2_4_1h <- read.delim('data/2_4_1h_UPEC_DNase_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_DNase_3_4_1h <- read.delim('data/3_4_1h_UPEC_DNase_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_DNase_4_4_1h <- read.delim('data/4_4_1h_UPEC_DNase_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_DNase_5_4_1h <- read.delim('data/5_4_1h_UPEC_DNase_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))

# UPEC + Intact NETs
UPEC_NETs_1_5_1h <- read.delim('data/1_5_1h_NETs_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_NETs_2_5_1h <- read.delim('data/2_5_1h_NETs_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_NETs_3_5_1h <- read.delim('data/3_5_1h_NETs_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_NETs_4_5_1h <- read.delim('data/4_5_1h_NETs_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_NETs_5_5_1h <- read.delim('data/5_5_1h_NETs_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))

# UPEC + Digested NETs
UPEC_DigestNETs_1_6_1h <- read.delim('data/1_6_1h_DigestedNETs_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_DigestNETs_2_6_1h <- read.delim('data/2_6_1h_DigestedNETs_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_DigestNETs_3_6_1h <- read.delim('data/3_6_1h_DigestedNETs_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_DigestNETs_4_6_1h <- read.delim('data/4_6_1h_DigestedNETs_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))
UPEC_DigestNETs_5_6_1h <- read.delim('data/5_6_1h_DigestedNETs_UPEC_featureCounts_Counts.tabular', header=T, col.names = c('GeneID', 'GeneCount'))

# Import gene lengths and assembly annotations
GeneID_Lenght <- read.delim('data/uti89_gene_feature_lenghts.tabular', header=T, col.names = c('GeneID', 'GeneLenght'))
UTI89_assembly <- rtracklayer::import('data/genomic.gtf') %>% base::as.data.frame()

# 3. Data clean

# Isolate count vectors by removing GeneID column
UPEC_1_2_1h <- UPEC_1_2_1h[,-1]; UPEC_2_2_1h <- UPEC_2_2_1h[,-1]; UPEC_3_2_1h <- UPEC_3_2_1h[,-1]; UPEC_4_2_1h <- UPEC_4_2_1h[,-1]; UPEC_5_2_1h <- UPEC_5_2_1h[,-1]
UPEC_PMA_1_3_1h <- UPEC_PMA_1_3_1h[,-1]; UPEC_PMA_2_3_1h <- UPEC_PMA_2_3_1h[,-1]; UPEC_PMA_3_3_1h <- UPEC_PMA_3_3_1h[,-1]; UPEC_PMA_4_3_1h <- UPEC_PMA_4_3_1h[,-1]; UPEC_PMA_5_3_1h <- UPEC_PMA_5_3_1h[,-1]
UPEC_DNase_1_4_1h <- UPEC_DNase_1_4_1h[,-1]; UPEC_DNase_2_4_1h <- UPEC_DNase_2_4_1h[,-1]; UPEC_DNase_3_4_1h <- UPEC_DNase_3_4_1h[,-1]; UPEC_DNase_4_4_1h <- UPEC_DNase_4_4_1h[,-1]; UPEC_DNase_5_4_1h <- UPEC_DNase_5_4_1h[,-1]
UPEC_NETs_1_5_1h <- UPEC_NETs_1_5_1h[,-1]; UPEC_NETs_2_5_1h <- UPEC_NETs_2_5_1h[,-1]; UPEC_NETs_3_5_1h <- UPEC_NETs_3_5_1h[,-1]; UPEC_NETs_4_5_1h <- UPEC_NETs_4_5_1h[,-1]; UPEC_NETs_5_5_1h <- UPEC_NETs_5_5_1h[,-1]
UPEC_DigestNETs_1_6_1h <- UPEC_DigestNETs_1_6_1h[,-1]; UPEC_DigestNETs_2_6_1h <- UPEC_DigestNETs_2_6_1h[,-1]; UPEC_DigestNETs_3_6_1h <- UPEC_DigestNETs_3_6_1h[,-1]; UPEC_DigestNETs_4_6_1h <- UPEC_DigestNETs_4_6_1h[,-1]; UPEC_DigestNETs_5_6_1h <- UPEC_DigestNETs_5_6_1h[,-1]

GeneID <- GeneID_Lenght[,-2]
GeneID_Lenght <- GeneID_Lenght[,-1]

# Construct primary Count Matrix
df_UPEC_NETs <- data.frame(UPEC_1_2_1h, UPEC_2_2_1h, UPEC_3_2_1h, UPEC_4_2_1h, UPEC_5_2_1h, 
                           UPEC_PMA_1_3_1h, UPEC_PMA_2_3_1h, UPEC_PMA_3_3_1h, UPEC_PMA_4_3_1h, UPEC_PMA_5_3_1h,
                           UPEC_DNase_1_4_1h, UPEC_DNase_2_4_1h, UPEC_DNase_3_4_1h, UPEC_DNase_4_4_1h, UPEC_DNase_5_4_1h,
                           UPEC_NETs_1_5_1h, UPEC_NETs_2_5_1h, UPEC_NETs_3_5_1h, UPEC_NETs_4_5_1h, UPEC_NETs_5_5_1h,
                           UPEC_DigestNETs_1_6_1h, UPEC_DigestNETs_2_6_1h, UPEC_DigestNETs_3_6_1h, UPEC_DigestNETs_4_6_1h, UPEC_DigestNETs_5_6_1h,
                           row.names = GeneID)

# Define Metadata
Conditions <- c(rep('Unstim', 5), rep('PMA', 5), rep('DNase', 5), rep('NETs', 5), rep('DigestNETs', 5))
Replicates <- rep(c('1', '2', '3', '4', '5'), 5)
md_UPEC_NETs <- data.frame(Replicates, Conditions)
rownames(md_UPEC_NETs) <- colnames(df_UPEC_NETs)

# Refine Genome Annotations
UTI89_assembly <- tibble(UTI89_assembly)
UTI89_genes <- UTI89_assembly %>% filter(type == 'gene') %>% select(-score, -phase, -transcript_id)
UTI89_cds_products <- UTI89_assembly %>% filter(type == 'CDS')

# Filter out ribosomal RNA
ribosomal_genes <- UTI89_genes[str_detect(UTI89_genes$gene_biotype, 'rRNA'),] %>% pull(gene_id)
df_UPEC_NETs <- df_UPEC_NETs[!(row.names(df_UPEC_NETs) %in% ribosomal_genes), ]

# 4. DESeq2 analysis

# DESeqDataSet
dds_UPEC_NETs <- DESeqDataSetFromMatrix(countData = df_UPEC_NETs,
                                        colData = md_UPEC_NETs,
                                        design = ~Conditions)

# Set 'Unstimulated' as the baseline
dds_UPEC_NETs$Conditions <- relevel(dds_UPEC_NETs$Conditions, ref = 'Unstim')

# Normalization
dds_UPEC_NETs <- estimateSizeFactors(dds_UPEC_NETs)
normalized_counts_UPEC_NETs <- counts(dds_UPEC_NETs, normalized=TRUE)

# VST for visualization
vsd_UPEC_NETs <- vst(dds_UPEC_NETs, blind=TRUE)
vsd_mat_UPEC_NETs <- assay(vsd_UPEC_NETs)

# 5. Visualization

# Sample Correlation Heatmap
vsd_cor_UPEC_NETs <- cor(vsd_mat_UPEC_NETs) 
pheatmap(vsd_cor_UPEC_NETs, annotation = select(md_UPEC_NETs, Conditions), main = "Sample Correlation")

# PCA Plot
pca_data <- plotPCA(vsd_UPEC_NETs, intgroup='Conditions', returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color=Conditions)) +
  geom_point(size=4) +
  theme_bw() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# 6. Differental expression

dds_UPEC_NETs <- DESeq(dds_UPEC_NETs)

# Extracting results - Intact NETs vs Unstimulated
de_NETs_Unstim <- results(dds_UPEC_NETs, contrast = c('Conditions', 'NETs', 'Unstim'),
                        alpha = 0.05, lfcThreshold = 0.5)

# Log2FoldChange Shrinkage
de_NETs_Unstim <- lfcShrink(dds_UPEC_NETs, coef = 4, res = de_NETs_Unstim, type = 'apeglm')

# Merge with gene annotations
UTI89_merge <- merge(UTI89_genes, UTI89_cds_products, by = 'gene_id', all.x = TRUE) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  select(-ends_with('.y'))

de_NETs_Unstim_all <- data.frame(de_NETs_Unstim) %>% 
  rownames_to_column(var = 'gene_id') %>%
  left_join(UTI89_merge[, c('gene_id', 'gene.x', 'product')], by = 'gene_id')

# Filter for significant results
de_NETs_Unstim_sig <- subset(de_NETs_Unstim_all, padj < 0.05) %>% arrange(padj)

# 7. Further visualization

# Volcano Plot with Custom Expressions
de_NETs_Unstim_all <- de_NETs_Unstim_all %>%
  mutate(Expression = case_when(
    log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
    TRUE ~ "Unchanged"))

ggplot(de_NETs_Unstim_all, aes(x = log2FoldChange, y = -log10(padj), color = Expression)) +
  geom_point() +
  xlab('log2 fold change') + ylab('-log10 adjusted p-value') +
  ggtitle('NETs vs Unstimulated') +
  theme_bw()

# 8. Export results
write_xlsx(de_NETs_Unstim_sig, 'results/NETs_vs_Unstimulated_sig.xlsx')
