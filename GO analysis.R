# Project: Gene Ontology (GO) Enrichment Analysis
# Description: Functional enrichment analysis of DEGs from NETs and Mutant 
#              experiments using clusterProfiler and E. coli K12 database.
# Author: Dora Cerina

# Setup and libraries

# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.EcK12.eg.db"))

library(ggplot2)
library(GOplot)
library(clusterProfiler)
library(org.EcK12.eg.db)
library(dplyr)
library(gridExtra)
library(readxl)

# Check available columns in the E. coli K12 database
columns(org.EcK12.eg.db) 

# GO analysis: NETs vs Unstimulated

Unstimulated_vs_NETs_sig <- read_excel("results/Unstimulated_vs_NETs_significant.xlsx")

upregulated_NETs <- subset(Unstimulated_vs_NETs_sig, log2FoldChange > 0)
downregulated_NETs <- subset(Unstimulated_vs_NETs_sig, log2FoldChange < 0)

go_up_NETs <- upregulated_NETs$gene.x
go_down_NETs <- downregulated_NETs$gene.x

go_up_enrich <- enrichGO(gene = go_up_NETs,
                          OrgDb = org.EcK12.eg.db,  # Specify the organism database
                          keyType = "ALIAS", # Gene ID type
                          ont = "BP",           # Ontology type: "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                          pvalueCutoff = 0.05,
                          readable = TRUE)

go_down_enrich <- enrichGO(gene = go_down_NETs,
                         OrgDb = org.EcK12.eg.db,  # Specify the organism database
                         keyType = "ALIAS", # Gene ID type
                         ont = "BP",           # Ontology type: "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                         pvalueCutoff = 0.05,
                         readable = TRUE)

go_up_enrich_simp <- simplify(go_up_enrich, cutoff = 0.6, by='p.adjust', select_fun=min)
go_down_enrich_simp <- simplify(go_down_enrich, cutoff = 0.6, by='p.adjust', select_fun=min)

df_go_up_simp <- go_up_enrich_simp@result
df_go_down_simp <- go_down_enrich_simp@result

df_go_up_simp_unique <- df_go_up_simp %>%
  distinct(geneID, .keep_all = TRUE)

df_go_down_simp_unique <- df_go_down_simp %>%
  distinct(geneID, .keep_all = TRUE)

#Log transform the p values
df_go_up_simp_unique <- df_go_up_simp_unique %>%
  mutate(log_pvalue = -log(p.adjust, base = 10))

df_go_down_simp_unique <- df_go_down_simp_unique %>%
  mutate(log_pvalue = -log(p.adjust, base = 10))

#Adding coloum with the Up or down expression
df_go_up_simp_unique$Expression <- c('Up-regulated')
df_go_down_simp_unique$Expression <- c('Down-regulated')

#Top10 terms
df_go_up_simp_unique <- df_go_up_simp_unique[1:10,]
df_go_down_simp_unique <- df_go_down_simp_unique[1:10,]

cols <- c("Up-regulated" = "firebrick", "Down-regulated" = "steelblue") 

ggplot(data = df_go_up_simp_unique,
       aes(x=log_pvalue, 
           y=reorder(Description, log_pvalue),
           fill = Expression)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(0, 10)) +
  theme_bw() +
  theme(legend.position = 'top', # Sets the position of the legend.
        plot.title = element_text (size = rel(1), hjust = 1),
        axis.title = element_text (size = rel(1))) +
  labs(x = '-log10 adjusted p-value', y = 'GO-term') 



# GO analysis: WT vs DEL32 (DanI mutant)

WTvsDEL32_sig <-  read_excel("results/WT_vs_DEL32_sig.xlsx")

up_WTvsDEL32 <- subset(WTvsDEL32_sig, log2FoldChange > 0)
down_WTvsDEL32 <- subset(WTvsDEL32_sig, log2FoldChange < 0)

go_up_DEL32<- up_WTvsDEL32$gene.x
go_down_DEL32 <- down_WTvsDEL32$gene.x

go_up_enrich <- enrichGO(gene = go_up_DEL32,
                         OrgDb = org.EcK12.eg.db,  # Specify the organism database
                         keyType = "ALIAS", # Gene ID type
                         ont = "BP",           # Ontology type: "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                         pvalueCutoff = 0.05,
                         readable = TRUE)

go_down_enrich <- enrichGO(gene = go_down_DEL32,
                           OrgDb = org.EcK12.eg.db,  # Specify the organism database
                           keyType = "ALIAS", # Gene ID type
                           ont = "BP",           # Ontology type: "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                           pvalueCutoff = 0.05,
                           readable = TRUE)

go_up_enrich_simp <- simplify(go_up_enrich, cutoff = 0.6, by='p.adjust', select_fun=min)
go_down_enrich_simp <- simplify(go_down_enrich, cutoff = 0.6, by='p.adjust', select_fun=min)

df_go_up_simp <- go_up_enrich_simp@result
df_go_down_simp <- go_down_enrich_simp@result

df_go_up_simp_unique <- df_go_up_simp %>%
  distinct(geneID, .keep_all = TRUE)

df_go_down_simp_unique <- df_go_down_simp %>%
  distinct(geneID, .keep_all = TRUE)

#Log transform the p values
df_go_up_simp_unique <- df_go_up_simp_unique %>%
  mutate(log_pvalue = -log(p.adjust, base = 10))

df_go_down_simp_unique <- df_go_down_simp_unique %>%
  mutate(log_pvalue = -log(p.adjust, base = 10))

#Adding coloum with the Up or down expression
df_go_up_simp_unique$Expression <- c('Up-regulated')
df_go_down_simp_unique$Expression <- c('Down-regulated')

#Top10 terms
df_go_up_simp_unique <- df_go_up_simp_unique[1:10,]
df_go_down_simp_unique <- df_go_down_simp_unique[1:10,]

cols <- c("Up-regulated" = "firebrick", "Down-regulated" = "steelblue") 

ggplot(data = df_go_up_simp_unique,
       aes(x=log_pvalue, 
           y=reorder(Description, log_pvalue),
           fill = Expression)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(0, 15.5)) +
  theme_bw() +
  theme(legend.position = 'top', # Sets the position of the legend.
        plot.title = element_text (size = rel(1), hjust = 1),
        axis.title = element_text (size = rel(1))) +
  labs(x = '-log10 adjusted p-value', y = 'GO-term') 
