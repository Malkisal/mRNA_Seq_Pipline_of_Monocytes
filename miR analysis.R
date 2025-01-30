# My Script 

#activate required libraries
library(tidyverse)
library(DESeq2)
library(readxl)
library(ggplot2)
library(ggrepel)
library(ggpubr)

#read in the gene counts and the sample information files
coldata <- read_excel("coldata miRs.xlsx")
gene_count <- read_excel("Readcount Only.xlsx")


# convert the count data into a data frame.
gene_count <- as.data.frame(gene_count)


# need first to remove dublicate from the gene name columns, as it has to be unique.
gene_count2 <- gene_count [!duplicated(gene_count$sRNA.readcount),]


# convert the gene names column to row names and deltete the column woth gene name
gene_count3 <- gene_count2[,-1]
# take gene name from gene_count2 and add it into gene_coun 3 in the row name
rownames(gene_count3) <- gene_count2 [,1]

# finally we can run the dds for gene_count3
dds_miRs <- DESeqDataSetFromMatrix(
  countData = gene_count3,
  colData = coldata,
  design = ~ Condition)
# Filtration of count gene with very low counts.

low_expressed <- counts(dds_miRs) >10
keep <- rowSums(low_expressed) == ncol(dds_miRs)

Keep_dds_miRs <- rowSums(counts(dds_miRs)) >= 10
dds_miRs <- dds_miRs [keep,]

# PCA for gene counts
# first we need to creeate a transform tool like the Log to reduce the variants


# Error happend! To resolve it:

# Check the number of rows in your dataset
nrow_dds_miRs <- nrow(dds_miRs)

# Use varianceStabilizingTransformation directly
vsd <- varianceStabilizingTransformation(dds_miRs)


# how to see things in your files
unique(dds_miRs$Condition)


# Plot PCA:
plotPCA (vsd,intgroup = "Participant",ntop=500)


# convert Participant to Factor:
dds_miRs$ Participant<- factor(dds_miRs$Participant)


# to run the differential analysis:
#model comparing condition without adjustment
de_miRs <- DESeq(dds_miRs)
res_ctrl_vs_nLDL <- results(de_miRs, 
                            contrast = c("Condition", "CTRL","nLDL"),tidy=TRUE) %>% arrange(pvalue)

# to see the rsults file:

head (res_ctrl_vs_nLDL)

# call result:

res_ctrl_vs_oxLDL<- results(de_miRs, contrast = c("Condition", "CTRL", "oxLDL"),tidy=TRUE) %>%
  arrange(pvalue)

res_ctrl_vs_oxLDL

head (res_ctrl_vs_oxLDL)

res_nLDL_vs_oxLDL<- results(de_miRs, contrast = c("Condition", "nLDL", "oxLDL"),tidy=TRUE) %>%
  arrange(pvalue)

head (res_nLDL_vs_oxLDL)

resultsNames(de_miRs)

# to export table:

write.csv(res_ctrl_vs_nLDL,'CTRL vs nLDL_NMVs No correction.csv')
write.csv(res_ctrl_vs_oxLDL,'CTRL vs oxLDL_NMVs No correction.csv')
write.csv(res_nLDL_vs_oxLDL,'nLDL_NMVs vs oxLDL_NMVs No correction.csv')

#model comparing condition with adjustment for participant
dds_miRs2 <- dds_miRs
design(dds_miRs2) <- ~ Condition + Participant

# to see the new design of the new file (2)
design(dds_miRs2)

# Convert character variables to factors
# For example, if "condition" is a character variable
dds_miRs$ Gender<- factor(dds_miRs$Gender)

#model comparing condition with adjustment for Gender
dds_miRs3 <- dds_miRs
design(dds_miRs3) <- ~ Condition + Participant



# to rerun the deseq after changing the design (3):
de_miRs3 <- DESeq(dds_miRs3)
design(de_miRs3)

# call results:
res_ctrl_vs_nLDL3 <- results(de_miRs3, contrast = c("Condition", "CTRL", "nLDL"),tidy=TRUE) %>%
  arrange(pvalue)

# to see the results file:
head(res_ctrl_vs_nLDL3)

res_ctrl_vs_nLDL3

# to create file for sig. only:
sig_res_ctrl_vs_nLDL3 <- subset(res_ctrl_vs_nLDL3, pvalue <0.05)

write.csv(sig_res_ctrl_vs_nLDL3, 'sig <0.1 CTRL vs nLDL_Participant Correction.csv')

# call result CTRL vs oxLDL:
res_ctrl_vs_oxLDL3 <- results(de_miRs3, contrast = c("Condition", "CTRL", "oxLDL"),tidy=TRUE) %>%
  arrange(pvalue)

head(res_ctrl_vs_oxLDL3)

# to create file for sig. only:
sig_res_ctrl_vs_oxLDL3 <- subset(res_ctrl_vs_oxLDL3, pvalue <0.1)

write.csv(sig_res_ctrl_vs_oxLDL3, 'sig<0.1 CTRL vs oxLDL_ Participant Correction.csv')

# call result nLDL vs oxLDL:
res_nLDL_vs_oxLDL3<- results(de_miRs3, contrast = c("Condition", "nLDL", "oxLDL"),tidy=TRUE) %>%
  arrange(pvalue)

head(res_nLDL_vs_oxLDL3)


# to create file for sig. only:
sig_res_nLDL_vs_oxLDL3 <- subset(res_nLDL_vs_oxLDL3, pvalue <0.05)

write.csv(sig_res_nLDL_vs_oxLDL3, 'sig<0.1 nLDL vs oxLDL_Participant Correction.csv')

#volcano

#create a new coloumn called diffexpressed to categorise the expression to up, down and no

#ctrl vs nLDL
res_ctrl_vs_nLDL3$diffexpressed <- "NO"
res_ctrl_vs_nLDL3$diffexpressed[res_ctrl_vs_nLDL3$log2FoldChange > 0.0 & res_ctrl_vs_nLDL3$pvalue < 0.05] <- "UP"
res_ctrl_vs_nLDL3$diffexpressed[res_ctrl_vs_nLDL3$log2FoldChange < -0.0 & res_ctrl_vs_nLDL3$pvalue < 0.05] <- "DOWN"
head(res_ctrl_vs_nLDL3)
unique(res_ctrl_vs_nLDL3$diffexpressed)
table(res_ctrl_vs_nLDL3$diffexpressed) 

#ctrl vs oxLDL
res_ctrl_vs_oxLDL3$diffexpressed <- "NO"
res_ctrl_vs_oxLDL3$diffexpressed[res_ctrl_vs_oxLDL3$log2FoldChange > 0.0 & res_ctrl_vs_oxLDL3$pvalue < 0.05] <- "UP"
res_ctrl_vs_oxLDL3$diffexpressed[res_ctrl_vs_oxLDL3$log2FoldChange < -0.0 & res_ctrl_vs_oxLDL3$pvalue < 0.05] <- "DOWN"
head(res_ctrl_vs_oxLDL3)
unique(res_ctrl_vs_oxLDL3$diffexpressed)
table(res_ctrl_vs_oxLDL3$diffexpressed)

#nLDL vs oxLDL
res_nLDL_vs_oxLDL3$diffexpressed <- "NO"
res_nLDL_vs_oxLDL3$diffexpressed[res_nLDL_vs_oxLDL3$log2FoldChange > 0.0 & res_nLDL_vs_oxLDL3$pvalue < 0.05] <- "UP"
res_nLDL_vs_oxLDL3$diffexpressed[res_nLDL_vs_oxLDL3$log2FoldChange < -0.0 & res_nLDL_vs_oxLDL3$pvalue < 0.05] <- "DOWN"
head(res_nLDL_vs_oxLDL3)
unique(res_nLDL_vs_oxLDL3$diffexpressed)
table(res_nLDL_vs_oxLDL3$diffexpressed)



#Create a new column "Top" that includes the names of the top differentially expressed genes

#ctrl vs nLDL
res_ctrl_vs_nLDL3$Top <- ifelse(res_ctrl_vs_nLDL3$row %in% head(arrange(res_ctrl_vs_nLDL3,pvalue),20)$row & res_ctrl_vs_nLDL3$diffexpressed != "NO", res_ctrl_vs_nLDL3$row, NA)
head(res_ctrl_vs_nLDL3)

#ctrl vs oxLDL
res_ctrl_vs_oxLDL3$Top <- ifelse(res_ctrl_vs_oxLDL3$row %in% head(arrange(res_ctrl_vs_oxLDL3,pvalue),20)$row & res_ctrl_vs_oxLDL3$diffexpressed != "NO", res_ctrl_vs_oxLDL3$row, NA)
head(res_ctrl_vs_oxLDL3)

#nLDL vs oxLDL
res_nLDL_vs_oxLDL3$Top <- ifelse(res_nLDL_vs_oxLDL3$row %in% head(arrange(res_nLDL_vs_oxLDL3,pvalue),20)$row & res_nLDL_vs_oxLDL3$diffexpressed != "NO", res_nLDL_vs_oxLDL3$row, NA)
head(res_nLDL_vs_oxLDL3)

#plot

table(res_ctrl_vs_nLDL3$diffexpressed)
ctrl_vs_nLDL3 <- ggplot(res_ctrl_vs_nLDL3, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label=Top))+
  geom_vline(xintercept = 0, col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("aquamarine4", "grey", "sienna3"),
                     labels = c("DOWN", "NO", "UP")) +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))+
  ggtitle('Ctrl vs nLDL-NMVs') +geom_text_repel(max.overlaps = Inf, size=2.5, fontface="bold")+
  theme(legend.title = element_blank(), plot.title.position = "right", plot.title = element_text(hjust = 0.5))+theme_classic()


table(res_ctrl_vs_oxLDL3$diffexpressed)
ctrl_vs_oxLDL3 <- ggplot(res_ctrl_vs_oxLDL3, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label=Top))+
  geom_vline(xintercept = 0, col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("aquamarine4", "grey", "sienna3"),
                     labels = c("DOWN", "NO", "UP")) +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))+
  ggtitle('Ctrl vs oxLDL-NMVs') +geom_text_repel(max.overlaps = Inf, size=2.5, fontface="bold")+
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+theme_classic()

table(res_nLDL_vs_oxLDL3$diffexpressed)
nLDL_vs_oxLDL3 <- ggplot(res_nLDL_vs_oxLDL3, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label=Top))+
  geom_vline(xintercept = 0, col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("aquamarine4", "grey", "sienna3"),
                     labels = c("DOWN", "NO", "UP")) +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))+
  ggtitle('nLDL-NMVs vs oxLDL-NMVs') +geom_text_repel(max.overlaps = Inf, size=2.5, fontface="bold")+
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+theme_classic()

ggarrange(ctrl_vs_nLDL3, ctrl_vs_oxLDL3,nLDL_vs_oxLDL3, nrow=1, common.legend = TRUE) %>% 
  annotate_figure(top = text_grob("Model: Condition with adjustment of Participant",
                                  
                                  face = "bold", size = 14))
# to Save it:
ggsave("Volcano plot condition Plus Part (miRs).tiff", width = 15, height=6)

# Heatmap:

install.packages("pheatmap")
library(pheatmap)

gene_expression <- read.csv("Readcount Only.csv", row.names = 1)  

# Assuming gene names are in the first column
gene_expression <- t(gene_expression)

gene_expression <- log2(gene_expression + 1)  # Adding 1 to avoid log(0)


pheatmap(gene_expression,
         cluster_rows = TRUE,    # Cluster rows (genes)
         cluster_cols = TRUE,    # Cluster columns (samples)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette
         main = "Gene Expression Heatmap",  # Title of the heatmap
         fontsize = 8            # Font size
)


#Enrichment analysis

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")



library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(tidyverse)
library(ggpubr)



genes_to_test <- sig_res_ctrl_vs_nLDL3[,1]
head(genes_to_test,20)

GO_BP <- enrichGO (gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_CC <- enrichGO (gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC")
GO_MF <- enrichGO (gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF")

GO_BP_df <- as.data.frame(GO_BP)
GO_CC_df <- as.data.frame(GO_CC)
GO_MF_df <- as.data.frame(GO_MF)

# test my filtered genes.
genes_to_test <- res_ctrl_vs_nLDL3[,1]
head(genes_to_test,20)

GO_BP <- enrichGO (gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_CC <- enrichGO (gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC")
GO_MF <- enrichGO (gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF")

GO_BP_df <- as.data.frame(GO_BP)
GO_CC_df <- as.data.frame(GO_CC)
GO_MF_df <- as.data.frame(GO_MF)

# plot results:

BP_plot_CTRL_vs_nLDL <- head(GO_BP_CTRL_vs_nLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Biological Process")
CC_plot_CTRL_vs_nLDL <- head(GO_CC_CTRL_vs_nLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Cellular Component")
MF_plot_CTRL_vs_nLDL <- head(GO_MF_CTRL_vs_nLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Molecular Function")

ggarrange(BP_plot, CC_plot,MF_plot, ncol=1, common.legend = TRUE) %>%
  annotate_figure(top = text_grob("Control vs nLDL MMVs" ,
                                  face = "bold", size = 14))
#my dot plots of GO


# Plot Biological Process (BP)
dotplot(GO_BP, title = "Biological Process", showCategory = 20)

# Plot Cellular Component (CC)
dotplot(GO_CC, title = "Cellular Component", showCategory = 20)

# Plot Molecular Function (MF)
dotplot(GO_MF, title = "Molecular Function", showCategory = 20)

# Combine the plots
install.packages("clusterProfiler")
install.packages("cowplot")

library(clusterProfiler)
library(cowplot)
combined_plot <- plot_grid(BP_plot, CC_plot, MF_plot, ncol = 1)

# Display the combined plot
print(combined_plot)




#KEGG

entrez_genes <- mapIds(org.Hs.eg.db, genes_to_test, 'ENTREZID', 'SYMBOL')
kegg <- enrichKEGG(gene = genes_to_test, organism = 'hsa', keyType="kegg", pvalueCutoff = 0.05)
kegg_df <- as.data.frame(kegg)

#Plot KEGG
dotplot(kegg_3, showCategory=20)

#====================================================================================

# CTRL vs oxLDL Erich.

genes_ctrl_vs_oxLDL3 <- res_ctrl_vs_oxLDL3[,1]
head(genes_ctrl_vs_oxLDL3,20)

GO_BP_CTRL_vs_oxLDL <- enrichGO (gene = genes_ctrl_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO_CC_CTRL_vs_oxLDL <- enrichGO (gene = genes_ctrl_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC")
GO_MF_CTRL_vs_oxLDL <- enrichGO (gene = genes_ctrl_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF")

GO_BP_CTRL_vs_oxLDL_df <- as.data.frame(GO_BP_CTRL_vs_oxLDL)
GO_CC_CTRL_vs_oxLDL_df <- as.data.frame(GO_CC_CTRL_vs_oxLDL)
GO_MF_CTRL_vs_oxLDL_df <- as.data.frame(GO_MF_CTRL_vs_oxLDL)


BP_plot_CTRL_vs_oxLDL <- head(GO_BP_CTRL_vs_oxLDL_df,20) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Biological Process")
CC_plot_CTRL_vs_oxLDL <- head(GO_CC_CTRL_vs_oxLDL_df,20) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Cellular Component")
MF_plot_CTRL_vs_oxLDL <- head(GO_MF_CTRL_vs_oxLDL_df,20) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Molecular Function")

ggarrange(BP_plot, CC_plot,MF_plot, ncol=1, common.legend = TRUE) %>%
  annotate_figure(top = text_grob("Control vs nLDL MMVs" ,
                                  face = "bold", size = 14))
#my dot plots of GO


# Plot Biological Process (BP)
dotplot(GO_BP, title = "Biological Process", showCategory = 20)

# Plot Cellular Component (CC)
dotplot(GO_CC, title = "Cellular Component", showCategory = 20)

# Plot Molecular Function (MF)
dotplot(GO_MF, title = "Molecular Function", showCategory = 20)

# Combine the plots
install.packages("clusterProfiler")
install.packages("cowplot")

library(clusterProfiler)
library(cowplot)
combined_plot <- plot_grid(BP_plot, CC_plot, MF_plot, ncol = 1)

# Display the combined plot
print(combined_plot)




#KEGG

entrez_genes <- mapIds(org.Hs.eg.db, genes_to_test, 'ENTREZID', 'SYMBOL')
kegg_CTRL_vs_oxLDL <- enrichKEGG(gene = entrez_genes, organism = 'hsa', keyType="kegg", pvalueCutoff = 0.05)
kegg_kegg_CTRL_vs_oxLDL_df <- as.data.frame(kegg)

#Plot KEGG
dotplot(kegg_3, showCategory=20)


#================================================================================

# nLDL vs oxLDL Enrich.

genes_nLDL_vs_oxLDL3 <- res_nLDL_vs_oxLDL3[,1]
head(genes_nLDL_vs_oxLDL3,20)

GO_BP_nLDL_vs_oxLDL <- enrichGO (gene = genes_nLDL_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_CC_nLDL_vs_oxLDL <- enrichGO (gene = genes_nLDL_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC")
GO_MF_nLDL_vs_oxLDL <- enrichGO (gene = genes_nLDL_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF")

GO_BP_nLDL_vs_oxLDL_df <- as.data.frame(GO_BP_nLDL_vs_oxLDL)
GO_CC_nLDL_vs_oxLDL_df <- as.data.frame(GO_CC_nLDL_vs_oxLDL)
GO_MF_nLDL_vs_oxLDL_df <- as.data.frame(GO_MF_nLDL_vs_oxLDL)


BP_plot_nLDL_vs_oxLDL <- head(GO_BP_nLDL_vs_oxLDL_df,20) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Biological Process")
CC_plot_nLDL_vs_oxLDL <- head(GO_CC_nLDL_vs_oxLDL_df,20) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Cellular Component")
MF_plot_nLDL_vs_oxLDL <- head(GO_MF_nLDL_vs_oxLDL_df,20) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Molecular Function")

ggarrange(BP_plot, CC_plot,MF_plot, ncol=1, common.legend = TRUE) %>%
  annotate_figure(top = text_grob("Control vs nLDL MMVs" ,
                                  face = "bold", size = 14))
#my dot plots of GO


# Plot Biological Process (BP)
dotplot(GO_BP, title = "Biological Process", showCategory = 20)

# Plot Cellular Component (CC)
dotplot(GO_CC, title = "Cellular Component", showCategory = 20)

# Plot Molecular Function (MF)
dotplot(GO_MF, title = "Molecular Function", showCategory = 20)

# Combine the plots
install.packages("clusterProfiler")
install.packages("cowplot")

library(clusterProfiler)
library(cowplot)
combined_plot <- plot_grid(BP_plot, CC_plot, MF_plot, ncol = 1)

# Display the combined plot
print(combined_plot)




#KEGG

entrez_genes <- mapIds(org.Hs.eg.db, genes_nLDL_vs_oxLDL3, 'ENTREZID', 'SYMBOL')
kegg_nLDL_vs_oxLDL <- enrichKEGG(gene = entrez_genes, organism = 'hsa', keyType="kegg", pvalueCutoff = 0.05)
kegg_nLDL_vs_oxLDL_df3 <- as.data.frame(kegg_nLDL_vs_oxLDL)

#Plot KEGG
dotplot(kegg_nLDL_vs_oxLDL, showCategory=20)


