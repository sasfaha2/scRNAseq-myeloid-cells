library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cluster)
library(reshape2)
library(clusterProfiler)
library(org.Mm.eg.db)

### Differential Expression Analysis and Visualization
# Shared mono and macro markers
#MonoMac.markers <- FindMarkers(hsmm, ident.1 = c("Inflammatory Monocytes", "Inflammatory Macrophages"), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# log10 transform and make gene name column
log10(MonoMac.markers$p_val_adj)*-1 -> MonoMac.markers$p_val_adj_neglog10
rownames(MonoMac.markers) -> MonoMac.markers$Name
# label and color when log2 FC > 1 and p_val_adj_neglog10 > 2 (FDR > 1%)
ggplot(MonoMac.markers, aes(x = avg_log2FC, y = p_val_adj_neglog10)) +
  geom_point(color = dplyr::case_when
             (MonoMac.markers$avg_log2FC > 1 &
               MonoMac.markers$p_val_adj_neglog10 > 2 ~ "red",
               TRUE ~ "black"), size = 3, alpha = 0.8) +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = 1) +
  ggtitle("Inflammatory Monocyte & Macrophage Markers") + 
  ylab("-log10 Adjusted p Value") +
  xlab("Average log2FC") +
  theme(plot.title = element_text(size=25, hjust = 0.5), text = element_text(size = 20), axis.text = element_text(size = 15)) +
  geom_text_repel(data = . %>% mutate(label = ifelse(rownames(MonoMac.markers) %in% rownames(MonoMac.markers[(MonoMac.markers$avg_log2FC > 1) &
                                    MonoMac.markers$p_val_adj_neglog10 > 2,]), Name, "")), aes(avg_log2FC, p_val_adj_neglog10, label = label), max.overlaps = 50, size = 5, box.padding = .5, force = 2)

### Representative GO enrichment Analysis and Visualization
# convert symbol to ensembl IDs
ensembl <- bitr(rownames(hsmm), fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = "org.Mm.eg.db")
# identify top genes by condition - only inflammatory cells
# get DEG by comparing inflammatory conditions to health independently
hsmm <- SetIdent(hsmm, value = "condition")
SetIdent(hsmm, value = "condition") %>%
  subset(hsmm, idents = c("UC", "Healthy")) %>%
   FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) -> condition.markers
SetIdent(hsmm, value = "condition") %>%
  subset(hsmm, idents = c("Control", "DSS", "Oxa")) %>% 
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% 
  rbind(., condition.markers) -> condition.markers
# filter for top genes by group
condition_ids <- condition.markers %>% group_by(cluster) %>% 
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  dplyr::filter(row_number() <= 50 & cluster %in% c("UC", "DSS", "Oxa"))
p <- as.character(unique(condition_ids$cluster))
condition_ids$cluster <- factor(condition_ids$cluster)
# make list of genes and their ensemble IDs
condition_ids <- sapply(levels(as.factor(condition_ids$cluster)), function(x) {ensembl[match((condition_ids[condition_ids$cluster==x,] %>% .$gene), ensembl$SYMBOL), "ENTREZID"]})
condition_ids <- split(condition_ids, rep(1:ncol(condition_ids), each = nrow(condition_ids)))
names(condition_ids) <- p
condition_mods <- compareCluster(condition_ids, fun = "enrichGO", OrgDb = org.Mm.eg.db, ont = "BP", minGSSize = 50, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
p1 <- dotplot(condition_mods, showCategory = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15), text = element_text(size = 15), axis.title.x = element_text(size = 15))
p1

### Human and Mouse Cluster Similarity Analysis
# Identify human cluster marker genes
SetIdent(hsmm, value = "condition") %>%
  subset(., idents = c("Inflamed", "Non-inflamed", "Healthy")) %>%
  SetIdent(value = "ids") %>% FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) -condition.markers
# Identify mouse cluster marker genes
SetIdent(hsmm, value = "condition") %>%
  subset(., idents = c("control", "dss", "oxa")) %>%
  SetIdent(value = "cell_ids") %>% FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% rbind(., condition.markers) -condition.markers
# 78 is minimum number of significant markers for a cluster
condition_sim <- condition.markers %>% group_by(cluster) %>% dplyr::arrange(desc(avg_log2FC), .by_group = TRUE) %>% dplyr::filter(row_number() <= 78)
sim <- matrix(NA, nrow = 13, ncol = 13)
colnames(sim) <- levels(condition_sim$cluster)
rownames(sim) <- levels(condition_sim$cluster)
for (i in 1:13){
  for (l in 1:13){
# determine shared cluster markers between all clusters out of 78
    sim[i,l] <- length(which((condition_sim %>% dplyr::filter(cluster == rownames(sim)[i]) %>% .$gene) %in% (condition_sim %>% dplyr::filter(cluster == rownames(sim)[l]) %>% .$gene)) == TRUE)/78
  }
}
heatmap(sqrt(sim), scale = "none", cexRow = 2, cexCol = 2, margins = c(25,25))
