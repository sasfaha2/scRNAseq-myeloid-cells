library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cluster)
library(reshape2)
library(SingleR)
library(celldex)
library(SeuratWrappers)
library(monocle3)
library(clusterProfiler)
library(org.Mm.eg.db)

### Annotation
## SingleR Annotation
# Import reference dataset
ref.data <- ImmGenData()
pred.immune <- SingleR(test=mimmune@assays$RNA@data, ref = ref.data, assay.type.test = 1, labels = ref.data$label.main)

## Hand picked markers
mimmune$cell_class <- mimmune$cell_ids
levels(mimmune$cell_class) <- c("Neutrophils", "Neutrophils", "T Cells", "Neutrophils", "Monocytes & Macrophages", "Monocytes & Macrophages", "T Cells", "Monocytes & Macrophages", "B Cells", "Monocytes & Macrophages", "Neutrophils", "T Cells", "Monocytes & Macrophages")
SetIdent(mimmune, value = "cell_class") %>%
  DotPlot(., features = c("Ly6g", "Cxcr2", "Hdc", "Cd3d", "Cd3g", "Trbc2", "Ccr2", "Csf1r", "Lyz2", "Igkc", "Cd79a", "Ms4a1"), dot.scale = 8, scale.max = 100, scale.min = 0) + RotatedAxis() + scale_color_viridis_c()

## AB vs GD TCR Scores
# No joining or diversity genes; one j gene for delta
ab.markers <- c(
  rownames(mimmune)[which(substring(rownames(mimmune),1,4) %in% "Trav")],
  "Trac",
  rownames(mimmune)[which(substring(rownames(mimmune),1,3) %in% "Trb")])
gb.markers <- c(
  rownames(mimmune)[which(substring(rownames(mimmune),1,4) %in% "Tcrg")],
  "Trdv1", "Trdv2-2", "Trdv4", "Trdj1", "Trdc", "Trdv5")
# assign count data for markers to new df
ab.markers <- as.data.frame(t(cbind(mimmune@assays$RNA@scale.data[ab.markers,])))
gb.markers <- as.data.frame(t(cbind(mimmune@assays$RNA@scale.data[gb.markers,])))
# each marker is given an equal weight for each cell
# for loop: (value-min)/(max-min)
# marker ranks are summed across genes
for(i in 1:ncol(ab.markers)){
  #min is negative value, adding the abs(min)
  #lowest value is now 0
  ab.markers[,i] <- (ab.markers[,i] + abs(min(ab.markers[,i])))
  #(value-min)/(max-min); now value/max
  ab.markers[,i] <- (ab.markers[,i]/(max(ab.markers[,i])))
}
for(i in 1:ncol(gb.markers)){
  gb.markers[,i] <- (gb.markers[,i] + abs(min(gb.markers[,i])))
  gb.markers[,i] <- (gb.markers[,i]/(max(gb.markers[,i])))
}
# sum marker scores
ab.markers <- as.data.frame(rowSums(ab.markers))
gb.markers <- as.data.frame(rowSums(gb.markers))
colnames(ab.markers) <- "ab_score"
colnames(gb.markers) <- "gd_score"
# append to meta data
mimmune$ab_score <- ab.markers
mimmune$gd_score <- gb.markers

### Differential Expression Analysis and Visualization
# DSS vs Oxa
DSSOxa.markers <- SetIdent(mimmune, value = "orig.ident") %>%
  FindMarkers(., ident.1 = "dss", ident.2 = "oxa", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
# log10 transform and make gene name column
log10(DSSOxa.markers$p_val_adj)*-1 -> DSSOxa.markers$p_val_adj_neglog10
rownames(DSSOxa.markers) -> DSSOxa.markers$Name
# label and color when log2 FC > 1 and p_val_adj_neglog10 > 2 (FDR 1%)
ggplot(DSSOxa.markers, aes(x = avg_log2FC, y = p_val_adj_neglog10)) +
  geom_point(color = dplyr::case_when
             ((DSSOxa.markers$avg_log2FC > 1 |
                 DSSOxa.markers$avg_log2FC < -1) &
                 DSSOxa.markers$p_val_adj_neglog10 > 2 ~ "red",
               TRUE ~ "black"), size = 3, alpha = 0.8) +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  ggtitle("Oxa vs. DSS Differential Gene Expression") + 
  ylab("-log10 Adjusted p Value") +
  xlab("Average log2FC") +
  theme(plot.title = element_text(size=25, hjust = 0.5), text = element_text(size = 20), axis.text = element_text(size = 15)) +
  geom_text_repel(data = subset(DSSOxa.markers,
                                (DSSOxa.markers$avg_log2FC > 1 |
                                   DSSOxa.markers$avg_log2FC < -1) &
                                  DSSOxa.markers$p_val_adj_neglog10 > 2), aes(avg_log2FC, p_val_adj_neglog10, label = Name), 
                  max.overlaps = 50, size = 5, box.padding = .5, force = 2)

### Representative GO enrichment Analysis and Visualization
# convert symbol to ensembl IDs
ensembl <- bitr(rownames(mimmune), fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = "org.Mm.eg.db")
# filter for top genes by group
cluster_ids <- mimmune.markers %>%
  group_by(cluster) %>% 
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  dplyr::filter(row_number() <= 50)
p <- as.character(unique(cluster_ids$cluster))
# make list of genes and their ensemble IDs
cluster_ids <- sapply(levels(as.factor(cluster_ids$cluster)), function(x) {ensembl[match((cluster_ids[cluster_ids$cluster==x,] %>% .$gene), ensembl$SYMBOL), "ENTREZID"]})
cluster_ids <- split(cluster_ids, rep(1:ncol(cluster_ids), each = nrow(cluster_ids)))
names(cluster_ids) <- p
condition_mods <- compareCluster(cluster_ids, fun = "enrichGO", OrgDb = org.Mm.eg.db, ont = "BP", minGSSize = 50, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
p1 <- dotplot(condition_mods, showCategory = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15), text = element_text(size = 15), axis.title.x = element_text(size = 15))
p1

### Pseudotime Trajectory Analysis of Neutrophils
# Generate Monocle3 object and cluster
SeuratWrappers::as.cell_data_set(subset(mimmune, ident = c("Immature Neutrophils", "Mature Neutrophils", "Ly6gLo Neutrophils"))) -mimmune.cds
mimmune.cds <- cluster_cells(cds = mimmune.cds, reduction_method = "UMAP")
mimmune.cds <- learn_graph(mimmune.cds, use_partition = TRUE)
# CRITICAL: ident for root cells needs to be correct
root.mimmune <- WhichCells(mimmune, ident = c("Immature Neutrophils", "Mature Neutrophils", "Ly6gLo Neutrophils"), expression = Camp 2)
mimmune.cds <- order_cells(mimmune.cds, reduction_method = "UMAP", root_cells = root.mimmune)

### Ly6c2 Correlation Analysis
matrix <- subset(mimmune, ident = c("Itgax+ Macrophages", "Itgax- Macrophages", "Ly6cHi Monocytes", "Resident Monocytes"))
matrix <- as.matrix(matrix@assays$RNA@scale.data)[c(VariableFeatures(mimmune), "Sell"),]
gene <- as.numeric(matrix["Ly6c2",])
correlations <- as.data.frame(apply(matrix,1,function(x){cor(gene,x)}))
cor <- as.data.frame(correlations[order(correlations$correlations, decreasing = T),])
colnames(cor) <- "correlations"
# sort by strength of correlation
rownames(cor) <- rownames(correlations)[order(correlations$correlations, decreasing = T)]
cor$rownames <- rownames(cor)
cor$rownames <- factor(cor$rownames, levels = rownames(cor))
