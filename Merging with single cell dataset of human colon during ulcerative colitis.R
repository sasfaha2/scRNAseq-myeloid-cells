library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cluster)
library(reshape2)
library(harmony)

### Human Data: Harmony Batch Correction
uc <- RunHarmony(uc, "patient", plot_convergence = TRUE)

### Human Data: Gene Ortholog Conversion
mouse_homologs <- read.delim("../HOM_MouseHumanSequence.rpt.txt", header = TRUE, sep = "\t")
convert_mouse_to_human <- function(gene_list){
 output = c()
 for(gene in gene_list){
   class_key = (mouse_homologs %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
   if(!identical(class_key, integer(0)) ){
     human_genes = (mouse_homologs %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
     for(human_gene in human_genes){
       output = append(output,human_gene)
     }
   }
 }
 return (output)
>}

### Human and Mouse Integration
hsmm <- merge(uc, mimmune)
hsmm.list <- SplitObject(hsmm, split.by = "species")
# normalize and identify variable features for each dataset independently
hsmm.list <- lapply(X = hsmm.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nrow(hsmm), verbose = FALSE)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = hsmm.list)
immune.anchors <- FindIntegrationAnchors(object.list = hsmm.list, anchor.features = features)
# this command creates an 'integrated' data assay
hsmm <- IntegrateData(anchorset = immune.anchors, features.to.integrate = rownames(hsmm))
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(hsmm) <- "integrated"

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
