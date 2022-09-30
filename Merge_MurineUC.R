library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cluster)
library(reshape2)
library(Matrix)
library(harmony)

### Human Data: Import
uc.sparse <- readMM(file = "./SCP259/gene_sorted-Imm.matrix.mtx")
uc.cells <- read.csv("./SCP259/Imm.barcodes2.tsv", header = F)
uc.genes <- read.csv("./SCP259/Imm.genes.tsv", header = F)
dimnames(uc.sparse) <- list(uc.genes$V1,uc.cells$V1)
# tab delimited file requires read.delim
uc.meta <- read.delim("./SCP259/all.meta2.txt", header = TRUE, sep = "\t")
# remove first row of non-cell data
uc.meta <- uc.meta[-1,]
# trim non-myeloid cells
myeloid <- c("CD69+ Mast", "CD69- Mast", "Cycling Monocytes", "DC1", "DC2", "Inflammatory Monocytes", "Macrophages")
uc.meta <- uc.meta[which(uc.meta$Cluster %in% myeloid),]
uc.sparse <- uc.sparse[,which(colnames(uc.sparse) %in% uc.meta$NAME)]
# or trim non-immune cells (only necessary in meta)
immune <- c("Plasma", "CD8+ IELs", "Tregs", "Cycling T", "CD8+ IL17+", "CD8+ LP", "CD4+ Activated Fos-hi", "CD4+ Activated Fos-lo", "Cycling B", "CD4+ PD1+", "Follicular", "CD4+ Memory", "GC", "ILCs", "NKs", "MT-hi", "CD69+ Mast", "CD69- Mast", "Cycling Monocytes", "DC1", "DC2", "Inflammatory Monocytes", "Macrophages")
uc.meta <- uc.meta[which(uc.meta$Cluster %in% immune),]
# generate seurat object
uc <- CreateSeuratObject(counts = uc.sparse, project = "uc", min.cells = 3, min.features = 200)
# add meta data to seurat object
uc$cell_type <- uc.meta$Cluster
uc$disease <- uc.meta$Health
uc$patient <- uc.meta$Subject
# remove extras
rm(uc.data, uc.genes, uc.cells, uc.meta, uc.sparse, myeloid, immune, evens)
# 27310 cells, 18671 genes
dim(uc)
save(list = "uc", file = "./UC Raw Data.RData")

### Human Data: Harmony Batch Correction
# Visualize by patient prior to batch correction
DimPlot(uc, reduction = 'pca', group.by = "patient")
uc <- RunHarmony(uc, "patient", plot_convergence = TRUE)
# Visualize by patient after to batch correction
DimPlot(uc, reduction = 'harmony', group.by="patient")

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

### Human and Mouse Merging
# extract variable features for uc and mimmune
var_genes_human <- VariableFeatures(uc)
var_genes_mouse <- VariableFeatures(mimmune)
var_genes <- unique(c(var_genes_human, var_genes_mouse))
# add species to meta
uc$species <- "human"
mimmune$species <- "mouse"
# remove t, b cells and neutros
mimmune <- subset(mimmune, idents = c("Ly6cHi Monocytes", "Itgax- Macrophages", "Itgax+ Macrophages", "Resident Monocytes", "Alox15+ Macrophages"))
# remove mast cells
uc <- subset(uc, idents = c("Cd14Hi/Cd16Lo Macrophages", "Inflammatory Monocytes", "DC2", "Inflammatory Macrophages", "pre-DC2", "Nonclassical Monocytes", "Inflammatory Monocytes-2", "Cd14Hi/Cd16Lo Macrophages", "DC1", "pre-DC1", "Cycling Monocytes"))
# remove non-variable genes
uc <- uc[var_genes,]
mimmune <- mimmune[var_genes,]
# combine the two human and mouse objects
hsmm <- merge(uc, mimmune)

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
