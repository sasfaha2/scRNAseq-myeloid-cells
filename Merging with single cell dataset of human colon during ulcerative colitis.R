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
