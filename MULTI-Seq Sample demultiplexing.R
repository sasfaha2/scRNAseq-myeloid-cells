library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cluster)
library(reshape2)
library(deMULTIplex)

# Demultiplexing
bar_control <- as.data.frame(read.csv(file="../Control_barTable.csv"))
bar_dss <- as.data.frame(read.csv(file="../DSS_barTable.csv"))
bar_oxa <- as.data.frame(read.csv(file="../Oxa_barTable.csv"))

# making compatible with MULTI-seq format
rownames(bar_control) <- bar_control$X
rownames(bar_dss) <- bar_dss$X
rownames(bar_oxa) <- bar_oxa$X
bar_control <- bar_control[,2:4]
bar_dss <- bar_dss[,2:4]
bar_oxa <- bar_oxa[,2:4]

# https://satijalab.org/seurat/archive/v3.1/hashing_vignette.html
# Select cell barcodes detected by both RNA and MULTI-seq
rownames(bar_control) <- paste("control_", rownames(bar_control), "-1", sep = "")
rownames(bar_dss) <- paste("dss_", rownames(bar_dss), "-1", sep = "")
rownames(bar_oxa) <- paste("oxa_", rownames(bar_oxa), "-1", sep = "")
joint.bcs <- intersect(colnames(myeloid), c(rownames(bar_control), rownames(bar_dss), rownames(bar_oxa)))
ncol(myeloid) == length(joint.bcs)
nrow(rbind(bar_control, bar_dss, bar_oxa)) == length(joint.bcs)
bar <- rbind(bar_control, bar_dss, bar_oxa)[joint.bcs,]

# Add multi-seq data as a new assay independent from RNA
myeloid[["HTO"]] <- CreateAssayObject(counts = t(bar))

# Normalize multi-seq data, here we use centered log-ratio (CLR) transformation
# split by condition
myeloid <- SplitObject(myeloid, split.by = "orig.ident")
myeloid$control <- NormalizeData(myeloid$control, assay = "HTO", normalization.method = "CLR")
myeloid$dss <- NormalizeData(myeloid$dss, assay = "HTO", normalization.method = "CLR")
myeloid$oxa <- NormalizeData(myeloid$oxa, assay = "HTO", normalization.method = "CLR")
myeloid$control <- HTODemux(myeloid$control, assay = "HTO", positive.quantile = 0.99, kfunc = "clara")
myeloid$dss <- HTODemux(myeloid$dss, assay = "HTO", positive.quantile = 0.99, kfunc = "clara")
myeloid$oxa <- HTODemux(myeloid$oxa, assay = "HTO", positive.quantile = 0.99, kfunc = "clara")
myeloid <- merge(myeloid$control, y = c(myeloid$dss, myeloid$oxa), project = "myeloid")

# remove neg cells
myeloid <- SetIdent(myeloid, value = "HTO_classification.global") %>% subset(., idents = "Negative", invert = TRUE)

# adjust hash (sacrifices doublet sensitivity)
myeloid$hash.ID[which(myeloid$HTO_classification.global == "Doublet")] <- myeloid$HTO_maxID[which(myeloid$HTO_classification.global == "Doublet")]
