library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cluster)
library(reshape2)

### Quality Control, First Round (showing control only)
control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^mt-")
# estimate 1.5*IQR for number of UMIs expressed; cells 1.5*IQR to be removed
control.iqr <- as.vector(quantile(control$nFeature_RNA)[4] - quantile(control$nFeature_RNA)[2])
control.upper <- as.vector((control.iqr * 1.5) + quantile(control$nFeature_RNA)[4])
control <- subset(control, subset = nFeature_RNA 200 & nFeature_RNA < control.upper)
# remove cells with 10% mito transcripts
control <- subset(control, subset = percent.mt < 10)
# remove essential 0s for genes (genes that are 0 across cells)
control.zero <- which(!as.vector(rowSums(control) == 0))
# merge control, dss and oxa into a single dataset
myeloid <- merge(control, y = c(dss, oxa), add.cell.ids = c("control", "dss", "oxa"), project = "myeloid")

### Preprocessing
myeloid <- NormalizeData(myeloid, normalization.method = "LogNormalize", scale.factor = 1e4)
myeloid <- FindVariableFeatures(myeloid, selection.method = 'vst', nfeatures = 2000)
myeloid <- ScaleData(myeloid, vars.to.regress = c("percent.mt"), features = rownames(myeloid))
myeloid <- RunPCA(myeloid, features = VariableFeatures(object = myeloid))

### Quality Control, Second Round
# Subset epithelial and immune cells; remove low quality, doublet and erythroid clusters (7, 8, 14)
mepi <- subset(myeloid, idents = c("0", "2", "4", "9", "13", "15", '16', "17"))
mimmune <- subset(myeloid, idents = c("1", "3", "5", "6", "10", "11", "12", "18", "19", "20"))
