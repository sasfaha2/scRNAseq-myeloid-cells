library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cluster)
library(reshape2)
library(SoupX)

# SoupX Ambient RNA Removal in Immune Cells (showing control only)
control <- Read10X(data.dir = "../Data/Control Filtered/")
mimmune_cells <- colnames(mimmune)
# clean up cell names
# use substring() to remove "control_", "dss_" and "oxa_"
control <- control[,which(colnames(control) %in% substring(mimmune_cells[1:904],9,26))]
# create soup object
soup.control <- SoupChannel(raw.cont, control)
# perform following preprocessing and clustering
# transfering clustering info
soup.control <- setClusters(soup.control, setNames(control$seurat_clusters, colnames(control)))
# autoestcont method
soup.control <- autoEstCont(soup.control)
