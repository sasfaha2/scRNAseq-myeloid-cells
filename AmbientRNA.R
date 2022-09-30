library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cluster)
library(reshape2)
library(SoupX)

### Assessing Ambient RNA Contamination (showing DSS vs Oxa only)
load("../epithelial_degbycondition.RData")
load("../immune_degbycondition.RData")
# assign gene name column
rownames(DSSOxa.markers_imm) -> DSSOxa.markers_imm$Name
rownames(DSSOxa.markers_epi) -> DSSOxa.markers_epi$Name
# Which DEG shared by immune and goblet cells (log2FC > 1)
immepi_dssoxa <- DSSOxa.markers_imm$Name[which(DSSOxa.markers_imm$avg_log2FC > 1 | DSSOxa.markers_imm$avg_log2FC < -1)]
immepi_dssoxa <- immepi_dssoxa[which(immepi_dssoxa %in% DSSOxa.markers_epi$Name[which(DSSOxa.markers_epi$avg_log2FC > 1 | DSSOxa.markers_epi$avg_log2FC < -1)])]
# create dataframe of log2FCs for immepi_dssoxa genes
FCimmepi_dssoxa <- data.frame(matrix(nrow = length(immepi_dssoxa), ncol = 2))
rownames(FCimmepi_dssoxa) <- gobimmune
colnames(FCimmepi_dssoxa) <- c("Immune_log2FC", "Gob_log2FC")
FCimmepi_dssoxa$Immune_log2FC <- DSSOxa.markers_imm[immepi_dssoxa,"avg_log2FC"]
FCimmepi_dssoxa$Gob_log2FC <- DSSOxa.markers_epi[immepi_dssoxa,"avg_log2FC"]
# line of best fit calc
coeff <- lm(FCimmepi_dssoxa$Gob_log2FC ~ FCimmepi_dssoxa$Immune_log2FC)
coeff <- coefficients(coeff)
# Visualize relationship
ggplot(FCimmepi_dssoxa, aes(x = Immune_log2FC, y = Gob_log2FC)) +
  geom_point(color = "red", size = 5, alpha = 0.8) +
  ggtitle("Oxa vs. DSS Immune and Goblet DEGs") + 
  theme(plot.title = element_text(size=45, hjust = 0.5), 
                                  text = element_text(size = 30),
                                  axis.text = element_text(size = 20)) +
  geom_text_repel(data = FCimmepi_dssoxa, aes(Immune_log2FC, Gob_log2FC, label = rownames(FCimmepi_dssoxa)), max.overlaps = 50, size = 9) +
  geom_abline(intercept = round(coeff[1],1), slope = round(coeff[2],1), lty = 2, size = 2) +
  geom_hline(yintercept = 0, lty = 1, size = 1) +
  geom_vline(xintercept = 0, lty = 1, size = 1) +
  annotate("text", x = 1.5, y= -0.75, col = "red", label = paste("r^2:", round(cor(FCimmepi_dssoxa$Immune_log2FC,FCimmepi_dssoxa$Gob_log2FC)^2,2)), size = 12)
# pearson's r
cor.test(FCimmepi_dssoxa$Immune_log2FC,FCimmepi_dssoxa$Gob_log2FC)

### SoupX Ambient RNA Removal in Immune Cells (showing control only)
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
