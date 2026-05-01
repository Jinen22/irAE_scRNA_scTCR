# Figure 2.R
# This script generates the main figures for immune cell analysis (Figure 2).
# Input data should be placed in the 'data/' directory.
# Download required files from NGDC (OMIX:OMIX016674) and put them in 'data/'.

### 1. Set orders and color palettes
# Cell type order (original labels)
cellorder <- c("Naive T", "CD4 CTLA4+", "Treg", "Trm", "CTL", "Memory T", 
               "Proliferating T","NKT","NK", "B cell", 
               "Macrophage", "Monocyte","Neutrophil", "DC")

# Cell type order with numeric prefix for plotting
cellorder <- c("1 Naive T", "2 CD4 CTLA4+", "3 Treg", "4 Trm", "5 CTL", "6 Memory T", 
               "7 Proliferating T", "8 NKT","9 NK", 
               "10 B cell", "11 Macrophage", "12 Monocyte","13 Neutrophil", "14 DC")

# Group order
grouporder <- c("Saline", "PD1")

# Tissue order
tissueorder <- c('Blood','Gut','Heart','Liver','Ln','Lung','Skin','Thyroid','Tumor')

### 2. Define color palettes
# Group colors
mypal_group <- c("#E99D22", "#3A5AA5")
names(mypal_group) <- grouporder

# Tissue colors
library(ggsci)
mypal_tissue <- pal_npg()(length(tissueorder))
names(mypal_tissue) <- tissueorder

# Cell type colors (enough colors for all cell types)
mypal_cell <- c("#3C5488FF", "#F39B7FFF","#377DB8", "#4DAF4A", "#E51A1D", "#A4572B", "#FF7F00", 
                "#984EA3", "#F781BE", "#33A9CF",
                "#F4CA18", "#4C856B", "#ACACCA", "#CBC7BA", "#AAD47A","#B78D91", "#CA449E")[1:length(cellorder)]
names(mypal_cell) <- cellorder



#' Fig2
#' Immune cell DimPlot on UMAP
#-----
# Set working directory to the location of this script (or use project root)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # for RStudio; remove if not using
# Alternatively, just assume working directory is the repository root
# All file paths are relative to the repository root

library(Seurat)
library(ggpubr)  # for theme_classic etc.

# Load the annotated Seurat object from the 'data' folder
seuObject <- readRDS("data/Mouse_anit-PD1_irAEs_scRNA_immune_cell_label.rds")
table(seuObject$label.sum)

### 1. Modify labels for plotting
df.meta <- seuObject@meta.data

df.meta$label <- factor(df.meta$label.sum, levels = cellorder)
df.meta <- df.meta[order(df.meta$label), ]

cellnum <- as.data.frame(table(df.meta$label))
test <- paste(rep(1:14, cellnum$Freq), df.meta$label, sep = " ")
table(test)
df.meta$groupnum1 <- test

df.meta <- df.meta[colnames(seuObject), ]
seuObject$groupnum1 <- df.meta$groupnum1

# Simplify tissue names by removing suffixes
seuObject$tissue <- unlist(lapply(strsplit(seuObject$tissue, split = "_"), function(x) x[1]))

### 2. Redefine order and colors (consistent with above)
cellorder <- c("1 Naive T", "2 CD4 CTLA4+", "3 Treg", "4 Trm", "5 CTL", "6 Memory T", 
               "7 Proliferating T", "8 NKT","9 NK", 
               "10 B cell", "11 Macrophage", "12 Monocyte","13 Neutrophil", "14 DC")
grouporder <- c("Saline", "PD1")
tissueorder <- c('Blood','Gut','Heart','Liver','Ln','Lung','Skin','Thyroid','Tumor')

# Cell type colors
mypal_cell <- c("#3C5488FF", "#F39B7FFF","#377DB8", "#4DAF4A", "#E51A1D", "#A4572B", "#FF7F00", 
                "#984EA3", "#F781BE", "#33A9CF",
                "#F4CA18", "#4C856B", "#ACACCA", "#CBC7BA", "#AAD47A","#B78D91", "#CA449E")[1:length(cellorder)]
names(mypal_cell) <- cellorder

# Group colors
mypal_group <- c("#E99D22", "#3A5AA5")
names(mypal_group) <- grouporder

# Tissue colors
mypal_tissue <- pal_npg()(length(tissueorder))
names(mypal_tissue) <- tissueorder

Idents(seuObject) <- factor(seuObject$groupnum1, levels = cellorder)

## 1.1 UMAP colored by cell type
DimPlot(seuObject, label = T, repel = T, cols = mypal_cell, pt.size = 3, reduction = "umap.rpca",
        label.size = 5, raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16)) 
ggsave(filename = "Fig2_All_cell_dimplot_celltype_figure_legend.pdf", width = 6.4, height = 5)


## 1.2 UMAP colored by tissue
DimPlot(seuObject, cols = mypal_tissue, group.by = "tissue", reduction = "umap.rpca") +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig2_All_cell_dimplot_tissue.pdf", width = 6.2, height = 5)

## 1.3 UMAP colored by group (treatment)
DimPlot(seuObject, cols = mypal_group, group.by = "orig.ident", reduction = "umap.rpca") +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig2_All_cell_dimplot_group.pdf", width = 6.1, height = 5)


### Save the modified Seurat object for downstream analyses
saveRDS(seuObject, file = "data/Mouse_anit-PD1_irAEs_scRNA_immune_cell_label.rds")
#-----





#' Fig2
#' Cell type definition dotplot
#' Dotplot showing marker gene expression per cluster
#-----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # ensure working directory

library(Seurat)

seuObject <- readRDS("data/Mouse_anit-PD1_irAEs_scRNA_immune_cell_label.rds")

## 1.1 Prepare data and gene list
# Cell order (without numeric prefix) for dotplot
cellorder <-c("NK","NKT","CTL", "Trm","Proliferating T", "Memory T", "Naive T", "CD4 CTLA4+", "Treg",
              "B cell","Neutrophil", "Monocyte","Macrophage", "DC")
Idents(seuObject) <- factor(seuObject$label.sum, levels = rev(cellorder))

# Marker genes for each cell type
genes <- list("NK" = c("Klrd1", 'Klrk1', "Klrb1c",'Ncr1'),
              "NKT"= c("Cd3d", "Cd3e", "Cd3g"),
              "CTL" = c("Cd8a", "Cd8b1", "Gzma", "Gzmb", 'Gzmk', 'Ifng', 'Fasl'),
              "Trm" = c("Itgae", "Itga1", "Actn2", "Cd7"),
              "Proliferating T"= c("Mki67",'Rrm2','Ccna2','Kif11'),
              "Memory T" = c("Il7r", 'Il18r1','Cxcr6', 'Rora','Icos'),
              "Naive T" = c("Ccr7", "Sell", "Lef1", "Tcf7"),
              "CD4 CTLA4+"= c("Cd4","Ctla4", 'Cd5','Cd28'),
              "Treg" = c("Foxp3", 'Tnfrsf4','Il2ra','Ikzf2'),
              "B cell" = c("Cd79a", "Cd79b","Cd19", "Ms4a1"),
              "Neutrophil"= c("Csf3r", "S100a8", "S100a9", "Mmp9"),
              "Monocyte" = c('Cd300a', 'Mgst1', 'Ecm1','Fn1'),
              "Macrophage" = c('Cd63','Mrc1','Stab1','Itgb5','C1qb'),
              "DC"= c("Plbd1", "Ifitm3", 'Plac8', 'Cst3'))

DotPlot(seuObject, features = genes, scale = T, cols = c("#4393C3", "#D6604D")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"))
ggsave(filename = "Fig2_cell_type_dotplot.pdf", width = 18, height = 6)
#-----




#' Fig2
#' Cell type proportion barplot
#' Compare proportions between PD1 and Saline groups
#' (1) For all cells combined
#' (2) By tissue, showing log2 fold change of PD1/Saline ratios
#' Note: Skin has very few immune cells, causing extreme ratio changes; 
#'       log2FC values >2 or <-2 are capped at ±2 for visualization.
#-----
library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### 1. Data input
seuObject <- readRDS("data/Mouse_anit-PD1_irAEs_scRNA_immune_cell_label.rds")

### 1.2 Barplot of cell fraction by group and tissue
cellorder <- c("1 Naive T", "2 CD4 CTLA4+", "3 Treg", "4 Trm", "5 CTL", "6 Memory T", 
               "7 Proliferating T", "8 NKT","9 NK", 
               "10 B cell", "11 Macrophage", "12 Monocyte","13 Neutrophil", "14 DC")
mypal_cell <- c("#3C5488FF", "#F39B7FFF","#377DB8", "#4DAF4A", "#E51A1D", "#A4572B", "#FF7F00", 
                "#984EA3", "#F781BE", "#33A9CF",
                "#F4CA18", "#4C856B", "#ACACCA", "#CBC7BA", "#AAD47A","#B78D91", "#CA449E")[1:length(cellorder)]
names(mypal_cell) <- cellorder

tissueorder <- c('Blood','Gut','Heart','Liver','Ln','Lung','Skin','Thyroid','Tumor')
df.all <- data.frame("Group" = seuObject$orig.ident, 
                     "Cluster" = seuObject$groupnum1,
                     "Tissue" = seuObject$tissue)
p.plot <- list()
for (i in tissueorder) {
  
  # Filter for the specific tissue
  df <- df.all[df.all$Tissue == i, 1:2]
  df$Cluster <- factor(df$Cluster, levels = cellorder) 
  
  df.group <- group_by(df, Group, Cluster)
  ratio <- summarise(df.group, num = n())
  # Calculate proportion within each group
  ratio <- group_by(ratio, Group) %>% mutate(ratio.group = num / sum(num))
  cell.Saline <- ratio[ratio$Group == "Saline", c("Cluster","ratio.group")]
  cell.pd1 <- ratio[ratio$Group == "PD1", c("Cluster","ratio.group")]
  # Exclude cell types not present in both groups (cannot calculate fold change)
  cell.use <- intersect(cell.Saline$Cluster, cell.pd1$Cluster)
  rownames(cell.Saline) <- as.character(cell.Saline$Cluster)
  rownames(cell.pd1) <- as.character(cell.pd1$Cluster)
  cell.Saline <- cell.Saline[as.character(cell.use), "ratio.group"]
  cell.pd1 <- cell.pd1[as.character(cell.use), "ratio.group"]
  
  df.plot <- data.frame(FC = (cell.pd1+0.01) / (cell.Saline+0.01),
                        celltype = cell.use)
  df.plot$ratio.group <- log2(df.plot$ratio.group)
  # Cap extreme fold changes
  df.plot$ratio.group[df.plot$ratio.group > 2] <- 2
  df.plot$ratio.group[df.plot$ratio.group < -2] <- -2
  colnames(df.plot) <- c("FC", "celltype")
  df.plot <- df.plot[order(df.plot$FC, decreasing = T), ]
  df.plot$celltype <- factor(df.plot$celltype, levels = df.plot$celltype)
  mypal_group <- mypal_cell[unique(ratio$Cluster)]
  
  p1 <- ggplot(df.plot, aes(x=celltype, y= FC, fill=celltype)) +
    geom_bar(stat = "identity") +
    # Add labels: right-aligned for positive, left-aligned for negative
    geom_text(data = subset(df.plot, FC < 0),
              aes(x=celltype, y= 0.3, label= paste0(" ", celltype)),
              size = 3, 
              hjust = "inward" ) +  
    geom_text(data = subset(df.plot, FC > 0),
              aes(x=celltype, y= -0.1, label=celltype),
              size = 3, hjust = "outward") +
    
    theme_bw() + 
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    theme(panel.border = element_rect(linewidth = 0.6)) + 
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    coord_flip() +
    scale_fill_manual(values = mypal_group) +
    ylab("Log2 (Fold change)") +
    xlab("Cell types") +
    guides(fill="none") +
    ggtitle(i) +
    ylim(-2, 2) +
    geom_hline(yintercept = 0,  linetype="longdash", lwd = 0.5, color="grey40")
  p.plot[[i]] <- p1
}

p.plot[[1]] + p.plot[[2]] + p.plot[[3]] + p.plot[[4]] + p.plot[[5]] + p.plot[[6]] + p.plot[[7]] + p.plot[[8]] + p.plot[[9]] +
  plot_layout(ncol = 5)
ggsave("All_immune_cell_1th_2th_tissue_log2FC_barplot.pdf", width = 15, height = 6)



### 3. Combined tissue (overall fold change barplot)
cellorder <- c("1 Naive T", "2 CD4 CTLA4+", "3 Treg", "4 Trm", "5 CTL", "6 Memory T", 
               "7 Proliferating T", "8 NKT","9 NK", 
               "10 B cell", "11 Macrophage", "12 Monocyte","13 Neutrophil", "14 DC")
mypal_cell <- c("#3C5488FF", "#F39B7FFF","#377DB8", "#4DAF4A", "#E51A1D", "#A4572B", "#FF7F00", 
                "#984EA3", "#F781BE", "#33A9CF",
                "#F4CA18", "#4C856B", "#ACACCA", "#CBC7BA", "#AAD47A","#B78D91", "#CA449E")[1:length(cellorder)]
names(mypal_cell) <- cellorder

df <- data.frame("Group" = seuObject$orig.ident, 
                 "Cluster" = seuObject$groupnum1)
df$Cluster <- factor(df$Cluster, levels = cellorder)

df.group <- group_by(df, Group, Cluster)
ratio <- summarise(df.group, num = n())
ratio <- group_by(ratio, Group) %>% mutate(ratio.group = num / sum(num))
cell.Saline <- ratio[ratio$Group == "Saline", c("Cluster","ratio.group")]
cell.pd1 <- ratio[ratio$Group == "PD1", c("Cluster","ratio.group")]
cell.use <- intersect(cell.Saline$Cluster, cell.pd1$Cluster)
rownames(cell.Saline) <- cell.Saline$Cluster
rownames(cell.pd1) <- cell.pd1$Cluster
cell.Saline <- cell.Saline[as.character(cell.use), "ratio.group"]
cell.pd1 <- cell.pd1[as.character(cell.use), "ratio.group"]

df.plot <- data.frame(FC = (cell.pd1+0.01) / (cell.Saline+0.01),
                      celltype = cell.use)
colnames(df.plot) <- c("FC", "celltype")
df.plot$FC <- log2(df.plot$FC)

df.plot <- df.plot[order(df.plot$FC, decreasing = T), ]
df.plot$celltype <- factor(df.plot$celltype, levels = df.plot$celltype)
mypal_group <- mypal_cell[unique(ratio$Cluster)]

ggplot(df.plot, aes(x=celltype, y= FC, fill=celltype)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5)) +
  coord_flip() +
  scale_fill_manual(values = mypal_group) +
  ylab("Log2 (Fold change)") +
  xlab("Cell types") +
  geom_hline(yintercept = 0,  linetype="longdash", lwd = 0.5, color="grey40")
ggsave("All_immune_cell_All_tissue_log2FC_barplot.pdf", width = 5, height = 3)
#-----