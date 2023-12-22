# set up library using dir already created
lib <- c("/igm/home/cnh008/r_libs/Seurat_v3.1.5")
# install.packages defaults to first element of .libPaths
.libPaths(c(lib, .libPaths()))
# install packages
# seurat <- c("https://urldefense.com/v3/__https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.1.5.tar.gz__;!!AU3bcTlGKuA!CH8sbrEHBI-WtPYCom1QLSp2g4_iJ2W4KbjNQaMaMa1sJBJMD5gYbQJJiujVOkV8fGKjeQU2IZaq2bbn3jB_ZFsz3ZYrXHOk2lASnb2dpxY$ ")
# install.packages(seurat, repos=NULL, type="source")
# install.packages(c("ggplot2",
#                    "patchwork",
#                    "tidyverse",
#                    "pheatmap"))
# load packages
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(pheatmap)

# set working dir
setwd("/igm/home/cnh008/projects/Baskin_lab/REVISIONS/GSE156498/")

# read in data downloaded from GEO (https://urldefense.com/v3/__https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156498__;!!AU3bcTlGKuA!CH8sbrEHBI-WtPYCom1QLSp2g4_iJ2W4KbjNQaMaMa1sJBJMD5gYbQJJiujVOkV8fGKjeQU2IZaq2bbn3jB_ZFsz3ZYrXHOk2lASCpFFraw$ )
# WT
counts <- Read10X("./raw/GSM4732631_WT_strained/")
WT <- CreateSeuratObject(counts = counts,
                         assay = "RNA", project = "WT")

# D51
counts <- Read10X("./raw/GSM4732632_D51_strained/")
D51 <- CreateSeuratObject(counts = counts,
                          assay = "RNA", project = "D51")

#  For each sample, single-nucleus transcriptomes with fewer than 200 or more than 4,000 genes or with more than 15,000 UMIs were further filtered out from the analysis. 
WT <- subset(WT, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 15000)
D51 <- subset(D51, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 15000)

# Subsequently, data were log normalized and scaled, and principle component analysis was performed using the top 2,000 genes that showed highly variable expression in the integrated dataset. Cell clusters were called using the first 15 principle components under a clustering resolution of 0.6. Dimensional reduction was performed by UMAP using the first 15 principle components. 
data.list <- list("WT" = WT,
                  "D51" = D51)
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, 
                            selection.method = "vst", 
                            nfeatures = 2000)
  x <- ScaleData(x,
                 features = rownames(x))
})
features <- SelectIntegrationFeatures(object.list = data.list)
anchors <- FindIntegrationAnchors(object.list = data.list,
                                  anchor.features = features)
data <- IntegrateData(anchorset = anchors,
                      features.to.integrate = rownames(data@assays$RNA@counts))
DefaultAssay(data) <- "integrated"
data <- ScaleData(data, 
                  verbose = FALSE)
data <- RunPCA(data, 
               npcs = 30, 
               verbose = FALSE)
data <- RunUMAP(data, 
                reduction = "pca", 
                dims = 1:15)
data <- FindNeighbors(data, 
                      reduction = "pca", 
                      dims = 1:15)
data <- FindClusters(data, 
                     resolution = 0.6)
saveRDS(data, "./integrated.rds")

# Marker genes were unbiasedly analyzed using Seurat function FindAllMarkers. 
Idents(data) <- "seurat_clusters"
markers <- FindAllMarkers(data, 
                          assay = "RNA",
                          only.pos = TRUE)
top50 <- markers %>%
  group_by(cluster) %>% 
  slice_max(avg_logFC, n = 50) %>%
  as.data.frame()
write.csv(top50, "./cluster_markers_top50.csv")
# make gene heatmap
top5_genes <- top50 %>%
  group_by(cluster) %>% 
  slice_max(avg_logFC, n = 5) %>%
  as.data.frame()
genes <- unique(top5_genes$gene)
DefaultAssay(data) <- "RNA"
heatmap <- DoHeatmap(data,
                     features = genes,
                     group.by = "seurat_clusters",
                     size = 4)
ggsave(heatmap,
       filename = "./markers_top5_heatmap.jpeg",
       width = 14,
       height = 14,
       units = "in",
       device = "jpeg")




labels <- matrix(nrow = nrow(test@meta.data),
                 ncol = 2) %>%
  as.data.frame()
colnames(labels) <- c("barcode", "celltype")
labels$barcode <- rownames(test@meta.data)
labels$celltype <- test$celltype
data <- AddMetaData(data, metadata = labels$celltype, col.name = "celltype")
Idents(data) <- rownames(data@meta.data)
celltype <- as.character(labels$celltype)
names(celltype) <- levels(data)
data <- RenameIdents(data, celltype)
data$celltype <- Idents(data)







# Next, cluster identities were annotated based on differentially expressed genes, as well as expression of known marker genes for different cell types as described in the manuscript.
genes <- c("Ckm", "Myh2", "Myh1", "Myh4", "Chrne", "Col22a1", "Pax7",
           "Megf10", "Myh3", "Myh11", "Pecam1", "Pdgfra", "Mkx", "Adgre1")
plot <- VlnPlot(data,
                features = genes,
                assay = "RNA",
                group.by = "seurat_clusters",
                pt.size = 0)
ggsave(plot,
       filename = "./markers_violinplot.jpeg",
       width = 14,
       height = 14,
       units = "in",
       device = "jpeg")
genes2 <- c("Ckm", "Myh2", "Myh1", "Myh4", "Chrne")
plot2 <- VlnPlot(data,
                 features = genes2,
                 group.by = "seurat_clusters",
                 assay = "RNA",
                 pt.size = 0,
                 ncol = 1)
ggsave(plot2,
       filename = "./markers_violinplot2.jpeg",
       width = 5,
       height = 14,
       units = "in",
       device = "jpeg")
map <- DoHeatmap(data,
                 features = genes)

# rename idents and make UMAPs with clusters and celltypes
Idents(data) <- "seurat_clusters"
p1 <- DimPlot(data) &
  ggtitle("Seurat Clusters")
cellids <- list("0" = "IIx",
                "1" = "IIb",
                "2" = "IIa",
                "3" = "IIb",
                "4" = "FAP",
                "5" = "EC",
                "6" = "IIx_b",
                "7" = "IIb",
                "8" = "MPH",
                "9" = "MuSC",
                "10" = "SMC",
                "11" = "Myob_RegMyon",
                "12" = "IIx_b",
                "13" = "TC",
                "14" = "MuSC",
                "15" = "FAP",
                "16" = "NMJ")
sum(names(cellids) == levels(data))
data <- RenameIdents(data, cellids)
data$celltype <- Idents(data)
Idents(data) <- "celltype"
p2 <- DimPlot(data) &
  ggtitle("Annotated w/ ViolinPlot Markers")
plot <- p1 + p2
ggsave(plot,
       filename = "./umaps.jpeg",
       width = 14,
       height = 5)
data$orig.ident <- factor(data$orig.ident,
                          levels = c("WT", "D51"))
plot <- DimPlot(data,
                split.by = "orig.ident")
ggsave(plot,
       filename = "./celltype_umap_split_by_sample.jpeg",
       width = 14,
       height = 5)
saveRDS(data, "./integrated.rds")

# scale all genes
DefaultAssay(data) <- "RNA"
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)

# condense data into average z-score (data@assays$integratedRNA@scale.data)
Idents(data) <- "celltype"
for (cell in unique(data$celltype)){
  sub <- subset(data, idents = cell)
  Idents(sub) <- "orig.ident"
  zscores <- AverageExpression(sub,
                               assays = "RNA",
                               slot = "scale.data") %>%
    as.data.frame()
  colnames(zscores) <- gsub("RNA",
                            "",
                            colnames(zscores))
  write.csv(zscores,
            file = paste0("./average_heatmaps/", cell, "_average_zscores.csv"))
}
files <- list.files("./average_heatmaps/",
                    pattern = "_zscores.csv",
                    full.names = TRUE)
names(files) <- basename(files) %>%
  gsub("_average_zscores.csv", "", .)
zscores <- lapply(files, read.csv)
for (cell in names(zscores)){
  df <- zscores[[cell]]
  rownames(df) <- df$X
  df <- df[,-c(1)]
  colnames(df) <- paste0(cell, "_",
                       colnames(df))
  zscores[[cell]] <- df
}

# make heatmap from average zscores for each cell type
genes_list <- read.csv("../../gene_list.csv")
mediator <- str_to_title(genes_list$mediator)
control <- genes_list$control[1:12]
for (cell in names(zscores)){
  df <- zscores[[cell]]
  colnames(df) <- gsub(paste0(cell, "_"), "", colnames(df))
  # make mediator heatmap
  sub <- df[mediator,]
  # set up colors
  paletteLength <- 50
  myColor <- colorRampPalette(c("red", "white", "blue"))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(sub), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(sub)/paletteLength, max(sub), length.out=floor(paletteLength/2)))
  heatmap <- pheatmap(sub,
                      color = myColor,
                      breaks = myBreaks,
                      cluster_cols = FALSE,
                      cluster_rows = FALSE)
  ggsave(heatmap,
         filename = paste0(".//average_heatmaps/", cell, "_mediator_average_zscores_heatmap.jpeg"), 
         height = 10, 
         width = 8)
  write.csv(sub,
            paste0("./average_heatmaps/", cell, "_mediator_average_zscores.csv"))
  # make control heatmap
  sub <- df[control,]
  sub <- sub[-c(grep("NA", rownames(sub))),]
  # set up colors
  paletteLength <- 50
  myColor <- colorRampPalette(c("red", "white", "blue"))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(sub), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(sub)/paletteLength, max(sub), length.out=floor(paletteLength/2)))
  heatmap <- pheatmap(sub,
                      color = myColor,
                      breaks = myBreaks,
                      cluster_cols = FALSE,
                      cluster_rows = FALSE)
  ggsave(heatmap,
         filename = paste0("./average_heatmaps/", cell, "_control_average_zscores_heatmap.jpeg"), 
         height = 10, 
         width = 8)
  write.csv(sub,
            paste0("./average_heatmaps/", cell, "_control_average_zscores.csv"))
}

genes_list <- read.csv("../../gene_list.csv")
mediator <- str_to_title(genes_list$mediator)
control <- genes_list$control[1:12]
for (cell in unique(data$celltype)){
  Idents(data) <- "celltype"
  sub <-subset(data, idents = cell)
  Idents(sub) <- "orig.ident"
  sub_small <- subset(sub, downsample = 300)
  heatmap <- DoHeatmap(sub_small,
                       features = mediator,
                       group.by = "orig.ident",
                       size = 4,
                       assay = "RNA",
                       disp.min = -1,
                       disp.max = 1) &
    scale_fill_gradient2(high = "red",
                         mid = "white",
                         low = "blue",
                         midpoint = 0) &
    ggtitle("Mediator Genes")
  ggsave(heatmap,
         filename = paste0("./sc_heatmaps/", cell, "_mediator_sc_heatmap.jpeg"), 
         height = 10, 
         width = 8)
  heatmap <- DoHeatmap(sub_small,
                       features = control,
                       group.by = "orig.ident",
                       size = 4,
                       assay = "RNA",
                       disp.min = -1,
                       disp.max = 1) &
    scale_fill_gradient2(high = "red",
                         mid = "white",
                         low = "blue",
                         midpoint = 0) &
    ggtitle("Control Genes")
  ggsave(heatmap,
         filename = paste0("./sc_heatmaps/", cell, "_control_sc_heatmap.jpeg"), 
         height = 10, 
         width = 8)
}

# full data
for (cell in unique(data$celltype)){
  Idents(data) <- "celltype"
  sub <-subset(data, idents = cell)
  Idents(sub) <- "orig.ident"
  heatmap <- DoHeatmap(sub,
                       features = mediator,
                       group.by = "orig.ident",
                       size = 4,
                       assay = "RNA",
                       disp.min = -1,
                       disp.max = 1) &
    scale_fill_gradient2(high = "red",
                         mid = "white",
                         low = "blue",
                         midpoint = 0) &
    ggtitle("Mediator Genes")
  ggsave(heatmap,
         filename = paste0("./sc_heatmaps/", cell, "_mediator_sc_heatmap_full_data.jpeg"), 
         height = 10, 
         width = 8)
  heatmap <- DoHeatmap(sub,
                       features = control,
                       group.by = "orig.ident",
                       size = 4,
                       assay = "RNA",
                       disp.min = -1,
                       disp.max = 1) &
    scale_fill_gradient2(high = "red",
                         mid = "white",
                         low = "blue",
                         midpoint = 0) &
    ggtitle("Control Genes")
  ggsave(heatmap,
         filename = paste0("./sc_heatmaps/", cell, "_control_sc_heatmap_full_data.jpeg"), 
         height = 10, 
         width = 8)
}


# making new plots for REVISIONS
data <- readRDS("./integrated.rds")
genes_list <- read.csv("../../gene_list.csv")
mediator <- str_to_title(genes_list$mediator)
control <- genes_list$control[1:12]
control <- control[control %in% rownames(data@assays$RNA@scale.data)]

# make dotplots to replace Figure 7E
for (cell in unique(data$celltype)){
  Idents(data) <- "celltype"
  sub <- subset(data, idents = cell)
  min <- -4
  max <- 4
  p <- DotPlot(sub,
               features = mediator,
               assay = "RNA",
               group.by = "orig.ident",
               scale = FALSE)  &
    scale_color_gradient2(high = "darkred",
                          mid = "white",
                          low = "darkblue",
                          midpoint = 0,
                          limits = c(min, max),
                          breaks = seq(min, max, (max-min)/4)) &
    scale_y_discrete(position = "right") &
    xlab("") &
    coord_flip() &
    ggtitle(paste(cell)) &
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 8),
          plot.margin=unit(c(0.1,0.1,3.5,0.1),"cm"),
          legend.position = c(-0.55, -0.18),
          legend.box = "vertical",
          legend.justification = "left") &
    guides(color = guide_colorbar(title = 'Average Expression',
                                  direction = "horizontal",
                                  title.position = "top"),
           size = guide_legend(nrow = 1,
                               title = 'Percent Expressed'))
  ggsave(p,
         filename = paste0("./dotplots/", cell,"_Figure7E_mediator.jpeg"),
         width = 1.75,
         height = 6)
}

# Control dotplots
for (cell in unique(data$celltype)){
  Idents(data) <- "celltype"
  sub <- subset(data, idents = cell)
  min <- -4
  max <- 4
  p <- DotPlot(sub,
               features = control,
               assay = "RNA",
               group.by = "orig.ident",
               scale = FALSE)  &
    scale_color_gradient2(high = "darkred",
                          mid = "white",
                          low = "darkblue",
                          midpoint = 0,
                          limits = c(min, max),
                          breaks = seq(min, max, (max-min)/4)) &
    scale_y_discrete(position = "right") &
    xlab("") &
    coord_flip() &
    ggtitle(paste(cell)) &
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 8),
          plot.margin=unit(c(0.1,0.1,3.5,0.1),"cm"),
          legend.position = c(-0.55, -0.18),
          legend.box = "vertical",
          legend.justification = "left") &
    guides(color = guide_colorbar(title = 'Average Expression',
                                  direction = "horizontal",
                                  title.position = "top"),
           size = guide_legend(nrow = 1,
                               title = 'Percent Expressed'))
  ggsave(p,
         filename = paste0("./dotplots/", cell,"_Figure7E_control.jpeg"),
         width = 1.75,
         height = 6)
}

# all subset of celltypes combined mediator genes
data$celltype <- gsub("Myob_RegMyon", "Myoblasts\nReg.Myonuclei", data$celltype)
data$orig_cell <- paste0(data$celltype, "\n", data$orig.ident)
Idents(data) <- "orig_cell"
sub <- subset(data,
              idents = c("MuSC\nWT", "MuSC\nD51",
                         "Myoblasts\nReg.Myonuclei\nWT", "Myoblasts\nReg.Myonuclei\nD51",
                         "IIa\nWT", "IIa\nD51",
                         "IIx\nWT", "IIx\nD51",
                         "IIb\nWT", "IIb\nD51",
                         "IIx_b\nWT", "IIx_b\nD51"))
sub$orig_cell <- factor(sub$orig_cell,
                        levels = c("MuSC\nWT", "MuSC\nD51",
                                   "Myoblasts\nReg.Myonuclei\nWT", "Myoblasts\nReg.Myonuclei\nD51",
                                   "IIa\nWT", "IIa\nD51",
                                   "IIx\nWT", "IIx\nD51",
                                   "IIb\nWT", "IIb\nD51",
                                   "IIx_b\nWT", "IIx_b\nD51"))
min <- -2
max <- 3
p <- DotPlot(sub,
             features = mediator,
             assay = "RNA",
             group.by = "orig_cell",
             scale = TRUE)  &
  scale_color_gradient2(high = "darkred",
                        mid = "white",
                        low = "darkblue",
                        midpoint = 0,
                        limits = c(min, max),
                        breaks = seq(min, max, (max-min)/4)) &
  scale_y_discrete(position = "right") &
  xlab("") &
  coord_flip() &
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 8),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.box = "horizontal",
        legend.position = "bottom") &
  guides(color = guide_colorbar(title = 'Scaled Expression',
                                direction = "horizontal"),
         size = guide_legend(nrow = 1,
                             title = 'Percent Expressed'))
ggsave(p,
       filename = paste0("./dotplots/Figure7E_mediator.jpeg"),
       width = 7,
       height = 6)

# all subset of celltypes combined control genes
min <- -2
max <- 3
p <- DotPlot(sub,
             features = control,
             assay = "RNA",
             group.by = "orig_cell",
             scale = TRUE)  &
  scale_color_gradient2(high = "darkred",
                        mid = "white",
                        low = "darkblue",
                        midpoint = 0,
                        limits = c(min, max),
                        breaks = seq(min, max, (max-min)/4)) &
  scale_y_discrete(position = "right") &
  xlab("") &
  coord_flip() &
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 8),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.box = "horizontal",
        legend.position = "bottom") &
  guides(color = guide_colorbar(title = 'Scaled Expression',
                                direction = "horizontal"),
         size = guide_legend(nrow = 1,
                             title = 'Percent Expressed'))
ggsave(p,
       filename = paste0("./dotplots/Figure7E_control.jpeg"),
       width = 7,
       height = 6)


# all celltypes combined mediator genes
Idents(data) <- "orig_cell"
data$orig_cell <- factor(data$orig_cell,
                         levels = c("MuSC\nWT", "MuSC\nD51",
                                    "Myoblasts\nReg.Myonuclei\nWT", "Myoblasts\nReg.Myonuclei\nD51",
                                    "IIa\nWT", "IIa\nD51",
                                    "IIx\nWT", "IIx\nD51",
                                    "IIb\nWT", "IIb\nD51",
                                    "IIx_b\nWT", "IIx_b\nD51",
                                    "SMC\nWT", "SMC\nD51",
                                    "NMJ\nWT", "NMJ\nD51",
                                    "Tenocytes\nWT", "Tenocytes\nD51",
                                    "FAP\nWT", "FAP\nD51",
                                    "EC\nWT", "EC\nD51",
                                    "Macrophages\nWT", "Macrophages\nD51"))
min <- -3
max <- 3
p <- DotPlot(data,
             features = mediator,
             assay = "RNA",
             group.by = "orig_cell",
             scale = TRUE)  &
  scale_color_gradient2(high = "darkred",
                        mid = "white",
                        low = "darkblue",
                        midpoint = 0,
                        limits = c(min, max),
                        breaks = seq(min, max, (max-min)/4)) &
  scale_y_discrete(position = "right") &
  xlab("") &
  coord_flip() &
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 8),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.box = "horizontal",
        legend.position = "bottom") &
  guides(color = guide_colorbar(title = 'Scaled Expression',
                                direction = "horizontal"),
         size = guide_legend(nrow = 1,
                             title = 'Percent Expressed'))
ggsave(p,
       filename = paste0("./dotplots/Figure7E_mediator_all_cells.jpeg"),
       width = 9,
       height = 6)

# all subset of celltypes combined control genes
min <- -1
max <- 3
p <- DotPlot(data,
             features = control,
             assay = "RNA",
             group.by = "orig_cell",
             scale = TRUE)  &
  scale_color_gradient2(high = "darkred",
                        mid = "white",
                        low = "darkblue",
                        midpoint = 0,
                        limits = c(min, max),
                        breaks = seq(min, max, (max-min)/4)) &
  scale_y_discrete(position = "right") &
  xlab("") &
  coord_flip() &
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 8),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.box = "horizontal",
        legend.position = "bottom") &
  guides(color = guide_colorbar(title = 'Scaled Expression',
                                direction = "horizontal"),
         size = guide_legend(nrow = 1,
                             title = 'Percent Expressed'))
ggsave(p,
       filename = paste0("./dotplots/Figure7E_control_all_cells.jpeg"),
       width = 7,
       height = 6)
