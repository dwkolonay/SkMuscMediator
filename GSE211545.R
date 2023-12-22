# This uses standard (i.e. NOT SCTransform method)

# install packages
## set up library using dir already created
#lib <- c("/igm/home/cnh008/r_libs/multiome")
## install.packages defaults to first element of .libPaths
#.libPaths(c(lib, .libPaths()))
##install packages
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.17",
#                     lib = lib)
#BiocManager::install(c("EnsDb.Mmusculus.v79",
#                       "BSgenome.Mmusculus.UCSC.mm10",
#                       "JASPAR2020",
#                       "TFBSTools",
#                       "GenomicRanges",
#                       "motifmatchr",
#                       'AnnotationFilter',
#                       'BiocGenerics',
#                       'GenomicFeatures',
#                       'IRanges',
#                       'Rsamtools',
#                       'S4Vectors',
#                       'ggbio',
#                       'AnnotationDbi'),
#                     lib = lib,
#                     lib.loc = lib,
#                     force = TRUE)
#install.packages(c("SoupX",
#                   "Seurat",
#                   "cowplot",
#                   "Signac",
#                   "ggplot2",
#                   "patchwork",
#                   "dplyr",
#                   "tidyverse",
#                   "BiocManager",
#                   "patchwork",
#                   "harmony",
#                   "future",
#                   "sctransform"))
#setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
#install.packages("Signac")

# set up library using dir already created
lib <- c("/igm/home/cnh008/r_libs/multiome")
# install.packages defaults to first element of .libPaths
.libPaths(c(lib, .libPaths()))
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(pheatmap)

setwd("/igm/home/cnh008/projects/Baskin_lab/REVISIONS/")

# read in data downloaded from GEO (https://urldefense.com/v3/__https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211543__;!!AU3bcTlGKuA!CH8sbrEHBI-WtPYCom1QLSp2g4_iJ2W4KbjNQaMaMa1sJBJMD5gYbQJJiujVOkV8fGKjeQU2IZaq2bbn3jB_ZFsz3ZYrXHOk2lAS1ok4ABA$ ):
# read in E14 data
counts <- Read10X("../GSE211545/raw_GEO/GSM6475260_E14_rna")
E14 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "E14"
)

# read in E18 data
counts <- Read10X("../GSE211545/raw_GEO/GSM6475261_E18_rna")
E18 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "E18"
)

# read in P5 data
counts <- Read10X("../GSE211545/raw_GEO/GSM6475262_P5_rna")
P5 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "P5"
)

# read in Adult data
counts <- Read10X("../GSE211545/raw_GEO/GSM6475263_Adult_rna")
Adult <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "Adult"
)

DefaultAssay(E14) <- 'RNA'
DefaultAssay(E18) <- 'RNA'
DefaultAssay(P5) <- 'RNA'
DefaultAssay(Adult) <- 'RNA'

E14$orig.ident <- 'E14'
E18$orig.ident <- 'E18'
P5$orig.ident <- 'P5'
Adult$orig.ident <- 'Adult'

mergeData <- merge(E14, y=c(E18, P5, Adult), add.cell.ids = c("E14", "E18", "P5", "Adult"), project = "Full Merge Project")

mergeData[["percent.mt"]] <- PercentageFeatureSet(mergeData, pattern = "^mt-")

mergeData <- subset(
  x = mergeData,
  subset = nCount_RNA < 25000 &
    nCount_RNA > 400 &
    percent.mt < 20
)

DefaultAssay(mergeData) <- 'RNA'

saveRDS(mergeData, "./GSE211545/merged.rds")

so_list_normalized<- SplitObject(object = mergeData, split.by = "orig.ident")

for (i in 1:length(x = so_list_normalized)) {
  so_list_normalized[[i]] <- NormalizeData(object = so_list_normalized[[i]], 
                                           verbose = FALSE)
  so_list_normalized[[i]] <- FindVariableFeatures(object = so_list_normalized[[i]], 
                                                  selection.method = "vst", 
                                                  nfeatures = 2000, 
                                                  verbose = FALSE)
  so_list_normalized[[i]] <- PercentageFeatureSet(object = so_list_normalized[[i]], 
                                                  pattern = "^mt-", 
                                                  col.name = "percent.mt")
  so_list_normalized[[i]]<- ScaleData(object = so_list_normalized[[i]],
                                      features = rownames(so_list_normalized[[i]]))
  so_list_normalized[[i]]<- RunPCA(object = so_list_normalized[[i]], npcs =30)
  so_list_normalized[[i]]<- RunUMAP(so_list_normalized[[i]], dims = 1:30)
}

#Set default assay to RNA for all objects in list 
for (i in seq_len(length(so_list_normalized))) {
  DefaultAssay(so_list_normalized[[i]]) <- "RNA"
}
#integrate RNA using rpca
features <- SelectIntegrationFeatures(so_list_normalized,
                                      nfeatures = 3000)
integration_anchors <- FindIntegrationAnchors(object.list = so_list_normalized,
                                              normalization.method = "LogNormalize",
                                              anchor.features = features,
                                              dims = 1:30,
                                              reduction = "rpca",
                                              k.anchor = 20)
so <- IntegrateData(anchorset = integration_anchors,
                    normalization.method = "LogNormalize",
                    features.to.integrate = rownames(so_list_normalized[[1]]@assays$RNA@counts),
                    new.assay.name = "integratedRNA",
                    dims = 1:30)
DefaultAssay(so) <- "integratedRNA"
# Run the standard workflow for visualization and clustering
so <- ScaleData(so, features = rownames(so@assays$integratedRNA@data))
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30)
so <- RunUMAP(so, reduction = "pca", reduction.name = "umap.RNA",dims = 1:30)
so <- FindClusters(so, algorithm = 3, verbose = FALSE)
DimPlot(so, reduction = "umap.RNA")
saveRDS(so, "./GSE211545/integrated.rds")

# use previously annotated iteration to annotate data:
# note: this data was annotated using marker genes from manuscript described below
ref <- readRDS("../GSE211545/rna/outs/objects/integrated.rds")
labels <- matrix(nrow = nrow(ref@meta.data),
                 ncol = 2) %>%
  as.data.frame()
colnames(labels) <- c("barcode", "celltype")
labels$barcode <- rownames(ref@meta.data)
labels$celltype <- ref$celltypev2
so <- AddMetaData(so, metadata = labels$celltype, col.name = "celltype")
Idents(so) <- rownames(so@meta.data)
celltype <- as.character(labels$celltype)
names(celltype) <- levels(so)
so <- RenameIdents(so, celltype)
so$celltype <- Idents(so)
# make some marker plots
DefaultAssay(so) <- "RNA"
plot <- VlnPlot(so,
                features = c("Ttn", "Pdgfra", "Pax7", "Mymk", "Col1a1", "Flt1",
                             "Mrc1", "Myh11", "Pparg"),
                pt.size = 0,
                group.by = "seurat_clusters",
                stack = TRUE,
                fill.by = "ident")
ggsave(plot,
       filename = "./GSE211545/celltype_markers_vlnplot_by_clusters.jpeg",
       height = 10,
       width = 10)
# make heatmap with markers from manuscript for each cell type
genes <- c("Ttn", "Trdn", "Mylk4", "Neb", "Dmd",
           "Tmeff2", "Ebf2", "Abca8a", "Col3a1", "Fap",
           "Pax7", "Pde1c", "Map2", "Fgf12", "Hs6st2",
           "Myo16", "Megf10", "Tmem178b", "Slc24a3", "Dpp6",
           "Col11a1", "Col1a1", "Col1a2", "Itgbl1", "Mkx",
           "Flt1", "Ptprb", "Prex2", "Pecam1", "Mecom",
           "F13a1", "Mctp1", "Mrc1", "Slc9a9", "Rbpj",
           "Myh11", "Itga1", "Cacnb2", "Dgkb", "Kcnab1",
           "Pparg", "Cidec", "Prkar2b", "Pde3b", "Acsl1")
DefaultAssay(so) <- "integratedRNA"
# make heatmap with just clusters 0:10 so you can see them
cells <- rownames(so@meta.data[so@meta.data$seurat_clusters %in% c(0:10),])
heatmap <- DoHeatmap(so,
                     features = genes,
                     cells = cells,
                     group.by = "seurat_clusters",
                     size = 4)
ggsave(heatmap,
       filename = "./GSE211545/celltype_markers_heatmap_by_clusters_0_to_10.jpeg",
       height = 8,
       width = 10)
# make heatmap with just clusters 11:27 so you can see them
cells <- rownames(so@meta.data[so@meta.data$seurat_clusters %in% c(11:27),])
heatmap <- DoHeatmap(so,
                     features = genes,
                     cells = cells,
                     group.by = "seurat_clusters",
                     size = 4)
ggsave(heatmap,
       filename = "./GSE211545/celltype_markers_heatmap_by_clusters_11_to_27.jpeg",
       height = 8,
       width = 10)
Idents(so) <- "seurat_clusters"
so.small <- subset(so, downsample = 300)
heatmap <- DoHeatmap(so.small,
                     features = genes,
                     group.by = "seurat_clusters",
                     size = 4)
ggsave(heatmap,
       filename = "./GSE211545/celltype_markers_heatmap_by_clusters.jpeg",
       height = 8,
       width = 10)

# make UMAPs with clusters and celltypes
Idents(so) <- "seurat_clusters"
p1 <- DimPlot(so) &
  ggtitle("Seurat Clusters")
Idents(so) <- "celltype"
p2 <- DimPlot(so) &
  ggtitle("Annotated w/ ViolinPlot Markers")
plot <- p1 + p2
ggsave(plot,
       filename = "./GSE211545/GSE211545_umaps.jpeg",
       width = 21,
       height = 5)

# make heatmap with celltypev2
Idents(so) <- "celltype"
so.small <- subset(so, downsample = 300)
so.small$celltypev2 <- factor(so.small$celltype,
                              levels = c("Myonuclei",
                                         "Mesenchymal",
                                         "MuSCs",
                                         "Myoblasts",
                                         "Tenocytes",
                                         "Endothelial",
                                         "Macrophages",
                                         "Smooth muscle",
                                         "Adipocytes"))
heatmap <- DoHeatmap(so.small,
                     features = genes,
                     group.by = "celltype",
                     size = 4)
ggsave(heatmap,
       filename = "./GSE211545/celltype_markers_heatmap_by_celltype.jpeg",
       height = 8,
       width = 10)

# make umap labeled by celltype split by sample
Idents(so) <- "celltype"
so$orig.ident <- factor(so$orig.ident,
                        levels = c("E14", "E18", "P5", "Adult"))
plot <- DimPlot(so,
                split.by = "orig.ident",
                pt.size = 0.7)
ggsave(plot,
       filename = "./GSE211545/umap_split_by_sample.jpeg",
       width = 28,
       height = 5)

# make umap labeled by sample
Idents(so) <- "orig.ident"
plot <- DimPlot(so)
ggsave(plot,
       filename = "./GSE211545/umap_samples.jpeg",
       width = 7,
       height = 5)

# save integrated data again
saveRDS(so, "./GSE211545/integrated.rds")

# condense data into average z-score (so@assays$integratedRNA@scale.data)
Idents(so) <- "celltype"
for (cell in unique(so$celltype)){
  sub <- subset(so, idents = cell)
  zscores <- AverageExpression(sub,
                               assays = "integratedRNA",
                               slot = "scale.data",
                               group.by = "orig.ident") %>%
    as.data.frame()
  write.csv(zscores,
            file = paste0("./GSE211545/", cell, "_average_zscores.csv"))
}
files <- list.files("./GSE211545/",
                    pattern = "_zscores.csv",
                    full.names = TRUE)
names(files) <- basename(files) %>%
  gsub("_average_zscores.csv", "", .)
zscores <- lapply(files, read.csv)
for (cell in names(zscores)){
  df <- zscores[[cell]]
  colnames(df) <- gsub("integratedRNA.",
                       paste0(cell, "_"),
                       colnames(df))
  zscores[[cell]] <- df
}
test <- Reduce(merge, zscores)
write.csv(test, "./GSE211545/all_average_zscores.csv")

# make heatmap from average zscores for each cell type
genes_list <- read.csv("../gene_list.csv")
mediator <- str_to_title(genes_list$mediator)
control <- genes_list$control[1:12]
for (cell in names(zscores)){
  df <- zscores[[cell]]
  colnames(df) <- gsub(paste0(cell, "_"), "", colnames(df))
  rownames(df) <- df$X
  df <- df[,-c(1)]
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
                      cluster_rows = FALSE,
                      scale = "row")
  ggsave(heatmap,
         filename = paste0("./GSE211545/average_heatmaps/", cell, "_mediator_average_zscores_heatmap.jpeg"), 
         height = 10, 
         width = 8)
  write.csv(sub,
            paste0("./GSE211545/average_heatmaps/", cell, "_mediator_average_zscores.csv"))
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
         filename = paste0("./GSE211545/average_heatmaps/", cell, "_control_average_zscores_heatmap.jpeg"), 
         height = 10, 
         width = 8)
  write.csv(sub,
            paste0("./GSE211545/average_heatmaps/", cell, "_control_average_zscores.csv"))
}

# make single cell heatmaps with gene lists
genes_list <- read.csv("../gene_list.csv")
mediator <- str_to_title(genes_list$mediator)
control <- genes_list$control[1:12]
for (cell in unique(so$celltype)){
  Idents(so) <- "celltype"
  sub <-subset(so, idents = cell)
  Idents(sub) <- "orig.ident"
  sub_small <- subset(sub, downsample = 300)
  heatmap <- DoHeatmap(sub_small,
                       features = mediator,
                       group.by = "orig.ident",
                       size = 4,
                       assay = "integratedRNA",
                       disp.min = -1,
                       disp.max = 1) &
    scale_fill_gradient2(high = "red",
                          mid = "white",
                          low = "blue",
                          midpoint = 0) &
    ggtitle("Mediator Genes")
  ggsave(heatmap,
         filename = paste0("./GSE211545/sc_heatmaps/", cell, "_mediator_sc_heatmap.jpeg"), 
         height = 10, 
         width = 8)
  heatmap <- DoHeatmap(sub_small,
                       features = control,
                       group.by = "orig.ident",
                       size = 4,
                       assay = "integratedRNA",
                       disp.min = -1,
                       disp.max = 1) &
    scale_fill_gradient2(high = "red",
                         mid = "white",
                         low = "blue",
                         midpoint = 0) &
    ggtitle("Control Genes")
  ggsave(heatmap,
         filename = paste0("./GSE211545/sc_heatmaps/", cell, "_control_sc_heatmap.jpeg"), 
         height = 10, 
         width = 8)
}

# make single cell heatmaps (NOT DOWNSAMPLED)
for (cell in unique(so$celltype)){
  Idents(so) <- "celltype"
  sub <-subset(so, idents = cell)
  Idents(sub) <- "orig.ident"
  heatmap <- DoHeatmap(sub,
                       features = mediator,
                       group.by = "orig.ident",
                       size = 4,
                       assay = "integratedRNA",
                       disp.min = -1,
                       disp.max = 1) &
    scale_fill_gradient2(high = "red",
                         mid = "white",
                         low = "blue",
                         midpoint = 0) &
    ggtitle("Mediator Genes")
  ggsave(heatmap,
         filename = paste0("./GSE211545/sc_heatmaps/", cell, "_mediator_sc_heatmap_all_data.jpeg"), 
         height = 10, 
         width = 8)
  heatmap <- DoHeatmap(sub,
                       features = control,
                       group.by = "orig.ident",
                       size = 4,
                       assay = "integratedRNA",
                       disp.min = -1,
                       disp.max = 1) &
    scale_fill_gradient2(high = "red",
                         mid = "white",
                         low = "blue",
                         midpoint = 0) &
    ggtitle("Control Genes")
  ggsave(heatmap,
         filename = paste0("./GSE211545/sc_heatmaps/", cell, "_control_sc_heatmap_all_data.jpeg"), 
         height = 10, 
         width = 8)
}

# make dotplots to replace heatmap
setwd("/igm/home/cnh008/projects/Baskin_lab/REVISIONS")
so <- readRDS("./GSE211545/integrated.rds")

# Figure 3A
p <- DotPlot(so,
             features = rev(mediator),
             assay = "RNA",
             group.by = "orig.ident",
             scale = TRUE,
             cols = "RdBu",
             col.min = -1,
             col.max = 1)  &
  scale_y_discrete(position = "right") &
  coord_flip() &
  theme(legend.position="bottom",
        legend.box = "vertical",
        legend.justification = "left",
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 11)) &
  guides(color = guide_colorbar(title = 'Scaled Expression'))
ggsave(p,
       filename = "./GSE211545/Figure3A.jpeg",
       width = 3.9,
       height = 7)

# Supplemental Figure 3A
p <- DotPlot(so,
             features = rev(control),
             assay = "RNA",
             group.by = "orig.ident",
             scale = TRUE,
             cols = "RdBu",
             col.min = -1,
             col.max = 1)  &
  scale_y_discrete(position = "right") &
  coord_flip() &
  theme(legend.position="bottom",
        legend.box = "vertical",
        legend.justification = "left",
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 11)) &
  guides(color = guide_colorbar(title = 'Scaled Expression'))
ggsave(p,
       filename = "./GSE211545/SuppFigure3A.jpeg",
       width = 3.75,
       height = 5)
  

# Supplemental Figure 3B
so$sample <- paste(so$orig.ident)
so$sample <- gsub("Adult", "P40", so$sample)
so$time_cell <- paste(so$celltype, "\n", so$sample)
so$time_cell <- factor(so$time_cell,
                       levels = c("Myonuclei \n E14",
                                  "Myonuclei \n E18",
                                  "Myonuclei \n P5",
                                  "Myonuclei \n P40",
                                  "Myoblasts \n E14",
                                  "Myoblasts \n E18",
                                  "Myoblasts \n P5",
                                  "MuSCs \n E14",
                                  "MuSCs \n E18",
                                  "MuSCs \n P5",
                                  "MuSCs \n P40",
                                  "Smooth muscle \n E14",
                                  "Smooth muscle \n E18",
                                  "Smooth muscle \n P5",
                                  "Smooth muscle \n P40",
                                  "Adipocytes \n E18",
                                  "Adipocytes \n P5",
                                  "Adipocytes \n P40",
                                  "Tenocytes \n E14",
                                  "Tenocytes \n E18",
                                  "Tenocytes \n P5",
                                  "Tenocytes \n P40",
                                  "Mesenchymal \n E14",
                                  "Mesenchymal \n E18",
                                  "Mesenchymal \n P5",
                                  "Mesenchymal \n P40",
                                  "Endothelial \n E14",
                                  "Endothelial \n E18",
                                  "Endothelial \n P5",
                                  "Endothelial \n P40",
                                  "Macrophages \n E14",
                                  "Macrophages \n E18",
                                  "Macrophages \n P5",
                                  "Macrophages \n P40"))
p <- DotPlot(so,
             features = rev(mediator),
             assay = "RNA",
             group.by = "time_cell",
             scale = TRUE,
             cols = "RdBu",
             col.min = -1,
             col.max = 1)  &
  scale_y_discrete(position = "right") &
  coord_flip() &
  theme(legend.position="bottom",
        legend.box = "horizontal",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 10)) &
  guides(color = guide_colorbar(title = 'Scaled Expression'))
ggsave(p,
       filename = "./GSE211545/SuppFigure3B.jpeg",
       width = 11,
       height = 7)
