library(Seurat)
library(SeuratData)
library(patchwork)

# install dataset
InstallData("ifnb")

# load dataset
LoadData("ifnb")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

# Find integration ancher
immune.anchors <- FindIntegrationAnchors(
  object.list = ifnb.list, anchor.features = features
)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

library(magrittr)  # Load the magrittr package for the %>% operator

# Now you can use the pipe operator
immune.combined <- immune.combined %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(
  immune.combined,
  ident.1 = 5,
  ident.2 = c(0, 3),
  min.pct = 0.25
)
head(cluster5.markers, n = 5)

# トップの遺伝子を選択
top_genes <- head(row.names(cluster5.markers))

# バイオリンプロットを作成
VlnPlot(immune.combined, features = top_genes, group.by = "seurat_clusters")


cluster0.markers <- FindMarkers(
  immune.combined,
  ident.1 = 0,
  logfc.threshold = 0.25,
  test.use = "roc",
  only.pos = TRUE
)
head(cluster0.markers, n = 5)


VlnPlot(
  immune.combined,
  features = c("FTL", "CD63"),
  pt.size = 0.1,
  group.by = "seurat_clusters"
)


# Identify markers for every cluster compared to all remaining cells, report only the positive ones
cluster.markers <- FindAllMarkers(
  immune.combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
library(dplyr)
cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>%
  print(n = 28)


plots <- VlnPlot(
  immune.combined,
  features = c("CCL8"),
  group.by = "seurat_clusters",
  pt.size = 0.1,
  combine = FALSE
)



marker_genes <- cluster.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 2) %>%
  pull(gene)

unique_marker_genes <- unique(marker_genes)

DotPlot(
  immune.combined,
  features = unique_marker_genes,
  group.by = "seurat_clusters"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(immune.combined) <- "RNA"


# multtestパッケージをインストール
BiocManager::install('multtest')

# metapパッケージをインストール
install.packages('metap')


cluster6.markers <- FindConservedMarkers(
  immune.combined,
  ident.1 = 6,
  grouping.var = "stim"
)
head(cluster6.markers)

FeaturePlot(
  immune.combined,
  features = c(
    "CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"
  ),
  min.cutoff = "q9"
)


immune.combined <- RenameIdents(
  immune.combined,
  `0` = "CD14 Mono",
  `1` = "CD4 Naive T",
  `2` = "CD4 Memory T",
  `3` = "CD16 Mono",
  `4` = "B",
  `5` = "CD8 T",
  `6` = "NK",
  `7` = "T activated",
  `8` = "DC",
  `9` = "B Activated",
  `10` = "Mk",
  `11` = "pDC",
  `12` = "Eryth",
  `13` = "Mono/Mk Doublets"
)
DimPlot(immune.combined, label = TRUE)
