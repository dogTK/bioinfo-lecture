#Chpater30のおさらい
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)

BiocManager::install("multtest")
install.packages("metap")

LoadData("ifnb")
ifnb <- UpdateSeuratObject(ifnb)

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

# Seuratオブジェクト内のクラスター情報を確認
unique(Idents(immune.combined))

#Chapter31ここから
# find all markers distinguishing cluster 5 from clusters 0 and 3
# 例: IMMUNE_CTRL と IMMUNE_STIM の間でマーカーを比較
cluster5_markers <- FindMarkers(
  immune.combined,
  ident.1 = "IMMUNE_CTRL",
  ident.2 = "IMMUNE_STIM",
  min.pct = 0.25
)
head(cluster5_markers, n = 5)

# トップの遺伝子を選択
#トップ遺伝子がscRNA-seq本とちょっと違った。
top_genes <- head(row.names(cluster5_markers))

#トップ遺伝子が異なっていたため、ちゃんと同じようなバイオリンプロットが描けるか、遺伝子名を指定して行った。
# バイオリンプロットを作成
VlnPlot(immune.combined, features = top_genes, group.by = "RNA")
# 複数の遺伝子の発現をバイオリンプロットで表示
VlnPlot(object = immune.combined, features = c("CD2", "CD8A", "CD8B", "HLA-DRA", "TMSB4X", "CD3D"), group.by = "seurat_annotations")
# RNAアッセイを使用してバイオリンプロットを描画
VlnPlot(object = immune.combined, features = c("CD2", "CD8A", "CD8B", "HLA-DRA", "TMSB4X", "CD3D"), group.by = "seurat_annotations", assay = "RNA")

# トップの遺伝子を選択
top_genes <- head(row.names(cluster5_markers), n=5)

# ドットプロットを作成
DotPlot(immune.combined, features = top_genes, group.by = "seurat_clusters")
# ドットプロットを作成(遺伝子指定)
DotPlot(immune.combined, features = c("CD2", "CD8A", "CD8B", "HLA-DRA", "TMSB4X", "CD3D"), group.by = "seurat_annotations", assay = "RNA")

# クラスター情報を seurat_annotations に基づいて変更
Idents(immune.combined) <- "seurat_annotations"

# その後に "CD14 Mono" と他のクラスターの比較を行う
#細胞の名前がすでにアノテーション済みなので、ここでは具体的なクラスター名を指定
cluster_CD14Mono_markers <- FindMarkers(
  immune.combined,
  ident.1 = "CD14 Mono",     # CD14 Mono クラスターを指定
  ident.2 = NULL,            # 他のクラスターを比較対象
  logfc.threshold = 0.25,
  test.use = "roc",
  only.pos = TRUE
)

# 結果を確認
head(cluster_CD14Mono_markers, n = 5)

#assayを指定しておかないと別のデータが参照される。
VlnPlot(
  immune.combined,
  features = c("FTL", "CD63"),
  pt.size = 0.1,
  group.by = "seurat_annotations",
  assay = "RNA"
)

# Identify markers for every cluster compared to all remaining cells, report only the positive ones
cluster.markers <- FindAllMarkers(
  immune.combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>%
  print(n = 28)

VlnPlot(
  immune.combined,
  features = c("CCL8"),
  group.by = "seurat_annotations",
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
  group.by = "seurat_annotations"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# 結果を確認
head(conserved_markers)

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(immune.combined) <- "RNA"
cluster6.markers <- FindConservedMarkers(
  immune.combined,
  ident.1 = "NK",
  grouping.var = "stim"
)
head(cluster6.markers)

DefaultAssay(immune.combined) <- "integrated"

library(magrittr)  # Load the magrittr package for the %>% operator

# Now you can use the pipe operator
immune.combined <- immune.combined %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

FeaturePlot(
  immune.combined,
  features = c(
    "CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"
  ),
  assay = "RNA",
  min.cutoff = "q9"
)

