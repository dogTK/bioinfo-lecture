install.packages("dplyr")
install.packages("Seurat") # Docker imageを使ってる方はv5がインストールされるためスキップ推奨
install.packages("patchwork")

library(dplyr)
library(Seurat)
library(patchwork)

pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(
  counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200
)
pbmc


# MT-で始まる全ての遺伝子の集合をミトコンドリア遺伝子の集合として使用する
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)

# LogNormalize
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize", scale.factor = 10000
)

# pbmcに再度入れる
pbmc <- NormalizeData(pbmc)


# FindVariableFeatures関数を使ってデータを用意
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# top10の遺伝子名を記述のため、headを取得
top10 <- head(VariableFeatures(pbmc), 10)

# プロット
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

GetAssayData(object = pbmc, slot = "scale.data")


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
##PCAの結果を確認します。
##上での実行結果はpbmc[["pca"]]に格納されています。
##各主成分での特徴量を表示
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

##第1、2主成分での特徴量とその係数を表示
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

##各細胞を第1主成分、第2主成分でプロット
DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 10, balanced = TRUE)


pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

# セルのクラスタリングを行う
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAPの実行
pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap")
