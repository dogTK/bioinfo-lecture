# import library
library(dplyr)
library(Seurat)
library(patchwork)
# データのパス指定
pbmc.data = Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/");
# Seurat-Objectを作成
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# MT-で始まる遺伝子の集合をミトコンドリア由来の遺伝子としてみなす
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
# ヴァイオリンプロット図示
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNAの数-ミトコンドリア遺伝子の数（割合？）の散布図
plot1 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# RNAの数-
plot2 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# 両方のグラフを合わせたものを表示
plot1 + plot2


# 対数変換（pbmcデータ全体で正規化メソッドを作る）
# その他正規化メソッド:Center-log-ratio, Relative-counts
# ↑後でよく調べておく
pbmc = NormalizeData(pbmc, normalization.method = "LogNormalize", scale.facto=10000)
# pbmcに入れる
pbmc = NormalizeData(pbmc)


# 外れ値を見つける
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# top10の遺伝子名を取得
top10_val = head(VariableFeatures(pbmc), 10)

# プロット
val_plot1 = VariableFeaturePlot(pbmc)
val_plot2 = LabelPoints(plot = val_plot1, points = top10_val, repel = TRUE, xnudge = 0, ynudge = 0)
val_plot1 + val_plot2

# 主成分分析準備
all.genes = rownames(pbmc)
pbmc      = ScaleData(pbmc, features = all.genes)
# 主成分分析実行
pbmc      = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# 出力
print( pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# 寄与率図示
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# Dimplot
DimPlot(pbmc, reduction = "pca")
# Dim-Heatmap
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# JackStrawPlot(3m40s程度かかる)
pbmc = JackStraw(pbmc, num.replicate = 100)
pbmc = ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

# ElbowPlot
ElbowPlot(pbmc)

# クラスタリング
pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)
pbmc = RunUMAP(pbmc, dims = 1:10)

# Dim-plot
DimPlot(pbmc, reduction = "umap")