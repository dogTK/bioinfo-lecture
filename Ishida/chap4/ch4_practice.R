install.packages("dplyr")
install.packages("Seurat")
install.packages("patchwork")
#install.packages("Matrix")

library(dplyr)
library(Seurat)
library(patchwork)
sessionInfo()
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(
  counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200
)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#バイオリンプロット：分布を視覚的に見やすくする手法。外れ値が見やすくなるため有用
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#代入しているので注意
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)


pbmc <- NormalizeData(
  pbmc, 
  normalization.method = "LogNormalize", scale.factor = 10000
)

pbmc <- NormalizeData(pbmc)

#vstは変動係数の分散版(V/μ)を行った後に対数変換(log(x+1))することで行われる
#vst(Variance Stabilizing Transformation)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#rownames関数では名前そのままに行名を返す
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]])

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap")
