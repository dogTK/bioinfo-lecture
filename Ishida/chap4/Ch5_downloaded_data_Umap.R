library(dplyr)
library(Seurat)
library(patchwork)

pbmc.data <- Read10X(data.dir= "./GSE149689/")

pbmc <- CreateSeuratObject(
  counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200
)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- subset(
  pbmc, nFeature_RNA >= 200 & nFeature_RNA <= 5000 & percent.mt < 15
)

pbmc <- pbmc %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = VariableFeatures(object = pbmc)) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5)

DimPlot(pbmc, reduction = "umap", label = T)




if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("SingleCellExperiment")


library(SingleR)
library(celldex)
library(SingleCellExperiment)


ref <- celldex::HumanPrimaryCellAtlasData()


results <- SingleR(
  test = as.SingleCellExperiment(pbmc), ref = ref, labels = ref$label.main
)


pbmc$singlr_labels <- results$labels


DimPlot(pbmc, reduction = "umap", group.by = "singlr_labels", label = TRUE)

FeaturePlot(pbmc, features = c("CD14", "CD3E"))
