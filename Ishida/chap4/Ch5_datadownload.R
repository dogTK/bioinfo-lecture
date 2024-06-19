library(dplyr)
library(Seurat)
library(patchwork)


pbmc_v2.data <- Read10X(data.dir = "./GSE139386_RAW/")

pbmc_v2 <- CreateSeuratObject(
  counts = pbmc_v2.data, 
  project = "pbmc_v23k",
  min.cells = 3, 
  min.features = 200
)

pbmc_v2[["percent.mt"]] <- PercentageFeatureSet(pbmc_v2, pattern = "^MT-")
VlnPlot(
  pbmc_v2,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3
)
