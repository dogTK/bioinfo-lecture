#start from the last line!!!!!!!!!!!
#start from the last line!!!!!!!!!!!
#start from the last line!!!!!!!!!!!
#start from the last line!!!!!!!!!!!
#start from the last line!!!!!!!!!!!
#start from the last line!!!!!!!!!!!
#start from the last line!!!!!!!!!!!
#start from the last line!!!!!!!!!!!
install.packages("dplyr")
install.packages("devtools")
devtools::install_github("satijalab/seurat", ref = "tags/v4.3.0")
install.packages("patchwork")


library(dplyr)
library(Seurat)
library(patchwork)

# Load the dataset
pbmc.data <- Read10X(data.dir = "./GSE149689/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(
  counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200
)

# Calculate percent mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Quality control
pbmc <- subset(
  pbmc, nFeature_RNA >= 200 & nFeature_RNA <= 5000 & percent.mt < 15
)

"""
pbmc <- pbmc %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = VariableFeatures(object = pbmc)) %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:20) %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 0.5)
"""

pbmc <-  NormalizeData(pbmc) 
pbmc <-FindVariableFeatures(pbmc,selection.method = "vst", nfeatures = 2000) 
pbmc <-ScaleData(pbmc,features = VariableFeatures(object = pbmc))
pbmc <-RunPCA(pbmc) 
pbmc <-RunUMAP(pbmc,dims = 1:20) 
pbmc <-FindNeighbors(pbmc,dims = 1:10)
pbmc <-FindClusters(pbmc,resolution = 0.5)

# Visualization
DimPlot(pbmc, reduction = "umap", label = T)

devtools::install_version("dbplyr", version = "2.3.4")

# Check if BiocManager is installed, if not, install it
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required packages using BiocManager
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("SingleCellExperiment")

# Load required libraries
library(SingleR)
library(celldex)
library(SingleCellExperiment)

# Load reference data from celldex
ref <- celldex::HumanPrimaryCellAtlasData()

# Run SingleR to infer cell types of pbmc dataset using reference data
results <- SingleR(
  test = as.SingleCellExperiment(pbmc), ref = ref, labels = ref$label.main
)
print(results)
# Add inferred cell type labels to pbmc object
pbmc$singlr_labels <- results$labels

# Visualize cell types in a UMAP plot with labels
DimPlot(pbmc, reduction = "umap", group.by = "singlr_labels", label = TRUE)

FeaturePlot(pbmc, features = c("CD14", "CD3E"))
saveRDS(pbmc, file = "pbmc_seurat_object.rds")


#Start FROM HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pbmc <- readRDS("pbmc_seurat_object.rds")
