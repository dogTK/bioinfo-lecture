# library
# ------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
options(Seurat.object.assay.version = "v3")
# ------------------------

# Version Management
# ------------------------
#devtools::install_github("satijalab/seurat", ref = "tags/v4.3.0")
#pkgbuild::check_build_tools(debug = TRUE)
#packageVersion("Seurat.object.assay")
# ------------------------

# Load the data
# ------------------------
KD1.data <- Read10X(data.dir = "./KD_GSM6042979_RAW/")
# ------------------------


# Create SeuratObject
# ------------------------
## KD1
# ---------
KD1 <- CreateSeuratObject(counts = KD1.data, project = "KD1", min.cells = 3, min.features = 200)
KD1 [["percent.mt"]] <- PercentageFeatureSet(KD1, pattern = "^MT-")
# VlnPlot(KD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ---------


# Celldex
# ------------------------
ref <- celldex::HumanPrimaryCellAtlasData()
# ------------------------

# Extract  subset
# ------------------------
KD1_sub <- subset(KD1,subset=CX3CR1>2)
KD1_sub <- SingleR(test = as.SingleCellExperiment(KD1_sub), ref = ref, labels = ref$label.main)
KD1_sub$singlr_labels <- KD1_sub$labels
KD1_sub$KDandHC <- KD1_sub$labels
KD1_sub$singlr_labels
# ------------------------
# singlr_labels="Monocyte"


cell.anchors <- FindIntegrationAnchors(object.list = cell, dims = 1:30)
#object <- JoinLayers(object = object, layers = layer)
#cell <- JoinLayers(CreateSeuratObject(cell))
cell.int <- IntegrateData(anchorset = cell.anchors)

# 
# 








# Extract Subset
Monocyte <- subset(cell, idents = "Monocyte")

## Normalization
# ---------
KD1 <- NormalizeData(KD1)
# ---------
# ------------------------
