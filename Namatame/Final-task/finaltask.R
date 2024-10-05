# library
# ------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
# ------------------------


# Load the data
# ------------------------
KD1.data <- Read10X(data.dir = "./KD_GSM6042979_RAW/")
KD2.data <- Read10X(data.dir = "./KD_GSM6042980_RAW/")
KD3.data <- Read10X(data.dir = "./KD_GSM6470208_RAW/")
HC1.data <- Read10X(data.dir = "./HC1_GSM4509015_RAW/")
HC2.data <- Read10X(data.dir = "./HC2_GSM4509023_RAW/")
HC4.data <- Read10X(data.dir = "./HC4_GSM4509029_RAW/")
# ------------------------


# Create SeuratObject
# ------------------------
## KD1
# ---------
KD1 <- CreateSeuratObject(counts = KD1.data, project = "KD1", min.cells = 3, min.features = 200)
KD1 [["percent.mt"]] <- PercentageFeatureSet(KD1, pattern = "^MT-")
# VlnPlot(KD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ---------

## KD2
# ---------
KD2 <- CreateSeuratObject(counts = KD2.data, project = "KD2", min.cells = 3, min.features = 200)
KD2 [["percent.mt"]] <- PercentageFeatureSet(KD2, pattern = "^MT-")
# VlnPlot(KD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ---------

## KD3
# ---------
KD3 <- CreateSeuratObject(counts = KD3.data, project = "KD3", min.cells = 3, min.features = 200)
KD3 [["percent.mt"]] <- PercentageFeatureSet(KD3, pattern = "^MT-")
# VlnPlot(KD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ---------

## HC1
# ---------
HC1 <- CreateSeuratObject(counts = HC1.data, project = "HC1", min.cells = 3, min.features = 200)
HC1 [["percent.mt"]] <- PercentageFeatureSet(HC1, pattern = "^MT-")
# VlnPlot(HC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ---------

## HC2
# ---------
HC2 <- CreateSeuratObject(counts = HC2.data, project = "HC2", min.cells = 3, min.features = 200)
HC2 [["percent.mt"]] <- PercentageFeatureSet(HC2, pattern = "^MT-")
# VlnPlot(HC2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ---------

## HC4
# ---------
HC4 <- CreateSeuratObject(counts = HC4.data, project = "HC4", min.cells = 3, min.features = 200)
HC4 [["percent.mt"]] <- PercentageFeatureSet(HC4, pattern = "^MT-")
# VlnPlot(HC4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ---------
# ------------------------


# All Data Merge
# ------------------------
KDALL <- merge(KD1, y = c(KD2, KD3), add.cell.ids = c("KD1", "KD2", "KD3"), project = "KDALL")
HCALL <- merge(HC1, y = c(HC2, HC4), add.cell.ids = c("HC1", "HC2", "HC4"), project = "HCALL")
cell<- c(KDALL,HCALL)
for (i in 1:length(cell)) {
  cell[[i]] <- NormalizeData(cell[[i]], verbose = F)
  cell[[i]] <- FindVariableFeatures(cell[[i]], selection.method = "vst", 
                                    nfeature = 2000, verbose = F)
}

cell.anchors <- FindIntegrationAnchors(object.list = cell, dims = 1:30)
cell.int <- IntegrateData(anchorset = cell.anchors)

# Load reference data from celldex
ref <- celldex::HumanPrimaryCellAtlasData()
# Run SingleR to infer cell types of pbmc dataset using reference data
results <- SingleR(test = as.SingleCellExperiment(cell.int1), ref = ref, labels = ref$label.main)
# Add inferred cell type labels to pbmc object
cell.int1$singlr_labels <- results$labels
cell.int1$KDandHC <- results$labels
# ------------------------


# QC & Normalization
# ------------------------
## QC
# ---------
KD1 <- subset(KD1, subset = percent.mt < 20)
KD2 <- subset(KD2, subset = percent.mt < 20)
KD3 <- subset(KD3, subset = percent.mt < 20)
HC1 <- subset(HC1, subset = percent.mt < 20)
HC2 <- subset(HC2, subset = percent.mt < 20)
HC4 <- subset(HC4, subset = percent.mt < 20)
# ---------

## Normalization
# ---------
KD1 <- NormalizeData(KD1)
KD2 <- NormalizeData(KD2)
KD3 <- NormalizeData(KD3)
HC1 <- NormalizeData(HC1)
HC2 <- NormalizeData(HC2)
HC4 <- NormalizeData(HC4)
# ---------
# ------------------------





# Separate into Subsets
# ------------------------
## KD1




# ------------------------
