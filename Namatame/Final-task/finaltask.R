# library
# ------------------------
library(dplyr)
library(Seurat)
library(patchwork)
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

## KD1
# ------------------------
KD1 <- CreateSeuratObject(counts = KD1.data, project = "KD1", min.cells = 3, min.features = 200)
KD1 [["percent.mt"]] <- PercentageFeatureSet(KD1, pattern = "^MT-")
# VlnPlot(KD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ------------------------

## KD2
# ------------------------
KD2 <- CreateSeuratObject(counts = KD2.data, project = "KD2", min.cells = 3, min.features = 200)
KD2 [["percent.mt"]] <- PercentageFeatureSet(KD2, pattern = "^MT-")
# VlnPlot(KD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
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
# ------------------------


# QC & Normalization

## KD1
# ------------------------
KD1 <- subset(KD1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 300 & nCount_RNA < 7500)
#VlnPlot(KD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ------------------------





# Separate into Subsets

## KD1
