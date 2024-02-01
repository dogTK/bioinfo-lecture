# library
# ------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(clusterProfiler)
library(celldex)
library(ggplot2)
library(org.Hs.eg.db)
#library(igraph)
#options(Seurat.object.assay.version = "v3")
#devtools::install_version("dbplyr", version = "2.3.4")
#packageVersion("dbplyr")
# ------------------------

#install.packages("igraph")
#BiocManager::install("clusterProfiler")
#install.packages("ggplot2")
#install.packages("org.Hs.eg.db")
#install.packages("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")


# Load data
# ------------------------
KD1.data <- Read10X(data.dir = "./KD_GSM6042979_RAW/")
KD2.data <- Read10X(data.dir = "./KD_GSM6042980_RAW/")
KD3.data <- Read10X(data.dir = "./KD_GSM6470208_RAW/")
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
# ------------------------

# Normalization
# ------------------------
KD1 <- NormalizeData(KD1, normalization.method = "LogNormalize", scale.factor = 10000)
KD2 <- NormalizeData(KD2, normalization.method = "LogNormalize", scale.factor = 10000)
KD3 <- NormalizeData(KD3, normalization.method = "LogNormalize", scale.factor = 10000)
# ------------------------


# Merge
# ------------------------
ALL_KD <- merge(KD1, y = c(KD2, KD3), add.cell.ids = c("KD1", "KD2", "KD3"), project = "ALL_KD",merge.data = TRUE)
table(ALL_KD$orig.ident)
# ------------------------


# Celldex
# ------------------------
ref <- celldex::HumanPrimaryCellAtlasData()
# ------------------------


# Extract  subset
# ------------------------
ALL_KD_sub <- subset(ALL_KD,subset = CX3CR1>2)
results <- SingleR(test = as.SingleCellExperiment(ALL_KD_sub), ref = ref, labels = ref$label.main)
ALL_KD_sub$singlr_labels <- results$labels
#ALL_KD_sub$KDandHC <- results$labels
#ALL_KD_sub$singlr_labels <- results$labels
# ------------------------
class(ALL_KD_sub)

# Extract Subset
# ------------------------
## ---------
# Extract data with the label Monocyte
Monocyte <- subset(ALL_KD_sub,subset = singlr_labels == "Monocyte")
Monocyte$singlr_labels
# ------------------------


#preprocessing
# ------------------------
Monocyte <- Monocyte %>%
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10) %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 0.5)
# ------------------------


# DEG analyze
# ------------------------
deg <- FindAllMarkers(object = Monocyte, only.pos = T)
## Aggregate number of DEGs per cluster
table(deg$cluster)
# ------------------------


# dotplot
# ------------------------
## Vectorize DEGs by clusters and compile them into a list
gene_list <- split(deg$gene, deg$cluster)

## Conversion from gene name to gene ID
for(i in names(gene_list)){
  gene_list[[i]] <- bitr(gene_list[[i]], 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = "org.Hs.eg.db")
}

## Number of genes remaining
length(unlist(gene_list))/2

for(i in names(gene_list)){
  gene_list[[i]] <- gene_list[[i]]$ENTREZID
}

## Gene Ontology
GOBP <- compareCluster(geneClusters = gene_list, 
                       fun = "enrichGO", 
                       ont="BP",
                       OrgDb = "org.Hs.eg.db")

dotplot(GOBP)
# ------------------------
