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
options(Seurat.object.assay.version = "v3")
#devtools::install_version("dbplyr", version = "2.3.4")
packageVersion("dbplyr")
# ------------------------


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

#preprocessing
# ------------------------
ALL_KD <- ALL_KD %>%
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10) %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 0.5)
# ------------------------

# Celldex
# ------------------------
ref <- celldex::HumanPrimaryCellAtlasData()
# ------------------------


# Cell labeling
# ------------------------
results <- SingleR(test = as.SingleCellExperiment(ALL_KD), ref = ref, labels = ref$label.main)
ALL_KD$singlr_labels <- results$labels
# ------------------------


# Check labels without duplicates
# ------------------------
unique_labels <- unique(ALL_KD$singlr_labels)
print(unique_labels)
# ------------------------


# Dimplot
# ------------------------
DimPlot(ALL_KD, reduction = "umap", group.by = 'singlr_labels' ,label = TRUE)
# ------------------------


# Get list of genes in ALL_KD
# ------------------------
gene_li <- row.names(ALL_KD)
print(gene_li)
write.table(gene_li, file = "gene_li.txt", quote = FALSE, row.names = FALSE)
# ------------------------


# Feature Plot
# ------------------------
FeaturePlot(ALL_KD, features = c("CX3CR1"))
FeaturePlot(ALL_KD, features = c("UQCRHL"))
FeaturePlot(ALL_KD, features = c("ACTB"))
FeaturePlot(ALL_KD, features = c("GAPDH"))
FeaturePlot(ALL_KD, features = c("IFNG"))
FeaturePlot(ALL_KD, features = c("TNF"))
# ------------------------


# Extract Subset
# ------------------------
## ---------
#  Extract subset by CX3CR1 expression level
ALL_KD_sub <- subset(ALL_KD,subset = CX3CR1>2)
class(ALL_KD_sub)

# Extract data with the label Monocyte
Monocyte <- subset(ALL_KD_sub,subset = singlr_labels == "Monocyte")
Monocyte$singlr_labels

# Dimplot of Monocyte
#DimPlot(Monocyte,reduction = "umap", group.by = 'singlr_labels' ,label = TRUE)
# No particularly useful findings were found.
# ------------------------


# DEG analyze
# ------------------------
deg <- FindAllMarkers(object = Monocyte, only.pos = T)
## Aggregate number of DEGs per cluster
table(deg$cluster)
# ------------------------


# Plot result
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

## Molded for clusterProfiler
for(i in names(gene_list)){
  gene_list[[i]] <- gene_list[[i]]$ENTREZID
}

## Gene Ontology
GOBP <- compareCluster(geneClusters = gene_list, 
                       fun = "enrichGO", 
                       ont="BP",
                       OrgDb = "org.Hs.eg.db")

## dotplot
## Link : 
#dotplot(GOBP)

## barplot
GOBP_bar <- enrichGO(unlist(gene_list), OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
barplot(GOBP_bar)
### Specify showCategory as 20.
par(mar = c(5, 4, 4, 2),cex.axis=0.3, cex.lab=1)
barplot(GOBP_bar,showCategory=10)

## Fixed to original margins
par(mar=c(5, 4, 4, 2), cex.axis=1, cex.lab=1)

## KEGG Pathway Enrichment Analysis
#enrich_result <- compareCluster(geneClusters = gene_list,
#                               fun = "enrichKEGG",
#                                organism = "hsa")
#enrich_result <- enrichKEGG(gene = gene_list, organism = 'hsa', pvalueCutoff = 0.05)
## barplot
#barplot(enrich_result, showCategory = 20)
# ------------------------
