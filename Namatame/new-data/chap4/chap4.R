install.packages("dplyr")
install.packages("Seurat")
install.packages("patchwork")

library(dplyr)
library(Seurat)
library(patchwork)

## Load the dataset
hc1.data <- read.csv(file="GSM5397337_LCWE_count_matrix.csv.gz")
hc2.data <- read.csv(file="GSM5397338_PBS_count_matrix.csv.gz")

## Initialize the Seurat object with the raw (non-normalized data).
hc1 <- CreateSeuratObject(counts = hc1.data, project = "hc1", min.cells = 3, min.features = 200)
hc2 <- CreateSeuratObject(counts = hc2.data, project = "hc2", min.cells = 3, min.features = 200)


###QC and selecting cells for further analysis
## The [[ operator can add columns to object metadata. This is a great place to stash QC stats
hc1[["percent.mt"]] <- PercentageFeatureSet(hc1, pattern = "^mt-")
hc2[["percent.mt"]] <- PercentageFeatureSet(hc2, pattern = "^mt-")

## Visualize QC metrics as a violin plot
options(repr.plot.width = 8, repr.plot.height = 5) #size change
VlnPlot(hc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)
VlnPlot(hc2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(hc1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


##select cells
hc1 <- subset(hc1, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 10)
hc2 <- subset(hc2, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt <10)


###integration 
###normalization
##Identification of highly variable features (feature selection)
Tcell <- c(hc1,hc2)
for (i in 1:length(Tcell)) {
  Tcell[[i]] <- NormalizeData(Tcell[[i]], verbose = F)
  Tcell[[i]] <- FindVariableFeatures(Tcell[[i]], selection.method = "vst", 
                                     nfeature = 2000, verbose = F)
}

Tcell.anchors <- FindIntegrationAnchors(object.list = Tcell, dims = 1:30)
Tcell.int <- IntegrateData(anchorset = Tcell.anchors, dims = 1:30)

### Scaling and clustering
##Scaling the data
Tcell.int <- ScaleData(Tcell.int, verbose = T) 
#Perform linear dimensional reduction
Tcell.int <- RunPCA(Tcell.int, npcs = 30, verbose = T)

#Elbow Plot
ElbowPlot(Tcell.int)

##Cluster the cells
Tcell.int1<- FindNeighbors(Tcell.int, dims = 1:15)
Tcell.int1 <- FindClusters(Tcell.int1, resolution = 0.5)

#UMAP
Tcell.int1 <- RunUMAP(Tcell.int1, dims = 1:15)
DimPlot(Tcell.int1, reduction = "umap", label = TRUE)

# Visualization
p1 <- DimPlot(Tcell.int1, reduction = "umap", group.by="orig.ident")
p2 <- DimPlot(Tcell.int1, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


DimPlot(Tcell.int1, reduction = "umap", split.by = "orig.ident")


# Check if BiocManager is installed, if not, install it
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required packages using BiocManager
BiocManager::install("SingleR", force = TRUE)
BiocManager::install("celldex", force = TRUE)
BiocManager::install("SingleCellExperiment", force = TRUE)
install.packages("devtools")
devtools::install_version("dbplyr", version = "2.3.4")

# Load required libraries
library(SingleR)
library(celldex)
library(SingleCellExperiment)

# Load reference data from celldex
ref <- celldex::MouseRNAseqData()

# Run SingleR to infer cell types of pbmc dataset using reference data
results <- SingleR(test = as.SingleCellExperiment(Tcell.int1), ref = ref, labels = ref$label.main)

# Add inferred cell type labels to pbmc object
Tcell.int1$singlr_labels <- results$labels

# Visualize cell types in a UMAP plot with labels
DimPlot(Tcell.int1, reduction = 'umap', group.by = 'singlr_labels', label = TRUE)

# Visualize expression of selected genes split by stimulation status
FeaturePlot(Tcell.int1, features = c("Cd3e","Il1r1","Bhlhe40","Tnf"), split.by = "orig.ident", group.by = 'singlr_labels', max.cutoff = 3, cols = c("grey", "red"))


Idents(Tcell.int1) <- factor(Idents(Tcell.int1))
markers.to.plot <- c("Cd3e", "Bhlhe40", "Tnf")
DotPlot(Tcell.int1, features = markers.to.plot, split.by = "orig.ident", group.by = "singlr_labels") +
  RotatedAxis()


# Generate violin plots for specified genes, separated by stimulation status and grouped by cell type
VlnPlot(Tcell.int1, features = c("Cd3e", "Bhlhe40", "Tnf"), pt.size = 0.1, split.by = "orig.ident",group.by = "singlr_labels", combine = FALSE)
