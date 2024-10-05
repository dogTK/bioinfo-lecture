if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("satijalab/seurat-data")

library(Seurat)
library(SeuratData)
library(patchwork)


#installation problem
options(timeout=300)  # タイムアウトを300秒（5分）に設定
download.file("http://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.1.0.tar.gz", destfile = "ifnb.SeuratData_3.1.0.tar.gz")

install.packages("ifnb.SeuratData_3.1.0.tar.gz", repos = NULL, type = "source")

# install dataset
InstallData("ifnb")
# load dataset
ifnb <- LoadData("ifnb")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

head(ifnb@meta.data)


# Find integration ancher
immune.anchors <- FindIntegrationAnchors(
  object.list = ifnb.list, anchor.features = features
)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

head(ifnb@meta.data)
DefaultAssay(immune.combined) <- "integrated"


library(magrittr)  # Load the magrittr package for the %>% operator

# Now you can use the pipe operator
immune.combined <- immune.combined %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(immune.combined, reduction = "umap", split.by = "stim")

immune.combined <- RenameIdents(immune.combined,
                                `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated",
                                `8` = "DC", `9` = "B Activated", `10` = "Mk", `11` = "pDC", `12` = "Eryth",
                                `13` = "Mono/Mk Doublets", `14` = "HSPC"
)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(
  immune.combined,
  reduction = "umap",
  group.by = "seurat_annotations",
  label = TRUE,
  repel = TRUE
)
p1 + p2


Idents(immune.combined) <- factor(
  Idents(immune.combined),
  levels = c(
    "HSPC", "Mono/Mk Doublets", "pDC", "Eryth", "Mk", "DC", "CD14 Mono",
    "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated",
    "CD4 Naive T", "CD4 Memory T"
  )
)
markers.to.plot <- c(
  "CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7",
  "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1",
  "CCL2", "S100A9", "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2", "HBB",
  "TSPAN13", "IL3RA", "IGJ", "PRSS57"
)
DotPlot(
  immune.combined,
  features = markers.to.plot,
  cols = c("blue", "red"),
  dot.scale = 3,
  split.by = "stim"
) + RotatedAxis()

# Load visualization libraries
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Extract and process data for "CD4 Naive T" cells
t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "stim"
avg.t.cells <- as.data.frame(
  log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
)
avg.t.cells$gene <- rownames(avg.t.cells)

# Extract and process data for "CD14 Mono" cells
cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- as.data.frame(
  log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA)
)
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

# Specify genes to be highlighted in plots
genes.to.label <- c(
  "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8"
)

# Create scatter plots and label specified genes
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) +
  geom_point() +
  ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) +
  geom_point() +
  ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)

# Display plots side by side
p1 + p2


# Create a new identifier combining cell type and stimulation status
immune.combined$celltype.stim <- paste(
  Idents(immune.combined),
  immune.combined$stim,
  sep = "_"
)

# Backup current cell type identifiers and update to the new combined identifier
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"

# Identify genes differentially expressed between stimulated and control B cells
b.interferon.response <- FindMarkers(
  immune.combined,
  ident.1 = "B_STIM",
  ident.2 = "B_CTRL",
  verbose = FALSE
)

# Visualize expression of selected genes split by stimulation status
FeaturePlot(
  immune.combined,
  features = c("CD3D", "GNLY", "IFI6"),
  split.by = "stim",
  max.cutoff = 3,
  cols = c("grey", "red")
)

# Generate violin plots for specified genes, separated by stimulation status and grouped by cell type
plots <- VlnPlot(
  immune.combined,
  features = c("LYZ", "ISG15", "CXCL10"),
  split.by = "stim",
  group.by = "celltype",
  pt.size = 0,
  combine = TRUE
)

# Display the generated plots in a single column layout
wrap_plots(plots = plots, ncol = 1)




