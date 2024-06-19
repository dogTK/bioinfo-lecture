#install.packages("SeuratData")
#remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
#sink("console_output_0527.txt")
#print("Hello, world!")
#summary(mtcars)

#remove.packages("SeuratObject")
remove.packages("Seurat")
#pkgbuild::check_build_tools(debug = TRUE)
devtools::install_github("satijalab/seurat", ref = "tags/v4.3.0")

#devtools::install_github("satijalab/seurat-object", ref = "release/4.0.0")
install.packages("Seuratobject")
packageVersion("Seurat")
packageVersion("SeuratObject")

library(Seurat)
library(SeuratData)
library(patchwork)
#RemoveData("ifnb")
options(timeout=600) 

packageVersion("Seurat")
packageVersion("SeuratData")



InstallData("ifnb")
LoadData("ifnb")
head(ifnb@meta.data)

ifnb <- UpdateSeuratObject(ifnb)

ifnb.list <- SplitObject(ifnb, split.by = "stim")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
} )

features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(
  object.list =ifnb.list, anchor.features = features
  )

immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

library(magrittr)

immune.combined <- immune.combined %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.5)


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
  dot.scale = 8,
  split.by = "stim"
) + RotatedAxis()



library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "stim"
avg.t.cells <- as.data.frame(
  log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
)
avg.t.cells$gene <- rownames(avg.t.cells)


cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- as.data.frame(
  log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA)
)
avg.cd14.mono$gene <- rownames(avg.cd14.mono)


genes.to.label <- c(
  "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8"
)


p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) +
  geom_point() +
  ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) +
  geom_point() +
  ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)


p1 + p2



immune.combined$celltype.stim <- paste(
  Idents(immune.combined),
  immune.combined$stim,
  sep = "_"
)


immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"


b.interferon.response <- FindMarkers(
  immune.combined,
  ident.1 = "B_STIM",
  ident.2 = "B_CTRL",
  verbose = FALSE
)


FeaturePlot(
  immune.combined,
  features = c("CD3D", "GNLY", "IFI6"),
  split.by = "stim",
  max.cutoff = 3,
  cols = c("grey", "red")
)



plots <- VlnPlot(
  immune.combined,
  features = c("LYZ", "ISG15", "CXCL10"),
  split.by = "stim",
  group.by = "celltype",
  pt.size = 0,
  combine = FALSE
)


wrap_plots(plots = plots, ncol = 1)



