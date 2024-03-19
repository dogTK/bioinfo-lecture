# Load package
library(Seurat)
library(dplyr)


# data processing
pbmc.data <- Read10X(data.dir = "/Users/idrc/Desktop/singlecell_test/scRNA_test")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


pbmc <- pbmc %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10) %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 0.5)

#AddModuleScore
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$Cell_type <- pbmc@active.ident
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pbmc <- AddModuleScore(object = pbmc,
                       features = list(c("CD8A","GZMA")),
                       name = "CD8_Tcell_score")
#Feature Plotting Scores
FeaturePlot(pbmc,"CD8_Tcell_score1")

#multi geneset
pbmc <- AddModuleScore(object = pbmc,
                       features = list(c("PTPRC","CD3E","CD3D","CD8A","GZMA")),
                       name = "CD8_Tcell_score")
FeaturePlot(pbmc,"CD8_Tcell_score1")		       

pbmc <- AddModuleScore(object = pbmc,
                       features = list(c("CD8A","GZMA"),
                                       c("MS4A1","CD19"),
                                       c("CD14","CD68")),
                       name = c("CD8_Tcell_score","Bcell_score","Mono_score"))

FeaturePlot(pbmc, features = c("CD8_Tcell_score1","Bcell_score2","Mono_score3"), ncol = 3)

