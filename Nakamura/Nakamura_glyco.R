# Load package
library(Seurat)
library(dplyr)
library(UCell)

# data processing
pbmc.data <- Read10X(data.dir = "/path/to/scRNA_test")
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
# Loading glyco list
glyco_list <- read_tsv('GlycogenesList.tsv', col_names = c("Gene", "Pathway"))

# Calculate module score per pathway
for (i in 1:nrow(glyco_list)) {
# Get gene list for current pathway
  genes <- unlist(strsplit(glyco_list$Genes[i], ","))
  
# Calculate module score using AddModuleScore_UCell function
  pbmc <- AddModuleScore_UCell(obj = pbmc, features = list(Pathway = genes), name = paste0("pathway_", glyco_list$Pathway[i]))
}
