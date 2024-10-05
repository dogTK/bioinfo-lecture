#Chapter27
setwd("~/Desktop/R")

install.packages("dplyr")
install.packages("Seurat") # Docker imageを使ってる方はv5がインストールされるためスキップ推奨
install.packages("patchwork")
library(dplyr)
library(Seurat)
library(patchwork)

# Load the dataset
#"./GSE149689/"は、.が現在のディレクトリを表している
pbmc.data <- Read10X(data.dir = "./GSE149689/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(
  counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200
)

# Calculate percent mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Quality control
pbmc <- subset(
  pbmc, nFeature_RNA >= 200 & nFeature_RNA <= 5000 & percent.mt < 15
)

pbmc <- pbmc %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = VariableFeatures(object = pbmc)) %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:20) %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 0.5)

# Visualization
DimPlot(pbmc, reduction = "umap", label = T)

# Check if BiocManager is installed, if not, install it
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required packages using BiocManager
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("SingleCellExperiment")

# Load required libraries
# SingleR version 2.6.0
# celldex version 1.14.0
# SingleCellExperiment version 1.26.0
library(SingleR)
#出力
#Loading required package: SummarizedExperiment
#Loading required package: MatrixGenerics
#Loading required package: matrixStats

#Attaching package: ‘matrixStats’

#The following object is masked from ‘package:dplyr’:
  
  #count


#Attaching package: ‘MatrixGenerics’

#The following objects are masked from ‘package:matrixStats’:
  
  #colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse, colCounts, colCummaxs,
#colCummins, colCumprods, colCumsums, colDiffs, colIQRDiffs, colIQRs, colLogSumExps,
#colMadDiffs, colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds, colSums2,
#colTabulates, colVarDiffs, colVars, colWeightedMads, colWeightedMeans,
#colWeightedMedians, colWeightedSds, colWeightedVars, rowAlls, rowAnyNAs, rowAnys,
#rowAvgsPerColSet, rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps, rowMadDiffs, rowMads,
#rowMaxs, rowMeans2, rowMedians, rowMins, rowOrderStats, rowProds, rowQuantiles,
#rowRanges, rowRanks, rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs,
#rowVars, rowWeightedMads, rowWeightedMeans, rowWeightedMedians, rowWeightedSds,
#rowWeightedVars

#Loading required package: GenomicRanges
#Loading required package: stats4
#Loading required package: BiocGenerics

#Attaching package: ‘BiocGenerics’

#The following object is masked from ‘package:SeuratObject’:
  
  #intersect

#The following objects are masked from ‘package:dplyr’:
  
  #combine, intersect, setdiff, union

#The following objects are masked from ‘package:stats’:
  
  #IQR, mad, sd, var, xtabs

#The following objects are masked from ‘package:base’:
  
  #anyDuplicated, aperm, append, as.data.frame, basename, cbind, colnames, dirname,
#do.call, duplicated, eval, evalq, Filter, Find, get, grep, grepl, intersect,
#is.unsorted, lapply, Map, mapply, match, mget, order, paste, pmax, pmax.int, pmin,
#pmin.int, Position, rank, rbind, Reduce, rownames, sapply, setdiff, table, tapply,
#union, unique, unsplit, which.max, which.min

#Loading required package: S4Vectors

#Attaching package: ‘S4Vectors’

#The following objects are masked from ‘package:dplyr’:
  
  #first, rename

#The following object is masked from ‘package:utils’:
  
  #findMatches

#The following objects are masked from ‘package:base’:
  
  #expand.grid, I, unname

#Loading required package: IRanges

#Attaching package: ‘IRanges’

#The following object is masked from ‘package:sp’:
  
  #%over%
  
  #The following objects are masked from ‘package:dplyr’:
  
  #collapse, desc, slice

#Loading required package: GenomeInfoDb
#Loading required package: Biobase
#Welcome to Bioconductor

#Vignettes contain introductory material; view with 'browseVignettes()'. To cite
#Bioconductor, see 'citation("Biobase")', and for packages 'citation("pkgname")'.


#Attaching package: ‘Biobase’

#The following object is masked from ‘package:MatrixGenerics’:
  
  #rowMedians

#The following objects are masked from ‘package:matrixStats’:
  
  #anyMissing, rowMedians


#Attaching package: ‘SummarizedExperiment’

#The following object is masked from ‘package:Seurat’:
  
  #Assays

#The following object is masked from ‘package:SeuratObject’:
  
  #Assays

library(celldex)
#出力
#Attaching package: ‘celldex’

#The following objects are masked from ‘package:SingleR’:
  
#BlueprintEncodeData, DatabaseImmuneCellExpressionData, HumanPrimaryCellAtlasData,
#ImmGenData, MonacoImmuneData, MouseRNAseqData, NovershternHematopoieticData

library(SingleCellExperiment)

# Load reference data from celldex
ref <- celldex::HumanPrimaryCellAtlasData()

# Run SingleR to infer cell types of pbmc dataset using reference data
results <- SingleR(
  test = as.SingleCellExperiment(pbmc), ref = ref, labels = ref$label.main
)

# Add inferred cell type labels to pbmc object
pbmc$singlr_labels <- results$labels

# Visualize cell types in a UMAP plot with labels
DimPlot(pbmc, reduction = "umap", group.by = "singlr_labels", label = TRUE)