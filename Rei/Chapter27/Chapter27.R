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
#出力
  #Normalizing layer: counts
  #Performing log-normalization
  # 0%   10   20   30   40   50   60   70   80   90   100%
  # [----|----|----|----|----|----|----|----|----|----|
  # **************************************************|
     #Finding variable features for layer counts
  #Calculating gene variances
  # 0%   10   20   30   40   50   60   70   80   90   100%
  # [----|----|----|----|----|----|----|----|----|----|
  # **************************************************|
  # Calculating feature variances of standardized and clipped values
  # 0%   10   20   30   40   50   60   70   80   90   100%
  # [----|----|----|----|----|----|----|----|----|----|
  # **************************************************|
  #       Centering and scaling data matrix
  #       |=========================================================================================| 100%
  #       PC_ 1 
  #       Positive:  IL32, CCL5, CD3E, GZMM, TRBC2, CD69, CST7, CD247, CTSW, RORA 
  #       NKG7, GZMA, LCK, CD7, SPOCK2, GZMH, C12orf75, CD3D, HOPX, IFITM1 
  #       PRF1, CD2, GZMB, KLRD1, TRBC1, FGFBP2, CD3G, GNLY, ISG20, TRAC 
  #       Negative:  FCN1, LYZ, AIF1, SERPINA1, CTSS, S100A9, LST1, TYMP, CYBB, SPI1 
  #       GRN, CSTA, CST3, VCAN, S100A8, MNDA, CD14, AC020656.1, MS4A6A, NCF2 
  #       BRI3, PSAP, FGL2, S100A11, MPEG1, TKT, FTL, CFD, CFP, MAFB 
  #       PC_ 2 
  #       Positive:  HSP90AA1, NKG7, S100A4, GZMA, IL32, CST7, GZMM, GZMH, PRF1, PPIB 
  #       ANXA1, GZMB, CD247, IFITM1, FGFBP2, HSPA5, GNLY, CD3E, KLRD1, CTSW 
  #       HOPX, CD7, S100A10, CALR, PSME2, LCK, S100A6, TRBC2, RORA, DDIT4 
  #       Negative:  CAVIN2, GNG11, PF4, TUBB1, GP9, ACRBP, HIST1H2AC, CMTM5, TMEM40, PPBP 
  #       RGS18, NRGN, PTCRA, CLEC1B, TSC22D1, MTURN, PRKAR2B, TREML1, C2orf88, SPARC 
  #       RUFY1, ESAM, LIMS1, HIST1H3H, MPP1, BEX3, HIST1H2BJ, RGS10, ITGA2B, PGRMC1 
  #       PC_ 3 
  #       Positive:  CCL5, NKG7, CST7, GZMA, CTSW, PRF1, GZMH, GZMB, FGFBP2, GZMM 
  #       GNLY, KLRD1, S100A4, IL32, ANXA1, CD247, HOPX, CD7, ACTB, SPON2 
  #       GAPDH, C12orf75, CCL4, MATK, CLIC3, CD3E, CD63, LCK, ID2, ACTG1 
  #       Negative:  CD79A, MS4A1, IGHM, LINC00926, TNFRSF13C, IGHD, BANK1, CD79B, RALGPS2, AFF3 
  #       CD22, VPREB3, TCL1A, FAM129C, BLK, FCER2, PAX5, FCRL1, IGKC, TCF4 
  #       GNG7, BCL11A, HLA-DQB1, HLA-DRA, P2RX5, SPIB, CD83, CD74, ADAM28, CD40 
  #       PC_ 4 
  #       Positive:  S100A12, IL7R, SLC25A37, S100A8, MAL, ALAS2, HBA2, HBB, IL1R2, IRS2 
  #       RCAN3, MCEMP1, VCAN, HBA1, AHSP, GCA, SELL, PLBD1, CA1, LEF1 
  #       SLC4A1, CYP1B1, CXCL8, SOCS3, HBM, S100A9, SLC2A3, MEGF9, TRABD2A, CKAP4 
  #       Negative:  HLA-DPA1, HLA-DPB1, FCGR3A, HLA-DRB1, RHOC, HLA-DQA1, CD74, HLA-DRA, MTSS1, CTSC 
  #       ABI3, GZMB, CDKN1C, HLA-DRB5, HLA-DQB1, LY6E, HLA-DMA, PRF1, HES4, FGFBP2 
  #       CSF1R, KLRD1, GNLY, SPON2, PSME2, CALR, HLA-DQA2, CTSW, KLRF1, CD79B 
  #       PC_ 5 
  #       Positive:  CDKN1C, IL7R, CSF1R, MAL, TCF7L2, CD4, RCAN3, NEURL1, TRAC, CAMK1 
  #       BATF3, CD3D, CTSL, ZNF703, HES4, LEF1, CAMK4, TRAT1, CKB, CD3G 
  #       LRRC25, MS4A7, CD28, LTB, TCF7, SIGLEC10, CD5, TRABD2A, C1QA, COTL1 
  #       Negative:  GZMB, KLRF1, CCL4, SPON2, PRF1, FGFBP2, KLRD1, GNLY, ALOX5AP, CXCL8 
  #       CCL3, S100A12, NKG7, CLIC3, ITGAM, TRDC, S1PR5, NFKBIA, CST7, MCEMP1 
  #       CTSW, PLAC8, TTC38, NCF1, IGFBP7, ADGRG1, IER2, DYSF, HOPX, MYOM2 
  #       Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
  #       To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
  #       This message will be shown once per session
  #       00:55:19 UMAP embedding parameters a = 0.9922 b = 1.112
  #       00:55:19 Read 59176 rows and found 20 numeric columns
  #       00:55:19 Using Annoy for neighbor search, n_neighbors = 30
  #       00:55:19 Building Annoy index with metric = cosine, n_trees = 50
  #       0%   10   20   30   40   50   60   70   80   90   100%
  #       [----|----|----|----|----|----|----|----|----|----|
  #       **************************************************|
  #          00:55:26 Writing NN index file to temp file /var/folders/8t/mclc_3wx2b784wj_ns28_vmm0000gn/T//Rtmpb4MnGp/file7306d36f915
  #          00:55:26 Searching Annoy index using 1 thread, search_k = 3000
  #          00:55:49 Annoy recall = 100%
  #          00:55:49 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
  #          00:55:51 Initializing from normalized Laplacian + noise (using RSpectra)
  #          00:55:54 Commencing optimization for 200 epochs, with 2548648 positive edges
  #          Using method 'umap'
  #          0%   10   20   30   40   50   60   70   80   90   100%
  #          [----|----|----|----|----|----|----|----|----|----|
  #          **************************************************|
  #             00:56:25 Optimization finished
  #             Computing nearest neighbor graph
  #             Computing SNN
  #             Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
               
  #             Number of nodes: 59176
  #             Number of edges: 1862658
               
  #             Running Louvain algorithm...
  #             0%   10   20   30   40   50   60   70   80   90   100%
  #             [----|----|----|----|----|----|----|----|----|----|
  #             **************************************************|
  #             Maximum modularity in 10 random starts: 0.9291
  #             Number of communities: 21
  #             Elapsed time: 19 seconds

# Visualization
DimPlot(pbmc, reduction = "umap", label = T)