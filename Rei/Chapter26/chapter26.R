#Chapter26
#必要なパッケージのインストールとライブラリの読み込み
#Rのバージョン4.4.1
#dplyrのバージョン1.1.4
#Seuratのバージョン5.1.0
#patchworkのバージョン1.3.0
install.packages("dplyr")
install.packages("Seurat") # Docker imageを使ってる方はv5がインストールされるためスキップ推奨
install.packages("patchwork")

library(dplyr)
#出力
  #Attaching package: ‘dplyr’

  #The following objects are masked from ‘package:stats’:
  
    #filter, lag

  #The following objects are masked from ‘package:base’:
  
    #intersect, setdiff, setequal, union
library(Seurat)
#出力
  #Loading required package: SeuratObject
  #Loading required package: sp
  #‘SeuratObject’ was built under R 4.4.0 but the current version is 4.4.1; it is recomended
  #that you reinstall ‘SeuratObject’ as the ABI for R may have changed

  #Attaching package: ‘SeuratObject’

  #The following objects are masked from ‘package:base’:
  
    #intersect, t

library(patchwork)

setwd("~/Desktop/R")

# サンプルを読み出す
pbmc_v2.data <- Read10X(data.dir = "./GSE139386_RAW/")

# バイオリンプロットを表示する
pbmc_v2 <- CreateSeuratObject(
  counts = pbmc_v2.data,
  project = "pbmc_v23k",
  min.cells = 3,
  min.features = 200
)
#出力
  #Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
pbmc_v2[["percent.mt"]] <- PercentageFeatureSet(pbmc_v2, pattern = "^MT-")
VlnPlot(
  pbmc_v2,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3
)
#出力
  #Warning: Default search for "data" layer in "RNA" assay yielded no results; utilizing "counts" layer instead.