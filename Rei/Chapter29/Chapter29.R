#Chapter29
#下記のディレクトリ構造になるようにデスクトップにdataフォルダを作成する。
# pbmc4k
"./data/pbmc4k/filtered_gene_bc_matrices/GRCh38/"

# pbmc8k
"./data/pbmc8k/filtered_gene_bc_matrices/GRCh38/"

# pbmc3k
"./data/pbmc3k/filtered_gene_bc_matrices/hg19/"

#ディレクトリの設定
setwd("~/Desktop/data")

#libraryの読み込み
library(dplyr)
library(Seurat)
library(patchwork)

#データ解析を行うために、Seuratオブジェクトを作成する。
#Seuratオブジェクトは、scRNA-seqデータのすべての解析結果をひとつにまとめるため、ユーザーは個別にデータを保持したり操作する必要がなく、統一的に管理しやすい構造
pbmc4k.data <- Read10X(data.dir = "./data/pbmc4k/filtered_gene_bc_matrices/GRCh38/")
pbmc4k <- CreateSeuratObject(counts = pbmc4k.data, project = "PBMC4K")
pbmc4k

#pbmc4kの出力結果
#An object of class Seurat 
#33694 features across 4340 samples within 1 assay 
#Active assay: RNA (33694 features, 0 variable features)
#1 layer present: counts

pbmc8k.data <- Read10X(data.dir = "./data/pbmc8k/filtered_gene_bc_matrices/GRCh38/")
pbmc8k <- CreateSeuratObject(counts = pbmc8k.data, project = "PBMC8K")
pbmc8k

#pbmc8kの出力結果
#An object of class Seurat 
#33694 features across 8381 samples within 1 assay 
#Active assay: RNA (33694 features, 0 variable features)
#1 layer present: counts

#merge関数を使って、pbmc4kとpbmc8kのデータセットを統合する。
pbmc.combined <- merge(
  pbmc4k,
  y = pbmc8k,
  add.cell.ids = c("4K", "8K"),
  project = "PBMC12K"
)
pbmc.combined

#mergeさせたデータセットは,pbmc.combinedに入っている。
#pbmc.combinedの出力結果
#An object of class Seurat 
#33694 features across 12721 samples within 1 assay 
#Active assay: RNA (33694 features, 0 variable features)
#2 layers present: counts.PBMC4K, counts.PBMC8K

#head関数は、データの先頭部分を表示
#tail関数は、データの末尾部分を表示
#デフォルトでは、６行表示する
#head関数で4K_から始まる配列、tail関数8K_から始まる配列が表示されればOK
head(colnames(pbmc.combined))
tail(colnames(pbmc.combined))

#pbmc4k と pbmc8k の数を集計
#集計するときはtable関数を使う
table(pbmc.combined$orig.ident)

#結果の出力
#PBMC4K=4340, PBMC8K=8381

#3つ目のデータセットを読み込む。
pbmc3k.data <- Read10X(
  data.dir = "./data/pbmc3k/filtered_gene_bc_matrices/hg19/"
)
pbmc3k <- CreateSeuratObject(counts = pbmc3k.data, project = "PBMC3K")
pbmc3k

#pbmc3kの出力結果。(scRNA-seq本誤植あり、本文中にはpbmc4kの出力と記載されているがpbmc3k)
#An object of class Seurat 
#32738 features across 2700 samples within 1 assay 
#Active assay: RNA (32738 features, 0 variable features)
#1 layer present: counts

#y = c(pbmc4k, pbmc8k):
#ここで pbmc4k（4,000細胞）と pbmc8k（8,000細胞）の2つのSeuratオブジェクトを指定し、これらを基準オブジェクト pbmc3k に統合
#add.cell.ids = c("3K", "4K", "8K"):
#細胞ごとに元データセットを識別できるように、それぞれのデータセットにラベル
#project = "PBMC15K":
#統合後のSeuratオブジェクトにプロジェクト名として "PBMC15K" を指定
pbmc.big <- merge(
  pbmc3k,
  y = c(pbmc4k, pbmc8k),
  add.cell.ids = c("3K", "4K", "8K"),
  project = "PBMC15K"
)
pbmc.big

#pbmc.bigの出力結果
#An object of class Seurat 
#36116 features across 15421 samples within 1 assay 
#Active assay: RNA (36116 features, 0 variable features)
#3 layers present: counts.PBMC3K, counts.PBMC4K, counts.PBMC8K

#unique関数
#ベクトルやデータフレーム、リスト内の重複している要素を取り除き、一意の要素だけを返すために使用。
#colnames(pbmc.big)=>pbmc.big オブジェクトの列名（細胞名）を取得
#strsplit(colnames(pbmc.big), split = "_")=>"_"で分割(例："3K_cell001"なら、"3K"と"cell001")
#sapply(..., FUN = "[", 1)=>sapply 関数を使って、分割された各リストから最初の要素（1番目の要素）を抽出,FUN = "[" はリストから特定のインデックスの要素を取り出すためのインデックス操作を意味する。1 が指定されているので、各リストの1番目の要素（この場合、"3K", "4K" など）を取り出す。
unique(sapply(X = strsplit(colnames(pbmc.big), split = "_"), FUN = "[", 1))
#uniqueの出力結果
#[1]"3K" "4K" "8K"

table(pbmc.big$orig.ident)
#tableの出力結果
#PBMC3K=2700, PBMC4K=4340, PBMC8K=8381

#pbmc4kとpbmc8kをNormalizeData関数によって正規化した後にそれぞれを結合
pbmc4k <- NormalizeData(pbmc4k)
pbmc8k <- NormalizeData(pbmc8k)
pbmc.normalized <- merge(
  pbmc4k,
  y = pbmc8k,
  add.cell.ids = c("4K", "8K"),
  project = "PBMC12K",
  merge.data = TRUE
)

#GetAssayData 関数は特定のアッセイの特定のデータを閲覧することができる関数でseuratオブジェクトの中のデータを見ることができる。
GetAssayData(pbmc.normalized)[1:10, 1:15]



