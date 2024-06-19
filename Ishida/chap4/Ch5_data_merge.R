


pbmc4k.data <- Read10X(data.dir = "./data/pbmc4k/filtered_gene_bc_matrices/GRCh38/")
pbmc4k <- CreateSeuratObject(counts = pbmc4k.data, project = "PBMC4K")
pbmc4k

pbmc8k.data <- Read10X(data.dir = "./data/pbmc8k/filtered_gene_bc_matrices/GRCh38/")
pbmc8k <- CreateSeuratObject(counts = pbmc8k.data, project = "PBMC8K")
pbmc8k

pbmc.combined <- merge(
  pbmc4k,
  y = pbmc8k,
  add.cell.ids = c("4K", "8K"),
  project = "PBMC12K"
)
pbmc.combined

head(colnames(pbmc.combined))
tail(colnames(pbmc.combined))

table(pbmc.combined$orig.ident)

pbmc3k.data <- Read10X(
  data.dir = "./data/pbmc3k/filtered_gene_bc_matrices/hg19/"
)
pbmc3k <- CreateSeuratObject(counts = pbmc3k.data, project = "PBMC3K")
pbmc3k

pbmc.big <- merge(
  pbmc3k,
  y = c(pbmc4k, pbmc8k),
  add.cell.ids = c("3K", "4K", "8K"),
  project = "PBMC15K"
)
pbmc.big

unique(sapply(X = strsplit(colnames(pbmc.big), split = "_"), FUN = "[", 1))

table(pbmc.big$orig.ident)

pbmc4k <- NormalizeData(pbmc4k)
pbmc8k <- NormalizeData(pbmc8k)
pbmc.normalized <- merge(
  pbmc4k,
  y = pbmc8k,
  add.cell.ids = c("4K", "8K"),
  project = "PBMC12K",
  merge.data = TRUE
)

pbmc.normalized

Assays(pbmc.normalized)

data_matrix <- GetAssayData(object = pbmc.normalized, assay = "RNA", layer = "data")


data_matrix[1:10, 1:15]


GetAssayData(pbmc.normalized)[1:10, 1:15]
