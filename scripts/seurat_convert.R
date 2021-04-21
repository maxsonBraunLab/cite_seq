

library(Seurat)
is_installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

#seurat abandoned loomR for seuratDisk
#https://mojaveazure.github.io/seurat-disk/
if (!is_installed("remotes")) {
  install.packages("remotes",repos = "http://cran.us.r-project.org")
}
if (!is_installed('SeuratDisk')){
	remotes::install_github("mojaveazure/seurat-disk")	
}
library(SeuratDisk)

#snakemake parameters
output_file <- snakemake@output[["out_file"]]
seurat_file <- snakemake@input[["seurat_file"]]
embedding_file <- snakemake@output[["embedding_file"]]

seuratObj = readRDS(seurat_file)


#save as loom
as.loom(seuratObj[["integrated"]], filename = output_file)

#save umap embedding
embeddings <- Embeddings(object = seuratObj[["integrated"]], reduction = 'umap')
write.table(embeddings, file = embedding_file, sep = "\t", row.names=T, quote=F,col.names = F)

