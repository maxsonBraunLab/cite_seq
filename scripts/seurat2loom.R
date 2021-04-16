library(loomR)
library(Seurat)


output_loom <- snakemake@output[["out_loom"]]#paste(getwd(),snakemake@output[["out_loom"]],sep="/")
seurat_file <- snakemake@input[["seurat_file"]]#paste(getwd(),snakemake@input[["seurat_file"]],sep="/")

seuratObj = readRDS(seurat_file)

#has NA is going to break loom 
#correct any NA to unclassified, could potentially break if column is made of numeric
for (colum in colnames(seuratObj@meta.data)){
    if(any(is.na(seuratObj@meta.data[[colum]]))){
        seuratObj@meta.data[[colum]][which(is.na(seuratObj@meta.data[[colum]]))] <- "unclassified"
    }
}


as.loom(seuratObj, assay = "integrated", filename = output_loom)