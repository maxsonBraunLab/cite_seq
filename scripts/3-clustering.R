## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)


## ----setup_libs---------------------------------------------------------------
library(Seurat)
library(data.table)
library(ggplot2)
library(knitr)
library(networkD3)
library(patchwork)
library(dplyr)
library(clustree)
library(ggthemes)
library(viridis)
library(yaml)
library(tibble)
library(future)


## ----setup_variables----------------------------------------------------------

# Name of the single cell project
projectName <- snakemake@config$projectName

doc_title <- paste(snakemake@config$title, "- clustering")
author_list <- paste(snakemake@config$authors, collapse = ", ")

# Clustering Resolution

clusRes <- lapply(names(snakemake@config$resolutions), function(x) {snakemake@config$resolutions[[x]]}) # resolution values (floats)
names(clusRes) <- names(snakemake@config$resolutions) # map sample name to resolution

# Number of dimensions to reduce on
nDims <- snakemake@config$nDims

# check if resolutions > 2
if (any(clusRes > 2)) {
  stop("Cluster resolution is > 2")
}

# define cores
plan("multiprocess", workers = snakemake@threads) # no need to set seed despite random number generator warnings.

# Create plot that explains variance by Principle Components (PCs) to find the dimensionality
percentExplained <- function( seuratObject ){
  return(cumsum(seuratObject@reductions$pca@stdev^2)/seuratObject@reductions$pca@misc$total.variance)
}

plotPercentExplained <- function( seuratObject , ndims = nDims) {
    data.use <- percentExplained( seuratObject )
    plot <- ggplot(data = data.frame(dims = 1:ndims, var = data.use[1:ndims])) +
            geom_point(mapping = aes_string(x = 'dims', y = 'var')) + labs( y = "% variance explained", x = "PC") +
            ggtitle(sprintf("Variance Explained by Principle Components in %s", sample)) +
            theme_clean()
    return(plot)
}


## ----setup_inherent_variables-------------------------------------------------

# Load Seurat objects
exptsList <- readRDS(snakemake@input[[1]])

# Clustered variables list. Each SO contains cluster information.
clusExptsList <- list()


## ----id_umap_vars-------------------------------------------------------------
# Perform linear dimensional reduction (UMAP) on myso object
plotList <- list()

for (sample in names(exptsList)) {
  
  myso <- exptsList[[sample]]
  
  sampleAssay <- ifelse(sample=='integrated', "Integrated_snn_res." ,"SCT_snn_res.")
  
  plotList[[sample]][["plotPercentExplained"]] <- plotPercentExplained(myso, ndims = nDims)
  
  # Cluster the cells over a range of resolutions
  myso <- FindNeighbors(object = myso, dims = 1:nDims, verbose = FALSE)
  myso <- myso %>%
    FindClusters(resolution = seq(0.1, 2, 0.1), verbose = FALSE) %>%
    RunUMAP(dims = 1:nDims, verbose = FALSE)
  
  # Add UMAP embedding
  embeddings <- Embeddings(myso)
  myso@meta.data$PC_1 <- embeddings[,1]
  myso@meta.data$PC_2 <- embeddings[,2]
  
  # Clustree visualizations to determine best cluster resolution
  plotList[[sample]][["clustree"]] <- clustree(myso, prefix = sampleAssay) 
  
  plotList[[sample]][["clustreeOverlay"]] <- clustree_overlay(x = myso[[]],
                   prefix = sampleAssay,
                   x_value = "PC_1",
                   y_value = "PC_2") 
  
}

rm(myso)



## ----results='asis'-----------------------------------------------------------
# for (i in names(exptsList)) {
#   cat(paste( "##", toupper(i), "\n\n" ))

#   cat("### Percent Variance Explained by PC\n")
#   plot(plotList[[i]][["plotPercentExplained"]])
#   cat("\n\n---\n\n")

#   cat("### Clustree Hierarchy\n")
#   plot(plotList[[i]][["clustree"]])
#   cat("\n\n---\n\n")
 
#   cat("### Clustree Overlay\n")
#   plot(plotList[[i]][["clustreeOverlay"]])
#   cat("\n\n---\n\n")
# }


## ----umap---------------------------------------------------------------------

for (sample in names(exptsList)) {
  
  paste(sprintf("UMAP Clustering the %s object", sample))
  
  myso <- exptsList[[sample]]
  
  # Cluster the cells over a range of resolutions
  myso <- myso %>%
    FindNeighbors(dims = 1:nDims, verbose = FALSE) %>%
    FindClusters(resolution = clusRes[[sample]], verbose = FALSE) %>%
    RunUMAP(dims = 1:nDims, verbose = FALSE)
  
  # Plot and save UMAP clusters
  Idents(myso) <- myso$seurat_clusters
  plotList[[sample]][["umap"]] <- DimPlot(object = myso, reduction = 'umap', label = TRUE, label.size = 7) 
  
  # Save clustering
  clusExptsList[[sample]] <- myso
}

rm(myso)


## ----results='asis'-----------------------------------------------------------
# for (i in names(exptsList)) {
#   cat(paste( "##", toupper(i), "\n" ))

#   cat("\n")
#   plot(plotList[[i]][["umap"]])
#   cat("\n\n")
# }

## ----save_data----------------------------------------------------------------
saveRDS(clusExptsList, file = snakemake@output[[1]])

## ----session_info-------------------------------------------------------------
sessionInfo()

