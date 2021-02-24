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
library(circlize)
library(ComplexHeatmap)
library(ggthemes)
library(viridis)
library(yaml)
library(tibble)


## ----setup_variables----------------------------------------------------------

config_file <- read_yaml("config.yaml")

# Samples to process named by their treatment conditions
samples2process <- config_file$samples2process

# Name of the single cell project
projectName <- config_file$projectName

# Date that Seurat object with pre-processed data was generated
importRDS <- list.files("data/rda", pattern="*integrated*", full.names=TRUE)[1] # this is temporary workaround. This will change to the output of preprocess once in snakemake format.

# Clustering Resolution
clusRes <- list()
clusRes[[samples2process[1]]] <- 0.9
clusRes[[samples2process[2]]] <- 0.8
clusRes[[samples2process[3]]] <- 1.0
clusRes[[samples2process[4]]] <- 0.6
clusRes[['integrated']] <- 0.5

# Number of dimensions to reduce on
nDims <- 150


## ----setup_inherent_variables-------------------------------------------------

# Specify directory paths
directory = list(raw = "./data/raw",
                 rda = "./data/rda")

# Load Seurat objects
exptsList <- readRDS(sprintf("%s/integrated.%s.%s.rds", directory$rda, projectName,importDate))

# Clustered variables list
clusExptsList <- list()


## ----id_umap_vars-------------------------------------------------------------
# Perform linear dimensional reduction (UMAP) on myso object
plotList <- list()

for (sample in names(exptsList)) {
  
  myso <- exptsList[[sample]]
  
  sampleAssay <- ifelse(sample=='integrated', "Integrated_snn_res." ,"SCT_snn_res.")
  
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
  
  plotList[[sample]][["plotPercentExplained"]] <- plotPercentExplained(myso)
  
  # Cluster the cells over a range of resolutions
  myso <- FindNeighbors(object = myso, dims = 1:nDims, verbose = FALSE)
  myso <- myso %>%
    FindClusters(resolution = seq(0,2,0.1), verbose = FALSE) %>%
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


## ----results='asis'-----------------------------------------------------------
for (sample in samples2process) {
  cat( paste(rep("#", 2), toupper(sample)) )
  plot(plotList[sample][["plotPercentExplained"]])
}


## -----------------------------------------------------------------------------
plot(plotList[[names(exptsList)[i]]][["plotPercentExplained"]])


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["clustree"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["clustreeOverlay"]]


## -----------------------------------------------------------------------------
plot(plotList[[names(exptsList)[i]]][["plotPercentExplained"]])


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["clustree"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["clustreeOverlay"]]


## -----------------------------------------------------------------------------
plot(plotList[[names(exptsList)[i]]][["plotPercentExplained"]])


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["clustree"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["clustreeOverlay"]]


## -----------------------------------------------------------------------------
plot(plotList[[names(exptsList)[i]]][["plotPercentExplained"]])


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["clustree"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["clustreeOverlay"]]


## -----------------------------------------------------------------------------
plot(plotList[[names(exptsList)[i]]][["plotPercentExplained"]])


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["clustree"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["clustreeOverlay"]]


## ----umap---------------------------------------------------------------------
plostList <- list()

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


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["umap"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["umap"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["umap"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["umap"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["umap"]]


## ----save_data----------------------------------------------------------------
file2save <- sprintf("clustered.%s.%s.rds", projectName, Sys.Date())
print(sprintf("Saving clustered data by individual samples in %s", file2save))
saveRDS(clusExptsList, file = paste(directory$rda, file2save, sep = "/"))


## ----session_info-------------------------------------------------------------
sessionInfo()

