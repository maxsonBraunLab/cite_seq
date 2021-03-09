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
library(plotly)
library(RColorBrewer)
library(parallel)


## ----setup_variables----------------------------------------------------------

# stderr to log file
f <- file("logs/cluster.err", "w")
sink(f, type = "message")

# Name of the single cell project
projectName <- snakemake@config$projectName

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

# Find all DE genes for each cluster. This will take a long time.
doFindAllMarkers <- snakemake@config$findAllMarkers



## ----setup_inherent_variables-------------------------------------------------

# Load Seurat objects
exptsList <- readRDS(snakemake@input[[1]])



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

plan("sequential")

## ----umap---------------------------------------------------------------------

for (sample in names(exptsList)) {
  
  paste(sprintf("UMAP Clustering the %s object", sample))
  
  myso <- exptsList[[sample]]
  
  # Cluster the cells over a range of resolutions
  myso <- myso %>%
    FindNeighbors(dims = 1:nDims, verbose = FALSE) %>%
    FindClusters(resolution = clusRes[[sample]], verbose = FALSE) %>%
    RunUMAP(dims = 1:nDims, verbose = FALSE, n.components = 3)
  
  # Plot and save UMAP clusters
  Idents(myso) <- myso$seurat_clusters
  plotList[[sample]][["umap"]] <- DimPlot(object = myso, reduction = 'umap', label = TRUE, label.size = 7) 
  
  # Save clustering
  exptsList[[sample]] <- myso
}

rm(myso)

## -----------------------------------------------------------------------------

# define umap output directory
if (!dir.exists("data/umaps")) {
  dir.create("data/umaps")
}

for (i in names(exptsList)) {

  myso <- exptsList[[i]]

  cell_meta <- merge(myso@reductions$umap@cell.embeddings, myso@meta.data, by = 0, all = TRUE)
  cell_meta <- cell_meta %>% column_to_rownames("Row.names")
  col_scheme <- colorRampPalette(brewer.pal(8, snakemake@config$rcolorbrewer_palette))(length(unique(cell_meta$seurat_clusters)))

  p <- plot_ly(cell_meta,
          x = ~UMAP_1,
          y = ~UMAP_2,
          z = ~UMAP_3,
          size = 1,
          color = ~seurat_clusters,
          colors = col_scheme,
          # hover text
          text = ~paste("Cluster:", seurat_clusters, "<br>nFeature_ADT:", nFeature_ADT, "<br>nCount_SCT:", nCount_SCT)) %>% 
    add_markers() %>%
    layout(title = paste(toupper(i), "UMAP Clusters"),
           xaxis = list(title = "UMAP_1"),
           yaxis = list(title = "UMAP_2"),
           zaxis = list(title = "UMAP_3"))

  # save and print where they are located.
  out_file <- paste0("data/umaps/", tolower(i), "_umap.html")
  print(paste("Sample", i, "UMAP file is located at:", out_file))
  htmlwidgets::saveWidget(p, out_file)
}

rm(myso)


## -----UMAP-by-identity--------------------------------------------------------

if (!dir.exists("data/umaps")) {
  dir.create("data/umaps")
}

myso <- exptsList[["integrated"]]

cell_meta <- merge(myso@reductions$umap@cell.embeddings, myso@meta.data, by = 0, all = TRUE)
cell_meta <- cell_meta %>% column_to_rownames("Row.names")
col_scheme <- colorRampPalette(brewer.pal(8, snakemake@config$rcolorbrewer_palette))(length(unique(cell_meta$seurat_clusters)))

p <- plot_ly(cell_meta,
        x = ~UMAP_1,
        y = ~UMAP_2,
        z = ~UMAP_3,
        size =1,
        color = ~orig.ident,
        colors = col_scheme,
        # hover text
        text = ~paste("Cluster:", seurat_clusters, "<br>nFeature_ADT:", nFeature_ADT, "<br>nCount_SCT:", nCount_SCT)) %>% 
  add_markers() %>%
  layout(title = "Integrated UMAP by Sample Identity",
         xaxis = list(title = "UMAP_1"),
         yaxis = list(title = "UMAP_2"),
         zaxis = list(title = "UMAP_3"))

out_file <- "data/umaps/integrated_identity_umap.html"
print(paste("Integrated UMAP colored by sample identity is here:", out_file))
htmlwidgets::saveWidget(p, out_file)

rm(cell_meta, col_scheme, p, myso, out_file)



## ----cluster_markers----------------------------------------------------------

# parallelize FindAllMarkers. Share snakemake@config$cores across all samples.
if (doFindAllMarkers) {

  if (!dir.exists("data/markers")) {
    dir.create("data/markers")
  }

  allSampleMarkers <- mclapply(names(exptsList), function(x) {
    a <- FindAllMarkers(exptsList[[x]], max.cells.per.ident = 100)
    a
  }, mc.cores = floor(snakemake@config$cores / length(names(exptsList))) )

  names(allSampleMarkers) <- names(exptsList)
}
# mclapply returns results in the right order. https://stackoverflow.com/questions/14697901/is-mclapply-guaranteed-to-return-its-results-in-order

# make plots and export tables
if (doFindAllMarkers) {
  plotList <- list()

  for (sample in names(allSampleMarkers)) {
    # j is all DE genes per cluster for an individual sample
    j <- allSampleMarkers[[sample]]
    top <- j %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
    plotList[[sample]][['topMarkerHeatmap']] <- DoHeatmap(exptsList[[sample]], features = top$gene) + NoLegend()
    # export allMarkers for each sample and heatmap
    write.table(j, file = paste0("data/markers/", sample, "-markers.txt"), row.names = TRUE, sep = "\t", quote = FALSE)
    png(paste0("data/markers/", sample, "_heatmap.png"), width = 16, height = 9, units = "in", res = 300)
    plotList[[sample]][['topMarkerHeatmap']]
    dev.off()
  }
}

## ----save_data----------------------------------------------------------------
print(  sprintf("Saving preprocessed individual samples and the integrated object in %s", snakemake@output[[1]]) )
saveRDS(exptsList, file = snakemake@output[[1]])

## ----session_info-------------------------------------------------------------
sessionInfo()

