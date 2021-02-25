## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)


## ----setup_libs---------------------------------------------------------------
library(dplyr)
library(Seurat)
library(CHETAH)
library(SingleCellExperiment)
library(circlize)
library(ComplexHeatmap)


## ----setup_variables, include=TRUE--------------------------------------------

# Samples to process named by their treatment conditions
samples2process = c("dmso", "quzi", "ory", "combo")

# Name of the single cell project
projectName <- "flt3_aml"

# Date that Seurat object with pre-processed data was generated
importDate <- "2021-02-17"

# The FindAllMarkers step takes a long time
# Change to TRUE if you want to perform this step
doFindAllMarkers <- FALSE


## ----setup_inherent_variables, include=TRUE-----------------------------------

# Specify directory paths
directory = list(raw = "./data/raw",
                 rda = "./data/rda",
                 reference = "./data/reference",
                 clusterMarkers = "./results/cluster_markers")

# Load Seurat objects
exptsList <- readRDS(sprintf("%s/clustered.%s.%s.rds", directory$rda, projectName,importDate))

# Load CHETAH single-cell reference (https://figshare.com/s/aaf026376912366f81b6)
load(sprintf("%s/CHETAH_TME_reference.Rdata", directory$reference), CHETAH_TME_reference <- new.env())
chetahRef <- CHETAH_TME_reference$reference


## ----ref_set_up---------------------------------------------------------------
# Set up reference using Bioconductor Vignette (Step #7)
# http://www.bioconductor.org/packages/devel/bioc/vignettes/CHETAH/inst/doc/CHETAH_introduction.html 

# Normalize
assay(chetahRef, "counts") <- apply(assay(chetahRef, "counts"), 2, function(column) log2((column/sum(column) * 100000) + 1))

# Discard Ribosomal proteins, which have a high drop-out rate
ribo <- read.table("./data/reference/ribosomal.txt", header = FALSE, sep = '\t')
chetahRef <- chetahRef[!rownames(chetahRef) %in% ribo[,1], ]

# Correlate reference database cell type identifiers by genes with the highest fold-change 
# If the reference is good, all types will correlate poorly or even better, will anti-correlate.
CorrelateReference(ref_cells = chetahRef)

# Uses the reference to classify the reference itself
# If CHETAH works well with the reference, there should be almost no mix-ups in the classification
ClassifyReference(ref_cells = chetahRef)


## ----CHETAH-------------------------------------------------------------------
plotList <- list()

for (sample in names(exptsList)) {
  myso <- exptsList[[sample]]
  
  # Load the Seurat object as an SCE object, which is compatible with CHETAH
  mysoSCE <- as.SingleCellExperiment(myso)
  
  # Run the CHETAH classifier
  mysoSCE <- CHETAHclassifier(input = mysoSCE,
                              ref_cells = chetahRef,
                              thresh = 0.05)
  
  # Plot the cell type assignments and color the non-intermediary cell types
  plotList[[sample]][["CHETAH"]] <- PlotCHETAH(mysoSCE, redD='UMAP', return = TRUE)
  
  # Plot the cell type assignments and color the intermediary cell types
  plotList[[sample]][["CHETAHInts"]] <-PlotCHETAH(mysoSCE, redD='UMAP', interm = TRUE, return = TRUE)
  
  # Add the CHETAH results to the Seurat object metadata
  myso <- AddMetaData(object = myso, metadata=mysoSCE$celltype_CHETAH, col.name="cellTypeChetah")
  
  # Identify the main cell type for each cluster and assign it to the cluster
  nClusters <- max(as.numeric(unlist(unique(myso[['seurat_clusters']]))))
  seuratClustersChetahDF <- data.frame(seuratClustersChetah=character(0))

  # Identify the most prevalent cell type in each cluster and store in new metadat column
  for (i in seq(0,nClusters - 1)) {
    clusterCells <- myso[[]][(myso[['seurat_clusters']] == i), ]
    cellNames <- rownames(clusterCells)
    maxCellType <- names(sort(summary(as.factor(clusterCells$cellTypeChetah), decreasing=T)[1]))
    clusterDF <- data.frame(seuratClustersChetah=rep(maxCellType,length(cellNames)))
    rownames(clusterDF) <- cellNames
    seuratClustersChetahDF <- rbind(seuratClustersChetahDF,clusterDF)
  }
  
   myso <- AddMetaData(object = myso, metadata=seuratClustersChetahDF, col.name="seuratClustersChetah")
   
   plotList[[sample]][['seuratClustersChetah']] <- DimPlot(myso, group.by = 'seuratClustersChetah')
   
   exptsList[[sample]] <- myso
}


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CHETAH"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CHETAHInts"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['seuratClustersChetah']] + labs(title = "UMAP Clusters Annotated by CHETAH Cell Types")


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CHETAH"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CHETAHInts"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['seuratClustersChetah']] + labs(title = "UMAP Clusters Annotated by CHETAH Cell Types")


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CHETAH"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CHETAHInts"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['seuratClustersChetah']] + labs(title = "UMAP Clusters Annotated by CHETAH Cell Types")


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CHETAH"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CHETAHInts"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['seuratClustersChetah']] + labs(title = "UMAP Clusters Annotated by CHETAH Cell Types")


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CHETAH"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CHETAHInts"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['seuratClustersChetah']] + labs(title = "UMAP Clusters Annotated by CHETAH Cell Types")


## ----stage_markers------------------------------------------------------------

for (sample in names(exptsList)) {
  
  myso <- exptsList[[sample]]
  
  stageMarkers <- list(early = c("SOX4"),
                       mid = c("CEBPD"),
                       late = c("CD14"))
  
  for (smarker in stageMarkers) {
    plotList[[sample]][[smarker]] <- FeaturePlot(myso, features = smarker)
  }

}


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["SOX4"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CEBPD"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CD14"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["SOX4"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CEBPD"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CD14"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["SOX4"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CEBPD"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CD14"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["SOX4"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CEBPD"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CD14"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["SOX4"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CEBPD"]]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][["CD14"]]


## ----vG_markers---------------------------------------------------------------

# Load cell type markers reported by van Galen et al. 2019 https://doi.org/10.1016/j.cell.2019.01.031
undiff <- c("MEIS1", "EGR1", "MSI2", "CD38", "CD34", "PROM1", "EGFL7")
gmp <- c("MPO", "ELANE", "CTSG", "AZU1", "LYST")
promono <- c("LYST", "LYZ", "CEBPD", "MNDA")
mono <- c("FCER1G", "FCN1", "CD14", "C5AR1")
cDC <- c("C5AR1", "CLEC4A", "CLEC10A")
pDC <- c("FCER1A", "CLEC4C", "PTPRS", "IRF8", "TCF4")
early_ery <- c("CSF1", "KIT", "HBB", "HBD")
late_ery <- c("CSF1", "KIT", "HBB", "HBD", "GYPA")
pro_b <- c("CD24", "EBF1", "MME", "VPREB1", "PAX5")
b <- c("PAX5", "CD19", "CD79A", "MS4A1", "BANK1")
plasma <- c("MZB1", "IGLL5", "JCHAIN")
t <- c("CD3D", "CD3G", "IL32", "IL7R", "TCF7", "CCL5", "GZMK", "CD8A", "KLRB1")
ctl <- c("KLRB1")
nk <- c("GZMB", "NCAM1")
blasts <- c("ITGAM", "ANPEP", "CD14", "CD33", "FCGR1A", "IL3RA")

markerCats <- factor(c(rep("", 7), 
                       rep("Myeloid", 18), 
                       rep("Erythroid", 5), 
                       rep("Lymphoid", 23), 
                       rep("Myeloid Blasts", 6)),
                     levels=c("", "Myeloid", "Erythroid", "Lymphoid", "Myeloid Blasts"))

cellTypeMarkers <- c(undiff, gmp, promono, mono, cDC, pDC, early_ery, late_ery, pro_b, b, plasma, t, ctl, nk, blasts)
cellTypeMarkers <- cellTypeMarkers[!duplicated(cellTypeMarkers)]


## ----vg_heatmap_dmso, message=FALSE-------------------------------------------
for (sample in names(exptsList)) {
  myso <- exptsList[[sample]]
  
  Idents(myso) <- "seurat_clusters"
  clusterAvg <- as.data.frame(AverageExpression(object = myso, assays = 'SCT'))
  clusterAvg <- subset(clusterAvg, rownames(clusterAvg) %in% cellTypeMarkers)
  
  missing <- cellTypeMarkers[!(cellTypeMarkers %in% rownames(clusterAvg))]
  
  print(sprintf("The following markers are not present in %s and are set to zero: %s",sample, paste(missing, collapse=", ")))
   
  for (marker in missing) {
    clusterAvg[nrow(clusterAvg) + 1, ] <- 0
    rownames(clusterAvg)[nrow(clusterAvg)] <- marker
  }
  
  idx <- match(cellTypeMarkers, rownames(clusterAvg))
  idx <- idx[!is.na(idx)]
  clusterAvg <- data.matrix(clusterAvg[idx,,drop=FALSE])
  
  colorFunc = colorRamp2(c(0, 2, 6), c("light blue", "white", "dark red"))
  invisible(colorFunc(seq(-6, 6)))
  
  plotList[[sample]][['vanGalenHeatmap']] <- Heatmap(clusterAvg,
          na_col = "black",
          col=colorFunc,
          name="Cell Type Markers",
          column_title = "Seurat Clusters",
          cluster_rows = FALSE,
          column_km = 6,
          column_km_repeats = 100,
          row_split = markerCats,
          row_title_rot = 0,
        row_names_gp = grid::gpar(fontsize = 4))

}



## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['vanGalenHeatmap']]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['vanGalenHeatmap']]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['vanGalenHeatmap']]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['vanGalenHeatmap']]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['vanGalenHeatmap']]


## ----cluster_markers----------------------------------------------------------
if (doFindAllMarkers) {
  
  for (sample in names(exptsList)) {
  
    myso <- exptsList[[sample]]
    
    allMarkers <- FindAllMarkers(myso, logfc.threshold = 0.25)
    
    top <- allMarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
    
    plotList[[sample]][['topMarkerHeatmap']] <- DoHeatmap(myso, features = top$gene) + NoLegend()
    
    write.table(allMarkers,
              file = sprintf("%s/%s_top_cluster_markers.tsv",directory$clusterMarkers,sample),
              row.names = FALSE,
              sep = "\t")
  }
  
}


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['topMarkerHeatmap']] + theme(text = element_text(size = 4, face="bold"))


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['topMarkerHeatmap']] + theme(text = element_text(size = 4, face="bold"))


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['topMarkerHeatmap']] + theme(text = element_text(size = 4, face="bold"))


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['topMarkerHeatmap']] + theme(text = element_text(size = 4, face="bold"))


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['topMarkerHeatmap']] + theme(text = element_text(size = 4, face="bold"))


## ----prelimClusterID----------------------------------------------------------

# Assign labels 
labelsList <- list(dmso = c(cluster_0 = "Myeloid Blast", cluster_1 = "Myeloid Blast", cluster_2 = "Monocyte",
                          cluster_3 = "Monocyte", cluster_4 = "T Cell", cluster_5 = "Monocyte",
                          cluster_6 = "GMP", cluster_7 = "Undifferentiated", cluster_8 = "T Cell",
                          cluster_9 = "Monocyte", cluster_10 = "B Cell"),
                  quzi = c(cluster_0 = "Monocyte", cluster_1 = "Monocyte", cluster_2 = "Myeloid Blast",
                          cluster_3 = "T Cell", cluster_4 = "Undifferentiated", cluster_5 = "Monocyte",
                          cluster_6 = "T Cell", cluster_7 = "GMP", cluster_8 = "Monocyte",
                          cluster_9 = "B Cell", cluster_10 = "Monocyte"),
                  ory = c(cluster_0 = "Monocyte", cluster_1 = "Monocyte", cluster_2 = "T Cell",
                          cluster_3 = "Undifferentiated", cluster_4 = "Myeloid Blast", cluster_5 = "GMP",
                          cluster_6 = "Myeloid Blast", cluster_7 = "Monocyte", cluster_8 = "Monocyte",
                          cluster_9 = "Monocyte", cluster_10 = "Monocyte", cluster_11 = "T Cell",
                          cluster_12 = "B Cell", cluster_13 = "Monocyte"),
                  combo = c(cluster_0 = "Monocyte", cluster_1 = "Monocyte", cluster_2 = "Myeloid Blast",
                            cluster_3 = "Undifferentiated", cluster_4 = "T Cell", cluster_5 = "Myeloid Blast",
                            cluster_6 = "GMP", cluster_7 = "Monocyte", cluster_8 = "Monocyte",
                            cluster_9 = "T Cell", cluster_10 = "B Cell", cluster_11 = "Monocyte"),
                  integrated = c(cluster_0 = "Monocyte", cluster_1 = "Myeloid Blast", cluster_2 = "Monocyte",
                            cluster_3 = "Monocyte", cluster_4 = "Myeloid Blast", cluster_5 = "T Cell",
                            cluster_6 = "Undifferentiated", cluster_7 = "GMP", cluster_8 = "Monocyte",
                            cluster_9 = "Monocyte", cluster_10 = "T Cell", cluster_11 = "T Cell",
                            cluster_12 = "B Cell"))

for (sample in names(exptsList)) {
  myso <- exptsList[[sample]]

  Idents(myso) <- "seurat_clusters"
  newClusterIDs <- labelsList[[sample]]
  names(newClusterIDs) <- levels(myso)
  
  myso <- RenameIdents(myso, newClusterIDs)
  myso[['manualClusterIDs']] <- Idents(myso)
  
  plotList[[sample]][['manualClusterAnnotation']] <- DimPlot(myso, reduction = "umap", pt.size = 0.5)
  
  exptsList[[sample]] <- myso
}


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['manualClusterAnnotation']]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['manualClusterAnnotation']]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['manualClusterAnnotation']]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['manualClusterAnnotation']]


## -----------------------------------------------------------------------------
plotList[[names(exptsList)[i]]][['manualClusterAnnotation']]


## ----save_data, include=TRUE--------------------------------------------------
file2save <- sprintf("annotated_clusters.%s.%s.rds", projectName, Sys.Date())
print(sprintf("Saving clustered data of all samples in %s", file2save))
saveRDS(exptsList, file = paste(directory$rda, file2save, sep = "/"))


## ----session_info, include=TRUE-----------------------------------------------
sessionInfo()

