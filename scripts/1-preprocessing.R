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
library(ggthemes)
library(rlist)
library(cowplot)
library(viridis)
library(yaml)
library(tibble)

## ----setup_variables, warnings=F----------------------------------------------

# stderr to log file
f <- file("logs/preprocessing.err", "w")
sink(f, type = "message")

# Samples to process named by their treatment conditions
samples2process <- snakemake@config$samples2process

# Name of the single cell project
projectName <- snakemake@config$projectName

# Determines whether or not to regress by cell cycle in scaling step 
cellCycleRegression <- snakemake@config$cellCycleRegression

# List which quantiles to map onto the summary figures
metadataQuants<- unlist(snakemake@config$metadataQuants)

# Set filtering criterion
# Mitochondrial gene expression ratio per cell
percentMitoFilt = snakemake@config$percentMitoFilt

# Minimum and maximum counts per cell
nCountMinFilt = snakemake@config$nCountMinFilt
nCountMaxFilt = snakemake@config$nCountMaxFilt

# Set how many principle components to calculate
nPCs <- snakemake@config$nPCs

# visualize cell cycle states versus expression level of these cell cycle genes.
cc_genes <- snakemake@config$cc_genes

doc_title <- paste(snakemake@config$title, '- preprocessing')
author_list <- paste(snakemake@config$authors, collapse = ", ")


## ----setup_inherent_variables-------------------------------------------------

# List to store Seurat objects for each sample
exptsList <- list()

# List to store CITE-seq antibodies sequenced in each sample
antibodies <- list()

# List to store sample metadata 
metadataList <- list()

# List of data frames to store quantile information for metadata features
metadataQuantList <- list()

# List of data frames to store number of cells within each quantile
metadataQuantCellCountsList <- list()

# Settings list to add to misc slot in Seurat objects
settingsList <- list()

# list of plots to see effects of cell cycle regression
regressionPlots <- list()

# fx to print out metadata features more efficiently in table format
summary_table <- function(samples2process, summVar, indentation, col.names) {
  for (i in c(samples2process)) {
    temp_counts <- cbind(metadataQuantList[[i]][summVar], metadataQuantCellCountsList[[i]][summVar]) %>% rownames_to_column("Percentile")
    cat(paste( paste(rep("#", indentation), collapse = "") , toupper(i), "\n") )
    print(knitr::kable(temp_counts, col.names = col.names))
    cat("\n")
  }
}


## ----cellranger_data----------------------------------------------------------

# Get a list of directories with data

for (sample in samples2process) {
  
  # Load in CellRanger Data
  myData = Read10X(data.dir = paste0("data/raw/", sample, "/outs/filtered_feature_bc_matrix"))
  
  # Gene Expression matrix
  gex <- myData$`Gene Expression`
  
  # Antibody matrix 
  # Replace _ with - to suppress seurat error
  adt <- myData$`Antibody Capture`[,]
  adt_abs <- rownames(adt)
  rownames(adt) <- gsub("_", "-", rownames(adt))
  
  # Save antibodies to a list
  for (antibody in adt_abs) {
    if (!(antibody %in% antibodies)) {
      antibodies <- append(antibodies, antibody)
    }
  }

  # Consider the counts per gene
  print("***********************************************")
  
  # Summarize Features
  print(sprintf("CellRanger output for %s treatment: %i features in %i cells", sample,
                dim(gex)[1], dim(gex)[2]))
  minFeat <- min(Matrix::rowSums(gex))
  print(sprintf("Feature counts range %i to %i", minFeat, max(Matrix::rowSums(gex))))
  print(sprintf("Features with %i counts: %i", minFeat,
                sum(Matrix::rowSums(gex) == minFeat)))
  print(sprintf("Features with 1-3 counts: %i", 
                sum(Matrix::rowSums(gex) > 0 & Matrix::rowSums(gex) < 4)))
  minCells <- min(Matrix::colSums(gex))
  print(sprintf("Cell feature counts range %i to %i", minCells,
                max(Matrix::colSums(gex))))
  print(sprintf("Cells with %i counts: %i", minCells,
                sum(Matrix::colSums(gex) == minCells)))

  # Summarize antibody counts
  paste("Includes the following antibodies:", paste(sapply(rownames(adt), paste), collapse=", "))
  print(sprintf("CellRanger antibody output for %s treatment: %i features in %i cells", sample,
                dim(adt)[1], dim(adt)[2]))
  minAb <- min(Matrix::rowSums(adt))
  print(sprintf("Antibody counts range %i to %i", minAb, max(Matrix::rowSums(adt))))
  print(sprintf("Antibody with %i counts: %i", minAb,
                sum(Matrix::rowSums(adt) == minAb)))
  minAbCells <- min(Matrix::colSums(adt))
  print(sprintf("Cell antibody counts range %i to %i", minAbCells,
                max(Matrix::colSums(adt))))
  print(sprintf("Cells with %i counts: %i", minAbCells,
                sum(Matrix::colSums(adt) == minAbCells)))

  # Create Seurat Object
  # Set thresholds for cell, feature inclusion
  min.cells = 3
  min.features = 200
  settingsList[["step_CreateSeuratObject"]][["min.cells"]] <- min.cells
  settingsList[["step_CreateSeuratObject"]][["min.features"]] <- min.features

  # Initialize the Seurat object with the raw (non-normalized data).
  myso <- CreateSeuratObject(counts = gex, project = sample,
                             min.cells = min.cells, min.features = min.features)
  
  print(sprintf("%s treatment: %i cells in ADT matrix that are not present in gene expression matrix", sample, sum(!(colnames(adt) %in% colnames(gex)))))
  print(sprintf("%s treatment: %i cells in gene expression matrix that are not present in ADT matrix", sample, sum(!(colnames(gex) %in% colnames(adt)))))
  
  myso[["ADT"]] <- CreateAssayObject(counts = adt[, colnames(myso)])
  
  exptsList[[sample]] <- myso
}


## ----metadata-----------------------------------------------------------------

# Convert metadata quantile levels to percent format
percentFormat <- function(x) {
  return(paste(round(100*x, 2), "%", sep=""))
}
metadataQuantPercent <- sapply(metadataQuants, percentFormat)

for (sample in names(exptsList)) {
  myso <- exptsList[[sample]]
  
  # Quantify % mitochondrial features
  myso[["percentMito"]] <- PercentageFeatureSet(myso, pattern = "^MT-")
  
  # Add number of genes per unique molecular identifier (UMI) for each cell to metadata
  myso[['log10GenesPerUMI']] <- log10(myso[['nFeature_RNA']]) / log10(myso[['nCount_RNA']])
  
  # Add number of antibodies per UMI for each cell to metadata
  myso[['log10FeaturesPerUMI_ADT']] <- log10(myso[['nFeature_ADT']]) / log10(myso[['nCount_ADT']])
  
  # Extract metadata
  metadata <- myso@meta.data
  rownames(metadata) <- paste0(sample,'_',rownames(metadata))
  metadata$cells <- rownames(metadata)
  
  # Add to the metadata list
  metadataList[[sample]] <- metadata
  
  # Add to the Seurat object list
  exptsList[[sample]] <- myso
  
  # Calculate the quantiles in the metadata
  metadataNumeric <- which(sapply(metadata, is.numeric))
  metadataQuantList[[sample]] <- data.frame(sapply(metadataNumeric, function(y){
     quantile(x=unlist(metadata[,y]),
              probs = metadataQuants,
              na.rm = TRUE)}))
  
  # Calculate the number of cells in each quantile
  metadataQuantCellCounts <- data.frame(row.names = metadataQuantPercent)
  
  # Loop through the variables in the quantiles list
  for (metadataCol in colnames(metadataQuantList[[sample]])) {
    colVector <- list()
    
    # Calculate number of cells that are within the quantile
    for (level in metadataQuantPercent) {
      elem <- sum(metadata[[metadataCol]] < metadataQuantList[[sample]][level,metadataCol])
      colVector <- list.append(colVector, elem)
    }
    metadataQuantCellCounts[[metadataCol]] <- unlist(colVector)
  }
  
  # Save a column with the total number of cells
  metadataQuantCellCounts[['nCellsTotal']] <- rep(length(metadata$cells),length(metadataQuantPercent))
  
  # Add to the list with the other samples
  metadataQuantCellCountsList[[sample]] <- metadataQuantCellCounts

}

# Integrate the metadata information
integratedMetadata <- bind_rows(metadataList)
  
# Find quantiles of numeric metadata columns
# Mapped as dashed blue lines in the QC figures
integratedMetadataNumeric <- which(sapply(integratedMetadata, is.numeric))
integratedMetadataQuants <- sapply(integratedMetadataNumeric, function(y){
   quantile(x=unlist(integratedMetadata[,y]),
            metadataQuants,
            na.rm = TRUE)
})

## ----filter-------------------------------------------------------------------
preProcExptsList <- list()
for (sample in names(exptsList)) {
  myso <- exptsList[[sample]]
  # Filter high mitochondria cells
  settingsList[["step_filter"]][["percentMitoFilt"]] <- percentMitoFilt
  nMitoFilt <- sum(myso@meta.data$percentMito > percentMitoFilt)
  print(sprintf("%s:  Filter cells > %i%% mitochondria: %i/%i (%.1f%%)", sample, percentMitoFilt, nMitoFilt, ncol(myso), nMitoFilt/ncol(myso)*100))
  # Filter high and low RNA count cells
  settingsList[["step_filter"]][["nCountMinFilt"]] <- nCountMinFilt
  settingsList[["step_filter"]][["nCountMaxFilt"]] <- nCountMaxFilt
  nFeatFilt <- sum(myso@meta.data$nCount_RNA < nCountMinFilt) + sum(myso@meta.data$nCount_RNA > nCountMaxFilt)
  print(sprintf("%s:  Filter cells with %i > nCount > %i: %i/%i (%.1f%%)",
                sample, nCountMinFilt, nCountMaxFilt, nFeatFilt, ncol(myso), nFeatFilt/ncol(myso)*100))
  # Add settings info to the object under "misc"
  myso@misc$settings <- settingsList
  myso <- subset(myso, subset = nCount_RNA > nCountMinFilt & nCount_RNA < nCountMaxFilt & percentMito < percentMitoFilt)
  print(sprintf("Filtered %s:  %i cells, %i features remaining", sample, dim(myso)[1], dim(myso)[1]))
  preProcExptsList[[sample]] <- myso
  rm(myso) # remove seurat object from workspace.
}


## ---- cc_diff_regression------------------------------------------------------
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes

for (sample in samples2process) {
  paste(sprintf("Normalize & Scaling the %s object", sample))
  myso <- preProcExptsList[[sample]]
  # Assign cell cycle scores
  myso <- CellCycleScoring(myso, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
  # Per Seurat "Alternate Workflow" to reduce loss of relevant cell cycle information (eg, during hematopoiesis)
  myso$CC.Difference <- myso$S.Score - myso$G2M.Score
  # Analyze without cell cycle regression
  mysoNoReg <- SCTransform(myso, verbose = FALSE)
  mysoNoReg <- RunPCA(mysoNoReg, npcs = nPCs, verbose = FALSE)
  # Analyze with cell cycle regression
  mysoWithReg <- SCTransform(myso,
                             vars.to.regress = "CC.Difference",
                             verbose = FALSE)
  mysoWithReg <- RunPCA(mysoWithReg, npcs = nPCs, verbose = FALSE)
  # Ridge plot of cell cycle genes
  regressionPlots[[sample]][["ridgePlotCC"]] <- RidgePlot(mysoWithReg, features = cc_genes, group.by = "Phase", ncol = length(snakemake@config$cc_genes) / 2)
  # Plot PCAs Pre- and Post- Cell Cycle Regression
  regressionPlots[[sample]][['pcaPreCCR']] <- DimPlot(mysoNoReg, reduction = "pca", group.by = "Phase", split.by = "Phase")
  regressionPlots[[sample]][['pcaPostCCR']] <- DimPlot(mysoWithReg, reduction = "pca", group.by = "Phase", split.by = "Phase")
  # Save Seurat Object with cell cycle regression if indicated above
  if (cellCycleRegression) {
    preProcExptsList[[sample]] <- mysoWithReg
  } else{
    preProcExptsList[[sample]] <- mysoNoReg
  }
  print(paste("Finished regressing", sample))
  rm(myso)
  rm(mysoWithReg)
  rm(mysoNoReg)
}

## ----save_data----------------------------------------------------------------
print(sprintf("Saving preprocessed data in individual objects in %s", snakemake@output[[1]]))
saveRDS(preProcExptsList, file = snakemake@output[[1]])

## ----session_info-------------------------------------------------------------
sessionInfo()

