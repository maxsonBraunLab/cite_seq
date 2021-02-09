# cite_seq

Process CITE-Seq data using Seurat. Steps include pre-processing, individual clustering, sample integration, and group clustering.



---

Developer to-do list:

  * Develop and test the 4 scripts iteratively.
  * Define a common environment (SnakeMake?)
  * Give constructive feedback!

---



# Input

1. Raw output of cellranger should be dropped into the `data/raw` folder and fit the following structure: `data/raw/{sample}/filtered_feature_bc_matrix`. The folder should contain the barcodes, features, and matrix files. 
2. Properly configured `config.yaml` file for data filtration (% mitochondria, nCount_RNA, nCount_ADT). 

# Output

## Pre-processing

## Individual Clustering

## Sample Integration

## Group Clustering

## DE genes and RDS files

Describe each sub-header in detail.

# Reference

[CITE-Seq](https://www.nature.com/articles/nmeth.4380)

