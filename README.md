# cite_seq

Process CITE-Seq data using Seurat. Steps include pre-processing, individual clustering, sample integration, and group clustering.



---

Developer to-do list:

  * Develop and test the 4 scripts iteratively.
  * Please develop code generically - no hard coding or use analysis-specific files. 
  * Define a common environment (SnakeMake?)
  * Give constructive feedback!

---



# Input

1. Raw output of cellranger should be symlinked into the `data/raw` folder and fit the following structure: `data/raw/{sample}/outs/filtered_feature_bc_matrix/`. The folder should contain the barcodes, features, and matrix files. 
2. Properly configured `config.yaml` file for analysis-specific metadata and data filtration parameters (% mitochondria, nCount_RNA, nCount_ADT). 

# Output

One HTML report for each of the following steps:

## 1 - pre-processing

## 2 - individual Clustering

## 3 - sample Integration

## 4 - group Clustering

## 5 - cluster annotation

## Miscellaneous

​	* Differential genes by cluster table

# Reference

[CITE-Seq](https://www.nature.com/articles/nmeth.4380)