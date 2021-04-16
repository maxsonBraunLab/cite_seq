import loompy
import numpy as np
import os
import datetime

curr_dir = os.getcwd()
begin_time = datetime.datetime.now().timestamp()
sys.stderr.write("beginning scvelo!")
seurat_loom = snakemake.input.seurat_loom
subset_CB = snakemake.params.subset_CB
velocity_loom = snakemake.input.subset_CB
loom_out = snakemake.output.loom_out
out_dir = os.path.dirname(loom_out)
indices = snakemake.params.indices

sample_id = snakemake.params.seurat_batch
sample_id = sample_id.replace(".","_")
cluster_identity = snakemake.params.seurat_cluster
cluster_identity = cluster_identity.replace(".","_")
#walkthrough
#https://colab.research.google.com/github/theislab/scvelo_notebooks/blob/master/VelocityBasics.ipynb#scrollTo=iHl8jdCUd1j8

#scvelo documentation
#https://readthedocs.org/projects/scvelo/downloads/pdf/latest/


# corresponding number for the seurat sample to the velocity sample
base_number = indices

### intermediate step of making a loom file to match seurat object
#load loom objects
ds = loompy.connect(seurat_loom, mode ='r') #seurat object to loom ... in r 'as.loom(seuratObj, filename = "seuratObj.loom") 
merged = loompy.connect(velocity_loom, mode = 'r')

#seurat loom meta data
seurat_cells = ds.ca["CellID"]						#Cell ID (barcodes)
umap_coord = ds.ca["umap_cell_embeddings"]		 #umap coordinates  
cluster_ID = ds.ca[cluster_identity]			#cluster id  such as "seurat_clusters" or "integrated_snn_res.0.5"
sample_ids = ds.ca[sample_id]
#make copy of merged loom file
view = merged.view[:,:]

#close and tidy

merged.close()
#make corrected cell barcodes to match seurat objects
view.ca['CellID'] = np.array([s.split(":")[1].replace("x","_") + base_number for s in view.ca['CellID']])

#filter to keep seurat cells
view = view.view[:, np.isin(view.ca.CellID, seurat_cells)]

#add all keys
for ca in ds.ca.keys():
	curr_dict = {}
	for i,cb in enumerate(seurat_cells):
		curr_dict[cb] = ds.ca[ca][i]
	view.ca[ca] = np.array([curr_dict[j] for j in view.ca['CellID']])
	
ds.close()
cell_cluster_dict = {}
for i in range(len(seurat_cells)):
	cell_cluster_dict[seurat_cells[i]] = cluster_ID[i]

# dictionary of cell id to umap coord
cell_umap_dict = {}
for i in range(len(seurat_cells)):
	cell_umap_dict[seurat_cells[i]] = umap_coord[i]

# dictionary of cell id to umap coord
cell_sample_dict = {}
for i in range(len(seurat_cells)):
	cell_sample_dict[seurat_cells[i]] = sample_ids[i]

#add meta data array to velocyto loom.
# every cell gets assigned an cluster
view.ca['cluster'] = np.array([cell_cluster_dict[i] for i in view.ca['CellID']])
# every cell gets assigned umap_coordinates
view.ca['umap'] = np.array([cell_umap_dict[i] for i in view.ca['CellID']])
# every cell gets assigned a sample name
view.ca['samples'] = np.array([cell_sample_dict[i] for i in view.ca['CellID']])
#create filtered loom object
loompy.create(loom_out,view.layers,view.ra,view.ca)