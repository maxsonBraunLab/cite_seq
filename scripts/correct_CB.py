
import loompy
import numpy as np
import os
import datetime
#import sys

# Hardcoded to enforce that the cellbarcode IDs match the seurat objects and the merged loom file.
# setup to change the basename attached to the cellbarcodes with a colon in velocity such as:
# '[basename(_#)]:[barcode]x' to seurat format '[barcode]_#'
# 'CEL0024_A_1:TTTACGTCAACCAATCx'  -> 'TTTACGTCAACCAATC_1'

## consult the following github source
#https://github.com/velocyto-team/velocyto.py/blob/master/velocyto/analysis.py

#logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
#%matplotlib inline
#plt.rcParams['pdf.fonttype'] = 42

#access snakemake variables
velocity_loom = snakemake.input.velocity_loom
seurat_loom = snakemake.input.seurat_loom
out_file = snakemake.output.out_file
cluster_identity = snakemake.params.seurat_cluster
sample_id = snakemake.params.seurat_sample
sample_id = sample_id.replace(".","_") #remove any periods and replace with underscores
curr_cwd = os.getcwd()

#load loom objects
ds = loompy.connect(seurat_loom, mode ='r') #seurat object to loom ... in r 'as.loom(seuratObj, filename = "seuratObj.loom") 
merged = loompy.connect(velocity_loom, mode = 'r')

#seurat loom meta data
seurat_cells = ds.ca["CellID"]						#Cell ID (barcodes)
umap_coord = ds.ca["umap_cell_embeddings"]		 #umap coordinates  
cluster_ID = ds.ca[cluster_identity]			#cluster id  such as "seurat_clusters" or "integrated_snn_res.0.5"
sample_ids = ds.ca[sample_id]

#lookup the corresponding sample names to the posfix numbers of the seurat object
lookup = dict()
with open('input/paths.tsv') as file_in:
	 for i,line in enumerate(file_in):
			if i == 0:
				continue #skip header line
			else:
				location = line.split("\t")[0]
				uniq	=	str(line.split("\t")[1]).rstrip()
				basename = location.split("/")[-1]
				lookup[basename] = uniq

#make copy of merged loom file
view = merged.view[:,:]
#close and tidy

merged.close()

#make corrected cell barcodes to match seurat objects
view.ca['CellID'] = np.array([s.split(":")[1].replace("x","_") + lookup[s.split(":")[0]] for s in view.ca['CellID']])


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
output_file = os.path.join(curr_cwd,out_file)
loompy.create(output_file,view.layers,view.ra,view.ca)
