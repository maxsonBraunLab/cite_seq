import loompy
import numpy as np
import os
import datetime
import pandas as pd
import re


def main():
#merge and subset velocity loom with seurat loom
### intermediate step of making a loom file to match seurat object

    
    #load snakemake and other data
    
    begin_time = datetime.datetime.now().timestamp()
    sys.stderr.write("beginning cb correction!\n")
    seurat_loom = snakemake.input.seurat_convert
    subset_CB = snakemake.params.subset_CB
    velocity_file = snakemake.input.subset_CB
    loom_out = snakemake.output.loom_out
    out_dir = os.path.dirname(loom_out)
    seurat_embedding = snakemake.input.seurat_embedding

    cluster_identity = snakemake.params.seurat_cluster
    #cluster_identity = cluster_identity.replace(".","_")



    
    
    #load loom objects
    ds = loompy.connect(seurat_loom, mode ='r') #seurat object to loom ... in r 'as.loom(seuratObj, filename = "seuratObj.loom") 
    velocity_loom = loompy.connect(velocity_file, mode = 'r')

    #seurat loom meta data
    seurat_cells = ds.ca["CellID"]						#Cell ID (barcodes)
    
    try:
        cluster_ID = ds.ca[cluster_identity]			#cluster id  such as "seurat_clusters" or "integrated_snn_res.0.5"
    except:
        sys.stderr.write("The following meta column is not in data: {}\n ".format(cluster_identity))
    
    
    #sample names
    sample_ids = ds.ca["orig.ident"]

    ##seurat correction of cell barcode and integrated post fix
   
    ## we could find it by looping through variations in the postfix
    ## or/else we guess that it has an underscore so something like this: seurat_correct = "-1_" 
    
    #number of nucleotides in cell barcode
    dummy_cell = ds.ca["CellID"][0]
    number_of_nucleotides = len(re.findall(r'[A|C|G|T]',dummy_cell))
    #get the set of different post fixes
    all_postfix = list(set([x[number_of_nucleotides:] for x in ds.ca["CellID"]]))
    seurat_correct = []
    #find the commonality if more than one different set
    if len(all_postfix) > 1:
        for i,letter in enumerate(all_postfix[0]):
            if letter != all_postfix[1][i]:
                break
            else:
                seurat_correct.append(letter)
        seurat_correct = ''.join([str(elem) for elem in seurat_correct])
    # otherwise make a guess with an underscore
    else:
        seurat_correct = dummy_cell[number_of_nucleotides:].split("_")[0]
        seurat_correct = "{}_".format(seurat_correct)
    ## corresponding number for the seurat sample to the velocity sample
    base_number = dict(zip([y for y in ds.ca['orig.ident']],[x.split(seurat_correct)[1] for x in ds.ca['CellID']]))

    #make new loom (copy of velocity loom file)
    view = velocity_loom.view[:,:]

    #close and tidy
    velocity_loom.close()



    #make corrected cell barcodes to match seurat objects
    view.ca['CellID'] = np.array([s.split(":")[1].replace("x",seurat_correct) + base_number[subset_CB] for s in view.ca['CellID']])

    #filter to keep seurat cells
    view = view.view[:, np.isin(view.ca.CellID, seurat_cells)]

    #umap
    umap_coord = None
    try:
        umap_coord = ds.ca["umap_cell_embeddings"]		 #umap coordinates
    except:
        #read embedding df
        embed_df = pd.read_csv(seurat_embedding, sep='\t', index_col=0, header = None)
        #mangle into array
        umap_coord = np.array([list(x) for x in list(zip(embed_df[1],embed_df[2]))])
            
    
    #add all keys to new loom
    for ca in ds.ca.keys():
        curr_dict = {}
        for i,cb in enumerate(seurat_cells):
            curr_dict[cb] = ds.ca[ca][i]
        view.ca[ca] = np.array([curr_dict[j] for j in view.ca['CellID']])
        
    #close seurat loom    
    ds.close()
    
    # make dictionary for cluster ID
    cell_cluster_dict = {}
    for i in range(len(seurat_cells)):
        cell_cluster_dict[seurat_cells[i]] = cluster_ID[i]

    # dictionary of cell id to umap coord
    cell_umap_dict = {}
    for i in range(len(seurat_cells)):
        cell_umap_dict[seurat_cells[i]] = umap_coord[i]

    # dictionary of cell id to sample ID
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
    
if __name__ == '__main__':
    main()