import scanpy
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
import os
import datetime
import scvelo as scv
import sys

curr_dir = os.getcwd()
begin_time = datetime.datetime.now().timestamp()


subset_CB = snakemake.params.subset_CB
sys.stderr.write("beginning scvelo on sample: {}!\n".format(subset_CB))
velocity_loom = snakemake.input.subset_CB
genes_of_interest = snakemake.params.genes
out_object = snakemake.output.out_object
out_dir = os.path.dirname(out_object)


#walkthrough
#https://colab.research.google.com/github/theislab/scvelo_notebooks/blob/master/VelocityBasics.ipynb#scrollTo=iHl8jdCUd1j8

#scvelo documentation
#https://readthedocs.org/projects/scvelo/downloads/pdf/latest/

adata = scv.read(velocity_loom)

os.chdir(out_dir)
#matplotlib settings to 'upgraded' images
scv.set_figure_params('scvelo')
scv.logging.print_version()
#proportions of spliced/unspliced counts 
prop_plot = scv.pl.proportions(adata, show = False)


#filter genes, normalize per cell, filter genes dispersion, and scale (log1p)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=None)
#first and second order moments (means and uncentered variances) computed among nearest neighbors in PCA space, computes: pca and neighbors
scv.pp.moments(adata, n_neighbors=30,n_pcs=30)
#default mode for velocity is stochastic,  mode = 'dynamical' and mode = "deterministic" are also available.   see https://scvelo.readthedocs.io/about.html
scv.tl.velocity(adata,  mode = 'dynamical')
#transition probabilties calculated by cosine correlation between the potential cell-to-cell transitions
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'cluster', save = "scvelo_stream.png")
scv.pl.velocity_embedding(adata, basis='umap', color = 'cluster', arrow_length=3, arrow_size=2, dpi=120, save = "scvelo_embedding.png")
scv.pl.velocity_embedding_grid(adata,basis='umap', color = 'cluster', save = "scvelo_grid.png")

# timestamp
plots_time = datetime.datetime.now().timestamp()
sys.stderr.write("finished plots: " + str(round((plots_time-begin_time)/60/60,2)) + " hours\n")

# velocity confidence 
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save = "scatter_confidence.png")


df = adata.obs.groupby('cluster')[keys].mean().T
#df.style.background_gradient(cmap='coolwarm', axis=1)
df.to_csv("velo_confidence_cluster.tsv",sep="\t")


for gene in genes_of_interest:
	try:
		scv.pl.velocity(adata,str(gene), dpi = 120, figsize = (7,5),color = 'cluster',legend_loc = 'best',save = "scatter_gene_{}.png".format(gene))
	except:
		sys.stderr.write("{} not included".format(gene))

almost_time = datetime.datetime.now().timestamp()
sys.stderr.write("almost finished in: " + str(round((almost_time-begin_time)/60/60,2)) + " hours\n")

#save plot proportions
fig = prop_plot.get_figure()
fig.savefig('figures/proportions.png')

#meta_data = scv.get_df(adata, keys=None, layer=None)
adata.obs.to_csv("scvelo_obs.tsv",sep="\t")

os.chdir(curr_dir)

adata.write_h5ad(out_object)



#completed timestamp
end_time = datetime.datetime.now().timestamp()
sys.stderr.write("finished in: " + str(round((end_time - begin_time)/60/60,2)) + " hours")