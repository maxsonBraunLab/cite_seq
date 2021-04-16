import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
import os
import datetime
import scvelo as scv
import scanpy
import sys
import pandas as pd
from scipy import stats
from statannot import add_stat_annotation

curr_dir = os.getcwd()
begin_time = datetime.datetime.now().timestamp()

sys.stderr.write("beginning analyze!")

in_object         = snakemake.input.in_object

out_object        = snakemake.output.out_object
out_dir           = os.path.dirname(out_object)

genes_of_interest = snakemake.params.genes

condition = snakemake.params.seurat_status

if snakemake.params.color_dict == "None":
    color_dict = None
else:
    color_dict = snakemake.params.color_dict

adata = scv.read(in_object)
os.chdir(out_dir)

#loop through the clusters
#analyze the condition

unique_clusters = pd.unique(adata.obs["cluster"])


keys = 'velocity_length', 'velocity_confidence'

sys.stderr.write("{}".format(condition))
if condition == "None":
	os.chdir(curr_dir)
	adata.obs.to_csv(out_object,sep="\t")
else:
	for clust in unique_clusters:
		adata_subset = adata[adata.obs["cluster"].isin([str(clust)])]
		for gene in genes_of_interest:
			try:	
				scv.pl.velocity(adata_subset,str(gene), dpi = 120, figsize = (7,5), color = condition,legend_loc = 'best',save = "scatter_gene_cluster_{}_{}.png".format(str(clust),gene))
			except:
				sys.stderr.write("{} not included in {}, ".format(gene,str(clust)))
		sys.stderr.write("\n")
		try:
			scv.pl.scatter(adata_subset, c=keys, cmap='coolwarm', perc=[5, 95], save = "scatter_confidence_{}.png".format(str(clust)))
			scv.pl.velocity_embedding_stream(adata_subset, basis='umap', color = condition, palette = color_dict,save = "scvelo_stream_{}.png".format(str(clust)))
			#scv.pl.velocity_embedding(adata_subset, basis='umap', color = 'Condition', palette = color_dict,arrow_length=0, arrow_size=0, dpi=120, save = "scvelo_embedding_{}.png".format(str(clust)))
			scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = condition, save = "velocity_length_{}.png".format(str(clust)), palette = color_dict)
			scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = "orig_ident", save = "orig_velocity_length_{}.png".format(str(clust)), rotation = 90)
		except:
			sys.stderr.write("Error in one of the plots\n")
		
		Count = adata_subset.obs.groupby(condition)['velocity_length'].count()
		Max = adata_subset.obs["velocity_length"].max()
		Means = adata_subset.obs.groupby(condition)["velocity_length"].mean()
		Median = adata_subset.obs.groupby(condition)["velocity_length"].median()
		if color_dict is not None:
			Order_plot = [key for key in color_dict.keys()]
			p_plot = scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = condition, show = False, inner = "box", size = 0, order = Order_plot, palette = color_dict)
		
			fig = add_stat_annotation(p_plot, data = adata_subset.obs, x=condition, y = "velocity_length", box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
		
			add_fig = fig[0]
			for i,x in enumerate(Count):
				add_fig.text(Order_plot.index(Count.index.values[i]),Max+1,"{}".format(x))
			add_fig.get_figure().savefig("figures/stats_velocity_length_{}.png".format(str(clust)))
	

	os.chdir(curr_dir)
	adata.obs.to_csv(out_object,sep="\t")