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
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats.kde import gaussian_kde
from functools import reduce
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import pearsonr
from scipy import spatial
from velocity_fx import get_top_n_genes_per_cluster,cos_sim,add_entropy,scatter_velocity,heatmap_velocity,plot_gene_velocity,get_velocity_table







def main():

	# set parameters and load objects
	curr_dir = os.getcwd()
	begin_time = datetime.datetime.now().timestamp()
	
	sys.stderr.write("beginning analyze!\n")
	sys.stderr.write("loading parameters and variables\n")
	in_object		 = snakemake.input.in_object
	out_object		= snakemake.output.out_object
	out_dir		   = os.path.dirname(out_object)
	adata_out		 = snakemake.output.adata_out
	genes_of_interest = snakemake.params.genes

	condition = snakemake.params.seurat_status
	
	
	adata = scv.read(in_object)
	Order_plot,Color_hex = None,None
	
	if (snakemake.params.color_hex == "None"):
		Color_hex = sns.color_palette(palette= None,n_colors = len(pd.unique(adata.obs[condition])))
	else:
		Color_hex		 = snakemake.params.color_hex	
	if (snakemake.params.order_plot == "None"):
		Order_plot = list(pd.unique(adata.obs[condition]))
	else:
		Order_plot		= snakemake.params.order_plot

	color_dict		= dict(zip(Order_plot,Color_hex))	
	
	sys.stderr.write("\n color hex: {}\n".format(color_dict))
	os.chdir(out_dir)
	
	#make directory for putting 'tsv' files
	#analysis
	tsv_file_path = os.path.join(os.getcwd(),"analysis")
	os.makedirs(tsv_file_path,exist_ok = True)

	
	#get the unique clusters
	unique_clusters = pd.unique(adata.obs["cluster"])
	
	#get top genes per cluster
	top_n_genes_cluster = get_top_n_genes_per_cluster(adata, n_genes = 5)
	
	

	sys.stderr.write("\n Applied condition:{}\n".format(condition))
	
	#get genes of interest and make sure they are in the object
	genes_of_interest = [gene for gene in genes_of_interest if gene in adata.var_names]
	
	#color for phase plot	
	phase_colors = {'G1':'red', 'G2M':'blue', 'S':'green'}		
	
	
	if condition == "None":
		os.chdir(curr_dir)
		sys.stderr.write("not starting analyses since there is no condition to analyze! check your config file\n")
		adata.obs.to_csv(out_object,sep="\t") 
	else:
		#combine cluster and condition into one meta data column
		adata.obs["cluster_cond"] = ["{}_{}".format(x,y) for x,y in zip(adata.obs["cluster"],adata.obs[condition])]

		#do some global plotting:
		for gene in genes_of_interest:
			scv.pl.velocity(adata,str(gene), dpi = 120, figsize = (7,5), color = "cluster",legend_loc = 'best',save = "scatter_gene_cluster_{}.png".format(gene))		
			scv.pl.velocity(adata,str(gene), dpi = 120, figsize = (7,5), color = "cluster_cond",legend_loc = 'best',save = "scatter_gene_cluster_cond_{}.png".format(gene))
		
		
		
		file_path ="{}/cluster_condition_velocity.tsv".format(tsv_file_path)
		if not os.path.exists(file_path):
			dataf = get_velocity_table(adata,groupby = 'cluster_cond')
			dataf.to_csv(file_path,sep="\t")
		
		#loop through each cluster
		sys.stderr.write("entering outer loop for cluster\n")
		for clust in unique_clusters:
			sys.stderr.write("Analyzing {}".format(clust))
			#make a dir for clust in figures and analyses
			os.makedirs(os.path.join(os.getcwd(),"figures/{}".format(str(clust).replace("/","."))),exist_ok = True)
			os.makedirs(os.path.join(tsv_file_path,"{}".format(str(clust).replace("/","."))),exist_ok=True)
			
			#subset the data to current cluster
			adata_subset = adata[adata.obs["cluster"].isin([str(clust)])]
			
			#for each gene of interest make plots!
			for gene in genes_of_interest:
				sys.stderr.write("On marker: {}".format(gene))
				
				#phase plot
				sys.stderr.write("\n attempting phase plot\n")
				scv.pl.velocity(adata_subset,str(gene), dpi = 120, figsize = (7,5), color = condition,palette = color_dict, legend_loc = 'best',save = "scatter_gene_cluster_{}_{}.png".format(gene,str(clust).replace("/",".")))
				
				
				#phase plot (color by CC phase)
				sys.stderr.write("\n attempting phase plot (phase plot)\n")
				scv.pl.velocity(adata_subset,str(gene), dpi = 120, figsize = (7,5), color = "Phase",legend_loc = 'best',save = "scatter_phase_{}_{}.png".format(gene,str(clust).replace("/",".")), palette = phase_colors)
				
				
				#entropy scatter
				sys.stderr.write("\n attempting orig ident scatter\n")
				#scatter_velocity(adata_subset, clust= str(clust), groupby = condition, color = 'orig.ident',cluster = 'cluster',gene = gene, k_cells = 10, save_path = "./figures")
				
				#heatmap velocity plot
				sys.stderr.write("\n attempting heatmap velocity\n")
				heatmap_velocity(adata_subset, clust = str(clust), gene = gene, nbins = 300, plot_max = True, save_path = "./figures" ,cmap = 'viridis')
				
				#scatter velocity plot
				# velocity or 'layer' on y-axis separated by groupby on x-axis 
				sys.stderr.write("\n attempting scatter velocity plot\n")
				add_args = { 'palette': color_dict}
				plot_gene_velocity(adata,gene = str(gene), groupby = condition, cluster = "cluster", clust = str(clust), layer = "velocity", stat = True, show = False, hue_order = Order_plot,**add_args)
			
			
			#plots for cluster
			sys.stderr.write("\n attempting scatter confidence\n")
			
			scv.pl.scatter(adata_subset, color=['velocity_length', 'velocity_confidence'], cmap='coolwarm', perc=[5, 95], save = "scatter_confidence_{}.png".format(str(clust).replace("/",".")))
			
			scv.pl.velocity_embedding_stream(adata_subset, basis='umap', color = condition, palette = color_dict,save = "scvelo_stream_{}.png".format(str(clust).replace("/",".")))
			
			scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = condition, save = "velocity_length_{}.png".format(str(clust).replace("/",".")), palette = color_dict)
				
			sys.stderr.write("\n attempting velocity length\n")		
			# calculate velocity length (not really interesting but why not)
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
		
			   
			#RANK velocity genes
			# tsv files
			# DE on velocity expression values
			sys.stderr.write("\n attempting DE velocity genes\n")
			scv.tl.rank_velocity_genes(adata_subset, groupby=condition, n_genes = adata_subset.n_vars)
			df = scv.DataFrame(adata_subset.uns['rank_velocity_genes']['names'])
			df_scores = scv.DataFrame(adata_subset.uns['rank_velocity_genes']['scores'])
			
			all_dicts = dict()
			for col in df.columns:
				all_dicts[col] = pd.DataFrame()
				all_dicts[col]["genes"] = df[col]
				all_dicts[col]["scores"] = df_scores[col]
			# after getting the velocity genes, this mangling produces a more readable format
			all_list = [value for key,value in all_dicts.items()] 
			df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['genes'],how='inner'), all_list)
			df_merged.columns = ["genes"] + list(df.columns)
			df_merged.to_csv("top_{}_{}_velocity_genes.tsv".format(clust.replace("/","."),condition) ,sep="\t")
			#save in a table
			df.to_csv("{}/top_{}_velocity_genes.tsv".format(tsv_file_path,clust.replace("/",".")), sep="\t")
			
			
   
			### additional testing specific to two conditions (ie. wild-type vs mutant)
			#adata_subset.layers['velocity'] is the velocity residual for each gene per cell. axis 0 = gene, axis 1 = cell
			valid_var = [np.isnan(adata_subset.layers['velocity'][:,x]).any() for x,i in enumerate(adata_subset.var_names)] #any velocities in which genes have NAN ?  False means there is no NAN and has a velocities for the cells in that feature
			valid_gene_names = [adata_subset.var_names[i] for i,x in enumerate(valid_var) if x == False]
			valid_index = [i for i,x in enumerate(valid_var) if x == False]
			valid_gene_index = dict(zip(valid_gene_names,valid_index)) 
			#list of conditions
			condish = list(adata_subset.obs[condition])
			#make a dictionary to hold values
			df2 = {"sum_values" : [], "mean_values" : [], "condition" : []}
			valid_velocities = adata_subset.layers['velocity'][:,[i for i,x in enumerate(valid_var) if x == False]] # -> [array of cells : array of velocities of valid genes]

			condition_dict = {}
			#fill in the dictionary
			for j,y in enumerate(condish): #for y in in conditions
				#append the mean/sum of the velocities of all valid genes. valid genes being ones that have velocities / spliced/unspliced
				df2["sum_values"].append(sum(valid_velocities[j,:]))
				df2["mean_values"].append(np.mean(valid_velocities[j,:]))
				df2["condition"].append(y)
				if y not in condition_dict:
					condition_dict[y] = [j]
				else:
					condition_dict[y].append(j)
			cond_mean_genes = {}
			
			for key,value in condition_dict.items():
				cond_mean_genes[key] = np.mean(valid_velocities[value,:], axis = 0)  #mean of each genes velocity across the  cells of this cluster
			cond_mean_genes["gene"] = valid_gene_names
			#turn ditioncary into a dataframe
			df = pd.DataFrame(df2)
			pplot = sns.violinplot(x = 'condition', y = 'sum_values', data = df, inner = "box")
			fig2 = add_stat_annotation(pplot, data = df, x='condition', y = "sum_values", box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
			fig2[0].get_figure().savefig("figures/violin_sum_{}_plot.png".format(str(clust).replace("/",".")))
			fig2[0].get_figure().clf()
			pplot.get_figure().clf()
			
			pplot = sns.violinplot(x = 'condition', y = 'mean_values', data = df, inner = "box", palette = color_dict,order = Order_plot)
			fig2 = add_stat_annotation(pplot, data = df, x='condition', y = "sum_values", box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
			fig2[0].get_figure().savefig("figures/violin_mean_{}_plot.png".format(str(clust).replace("/",".")))
			fig2[0].get_figure().clf()
			pplot.get_figure().clf()
			
			df3 = pd.DataFrame(cond_mean_genes)

			try:
				stat,p = stats.f_oneway(valid_velocities[condition_dict[list(condition_dict.keys())[0]],:],valid_velocities[condition_dict[list(condition_dict.keys())[1]],:])
				df3['anova'] = stat
				df3['p_val'] = p
				df3['p_val_Adjusted'] = sm.stats.multipletests(p,method = "fdr_bh")[1]
				
			except:
				sys.stderr.write("Error in Anova\n")
			
			#wilcoxon rank-sum test on velocity values
			scanpy.tl.rank_genes_groups(adata_subset, groupby = condition,layer = 'velocity', method = "wilcoxon", use_raw = False)
			results = adata_subset.uns['rank_genes_groups']
			groups = results['names'].dtype.names #get names of groupby (condition)
			keys = list(results.keys()) #params, names (genes), scores, p_val, ...
			keys = keys[1:] #remove params
			#results for DE for each condition
			results_readable = {group : pd.DataFrame({key:results[key][group] for key in keys}) for group in groups}
			for group in groups:
				results_readable[group].to_csv("{}/{}/wilcoxon_{}_vs_rest_velocity.tsv".format(tsv_file_path,clust.replace("/","."),group), sep="\t")
			
			#wilcoxon rank-sum test on total RNA
			scanpy.tl.rank_genes_groups(adata_subset, groupby = condition, method = "wilcoxon")
			RNA_results = adata_subset.uns['rank_genes_groups']
			RNA_results_readable = {group : pd.DataFrame({key:RNA_results[key][group] for key in keys}) for group in groups}
			
			# merge velocity and total RNA tests
			for group in groups:
				RNA_results_readable[group] = pd.merge(results_readable[group],RNA_results_readable[group],on = ['names'])
				RNA_results_readable[group].columns = [k.replace("x","Vel").replace("y","RNA") for k in list(RNA_results_readable[group].columns)]
				RNA_results_readable[group].to_csv("{}/{}/wilcoxon_{}_RNA_velocity_merged.tsv".format(tsv_file_path,clust.replace("/","."),group), sep="\t")
			# all genes 
			df3.to_csv("{}/{}/anova_mean_velocity_genes.tsv".format(tsv_file_path,clust.replace("/",".")), sep="\t")
			

	############################################## 
	# HEATMAPS SETUP
	##############################################
	order = sorted(np.unique(adata.obs['cluster']))
	
	top_genes = [gene for cluster,gene in top_n_genes_cluster.items()] #list of genes to plot for heatmap
	top_genes = [gene for group in top_genes for gene in group] #flatten list
	adata.obs['clusters'] = adata.obs['cluster'].cat.reorder_categories(list(order), ordered=True)

	################
	#   HEATMAPS   
	################ 
	scanpy.pl.heatmap(adata, var_names = top_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Mu", save = "scanpy_Mu_heatmap.png")
	scanpy.pl.heatmap(adata, var_names = top_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_heatmap.png")

	# latent time
	scv.tl.latent_time(adata)
	scv.pl.scatter(adata, color = 'latent_time',legend_loc = 'on data', color_map = "magma", save = "latent_time.png")
	
	start_cell = np.argmin(adata.obs['latent_time'])
	last_cell = np.argmax(adata.obs['latent_time'])
	
	# pseudotime 
	adata.uns['iroot'] = start_cell
	scanpy.tl.dpt(adata)
	
	#latent time and pseudotime plots
	scv.pl.scatter(adata, color=['dpt_pseudotime','latent_time'], cmap='gnuplot', save = "pseudo_latent_time.png")
	
	#PAGA
	adata.uns['neighbors']['distances'] = adata.obsp['distances']
	adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
	scv.tl.paga(adata, groups='cluster')
	scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "PAGA.png")
	
	#save objects
	os.chdir(curr_dir)
	adata.write_h5ad(adata_out)
	adata.obs.to_csv(out_object,sep="\t")
		
if __name__ == '__main__':
	main()