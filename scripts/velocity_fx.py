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
import itertools


#### FUNCTIONS ####
def get_top_n_genes_per_cluster(adata,n_genes = 10, cluster = "cluster"):
	"""
	#args: adata  = scanpy or scvelo anndata object
	#	  n_genes = number of genes to select from
	#	  cluster = obs column of the clusters or other categorical data
	#returns a dictionary of the top n genes per cluster according to each genes velocity 
	#keys are the cluster names and values are the list of top genes
	#example usage:
	#	top_genes = get_top_n_genes_per_cluster(adata, n_genes = 10)
	"""
	if "top_cluster_velocity" not in adata.uns:
		scanpy.tl.rank_genes_groups(adata, groupby = cluster,layer = 'velocity', method = "wilcoxon", use_raw = False,key_added = "top_cluster_velocity")
	results = adata.uns['top_cluster_velocity']
	groups = results['names'].dtype.names
	keys = list(results.keys())
	keys = keys[1:]
	all_dicts = {group : pd.DataFrame({"{}".format(key) : results[key][group] for key in keys}) for group in groups}
	for df in all_dicts.keys():
		all_dicts[df][cluster] = df
	top_n_genes_cluster = {key : list(value.nlargest(n_genes,"scores").names) for key,value in all_dicts.items()}
	return top_n_genes_cluster


def cos_sim(adata,cell1,cell2):
	#given two cell names and an object 
	#return the velocity cosine similarity
	
	#cell names to index
	cell1 = adata.obs_names.get_loc(cell1) 
	cell2 = adata.obs_names.get_loc(cell2)
	#velocity vectors for all valid genes (without NaN)
	cell1 = adata.layers['velocity'][cell1,:][~np.isnan(adata.layers['velocity'][cell1,:])]
	cell2 = adata.layers['velocity'][cell2,:][~np.isnan(adata.layers['velocity'][cell2,:])]
	#cosine simularity (inner product)
	return np.dot(cell1, cell2)/(np.linalg.norm(cell1)*np.linalg.norm(cell2))

def add_entropy(adata,dataframe, k = None):
	# takes a results dataframe and adds a 'entropy' column
	# entropy is the mean of cosine similarity of the velocity vectors between cells divided by the distance from
	# cell to each other cell
	# returns dataframe 
	
	#values of the 2 dimensions (velocity and gene count)
	velo_gene_dist = dataframe[['velocity','gene_count']]
	velo_gene_dist.index = dataframe['cell'].values
	#distance matrix of the 2 dimensions per cell
	distance_mat = pd.DataFrame(spatial.distance_matrix(velo_gene_dist.values,velo_gene_dist.values),index = velo_gene_dist.index, columns = velo_gene_dist.index)

	#get all values of cosine similarity divided by the distance between each cell
	cell_key = {}
	
	
	if k is None:
		k = len(distance_mat.index) - 1
	# for each cell in the distance matrix, calulate the cosine similarity for each other cell divided by the distance from the current cell
	for i,cell in enumerate(distance_mat.index):
		#sort the distance vector according to that cell
		sorted_dist = distance_mat.iloc[:,i].sort_values()
		#loop through the sorted vector of distances
		for j,dist in enumerate(sorted_dist):
			cell2 = sorted_dist.index[j] #cell name corresponding to the current distance 
			if cell == cell2:
				continue		   #skip if cells are the same
			elif cell not in cell_key:
				cell_key[cell] = [] #add cell to dictionary if not already in dict
			elif len(cell_key[cell]) == k:
				break			   #gathered results for k cells, break out of loop and onto next cell in outer loop
			else:
				cell_key[cell].append(cos_sim(adata,cell,cell2)/(dist+0.1)) #calculate result and add a little bit to dist in case it is really small
	# convert into a datafame for merging
	entropy = pd.DataFrame([np.mean(value) for key,value in cell_key.items()], index = cell_key.keys(),columns = ["entropy"])
	#return merged dataframe
	return dataframe.merge(entropy, left_on= 'cell', right_index = True)

def scatter_velocity(adata,clust, gene, colorby = 'entropy', groupby = "orig.ident", cluster = "cluster", order = None, save_path = None, cmap = 'RdBu',k_cells = None):
	#Function plot Heatmap of velocity by abundance
	##by gene and cluster
	## if save_path is None (default) then it will not save. otherwise save path should be the path to the directory to save
	#clust = "HSC"  #choose cluster
	#gene = "FOS" #top_n_genes_cluster[clust][12] #choose gene
	#cmap is the color parameter for the plot
	# k_cells is the number of nearest cells to calculater, default is all cells. runs a lot faster with 10 or so cells
	# and still gets the jist
	# can be used to color different parameters with 'colorby'
	if np.isnan(adata.layers['velocity'][:,adata.var_names.get_loc(gene)]).any():
		return "gene {} has no velocity".format(gene)
	
	
	
	#data of what we want to gather to make the plots
	data = {"condition" : [], "velocity" : [], "cluster" : [], "gene" : [], "sample" : [], "cell": [],"gene_count": []}
	if clust is None:
		
		n_len = len(list(adata.obs[groupby]))
		#dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in pd.unique(adata.obs[cluster])}
		data["condition"] = data["condition"] + list(adata.obs[groupby])
		data["cluster"] = data["cluster"] + list(adata.obs[cluster])
		data["velocity"] = data["velocity"] + list(adata.layers['velocity'][:,adata.var_names.get_loc(gene)])
		data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
		data["sample"] = data["sample"] + list(adata.obs['samples'])
		data["cell"] = data["cell"] + list(adata.obs.index)
		data["gene_count"] = data["gene_count"] + list(np.squeeze(np.asarray(adata.X[:,adata.var_names.get_loc(gene)].todense())))
		clust = set(data["cluster"])
	else:
		dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
	#gather results
		for key,value in dict_of_clusters.items():
			n_len = len(list(adata.obs[groupby][value]))
			data["condition"] = data["condition"] + list(adata.obs[groupby][value])
			data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
			data["velocity"] = data["velocity"] + list(adata.layers['velocity'][value,adata.var_names.get_loc(gene)])
			data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
			data["sample"] = data["sample"] + list(adata.obs["samples"][value])
			data["cell"] = data["cell"] + list(adata.obs.index[value])
			data["gene_count"] = data["gene_count"] + list(np.squeeze(np.asarray(adata.X[value,adata.var_names.get_loc(gene)].todense())))
	results = pd.DataFrame(data)
	#return results
	
	results = results.dropna(axis=0) #remove NaN
	cond_dict = {}
	cond_name = []
	vmin,vmax = None,None
	for condition in pd.unique(results['condition']):
		cond_name.append(condition)
		cond_dict[condition] = results[results["condition"] == condition]
		if colorby == 'entropy':
			cond_dict[condition] = add_entropy(adata,cond_dict[condition], k = k_cells)
		if vmax is None:
			vmax = max(cond_dict[condition][colorby])
			vmin = min(cond_dict[condition][colorby])
		else:
			vmax = max(vmax,max(cond_dict[condition][colorby]))
			vmin = min(vmin,min(cond_dict[condition][colorby]))
	#matplotlib setup
	fig, axes = plt.subplots(1,len(cond_dict.keys()), sharey=True)
	
	if order is not None:
		cond_name = order
	#counter
	last = None
	plots_capture = {}
	#plot counts and velocity
	for i,ax in enumerate(axes):
		plots_capture[i] = sns.scatterplot(x="gene_count",y="velocity", data = cond_dict[cond_name[i]], hue = colorby, ax=ax, palette=cmap)
		ax.set_title("{}: {} cells".format(cond_name[i],len(pd.unique(cond_dict[cond_name[i]]['cell']))))
		ax.set_xlabel('')	
		
		last = i
	axes[0].set_ylabel('velocity')
	if colorby == 'entropy':
		fig.suptitle("{}: {}\n color by velocity vector agreement\n mean(cosine similarity / cell distance)".format(clust,gene), y = 1.02,x = 0.58)
		for ax in axes:
			ax.get_legend().remove()
		norm = plt.Normalize(vmin, vmax)
		sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
		sm.set_array([])
		divider = make_axes_locatable(axes[last])
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cbar = fig.colorbar(sm, cax = cax, ticks = [0,vmax])
		cbar.ax.set_yticklabels(['low', 'high'])
	else:
		fig.suptitle("{}: {}\n color by {}".format(clust,gene,colorby), y = 1.02,x = 0.58)
		# put label on right side
		for i,ax in enumerate(axes):
			plots_capture[i].legend(bbox_to_anchor = (1.05,1), loc = 2, borderaxespad=0)
	
	#set x label
	fig.text(0.5, 0, 'abundance (log normalized mRNA counts)', ha='center')
	
	
	#fig.tight_layout()
	if save_path is not None:
		if isinstance(clust,set):
			clust = "all_clusters"
		path_to_save = os.path.join("{}/{}/velocity_{}_by_{}_scatter.png".format(save_path,clust,gene,colorby))
		sys.stderr.write("saving to, {}: ".format(path_to_save))
		#os.makedirs(os.path.dirname(path_to_save),exist_ok = True)
		fig.savefig(path_to_save,bbox_inches='tight')			
	else:
		return fig.tight_layout()

def heatmap_velocity(adata,clust, gene, groupby = "orig.ident", cluster = "cluster" ,save_path = None, nbins = 300, bandwidth = None, plot_max = True,cmap = 'viridis',**kwargs):
	#Function plot Heatmap of velocity by abundance
	##by gene and cluster
	## if save_path is None (default) then it will not save. otherwise save path should be the path to the directory to save
	#clust = "HSC"  #choose cluster
	#gene = "FOS" #top_n_genes_cluster[clust][12] #choose gene
	#nbins is the number to scale the plot (higher number for more bins thus more pixels and smoother image)
	# plot_max (default = True) finds the maximum value for earh plot and caps the scale at the maximum for both plots.  otherwise each plot has its own scale
	# **kwargs any additional arguments for matplotlib pcolormesh
	
	
	"""
	required libraries:
	import scvelo as scv
	import matplotlib.pyplot as plt
	import numpy as np
	import pandas as pd
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from scipy.stats.kde import gaussian_kde
	import matplotlib.colors as colors
	"""
	#import gaussiane kde function
	from scipy.stats.kde import gaussian_kde
	import matplotlib.colors as colors
	
	
	if np.isnan(adata.layers['velocity'][:,adata.var_names.get_loc(gene)]).any():
		return "gene {} has no velocity".format(gene)
	
	  
	
	#get the data from the adata to work
	data = {"condition" : [], "velocity" : [], "cluster" : [], "gene" : [], "sample" : [], "cell": [],"gene_count": []}
	if clust is None:
		#if not specifying clust then plot all clusters
		#usually not very interesting, not recommended
		n_len = len(list(adata.obs[groupby]))
		#dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in pd.unique(adata.obs[cluster])}
		data["condition"] = data["condition"] + list(adata.obs[groupby])
		data["cluster"] = data["cluster"] + list(adata.obs[cluster])
		data["velocity"] = data["velocity"] + list(adata.layers['velocity'][:,adata.var_names.get_loc(gene)])
		data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
		data["sample"] = data["sample"] + list(adata.obs["orig.ident"])
		data["cell"] = data["cell"] + list(adata.obs.index)
		data["gene_count"] = data["gene_count"] + list(np.squeeze(np.asarray(adata.X[:,adata.var_names.get_loc(gene)].todense())))
		clust = set(data["cluster"])
	else:
		#get a dictionary of values (cell barcode locations) for the cluster specified by 'clust'
		dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
	#gather results
		for key,value in dict_of_clusters.items():
			n_len = len(list(adata.obs[groupby][value]))
			data["condition"] = data["condition"] + list(adata.obs[groupby][value])
			data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
			data["velocity"] = data["velocity"] + list(adata.layers['velocity'][value,adata.var_names.get_loc(gene)])
			data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
			data["sample"] = data["sample"] + list(adata.obs["orig.ident"][value])
			data["cell"] = data["cell"] + list(adata.obs.index[value])
			data["gene_count"] = data["gene_count"] + list(np.squeeze(np.asarray(adata.X[value,adata.var_names.get_loc(gene)].todense())))
	#convert results into dataframe for easier plotting
	results = pd.DataFrame(data)
	
	
	results = results.dropna(axis=0) #remove NaN
	#distingiush conditions
	_xmin, _xmax = None,None
	_ymin, _ymax = None,None
	cond_dict = {}
	cond_name = []
	for condition in pd.unique(results['condition']):
		cond_name.append(condition)
		cond_dict[condition] = results[results["condition"] == condition]
		if _xmin is None:
			_xmin,_xmax = np.amin(cond_dict[condition]["gene_count"]),np.amax(cond_dict[condition]["gene_count"])
			_ymin,_ymax = np.amin(cond_dict[condition]["velocity"]),np.amax(cond_dict[condition]["velocity"])
		else:
			_xmin,_xmax = min(_xmin,np.amin(cond_dict[condition]["gene_count"])),max(_xmax,np.amax(cond_dict[condition]["gene_count"]))
			_ymin,_ymax = min(_ymin,np.amin(cond_dict[condition]["velocity"])),max(_ymax,np.amax(cond_dict[condition]["velocity"]))
   
  
	# create a grid starting with minimum and going to maximum separated evenly by the number of nbins
	xi, yi = np.mgrid[_xmin:_xmax:nbins*1j,_ymin:_ymax:nbins*1j]
	
	#min and max scale
	vmin,vmax = None,None
	
	 #get x and y, for each condition
	cond_zi = {}
	for condition,value in cond_dict.items():
		x, y  = value["gene_count"], value["velocity"]
		try:
			k_input = np.vstack([x,y])
			#HD gaussian kernel
			khd = gaussian_kde(k_input, bw_method = bandwidth) #/k_input.std(ddof=1))
			#scaling
			#root the length of x and y * 1j (to make it imaginary)
			#multi dimensional meshgrid
			#apply HD kernel to size of plot
			zi = khd(np.vstack([xi.flatten(), yi.flatten()]))
			cond_zi[condition] = zi
			if vmax is None:
				vmax = max(zi)
				vmin = min(zi)
			else:
				vmax = max(vmax,max(zi))
				vmin = min(vmin,min(zi))
		except:
			return "gene {} has no velocity HD singular matrix".format(gene)
	
	#matplotlib setup
	fig, axes = plt.subplots(1,len(cond_dict.keys()), sharey=True)
	last = None
	max_plot = None
	#plotting
	if plot_max:
		for i,ax in enumerate(axes):
		#plot onto the grid 
			if max(cond_zi[cond_name[i]]) == vmax:
				max_plots = ax.pcolormesh(xi, yi, cond_zi[cond_name[i]].reshape(xi.shape), vmax = cond_zi[cond_name[i]].max(), **kwargs)
			else:
				ax.pcolormesh(xi, yi, cond_zi[cond_name[i]].reshape(xi.shape), vmax = cond_zi[cond_name[i]].max(), **kwargs)
			last = i
			#set labels
			ax.set_title("{}: {} cells".format(cond_name[i],len(pd.unique(cond_dict[cond_name[i]]['cell']))))
			ax.set_xlabel('')
		#Create color bar
		divider = make_axes_locatable(axes[i])
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cbar = fig.colorbar(max_plots, cax=cax, ticks = [vmax/2,vmax])
		cbar.ax.set_yticklabels(['medium', 'high']) #and label them 

	# if not plotting max the plot two color bars
	else:
		for i,ax in enumerate(axes):
			#plot
			cur_plot = ax.pcolormesh(xi, yi, cond_zi[cond_name[i]].reshape(xi.shape), vmax = cond_zi[cond_name[i]].max(), **kwargs)
			#Create color bar
			plt.colorbar(cur_plot,ax = ax)
			#set labels
			ax.set_title("{}: {} cells".format(cond_name[i],len(pd.unique(cond_dict[cond_name[i]]['cell']))))
			ax.set_xlabel('')
		
	#set additional labels
	fig.text(0.5, 0, 'abundance (log normalized mRNA counts)', ha='center')
	fig.suptitle("{}: {} ".format(str(clust),gene), y = 1.02,x = 0.58)
   
	#save if path declared (default None)
	if save_path is not None:
		try:
			if isinstance(clust,set):
				clust = "all_clusters"
			path_to_save = os.path.join("{}/{}/velocity_{}_heatmap.png".format(save_path,clust.replace("/","."),gene))
			os.makedirs(os.path.dirname(path_to_save),exist_ok = True)
			fig.savefig(path_to_save,bbox_inches='tight')
		except:
			sys.stderr.write("did not save, check the save parameter be a string path: {}",format(save_path))
	else:
		#show the figure!
		return fig.tight_layout()

def plot_gene_velocity(adata,gene,groupby = "orig.ident",cluster = "cluster", clust = 'HSC',layer = "velocity", dir_base = "violin_plots", stat = False, show = False, hue_order = None,**kwargs):
	"""
	# seaborn violin plot 
	# if a list of genes: plots the mean of a given layer, x= 'groupy', y = mean(layer), subset by 'clust' in 'cluster' 
	# if a single gene: plots the values of the layer,
	# returns plot if show = True or saves plots in [current working directory]/figures/[dir_base]/[clust]/  
	# 
	# args: 
	#
	# adata  = scanpy or scvelo anndata object
	# gene = gene (str) or list of genes
	# groupby = (str)is the condition splitting the data
	# cluster = (str) obs column of the clusters or other categorical data
	# clust   = (str) a specific entry that is in 'cluster'
	# dir_base = (str) specifies an additional extension for the save location 
	# stat = (bool) perform stats across the groupby, requires hue_order defined in kwargs
	# show = (bool) if true returns plot (default = False)
	# **kwargs to be using in seaborn violinplot
	#returns nothing but saves an png image
	# 
	#example usage:  plot_gene_velocity(adata,gene = "MSI2", groupby = "Condition", cluster = "cluster", layer = "velocity",hue_order = Order_plot, palette = color_dict)
	"""
	if isinstance(gene,str):
		# if gene not in var names, return error
		if (gene not in adata.var_names):
			sys.stderr.write("{} not included in the object\n".format(gene))
			return None
		#dictionary of idx list of clust in cluster
		dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
		data = {groupby : [], layer : []}
		# gather data into a dictionary
		for key,value in dict_of_clusters.items():
			data[groupby] = data[groupby] + list(adata.obs[groupby][value])
			data[layer] = data[layer] + list(adata.layers[layer][value,adata.var_names.get_loc(gene)])
		# convert dictionary to pandas dataframe for ease of plot making
		data_f = pd.DataFrame(data)
		#sns.set(style="whitegrid")
		# plot data
		fig,ax = plt.subplots()
		pplot = sns.violinplot(x = groupby, y = layer, data = data_f,inner = "box", ax = ax,hue_order = hue_order, **kwargs)
		
		if (stat):
			#add statistic for groupby. easy to break if box_pairs parameter is not satisfied properly
			pplot = add_stat_annotation(pplot, data = data_f, x=groupby, y = layer, box_pairs = list(itertools.combinations(hue_order, 2)), test = 't-test_welch', text_format = "full", loc = 'outside')
			pplot = pplot[0]
			#raise plot title to make room for stats
			pplot.set_title("Cluster {}: {}".format(clust,gene), y = 1.2)
		else:
			pplot.set_title("Cluster {}: {}".format(clust,gene))
		if show:
			return pplot
		pplot.get_figure().savefig("figures/{}/violin_genes_{}_{}.png".format(clust,layer,gene),bbox_inches='tight')
		# clear figure to avoid plotting on the same figure
		pplot.get_figure().clf()

	elif isinstance(gene,list):
		genes = gene
		#dictionary of clusters being the keys and list of index locations being the values
		dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
		data = {groupby : [], layer : []}
		included_genes = []
		#determine genes that are valid (in adata object)
		for gene in genes:
			if gene not in adata.var_names:
				sys.stderr.write("{} not included in the object\n".format(gene))
				continue
			else:
				included_genes.append(gene)
		for key,value in dict_of_clusters.items():
			n_len = len(list(adata.obs[groupby][value]))
			data[groupby] = data[groupby] + list(adata.obs[groupby][value])
			#get layer of the cells to calculate
			results = adata.layers[layer][value,:] 
			#get the mean of each cell excluding not a number (NaN)
			results = np.nanmedian(results[:,[adata.var_names.get_loc(gene) for gene in included_genes]], axis = 1)
			data[layer] = data[layer] + list(results)
		data_f = pd.DataFrame(data)
		sns.set(style="whitegrid")
		pplot = sns.violinplot(x = groupby, y = layer, data = data_f,inner = "box", hue_order = hue_order,  **kwargs)
		if (stat):
			#add statistic for groupby. easy to break if box_pairs parameter is not satisfied properly
			pplot = add_stat_annotation(pplot, data = data_f, x=groupby, y = layer, box_pairs = list(itertools.combinations(hue_order, 2)), test = 't-test_welch', text_format = "full", loc = 'outside')
			pplot = pplot[0]
			#raise plot title to make room for stats
			pplot.set_title("{}: top {} DE genes".format(clust,len(genes)), y = 1.2)
		else:
			pplot.set_title("{}: top {} DE genes".format(clust,len(genes)))
		pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 50)
		pplot.set(ylabel = "median({})".format(layer))
		base_dir = "figures/{}/{}".format(dir_base,clust)
		os.makedirs(os.path.join(os.getcwd(),base_dir),exist_ok = True)
		if show:
			return pplot
		else:
			pplot.get_figure().savefig("{}/violin_genes_{}_{}_{}_multigene.png".format(base_dir,groupby,clust,layer),bbox_inches='tight')
			pplot.get_figure().clf()
	else:
		sys.stderr.write("{} not list nor str\n".format(gene))
		return
		
		
def get_velocity_table(adata, name = 'rank_velocity_genes',groupby = None):
	"""
	input: a scvelo or scanpy object
	
	Returns a merged pandas data frame of 'rank_velocity_genes'
	
	Must run 'scv.tl.rank_velocity_genes(adata)' before running
	or else it will run it
	"""
	if name not in adata.uns:
		scv.tl.rank_velocity_genes(adata, groupby=groupby, min_corr=.3, n_genes = adata.n_vars)
	df = scv.DataFrame(adata.uns[name]['names'])
	df_scores = scv.DataFrame(adata.uns[name]['scores'])
	all_dicts = dict()
	for col in df.columns:
		all_dicts[col] = pd.DataFrame()
		all_dicts[col]["genes"] = df[col]
		all_dicts[col]["scores"] = df_scores[col]
	all_list = [value for key,value in all_dicts.items()] 
	df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['genes'],how='inner'), all_list)
	df_merged.columns = ["genes"] + list(df.columns)
	return(df_merged)