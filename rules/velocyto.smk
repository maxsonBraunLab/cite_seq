rule samtools_sort:
	input:
		ancient("data/raw/{sample}/outs/possorted_genome_bam.bam")
	output:
		"data/raw/{sample}/outs/cellsorted_possorted_genome_bam.bam"
	conda:
		"../envs/samtools.yml"
	log:
		"logs/velocity/{sample}/samtools_sort.log"
	shell:
		"""
		samtools sort -t CB -O BAM -o {output} {input} > {log} 2>&1
		"""
		
rule velocyto:
	input:
		ancient("data/raw/{sample}/outs/cellsorted_possorted_genome_bam.bam")
	output:
		"data/raw/{sample}/velocyto/{sample}.loom"
	params:
		name = lambda wildcards: "data/raw/{sample}".format(sample = wildcards.sample),
		n_cores=config["cores"],
		gtf = config["GTF_ref"],
		repeat_mask = config["repeat_mask"],
		mem=config["sam_mem"]
	conda:
		"../envs/velocity.yml"
	log:
		"logs/velocity/{sample}/velocyto.log"	
	shell:
		"""
		velocyto run10x -m {params.repeat_mask} -@ {params.n_cores} {params.name} {params.gtf} > {log} 2>&1
		"""
		
rule seurat_convert:
	input:
		seurat_file=rules.cluster.output.rds
	output: 
		out_file="data/looms/seurat_integrated.loom",
		embedding_file = "data/embeddings/seurat_integrated_embeddings.tsv",
		seurat_sample_embedding = expand("data/embeddings/seurat_{sample}_embeddings.tsv", sample = SAMPLES)
	params:
		samples = SAMPLES
	conda:
		"../envs/seurat_convert.yml"
	log:
		"logs/velocity/seurat_convert.log"
	script:
		"../scripts/seurat_convert.R" 
		
rule CB_correct_ind:
	input:
		subset_CB      		 = "data/raw/{sample}/velocyto/{sample}.loom",
		seurat_convert    	 = "data/looms/seurat_integrated.loom",
		seurat_embedding     = "data/embeddings/seurat_integrated_embeddings.tsv",
		seurat_samp_embedding= "data/embeddings/seurat_{sample}_embeddings.tsv"
	output:
		loom_out       		 = "data/velocity/{sample}/velocity.loom"
	params:
		subset_CB      		 = lambda wc:"{}".format(wc.sample),
		seurat_cluster 		 = config["seurat_cluster"],
		seurat_batch   		 = config["seurat_batch"]
	conda:
		"../envs/scvelo_env.yaml"
	log:
		"logs/velocity/{sample}/CB_correct_ind.log"
	script:
		"../scripts/ind_cb_correction.py" 

rule loom_merge:
	input:  
		input_list = expand("data/velocity/{sample}/velocity.loom", sample = SAMPLES)
	output: 
		"data/looms/sorted_merged.loom"
	conda:
		"../envs/velocity.yml"
	log:
		"logs/velocity/loom_merge.log"	
	script:
		"../scripts/merge_looms.py" 
		
rule scvelo_ind:
	input:
		subset_CB				= ancient("data/velocity/{sample}/velocity.loom")
	output:
		out_object				= "data/velocity/{sample}/ind_scvelo_object.h5ad"
	params:
		subset_CB				= lambda wc:"{}".format(wc.sample),
		genes					= config["FeaturePlotList"],
		seurat_cluster			= config["seurat_cluster"],
		seurat_batch 			= config["seurat_batch"]
	conda:
		"../envs/scvelo_env.yaml"
	log:
		"logs/velocity/{sample}/scvelo_ind.log"	
	script:
		"../scripts/scvelo_ind.py" 

rule scvelo_batch:
	input:
		velocity_loom = "data/looms/sorted_merged.loom"
	output: 
		out_object="data/velocity/scvelo_object_batch.h5ad"
	params:
		seurat_cluster=config["seurat_cluster"],
		genes=config["FeaturePlotList"],
		n_cores=config["cores"]
	conda:
		"../envs/scvelo_env.yaml"
	log:
		"logs/velocity/scvelo_batch.log"	
	script:
		"../scripts/scvelo_batch.py" 


rule scvelo:
	input:
		velocity_loom = "data/looms/sorted_merged.loom"
	output: 
		out_object = "data/velocity/scvelo_object.h5ad"
	params:
		genes=config["FeaturePlotList"],
		n_cores=config["cores"]
	conda:
		"../envs/scvelo_env.yaml"
	log:
		"logs/velocity/scvelo.log"
	script:
		"../scripts/scvelo.py"
	

rule analyze_scvelo:
	input:
		in_object="data/velocity/scvelo_object.h5ad"
	output: 
		out_object="data/velocity/scvelo_obs.tsv",
		adata_out ="data/velocity/scvelo_analysis.h5ad"
	params:
		genes=config["FeaturePlotList"],
		seurat_status = config["seurat_status"],
		color_hex = config["color_hex"],
		order_plot = config["order_plot"]
	conda:
		"../envs/scvelo_env.yaml"
	log:
		"logs/velocity/analyze_scvelo.log"
	script:
		"../scripts/Analyze_Cluster_Condition.py" 


rule RNAvelocity:
	input:
		"data/velocity/scvelo_obs.tsv",
		expand("data/velocity/{sample}/ind_scvelo_object.h5ad", sample = SAMPLES),
		"data/velocity/scvelo_object_batch.h5ad",
		"data/velocity/scvelo_analysis.h5ad"
	output: 
		html="data/reports/5-RNAvelocity.html"
	params:
		sample_names = SAMPLES,
		script="scripts/5-RNAvelocity.Rmd",
		seurat=rules.cluster.output.rds,
		seurat_status=config["seurat_status"],
		seurat_cluster=config["seurat_cluster"],
		genes=config["FeaturePlotList"],
		base_line = config["baseline"],
		wave = config['project_name'],
		out_dir="data/velocity"
	conda:
		"../envs/seurat.yaml"
	log:
		"logs/velocity/RNAvelocity.log"
	shell: 
		"""
		Rscript -e 'rmarkdown::render(\"./{params.script}\", output_file = \"../{output.html}\", params=list(inputobs = \"../{input[0]}\", out_dir = \"../{params.out_dir}\", seurat=\"{params.seurat}\",contrast = \"{params.seurat_status}\",cluster = \"{params.seurat_cluster}\",genes=\"{params.genes}\",wave=\"{params.wave}\",samples=\"{params.sample_names}\",baseline=\"{params.base_line}\"))' > {log} 2>&1
		"""











			 
