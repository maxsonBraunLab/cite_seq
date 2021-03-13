""" Process CITE-Seq data """

""" 1 - preprocessing step includes QC metrics for data filtration """
""" 2 - integration will find anchors and integrate all samples into an integrated seurat object """
""" 3 - clustering will cluster individual samples plus integrated object """
""" 4 - annotate cell clusters by cell type """

import os
import sys
import time

configfile: "config.yaml"

SAMPLES, = glob_wildcards("data/raw/{sample}")

def message(m):
	sys.stderr.write("|--- {} \n".format(m))

for i in SAMPLES:
	message("Samples to process: {}".format(i))

if config["date"]:
	time_stamp = config["date"]
else:
	time_stamp = time.strftime("%Y-%m-%d")

message("time_stamp: {}".format(time_stamp))

# define output of findAllMarkers
findAllMarkersOutput = []
if config["findAllMarkers"]:
	findAllMarkersOutput = expand("data/markers/{sample}-markers.txt", sample = SAMPLES)
	findAllMarkersOutput.append("data/markers/integrated-markers.txt")

p = config["projectName"]

# snakemake -j 8 --use-conda --cluster-config cluster.yaml --profile slurm

rule all:
	input:
		# 1 - preprocessing
		os.getcwd() + "/" + "data/reports/1-preprocessing.html",
		os.getcwd() + "/" + "data/rda/preprocessed.{projectName}.{date}.rds".format(projectName = p, date = time_stamp),
		# 2 - integration
		os.getcwd() + "/" + "data/reports/2-integration.html",
		os.getcwd() + "/" + "data/rda/integrated.{projectName}.{date}.rds".format(projectName = p, date = time_stamp),
		# 3 - cluster
		os.getcwd() + "/" + "data/reports/3-cluster.html",
		os.getcwd() + "/" + "data/rda/cluster.{projectName}.{date}.rds".format(projectName = p, date = time_stamp),
		findAllMarkersOutput

rule preprocessing:
	input:
		rmd = "scripts/1-preprocessing.Rmd",
		samples = expand("data/raw/{sample}/outs/filtered_feature_bc_matrix", sample = SAMPLES)
	output:
		report = os.getcwd() + "/" + "data/reports/1-preprocessing.html",
		rds = os.getcwd() + "/" + "data/rda/preprocessed.{projectName}.{date}.rds".format(projectName = p, date = time_stamp)
	conda:
		"envs/seurat.yaml"
	shell:
		"""
		Rscript -e 'rmarkdown::render( "{input.rmd}", output_file = "{output.report}", knit_root_dir = getwd(), envir = new.env(), params = list(input_samples = "{input.samples}", output_rds = "{output.rds}" ))'
		"""

# knit rmarkdown report with multiple outputs workaround: https://github.com/snakemake/snakemake/issues/178
# how to use rmd params to specify input/output files: https://stackoverflow.com/questions/32479130/passing-parameters-to-r-markdown

rule integration:
	input:
		rmd = "scripts/2-integration.Rmd",
		rds = rules.preprocessing.output.rds
	output:
		report = os.getcwd() + "/" + "data/reports/2-integration.html",
		rds = os.getcwd() + "/" + "data/rda/integrated.{projectName}.{date}.rds".format(projectName = p, date = time_stamp)
	conda:
		"envs/seurat.yaml"
	shell:
		"""
		Rscript -e 'rmarkdown::render( "{input.rmd}", output_file = "{output.report}", knit_root_dir = getwd(), envir = new.env(), params = list(input_rds = "{input.rds}", output_rds = "{output.rds}" ))'
		"""

rule cluster:
	input:
		rmd = "scripts/3-clustering.Rmd",
		rds = rules.integration.output.rds
	output:
		report = os.getcwd() + "/" + "data/reports/3-cluster.html",
		rds = os.getcwd() + "/" + "data/rda/cluster.{projectName}.{date}.rds".format(projectName = p, date = time_stamp),
		# 3d_plots = expand("data/umaps/{sample}_umaps.html", sample = SAMPLES),
		# ab_plots = expand("data/markers/ADT-per-sample-{sample}.png", sample = SAMPLES),
		# de_ident = directory("data/markers/byIdent"),
		de_clust = findAllMarkersOutput
	conda:
		"envs/seurat.yaml"
	shell:
		"""
		Rscript -e 'rmarkdown::render( "{input.rmd}", output_file = "{output.report}", knit_root_dir = getwd(), envir = new.env(), params = list(input_rds = "{input.rds}", output_rds = "{output.rds}" ))'
		"""