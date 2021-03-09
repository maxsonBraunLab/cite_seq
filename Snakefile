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
		"data/reports/1-preprocessing.html",
		"data/rda/preprocessed.{projectName}.{date}.rds".format(projectName = p, date = time_stamp),
		# 2 - integration
		"data/reports/2-integration.html",
		"data/rda/integrated.{projectName}.{date}.rds".format(projectName = p, date = time_stamp),
		# 3 - cluster
		"data/reports/3-cluster.html",
		"data/rda/cluster.{projectName}.{date}.rds".format(projectName = p, date = time_stamp),
		findAllMarkersOutput

rule preprocessing:
	input:
		samples = expand("data/raw/{sample}/outs/filtered_feature_bc_matrix", sample = SAMPLES)
	output:
		"data/rda/preprocessed.{projectName}.{date}.rds".format(projectName = p, date = time_stamp)
	conda:
		"envs/seurat.yaml"
	script:
		"scripts/1-preprocessing.R"

rule preprocessing_report:
	input:
		samples = expand("data/raw/{sample}/outs/filtered_feature_bc_matrix", sample = SAMPLES)
	output:
		"data/reports/1-preprocessing.html"
	conda:
		"envs/seurat.yaml"
	script:
		"scripts/1-preprocessing.Rmd"

rule integration:
	input:
		rules.preprocessing.output
	output:
		"data/rda/integrated.{projectName}.{date}.rds".format(projectName = p, date = time_stamp)
	conda:
		"envs/seurat.yaml"
	script:
		"scripts/2-integration.R"

rule integration_report:
	input:
		rules.preprocessing.output
	output:
		"data/reports/2-integration.html"
	conda:
		"envs/seurat.yaml"
	script:
		"scripts/2-integration.Rmd"

rule cluster:
	input:
		rules.integration.output
	output:
		"data/rda/cluster.{projectName}.{date}.rds".format(projectName = p, date = time_stamp),
		findAllMarkersOutput
	conda:
		"envs/seurat.yaml"
	threads: config["cores"]
	script:
		"scripts/3-clustering.R"

rule cluster_report:
	input:
		rules.integration.output
	output:
		"data/reports/3-cluster.html"
	conda:
		"envs/seurat.yaml"
	threads: config["cores"]
	script:
		"scripts/3-clustering.Rmd"