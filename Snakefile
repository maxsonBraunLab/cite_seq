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

print("time_stamp:", time_stamp)

p = config["projectName"]

# snakemake -j 64 --use-conda --cluster-config cluster.yaml --profile slurm

rule all:
	input:
		# 1 - preprocessing
		"data/reports/1-preprocessing.html",
		"data/rda/preprocessed.{projectName}.{date}.rds".format(projectName = p, date = time_stamp)

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

# rule integration:
# 	input:
# 		rules.preprocessing.output
# 	output:
# 		"data/rda/integrated.{projectName}.{date}.rds".format(projectName = config["projectName"], date = time_stamp)
# 	conda:
# 		"envs/seurat.yaml"
# 	script:
# 		"scripts/2-preprocessing.R"

# rule integration_report:
# 	input:
# 		rules.preprocessing.output
# 	output:
# 		"data/rda/integrated.{projectName}.{date}.rds".format(projectName = config["projectName"], date = time_stamp)
# 	conda:
# 		"envs/seurat.yaml"
# 	script:
# 		"scripts/2-preprocessing.Rmd"
